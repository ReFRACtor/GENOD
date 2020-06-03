#!/usr/bin/env python

import os, sys, glob
import numpy as np
import subprocess as sub

sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils

class RTRefOD:
  def __init__(self, inFile, startWN, endWN, subset, profile, \
    rtPaths):
    """
    Radiative Transfer Reference Optical Depths

    - Build LBLRTM TAPE5s for each set of molecules of interest
      (full set of molecules and subsets) for a set of 4 profiles
      specified in an input netCDF and over all bands (25000-38000
      cm-1 broken up into 2000 cm-1 chunks for LBLRTM)
    - Build TAPE3 (binary line file) for bands
    - Run LBLRTM to generate ODInt files

    Generate optical depth tables (OD as a function of
    wavenumber, pressure, temperature, and band) for specified
    molecule set

    Input
      inFile -- string, netCDF file with profile specifications
      startWN -- float, starting wavenumber for band
      endWN -- float, ending wavenumber for band
      subset -- int, represents the zero-offset molecule number from
        JPL profiles (profile['molecules']) of the molecule to "keep",
        with subset=nMol being the full set of molecules
      profile -- dictionary from profile_extraction.singleProfile()
      rtPaths -- dictionary with LNFL and LBLRTM executable paths
        as well as the top-level directory for the AER Line File
    """

    self.ncFile = str(inFile)
    self.wn1 = float(startWN)
    self.wn2 = float(endWN)
    self.iSubset = int(subset)
    self.profile = dict(profile)
    self.nLev = profile['VMR'].shape[1]
    self.nProfMol = profile['molecules'].shape[0]
    self.pathLNFL = rtPaths['lnfl']
    self.pathLBL = rtPaths['lbl']
    self.pathLines = rtPaths['lines']

    # this is in version control
    self.lnflT5 = 'LNFL_TAPE5'
    utils.file_check(self.lnflT5)

    # HITRAN stuff; 'molecules.txt' is in version control
    htList = 'molecules.txt'
    utils.file_check(htList)

    # grab HITRAN metadata
    inDat = open(htList).read().splitlines()
    self.molNamesHT = [mol.split("_")[1] for mol in inDat]
    self.nMolMax = len(self.molNamesHT)

    # HDO is not active in the 25000-38000 cm-1 range
    # we are working on adding HCHO as a XS in the database
    # cannot find CHOHO (is it CHOCHO?)
    self.ignored = ['HDO', 'HCHO', 'CHOHO']

    # these guys are HITRAN molecules AND cross section species,
    # depending on the region. for the UV, they are only XS
    self.molXS = ['SO2', 'NO2']

    if subset == self.nProfMol:
      self.subStr = 'all_molecules'
    elif subset > self.nProfMol:
      print('Invalid molecule subset')
      sys.exit(1)
    else:
      self.subStr = self.profile['molecules'][subset].upper()
    # endif subset
  # end constructor

  def molIdx(self):
    """
    Map the indices of the molecules given in the input profiles to
    their names
    """

    profMol = list(self.profile['molecules'])

    self.iMol = {}
    for iMol, mol in enumerate(self.molNamesHT):
      if mol in self.molXS: continue
      if mol in profMol: self.iMol[mol.upper()] = profMol.index(mol)
    # end molecule loop
  # end molIdx()

  def xsIdx(self):
    """
    Cross section species do not have assigned numbers, so no mapping
    to HITRAN is necessary like the molecules. But this function
    automates the detection of cross section names and stores their
    indices from the profile['molecules'] array.
    """

    self.iXS = {}
    for iMol, mol in enumerate(self.profile['molecules']):
      if (mol not in self.molNamesHT and mol not in self.ignored) or \
          mol in self.molXS: self.iXS[mol.upper()] = iMol
    # end mol loop
    self.nXS = len(self.iXS.keys())
  # end xsIdx()

  def calcAlt(self):
    """
    Height/altitude (km) calculation with Hydrostatic Equation
    z = -(RT/g) * ln(p/p0)
    """

    # constants for Hydrostatic Equation
    R = 287.0 # dry air, J kg-1 K-1; Rydberg Constant
    g = 9.81 # m s-2; acceleration due to gravity

    pLev = self.profile['level_P']
    tLev = self.profile['level_T']
    pArg = pLev / self.profile['surface_P']
    thickness = -(R*tLev/g) * np.log(pArg)

    # because i was getting -0 that messes up LBLRTM
    self.profile['level_Z'] = np.abs(np.cumsum(thickness))/1000.0
  # end calcAlt()

  def profSubset(self):
    """
    Adjust profile VMRs to given subset -- i.e., zero out all of the
    VMRs except for the molecule of interest in a given subset
    """

    # save current VMR profiles
    cache = np.array(self.profile['VMR'])

    # need the array (row) index subset HITRAN molecule or XS species
    try:
      iCache = self.iMol[self.subStr]
    except:
      iCache = self.iXS[self.subStr]
    # end exception

    # overwrite the profile VMRs with just the subset
    self.profile['VMR'] = np.zeros_like(self.profile['VMR'])
    self.profile['VMR'][iCache] = cache[iCache]
  # end profSubset

  def profO2(self, iProfile):
    """
    It is assumed that O2_profiles/* were generated with the default
    configuration of o2_profiles.py and thus follow a specific naming
    convention. these should be in version control.

    NOTE: we are using the altitudes from the O2 profiles
    (standard atmosphere calculation). calcAlt() altitudes are not
    trusted.
    """

    import pandas as PD

    fileO2 = os.path.join(
      'O2_profiles', 'UV_O2_profile_{}.csv'.format(iProfile))
    utils.file_check(fileO2)

    dfO2 = PD.read_csv(fileO2)
    self.levelZ = dfO2['ALT']
    self.profileO2 = dfO2['O2'] * 1e6
    self.profile['molecules'] = \
      np.append((self.profile['molecules']), 'O2')
  # end profO2()

  def lblT5(self):
    """
    Generate an LBLRTM TAPE5 for a given band and subset of molecules

    See lblrtm_instructions.html that come with source code
    """

    def recordBlock(param, format='{:10.3E}', xs=False):
      """
      Records 3.3b, 3.6, and 3.8.2 are flexible in size, depending on
      the number of pressures or molecules being specified. The
      "record blocks" are pretty similarly constructed, so this
      function automates the process, given a parameter over which to
      loop

      Input
        param -- float array of pressures or densities for a profile

      Keywords
        format -- format string for a single value
          https://docs.python.org/3.4/library/string.html#formatexamples
        xs -- boolean; first line in XS profile is a set of 7 names,
          so this keyword takes that into account when constructing
          the record

      Output
        outRec -- string, format block record
      """

      nRecMol = 7 if xs else 8

      outRec = ''
      for iVal, val in enumerate(param):
        outRec += format.format(val)

        # eight molecules per line
        if ((iVal+1) % nRecMol) == 0: outRec += '\n'
      # end record36 loop

      # end of record, if nMol is not divisible by nRecMol
      if outRec[-1] != '\n': outRec += '\n'

      return outRec
    # end recordBlock()

    # organize the TAPE5s by profile timestamp in the working dir
    self.outDirT5 = os.path.join('LBL_TAPE5_dir', self.profile['time'])
    if not os.path.exists(self.outDirT5): os.makedirs(self.outDirT5)

    # problematic for fractional wavenumbers
    outT5 = 'TAPE5_{0:05.0f}-{1:05.0f}_{2:s}'.format(
      self.wn1, self.wn2, self.subStr)
    self.outT5 = os.path.join(self.outDirT5, outT5)

    pLevs = self.profile['level_P']
    zenith = 90-self.profile['obs_zenith']

    # record 1.1: simple description of calculation
    record11 = '$ OD computation for {}'.format(self.outT5)

    # record1.2: HI, F4: spectral line application
    # CN=1: all continua calculated, including Rayleigh extinction
    #   where applicable
    # OD=1, MG=1: optical depth computation, layer-by-layer
    # include lines and cross sections; use LBLATM
    record12 = ' HI=1 F4=1 CN=1 AE=0 EM=0 SC=0 FI=0 PL=0 ' + \
      'TS=0 AM=1 MG=1 LA=0 OD=1 XS=1'
    record12 += '{:>20s}'.format('1')

    # record 1.3 is kinda long...first, band limits
    record13 = '{0:10.3e}{1:10.3e}'.format(self.wn1, self.wn2)

    # concatenate (NOT append) 6 zeros in scientific notation
    # using defaults for SAMPLE, DVSET, ALFAL0, AVMASS,
    # DPTMIN, and DPTFAC params
    record13 += ''.join(['{:10.3e}'.format(0) for i in range(6)])

    # line rejection not recorded and output OD spectral resolution
    record13 += '{0:4s}{1:1d}{2:5s}{3:10.3e}'.format('', 0, '', 0.01)

    # records required with IATM=1 (we're using LBLATM to calculate
    # layer amounts for species -- we only have level amounts)
    # US Standard atmosphere, path type 2 (slant from H1 to H2), 2
    # pressure levels, no zero-filling, full printout, 7 molecules,
    # do not write to TAPE7
    record31 = '{0:5d}{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}{6:5d}'.format(
      0, 2, -self.nLev, 0, 0, self.nMolMax, 0)

    # record 3.2: observer pressure limits, nadir SZA
    record32 = '{0:10.3f}{1:10.3f}{2:10.3f}'.format(
      pLevs.max(), pLevs.min(), zenith)

    # record 3.3b -- list of pressures boundaries for calculations
    record33b = recordBlock(pLevs,format='{:10.3f}')[:-1]

    # record 3.4: user profile header for given molecule
    record34 = '{0:5d}{1:24s}'.format(-self.nLev, ' User profile')

    # record 3.7: number of XS species, user-provided profile
    record37 = '{0:5d}{1:5d}'.format(self.nXS, 0)

    # record 3.7.1: XS molecule name
    record371 = recordBlock(self.iXS.keys(), format='{:10s}', xs=True)
    record371 = record371[:-1]
    #record371 = ['{:10s}'.format(xs) for xs in self.iXS.keys()]
    #record371[70] = '\n'
    #record371 = ''.join(record371)

    # record 3.8: n pressure levels, pressure used for "height"
    record38 = '{0:5d}{1:5d} XS User Profile'.format(self.nLev, 1)

    # generate entire user profile (molecules and cross sections)
    # for given T and write to TAPE5
    records35_36, records381_382 = [], []
    for iLev in range(self.nLev):
      pLev = pLevs[iLev]
      xsAll = [self.profile['VMR'][self.iXS[xs], iLev]
        for xs in self.iXS.keys()]

      record35 = '{0:10.3E}{1:10.3E}{2:10.3E}{3:5s}AA\n'.format(
        self.levelZ[iLev], pLev, self.profile['level_T'][iLev], '')
      records35_36.append(record35)

      # record 3.6: provide VMR at a given level
      lblAll = np.repeat(0.0, self.nMolMax)
      for iHT in range(self.nMolMax):
        nameHT = self.molNamesHT[iHT]
        if nameHT not in self.profile['molecules'] or \
           nameHT in self.molXS: continue

        if iHT == 6:
          lblAll[iHT] = self.profileO2[iLev]
        else:
          lblAll[iHT] = self.profile['VMR'][self.iMol[nameHT], iLev]
        # endif O2

        # sanity check
        #print('{:4e}'.format(lblAll[iHT]), iHT+1, self.iMol[iHT])
      # end molecule loop

      # start building the string for record 3.6
      record36 = recordBlock(np.array(lblAll))
      records35_36.append(record36)

      # record 3.8.1: boundary pressure
      # plus VMR units for all specified XS species
      record381 = '{:<10.3f}'.format(pLev)
      record381 += ''.join(['A']*self.nXS) + '\n'
      records381_382.append(record381)

      # record 3.8.2: layer molecule VMR
      record382 = recordBlock(np.array(xsAll))
      records381_382.append(record382)
    # end level loop

    # strip off extra newlines
    records35_36 = ''.join(records35_36)[:-1]
    records381_382 = ''.join(records381_382)[:-1]

    # combine the record lists into single strings
    records = [record11, record12, record13, \
      record31, record32, record33b, \
      record34, records35_36, record37, record371, \
      record38, records381_382]

    # finally write the TAPE5
    outFP = open(self.outT5, 'w')
    for rec in records: outFP.write('{}\n'.format(rec))
    outFP.write('%%%%')
    outFP.close()
    print('Built {}'.format(self.outT5))
  # end lblT5()

  def runLNFL(self):
    """
    Run LNFL with the version-controlled LNFL_TAPE5 so we can produce
    a TAPE3 necessary for the LBLRTM runs in runLBL()
    """

    if os.path.exists('TAPE5'): os.remove('TAPE5')
    os.symlink(self.lnflT5, 'TAPE5')

    # stage files necessary for LNFL run
    # assuming the directory structure convention we establish in
    # common/build_models.getLineFile()
    lfpFiles = glob.glob(\
      os.path.join(self.pathLines, 'line_file')+'/*')
    for lfpFile, alias in zip(lfpFiles, ['TAPE1', 'TAPE2']):
      if os.path.islink(alias): os.unlink(alias)
      os.symlink(lfpFile, alias)
    # end line file loop

    # broadening and speed dependence parameters
    brdFiles = glob.glob(\
      os.path.join(self.pathLines, 'extra_brd_params') + '/*')
    for brdFile in brdFiles:
      if os.path.islink(alias): os.unlink(alias)
      os.symlink(brdFile, os.path.basename(brdFile))
    # end broadening file loop

    print('Running LNFL')
    sub.call([self.pathLNFL])
  # end runLNFL

  def runLBL(self):
    """
    Run LBLRTM on the TAPE5 generated with lblT5(), then
    save the OD files with their own unique names that contain band,
    subset, and layer
    """

    # stage files needed for these LBL runs
    # dependent on LBLRTM submodule added to repo
    # these likely do not need to be overwritten (assuming LNFL
    # was run and produced a TAPE3 for this project)
    if not os.path.islink('FSCDXS'):
      os.symlink('LBLRTM/cross-sections/FSCDXS', 'FSCDXS')
    if not os.path.islink('xs'):
      os.symlink('LBLRTM/cross-sections/xs', 'xs')
    if not os.path.exists('TAPE3'):
      print('TAPE3 not found -- consider running runLNFL()')

    if os.path.islink('TAPE5'): os.unlink('TAPE5')
    os.symlink(self.outT5, 'TAPE5')

    print('Running LBL')
    sub.call([self.pathLBL])

    # directory for all of the OD files
    self.outDirOD = self.outDirT5.replace('TAPE5', 'OD')
    if not os.path.exists(self.outDirOD): os.makedirs(self.outDirOD)

    odFiles = sorted(glob.glob('ODint_???'))
    for odFile in odFiles:
      outOD = os.path.join(self.outDirOD,
        os.path.basename(self.outT5)).replace('TAPE5', odFile)
      os.rename(odFile, outOD)
    # end OD file loop

    if os.path.exists('TAPE6'):
      outT6 = os.path.basename(self.outT5).replace('TAPE5', 'TAPE6')
      dirT6 = self.outDirOD.replace('OD', 'TAPE6')
      if not os.path.exists(dirT6): os.makedirs(dirT6)
      os.rename('TAPE6', os.path.join(dirT6, outT6))
    # end TAPE6
  # end runLBL
# end RTRefOD

class GENOD_netCDF:
  def __init__(self, odDir, subset, profiles, bands, totalOD=False):
    """
    Generate a netCDF for a given a subset of molecules and a set of
    ODInt files from LBLRTM that were generated with RTRefOD. The
    file will contain an (nLay x nWN x nProf) array of optical
    depths as well as pressure levels and spectral points
    (wavenumbers) used

    'all_molecule' subset ODs are the column-integrated ODs
    other, single-molecule subsets ODs are layer ODs

    Input
      odDir -- string, path to OD files generated with RTRefOD
      subset -- string, indicates what subset of molecules for which
        to create an OD netCDF
      profiles -- dictionary, output from
        profile_extraction.readProfiles()
      bands -- nBands x 2 array of starting and ending wavenumbers
        for each band used in the RTRefOD.runLBL()

    Keywords
      totalOD -- boolean, compute column-integrated optical depth
        at each layer (default is to just store OD for a given layer)
    """

    self.odDir = str(odDir); utils.file_check(odDir)
    self.subStr = str(subset)
    self.profiles = dict(profiles)
    self.bands = np.array(bands)
    self.totalOD = bool(totalOD)

    self.nBands = bands.shape[0]
    self.chunkSize = 5000
    self.compressLev = 4

    # put the netCDF files in a separate subdirectory
    self.ncDir = 'OD_netCDF'
    if not os.path.exists(self.ncDir): os.makedirs(self.ncDir)
  # end constructor

  def getProfP(self):
    """
    Not all of the pressures provided in the user profiles may have
    been necessarily used in the OD calculation, and we need to
    reverse the profiles, so this method extracts the profiles used
    and reverses the direction to TOA-to-Surface
    """

    # local module
    from profile_extraction import singleProfile

    newP = []
    for iProf, profile in enumerate(self.profiles['level_P']):
      profP = singleProfile(self.profiles, iProf)['level_P']
      newP.append(profP[::-1])
    # end profile loop

    self.profP = np.array(newP)

    # dimensions for netCDF file
    profShape = self.profP.shape
    self.nLay = profShape[1]-1
    self.nProf = profShape[0]
  # end getProfP()

  def getFilesOD(self):
    """
    For a molecule subset, gather all LBLRTM OD(int) files and
    organize by profile time
    """

    # there should be one subdirectory per profile time
    subDirs = self.profiles['time']

    self.odFiles = {}
    for iSub, subDir in enumerate(subDirs):
      search = 'ODint_*_*_{}'.format(self.subStr)
      search = os.path.join(self.odDir, subDirs[0], search)
      self.odFiles[subDir] = sorted(glob.glob(search))
    # end subdir loop

    keys = self.odFiles.keys()
    nFileOD = []
    for key in keys: nFileOD.append(len(self.odFiles[key]))

    # for now, all profiles are required to have the same number of
    # of layers; this may not always be the case, so more flexibility
    # here is on the TODO list
    if not np.all(np.array(nFileOD) == self.nLay * self.nBands):
      print('Profiles must have equal number of OD files, exiting')
      sys.exit(1)
    # endif nFileOD
  # end getFilesOD()

  def arrOD(self):
    """
    Read LBLRTM optical depth files and populate them intoan
    (nLay x nWN x nProf) array of ODs
    """

    # Git submodule
    import lblTools

    # spectral points will be the same for all profiles and layers
    wnAll = []

    # array will be constructed by appending arrays onto lists for
    # a given dimension, then converting to an array
    # not particularly efficient
    profDim = []
    for iProf, prof in enumerate(self.profiles['time']):
      layDim = []

      # LBL returns computations from the surface to TOA, but
      # we will calculate total OD from TOA to surface, so we
      # have to work backwards
      for iLay in range(self.nLay, 0, -1):
        bandDim = []
        sLay = str(iLay)

        for wn1, wn2 in self.bands:
          swn1, swn2 = str(wn1), str(wn2)

          odFile = 'ODint_{0:03d}_{1:05.0f}-{2:05.0f}_{3:s}'.format(
            iLay, wn1, wn2, self.subStr)
          odPath = os.path.join(self.odDir, prof, odFile)

          wn, od = lblTools.readOD(odPath, double=True)
          bandDim.extend(list(od))
          if iProf == 0 and iLay == 1: wnAll.extend(list(wn))
          tempStr = 'Storing ODs for {} '.format(self.subStr)
          tempStr += 'Profile {:d}, '.format(iProf+1)
          tempStr += 'Layer {:d}, '.format(iLay)
          tempStr += '{:.0f}-{:.0f} cm-1'.format(wn1, wn2)
          print(tempStr)
          #print(iProf, iLay, wn1, wn2, od.shape, len(bandDim), len(wnAll))
        # end band loop

        layDim.append(bandDim)
      # end layer loop
      profDim.append(layDim)
    # end profile loop

    self.allOD = np.array(profDim)
    self.allWN = np.array(wnAll)
    self.nWN = self.allWN.size

    # integrate OD along column if we're using all molecules
    if self.totalOD: self.allOD = self.allOD.cumsum(axis=1)
  # end arrOD()

  def writeNC(self):
    """
    Write the OD netCDF
    """

    import netCDF4 as nc

    ncFile = 'LBLRTM_OD_{}.nc'.format(self.subStr)
    ncFile = os.path.join(self.ncDir, ncFile)
    print('Building {}'.format(ncFile))

    npzDat = np.load('temp.npz')
    self.allOD = npzDat['od']
    self.allWN = npzDat['wn']
    self.nWN = self.allWN.size

    # switch around the dimensions of the OD array for easier viewing
    # in hdfview in mind
    ncOD = np.transpose(self.allOD, axes=(2,1,0))

    outFP = nc.Dataset(ncFile, 'w')
    outFP.set_fill_on()

    outFP.description = 'Optical depths for {} '.format(self.subStr)
    outFP.description += 'as a function of pressure and wavenumber'

    outFP.source = 'Profile: JPL; Optical Depths: AER'

    dimNames = ['wavenumber', 'layers', 'profiles', 'levels']
    dimVals = [self.nWN, self.nLay, self.nProf, self.nLay+1]
    dimOD = ('wavenumber', 'layers', 'profiles')

    for name, val in zip(dimNames, dimVals):
      outFP.createDimension(name, val)

    # use chunking for the OD array (stolen from ABSCO code)
    inDims = ncOD.shape
    chunksizes = list(inDims)
    chunksizes[0] = min(inDims[0], self.chunkSize)

    # now onto the variables
    # need to reverse the pressure levels because we go
    # TOA-to-surface in arrOD() calculations
    outVar = outFP.createVariable('P_level', float, \
      ('profiles', 'levels'), zlib=True, complevel=self.compressLev, \
      fill_value=np.nan)
    outVar[:] = self.profP[:, ::-1]
    outVar.units = 'mbar'
    outVar.long_name = 'Pressure Levels'
    outVar.valid_range = (0, 1050)
    outVar.description = 'JPL-provided layer boundary pressures'

    # optical depth values will depend on subset of molecules
    # all molecules: column-integrated OD; subsets: layer-OD
    if self.subStr == 'all_molecules':
      odDesc = 'TOA-to-layer path-integrated optical depths'
      odName = 'Path-integrated OD'
    else:
      odDesc = 'Layer-specific optical depths'
      odName = 'Layer OD'
    # endif subStr
    odName = 'LBLRTM-calculated ' + odName

    outVar = outFP.createVariable('LBLRTM_Optical_Depth', float, \
      dimOD, zlib=True, complevel=self.compressLev, \
      chunksizes=chunksizes, fill_value=np.nan)
    outVar[:] = np.ma.array(ncOD, mask=np.isnan(ncOD))
    outVar.units = 'unitless'
    outVar.long_name = odName
    outVar.valid_range = (0, 1e20)
    outVar.description = odDesc

    outVar = outFP.createVariable('Spectral_Grid', float, \
      ('wavenumber'), zlib=True, complevel=self.compressLev, \
      fill_value=np.nan)
    outVar[:] = self.allWN
    outVar.units = 'cm-1'
    outVar.long_name = 'Spectral Points'
    outVar.valid_range = (0, 50000)
    outVar.description = 'Spectral points corresponding to ODs'

    outFP.close()
  # end writeNC()
# end GENOD_netCDF
