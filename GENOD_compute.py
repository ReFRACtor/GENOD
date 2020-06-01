#!/usr/bin/env python

import os, sys, glob
import numpy as np
import subprocess as sub

sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils

class RTRefOD:
  def __init__(self, inFile, startWN, endWN, subset, profile, \
    exePaths):
    """
    Radiative Transfer Reference Optical Depths

    - Build LBLRTM TAPE5s for each set of molecules of interest
      (full set of molecules and subsets) for a set of 4 profiles
      specified in an input netCDF and over all bands (25000-38000
      cm-1 broken up into 2000 cm-1 chunks for LBLRTM)
    - Build TAPE3 (binary line file) for bands
    - Run LBLRTM to generate ODInt files
    - Use the ODInt files determine path-integrated optical depth and
      store in netCDF file

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
      exePaths -- dictionary with LNFL and LBLRTM executable paths
    """

    self.ncFile = str(inFile)
    self.wn1 = float(startWN)
    self.wn2 = float(endWN)
    self.iSubset = int(subset)
    self.profile = dict(profile)
    self.nLev = profile['VMR'].shape[1]
    self.nProfMol = profile['molecules'].shape[0]
    self.pathLNFL = exePaths['lnfl']
    self.pathLBL = exePaths['lbl']

    # HITRAN stuff; 'molecules.txt' is in version control
    htList = 'molecules.txt'

    paths = [self.ncFile, self.pathLNFL, self.pathLBL, htList]
    for path in paths: utils.file_check(path)

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

    def recordBlock(param, format='{:10.3E}'):
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

      Output
        outRec -- string, format block record
      """

      outRec = ''
      for iVal, val in enumerate(param):
        outRec += format.format(val)

        # eight molecules per line
        if ((iVal+1) % 8) == 0: outRec += '\n'
      # end record36 loop

      # end of record, if nMol is not divisible by 8
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
    record371 = ['{:10s}'.format(xs) for xs in self.iXS.keys()]
    record371 = ''.join(record371)

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

    
  # end runLNFL

  def runLBL(self):
    """
    Run LBLRTM on the TAPE5 generated with lblT5(), then
    save the OD files with their own unique names that contain band,
    subset, and layer
    """

    # stage files needed for these LBL runs
    # dependent on LBLRTM submodule added to repo
    if not os.path.islink('FSCDXS'):
      os.symlink('LBLRTM/cross_sections/FSCDXS', 'FSCDXS')
    if not os.path.islink('xs'):
      os.symlink('LBLRTM/cross_sections/xs', 'xs')
    if not os.path.islink('TAPE3'):
      os.symlink('LNFL_Out/TAPE3', 'TAPE3')

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
  # end runLBL
# end LBLOD
