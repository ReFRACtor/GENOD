#!/usr/bin/env python

import os, sys, argparse
import netCDF4 as nc
import datetime as DT

# submodule
sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils

# subdir in ABSCO library (https://github.com/ReFRACtor/ABSCO)
# didn't think it was worth it to make this a submodule, because
# this script can be executed once and we have everything we need
ABSCO = '/Users/rpernak/Work/RC/ABSCO'; utils.file_check(ABSCO)
sys.path.append(os.path.join(ABSCO, 'VMR'))
from standard_atm_profiles import vmrProfiles as VMR

# local module
from profile_extraction import readProfiles, singleProfile

class vmrProfO2:
  def __init__(self, inDict):
    """
    Read atmospheric specifications, write pressure ASCII files for
    each profile, select standard atmosphere, then write a CSV with
    O2 profiles in them for the profile pressures
    """

    self.ncFile = inDict['nc_file']
    self.molCSV = inDict['csvMol']
    self.xsCSV = inDict['csvXS']
    inFiles = [self.ncFile, self.molCSV, self.xsCSV]
    for inFile in inFiles: utils.file_check(inFile)

    self.outDir = inDict['out_dir']
    if not os.path.exists(self.outDir): os.makedirs(self.outDir)

  # end constructor

  def readNC(self):
    """
    Read netCDF and keep pressures, latitudes, and times, all of
    which are relevant to standard atmosphere selection
    """

    profiles = readProfiles(self.ncFile)
    self.levP = profiles['level_P']
    self.lats = profiles['lat']
    self.time = [time.split('.')[0] for time in profiles['time']]
    self.timeDT = [DT.datetime.strptime(
      time, '%Y-%m-%dT%H:%M:%S') for time in self.time]

  # end readNC()

  def profP(self):
    """
    This process is necessary for vmrProfiles() -- it needs a list of
    pressures so it can interpolate from standard atmosphere
    pressures onto that grid

    Also construct file names for O2 profiles, which will follow the
    convention used for the pressure profiles
    """

    self.pFiles, self.o2Files = [], []
    for iProf, profP in enumerate(self.levP):
      outFile = '{}_{:d}.txt'.format('UV_pressure_profile', iProf)
      outFile = os.path.join(self.outDir, outFile)
      outFP = open(outFile, 'w')
      for p in profP: outFP.write('{:f}\n'.format(p))
      outFP.close()
      print('Wrote {}'.format(outFile))
      self.pFiles.append(outFile)

      outCSV = outFile.replace('pressure', 'O2')
      outCSV = outCSV.replace('txt', 'csv')
      self.o2Files.append(outCSV)
    # end profile loop
  # end profP

  def selectAtm(self):
    """
    LBLRTM Standard Atmosphere selection based on location and time
    of profile
    """

    # standard atmosphere that should be used in LBL
    # (record 3.1 in LBLRTM instructions)
    # WARNING: not particularly robust, and ad-hoc winter/summer def
    summer = [3, 4, 5, 6, 7, 8]
    winter = [1, 2, 9, 10, 11, 12]

    self.iStdAtm = []
    for lat, time in zip(self.lats, self.timeDT):
      absLat = abs(lat)
      month = time.month
      if absLat <= 15.0:
        iSA = 1
      elif absLat > 15.0 and absLat <= 60:
        iSA = 3 if month in winter else 2
      elif absLat > 60.0:
        iSA = 5 if month in winter else 4
      # endif lat
      self.iStdAtm.append(iSA)
    # end lat/time loop
  # end selectAtm()
# end vmrProfO2

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Generate standard atmosphere O2 profiles for ' + \
      'atmospheric conditions given in netCDF file.')
  parser.add_argument('--nc_file', '-n', type=str, \
    default='uv_benchmark_scenes.nc', help='UV profiles.')
  parser.add_argument('--out_dir', '-d', type=str, \
    default='O2_profiles', \
    help='Directory into which output files will be written.')
  parser.add_argument('--csvMol', '-cm', type=str, \
    default='{}/VMR/LBLATM_Standard_Profiles.csv'.format(ABSCO), \
    help='CSV file that contains LBLATM VMR profile blocks ' + \
    'for all 6 standard atmospheres.  The VMRs exist for ' + \
    'all HITRAN molecules.')
  parser.add_argument('--csvXS', '-cx', type=str, \
    default='{}/VMR/XS_LBLATM_Standard_Profiles.csv'.format(ABSCO), \
    help='CSV file that contains LBLATM VMR profile blocks ' + \
    'for all 6 standard atmospheres.  The VMRs exist for ' + \
    'all XS molecules (i.e., those where no HITRAN line ' + \
    'parameters exist.')
  args = parser.parse_args()

  vObj1 = vmrProfO2(vars(args))
  vObj1.readNC()
  vObj1.profP()
  vObj1.selectAtm()

  for pFile, oFile, i in \
    zip(vObj1.pFiles, vObj1.o2Files, vObj1.iStdAtm):
    vObj2 = VMR(vObj1.molCSV, vObj1.xsCSV, pFile, \
      outFile=oFile, stanAtm=i)
    vObj2.calcVMR()
  # end pFile loop
# end main()
