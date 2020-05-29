#!/usr/bin/env python

import os, sys, argparse

# Git submodules
sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils

# local modules
from GENOD_compute import RTRefOD as calcOD
from profile_extraction import readProfiles, singleProfile

parser = argparse.ArgumentParser(\
  description='Given profile specifications in a netCDF, ' + \
  'calculate optical depths with LBLRTM for different subsets of ' + \
  'molecules in the UV (25000-38000 cm-1).')
parser.add_argument('--nc_file', '-n', type=str, \
  default='uv_benchmark_scenes.nc', \
  help='netCDF file with profile specifications.')
parser.add_argument('--start_wn', '-wn1', type=float, default=25000, \
  help='Starting wavenumber of the entire region of interest')
parser.add_argument('--end_wn', '-wn2', type=float, default=38000, \
  help='Starting wavenumber of the entire region of interest')
args = parser.parse_args()

ncFile = args.nc_file; utils.file_check(ncFile)
profiles = readProfiles(ncFile, ppmv=True)
nProf = profiles['VMR'].shape[0]

# we have to work in 2000 cm-1 chunks for LBLRTM: figure out how many
# bands are needed for the given region
wnChunk = 2000
startWN, endWN = [], []
wn = args.start_wn
regEndWN = args.end_wn
while wn <= regEndWN:
  startWN.append(wn)
  tempWN = wn+wnChunk-1
  if tempWN >= regEndWN: tempWN = regEndWN
  endWN.append(tempWN)
  wn += wnChunk
# end while

# we'll do a separate object per band per profile
for iProf in range(nProf):
  profile = singleProfile(profiles, iProf)
  for wn1, wn2 in zip(startWN, endWN):
    for set in [27, 2, 6, 22, 26]:
      # full set of molecules and XS, then subset/single molecule
      # these indices "keep" the molecule corresponding to the index
      # and is used with profile['VMR']
      odObj = calcOD(ncFile, wn1, wn2, set, profile)
      odObj.molIdx()
      odObj.xsIdx()
      odObj.profO2(iProf)

      if set != odObj.nProfMol: odObj.profSubset()

      odObj.lblT5()
    # end subset loop
  # end band loop
# end profile loop
