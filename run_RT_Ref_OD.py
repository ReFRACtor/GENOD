#!/usr/bin/env python

import os, sys, argparse

# Git submodules
sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils
import build_models as BUILD

# local modules
from GENOD_compute import RTRefOD as calcOD
from profile_extraction import readProfiles, singleProfile

# defaults for AER systems
DEFLNFL = '/nas/project/rc/rc1//lblrtm_local_version/' + \
  'lblrtm_v12.9_linux_intel_dbl'
DEFLBL = '/nas/project/rc/rc1/lnfl_local_version/' + \
  'lnfl_v3.2_linux_intel_sgl'

parser = argparse.ArgumentParser(\
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
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
parser.add_argument('--lnfl_exe', '-lnfl', type=str, \
  help='Path to LNFL executable')
parser.add_argument('--lbl_exe', '-lbl', type=str, \
  help='Path to LBLRTM executable')
args = parser.parse_args()

ncFile = args.nc_file; utils.file_check(ncFile)
profiles = readProfiles(ncFile, ppmv=True)
nProf = profiles['VMR'].shape[0]

# build models if necessary; start with dictionary that is necessary
# for BUILD objects (see common/build_models.py)
buildDict = {'compiler': 'ifort', 'ini': None, 'lnfl_path': 'LNFL', \
  'lblrtm_path': 'LBLRTM', 'lines_path': 'AER_Line_File', \
  'record_id': 3837550, 'no_build': False, 'top_dir': os.getcwd()}

exeLNFL = args.lnfl_exe
if exeLNFL is None:
  lnflObj = submodules(vars(args), lnfl=True)
  lnflObj.build()
  exeLNFL = lnflObj.pathLNFL
# endif LNFL build

exeLBL = args.lbl_exe
if exeLBL is None:
  lblObj = submodules(vars(args), lnfl=True)
  lblObj.build()
  exeLBL = lnflObj.pathLBL
# endif LBL build

exeAll = {'lnfl': exeLNFL, 'lbl': exeLBL}
for exe in exeAll.keys(): utils.file_check(exeAll[exe])

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
      odObj = calcOD(ncFile, wn1, wn2, set, profile, exeAll)

      # set up for the model runs (create inputs)
      odObj.molIdx()
      odObj.xsIdx()
      odObj.profO2(iProf)

      if set != odObj.nProfMol: odObj.profSubset()

      odObj.lblT5()

      # run the model; the same TAPE3 will be used for all
      # profiles, subsets, and bands
      if not os.path.exists('TAPE3'): odObj.runLNFL
      odObj.runLBL()
    # end subset loop
  # end band loop
# end profile loop
