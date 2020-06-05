#!/usr/bin/env python

import os, sys, argparse
import numpy as np

# Git submodules
sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils
import build_models as BUILD

# local modules
from GENOD_compute import RTRefOD as calcOD
from GENOD_compute import GENOD_netCDF as GNC
from profile_extraction import readProfiles, singleProfile

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
parser.add_argument('--line_file_dir', '-lines', type=str, \
  help='Top-level directory that contains all AER Line File info.')
parser.add_argument('--lnfl_exe', '-lnfl', type=str, \
  help='Path to LNFL executable')
parser.add_argument('--lbl_exe', '-lbl', type=str, \
  help='Path to LBLRTM executable')
parser.add_argument('--only_lbl', '-rt', action='store_true', \
  help='Forego the LNFL and skip straight to radiative transfer ' + \
  'modeling with LBLRTM (so TAPE3 already exists).')
parser.add_argument('--only_netcdf', '-nc', action='store_true', \
  help='If LBLRTM has already been run and the OD files have ' + \
  'been stored to disk, this option can be used to just generate ' + \
  'the corresponding netCDF file.')
args = parser.parse_args()

ncFile = args.nc_file; utils.file_check(ncFile)
profiles = readProfiles(ncFile, ppmv=True)
nProf = profiles['VMR'].shape[0]

# build models if necessary; start with dictionary that is necessary
# for BUILD objects (see common/build_models.py)
# TODO: need more flexibility with these values
buildDict = {'compiler': 'ifort', 'ini': None, 'lnfl_path': 'LNFL', \
  'lblrtm_path': 'LBLRTM', 'lines_path': 'AER_Line_File', \
  'record_id': 3837550, 'no_build': False, 'top_dir': os.getcwd()}

onlyRT = args.only_lbl
onlyNC = args.only_netcdf
if not onlyRT or not onlyNC:
  exeLNFL = args.lnfl_exe
  if exeLNFL is None:
    lnflObj = BUILD.submodules(buildDict, lnfl=True)
    lnflObj.build()
    exeLNFL = lnflObj.pathLNFL

  # endif LNFL build

  utils.file_check(exeLNFL)

  lfpPath = args.line_file_dir
  if lfpPath is None:
    lfpObj = BUILD.submodules(buildDict, lines=True)
    lfpObj.getLineFile()
    lfpPath = buildDict['lines_path']
  # endif line file

  utils.file_check(buildDict['lines_path'])
else:
  exeLNFL = ''
  lfpPath = ''
# endif only_lbl/only_netcdf

exeLBL = args.lbl_exe
if not onlyNC:
  if exeLBL is None:
    lblObj = BUILD.submodules(buildDict, lbl=True)
    lblObj.build()
    exeLBL = lblObj.pathLBL
  # endif LBL build
  utils.file_check(exeLBL)
else:
  exeLBL = ''
# endif only_netcdf

rtAll = {'lnfl': exeLNFL, 'lbl': exeLBL, 'lines': lfpPath}

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

# we'll do a separate object per subset per band per profile
# TODO: flexibility with the subsets (also list available subsets)
#for set in [27, 2, 6, 22, 24, 26]:
for set in [2]:
  for iProf in range(nProf):
    profile = singleProfile(profiles, iProf)
    for wn1, wn2 in zip(startWN, endWN):
      # full set of molecules and XS, then subset/single molecule
      # these indices "keep" the molecule corresponding to the index
      # and is used with profile['VMR']
      odObj = calcOD(ncFile, wn1, wn2, set, profile, rtAll)
      if onlyNC: continue

      # set up for the model runs (create inputs)
      odObj.molIdx()
      odObj.xsIdx()
      odObj.profO2(iProf)

      if set != odObj.nProfMol: odObj.profSubset()

      odObj.lblT5()

      # run the model; the same TAPE3 will be used for all
      # profiles, subsets, and bands
      if not onlyRT and not os.path.exists('TAPE3'): odObj.runLNFL()

      odObj.runLBL()
    # end band loop
  # end profile loop

  totOD = True if set == odObj.nProfMol else False
  bandArr = np.array((startWN, endWN)).T
  gncObj = GNC('LBL_OD_dir', odObj.subStr, profiles, bandArr, \
    totalOD=totOD)
  gncObj.getProfP()
  gncObj.arrOD()
  gncObj.writeNC()
# end subset loop
