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
parser.add_argument('--compiler', '-c', type=str, default='ifort', \
  help='Name (not path!) of compiler to use when building LNFL ' + \
  'and LBL. This string is used in a search for the full path ' + \
  'to the executable.')
parser.add_argument('--zenodo_record_id', '-z', type=int, \
  default=3837550, help='Zenodo record ID for AER Line File')
parser.add_argument('--top_dir', '-t', type=str, \
  default=os.getcwd(), help='Top-level directory (used in build)')
parser.add_argument('--molecules', '-m', type=int, nargs='*', \
  default=[27, 2, 6, 22, 24, 26], \
  help='List of zero-offset indices that correspond to molecule ' + \
  'numbers inferred from nc_file Atmosphere.Absorber.name array. ' + \
  'All molecule index options can be displayed with --mol_list')
parser.add_argument('--mol_list', '-ml', action='store_true', \
  help='List all available molecules and their indices to be ' + \
  'used with --molecules argument, then exit. Zero-offset.')
parser.add_argument('--profiles', '-p', type=int, nargs='*', \
  default=[0,1,2,3], \
  help='Zero-offset indices of profiles to model.')
args = parser.parse_args()

ncFile = args.nc_file; utils.file_check(ncFile)
profiles = readProfiles(ncFile, ppmv=True)
nProf = profiles['VMR'].shape[0]

if args.mol_list:
  # it is assumed all molecules have the same molecule names array
  print('{0:<10s}{1:<10s}'.format('Index', 'Molecule'))
  for iMol, mol in enumerate(profiles['molecules'][0]):
    print('{0:<10d}{1:<10s}'.format(iMol, mol))
  print('{0:<10d}{1:<10s}'.format(iMol+1, 'All Molecules'))
  print('Exiting after listing molecule indices')
  sys.exit(0)
# end mol_list

# build models if necessary; start with dictionary that is necessary
# for BUILD objects (see common/build_models.py)
# not all keys have variable values since there is no
# configuration file to edit; "LNFL", "LBLRTM", and "AER_Line_File"
# are submodule paths automated in the clone; we want to build models
buildDict = {'compiler': args.compiler, 'ini': None, \
  'lnfl_path': 'LNFL', 'lblrtm_path': 'LBLRTM', \
  'lines_path': 'AER_Line_File', 'record_id': args.zenodo_record_id, \
  'no_build': False, 'top_dir': args.top_dir}

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
for subset in args.molecules:
  for iProf in range(nProf):
    if iProf not in args.profiles: continue
    profile = singleProfile(profiles, iProf)

    for wn1, wn2 in zip(startWN, endWN):
      # full subset of molecules and XS, then subset/single molecule
      # these indices "keep" the molecule corresponding to the index
      # and is used with profile['VMR']
      odObj = calcOD(ncFile, wn1, wn2, subset, profile, rtAll)
      if onlyNC: continue

      # subset up for the model runs (create inputs)
      odObj.molIdx()
      odObj.xsIdx()
      odObj.profO2(iProf)

      if subset != odObj.nProfMol: odObj.profSubset()

      odObj.lblT5()

      # run the model; the same TAPE3 will be used for all
      # profiles, subsets, and bands
      if not onlyRT and not os.path.exists('TAPE3'): odObj.runLNFL()

      odObj.runLBL()
    # end band loop
  # end profile loop

  totOD = True if subset == odObj.nProfMol else False
  bandArr = np.array((startWN, endWN)).T
  gncObj = GNC('LBL_OD_dir', odObj.subStr, profiles, bandArr, \
    totalOD=totOD)
  gncObj.getProfP()
  gncObj.arrOD()
  gncObj.writeNC()
# end subset loop
