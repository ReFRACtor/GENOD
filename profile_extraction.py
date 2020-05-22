#!/usr/bin/env python

import os, sys

sys.path.append('common')
import utils

def readProfiles(inFile):
  """
  inFile -- string, netCDF file with atmospheric profiles
  """

  import netCDF4 as nc
  import numpy as np

  atmDict = {}

  with nc.Dataset(inFile, 'r') as ncObj:
    atm = ncObj.groups['Atmosphere']
    scene = ncObj.groups['Scenario']
    absorb = atm.groups['Absorber']
    ground = atm.groups['Ground']

    # atmospheric specs
    atmDict['molecules'] = nc.chartostring(absorb.variables['name'][:])
    atmDict['VMR'] = absorb.variables['vmr'][:]
    atmDict['albedo'] = ground.variables['surface_albedo'][:]

    atmDict['level_P'] = atm.variables['pressure_levels'][:]
    atmDict['surface_P'] = atm.variables['surface_pressure'][:]
    atmDict['level_T'] = atm.variables['temperature_levels'][:]
    atmDict['surface_T'] = atm.variables['surface_temperature'][:]

    # scene specs
    atmDict['lat'] = scene.variables['latitude'][:]
    atmDict['lon'] = scene.variables['longitude'][:]
    atmDict['obs_zenith'] = scene.variables['observation_zenith'][:]
    atmDict['rel_azimuth'] = scene.variables['relative_azimuth'][:]
    atmDict['scatter_angle'] = scene.variables['scattering_angle'][:]
    atmDict['id'] = scene.variables['scene_id'][:]
    atmDict['solar_zenith'] = scene.variables['solar_zenith'][:]
    atmDict['surface_height'] = scene.variables['surface_height'][:]
    atmDict['time'] = scene.variables['time_string'][:]
  # endwith

  return atmDict
# end readProfiles()

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Read in JPL-provided netCDF with atmospheric ' + \
    'profiles')
  parser.add_argument('--profiles_nc', '-i', type=str, \
    default='uv_benchmark_scenes.nc', \
    help='netCDF file with atmospheric conditions for a number ' + \
      'of profiles.')
  args = parser.parse_args()

  ncFile = args.profiles_nc; utils.file_check(ncFile)
  atm = readProfiles(ncFile)
  profMol = atm['molecules'][0]

  inDat = open('molecules.txt').read().splitlines()
  molNames = [mol.split("_")[1] for mol in inDat]

  lnflRec3 = ['0']*47
  for iMol, mol in enumerate(molNames):
    if mol in profMol: lnflRec3[iMol] = '1'
    #print(mol, mol in profMol)

  print(''.join(lnflRec3))
# end main()
