#!/usr/bin/env python

import os, sys
import numpy as np

sys.path.append('common')
import utils

def readProfiles(inFile, ppmv=False):
  """
  inFile -- string, netCDF file with atmospheric profiles
  ppmv -- boolean, convert VMR profiles to ppmv units
  """

  import netCDF4 as nc

  print('Reading {}'.format(inFile))

  atmDict = {}

  with nc.Dataset(inFile, 'r') as ncObj:
    atm = ncObj.groups['Atmosphere']
    scene = ncObj.groups['Scenario']
    absorb = atm.groups['Absorber']
    ground = atm.groups['Ground']

    # atmospheric specs
    atmDict['molecules'] = nc.chartostring(absorb.variables['name'][:])
    atmDict['VMR'] = absorb.variables['vmr'][:] * 1e6 if ppmv else \
      absorb.variables['vmr'][:]
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
    atmDict['time'] = \
      nc.chartostring(scene.variables['time_string'][:])
  # endwith

  return atmDict
# end readProfiles()

def singleProfile(inDict, iProf, aboveSurface=True, offsetP=0.1):
  """
  Using the output from readProfiles() as a starting point, extract
  a single profile from the collection of profiles provided

  Input
    inDict -- dictionary from readProfiles()
    iProf -- zero-offset index of profile to extract

  Keywords
    aboveSurface -- boolean; only keep levels where pressure exceeds
      the given surface pressure
    offsetP -- float, "perturbation" to surface pressure so that
      pressure levels are between the surface and TOA (so P_O is
      P_surface-offsetP)
  """

  print('Extracting profile {}'.format(iProf+1))

  singleAtm = {}
  for field in inDict: singleAtm[field] = inDict[field][iProf]

  levP = singleAtm['level_P']
  nLev = levP.size

  if aboveSurface:
    iPosAlt = np.where(
      singleAtm['level_P'] < singleAtm['surface_P'])[0]
    if iPosAlt.size == 0:
      print('Found no levels above the surface')

      # need to do something else here
      return singleAtm
    elif iPosAlt.size == nLev-1:
      levP[levP.argmax()] = singleAtm['surface_P']-0.1
    else:
      levP = singleAtm['level_P'][iPosAlt]
    # endif 0

    singleAtm['level_P'] = np.array(levP)
    singleAtm['level_T'] = singleAtm['level_T'][iPosAlt]
    singleAtm['VMR'] = singleAtm['VMR'][:,iPosAlt]
  # end surface

  return singleAtm
# end singleProfile

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Read in JPL-provided netCDF with atmospheric ' + \
    'profiles')
  parser.add_argument('--profiles_nc', '-p', type=str, \
    default='uv_benchmark_scenes.nc', \
    help='netCDF file with atmospheric conditions for a number ' + \
      'of profiles.')
  parser.add_argument('--index', '-i', type=int, \
    help='Zero-offset index for profile of interest.')
  args = parser.parse_args()

  ncFile = args.profiles_nc; utils.file_check(ncFile)
  atm = readProfiles(ncFile)
  iProf = args.index
  if iProf is not None: atm = singleProfile(atm, iProf)

# end main()
