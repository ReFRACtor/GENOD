#!/usr/bin/env python

import os, sys, glob, argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plot

# for multi-page PDF files
from matplotlib.backends.backend_pdf import PdfPages

sys.path.append('common')
import utils

parser = argparse.ArgumentParser(\
  description='Plot OD spectra given netCDF output file from ' + \
    'GENOD_compute.GENOD_netCDF.writeNC().')
parser.add_argument('--netcdf_file', '-i', type=str, \
  default='OD_netCDF/LBLRTM_OD_all_molecules.nc', \
  help='Output file from GENOD_compute.GENOD_netCDF.writeNC(). ' + \
  'Each layer will have its spectrum plotted in a separate page ' + \
  'of a PDF file.')
parser.add_argument('--profile_number', '-p', type=int, default=0, \
  help='Number of profile for which to plot spectra. Profile ' + \
  'numbers corresond to profiles provided in input to ' + \
  'GENOD_compute.RTRefOD class.')
parser.add_argument('--outfile', '-o', type=str, default='OD.pdf', \
  help='Name of PDF that will contain figures.')
args = parser.parse_args()

ncFile = args.netcdf_file; utils.file_check(ncFile)
outFile = args.outfile
iProf = args.profile_number

# grab the data from the netCDF
with nc.Dataset(ncFile, 'r') as ncObj:
  # transform the spectrum for easier looping
  od = ncObj.variables['LBLRTM_Optical_Depth'][:, :, iProf].T
  wn = ncObj.variables['Spectral_Grid'][:]
# endwith

pdf = PdfPages(outFile)

for iLay, layerOD in enumerate(od):
  lay = iLay + 1
  print('Layer {}'.format(lay))

  plot.plot(wn, layerOD)
  plot.xlabel('Wavenumber')
  plot.ylabel(r'$\tau$')
  plot.title('{}, Layer {}'.format(os.path.basename(ncFile), lay))
  pdf.savefig()
  plot.close()
# end odFile loop

pdf.close()

print('Wrote figures to {}'.format(outFile))
