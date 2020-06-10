# GENOD <a name='intro'></a>

Objectives for this library are similar to those in [ABSCO](https://github.com/ReFRACtor/ABSCO) code, but instead of absorption coefficients for a flexible pressure and temperature grid, it GENerates an Optical Depth lookup table given a set of atmospheric conditions.

To clone the repository:

```
git clone --recursive git@github.com:pernak18/GENOD.git
```

Note the `--recursive` keyword. This repository is dependent on a number of other repositories and thus needs submodules also cloned. If the `--recursive` keyword is not used, the submodules can still be cloned with:

```
git submodule init
git submodule update
```

# Objective <a name='objective'></a>

For this project, a set of profiles was provided (`uv_benchmark_scenes.nc`) that contained profiles for all molecules of interest in the 25000-38000 cm<sup>-1</sup> spectral range. Optical depths were first calculated for all molecules in the specifications, then the ODs for specific molecules ("subsets") were computed. The "all molecules" run contains cumulative ODs at each layer, starting at the top of the atmosphere (TOA), while individual molecule subsets contain per-layer ODs.

When all molecules are "on", LBLRTM uses the MT_CKD continuum model along with the line parameters that are available in the UV and Visible spectral domains. Water vapor self- and foreign-broadened, carbon dioxide, ozone, oxygen, and nitrogen continua are all selected and not scaled with a multiplicative factor. Rayleigh extinction is also included.

For Ozone-only experiments, its volume mixing ratio profile is used along with its continuum. For all other subsets, no continua are used.

# Dependencies <a name='dependencies'></a>

As mentioned in [the introduction](#intro), there are a number of dependencies associated with this repository, all of which are publicly available (mostly through GitHub):

- [LNFL](https://github.com/AER-RC/LNFL) (v3.2)
- [LBLRTM](https://github.com/AER-RC/LBLRTM) (v12.9)
- [common](https://github.com/AER-RC/common)
- [AER Line File](https://zenodo.org/record/3837550)

The AER Line File is too large for version control in GitHub and is archived in Zenodo.

Necessary Python libraries are mostly standard installs, but some additional packages are required and listed in `requirements.txt`.

# End-to-end GENOD Run <a name='e2e'></a>

`run_RT_Ref_OD.py` is the driver script that is equipped and intended to start from scratch -- i.e., out-of-box after a repository clone. It builds the models, unpacks the AER Line File, runs the models as necessary, and generates the netCDF lookup table.

## Molecule Subsets

Optical depths can be computed for any combination of molecules that are allowed by LBLRTM. What molecules are used is designated in either the LNFL `TAPE5` or LBLRTM `TAPE5`. The latter was chosen for this project -- any molecules that are not used are represented by profiles with zero-valued volume mixing ratios in the `TAPE5`. Subsets are defined with the `--molecules` keyword into `run_RT_Ref_OD.py` and can be any number of molecules. The `--mol_list` keyword can be used to list all of the available molecule indices:

```
./run_RT_Ref_OD.py -ml
Reading uv_benchmark_scenes.nc
Index     Molecule  
0         H2O       
1         CO2       
2         O3        
3         N2O       
4         CO        
5         CH4       
6         SO2       
7         NH3       
8         HNO3      
9         OCS       
10        N2        
11        HCN       
12        SF6       
13        HCOOH     
14        CCL4      
15        CFC11     
16        CFC12     
17        CFC22     
18        HDO       
19        CH3OH     
20        C2H4      
21        PAN       
22        NO2       
23        NO        
24        HCHO      
25        CHOHO     
26        BrO       
27        All Molecules
Exiting after listing molecule indices

./run_RT_Ref_OD.py -m 27 2 4 # all molecules, ozone, and CO experiments
```

## Model Builds <a name='building'></a>

LNFL and LBLRTM executables are built if their arguments (`lnfl_exe` and `lbl_exe`) in `run_RT_Ref_OD.py` are `None`. This is the default setting. They are built using the compiler defined with the `--compiler` keyword into `run_RT_Ref_OD.py` and `common/build_models.py`. The executables are built in their respective submodule directories -- `LNFL` and `LBLRTM`.

## Line File Extraction <a name='linefile'></a>

Similarly, the AER Line File is "built" with `common/build_models.py` -- the tarball is downloaded from Zenodo and extracted into the `AER_Line_File` subdirectory.

## Model Runs <a name='running'></a>

Once the models are built, LNFL is run once to produce a `TAPE3` (see [LNFL section](#lnfl) for some more details). The same binary line file is used for all iterations of LBLRTM.

LBLRTM is run for all molecule subsets, profiles, and bands (2000 cm<sup>-1</sup> chunks -- LBLRTM can only be run 2000 cm<sup>-1</sup> at a time). New `TAPE5` LBLRTM inputs are thus needed for every run. This is automated in `run_RT_Ref_OD.py`. A separate object is created for each iteration.

## Output <a name='output'></a>

Every LBLRTM run produces `ODint` binary files that contain optical depth spectra for every layer in the profile when the level-0 pressure is greater than the surface pressure that is provided. In this case, that is 63 layers. All `ODint` files are uniquely named -- they contain layer number, band, and molecule -- and are stored in a subdirectory for a given profile. This directory structure is assumed in generating the netCDF file that consolidates all profiles, wavenumbers, and layers into a single file.

# Standalone Runs <a name='standalone'></a>

`GENOD_compute.py` is compartmentalized such that LNFL, LBLRTM, and netCDF generation processes can be run separately, which is particularly useful when only a new netCDF (e.g., with different metadata) is desired and new radiative transfer calculations are not.

The rest of this section will be fleshed out when we run our end-to-end tests.

## Profile Extraction <a name='profiles'></a>

Profile data are stored into memory from the input netCDF file via the `profile_extraction.py` module.

## LNFL <a name='lnfl'></a>

Only one run of LNFL is necessary -- the "static" run. It should cover the entire spectral range or interest and all the molecules that will be used in any analysis. It will produce a `TAPE3` binary line parameter file, which is all that is needed. `GENOD_compute.py` contains a `runLNFL()` method that stages files for LNFL and runs it to produce the `TAPE3` then moves it to where it is needed for LBLRTM runs. Alternatively, the [AER-RC Docker Image](https://github.com/AER-RC/LNFL/packages/200491?version=v3.2) can be used to run LNFL with the static run:

```
% docker pull docker.pkg.github.com/aer-rc/lnfl/lnfl:latest
% docker tag docker.pkg.github.com/aer-rc/lnfl/lnfl lnfl
% time docker run --name lnfl --rm -v $GENOD/LNFL_TAPE5:/LNFL/TAPE5 -v $GENOD/LNFL_out:/LNFL/LNFL_Out lnfl
real	0m43.703s
user	0m0.034s
sys	0m0.014s
% ls -l LNFL_out/TAPE3
-rw-r--r--  1 rpernak  blue  79752 May 22 15:19 LNFL_out/TAPE3
```

If Docker is used to run LNFL, the `TAPE3` should then either be copied, moved, or linked into the working directory (i.e., the same directory as `run_RT_Ref_OD.py`).

The AER Line File only has spectral line parameters for O<sub>2</sub> in the UV and VIS parts of the spectrum, so the binary line file is not very large, even for the >25 molecules and 13000 wavenumbers in this study, which is why only a static run is necessary.

## LBLRTM <a name='lbl'></a>

```
% ./run_RT_Ref_OD.py -lbl LBLRTM/lblrtm_v12.9_linux_intel_dbl -rt
```

```
% ./run_RT_Ref_OD.py
Reading uv_benchmark_scenes.nc
Extracting profile 1
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_25000-26999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_27000-28999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_29000-30999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_31000-32999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_33000-34999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_35000-36999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:13:30.00010000Z/TAPE5_37000-38000_all_molecules
Extracting profile 2
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_25000-26999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_27000-28999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_29000-30999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_31000-32999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_33000-34999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_35000-36999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T00:56:26.00000000Z/TAPE5_37000-38000_all_molecules
Extracting profile 3
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_25000-26999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_27000-28999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_29000-30999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_31000-32999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_33000-34999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_35000-36999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:25:14.00000000Z/TAPE5_37000-38000_all_molecules
Extracting profile 4
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_25000-26999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_27000-28999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_29000-30999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_31000-32999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_33000-34999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_35000-36999_all_molecules
Building LBL_TAPE5_dir/2016-04-14T01:06:02.00003000Z/TAPE5_37000-38000_all_molecules
```

## Output netCDF <a name='outnc'></a>

If the ODs have been computed and the user only needs an output netCDF, they can use:

```
% ./run_RT_Ref_OD.py -lines AER_Line_File -lnfl LNFL/lnfl_v3.2_linux_intel_sgl -lbl LBLRTM/lblrtm_v12.9_linux_intel_dbl -nc
```

## O<sub>2</sub> <a name='o2'></a>

code is in version control, but it was run locally and its output was archived so that the code does not need to be run unless new O2 profiles are provided. if it does have to be re-run, the path to a cloned ABSCO library is needed.

used script from ABSCO. first ASCII pressure profiles:

```
% pwd
/Users/rpernak/Work/RC/GENOD
% ./o2_profiles.py
Reading uv_benchmark_scenes.nc
Wrote O2_profiles/UV_pressure_profile_0.txt
Wrote O2_profiles/UV_pressure_profile_1.txt
Wrote O2_profiles/UV_pressure_profile_2.txt
Wrote O2_profiles/UV_pressure_profile_3.txt
Wrote O2_profiles/UV_O2_profile_0.csv
Wrote O2_profiles/UV_O2_profile_1.csv
Wrote O2_profiles/UV_O2_profile_2.csv
Wrote O2_profiles/UV_O2_profile_3.csv
```

these files are in version control.
