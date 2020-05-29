# GENOD
Similar to ABSCO code, but GENerates an Optical Depth lookup table given a set of atmospheric conditions

# LNFL

I used the [AER-RC Docker Image]() to run LNFL with what I'm calling the "static" LNFL run -- the full spectral range for this project and all of the molecules for which we have line parameters:

```
% docker pull docker.pkg.github.com/aer-rc/lnfl/lnfl:latest
% docker tag docker.pkg.github.com/aer-rc/lnfl/lnfl lnfl
% time docker run --name lnfl --rm -v ~/Work/RC/GENOD/LNFL_TAPE5:/LNFL/TAPE5 -v ~/Work/RC/GENOD/LNFL_out:/LNFL/LNFL_Out lnfl
real	0m43.703s
user	0m0.034s
sys	0m0.014s
% ls -l LNFL_out/TAPE3
-rw-r--r--  1 rpernak  blue  79752 May 22 15:19 LNFL_out/TAPE3
```

# LBLRTM Inputs

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

## O<sub>2</sub>

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
