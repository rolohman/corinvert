# corinvert
Code for generating coregistered SLC stack, all possible coherence images, and inverting for factors influencing coherence as a function of time

steps:
1: make a directory called "data" in your working dir, and download all of the Sentinel-1 imagery that you are going to use.
2: Download appropriate dem and run fixImageXml.py with -f option for full pathname
3: in working directory, call stack processor with something like the following command:

stackSentinel.py -s /data/rlohman/Sentinel/Chile/T54/data -o /data/Sentinel/precise/ -a /data/Sentinel/aux_cal/ -w /data/rlohman/Sentinel/Chile/T54/ -d /data/rlo
hman/Sentinel/Chile/dem/demLat_S26_S23_Lon_W071_W067.dem.wgs84 -c 1  -n '1 2' -O 1 -m 20150724 -b '-25.47 -23 -70.6 -69' -r 1 -z 1 -W slc

In that example, I was only using swaths 1 and 2.

4: Run command, perhaps rerunning (after removing dir run_files) with different bounding boxes until the relevant dates are all going to be processed.
5: chmod+x and sequentially run all files in run_files, from working dir. (i.e., run_files/run_1_....)
6: Work on that paper you've been trying to finish, for a day or so.  CHeck in occasionally to make sure things worked.  coreg_pawns should all be around the same size, for instance.
7: Run doVH.py to repeat the resampling of the VH version using the same range and az offsets as in the VV version. This also renames merged/SLC to merged/SLC_VV
8: Manually rename merged/SLC to merged/SLC_VH, and create param.m file in working dir.  I think this can be empty, but I often populate it with pth/masterdate.

%now you should have all merged SLCS.  You can plot values at a pt with plot_all_onfly(xpt,ypt,1), or can run run_intcor_stack.m (change pol in first line) to make the skip 0 and skip 1 int and cor files, plus all vs. first date.

9:invert_dates_stack  -> change pol in first bit to run for vv or vh versions.  this makes the results_dates_VV or _VH dirs and populations with c0, cmin, rel_date and perm_date1_date2 files

%now plot_all_onfly should also show the synthetic cor and timeseries when you run it.
10: geocode all, then look_geo (downlooks by 4x4 to make next step more reasonable)

11: fit_expfun -> this fits exponential decays to time intervals between rain events that are currently hardwired.   When done, can plot with plot_fit_latlon.m
