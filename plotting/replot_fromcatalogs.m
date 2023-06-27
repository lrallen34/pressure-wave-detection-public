%%% Script for remaking plots from event catalogs, to be used when only
%%% updating plotting function plotevents_fromcatalog_simplified.m

catalogs_loc = '/home/disk/zathras/lukea41/pressure_event_catalogs_coef5_iterative_watersheds_allpairs/**';
datadir = '/home/disk/ivanova2/RPi_Pressure_Data/';
addpath(genpath('../'))

dirstruct = dir(catalogs_loc);
axisflag = 'centered';
minsensors = 3; % minimum number of sensors which must capture an event in order to plot it
coef = 5;
smoothParam = 10;

for curcontent = 1:length(dirstruct)
    if dirstruct(curcontent).isdir
        continue
    end
    
    curpath = [dirstruct(curcontent).folder '/' dirstruct(curcontent).name];
    load(curpath);
    pathcomponents = strsplit(dirstruct(curcontent).folder, '/');
    network = pathcomponents{end};
    threshloc = ['/home/disk/zathras/lukea41/pressure_wavelet_means/' network '_means_asof_20220310.txt'];
    
    sensor_ids = {struct_of_tables.sensor_id};
    name_of_catalog = dirstruct(curcontent).name(1:end-4);
    outdir_plots = ['/home/disk/zathras/lukea41/pressure_crosssensor_plots_' ...
        axisflag 'axes/' name_of_catalog '/'];
    makeneededdirs(outdir_plots)
    
    threshmat = readmatrix(threshloc);
    period_vec = minutes(threshmat(:, 1));
    mean_vec = threshmat(:, 2);
    threshold_vec = mean_vec * coef;

    % Identify events and get the pressure traces from each sensor
    % create cell arrays of the event tables, event traces, datetime vectors,
    % wavelet transforms, and full pressure traces
    wt_struct = struct();
    parfor cur_sensor = 1:length(sensor_ids)
        wt_struct(cur_sensor).sensor_id = sensor_ids{cur_sensor};
        sensor_dtvec = struct_of_tables(cur_sensor).datetime_vecs;
        if ~isempty(sensor_dtvec)
            datetime_span = [sensor_dtvec(1), sensor_dtvec(end)];
            wt_struct(cur_sensor).values = calc_wavelettransform(...
                datadir, datetime_span, sensor_ids{cur_sensor}, period_vec, smoothParam);
        end
    end

    plotevents_fromcatalog_simplified(eventcatalog, struct_of_tables, ...
        wt_struct, threshloc, outdir_plots, minsensors)

    close all
    
end