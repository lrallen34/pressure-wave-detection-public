% script for producing event catalogs
clear; clc;
addpath('calc\', 'cmocean\', 'plotting\', 'pressure-sensor-tools\')
%% Define variables
network = 'canada'; % define the network of sensors to process
% process one month of data at a time; include a 1 day buffer on both sides
% to account for edge effects in wavelet transform
yr = 2023; % define a year
mnth = 5; % define a month
datetime_span = [datetime(yr,mnth,1,0,0,0)-days(1), ...
    datetime(yr,mnth,1,0,0,0)+calmonths(1)+days(1)]; % leave as is; defined based on assigned year and month
datadir = 'path/to/data'; % FILLIN: define the data location
% after running calcthresh.m and outputting the mean wavelet powers, define
% threshloc as that output file
threshloc = 'path/to/calcthresh/output.txt'; % FILLIN

% the mean wavelet power file is then read in
threshmat = readmatrix(threshloc);
period_vec = minutes(threshmat(:, 1));
mean_vec = threshmat(:, 2);

coef = 5; % this is K/2, half the threshold in normalized wavelet power to 
% identify an event center, equal to the threshold in normalized wavelet 
% power to define the outer boundary of the wave event region
crosscorr_thresh = 0.65; % minimum optimized cross-correlation to pair events from two sensors together
maxdelay = 2; % max gap between events in two sensors to possibly pair them together
smoothParam = 10; % length of averaging period in seconds
name_of_catalog = [num2str(yr) num2str(mnth, '%02i') '_coef' ...
    num2str(coef) '_iterative_watersheds_allpairs']; % name of output .mat file

% define output directory for the .mat event catalog file
catalog_outdir = ['path/to/catalogs/' network '/']; %FILLIN
% define output directory for the .xlsx event spreadsheet file
spreadsheet_outdir = ['/path/to/spreadsheets/' network '/']; % FILLIN

% define the location of the file containing sensor names and lat/lon
% coordinates
sensorloc = [network '_sites.txt'];

% define output directory for plots
axisflag = 'centered';
outdir_plots = ['/path/to/plots/' ...
    axisflag 'axes/' name_of_catalog '/']; % FILLIN
minsensors = 3; % minimum number of sensors which must capture an event in order to plot it
makeneededdirs(catalog_outdir, spreadsheet_outdir, outdir_plots)

% read in the sensor list file
sensormat = readmatrix(sensorloc, 'OutputType', 'char');
sensor_ids = sensormat(:, 1);
alt_lats = str2double(sensormat(:, 2));
alt_lons = str2double(sensormat(:, 3));

%% Produce tables of events

[struct_of_tables, wt_struct] = ...
    determinecoherence_fullnetwork_iterative(datadir, datetime_span, ...
    sensor_ids, threshloc, coef, crosscorr_thresh, maxdelay, smoothParam);

%% check for two-way matches and assign IDs to each event

struct_of_tables = assign_eventIDs(struct_of_tables);

%% Go through each unique event ID and calculate the slowness vector

[wave_propertytable,sensors_with_event,timedeltas,modeled_timedeltas,...
    corr_coefs,sensor_pairs] = ...
    calc_slownessvector(struct_of_tables, alt_lats, alt_lons);

%% Output a catalog and spreadsheet of events

catalog_outloc = [catalog_outdir name_of_catalog '.mat'];
spreadsheet_outloc = [spreadsheet_outdir name_of_catalog '.xlsx'];

[eventcatalog, eventtable] = create_eventcatalog(struct_of_tables, ...
    sensor_ids, wave_propertytable, sensors_with_event, ...
    timedeltas, modeled_timedeltas, corr_coefs, sensor_pairs);

save(catalog_outloc, 'eventcatalog', 'struct_of_tables')
writetable(eventtable, spreadsheet_outloc)

%% Plot the events which appear in >= 3 sensors from the catalog

plotevents_fromcatalog_simplified(eventcatalog, struct_of_tables, ...
    wt_struct, threshloc, outdir_plots, minsensors)

close all
