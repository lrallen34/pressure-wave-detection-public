%%% Create a set of waves at different wave periods and phase velocities
%%% propagating through a 4-sensor network

% use sensors 004, 024, 026, and 034
sensor_lats = [43.634249, 43.7275, 43.699966, 43.781355];
sensor_lons = [-79.774896, -79.7866, -79.408338, -79.468298];
sensor_ids = {'EA_PiZero_004', 'eapizero024', 'eapizero026', 'eapizero034'};

wave_periods_sec = [120, 600, 3000, 7200]; % 2 min, 10 min, 50 min, 120 min
wave_amplitudes_hPa = [0.15, 0.45, 1.35, 3.24];
wave_speeds_ms1 = [4, 20, 60, 80];
wave_directions_deg = [5, 95, 185, 275];

wave_xspeeds = wave_speeds_ms1 .* sind(wave_directions_deg);
wave_yspeeds = wave_speeds_ms1 .* cosd(wave_directions_deg);

% use sensor 004 as an "anchor point" from which the time to reach the
% other sensors will be calculated
[dists, azimuths] = distance([sensor_lats(1), sensor_lons(1)], ...
    [sensor_lats(:), sensor_lons(:)], wgs84Ellipsoid("m"));
sensor_004_starts = datetime(2020,1,1,12,0,0):hours(24):datetime(2020,1,4,12,0,0);
sensor_004_ends = sensor_004_starts + hours(2);

xdists = dists .* sind(azimuths);
ydists = dists .* cosd(azimuths);
dist_vecs = [xdists ydists]';

timevec = zeros(4,4);
start_times = NaT(4,4);
end_time = NaT(4,4);
for event = 1:4
    timevec(:, event) = (1/wave_speeds_ms1(event)) * ...
            dists.* ...
            cosd(azimuths - ...
            wave_directions_deg(event));

    start_times(:, event) = sensor_004_starts(event) + round(seconds(timevec(:, event)));
    end_times(:, event) = sensor_004_ends(event) + round(seconds(timevec(:, event)));
end

noise_floor = 0.8/100; % hPa

%% Create the time series
datetime_vec = datetime(2020,1,1):seconds(1):datetime(2020,1,4,23,59,59);

for cursensor = 1:4
    TT_synth = ...
        create_synthetic_dataset(datetime_vec, start_times(cursensor, :), end_times(cursensor, :), ...
        wave_periods_sec, wave_amplitudes_hPa, noise_floor);

    % Output it to CSVs
    csv_directory = [pwd '/data/'];
    makeneededdirs(csv_directory)
    [h,m,s] = hms(TT_synth.Time);

    daily_starts = find(h == 0 & m == 0 & s == 0);
    if daily_starts(1) > 1
        daily_starts = [1; daily_starts];
    end

    daily_ends = find(h == 23 & m == 59 & s == 59);
    if daily_ends(end) < height(TT_synth)
        daily_ends = [daily_ends; height(TT_synth)];
    end

    for curday = 1:length(daily_starts)
        daily_TT = TT_synth(daily_starts(curday):daily_ends(curday), :);
        daily_TT.Time.Format = 'yyyyMMdd HHmmss.SSS';
        daily_str = char(daily_TT.Time(1), 'yyyyMMdd');
        filename = ['BMP388_TP_Log_' sensor_ids{cursensor} '_' daily_str '.csv'];
        writetimetable(daily_TT, [csv_directory filename], 'WriteVariableNames', false);
    end
end

%% Run the synthetic data through the processing algorithm
addpath(genpath(pwd))
datetime_span = [datetime(2020,1,1) datetime(2020,1,5)];
threshloc = 'C:/Users/lrallen3/pressure_wavelet_means/all_means_asof_20220310.txt';
coef = 5;
crosscorr_thresh = 0.65;
maxdelay = 2;
smoothParam = 10;

[struct_of_tables, wt_struct] = ...
    determinecoherence_fullnetwork_iterative('data/', datetime_span, ...
    sensor_ids, threshloc, coef, crosscorr_thresh, maxdelay, smoothParam);

struct_of_tables = assign_eventIDs(struct_of_tables);

[wave_propertytable,sensors_with_event,timedeltas,modeled_timedeltas,...
    corr_coefs,sensor_pairs] = ...
    calc_slownessvector(struct_of_tables, sensor_lats, sensor_lons);

[eventcatalog, eventtable] = create_eventcatalog(struct_of_tables, ...
    sensor_ids, wave_propertytable, sensors_with_event, ...
    timedeltas, modeled_timedeltas, corr_coefs, sensor_pairs);
