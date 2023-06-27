%%% Create a sequence of plots to illustrate the process of matching events
%%% in one sensor to events in another sensor
% define sensor IDs, timespan to look at, and locations of information on
% the event thresholds
addpath(genpath('../'))
sensor_id = 'eapizero025'; % 'primary' sensor
%datadir = '/home/disk/ivanova2/RPi_Pressure_Data/';
datadir = 'C:/Users/lrallen3/rpi_pressure_data/';
network = 'canada';
datetime_span = [datetime(2023,2,22,12,0,0), datetime(2023,2,23,12,0,0)];
plotting_span = [datetime(2023,2,23,1,30,0), datetime(2023,2,23,4,0,0)];
%threshloc = 'C:/Users/lrallen3/pressure_wavelet_means/all_means_asof_20220310.txt';
threshloc = ['C:/Users/lrallen3/pressure_wavelet_means/' network '_means_asof_20220310.txt'];
threshmat = readmatrix(threshloc);
period_vec = minutes(threshmat(:, 1));
mean_vec = threshmat(:, 2);
thresholds = 5.*mean_vec;
merge_regions = 0;

alt_sensors = {'EA_PiZero_004', 'eapizero023', 'eapizero024', 'eapizero034'}; % 'secondary' sensor
n_sensors = length(alt_sensors) + 1;

% calculate wavelet transform in primary sensor
[wavelet_transform, datetime_vec, ~, ~, primary_pressuretrace] = ...
    calc_wavelettransform(datadir, datetime_span, sensor_id, period_vec, 10);
wavelet_transform(isnan(wavelet_transform)) = 0;
wt_norm = abs(wavelet_transform) ./ repmat(mean_vec, 1, size(wavelet_transform, 2));
% calculate wavelet transform and identify events in secondary sensors
for cur_altsensor = 1:length(alt_sensors)
    [wt_alt{cur_altsensor}, alt_datetime_vec{cur_altsensor}, ~, ~, ...
        secondary_pressuretrace{cur_altsensor}] = ...
        calc_wavelettransform(datadir, datetime_span, alt_sensors{cur_altsensor}, period_vec, 10);
    wt_alt{cur_altsensor}(isnan(wt_alt{cur_altsensor})) = 0;
    wt_alt_norm{cur_altsensor} = ...
        abs(wt_alt{cur_altsensor}) ./ repmat(mean_vec, 1, size(wt_alt{cur_altsensor}, 2));
    [alt_events_table{cur_altsensor}, alt_omegaprime{cur_altsensor}, ...
        alt_omega{cur_altsensor}, altevent_centers{cur_altsensor}, ~, ~, ...
        secondary_extracted_series{cur_altsensor}, alt_dtvec{cur_altsensor}] = ...
        identifyevents_iterative(datadir, datetime_span, alt_sensors{cur_altsensor}, 10, ...
        period_vec, thresholds, 0);
end

% identify events in primary sensor
[events_table, omega_prime, omega, event_centers, omega_pixels, ...
    event_bounds, primary_extracted_series, prim_dtvec] = ...
    identifyevents_iterative(datadir, datetime_span, sensor_id, 10, ...
    period_vec, thresholds, 0);

%% First plot: original smoothed pressure traces from both sensors

f1 = figure;
f1.Position = [0 0 960 1080];
subplot(n_sensors,1,1)
plot(datetime_vec, primary_pressuretrace)
xlim(plotting_span)
ylabel('hPa')
grid on
title(['Full pressure trace for ' sensor_id])

for cur_altsensor = 1:length(alt_sensors)
    subplot(n_sensors,1,cur_altsensor+1)
    plot(alt_datetime_vec{cur_altsensor}, secondary_pressuretrace{cur_altsensor})
    xlim(plotting_span)
    ylabel('hPa')
    grid on
    title(['Full pressure trace for ' alt_sensors{cur_altsensor}], ...
        'Interpreter', 'none')
end

%% Second plot: reconstructed events for both sensors, unshifted
% invert wavelet transforms over the event regions
primary_event_trace = primary_extracted_series{1};

f2 = figure;
f2.Position = [0 0 960 1080];
curax = subplot(n_sensors,1,1);
plot(prim_dtvec, primary_event_trace, 'k', 'LineWidth', 2)
xlim(plotting_span)
ylabel('hPa')
grid on
title(['Extracted event in ' sensor_id])
curax.FontSize = 14;

for cur_altsensor = 1:length(alt_sensors)
    secondary_event_trace = secondary_extracted_series{cur_altsensor}{1};
    curax = subplot(n_sensors,1,cur_altsensor+1);
    plot(alt_dtvec{cur_altsensor}, secondary_event_trace, 'k', 'LineWidth', 2)
    xlim(plotting_span)
    ylabel('hPa')
    grid on

    % calculate unshifted correlation
    corrmat_unshiftedtrace = corrcoef(primary_event_trace, secondary_event_trace);
    corr_unshiftedtrace = corrmat_unshiftedtrace(2);
    title(['Extracted event in ' alt_sensors{cur_altsensor} ...
        '; Shift = 0 s, R = ' num2str(corr_unshiftedtrace)], 'Interpreter', ...
        'none')
    curax.FontSize = 14;
end
print(f2, '-r300', ['C:/Users/lrallen3/NCSU/writeups/pressurewaveletpaper/' ...
    'figures/reconstructedtraces_unshifted'], '-dpng')

%% Calculate optimal shift to maximize lagged correlation between the 
% reconstructed event traces
evstart = events_table.StartTime(1);
evend = events_table.EndTime(1);
for cur_altsensor = 1:length(alt_sensors) % for each secondary sensor
    secondary_event_trace = secondary_extracted_series{cur_altsensor}{1};
    alt_evstart = alt_events_table{cur_altsensor}.StartTime(1);
    alt_evend = alt_events_table{cur_altsensor}.EndTime(1);
    dtmin_rec = alt_datetime_vec{cur_altsensor}(1)-alt_evstart;
    dtmax_rec = alt_datetime_vec{cur_altsensor}(end)-alt_evend;
    deltatrange_rec = [dtmin_rec, dtmax_rec];
    earlieststart = min([evstart, alt_evstart]);
    latestend = max([evend, alt_evend]);
    % trim arrays to the duration of the event to simplify
    % cross-correlation calculation
    alt_rec_forcorr = secondary_event_trace( ...
        alt_datetime_vec{cur_altsensor} >= earlieststart & ...
        alt_datetime_vec{cur_altsensor} <= latestend);
    cur_rec_forcorr = primary_event_trace( ...
        datetime_vec >= earlieststart & ...
        datetime_vec <= latestend);

    % pad the arrays if necessary to make them the same length
    if max(datetime_vec) < latestend
        cur_rec_forcorr = [cur_rec_forcorr, zeros(1, seconds(latestend - max(datetime_vec))/10)];
    end
    if min(datetime_vec) > earlieststart
        cur_rec_forcorr = [zeros(1, seconds(min(datetime_vec) - earlieststart)/10), cur_rec_forcorr];
    end
    if max(alt_datetime_vec{cur_altsensor}) < latestend
        alt_rec_forcorr = [alt_rec_forcorr, zeros(1, seconds(latestend - max(alt_datetime_vec{cur_altsensor}))/10)];
    end
    if min(alt_datetime_vec{cur_altsensor}) > earlieststart
        alt_rec_forcorr = [zeros(1, seconds(min(alt_datetime_vec{cur_altsensor}) - earlieststart)/10), alt_rec_forcorr];
    end
    % calculate cross-correlations then get the maximum value and
    % corresponding delay
    [crosscorr_rec, lags] = xcorr(alt_rec_forcorr, cur_rec_forcorr, 'normalized');
    optimal_crosscorr_final{cur_altsensor} = max(crosscorr_rec);
    % in the event of a 'tie' regarding optimal delay time for
    % maximizing the lagged correlation, take an average
    delaytime{cur_altsensor} = mean(10*lags(crosscorr_rec == max(crosscorr_rec)));
end

%% Third plot: reconstructed events for both sensors, optimally shifted in time to maximize correlation
f3 = figure;
f3.Position = [0 0 1000 900];
curax = subplot(n_sensors,1,1);
plot(datetime_vec, primary_event_trace, 'k', 'LineWidth', 2)
xlim(plotting_span)
ylabel('hPa')
grid on
title(['Extracted event in ' sensor_id])
curax.FontSize = 12;

for cur_altsensor = 1:length(alt_sensors)
    curax = subplot(n_sensors,1,cur_altsensor+1);
    plot(alt_datetime_vec{cur_altsensor}, ...
        secondary_extracted_series{cur_altsensor}{1}, ...
        'LineWidth', 2, 'Color', [166 118 29]/255)
    hold on
    plot(alt_datetime_vec{cur_altsensor} - seconds(delaytime{cur_altsensor}), ...
        secondary_extracted_series{cur_altsensor}{1}, 'k', 'LineWidth', 2)
    xlim(plotting_span)
    ylabel('hPa')
    grid on
    title(['Extracted event in ' alt_sensors{cur_altsensor} ...
        '; Shift = ' char(duration(0, 0, -delaytime{cur_altsensor})) ...
        ', R = ' num2str(optimal_crosscorr_final{cur_altsensor})], ...
        'Interpreter', 'none')
    curax.FontSize = 12;
    legend('Original', 'Shifted')
end

print(f3, '-r300', 'C:/Users/lrallen3/NCSU/writeups/pressurewaveletpaper/figures/reconstructedtraces_shifted', '-dpng')
