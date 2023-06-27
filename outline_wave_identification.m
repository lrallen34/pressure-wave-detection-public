%%% Script to outline the process of identifying a wave event from pressure
%%% time trace data

% Before doing anything, add all of the subdirectories containing needed
% functions to the Matlab path
addpath(genpath(pwd)) % alternatively, running init.m also does this

%% (1) Reading in the data
% Need to define the path to the data, the sensor ID we want, and the
% datetime span we want to read data for
data_path = '/path/to/data'; %FILLIN
sensor_id = 'eapizero014';
datetime1 = datetime(2022,2,4,17,30,0);
datetime2 = datetime(2022,2,4,22,0,0);
datetime_span_narrow = [datetime1, datetime2];
pressure_timetable = pressure_sensor_csv2timetable(...
data_path, datetime_span_narrow, sensor_id);

% Plot to see that there are waves in this part of the data
f1 = figure;
f1.Position = f1.Position .* [1 1 2 1]; % make the figure window larger
plot(pressure_timetable.Datetime, pressure_timetable.Pressure)
t = title(['Pressure at sensor ' sensor_id]);
t.Interpreter = 'none';
ylabel('Pressure (hPa)')

%% (2) Calculating the wavelet transform
% The wavelet transform is a method for capturing transient wave events,
% localized both in time and in frequency. To avoid edge effects when
% calculating, we look at a broader datetime span.
datetime3 = datetime1-hours(6);
datetime4 = datetime2+hours(6);
datetime_span_broad = [datetime3, datetime4];

% the function calc_wavelettransform reads in the data, smooths the data
% using a defined sampling period, then calculates the Morlet wavelet
% transform over a defined range of periods
sampling_period = 10; % seconds
period_interval = [minutes(1), minutes(120)]; % duration values
[wavelet_transform, datetime_vec, period_mesh, ~, smoothed_pressuretrace] = ...
calc_wavelettransform(data_path, datetime_span_broad, sensor_id, period_interval, sampling_period);

% the wavelet transform is a complex array, so for plotting we take the
% absolute value (or wavelet 'power')
wavelet_power = abs(wavelet_transform);
f2 = figure;
f2.Position = f2.Position .* [1 1 2 1.5]; % make the figure window larger
subplot(2,1,1) % plot the smoothed pressure trace
plot(datetime_vec, smoothed_pressuretrace)
t=title(['Pressure at sensor ' sensor_id]); t.Interpreter = 'none';
ylabel('Pressure (hPa)')

subplot(2,1,2) % plot the wavelet power
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), wavelet_power);
p.EdgeColor = 'none';
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Wavelet power for sensor ' sensor_id]); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Wavelet Power (hPa)';
colormap(cmocean('rain'))
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')

%% (3) Normalize the wavelet power using the time-averaged values
% Pressure tends to vary more at longer scales, so to identify relatively
% wavy portions of the pressure trace, it is useful to divide by the 
% average wavelet energy at each period from the full dataset
means_file = 'ny_means.txt';
means_array = readmatrix(means_file);
period_vec = means_array(:, 1);
means_vec = means_array(:, 2);

normalized_power = wavelet_power ./ means_vec;

f3 = figure; % plot the normalized wavelet power
f3.Position = f3.Position .* [1 1 2 1]; % make the figure window larger
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), normalized_power);
p.EdgeColor = 'none';
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Normalized wavelet power for sensor ' sensor_id]); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Normalized Wavelet Power';
colormap(cmocean('rain'))
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')

% Now the wavy parts of the data at shorter periods stick out more

%% (4) Identify wave events and extract them from the pressure trace
% Peaks in the normalized wavelet power can be taken to represent wave
% events. The function identifyevents_iterative reads in the data,
% calculates the normalized wavelet power, and uses the following steps to
% define wave event regions in the time-period space:

% 1. get the connected region to the maximum in the normalized wavelet
% transform that exceeds coef
%   1a. find the 'watersheds' of the negative normalized wavelet transform,
%   and remove any of the watersheds within the connected region found in 1
%   which lie such that the maximum in the normalized wavelet transform is
%   completely outside the period range of that watershed
% 2. extend that region along the time axis until a local minimum in the
% wavelet transform is reached; this is the region for the event
% 3. invert the wavelet transform in the event region and remove the event
% from the signal
% 4. recalculate the wavelet transform and normalized wavelet transform
% 5. return to (1), unless the maximum value in the normalized wavelet
% transform does not exceed 2 or the number of iterations has reached
% max_events

% The inverted wavelet transform over the event region obtained in step 3
% represents a reconstruction of the wave event, which we can then use in
% analysis with multiple sensors to approximate wave velocity.

% multiply the time-averaged wavelet power by some coefficient to define a
% threshold for wave events
threshold_vec = 5 * means_vec;
[events_table, ~, event_points, ~, ~, event_bounds, reconstructed_series, datetime_vec] = ...
    identifyevents_iterative(data_path, datetime_span_broad, sensor_id, ...
    sampling_period, minutes(period_vec), threshold_vec, 0);

f4 = figure;
f4.Position = f4.Position .* [1 1 3 1.5]; % make the figure window larger
subplot(3,1,1) % plot the original pressure trace
plot(pressure_timetable.Datetime, pressure_timetable.Pressure)
t = title(['Pressure at sensor ' sensor_id]);
t.Interpreter = 'none';
ylabel('Pressure (hPa)')

subplot(3,1,2) % plot the normalized wavelet power with event regions outlined
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), normalized_power);
[~, plotStart] = min(abs(datetime_vec - datetime1));
[~, plotEnd] = min(abs(datetime_vec - datetime2));
xlim([plotStart plotEnd])
p.EdgeColor = 'none';
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Normalized wavelet power for sensor ' sensor_id ' with event regions']); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Normalized Wavelet Power';
colormap(cmocean('rain'))
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')
hold on
for cur_event = 1:length(event_bounds)
    cur_bounds = event_bounds{cur_event};
    curx = cur_bounds(:, 2);
    cury = period_vec(cur_bounds(:, 1));
    plot(curx, cury, 'm', 'LineWidth', 3)
end
hold off

subplot(3,1,3) % plot the extracted event trace for the 1st event in the 
% table (the one with the largest peak in normalized wavelet power) which
% is within the original datetime span
events_table_trimmed = events_table(events_table.CenterTime < datetime_span_narrow(2) ...
    & events_table.CenterTime > datetime_span_narrow(1), :);
reconstructed_series_trimmed = reconstructed_series(events_table.CenterTime < datetime_span_narrow(2) ...
    & events_table.CenterTime > datetime_span_narrow(1));
plot(datetime_vec, reconstructed_series_trimmed{1})
xlim([datetime1 datetime2])
t = title(['Extracted wave event centered at ' char(events_table_trimmed.CenterTime(1))]);
t.Interpreter = 'none';
ylabel('Pressure (hPa)')

%% Additional plots
% three-panel figure to show wavelet power, normalized wavelet power, and
% event regions

fe1 = figure;
fe1.Position = fe1.Position .* [1 1 3 1.5]; % make the figure window larger

[~, plotStart] = min(abs(datetime_vec - datetime1));
[~, plotEnd] = min(abs(datetime_vec - datetime2));

subplot(3,1,1) % plot the wavelet power
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), wavelet_power);
p.EdgeColor = 'none';
xlim([plotStart plotEnd])
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Wavelet power for sensor ' sensor_id]); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Wavelet Power (hPa)';
colormap(gca, cmocean('rain'))
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')

subplot(3,1,2) % plot the normalized wavelet power
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), normalized_power);
p.EdgeColor = 'none';
xlim([plotStart plotEnd])
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Normalized wavelet power for sensor ' sensor_id]); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Normalized Wavelet Power';
colormap(gca, cmocean('rain'))
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')

subplot(3,1,3) % plot the points exceeding the threshold and event region outlines
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), double(normalized_power > 5));
p.EdgeColor = 'none';
xlim([plotStart plotEnd])
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Regions where normalized power > 5 and event region outlines for sensor ' sensor_id]);
t.Interpreter = 'none';
c = colorbar('Ticks', [0.25 0.75], 'TickLabels', {'<= 5','> 5'});
c.Label.String = 'Normalized Wavelet Power';
colormap(gca, [0 0 0; 1 1 1])
set(gca,'YScale','log')
hold on
for cur_event = 1:length(event_bounds)
    cur_bounds = event_bounds{cur_event};
    curx = cur_bounds(:, 2);
    cury = period_vec(cur_bounds(:, 1));
    plot(curx, cury, 'm', 'LineWidth', 3)
end
hold off

% figure to illustrate imhmin transform
[wt_max, ind_max] = max(normalized_power, [], 'all', 'linear');
[row_max, col_max] = ind2sub(size(normalized_power), ind_max);
[omega_curevent, idx] = bwselect(normalized_power > 5, col_max, row_max);

normalized_power = wavelet_power ./ means_vec;
normalized_power(isnan(normalized_power)) = 0;
normalized_power_imhmin = imhmin(normalized_power, 7.5);
normalized_power(~omega_curevent) = 0;
normalized_power_imhmin(~omega_curevent) = 0;
fe2 = figure;
fe2.Position = fe2.Position .* [1 1 2 1.5]; % make the figure window larger

subplot(2,1,1)
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), normalized_power);
[~, plotStart] = min(abs(datetime_vec - datetime1));
[~, plotEnd] = min(abs(datetime_vec - datetime2));
xlim([plotStart plotEnd])
p.EdgeColor = 'none';
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Normalized wavelet power for sensor ' sensor_id]); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Normalized Wavelet Power';
colormap(cmocean('rain'))
clim([0 20])
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')

subplot(2,1,2)
p = pcolor(1:length(datetime_vec), minutes(period_mesh(:, 1)), normalized_power_imhmin);
[~, plotStart] = min(abs(datetime_vec - datetime1));
[~, plotEnd] = min(abs(datetime_vec - datetime2));
xlim([plotStart plotEnd])
p.EdgeColor = 'none';
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Normalized wavelet power for sensor ' sensor_id ' after imhmin']); t.Interpreter = 'none';
c = colorbar;
c.Label.String = 'Normalized Wavelet Power';
colormap(cmocean('rain'))
clim([0 20])
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')

% figure to show watersheds before and after imhmin transform
negative_wtnorm = -normalized_power;
negative_wtnorm_imhmax = imhmax(negative_wtnorm, 7.5);
wtnorm_watersheds_imhmax = watershed(negative_wtnorm_imhmax);
wtnorm_watersheds = watershed(negative_wtnorm);
wtnorm_watersheds_imhmax(~omega_curevent) = 0;
wtnorm_watersheds(~omega_curevent) = 0;

rgb_original = label2rgb(wtnorm_watersheds);
rgb_new = label2rgb(wtnorm_watersheds_imhmax);

% include lines at top and bottom of central watershed
central_watershed = wtnorm_watersheds(ind_max);
[centralwatershed_rows,~] = find(wtnorm_watersheds == central_watershed);
% include a 1 pixel buffer
centralwatershed_minperiod = period_vec(min(centralwatershed_rows)-1);
centralwatershed_maxperiod = period_vec(max(centralwatershed_rows)+1);

central_watershed_imhmax = wtnorm_watersheds_imhmax(ind_max);
[centralwatershed_rows_imhmax,~] = find(wtnorm_watersheds_imhmax == central_watershed_imhmax);
% include a 1 pixel buffer
centralwatershed_minperiod_imhmax = period_vec(min(centralwatershed_rows_imhmax)-1);
centralwatershed_maxperiod_imhmax = period_vec(max(centralwatershed_rows_imhmax)+1);

fe3 = figure;
fe3.Position = fe3.Position .* [1 1 2 1.5]; % make the figure window larger
zvals = zeros(size(normalized_power));

subplot(2,1,1)
surface(1:length(datetime_vec), minutes(period_mesh(:, 1)), zvals, rgb_original, 'EdgeColor', 'none');
[~, plotStart] = min(abs(datetime_vec - datetime1));
[~, plotEnd] = min(abs(datetime_vec - datetime2));
xlim([plotStart plotEnd])
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Watersheds within region exceeding event threshold for ' sensor_id]); t.Interpreter = 'none';
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')
hold on
plot(1:length(datetime_vec), repmat(centralwatershed_minperiod, 1, length(datetime_vec)), '--r')
plot(1:length(datetime_vec), repmat(centralwatershed_maxperiod, 1, length(datetime_vec)), '--r')
hold off

subplot(2,1,2)
surface(1:length(datetime_vec), minutes(period_mesh(:, 1)), zvals, rgb_new, 'EdgeColor', 'none');
[~, plotStart] = min(abs(datetime_vec - datetime1));
[~, plotEnd] = min(abs(datetime_vec - datetime2));
xlim([plotStart plotEnd])
% manually define x tick labels to make this work with older matlab
% versions
xt = xticks;
labelstrings = char(datetime_vec(xt), 'HH:MM:SS');
xticklabels(labelstrings)
ylabel('Wave Period (min.)')
t=title(['Watersheds within region exceeding event threshold for ' sensor_id ' after imhmin']); t.Interpreter = 'none';
% put the period axis on a log scale (equivalent to making the axis linear
% wrt frequency)
set(gca,'YScale','log')
hold on
plot(1:length(datetime_vec), repmat(centralwatershed_minperiod_imhmax, 1, length(datetime_vec)), '--r')
plot(1:length(datetime_vec), repmat(centralwatershed_maxperiod_imhmax, 1, length(datetime_vec)), '--r')
hold off
