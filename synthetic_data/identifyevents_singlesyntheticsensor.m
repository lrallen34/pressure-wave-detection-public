%%% Generate a synthetic pressure time series with embedded wave events,
%%% given a certain noise floor, then try to detect the wave events
addpath(genpath('../'))

%% Define parameters
noise_floor = 0.8/100; % hPa
[wave_amplitudes, wave_periods] = meshgrid(logspace(-2, 0, 7), 1./(linspace(1/7200, 1/120, 15)));
wave_amplitudes = wave_amplitudes(:);
wave_periods = wave_periods(:);
wave_starts = datetime(2020,1,1,3,30,0):hours(4):datetime(2020,1,1,3,30,0)+hours(4*length(wave_periods)-4);
wave_ends = wave_starts + hours(2); % make each event 2 hours long
datetime_vec = datetime(2020,1,1):seconds(1):wave_ends(end)+hours(4);

real_events_table = table(wave_starts(:), wave_ends(:), ...
    wave_amplitudes(:), wave_periods(:), 'VariableNames', ...
    {'Start_time', 'End_time', 'Amplitude_hPa', 'Period_sec'});

%% Create the dataset
TT_synth = ...
    create_synthetic_dataset(datetime_vec, wave_starts, wave_ends, ...
    wave_periods, wave_amplitudes, noise_floor);

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
    filename = ['BMP388_TP_Log_synth1_' daily_str '.csv'];
    writetimetable(daily_TT, [csv_directory filename], 'WriteVariableNames', false);
end

%% Try to identify the events
datetime_span = [datetime_vec(1), datetime_vec(end)];
sensor_id = 'synth1';
smoothParam = 10;
threshloc = 'C:/Users/lrallen3/pressure_wavelet_means/all_means_asof_20220310.txt';

for k = 10
    % the mean wavelet power file is then read in
    threshmat = readmatrix(threshloc);
    period_vec = minutes(threshmat(:, 1));
    mean_vec = threshmat(:, 2);
    threshold_vec = (k/2)*mean_vec;
    max_events = 0;

    [events_table, ~, ~, ~, ~, ~, reconstructed_series, subsampled_datetimes, ...
        original_wavelet_transform] = identifyevents_iterative(...
        csv_directory, datetime_span, sensor_id, smoothParam, period_vec, ...
        threshold_vec, max_events);
    events_table.Duration = events_table.EndTime - events_table.StartTime;

%% Check which events were detected
    successful_detection = zeros(height(real_events_table), 1);
    false_positive = zeros(height(events_table), 1);
    for real_event = 1:height(real_events_table)
        if any(...
                events_table.CenterTime > real_events_table.Start_time(real_event) ...
                & events_table.CenterTime < real_events_table.End_time(real_event) ...
                & seconds(events_table.MinPeriod) < real_events_table.Period_sec(real_event) ...
                & seconds(events_table.MaxPeriod) > real_events_table.Period_sec(real_event))
            successful_detection(real_event) = 1;
        end
    end
    
    for detected_event = 1:height(events_table)
        if ~any(...
                events_table.CenterTime(detected_event) > real_events_table.Start_time ...
                & events_table.CenterTime(detected_event) < real_events_table.End_time ...
                & seconds(events_table.MinPeriod(detected_event)) < real_events_table.Period_sec ...
                & seconds(events_table.MaxPeriod(detected_event)) > real_events_table.Period_sec)
            false_positive(detected_event) = 1;
        end
    end
    disp(['K = ' num2str(k) ': detected ' ...
        num2str(sum(successful_detection)) '/' ...
        num2str(height(real_events_table)) ' with ' ...
        num2str(sum(false_positive)) ' false positives'])

end

f = figure;
scatter(real_events_table.Period_sec(successful_detection == 1)/60, ...
    real_events_table.Amplitude_hPa(successful_detection == 1), 24, ...
    "blue", "filled")
hold on
scatter(real_events_table.Period_sec(successful_detection == 0)/60, ...
    real_events_table.Amplitude_hPa(successful_detection == 0), 24, ...
    "red", "x")
set(gca, "XScale", "log")
set(gca, "YScale", "log")
grid on
title("Synthetic event detection")
xlabel("Wave period (min)")
ylabel("Amplitude (hPa)")
legend("Detected", "Undetected")

print('-r300', f, 'synthetic_event_detection', '-dpng')

%% Is the amplitude-duration relationship an artifact of the wavelet processing?
if isduration(events_table.Duration)
    events_table.Duration = minutes(events_table.Duration);
end

f = figure;
scatter(events_table.Duration, events_table.Amplitude, "filled")
grid on
ylabel('Amplitude (hPa)')
xlabel('Event Duration (min)')
xlim([0 170])
title('Amplitude-Duration for detected synthetic events')
print('-r300', f, 'synthetic_event_amplitudevsduration', '-dpng')

%% Event example
wt_norm = abs(original_wavelet_transform) ./ mean_vec;

evnum = 52;
ev_lims = [datetime(2020,1,8,11,0,0), datetime(2020,1,8,14,0,0)];

f = figure;
f.Position = [1 1 1400 700];
tiledlayout("vertical")
nexttile
plot(TT_synth.Time, TT_synth.Var2)
xlim(ev_lims)
grid on
ylabel('Pressure (hPa)')
title('Original pressure trace')
tick_times = xticks;
tick_labels = xticklabels;

nexttile
[~, minind] = min(abs(subsampled_datetimes - ev_lims(1)));
[~, maxind] = min(abs(subsampled_datetimes - ev_lims(2)));
tick_inds = zeros(length(tick_times), 1);
for curtick = 1:length(tick_times)
    [~, tick_inds(curtick)] = min(abs(subsampled_datetimes - tick_times(curtick)));
end

p2 = pcolor(1:length(subsampled_datetimes), minutes(period_vec), abs(original_wavelet_transform));
p2.EdgeColor = 'none';
hold on
title('Wavelet power')
c = colorbar;
c.Label.String = 'Wavelet Power (hPa^2 s^{-1})';
clim([0 0.8])
ylabel('Wave Period (min.)')
colormap(cmocean('rain'))
% for curline = 1:length(event_bounds)
%     curpoints = event_bounds{curline};
%     curx = curpoints(:, 2);
%     cury = periods(curpoints(:, 1));
%     plot(gca, curx, minutes(cury), 'm', 'LineWidth', 3)
% end
contour(1:length(subsampled_datetimes), minutes(period_vec), wt_norm, [5 5], '--k')
contour(1:length(subsampled_datetimes), minutes(period_vec), wt_norm, [10 10], 'k')
xlim([minind maxind])
set(gca, 'YScale', 'log')
% manually define x tick labels to make this work with older matlab
% versions
xticks(tick_inds)
xticklabels(tick_labels)

nexttile
plot(subsampled_datetimes, reconstructed_series{evnum})
xlim(ev_lims)
grid on
ylabel('Extracted event (hPa)')
title('Extracted events')

print(f, '-r300', 'example_series_short', '-dpng')

%% Figure for paper
f = figure;
corder = colororder('default');
colororder([0 0 0; corder(1,:)])
f.Position = [200 200 1000 400];
f.PaperSize = f.Position(3:4)./96;

t = tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding', 'compact');

ax = nexttile(t, [1 2]);
plot(subsampled_datetimes, reconstructed_series{evnum}, '-k')
xlim(ev_lims)
grid on
ylabel('Extracted event (hPa)')
title({'(a) Example synthetic pressure wave event', '(period: 2 min, amplitude: 0.0464 hPa)'})
tick_times = xticks;
tick_labels = xticklabels;
yyaxis right
plot(TT_synth.Time, TT_synth.Var2, 'Color', corder(1,:))
ylabel('Total pressure (hPa)')
legend('Extracted event', 'Total pressure', 'location', 'southeast')
ax.FontSize = 10;

ax = nexttile(t, 5, [1, 2]);
[~, minind] = min(abs(subsampled_datetimes - ev_lims(1)));
[~, maxind] = min(abs(subsampled_datetimes - ev_lims(2)));
tick_inds = zeros(length(tick_times), 1);
for curtick = 1:length(tick_times)
    [~, tick_inds(curtick)] = min(abs(subsampled_datetimes - tick_times(curtick)));
end

p2 = pcolor(1:length(subsampled_datetimes), minutes(period_vec), abs(original_wavelet_transform));
p2.EdgeColor = 'none';
hold on
title('(b) Wavelet power')
c = colorbar;
c.Label.String = 'Wavelet Power (hPa^2 s^{-1})';
clim([0 0.1])
ylabel('Wave Period (min.)')
colormap(cmocean('rain'))
% for curline = 1:length(event_bounds)
%     curpoints = event_bounds{curline};
%     curx = curpoints(:, 2);
%     cury = periods(curpoints(:, 1));
%     plot(gca, curx, minutes(cury), 'm', 'LineWidth', 3)
% end
contour(1:length(subsampled_datetimes), minutes(period_vec), wt_norm, [5 5], '--k')
contour(1:length(subsampled_datetimes), minutes(period_vec), wt_norm, [10 10], 'k')
legend('', '$W = 5\left< |W(b,a)| \right> _b$', ...
    '$W = 10\left< |W(b,a)| \right> _b$', 'interpreter', 'latex', ...
    'Location', 'northwest')
xlim([minind maxind])
set(gca, 'YScale', 'log')
% manually define x tick labels to make this work with older matlab
% versions
xticks(tick_inds)
xticklabels(tick_labels)
ax.FontSize = 10;

ax = nexttile(t, 3, [2 2]);
scatter(real_events_table.Period_sec(successful_detection == 1)/60, ...
    real_events_table.Amplitude_hPa(successful_detection == 1), 48, ...
    "blue", "filled")
hold on
scatter(real_events_table.Period_sec(successful_detection == 0)/60, ...
    real_events_table.Amplitude_hPa(successful_detection == 0), 48, ...
    "red", "x")
set(gca, "XScale", "log")
set(gca, "YScale", "log")
grid on
title("(c) Synthetic event detection with K = 10")
xlabel("Wave period (min)")
ylabel("Amplitude (hPa)")
legend("Detected", "Undetected")
ax.FontSize = 10;

print(f, '-vector', 'C:/Users/lrallen3/NCSU/writeups/pressurewaveletpaper/figures/f06', '-dpdf')
