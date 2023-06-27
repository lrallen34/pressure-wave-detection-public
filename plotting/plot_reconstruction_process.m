function plot_reconstruction_process(dataloc, sensor_id, datetime_span, ...
    threshloc, imagedir, period_plotlims, varargin)
%PLOT_RECONSTRUCTION_PROCESS Creates figures illustrating the steps to
%identifying pressure wave events and outputs them to .png files.
%
%   INPUTS:
%
%       dataloc (str): location of pressure sensor .csv data files
%
%       sensor_id (str): ID of sensor
%
%       datetime_span (datetime): vector (length 2) of datetimes indicating
%       the start and end time of the span to be plotted
%
%       threshloc (str): location of .txt file with mean wavelet amplitudes
%       (produced using calcthresh.m)
%
%       imagedir (str): output image directory
%
%       period_plotlims (num): vector (length 2) of periods in minutes to
%       define y-limits of plots in time-period space
%
%       coef (optional, num): coefficient by which the mean wavelet
%       amplitudes are multiplied to define the event threshold (current
%       default is 5)
%
%   OUTPUTS:
%
%       none; images are output to subdirectories created in imagedir
%
%   Author: Luke Allen (lrallen3@ncsu.edu)
%

if numel(varargin) == 0
    coef = 5;
else 
    coef = varargin{1};
end

%omegadir = [imagedir '/omega/'];
omegaprimedir = [imagedir '/omega_prime/'];
wtdir = [imagedir '/wtamp/'];
normwtdir = [imagedir '/wt_norm/'];
recdir = [imagedir '/reconstructed/'];
makeneededdirs(imagedir, omegaprimedir, wtdir, normwtdir, recdir)

% read in the threshold table
threshmat = readmatrix(threshloc);
periods = minutes(threshmat(:, 1));
means = threshmat(:, 2);
thresholds = coef * means;

span_toread = [datetime_span(1)-hours(10), datetime_span(2)+hours(10)];
pressure_timetable = pressure_sensor_csv2timetable(dataloc, datetime_span, sensor_id);
smoothedtable = smoothdata(pressure_timetable, 'movmean', seconds(10), 'omitnan');
smoothedtable = smoothedtable(5:10:end, :);

[wt, datetime_vec, taumesh] = calc_wavelettransform(dataloc, span_toread, sensor_id, periods, 10);
[~, minind] = min(abs(datetime_vec - datetime_span(1)));
[~, maxind] = min(abs(datetime_vec - datetime_span(2)));
dtmesh = meshgrid(datetime_vec, taumesh(:, 1));

f1 = figure('visible', 'off');
ax = axes(f1, 'FontSize', 40);
f1.Position = [1 1 800 200];
p = pcolor(ax, dtmesh, minutes(taumesh), abs(wt));
p.EdgeColor = 'none';
title(['Wavelet power for ' sensor_id])
c = colorbar;
c.Label.String = 'Wavelet Power (hPa)';
ylabel('Wave Period (min.)')
colormap(cmocean('rain'))
clim([0 0.75])
ylim(period_plotlims)
xlim(datetime_span)
set(ax, 'YScale', 'log')
yticks(unique(round(logspace(log10(period_plotlims(1)), log10(period_plotlims(2)), 5))));
ax.FontSize = 14;
print('-r900', f1, [wtdir, sensor_id, char(datetime_span(1), 'yyyyMMdd_HHmmSS'), ...
    'to', char(datetime_span(2), 'yyyyMMdd_HHmmSS') '_wavelettransform'], '-dpng')

wt_norm = abs(wt) ./ means(1:size(wt,1));

[~, omega_prime, omega, event_centers, ~, event_bounds, ~, datetime_vec] = ...
    identifyevents_iterative(dataloc, span_toread, sensor_id, ...
    10, periods, thresholds, 0);

f3 = figure('visible', 'off');
ax3 = axes(f3, 'FontSize', 40);
f3.Position = [1 1 800 200];
p3 = pcolor(ax3, dtmesh, minutes(taumesh), double(omega_prime));
p3.EdgeColor = 'none';
colormap(cmocean('-greys'))
[rows,cols] = ind2sub(size(wt), event_centers);
hold on
scatter(cols, minutes(periods(rows)), 100, '.r')
for curline = 1:length(event_bounds)
    curpoints = event_bounds{curline};
    curx = curpoints(:, 2);
    cury = periods(curpoints(:, 1));
    plot(ax3, curx, minutes(cury), 'b', 'LineWidth', 3)
end
title('Points meeting threshold and reconstruction regions')
ylabel('Wave Period (min.)')
ylim(period_plotlims)
%xlim([minind maxind])
xlim(datetime_span)
set(ax3, 'YScale', 'log')
yticks(unique(round(logspace(log10(period_plotlims(1)), log10(period_plotlims(2)), 5))));
tick_times = xticks;
tick_labels = xticklabels;
legend('', 'Event Centers', 'Event Regions')
ax3.FontSize = 14;
print('-r900', f3, [omegaprimedir, sensor_id, ...
    char(datetime_span(1), 'yyyyMMdd_HHmmSS'), 'to', ...
    char(datetime_span(2), 'yyyyMMdd_HHmmSS') '_omegaprime'], '-dpng')

tick_inds = zeros(length(tick_times), 1);
for curtick = 1:length(tick_times)
    [~, tick_inds(curtick)] = min(abs(datetime_vec - tick_times(curtick)));
end

f2 = figure('visible', 'off');
ax2 = axes(f2);
f2.Position = [1 1 800 200];
%p2 = pcolor(ax2, dtmesh, minutes(taumesh), wt_norm);
p2 = pcolor(1:length(datetime_vec), minutes(periods), wt_norm);
p2.EdgeColor = 'none';
hold on
title(['Normalized wavelet power for ' sensor_id ' with event regions'])
c = colorbar;
c.Label.String = 'Norm. Wavelet Power';
clim([0 15])
ylabel('Wave Period (min.)')
colormap(cmocean('rain'))
for curline = 1:length(event_bounds)
    curpoints = event_bounds{curline};
    curx = curpoints(:, 2);
    cury = periods(curpoints(:, 1));
    plot(ax2, curx, minutes(cury), 'm', 'LineWidth', 3)
end
ylim(period_plotlims)
contour(ax2, 1:length(datetime_vec), minutes(periods), wt_norm, [5 5], '--k')
contour(ax2, 1:length(datetime_vec), minutes(periods), wt_norm, [10 10], 'k')
xlim([minind maxind])
set(ax2, 'YScale', 'log')
yticks(unique(round(logspace(log10(period_plotlims(1)), log10(period_plotlims(2)), 5))));
% manually define x tick labels to make this work with older matlab
% versions
xticks(tick_inds)
xticklabels(tick_labels)
ax2.FontSize = 14;
print('-r900', f2, [normwtdir, sensor_id, ...
    char(datetime_span(1), 'yyyyMMdd_HHmmSS'), 'to', ...
    char(datetime_span(2), 'yyyyMMdd_HHmmSS') '_normalizedwavelettransform'], '-dpng')

% reconstruct the events by inverting the wavelet transform over the
% reconstruction region (omega)
wt_reconstruct = wt;
wt_reconstruct(~logical(omega) | isnan(wt_reconstruct)) = 0;
pres_rec = icwt(wt_reconstruct, 'amor');

% plot the original time series with the reconstructed events
f5 = figure('visible', 'off');
f5.Position = [1 1 800 600];
ax5a = subplot(3,1,1);
plot(pressure_timetable.Datetime, pressure_timetable.Pressure, 'LineWidth', 2)
title(['Original pressure at ' sensor_id])
ylabel('Pressure (hPa)')
xlim(datetime_span)
grid on
ax5a.FontSize = 14;
ax5b = subplot(3,1,2);
plot(smoothedtable.Datetime, smoothedtable.Pressure, 'LineWidth', 2)
title(['10-second sampled pressure at ' sensor_id])
ylabel('Pressure (hPa)')
xlim(datetime_span)
grid on
ax5b.FontSize = 14;
ax5c = subplot(3,1,3);
plot(datetime_vec, pres_rec, 'k', 'LineWidth', 2)
title('Extracted Event')
ylabel('Pressure (hPa)')
xlim(datetime_span)
grid on
ax5c.FontSize = 14;
print('-r900', f5, [recdir, sensor_id, ...
    char(datetime_span(1), 'yyyyMMdd_HHmmSS'), 'to', ...
    char(datetime_span(2), 'yyyyMMdd_HHmmSS') '_eventreconstruction'], '-dpng')

close all

end

