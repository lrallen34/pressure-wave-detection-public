%%% Create a 4-panel plot showing two wave events from the Tonga volcanic
%%% shockwave
% Load in the event catalog
load('C:/Users/lrallen3/pressure_event_catalogs/all/tonga_shockwave_updated.mat')
figure_fontsize = 22;

% Choose the events to plot by hand and get the needed traces, plot
% endpoints, and sensor names
inbound_eventtrace = struct_of_tables(19).reconstructed_series{7};
inbound_fulltrace = struct_of_tables(19).fullTrace;
inbound_datetimes = struct_of_tables(19).datetime_vecs;
inbound_start = datetime(2022,1,15,15,15,0);
inbound_end = datetime(2022,1,15,16,15,0);

rebound_eventtrace = struct_of_tables(13).reconstructed_series{8};
rebound_fulltrace = struct_of_tables(13).fullTrace;
rebound_datetimes = struct_of_tables(13).datetime_vecs;
rebound_start = datetime(2022,1,16,6,45,0);
rebound_end = datetime(2022,1,16,7,45,0);

% Plot
f = figure('units','normalized','outerposition',[0 0 1 1]);
set(f, 'defaultAxesColorOrder', [0 0 0; 0 0.447 0.741])
ax = subplot(2,1,1);
plot(inbound_datetimes, inbound_eventtrace, 'LineWidth', 2)
grid on
xlim([inbound_start, inbound_end])
t=title('Inbound (ii) - EA_PiZero_002');
t.Interpreter = 'none';
ylabel('Extracted event (hPa)')
hold on
yyaxis right
plot(inbound_datetimes, inbound_fulltrace, 'LineWidth', 2)
ylabel('Total pressure (hPa)')
ax.FontSize = figure_fontsize;

% ax = subplot(2,2,2);
% plot(inbound_datetimes, inbound_fulltrace, 'LineWidth', 2)
% grid on
% xlim([inbound_start, inbound_end])
% t=title('Full pressure trace from EA_PiZero_002');
% t.Interpreter = 'none';
% ylabel('hPa')
% ax.FontSize = figure_fontsize;

ax = subplot(2,1,2);
plot(rebound_datetimes, rebound_eventtrace, 'LineWidth', 2)
grid on
xlim([rebound_start, rebound_end])
title('Rebound (ii) - eapizero018')
ylabel('Extracted event (hPa)')
hold on
yyaxis right
plot(rebound_datetimes, rebound_fulltrace, 'LineWidth', 2)
ylabel('Total pressure (hPa)')
ax.FontSize = figure_fontsize;

% ax = subplot(2,2,4);
% plot(rebound_datetimes, rebound_fulltrace, 'LineWidth', 2)
% grid on
% xlim([rebound_start, rebound_end])
% title('Full pressure trace from eapizero018')
% ylabel('hPa')
% ax.FontSize = figure_fontsize;

print('-r300', f, ['C:/Users/lrallen3/NCSU/writeups/pressurewaveletpaper/' ...
    'figures/tonga_examples_unedited'], '-dpng')
close(f)
