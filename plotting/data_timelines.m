%%% Script to generate timeline plot showing when each sensor had data available

% Should only need to change the network variable ('raleigh', 'ny', or
% 'canada')
network = 'raleigh';
% read in the list of sites for that network and get the IDs that have
% defined coordinates
sensorloc = ['../' network '_sites.txt'];
datadir = '/path/to/data/'; % FILLIN
sensormat = readmatrix(sensorloc, 'OutputType', 'char');
sensorlats = str2double(sensormat(:, 2));
all_ids = sensormat(~isnan(sensorlats), 1);

% plot
f = figure;
f.Position = f.Position .* [1 1 1 length(all_ids)/10];
hold on
for cur_sensor = 1:length(all_ids)
    sensor_id = all_ids{cur_sensor};
    datetime_spans = identifydatastreaks(datadir, sensor_id);
    ymatrix = repmat(cur_sensor, size(datetime_spans));
    
    plot(datetime_spans', ymatrix', 'Color', [0 0.4470 0.7410], 'LineWidth', 3)
end

% add in labels for each sensor
yticks(1:length(all_ids))
yticklabels(all_ids);
xlim([datetime(2020,1,1,0,0,0), datetime(2023,1,1,0,0,0)])
set(gca,'TickLabelInterpreter','none')
grid on
switch network
    case 'canada'
        desc = ' Toronto ';
    case 'ny'
        desc = ' NY/LI ';
    case 'raleigh'
        desc = ' Raleigh ';
    case 'all'
        desc = ' ';
    otherwise
        error('unsupported network string')
end
title(['Periods when each' desc 'sensor was active'])

% output to an image file and close the figure
print('-r300', f, ['sensor_data_timelines_' network], '-dpng')
close(f)
