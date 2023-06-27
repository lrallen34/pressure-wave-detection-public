%%% Script for plotting all pressure traces for a given network together
network = 'ny';
datetime_span = [datetime(2022,2,4,0,0,0), datetime(2022,2,5,0,0,0)];
datadir = '/path/to/data'; % FILLIN
outdir = ['/path/to/pressuretraces/' network]; % FILLIN: directory to output .png files to

sensorloc = [network '_sites.txt'];
sensormat = readmatrix(sensorloc, 'OutputType', 'char');
sensor_names = sensormat(:, 1);
legend_entries = sensor_names;
sensor_lats = str2double(sensormat(:, 2));
sensor_lons = str2double(sensormat(:, 3));

sensor_colormat = readtable('sensor_rgbcodes.txt');
sensor_colormat.Properties.VariableNames = {'sensor', 'red', 'green', 'blue'};


f = figure;
hold on
for cur_sensor = sensor_names'
    % in case no data are available, which will produce an error, use a
    % try-catch call on the reader function
    sensor_color = table2array(sensor_colormat(strcmp(sensor_colormat.sensor, cur_sensor{1}), 2:4))/255;
    try
        sensor_timetable = pressure_sensor_csv2timetable(datadir, datetime_span, cur_sensor{1});
    catch
        sensor_names(strcmp(sensor_names, cur_sensor{1})) = [];
        continue
    end
    plot(sensor_timetable.Datetime, sensor_timetable.Pressure, 'LineWidth', 2, 'color', sensor_color)
end
ylabel('Pressure (hPa)')
f.Position = [1 1 1920 1080];
title(['Pressure traces for ' network])
l=legend(sensor_names, 'Location', 'best');
l.Interpreter = 'none';
grid on
makeneededdirs(outdir)
print('-r300', f, [outdir ...
    '/' char(datetime_span(1), 'yyyymmddHHMM') '_' ...
    char(datetime_span(2), 'yyyymmddHHMM')], '-dpng')
