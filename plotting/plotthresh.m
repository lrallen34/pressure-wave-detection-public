function f = plotthresh(datetoendon, datapath, imageloc)

%%% PLOTTHRESH plots the average wavelet power for all of the data through
%%% a given data for each of the sensor networks, to be used in conjunction
%%% with calcthresh.m script
%
% INPUTS:
%   - datetoendon (datetime): the date through which the mean wavelet power
%   is calculated
%   - datapath (string): the path to the directory containing text files
%   with mean wavelet power values (output using calcthresh.m)
%   - imageloc (string): path to the directory to which the resulting
%   figure will be output
%
% OUTPUT:
%   - f (figure): handle to the output figure
%
% Written by Luke Allen (lrallen3@ncsu.edu), last updated 2023/06/14

asofstring = char(datetoendon, 'yyyyMMdd');
dataloc = dir([datapath '/*' asofstring '*.txt']);
dataloc = dataloc(~startsWith({dataloc.name}, 'all'));

array_to_plot = struct();
f = figure;
ax = axes();
hold on
for i = 1:length(dataloc)
    curfile = [dataloc(i).folder '/' dataloc(i).name];
    filenamearray = split(dataloc(i).name, '_');
    array_to_plot(i).Name = filenamearray{1};
    array_to_plot(i).Data = readmatrix(curfile);
    
    plot(array_to_plot(i).Data(:, 1), array_to_plot(i).Data(:, 2), ...
        'LineWidth', 2)
end

% calculate overall mean
overall_sums = (array_to_plot(1).Data(:,2) .* array_to_plot(1).Data(:,3)) ...
    + (array_to_plot(2).Data(:,2) .* array_to_plot(2).Data(:,3)) ...
    + (array_to_plot(3).Data(:,2) .* array_to_plot(3).Data(:,3));
totalobs = array_to_plot(1).Data(:,3) + array_to_plot(2).Data(:,3) ...
    + array_to_plot(3).Data(:,3);
overall_means = overall_sums ./ totalobs;
plot(array_to_plot(1).Data(:, 1), overall_means, 'LineWidth', 2)
legendnames = {'Toronto', 'New York', 'Raleigh', 'Combined'};

xlabel('Wave Period (minutes)');
ylabel('hPa');
title('Mean Wavelet Power')
grid on
legend(legendnames, 'Location', 'northwest')
ax.FontSize = 18;
print('-r300', f, [imageloc '/meanwavelets_asof_' asofstring], '-dpng')
