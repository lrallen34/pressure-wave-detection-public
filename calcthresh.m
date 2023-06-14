%%% Script for calculating and outputting mean wavelet transform by period
%%% across entire pressure data set for a given network
addpath(genpath(pwd)) % add subdirectories including dependencies to path
% should only need to change network and outloc
networks = {'ny', 'canada', 'raleigh'};
% define output location with the date on which this script is run
%datetoendon = datetime('yesterday');
datetoendon = datetime(2022,3,10);
stringtoendon = char(datetoendon, 'yyyymmdd');
smoothParam = 10; % duration in seconds to average pressure over
outdir = '/home/disk/zathras/lukea41/pressure_wavelet_means/';
imageloc = '/waveletpowermean_plots/';
makeneededdirs(outdir, imageloc)

for curnet = 1:length(networks)
    network = networks{curnet};
    
    % define location of .txt file with sensor IDs
    sensor_txt = ['/home/disk/zathras/lukea41/sensor_lists/' network '_sites.txt'];
    outloc = [outdir network '_means_asof_' stringtoendon '_smooth' ...
        num2str(smoothParam) '.txt'];
    % read in sensor IDs
    sensorMat = readmatrix(sensor_txt, 'OutputType', 'char');
    sensors = sensorMat(:,1);
    
    % location of pressure sensor CSV files
    datadir = '/home/disk/ivanova2/RPi_Pressure_Data/';
    
    % loop through each sensor; preallocate first
    Asensor = zeros(70, length(sensors));
    sensorobs = zeros(70, length(sensors));
    for cursensor = 1:length(sensors)
        % get datetime_spans corresponding to stretches of consecutive days
        % with data for the sensor
        streak_spans = identifydatastreaks(datadir, sensors{cursensor});
        
        % preallocate for this sensor
        totalobs = zeros(70, 1);
        totalA = zeros(70, 1);
        
        for cur_span = 1:size(streak_spans, 1)
            
            datetime_span = streak_spans(cur_span, :);
            if datetime_span(1) > datetoendon
                break
            end
            if datetime_span(2) > datetoendon
                datetime_span(2) = datetoendon + days(1) - seconds(1);
            end
            disp(['Processing streak for ' sensors{cursensor} ' from ' ...
                char(datetime_span(1)) ' to ' char(datetime_span(2))])
            
            % apply the wavelet transform to the streak
            [wt, ~, taumesh, coimesh] = calc_wavelettransform(datadir, datetime_span, ...
                sensors{cursensor}, [minutes(1) minutes(120)], 5);
            % the wavelet transform function will output empty arrays if there
            % is not enough data in the given timespan, so skip if this is the
            % case
            if isempty(wt)
                continue
            end
            
            % find the time average of the wavelet magnitude for each period
            A = mean(abs(wt), 2, 'omitnan');
            
            % get the number of valid observations at each period
            streakobs = sum(~isnan(wt), 2);
            totalobs(1:length(A)) = totalobs(1:length(A)) + streakobs;
            
            % multiply the time average by the number of valid observations
            Amult = streakobs .* A;
            % sum that multiplied value
            totalA(1:length(A)) = totalA(1:length(A)) + Amult;
        end
        % after going through all the data for the current sensor, put the
        % final summed values into arrays
        Asensor(1:length(A), cursensor) = totalA;
        sensorobs(1:length(A), cursensor) = totalobs;
    end
    
    % from each sensor's resulting data, get the overall mean at each period
    Aovrtotal = sum(Asensor, 2);
    obsovrtotal = sum(sensorobs, 2);
    % divide back by the number of observations to get the mean wavelet
    % magnitudes
    Aoverall = Aovrtotal ./ obsovrtotal;
    tau = taumesh(:, 1);
    
    % put into a table and output
    outtable = table(minutes(tau), Aoverall, obsovrtotal, ...
        'VariableNames', {'Period', 'MeanWaveletAmp', 'TotalObs'});
    writetable(outtable, outloc);
    
end

%% Output a plot of the mean wavelet power values
f = plotthresh(datetoendon, outdir, imageloc);
