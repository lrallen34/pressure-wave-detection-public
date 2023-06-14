function [event_table,sensors_with_event,timedeltas,modeled_timedeltas,...
    corr_coefs,sensor_pairs] = ...
    calc_slownessvector(struct_of_tables, latvec, lonvec)
%CALC_SLOWNESSVECTOR Calculates slowness vectors (essentially the inverse
%of velocity vectors) for each wave event in a struct_of_tables with wave
%event IDs (i.e., one which has been processed through assign_eventIDs.m).
%The function then outputs a table summarizing wave events and additional
%optional outputs with diagnostics on the calculation.
%
%   General form: [event_table,sensors_with_event,timedeltas,modeled_timedeltas,...
%   corr_coefs,sensor_pairs] =  calc_slownessvector(struct_of_tables, latvec, lonvec)
%
%   INPUTS:
%       -struct_of_tables (struct): structure array with tables of
%       information on the wave events in each pressure sensor within a
%       network produced by determinecoherence_fullnetwork_iterative.m.
%       Must include event IDs as assigned by assign_eventIDs.m
%       -latvec (numeric vector): vector of latitude values for the
%       locations of each pressure sensor. The order of the values must
%       align with the order of the pressure sensors in struct_of_tables
%       -lonvec (numeric vector): vector of longitude values for the
%       locations of each pressure sensor. The order of the values must
%       align with the order of the pressure sensors in struct_of_tables
%
%   OUTPUTS:
%       -event_table (table): table summarizing the wave events, represented by
%       unique event IDs from the struct_of_tables. Each wave event is on a
%       single row, with columns for the number of sensors, mean amplitude
%       in hPa, start and end times of the event, the central period of the
%       event, the wave velocity components, and various measures of
%       uncertainty/error in the slowness vector calculation
%       -sensors_with_event (cell): cell array in which each cell
%       represents an event ID. The elements within each cell are indices
%       corresponding to the sensors which captured the event (aligning
%       with the sensor indices in struct_of_tables)
%       -timedeltas (cell): cell array in which each cell represents an
%       event ID. The elements within each cell are 'observed' delay times
%       between the event passage in pairs of sensors, corresponding to
%       sensor_pairs.
%       -modeled_timedeltas (cell): cell array in which each cell represents an
%       event ID. The elements within each cell are the delay times
%       between the event passage in pairs of sensors which result from the
%       calculated slowness vector, corresponding to sensor_pairs.
%       -corr_coefs (cell): cell array in which each cell represents an
%       event ID. The elements within each cell are the maximum of the
%       cross-correlation function for the extracted wave traces in pairs
%       of sensors, corresponding to sensor_pairs.
%       -sensor_pairs (cell): cell array in which each cell represents an
%       event ID. The elements within each cell are each possible pairing
%       of sensors which captured a given wave event, given as a string.
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Dec 2022
%
% Reference: Del Pezzo and Giudicepietro (2002), Computers and Geosciences

% get an array of the unique event IDs in the struct_of_tables
allIDs = {struct_of_tables.eventIDs};
concatenatedIDs = {cat(1, allIDs{:})};
uniqueIDs = unique(concatenatedIDs{1});
% preallocating output arrays
slownessvec = cell(length(uniqueIDs), 1);
wavespeed = NaN(length(uniqueIDs), 1);
wavedirection = NaN(length(uniqueIDs), 1);
waveu = NaN(length(uniqueIDs), 1);
wavev = NaN(length(uniqueIDs), 1);
corr_coefs = cell(length(uniqueIDs), 1);
mean_corrs = NaN(length(uniqueIDs), 1);
sensor_pairs = cell(length(uniqueIDs), 1);
timedeltas = cell(length(uniqueIDs), 1);
N_sensors = zeros(length(uniqueIDs), 1);
modeled_timedeltas = cell(length(uniqueIDs), 1);
timedelta_RMSE = NaN(length(uniqueIDs), 1);
timedelta_NRMSE = NaN(length(uniqueIDs), 1);
error_x = NaN(length(uniqueIDs), 1);
error_y = NaN(length(uniqueIDs), 1);
sensors_with_event = cell(length(uniqueIDs), 1);
mean_amp = NaN(length(uniqueIDs), 1);
event_start = NaT(length(uniqueIDs), 1);
event_end = NaT(length(uniqueIDs), 1);
event_centerperiod = minutes(NaN(length(uniqueIDs), 1));
event_maxperiod = minutes(NaN(length(uniqueIDs), 1));
event_minperiod = minutes(NaN(length(uniqueIDs), 1));
event_minperiod_peakregion = minutes(NaN(length(uniqueIDs), 1));
event_maxperiod_peakregion = minutes(NaN(length(uniqueIDs), 1));

% loop through each event ID
for curevent = 1:length(uniqueIDs)
    curevid = uniqueIDs(curevent);
    sensors_with_event{curevent} = [];
    timedeltas{curevent} = [];
    distancematrix = [];
    sensor_pairs{curevent} = string([]);
    amplitudes_curevent = [];
    starts_curevent = [];
    ends_curevent = [];
    periods_curevent = [];
    periodmins_curevent = [];
    periodmaxs_curevent = [];
    periodmins_peakregion_curevent = [];
    periodmaxs_peakregion_curevent = [];
    for cursensor = 1:length(struct_of_tables)
        if any(struct_of_tables(cursensor).eventIDs == curevid)
            eventind_cursensor = find(struct_of_tables(cursensor).eventIDs == curevid);
            % for each previous sensor that had the event, get the event
            % traces to calculate delay times and get the distance vectors
            cursensor_trace = struct_of_tables(cursensor).reconstructed_series{eventind_cursensor};
            cursensor_datetimes = struct_of_tables(cursensor).datetime_vecs;
            cursensor_start = struct_of_tables(cursensor).events_tables.StartTime(eventind_cursensor);
            cursensor_end = struct_of_tables(cursensor).events_tables.EndTime(eventind_cursensor);
            
            amplitudes_curevent = [amplitudes_curevent; max(cursensor_trace) - min(cursensor_trace)];
            starts_curevent = [starts_curevent; cursensor_start];
            ends_curevent = [ends_curevent; cursensor_end];
            periods_curevent = [periods_curevent; ...
                struct_of_tables(cursensor).events_tables.CenterPeriod(eventind_cursensor)];
            periodmins_curevent = [periodmins_curevent; ...
                struct_of_tables(cursensor).events_tables.MinPeriod(eventind_cursensor)];
            periodmaxs_curevent = [periodmaxs_curevent; ...
                struct_of_tables(cursensor).events_tables.MaxPeriod(eventind_cursensor)];
            periodmins_peakregion_curevent = [periodmins_peakregion_curevent; ...
                struct_of_tables(cursensor).events_tables.MinPeriod_PeakRegion(eventind_cursensor)];
            periodmaxs_peakregion_curevent = [periodmaxs_peakregion_curevent; ...
                struct_of_tables(cursensor).events_tables.MaxPeriod_PeakRegion(eventind_cursensor)];
            for prevsensor = sensors_with_event{curevent}
                eventind_prevsensor = find(struct_of_tables(prevsensor).eventIDs == curevid);
                prevsensor_trace = struct_of_tables(prevsensor).reconstructed_series{eventind_prevsensor};
                prevsensor_datetimes = struct_of_tables(prevsensor).datetime_vecs;
                prevsensor_start = struct_of_tables(prevsensor).events_tables.StartTime(eventind_prevsensor);
                prevsensor_end = struct_of_tables(prevsensor).events_tables.EndTime(eventind_prevsensor);
                
                earlieststart = min(prevsensor_start, cursensor_start);
                latestend = max(prevsensor_end, cursensor_end);
                
                alt_rec_forcorr = prevsensor_trace(prevsensor_datetimes >= earlieststart & ...
                    prevsensor_datetimes <= latestend);
                cur_rec_forcorr = cursensor_trace(cursensor_datetimes >= earlieststart & ...
                    cursensor_datetimes <= latestend);
                
                % pad the arrays if necessary
                if max(cursensor_datetimes) < latestend
                    cur_rec_forcorr = [cur_rec_forcorr, zeros(1, seconds(latestend - max(cursensor_datetimes))/10)];
                end
                if min(cursensor_datetimes) > earlieststart
                    cur_rec_forcorr = [zeros(1, seconds(min(cursensor_datetimes) - earlieststart)/10), cur_rec_forcorr];
                end
                if max(prevsensor_datetimes) < latestend
                    alt_rec_forcorr = [alt_rec_forcorr, zeros(1, seconds(latestend - max(prevsensor_datetimes))/10)];
                end
                if min(prevsensor_datetimes) > earlieststart
                    alt_rec_forcorr = [zeros(1, seconds(min(prevsensor_datetimes) - earlieststart)/10), alt_rec_forcorr];
                end
                
                [crosscorr_rec, lags] = xcorr(alt_rec_forcorr, cur_rec_forcorr, 'normalized');
                % in the event of a 'tie' regarding optimal delay time for
                % maximizing the lagged correlation, take an average
                timedeltas{curevent} = [timedeltas{curevent}; ...
                    mean(10*lags(crosscorr_rec == max(crosscorr_rec)))];
                corr_coefs{curevent} = [corr_coefs{curevent}; max(crosscorr_rec)];
                sensor_pairs{curevent} = [sensor_pairs{curevent}; ...
                    [struct_of_tables(cursensor).sensor_id '-' struct_of_tables(prevsensor).sensor_id]];
                
                [dists_i, azimuths_i] = distance(latvec(cursensor),lonvec(cursensor),...
                    latvec(prevsensor),lonvec(prevsensor), referenceEllipsoid('GRS80'));
                xdists_i = dists_i .* sind(azimuths_i);
                ydists_i = dists_i .* cosd(azimuths_i);
                
                distancematrix = [distancematrix; [xdists_i,ydists_i]];
                
            end
            sensors_with_event{curevent} = [sensors_with_event{curevent}, cursensor];
            
        end
    end
    mean_corrs(curevent) = mean(corr_coefs{curevent});
    mean_amp(curevent) = mean(amplitudes_curevent);
    event_start(curevent) = min(starts_curevent);
    event_end(curevent) = max(ends_curevent);
    event_centerperiod(curevent) = mean(periods_curevent);
    event_minperiod(curevent) = min(periodmins_curevent);
    event_maxperiod(curevent) = max(periodmaxs_curevent);
    event_minperiod_peakregion(curevent) = min(periodmins_peakregion_curevent);
    event_maxperiod_peakregion(curevent) = max(periodmaxs_peakregion_curevent);
    N_sensors(curevent) = length(sensors_with_event{curevent});
    if length(timedeltas{curevent}) >= 3 % if we have enough sensors to do the slowness vector calculation
        % calculate the slowness vector (this will use a least squares
        % method for an overdetermined system of equations)
        slownessvec{curevent} = distancematrix \ timedeltas{curevent};
        wavespeed(curevent) = (sum(slownessvec{curevent}.^2))^(-0.5);
        wavedirection(curevent) = wrapTo360(atan2d(slownessvec{curevent}(1), slownessvec{curevent}(2)));
        waveu(curevent) = 1 / (slownessvec{curevent}(1));
        wavev(curevent) = 1 / (slownessvec{curevent}(2));
        
        % uncertainty estimate following Del Pezzo and Giudicepietro (2002)
        slowness_covmat = ((distancematrix' * distancematrix)^-1 * distancematrix') ...
            * eye(length(distancematrix))*20 * ...
            ((distancematrix' * distancematrix)^-1 * distancematrix')';
        error_x(curevent) = sqrt(slowness_covmat(1,1));
        error_y(curevent) = sqrt(slowness_covmat(2,2));
        
        % calculate what the delay times should be if the calculated wave
        % velocity is right
        modeled_timedeltas{curevent} = (1/wavespeed(curevent)) * ...
            sqrt(sum(distancematrix.^2, 2)).* ...
            cosd(atan2d(distancematrix(:,1), distancematrix(:,2)) - ...
            wavedirection(curevent));
        % from those, calculate the RMSE
        square_errors = (modeled_timedeltas{curevent} - timedeltas{curevent}) .^ 2;
        sum_square_errors = sum(square_errors);
        timedelta_RMSE(curevent) = sqrt(sum_square_errors / length(timedeltas{curevent}));
        % normalize the error by the delay time and calculate the root mean
        % square of that ('root mean square normalized error' or RMSNE)
        square_deltas = timedeltas{curevent} .^ 2;
        sum_square_deltas = sum(square_deltas);
        timedelta_NRMSE(curevent) = sqrt(sum_square_errors / sum_square_deltas);
    else  
        slownessvec{curevent} = [NaN; NaN];
    end
end

% put the outputs together into a table
event_table = table(uniqueIDs, N_sensors, mean_amp, event_start, ...
    event_end, event_centerperiod, event_minperiod, event_maxperiod, ...
    wavespeed, wavedirection, waveu, wavev, error_x, error_y, ...
    timedelta_RMSE, timedelta_NRMSE, mean_corrs, ...
    event_minperiod_peakregion, event_maxperiod_peakregion, ...
    'VariableNames', {'EventID', 'N_Sensors', 'Mean_Amplitude_hPa', ...
    'Event_Start', 'Event_End', 'Center_Period_min', 'Min_Period_min', ...
    'Max_Period_min', 'Wave_Speed_ms', 'Wave_Direction_deg', ...
    'WaveSpeed_x_ms', 'WaveSpeed_y_ms', 'Uncertainty_x', ...
    'Uncertainty_y', 'RMSE_seconds', 'NRMSE_unitless', ...
    'Mean_CrossCorr', 'MinPeriod_PeakRegion_min', ...
    'MaxPeriod_PeakRegion_min'});
    
end

