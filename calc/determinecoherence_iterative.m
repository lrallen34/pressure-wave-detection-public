function [full_eventstable, coherencetable, corrtable, delaytable, ...
    starttable, centertable, endtable, amptable] = ...
    determinecoherence_iterative(prim_eventtable, prim_eventtraces, ...
    prim_datetimevec, alt_sensors, alt_eventstables, alt_eventtraces, ...
    alt_datetimevecs, datadir, crosscorr_thresh, maxdelay)
%DETERMINECOHERENCE_iterative identifies events for a pressure time series 
%from a set of sensors(using identifyevents_iterative.m), determines 
%whether those events can also be identified in a set of 'alternate'
%pressure sensors with some time lag, then outputs a table of events in the
%original sensor with information on whether those events occur in each 
%alternate sensor and with what time lag
%
%   INPUTS:
%
%       prim_eventtable (table): table output by identifyevents_iterative.m
%       containing information on the detected wave events for the primary
%       sensor
%
%       prim_eventtraces (cell): cell array containing vectors
%       corresponding to the extracted waveforms of events in the primary
%       sensor, output by identifyevents_iterative.m
%
%       prim_datetimevec (datetime vector): datetimes corresponding to the
%       vectors in prim_eventtraces, output by identifyevents_iterative.m
%
%       alt_sensors: cell array (or a single string) of alternate sensor
%       IDs
%
%       alt_eventstables (cell): cell array of event tables, like
%       prim_eventtable, but for the 'alternate' sensors
%
%       alt_eventtraces (cell): cell array with nested cell arrays
%       containing vectors corresponding to the extracted waveforms of
%       events, like prim_eventtraces, but for the 'alternate' sensors
%
%       alt_datetimevecs (cell): cell array of datetime vectors
%       corresponding to the vectors in alt_eventtraces, like
%       prim_datetimevec, but for the 'alternate' sensors
%
%       datadir: string indicating the location of pressure sensor CSV data
%       files
%
%       crosscorr_thresh (num): the correlation threshold between sensors
%       to call an event 'coherent' (recommended values between 0.5 and
%       0.7)
%
%       maxdelay (num): the maximum allowable delay between sensors in
%       hours
%
%   OUTPUTS:
%       full_eventstable: a table containing information on events
%       identified from the input prim_wt (including the central
%       time, central period, and ranges on the time and period for each
%       event). The table also includes columns indicating whether an event
%       coherent with the event for each row was identified in each
%       alternative sensor, and the time lag between the event in the
%       original sensor and the event in the alternate sensor.
%
%       coherencetable, corrtable, delaytable, starttable, centertable,
%       endtable, amptable: tables with the same number of rows as new_eventstable,
%       containing one column for each sensor indicating whether events
%       were found coherent in that sensor, the correlation between the
%       event in the primary sensor and the event in that sensor, the delay
%       between the event in the primary sensor and the event in that
%       sensor, and the start, center, end time of events in that
%       sensor, and amplitude of wave events (in hPa),  respectively
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Aug 2022
%

% for each alternate sensor read in the data and get the
% wavelet transform, then for each event, determine whether it is present
% in the alternate sensor
% preallocating arrays to go into output table
eventappears = strings(height(prim_eventtable), length(alt_sensors));
delaytime = NaN(height(prim_eventtable), length(alt_sensors));
event_starts = NaT(height(prim_eventtable), length(alt_sensors));
alt_evcenters = NaT(height(prim_eventtable), length(alt_sensors));
event_ends = NaT(height(prim_eventtable), length(alt_sensors));
eventcorrelation = NaN(height(prim_eventtable), length(alt_sensors));
event_amplitude = NaN(height(prim_eventtable), length(alt_sensors));
for cur_sensor = 1:length(alt_sensors)

    alt_id = alt_sensors{cur_sensor};
    disp(['alt sensor: ' alt_id ' (' num2str(cur_sensor) '/' ...
        num2str(length(alt_sensors)) ')']);
    new_datetime_span = [min(prim_datetimevec)-hours(10), max(prim_datetimevec)+hours(10)];

    % make sure the data exist before proceeding
    dataindicator = datapresent(datadir, alt_id, new_datetime_span);
    if dataindicator == 0
        eventappears(:, cur_sensor) = "nodata";
        continue
    end

    % get events from the alternate sensor
    curalttable = alt_eventstables{cur_sensor};
    curalttraces = alt_eventtraces{cur_sensor};
    curalt_datetimevec = alt_datetimevecs{cur_sensor};
    
    % if no events were found in the alternate sensor, there is no need to
    % proceed for that sensor
    if isempty(curalttable)
        eventappears(:, cur_sensor) = "incoherent";
        continue
    end

    potentialmatch = possiblematchingevents(prim_eventtable, curalttable, maxdelay);

    for cur_event = 1:height(prim_eventtable)
        % to avoid edge effects, need 3*pi*tau*sqrt(2) on either side of the
        % event where tau is the max period -- go an extra 10 hours to account for
        % wave propagation speed
        cur_evstart = prim_eventtable.StartTime(cur_event);
        cur_evend = prim_eventtable.EndTime(cur_event);

        cur_rec = prim_eventtraces{cur_event};
        possible_alttraces = curalttraces(potentialmatch{cur_event});

        optimal_crosscorr_final = nan(1, length(possible_alttraces));
        possible_delay = nan(1, length(possible_alttraces));
        alt_evstart = curalttable.StartTime(potentialmatch{cur_event});
        alt_evend = curalttable.EndTime(potentialmatch{cur_event});
        alt_amplitudes = nan(1, length(possible_alttraces));
        for cur_alt_event = 1:length(possible_alttraces)
            % get the reconstructed event from the alternate sensor
            alt_rec = possible_alttraces{cur_alt_event};
            earlieststart = min(alt_evstart(cur_alt_event), cur_evstart);
            latestend = max(alt_evend(cur_alt_event), cur_evend);

            alt_rec_forcorr = alt_rec(curalt_datetimevec >= earlieststart & curalt_datetimevec <= latestend);
            cur_rec_forcorr = cur_rec(prim_datetimevec >= earlieststart & prim_datetimevec <= latestend);

            % pad the arrays if necessary
            if max(prim_datetimevec) < latestend
                cur_rec_forcorr = [cur_rec_forcorr, zeros(1, seconds(latestend - max(prim_datetimevec))/10)];
            end
            if min(prim_datetimevec) > earlieststart
                cur_rec_forcorr = [zeros(1, seconds(min(prim_datetimevec) - earlieststart)/10), cur_rec_forcorr];
            end
            if max(curalt_datetimevec) < latestend
                alt_rec_forcorr = [alt_rec_forcorr, zeros(1, seconds(latestend - max(curalt_datetimevec))/10)];
            end
            if min(curalt_datetimevec) > earlieststart
                alt_rec_forcorr = [zeros(1, seconds(min(curalt_datetimevec) - earlieststart)/10), alt_rec_forcorr];
            end

            [crosscorr_rec, lags] = xcorr(alt_rec_forcorr, cur_rec_forcorr, 'normalized');
            optimal_crosscorr_final(cur_alt_event) = max(crosscorr_rec);
            % in the event of a 'tie' regarding optimal delay time for
            % maximizing the lagged correlation, take an average
            possible_delay(cur_alt_event) = mean(10*lags(crosscorr_rec == max(crosscorr_rec)));
            alt_amplitudes(cur_alt_event) = max(alt_rec) - min(alt_rec);
        end
        % choose the event that has the best correlation with the original
        % event, if there are more than one
        if ~isempty(possible_alttraces)
            max_event_correlation = max(optimal_crosscorr_final);
            event_starts(cur_event, cur_sensor) = ...
                alt_evstart(optimal_crosscorr_final == max_event_correlation);
            alt_evcenters(cur_event, cur_sensor) = ...
                curalttable.CenterTime(potentialmatch{cur_event}...
                (optimal_crosscorr_final == max_event_correlation));
            event_ends(cur_event, cur_sensor) = ...
                alt_evend(optimal_crosscorr_final == max_event_correlation);
            eventcorrelation(cur_event, cur_sensor) = max_event_correlation;
            delaytime(cur_event, cur_sensor) = ...
                possible_delay(optimal_crosscorr_final == max_event_correlation);
            event_amplitude(cur_event, cur_sensor) = ...
                alt_amplitudes(optimal_crosscorr_final == max_event_correlation);
            % define the cross-correlation for a 'coherent' event
            if optimal_crosscorr_final < crosscorr_thresh
                eventappears(cur_event, cur_sensor) = "incoherent";
            else
                eventappears(cur_event, cur_sensor) = "coherent";
            end

        else % if there are no possible events to match, note as much
            eventappears(cur_event, cur_sensor) = "incoherent";

        end
    end
end

% make tables of whether events appear in each sensor, the event
% correlations, and the event delay times, and concatenate them to the
% events table
appearnames = cell(size(alt_sensors));
corrnames = cell(size(alt_sensors));
delaynames = cell(size(alt_sensors));
for cursensor = 1:length(alt_sensors)
    appearnames{cursensor} = ['appearsin_' alt_sensors{cursensor}];
    corrnames{cursensor} = ['correlationwith_' alt_sensors{cursensor}];
    delaynames{cursensor} = ['delayto_' alt_sensors{cursensor}];
end

coherencetable = array2table(eventappears, 'VariableNames', appearnames);
corrtable = array2table(eventcorrelation, 'VariableNames', corrnames);
delaytable = array2table(delaytime, 'VariableNames', delaynames');
starttable = array2table(event_starts, 'VariableNames', alt_sensors');
centertable = array2table(alt_evcenters, 'VariableNames', alt_sensors');
endtable = array2table(event_ends, 'VariableNames', alt_sensors');
amptable = array2table(event_amplitude, 'VariableNames', alt_sensors');
full_eventstable = [prim_eventtable, coherencetable, corrtable, delaytable];

end
