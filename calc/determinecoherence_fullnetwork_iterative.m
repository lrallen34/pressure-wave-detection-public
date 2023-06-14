function [struct_of_tables, wt_struct] = ...
    determinecoherence_fullnetwork_iterative(datadir, datetime_span, ...
    sensor_ids, threshloc, coef, crosscorr_thresh, maxdelay, smoothParam)
%DETERMINECOHERENCE_FULLNETWORK_iterative produces the tables output by
%determinecoherence_iterative.m for an entire network of pressure sensors, 
%iterating through each sensor as the 'primary' sensor
%
%   INPUTS:
%
%       datadir (str): path pointing to location of pressure sensor .csv
%       log files
%
%       datetime_span (datetime): vector (length 2) of datetimes indicating
%       start and end of period to be analyzed
%
%       sensor_ids (cell): cell array containing strings which correspond
%       to the IDs of sensors in the network
%
%       threshloc (str): path pointing to location of .txt file with
%       wavelet energy means by period for the sensor network
%
%       coef (scalar numeric): scalar indicating the amount to multiply the
%       wavelet power means by to get the threshold for connected regions 
%       to be considered events (double this number for event center 
%       threshold)
%
%       crosscorr_thresh (scalar numeric): value between 0 and 1 which
%       defines the correlation threshold between a coherent event and an 
%       incoherent event between two sensors (recommended values between
%       0.5 and 0.7)
%
%       maxdelay (scalar numeric): the maximum allowable delay time between
%       sensors in hours
%
%       smoothParam (scalar numeric): averaging period over which to smooth
%       original pressure trace in seconds
%
%   OUTPUTS:
%
%       struct_of_tables (struct): structure array containing each of the
%       output tables from determinecoherence.m using the given parameters.
%       each row of the structure array corresponds to using each of the
%       given sensor IDs as the primary sensor. all of the sensor IDs are
%       input as alternate sensors each time determinecoherence.m is ran.
%
%       wt_struct (struct): structure array of wavelet transform arrays for
%       each sensor
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Mar 2022
%

threshmat = readmatrix(threshloc);
period_vec = minutes(threshmat(:, 1));
mean_vec = threshmat(:, 2);
threshold_vec = mean_vec * coef;

struct_of_tables = struct();

% Identify events and get the pressure traces from each sensor
% create cell arrays of the event tables, event traces, datetime vectors,
% wavelet transforms, and full pressure traces
wt_struct = struct();
parfor cur_sensor = 1:length(sensor_ids)
    wt_struct(cur_sensor).sensor_id = sensor_ids{cur_sensor};
    struct_of_tables(cur_sensor).sensor_id = sensor_ids{cur_sensor};
    if datapresent(datadir, sensor_ids{cur_sensor}, datetime_span)
        disp(['identifying events for ' sensor_ids{cur_sensor}])
        [struct_of_tables(cur_sensor).events_tables, ~, ~, ~, ~, ...
            struct_of_tables(cur_sensor).event_bounds, ...
            struct_of_tables(cur_sensor).reconstructed_series, ...
            struct_of_tables(cur_sensor).datetime_vecs, ...
            wt_struct(cur_sensor).values, ...
            struct_of_tables(cur_sensor).fullTrace] = ...
            identifyevents_iterative(datadir,datetime_span,...
            sensor_ids{cur_sensor},smoothParam,period_vec,threshold_vec,0);
    end
end

% loop through each sensor and find coherent events
for cur_sensor = 1:length(sensor_ids)
    if datapresent(datadir, sensor_ids{cur_sensor}, datetime_span)
        
        disp(['matching events for ' sensor_ids{cur_sensor} ' as primary sensor'])
        prim_eventtable = struct_of_tables(cur_sensor).events_tables;
        prim_eventtraces = struct_of_tables(cur_sensor).reconstructed_series;
        prim_datetimevec = struct_of_tables(cur_sensor).datetime_vecs;
        % if data exists, but not enough data to process, skip the current
        % sensor
        if isempty(prim_datetimevec); continue; end

        [struct_of_tables(cur_sensor).table_allsensors, ...
            struct_of_tables(cur_sensor).coherencetable, ...
            struct_of_tables(cur_sensor).corrtable, ...
            struct_of_tables(cur_sensor).delaytable, ...
            struct_of_tables(cur_sensor).starttable, ...
            struct_of_tables(cur_sensor).centertable, ...
            struct_of_tables(cur_sensor).endtable, ...
            struct_of_tables(cur_sensor).amptable] = ...
            determinecoherence_iterative(...
            prim_eventtable, prim_eventtraces, prim_datetimevec, sensor_ids, ...
            {struct_of_tables.events_tables}, ...
            {struct_of_tables.reconstructed_series}, ...
            {struct_of_tables.datetime_vecs}, datadir, ...
            crosscorr_thresh, maxdelay);
    end
end

end

