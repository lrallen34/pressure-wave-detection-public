function [eventcatalog, eventtable] = create_eventcatalog(...
    struct_of_tables, sensor_ids, wave_propertytable, ...
    sensors_with_event, timedeltas, modeled_timedeltas, ...
    corr_coefs, sensor_pairs)
%CREATE_EVENTCATALOG creates a structure array 'event catalog' from a
%structure of tables output by assign_eventIDs.m, and vectors of sensor IDs,
%latitudes, and longitudes. the catalog includes information on which
%sensors had each event, the optimal 'central' or 'primary' sensor for the
%event, correlation coefficients between the central sensor and the other
%sensors for the event, and the propagation speed, direction, and error if
%at least three sensors caught the event.
%
%   INPUTS:
%
%       struct_of_tables (struct): structure output by assign_eventIDs.m
%       which contains a row for each sensor and fields with the sensor ID,
%       tables of events and how they correspond to the other sensors in
%       the array, and vectors indicating the event IDs assigned to each
%       event found for each sensor
%
%       sensor_ids (cell): cell of strings with each sensor ID
%
%       sensor_lats (numeric vector): latitudes of each sensor, in the same
%       order as sensor_ids
%
%       sensor_lons (numeric vector): longitudes of each sensor, in the
%       same order as sensor_ids
%       
%       (each of the following inputs are outputs from
%       calc_slownessvector.m)
%       wave_propertytable (table): table output from calc_slownessvector.m
%       with wave event modeled velocity components and error estimates
%
%       sensors_with_event (cell): each cell represents an event; the
%       elements within the cell are indices corresponding to the sensors
%       which captured the event (those indices align with
%       struct_of_tables)
%
%       timedeltas (cell): each cell represents an event; the elements
%       within the cell are observed delay times between each possible pair
%       of sensors which captured the event (those pairs are enumerated in
%       sensor_pairs)
%
%       modeled_timedeltas (cell): same as timedeltas, but with the delay
%       times which would be produced if the modeled wave velocity was
%       exactly right and constant
%
%       corr_coefs (cell): each cell represents an event; the elements
%       within the cell are max of the extracted event cross-correlations
%       between each possible pair of sensors which captured the event 
%       (those pairs are enumerated in sensor_pairs)
%
%       sensor_pairs (cell): each cell represents an event; the elements
%       within the cell are each possible pair of sensors which captured
%       the event, represented as a string array
%
%   OUTPUTS:
%
%       eventcatalog (struct): structure array with a row for each unique
%       event ID in struct_of_tables and fields for the event's central
%       sensor, alternate sensors, correlation coefficients, delay times,
%       wave propagation speed and direction, error in the delay times
%       given that speed/direction, and the allowed error to consider the
%       event coherent
%
%       eventtable (table): table with a row for each unique event ID in
%       struct_of_tables that corresponds to a coherent event, and columns
%       for the event's center date/time, center period, center sensor ID,
%       wave amplitude, number of secondary sensors, mean cross-correlation
%       with alternate sensors, modeled wave speed/direction, model delay
%       error, and the allowed error to consider the event coherent
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Mar 2022 (last updated Dec
%   2022)
%

% preallocate output arrays
eventcatalog = struct();
eventtable = table('Size', [0 17], ...
    'VariableNames', {'center_period', 'start_datetime', 'end_datetime', ...
    'min_period', 'max_period', 'amplitude_hPa', 'num_sensors', ...
    'mean_crosscorr', 'min_crosscorr', 'modeled_wavespeed_ms', ...
    'modeled_wavedir', 'delay_RMSE', 'delay_NRMSE', 'uncertainty_x', ...
    'uncertainty_y', 'minperiod_peakregion', 'maxperiod_peakregion'}, ...
    'VariableTypes', {'duration', 'datetime', 'datetime', 'duration', ...
    'duration', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'duration', ...
    'duration'});

%disp('Event summaries:')
n_coherent = 1;
allIDs = {struct_of_tables.eventIDs};
concatenatedIDs = {cat(1, allIDs{:})};
uniqueIDs = unique(concatenatedIDs{1});
% loop through each unique event
for curevent = 1:length(uniqueIDs)
    eventcatalog(curevent).event_id = curevent;
    eventcatalog(curevent).mean_amplitude = wave_propertytable.Mean_Amplitude_hPa(curevent);
    % assign values from table to catalog for the event
    eventcatalog(curevent).start_time = wave_propertytable.Event_Start(curevent);
    eventcatalog(curevent).end_time = wave_propertytable.Event_End(curevent);
    eventcatalog(curevent).center_period = wave_propertytable.Center_Period_min(curevent);
    eventcatalog(curevent).min_period = wave_propertytable.Min_Period_min(curevent);
    eventcatalog(curevent).max_period = wave_propertytable.Max_Period_min(curevent);
    eventcatalog(curevent).minperiod_peakregion = wave_propertytable.MinPeriod_PeakRegion_min(curevent);
    eventcatalog(curevent).maxperiod_peakregion = wave_propertytable.MaxPeriod_PeakRegion_min(curevent);
    eventcatalog(curevent).sensors_with_event = {sensor_ids{sensors_with_event{curevent}}};
    eventcatalog(curevent).sensor_pairs = sensor_pairs{curevent};
    eventcatalog(curevent).cross_correlations = corr_coefs{curevent};
    eventcatalog(curevent).delay_times = timedeltas{curevent};
    eventcatalog(curevent).wave_speed = wave_propertytable.Wave_Speed_ms(curevent);
    eventcatalog(curevent).wave_direction = wave_propertytable.Wave_Direction_deg(curevent);
    eventcatalog(curevent).modeled_delaytimes = modeled_timedeltas{curevent};
    eventcatalog(curevent).delay_error = wave_propertytable.RMSE_seconds(curevent);
    eventcatalog(curevent).delay_error_norm = wave_propertytable.NRMSE_unitless(curevent);
    eventcatalog(curevent).uncertainty = [wave_propertytable.Uncertainty_x(curevent), ...
        wave_propertytable.Uncertainty_y(curevent)];
    
    % for events appearing in at least 3 sensors, add them to an output
    % table
    if length(sensors_with_event{curevent}) >= 3 
        eventtable.start_datetime(n_coherent) = eventcatalog(curevent).start_time;
        eventtable.end_datetime(n_coherent) = eventcatalog(curevent).end_time;
        eventtable.center_period(n_coherent) = eventcatalog(curevent).center_period;
        eventtable.min_period(n_coherent) = eventcatalog(curevent).min_period;
        eventtable.max_period(n_coherent) = eventcatalog(curevent).max_period;
        eventtable.amplitude_hPa(n_coherent) = eventcatalog(curevent).mean_amplitude;
        eventtable.num_sensors(n_coherent) = length(eventcatalog(curevent).sensors_with_event);
        eventtable.mean_crosscorr(n_coherent) = wave_propertytable.Mean_CrossCorr(curevent);
        eventtable.min_crosscorr(n_coherent) = min(corr_coefs{curevent});
        eventtable.modeled_wavespeed_ms(n_coherent) = wave_propertytable.Wave_Speed_ms(curevent);
        eventtable.modeled_wavedir(n_coherent) = wave_propertytable.Wave_Direction_deg(curevent);
        eventtable.delay_RMSE(n_coherent) = wave_propertytable.RMSE_seconds(curevent);
        eventtable.delay_NRMSE(n_coherent) = wave_propertytable.NRMSE_unitless(curevent);
        eventtable.uncertainty_x(n_coherent) = wave_propertytable.Uncertainty_x(curevent);
        eventtable.uncertainty_y(n_coherent) = wave_propertytable.Uncertainty_y(curevent);
        eventtable.minperiod_peakregion(n_coherent) = wave_propertytable.MinPeriod_PeakRegion_min(curevent);
        eventtable.maxperiod_peakregion(n_coherent) = wave_propertytable.MaxPeriod_PeakRegion_min(curevent);
        
        n_coherent = n_coherent + 1;
        
    end
end

% sort events from earliest to latest according to their start time
eventtable = sortrows(eventtable, 'start_datetime');

end

