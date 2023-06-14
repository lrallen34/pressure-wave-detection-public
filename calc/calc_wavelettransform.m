function [wt, datetime_vec, taumesh, coimesh, smoothed_pressuretrace] = ...
    calc_wavelettransform(csv_directory, datetime_span, sensor_id, ...
    period_vec, smoothParam, varargin)
% CALC_WAVELETTRANSFORM Performs the wavelet transform of the pressure data
% for a given sensor over a given span of periods and datetimes. The data
% are smoothed using a moving average over a window defined by smoothparam.
% An analytic morlet wavelet is currently used.
%
%   INPUTS:
%   
%       csv_directory (str): the directory where pressure sensor data (CSV
%       files) are located
%
%       datetime_span (datetime): vector (length 2) of datetimes defining
%       the start and end time to perform the wavelet transform over
%
%       sensor_id (str): ID of the sensor for which to read and transform
%       the data
%
%       period_vec (num): vector of periods in seconds defining the
%       range of periods to include
%
%       smoothParam (num): duration of moving average window in seconds
%
%       outLoc (str, optional): location of output directory if saving wavelet
%       transform desired
%
%   OUTPUTS:
%
%       wt (complex double): wavelet transform array. Values outside the
%       cone of influence are set to NaN
%
%       datetime_vec (datetime): vector of datetimes corresponding to the
%       columns of wt
%
%       taumesh (duration): array of periods corresponding to each element
%       in wt
%
%       coimesh (duration): array of periods which define the cone of
%       influence in each column, corresponding to each element in wt
%
%       smoothed_pressuretrace (double): the vector of pressure values
%       smoothed using a moving average with window defifned by smoothParam
%
%   Written by Luke Allen (lrallen3@ncsu.edu)
%

presData = pressure_sensor_csv2timetable(csv_directory, datetime_span, sensor_id);
% start at the first time step evenly divisible by the smoothing parameter
secplace = second(presData.Datetime);
firstind = find(rem(secplace, smoothParam) == 0, 1);
presData = presData(firstind:end, :);


% this will not work if the timetable is too short, so if less than 3600
% rows (an hour of data), output empty arrays and warn
if height(presData) < 3600
    [wt, datetime_vec, taumesh, coimesh, smoothed_pressuretrace] = deal([]);
    warning('not enough data in this span, returning empty arrays...')
    return
end

smoothedtable = smoothdata(presData, 'movmean', seconds(smoothParam), 'omitnan');
smoothedtable = smoothedtable(ceil((smoothParam/2)):smoothParam:end, :);
smoothed_pressuretrace = smoothedtable.Pressure;

% perform wavelet transform
[wt, period_vec, coi] = cwt(smoothedtable.Pressure, 'amor', minutes(smoothParam/60), ...
    'PeriodLimits', [min(period_vec) max(period_vec)+minutes(1)]);

% change values along the edges of the data to NaN
[coimesh, taumesh] = meshgrid(coi, period_vec);
wt(taumesh >= coimesh) = NaN;
datetime_vec = smoothedtable.Datetime;

% check if save location input
if ~isempty(varargin) && ischar(varargin{1})
    makeneededdirs([varargin{1} '/' sensor_id])
    outfile = [varargin{1} '/' sensor_id '/' ...
        char(datetime_span(1), 'yyyymmddHH') 'to' ...
        char(datetime_span(2), 'yyyymmddHH') '_wtstruct.mat'];
    save(outfile, 'wt', 'datetime_vec', 'period_vec', '-mat');
elseif ~isempty(varargin) && ~ischar(varargin{1})
    error('output location should be a string')
end

end
