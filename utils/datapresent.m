function data_indicator = datapresent(dataloc, sensor_id, datetime_span)
% DATAPRESENT Check if pressure sensor data exist within a given time span
% for a given sensor ID
%   
%   INPUTS:
%   
%       dataloc (str): location of pressure data .csv files
%
%       sensor_id (str): ID of sensor to check
%
%       datetime_span (datetime): vector (length 2) of datetimes defining
%       start and end time of span to check
%
%   OUTPUTS:
%
%       data_indicator (num): 1 if any data exist for the given sensor_id
%       within the given datetime_span, 0 otherwise
%

data_indicator = 0;
[y(1), m(1), d(1)] = ymd(datetime_span(1));
[y(2), m(2), d(2)] = ymd(datetime_span(2));
daysvec = datetime(y(1), m(1), d(1)):days(1):datetime(y(2), m(2), d(2));
for curday = 1:length(daysvec)
    curdaystr = char(daysvec(curday), 'yyyymmdd');
    filestr = [dataloc '/*' sensor_id '_' curdaystr '.csv'];
    if ~isempty(dir(filestr))
        data_indicator = 1;
        break
    end
end
end

