function [ datetime_spans ] = identifydatastreaks( datadir, sensor_id )
%IDENTIFYDATASTREAKS Identifies stretches of consecutive days ('streaks')
% for which data are available from a given pressure sensor (sensor_id) in
% a given directory (datadir)
%
% INPUTS:
%
%   datadir (str): location of .csv files containing pressure sensor data
%
%   sensor_id (str): ID of pressure sensor
%
% OUTPUTS:
%
%   datetime_spans (datetime): n-by-2 array of datetimes representing the
%   days on which data were found for the input sensor_id. Each row
%   contains a start (00:00:00 on the first day) and end (23:59:59 on the
%   last day) to a streak of days with data.
%
% Written by Luke Allen (lrallen3@ncsu.edu)
%

sensorfiles = dir([datadir '*' sensor_id '*.csv']);

curfile = 1;
i = 1;
datetime_spans = datetime.empty(0,2);
while curfile <= length(sensorfiles)
    % extract the date and determine if the file for the same sensor
    % and the next day exists (while loop)
    fileExists = 1;
    curname = sensorfiles(curfile).name;
    datestart = regexp(curname, '\d{8}');
    streakstartstr = curname(datestart:datestart+7);
    streakstart = datetime(streakstartstr, 'InputFormat', 'yyyyMMdd');
    streaklen = 1;
    while fileExists == 1
        datestring = curname(datestart:datestart+7);
        nextday = datetime(datestring, 'InputFormat', 'yyyyMMdd') + days(1);
        nextstring = datestr(nextday, 'yyyymmdd');
        
        newname = strrep(curname, datestring, nextstring);
        if curfile == length(sensorfiles)
            fileExists = 0;
            curfile = curfile + 1;
        else
            nextname = sensorfiles(curfile+1).name;
            fileExists = strcmp(newname, nextname);
            curfile = curfile + 1;
            streaklen = streaklen + 1;
            curname = sensorfiles(curfile).name;
            datestart = regexp(curname, '\d{8}');
        end
    end
    % after while loop, extract the last date in the 'streak' of days
    % with data for the sensor
    
    % sensors 05 and 06 were moved; ignore dates before their moves
    if ~(strcmp(sensor_id, 'EA_PiZero_005') && nextday < datetime(2019,5,1)) ...
            && ~(strcmp(sensor_id, 'EA_PiZero_006') && nextday < datetime(2019,2,1))
        datetime_spans(i, :) = [streakstart, nextday-seconds(1)];
        i = i+1;
    end
end


end

