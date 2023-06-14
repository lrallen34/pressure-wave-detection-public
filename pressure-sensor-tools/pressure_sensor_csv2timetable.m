function TTout = pressure_sensor_csv2timetable(filePath, datetime_span, sensor_id)
%pressure_sensor_csv2timetable Read a pressure sensor CSV data log into a
%timetable
%   
%Inputs:
%   filePath - char - full file path to CSV file to be read in
%   datetime_span - array - 1x2 array of datetimes. The first value is the
%      start of the desired time span and the second is the end.
%   sensor_id - char - char of sensor id to use as a wildcard mask when 
%   searching the given directory for CSV files
%
%Outputs:
%   TTout - timetable - Data output as a timetable
%
%Notes: data are linearly interpolated and snapped to be on the second.
%
%Syntax:
%   TTout = pressure_sensor_csv2timetable(...
%      'C:\Users\mille\Desktop\PressureTestCase\BMP388_TP_Log_EA_PiZero_006_20201030.csv');
%
%   DateSpan = [datetime(2020,12,16,1,2,3),datetime(2020,12,18,11,12,33)];
%   TTout = pressure_sensor_csv2timetable(...
%      'C:\Users\mille\Desktop\PressureTestCase\, DateSpan, '003');
%
%Written By: Matthew Miller, Dec. 2020



% Determine if filePath is for a file or a directory and then build a
switch exist(filePath)
    case 2  % is a file
        list_of_files = string(filePath);
    case 7  % is a directory
        % build range of dates to pull files for
        date_list = dateshift(datetime_span(1),'start','day'):dateshift(datetime_span(2),'start','day');
        file_struct = struct([]);
        for date_mask = date_list
            % convert the datetime to a char to use as a wildcard mask
            date_mask_char = datestr(date_mask, 'yyyymmdd');
            % search the directory for the file(s)
            file_struct_temp = dir([filePath '*' sensor_id '_' date_mask_char '.csv']);
            % warn if no file is found.
            if isempty(file_struct_temp)
                warning('File for sensor id %s on date %s not found.\n', sensor_id, date_mask_char);
            end
            % concatenate file structs
            file_struct = [file_struct; file_struct_temp];
        end
        list_of_files = string();
        % build cell array of chars of full file paths
        for ii = 1:length(file_struct)
            list_of_files(ii) = fullfile(file_struct(ii).folder, file_struct(ii).name);
        end
                
    otherwise
        error('filePath does not exist or was not recognized');
end


TTout = timetable();
for csv_file_path = list_of_files
    % Parse file name to identify sensor type in order to know what data is
    % inside
    [~, fileName, ~] = fileparts(csv_file_path);
    [regexout] = regexp(fileName,'^(\w{6})_[TPH]{2,3}_Log_(\w*)_\d{8}','tokens');
    [sensorName, sensorID] = regexout{1}{:};
    
    switch sensorName
        case 'BME280'
            % Set readtimetableoptions
            opts = delimitedTextImportOptions(...
                'ExtraColumnsRule', 'ignore', ...
                'VariableNames', {'Datetime', 'Temperature', 'Pressure', 'Humidity'},...
                'VariableTypes', {'datetime', 'double', 'double', 'double'});
            opts = setvaropts(opts, 'Datetime', 'InputFormat', 'yyyyMMdd HHmmss.SSS');
        case 'BMP388'
            % Set readtimetableoptions
            opts = delimitedTextImportOptions(...
                'ExtraColumnsRule', 'ignore', ...
                'VariableNames', {'Datetime', 'Temperature', 'Pressure'},...
                'VariableTypes', {'datetime', 'double', 'double'});
            opts = setvaropts(opts, 'Datetime', 'InputFormat', 'yyyyMMdd HHmmss.SSS');
    end
    TTout_temp = readtimetable(csv_file_path, opts);
    
    % remove obviously unphysical pressure values
    TTout_temp(TTout_temp.Pressure > 1100 | TTout_temp.Pressure < 800, :) = [];
    
    % if there is a 'jump back' in time, the retiming will not work.
    % therefore, delete the subsequent rows until we get back to times
    % after the jump
    TTout_temp = TTout_temp(~isnat(TTout_temp.Datetime), :);
    while min(diff(TTout_temp.Datetime)) < 0
        jumpind = find(diff(TTout_temp.Datetime) < 0, 1, 'first');
        jumptime = TTout_temp.Datetime(jumpind);
        lastind_todelete = find(TTout_temp.Datetime < jumptime, 1, 'last');
        TTout_temp(jumpind:lastind_todelete, :) = [];
    end
    
    TTout = [TTout; TTout_temp];
end

% retime the table to make data fall exactly on the second
TTout = retime(TTout, 'secondly', 'linear');

% subset based on given time
TTout = TTout(timerange(datetime_span(1), datetime_span(2), 'closed'), :);

end

