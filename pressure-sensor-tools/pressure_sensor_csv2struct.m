function structOut = pressure_sensor_csv2struct(filePath)
%pressure_sensor_csv2struct Read a pressure sensor CSV data log into a
%structure
%   
%Inputs:
%   filePath - string - full file path to CSV file to be read in
%
%Outputs:
%   TTout - timetable - Data output as a struct
%
%Syntax:
%   structOut = pressure_sensor_csv2struct(...
%      'C:\Users\mille\Desktop\PressureTestCase\BMP388_TP_Log_EA_PiZero_006_20201030.csv');
%
%Written By: Matthew Miller, Dec. 2020

% Parse file name to identify sensor type in order to know what data is
% inside
[~, fileName, ~] = fileparts(filePath);
[regexout] = regexp(fileName,'^(\w{6})_[TPH]{2,3}_Log_(\w*)_\d{8}','tokens');
[sensorName, sensorID] = regexout{1}{:};

% use textscan to read the file and sort the output cells into structure
% fields
fid = fopen(filePath);
switch sensorName
    case 'BME280'
        C = textscan(fid,'%{yyyyMMdd HHmmss.SSS}D%f%f%f','Delimiter',',');
        structOut.datetime = C{1};
        structOut.temperature = C{2};
        structOut.pressure = C{3};
        structOut.humidity = C{4};
    case 'BMP388'
        C = textscan(fid,'%{yyyyMMdd HHmmss.SSS}D%f%f','Delimiter',',');
        structOut.datetime = C{1};
        structOut.temperature = C{2};
        structOut.pressure = C{3};
end
fclose(fid);
structOut.SensorID = sensorID;

end

