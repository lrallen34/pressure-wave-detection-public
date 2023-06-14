# Pressure Sensor Tools

## Reading the data
Use `pressure_sensor_csv2timetable.m` to read the data into Matlab as a timetable object. 
This function reads the raw CSV data then snaps all the data obsevations to be exactly on the sencond using linear interpolation. The function can either read a single file or return data from many files based on a user specified datetime span.

### Example Syntax
Single file syntax
```matlab
TTout = pressure_sensor_csv2timetable(...
      'C:\Users\mille\Desktop\PressureTestCase\BMP388_TP_Log_EA_PiZero_006_20201030.csv');
```

Datetime range syntax
```matlab
DateSpan = [datetime(2020,12,16,1,2,3),datetime(2020,12,18,11,12,33)];
TTout = pressure_sensor_csv2timetable(...
      'C:\Users\mille\Desktop\PressureTestCase\, DateSpan, '003');
```


## About the Sensor Data

The data are located at `/home/disk/ivanova2/RPi_Pressure_Data/`

### Sensor Data Format
#### File naming scheme

`<Sensor Type>_Log_<Sensor Name>_<Date>.csv`

e.g. `BME280_TPH_Log_EA_PiZero_002_20190324.csv`

Sensor Type - BME280_TPH or BMP388_TP (TPH indicates the data files contain temperature, pressure, and humidity data; TP indicates data files contain temperature and pressure data)

Sensor Name - The name of the sensor e.g. `EA_PiZero_002` or `eapizero025`

Date - `yyyymmdd` format

#### CSV Format Details

The data are in CSV format with no header line.

The time format is `yyyymmdd HHMMSS.sss`

The data are ordered: `datetime, temperature, pressure, humidity (if present)`

## About the Sensors
The pressure sensors are Bosch [BME280](https://www.adafruit.com/product/2652) and [BMP388](https://www.adafruit.com/product/3966) sensors. 
The sensors are attached to [Raspberry Pi Zero](https://www.raspberrypi.org/products/raspberry-pi-zero-w/) single board computers which run the [Raspbian OS](https://www.raspberrypi.org/documentation/raspbian/).
The data are collected at roughly 1 Hz.



