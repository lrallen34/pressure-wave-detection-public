function [Wave_azimuth,WaveSpeed, time_error] = waveFrontTriangulation(LatLon1,LatLon2,LatLon3,Time1,Time2,Time3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Locations as 2 element lat/lon arrays. e.g. Loc1=[35.808617, -78.768625];
% Times as datetime 
%
% Output in compass degrees and m/s
%
% Uses awful brute force approach.

D2_1 = distance(LatLon1(1),LatLon1(2),LatLon2(1),LatLon2(2),referenceEllipsoid('earth','m'));
D3_1 = distance(LatLon1(1),LatLon1(2),LatLon3(1),LatLon3(2),referenceEllipsoid('earth','m'));

A2_1 = azimuth(LatLon1(1),LatLon1(2),LatLon2(1),LatLon2(2),referenceEllipsoid('earth','m'));
A3_1 = azimuth(LatLon1(1),LatLon1(2),LatLon3(1),LatLon3(2),referenceEllipsoid('earth','m'));

T2_1 = seconds(Time2-Time1);
T3_1 = seconds(Time3-Time1);

A_range = linspace(0, 359.99, 1000);
A_range_speed = (D2_1./T2_1).*cosd(A_range);
A_range_speed(A_range_speed<0)=NaN;
Time_to_3 = (D3_1./A_range_speed).*cosd(A3_1-(A2_1+A_range));
[time_error,I]=nanmin(abs(T3_1-Time_to_3));%in seconds
Wave_azimuth = wrapTo360(A_range(I)+A2_1);
WaveSpeed = A_range_speed(I);
end

