# pressure-wave-identification
Matlab code for identifying gravity wave events from 1-Hz pressure sensor data.

This code ingests data from Bosch BME280 and BMP288 pressure sensors deployed by Environment Analytics using functions from Matt Miller's pressure-sensor-tools repository, then can be used to identify wave events using methods similar to Grivet-Talocia and Einaudi (1998) and Grivet-Talocia et al. (1999). An analytic morlet wavelet transform is applied to the data, thresholds are defined using the average wavelet power by period over the entire data set, then events are identified from regions in the time-wave period phase space which meet the threshold. Events are identified across all sensors individually, then events in different sensors with overlapping period ranges that occur within 2 hours of each other are checked to determine if they represent the same wave(s). The time lag for the event from the 'primary' sensor to the 'secondary' sensor is obtained by maximizing the cross-correlation between the 'primary' waveform and the 'alternate' waveform. If that cross-correlation exceeds some threshold (currently 0.65), the event is considered 'coherent' for the two sensors in question. Once wave events are linked for each possible combination of primary and secondary sensors in a given network, a plane wave model is fitted to the observed time lags. If the fitting error (RMSE) is less than 90 seconds, and normalized fitting error (NRMSE) is smaller than 0.1, the event is considered robust. The current processing status is summarized in [this spreadsheet](https://docs.google.com/spreadsheets/d/1GArBndze9BQ7TJqtc1_c58A32tpUe7rRBVTBL-DQ6_0/edit#gid=0). Output .mat and .csv files are located in the EA Google Drive [here](https://drive.google.com/drive/folders/1RgDo4aCpBK4_25GcXLAl6MzzAEctfWeS?usp=sharing).

Dependencies:

`pressure-sensor-tools`, a repository containing basic reader functions for the pressure sensor data output by Environment Analytics' pressure sensors (written by Matt Miller, mamille4@ncsu.edu).
`cmocean`, a set of perceptually uniform colormaps for plotting. Located on the Matlab file exchange [here](https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps). Reference: Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True colors of oceanography: Guidelines for effective and accurate colormap selection. Oceanography 29(3):9–13. http://dx.doi.org/10.5670/oceanog.2016.66.
The [Wavelet Toolbox](https://www.mathworks.com/products/wavelet.html) within Matlab.

# Procedure for processing pressure data

First, clone the repository.

1. Calculate the mean wavelet transform amplitude by wave period across each sensor network (preprocessing step; will not need to be done every time)

      To do this, run `calcthresh.m`. You'll need to adjust the `datetoendon` variable to define the date range over which the wavelet power is averaged; that range goes from the earliest date with data to the value of `datetoendon`. You may also need to adjust `smoothParam`; this variable defines the length of the interval over which the data are averaged before processing (e.g., for `smoothParam = 10`, 10-second averages of pressure are used). This script will also output a plot of the mean wavelet power by wave period to a given directory.
      
2. Edit variables in `masterrunner_produceeventcatalog.m` as necessary

      The script is designed to process one month of data at a time for a single network, since the process is pretty memory-intensive. Therefore, you'll need to edit the `network`, `yr`, and `mnth` values each time you process data. `datadir` points to where the raw pressure data are located, and `threshloc` points to a file with mean wavelet transform amplitudes by period (the output from step 1). `coef` describes the coefficient in the event region threshold. Local maxima in the wavelet transform which exceed `2*coef*(the mean wavelet transform for the given wave period)` are event centers, while the connected regions to event centers exceeding `coef*(the mean wavelet transform by wave period)`, extended along the time axis until a local minimum in the wavelet transform is reached, are the event regions. Grivet-Talocia et al. (1998, 1999) used a `coef = 1`; I am currently using `coef = 5` to only find the most robust events. You can also adjust `crosscorr_thresh`, the minimum optimized cross-correlation between the event traces for two sensors to match those event traces, and `maxdelay`, the maximum gap between events in 2 sensors to possibly pair them together (in hours). `catalog_outdir` is where .mat files containing a `struct` array of every event and a `struct` array of tables desribing the events in each sensor will be output to. `spreadsheet_outdir` is where a spreadsheet of robust events will be output to. `outdir_plots` is where .pngs of figures will be output; these include the full and event pressure traces for each event, the wavelet transforms for each event, and the wavelet transforms normalized by the mean wavelet transform amplitude by period for each event. `sensorloc` should point to where the text file containing the list of sensors with their coordinates is. Everything else should populate automatically or not need to be changed.

3. Run `masterrunner_produceeventcatalog.m` to process a single month of pressure data for one network of sensors

      The script runs the following functions:

      a) `determinecoherence_fullnetwork_iterative`, which identifies and extracts events in each individual sensor and determines the matching events for each possible combination of sensor pairs (order dependent, i.e., the function checks both sensor A-to-sensor B and sensor B-to-sensor A). This function also calculates the lag times between the event passage for every pair of sensors based on optimizing the cross-correlation between their respective extracted event traces.

      b) `assign_eventIDs`, which iterates through each sensor in the network and its events to define event IDs for tracking events across multiple sensors. The following steps are done in each iteration: if on the first sensor, simply assign a new ID to each event. Otherwise, check whether the event can be matched two-ways (i.e., with the ordering of sensors kept the same or reversed) with any events in previous sensors. If so, and there is only one event ID assigned to those matched events, use that ID; if not, assign a new ID to that event. If there are events with different IDs matched to the current event, give the current event the ID associated with the highest cross-correlation with the current sensor. If there are now multiple events for the current sensor with the same ID, the one with the highest cross-correlation to an event with that ID in any previous sensor gets assigned that ID, then the process is restarted for the other event with that ID excluded.

      c) `calc_slownessvector`, which runs a least-squares approximation to calculate the slowness vector (a vector containing the inverse of the components of the wave phase velocity) for each wave event that was captured by at least 3 sensors. It also calculates the RMSE and NRMSE in the least-squares approximations and outputs that information to a table.

      d) `create_eventcatalog` reorganizes the events into a new `struct` array, with one row for each unique event ID. The catalog array includes columns with information on each event, such as which sensors captured the event, which sensor was chosen as the central sensor, what the delay times from the central sensor to each other sensor were, etc. A table is also created, which only contains events that were captured by 4 or more sensors with sufficiently low error in the propagation model. The table summarizes the event start/end times, central sensor, mean event trace cross-correlation, propagation speed/direction, etc. The catalog and table are output to .mat and .csv files, respectively.

      e) `plotevents_fromcatalog` outputs plots (as .png files) of the full pressure trace, extracted wave trace, wavelet transform amplitude, and normalized wavelet transform amplitude for events which occurred in at least a given number of sensors.

      The outputs from `create_eventcatalog` and `plotevents_fromcatalog` are saved as described in step 2.

Relevant references:

Del Pezzo, E. and Giudicepietro, F.: Plane wave fitting method for a plane, small aperture, short period seismic array: a MATHCAD program, Computers & geosciences, 28, 59–64, 2002.

Grivet-Talocia, S., and F. Einaudi, 1998: Wavelet Analysis of a Microbarograph Network. IEEE Transactions on Geoscience and Remote Sensing, 36(2), 418-433.

Grivet-Talocia, S., Einaudi, F., Clark, W. L., Dennett, R. D., Nastrom, G. D., & VanZandt, T. E. (1999). A 4-yr Climatology of Pressure Disturbances Using a Barometer Network in Central Illinois, Monthly Weather Review, 127(7), 1613-1629. https://doi.org/10.1175/1520-0493(1999)127%3C1613:AYCOPD%3E2.0.CO;2
