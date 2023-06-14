function [events_table, omega_prime, omega, event_centers, omega_pixels, ...
    event_bounds, reconstructed_series, datetime_vec, original_wavelet_transform, ...
    original_fullTrace] = ...
    identifyevents_iterative(csv_directory, datetime_span, sensor_id, ...
    smoothParam, period_vec, threshold_vec, max_events)
%IDENTIFYEVENTS_iterative will produce a table with an entry for each
%event. Columns will include event start time, event end time, and the
%range of periods in the event. Inputs are the wavelet transform calculated
%from a pressure time series, a vector of durations representing periods, 
%and a vector of thresholds of equal length to the period vector 
%representing the minimum wavelet energies to define pressure events.
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
%       smoothParam (num): duration of moving average window in seconds
%
%       period_vec (num): vector of periods in seconds defining the
%       range of periods to include
%
%       threshold_vec (numeric): vector of wavelet amplitude thresholds
%       corresponding to the rows in wavelet_transform. These thresholds
%       are for the broader event regions, whereas the event center
%       thresholds will be double these values.
%
%       max_events (scalar int): defines the maximum number of events to
%       iterate through before terminating the process. If set to 0, all
%       events meeting the threshold will be identified and included
%
%   OUTPUTS:
%
%       events_table (table): table where each row corresponds to an event.
%       Columns indicate the borders of the event region in time-period
%       space (start_time, end_time, min_period, max_period), and the 
%       location of the event center in time-period space (center_time,
%       center_period).
%
%       omega_prime (logical): array of the same size as wavelet_transform
%       indicating where connected regions exceeding threshold_vec and
%       containing an event center are located
%
%       omega (logical): array of the same size as wavelet_transform
%       indicating where the final event regions are
%
%       event_centers (int): indices of event centers in the
%       wavelet_transform array space
%
%       omega_pixels (cell): cell array where each cell contains indices of
%       event regions in the wavelet_transform array space
%
%       event_bounds (cell): cell array of points which trace the outline
%       of the event regions in the wavelet_transform array space
%
%       reconstructed_series (cell): cell array of vectors representing the
%       extracted wave event traces where the indices line up with those in
%       event_table
%
%       datetime_vec (datetime vector): vector of datetimes corresponding
%       to the extracted and full pressure traces
%
%       original_wavelet_transform (complex array): the first iteration of
%       the wavelet transform calculated on the pressure trace
%
%       original_fulltrace (vector double): the original pressure trace
%
%   Author: Luke Allen (lrallen3@ncsu.edu)
%   Last updated 2023/06/13
%

% read in the full pressure trace and calculate the wavelet transform
[wavelet_transform, datetime_vec, ~, ~, fullTrace] = calc_wavelettransform(...
    csv_directory, datetime_span, sensor_id, period_vec, smoothParam);
original_wavelet_transform = wavelet_transform;
original_fullTrace = fullTrace;

% find the absolute maximum of the wavelet magnitude divided by the
% time-averaged wavelet magnitude (i.e., the threshold)
wt_norm = abs(wavelet_transform) ./ threshold_vec(1:size(wavelet_transform,1));
[wt_max, ind_max] = max(wt_norm, [], 'all', 'linear');
[row_max, col_max] = ind2sub(size(wt_norm), ind_max);

% get indices where the wavelet transform is NaN (i.e., outside the cone of
% influence)
naninds = find(isnan(wavelet_transform));
wavelet_transform(isnan(wavelet_transform)) = 0;

% iterate through the following process:
% 1. get the connected region to the maximum in the normalized wavelet
% transform that exceeds coef
%   1a. find the 'watersheds' of the negative normalized wavelet transform,
%   and remove any of the watersheds within the connected region found in 1
%   which lie such that the maximum in the normalized wavelet transform is
%   completely outside the period range of that watershed
% 2. extend that region along the time axis until a local minimum in the
% wavelet transform is reached; this is the region for the event
% 3. invert the wavelet transform in the event region and remove the event
% from the signal
% 4. recalculate the wavelet transform and normalized wavelet transform
% 5. return to (1), unless the maximum value in the normalized wavelet
% transform does not exceed 2 or the number of iterations has reached
% max_events
event_centers = [];
min_period = duration();
center_period = duration();
max_period = duration();
minperiod_peakregion = duration();
maxperiod_peakregion = duration();
event_start = datetime();
event_center = datetime();
event_end = datetime();
omega_prime = zeros(size(wt_norm));
omega = zeros(size(omega_prime));
omega_pixels = {};
reconstructed_series = {};
cur_event = 1;

% if there isn't enough data, calc_wavelettransform will return empty
% arrays, so need to include code to account for that case
if isempty(wt_max); wt_max = 0; end
while wt_max > 2 && (cur_event < max_events || max_events == 0)
    % step 1
    [omega_curevent, idx] = bwselect(wt_norm > 1, col_max, row_max);
    [~, peakregion_idx] = bwselect(wt_norm > 2, col_max, row_max);
    [peakregion_rows, ~] = ind2sub(size(wt_norm), peakregion_idx);
    
    % to for sure avoid infinite loops, make sure we havent already
    % identified an event at the current center index
    if omega(row_max,col_max) == 1
        % set all pixels of the preliminary event region to 0
        wavelet_transform(idx) = 0;
        wt_norm(idx) = 0;
        [wt_max, ind_max] = max(wt_norm, [], 'all', 'linear');
        [row_max, col_max] = ind2sub(size(wt_norm), ind_max);
        continue
    end
        
        
    event_centers = [event_centers, ind_max];

    % step 1a
    negative_wtnorm = -wt_norm; negative_wtnorm(isnan(negative_wtnorm)) = 0;
    negative_wtnorm = imhmax(negative_wtnorm, 1.5);
    wtnorm_watersheds = watershed(negative_wtnorm);
    wtnorm_watersheds(wt_norm < 1) = 0;
    watersheds_abovethresh = unique(wtnorm_watersheds(wtnorm_watersheds ~= 0))';
    % get the watershed which contains the event center and its period
    % range
    central_watershed = wtnorm_watersheds(ind_max);
%     % also force the connected region to the event center where the
%     % normalized wavelet transform exceeds 2*coef to be part of the same
%     % watershed
%     conn_exceeding2 = bwselect(wt_norm > 2, col_max, row_max);
%     wtnorm_watersheds(conn_exceeding2) = central_watershed;
    [centralwatershed_rows,~] = find(wtnorm_watersheds == central_watershed);
    % include a 1 pixel buffer
    centralwatershed_rows = [centralwatershed_rows; ...
        min(centralwatershed_rows)-1; max(centralwatershed_rows)+1];
    for cur_watershed = watersheds_abovethresh
        [watershed_rows,~] = find(wtnorm_watersheds == cur_watershed);
        % find which watersheds intersect with the initial event region
        intersecting_idx = intersect(find(wtnorm_watersheds == cur_watershed), idx);
        % if that watershed's period range has no overlap with the central 
        % watershed's period range, remove its points from the event region
        if ~isempty(intersecting_idx) && ~any(ismember(watershed_rows, centralwatershed_rows))
            omega_curevent(intersecting_idx) = 0;
            idx = setdiff(idx, intersecting_idx);
        end
    end
    omega_prime(idx) = 1;

    % step 2
    % extend the event region in each direction until the bounding box
    % runs into a local minimum in normalized wavelet amplitude
    localmins_row = islocalmin(abs(wavelet_transform), 2);
    localmins_col = islocalmin(abs(wavelet_transform), 1);
    for curcol = 1:size(wavelet_transform, 2)
        if abs(wavelet_transform(1, curcol)) < abs(wavelet_transform(2, curcol))
            localmins_col(1, curcol) = 1;
        end
        if abs(wavelet_transform(end, curcol)) < abs(wavelet_transform(end-1, curcol))
            localmins_col(end, curcol) = 1;
        end
    end
    for currow = 1:size(wavelet_transform, 1)
        if abs(wavelet_transform(currow, 1)) < abs(wavelet_transform(currow, 2))
            localmins_row(currow, 1) = 1;
        end
        if abs(wavelet_transform(currow, end)) < abs(wavelet_transform(currow, end-1))
            localmins_row(currow, end) = 1;
        end
    end
    localminsinds = find(localmins_row & localmins_col);
    [~, localmincols] = ind2sub(size(wt_norm), localminsinds);
    [omegaprime_rows, omegaprime_cols] = ind2sub(size(wt_norm), idx);
    startcol = localmincols(localmincols < min(omegaprime_cols));
    endcol = localmincols(localmincols > max(omegaprime_cols));
    if isempty(startcol);startcol=1;end
    if isempty(endcol);endcol=size(wavelet_transform,2);end
    omega(min(omegaprime_rows):max(omegaprime_rows), startcol(end):endcol(1)) = 1;
    omega_curevent(min(omegaprime_rows):max(omegaprime_rows), startcol(end):endcol(1)) = 1;
    [omega_rows, omega_cols] = meshgrid(min(omegaprime_rows):max(omegaprime_rows), startcol(end):endcol(1));
    omega_inds = sub2ind(size(omega), omega_rows(:), omega_cols(:));
    omega_inds = setdiff(omega_inds, naninds);
    omega_pixels{cur_event} = omega_inds;

    min_period(cur_event) = period_vec(min(omegaprime_rows));
    center_period(cur_event) = period_vec(row_max);
    max_period(cur_event) = period_vec(max(omegaprime_rows));
    event_start(cur_event) = datetime_vec(startcol(end));
    event_center(cur_event) = datetime_vec(col_max);
    event_end(cur_event) = datetime_vec(endcol(1));
    minperiod_peakregion(cur_event) = period_vec(min(peakregion_rows));
    maxperiod_peakregion(cur_event) = period_vec(max(peakregion_rows));

    % step 3
    wt_rec = wavelet_transform;
    wt_rec(~omega_curevent) = 0;
    reconstructed_series{cur_event} = icwt(wt_rec, 'amor');

    % step 4
    fullTrace = fullTrace - reconstructed_series{cur_event}';
    % rather than recalculate the whole wavelet transform, recalculate the
    % 24 hour period around the event (plus the length of the event) and
    % substitute it into the existing wavelet transform
    firstind_to_replace = max(startcol(end)-4320, 1);
    lastind_to_replace = min(endcol(1)+4320, length(fullTrace));
    [wt_replacement, tau_rep, coi_rep] = cwt(...
        fullTrace(firstind_to_replace:lastind_to_replace), 'amor', ...
        minutes(smoothParam/60), 'PeriodLimits', [min(period_vec) max(period_vec)+minutes(1)]);
    [coimesh_rep, taumesh_rep] = meshgrid(coi_rep, tau_rep);
    wt_replacement(taumesh_rep >= coimesh_rep) = NaN;
    wt_beforereplacement = wavelet_transform;
    wavelet_transform(1:size(wt_replacement,1), ...
        firstind_to_replace:lastind_to_replace) = wt_replacement;
    % avoid putting in any NaNs that weren't there before
    wavelet_transform(isnan(wavelet_transform)) = wt_beforereplacement(isnan(wavelet_transform));

    % step 5
    wt_norm = abs(wavelet_transform) ./ threshold_vec(1:size(wavelet_transform,1));
    [wt_max, ind_max] = max(wt_norm, [], 'all', 'linear');
    [row_max, col_max] = ind2sub(size(wt_norm), ind_max);

    cur_event = cur_event + 1;
end

% create a table of basic properties of each event
varnames = {'MinPeriod', 'CenterPeriod', 'MaxPeriod', ...
    'StartTime', 'CenterTime', 'EndTime', ...
    'MinPeriod_PeakRegion', 'MaxPeriod_PeakRegion'};

if isempty(omega_pixels)
    events_table = table([], [], [], [], [], [], [], [], 'VariableNames', varnames);
else
    events_table = table(min_period', center_period', max_period', ...
        event_start', event_center', event_end', minperiod_peakregion', ...
        maxperiod_peakregion', 'VariableNames', varnames);
end

% get boundaries of each event region
event_bounds = cell(height(events_table), 1);
for cur_region = 1:length(event_bounds)
    dummyarray = zeros(size(omega));
    dummyarray(omega_pixels{cur_region}) = 1;
    curbounds = bwboundaries(dummyarray);
    event_bounds{cur_region} = curbounds{1};
end

end

