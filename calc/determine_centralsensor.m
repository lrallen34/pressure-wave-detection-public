function [arrayofevents, arrayofeventnames, bestcenter] = determine_centralsensor(struct_of_tables)
%DETERMINE_CENTRALSENSOR Determines the best 'central' or 'primary' sensor
%for each event ID in a structure of tables output by assign_eventIDs.m.
%The process of determining this involves finding the sensor with the
%highest sum of optimal correlations for that event in the other sensors
%sharing the event ID.
%   
%   INPUTS:
%
%       struct_of_tables(struct): structure array output by
%       assign_eventIDs.m which contains fields for the
%       sensor ID, table of events in each sensor along with tables
%       indicating the corresponding event center times, start times, end
%       times, correlations with the other sensors, whether a reciprocal
%       event match was found in the other sensors, and vectors of the
%       event IDs corresponding to each row of the tables.
%
%   OUTPUTS:
%
%       arrayofevents (cell of numeric vectors): vectors corresponding to
%       each unique event ID which list the indices of the sensors found to
%       have that event
%
%       arrayofeventnames (cell of strings): comma-separated lists
%       corresponding to each unique event ID which list the IDs of the
%       sensors found to have that event
%
%       bestcenter (numeric vector): vector where each element corresponds
%       to a unique event ID which contains the index of the sensor found
%       to be the best 'central' sensor for that event
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Mar 2022
%

% loop through each event ID, determine which sensors have that event ID
allIDs = {struct_of_tables.eventIDs};
concatenatedIDs = {cat(1, allIDs{:})};
uniqueIDs = unique(concatenatedIDs{1});
arrayofevents = cell(length(uniqueIDs), 1);
arrayofeventnames = cell(length(uniqueIDs), 1);
bestcenter = NaN(length(uniqueIDs), 1);
for curevent = 1:length(uniqueIDs)
    curevid = uniqueIDs(curevent);
    sensors_with_event = [];
    names_with_event = [];
    for cursensor = 1:length(struct_of_tables)
        if any(struct_of_tables(cursensor).eventIDs == curevid)
            sensors_with_event = [sensors_with_event, cursensor];
            if isempty(names_with_event)
                names_with_event = [names_with_event, ...
                struct_of_tables(cursensor).sensor_id];
            else
                names_with_event = [names_with_event, ', ' ...
                struct_of_tables(cursensor).sensor_id];
            end
        end
    end
    
    i = 1;
    % once the sensors with the event have been determined, get each
    % sensor's correlation with the other sensors for this event
    corrsum = NaN(length(sensors_with_event), 1);
    nummatches = zeros(length(sensors_with_event), 1);
    for curcenter = sensors_with_event
        eventloc_inthissensor = find(struct_of_tables(curcenter).eventIDs == curevid);
        eventcorrs_curcenter = table2array(struct_of_tables(curcenter).corrtable(eventloc_inthissensor, :));
        nummatches(i) = length(eventcorrs_curcenter(eventcorrs_curcenter > 0.65));
        corrsum(i) = sum(eventcorrs_curcenter, 'omitnan');
        i = i+1;
    end
    arrayofevents{curevid} = sensors_with_event;
    arrayofeventnames{curevid} = names_with_event;
    % the sensor which matches with the most other sensors will be the
    % central sensor
    mostmatches = find(nummatches == max(nummatches));
    % in the event of a tie...
    if length(mostmatches) > 1
        corrsum = corrsum(mostmatches);
        sensors_with_mostmatches = sensors_with_event(mostmatches);
        % the sensor with highest sum of correlations will be the central sensor
        findbestcenters = sensors_with_mostmatches(corrsum == max(corrsum, [], 'omitnan'));
        % in event of two sensors being equally good, take the first one
        bestcenter(curevid) = findbestcenters(1);
    else
        bestcenter(curevid) = sensors_with_event(mostmatches);
    end
end
end

