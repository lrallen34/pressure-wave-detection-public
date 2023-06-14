function newStruct = assign_eventIDs(struct_of_tables)
%ASSIGN_EVENTIDS Determines whether events from a structure output by
%determinecoherence_fullnetwork.m are coherent in both directions (i.e., if
%swapping 'primary' and 'alternate' sensors). Then assigns event IDs by the
%following process:
% 1. Iterate through each sensor. On the first sensor, create an ID for
% each event.
% 2. In each subsequent sensor, if an event is coherent with a previous
% sensor, check whether that event is coherent when using the previous
% sensor as the primary sensor and the current sensor as the alternate
% sensor. If so, use that event ID. In the event of conflicting event IDs,
% use the one which is associated with the highest cross-correlation with
% the current sensor.
%
%   INPUTS:
%
%       struct_of_tables (struct): structure array output by
%       determinecoherence_fullnetwork.m which contains fields for the
%       sensor ID, table of events in each sensor along with tables
%       indicating the corresponding event center times, start times, end
%       times, and correlations with the other sensors
%
%   OUTPUTS:
%
%       newStruct (struct): structure array which is the same as
%       struct_of_tables but with two new fields added:
%       1. match_table: table indicating for each event whether or not the
%       event was reciprocally identified in the other sensor
%       2. eventIDs: vectors of event IDs with one element for each event
%       in the sensor's tables
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Mar 2022
%

newStruct = struct_of_tables;
alt_sensors = {struct_of_tables.sensor_id};
next_id = 1;
for sensor_n = 1:length(newStruct)
    cur_id = newStruct(sensor_n).sensor_id;
    
    if isempty(newStruct(sensor_n).centertable); continue; end
    current_sensor_event_centers = newStruct(sensor_n).centertable.(cur_id);
    match_table = table('Size', size(newStruct(sensor_n).centertable), ...
        'VariableTypes', repmat("logical", size(alt_sensors)));
    match_table.Properties.VariableNames = alt_sensors;
    eventIDs = zeros(height(match_table), 1);

    % loop through each event, check each other sensor's events
    allIDvecs = cell(1, height(newStruct(sensor_n).table_allsensors));
    ID_histories = cell(1, height(newStruct(sensor_n).table_allsensors));
    for cur_event = 1:height(newStruct(sensor_n).centertable)
        event_to_look_for = current_sensor_event_centers(cur_event);
        if isnat(event_to_look_for); continue; end
        for alt_sensor = 1:length(newStruct)
            % this finds where the current matched event is located in the
            % alternate sensor's event table
            if isempty(newStruct(alt_sensor).centertable)
                continue
            end
            alt_event_loc = find(newStruct(alt_sensor).centertable.(alt_sensors{alt_sensor}) == ...
                newStruct(sensor_n).centertable.(alt_sensors{alt_sensor})(cur_event));
            % this checks whether, when switching the alternate and primary
            % sensors, the same event is found in the original primary
            % sensor
            check = newStruct(alt_sensor).centertable.(cur_id)(alt_event_loc) == ...
                event_to_look_for & strcmp(newStruct(alt_sensor)...
                .coherencetable.(['appearsin_' cur_id])(alt_event_loc), ...
                'coherent');
            if check; match_table.(alt_sensors{alt_sensor})(cur_event) = 1; end
                
        end
        % if any sensors prior to the current sensor had the event, use
        % their event ID; otherwise, create a new one
        if sensor_n > 1 && any(table2array(match_table(cur_event,1:sensor_n-1)))
            sensors_with_event = find(table2array(match_table(cur_event,1:sensor_n-1)));
            IDvec = zeros(length(sensors_with_event), 1);
            sensor_j = 1;
            for curmatch = sensors_with_event
                IDvec(sensor_j) = ...
                    newStruct(curmatch).eventIDs(...
                    newStruct(curmatch).centertable.(alt_sensors{curmatch}) == ...
                newStruct(sensor_n).centertable...
                .(alt_sensors{curmatch})(cur_event));
                sensor_j = sensor_j + 1;
            end
            allIDvecs{cur_event} = IDvec;

            % in case different event IDs come up, base the event ID to go
            % with this on whichever one is associated with the highest correlation
            corrs_forthisevent = table2array(...
                newStruct(sensor_n).corrtable(...
                cur_event, sensors_with_event(sensors_with_event < sensor_n)));
            eventIDs(cur_event) = IDvec(...
                corrs_forthisevent == max(corrs_forthisevent));
            ID_histories{cur_event} = eventIDs(cur_event);

            % If this is a repeat of a previous event ID, have to choose one
            % to assign to the ID. The highest overall correlation between
            % the two events will get the respective event ID it
            % corresponds to, then if the other event can be corresponded
            % to another event ID it will. If not, that event gets assigned
            % a new ID.
            while any(eventIDs(1:cur_event-1) == eventIDs(cur_event))
                priorappearance = find(eventIDs(1:cur_event-1) == eventIDs(cur_event));
                sensors_with_priorevent = find(table2array(match_table(priorappearance,1:sensor_n-1)));
                corrs_forpriorappearance = table2array(...
                newStruct(sensor_n).corrtable(...
                priorappearance, sensors_with_priorevent(sensors_with_priorevent < sensor_n)));
                % in case of a tie, the earlier event gets the ID
                if max(corrs_forthisevent) > max(corrs_forpriorappearance)
                    % in this case, go back and check for other possible
                    % IDs for the previous event
                    IDs_prior = allIDvecs{priorappearance};
                    possible_newassignments = IDs_prior(IDs_prior ~= eventIDs(cur_event));
                    possible_newassignmentcorrs = corrs_forpriorappearance(IDs_prior ~= eventIDs(cur_event));
                    
                    % check that we haven't already tried any of the
                    % possible reassignments when we originally processed
                    % this event
                    possibleassignments_alreadyused = ismember(possible_newassignments, ID_histories{priorappearance});
                    possible_newassignments(possibleassignments_alreadyused) = [];
                    possible_newassignmentcorrs(possibleassignments_alreadyused) = [];

                    if isempty(possible_newassignments)
                        eventIDs(priorappearance) = next_id;
                        ID_histories{priorappearance} = [ID_histories{priorappearance}, next_id];
                        next_id = next_id + 1;
                    else
                        eventIDs(priorappearance) = ...
                            possible_newassignments(...
                            possible_newassignmentcorrs == ...
                            max(possible_newassignmentcorrs));
                        ID_histories{priorappearance} = [ID_histories{priorappearance}, eventIDs(priorappearance)];
                    end

                else
                    % in this case, have to reassign the current event
                    possible_newassignments = IDvec(IDvec ~= eventIDs(priorappearance));
                    possible_newassignmentcorrs = corrs_forthisevent(IDvec ~= eventIDs(priorappearance));

                    if isempty(possible_newassignments)
                        eventIDs(cur_event) = next_id;
                        ID_histories{cur_event} = [ID_histories{cur_event}, next_id];
                        next_id = next_id + 1;
                    else
                        eventIDs(cur_event) = ...
                            possible_newassignments(...
                            possible_newassignmentcorrs == ...
                            max(possible_newassignmentcorrs));
                    end
                end

            end

        else 
            eventIDs(cur_event) = next_id;
            ID_histories{cur_event} = next_id;
            next_id = next_id + 1;
        end
    end
    
    % before recording the event IDs to the struct array, check that there
    % are no conflicts
    [uniqueIDs, IDinds] = unique(eventIDs);
    while length(uniqueIDs) < length(eventIDs)
        % if there is a conflict, resolve it
        duplicate_indices = setdiff( 1:numel(eventIDs), IDinds );
        repeated_ID = eventIDs(duplicate_indices(1));
        % rule out the event that has the smallest maximum
        % cross-correlation with other sensors for that event
        n = 1;
        conflicted_events = find(eventIDs == repeated_ID);
        maxcorrs = zeros(length(conflicted_events), 1);
        for conflicting_event = conflicted_events'
           sensors_matched = find(table2array(match_table(conflicting_event,1:sensor_n-1)));
           curevent_IDvec = allIDvecs{conflicting_event};
           sensors_matching_ID = sensors_matched(curevent_IDvec == repeated_ID);
           crosscorrs_curevent = table2array(...
           newStruct(sensor_n).corrtable(conflicting_event, sensors_matching_ID));
           maxcorrs(n) = max(crosscorrs_curevent);
           n = n+1;
        end
        
        event_to_reprocess = conflicted_events(maxcorrs == min(maxcorrs));
        
        sensors_with_priorevent = find(table2array(match_table(event_to_reprocess,1:sensor_n-1)));
        corrs_forpriorappearance = table2array(...
            newStruct(sensor_n).corrtable(...
            priorappearance, sensors_with_priorevent(sensors_with_priorevent < sensor_n)));
        
        possible_newassignments = allIDvecs{event_to_reprocess};
        possible_newassignmentcorrs = corrs_forpriorappearance;
        
        % check that we haven't already tried any of the
        % possible reassignments when we originally processed
        % this event
        possibleassignments_alreadyused = ismember(possible_newassignments, ID_histories{event_to_reprocess});
        possible_newassignments(possibleassignments_alreadyused) = [];
        possible_newassignmentcorrs(possibleassignments_alreadyused) = [];
        
        if isempty(possible_newassignments)
            eventIDs(event_to_reprocess) = next_id;
            ID_histories{event_to_reprocess} = [ID_histories{event_to_reprocess}, next_id];
            next_id = next_id + 1;
        else
            eventIDs(event_to_reprocess) = ...
                possible_newassignments(...
                possible_newassignmentcorrs == ...
                max(possible_newassignmentcorrs));
            ID_histories{priorappearance} = [ID_histories{event_to_reprocess}, next_id];
        end
        [uniqueIDs, IDinds] = unique(eventIDs);
    end
    

    newStruct(sensor_n).matchtable = match_table;
    newStruct(sensor_n).eventIDs = eventIDs;
end
end

