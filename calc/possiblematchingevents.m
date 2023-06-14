function potentialmatch = possiblematchingevents(primary_event_table, ...
    secondary_event_table, maxdelay)
%POSSIBLEMATCHINGEVENTS Produces an array of events to check for matches
%between some primary sensor event table and some secondary sensor event
%table
%
%   INPUTS:
%       
%       primary_event_table and secondary_event_table (tables): tables of 
%       wave events produced by identifyevents.m
%
%       maxdelay (scalar numeric): the maximum possible gap between events
%       in the two sensors to pair them together, in hours
%
%   OUTPUTS:
%
%       potentialmatch (cell): cell array with each cell representing a row
%       of the primary_event_table and each element within the cell
%       representing a row of the secondary_event_table. the cell elements
%       indicate events in the secondary sensor which can be checked for
%       coherence with the event corresponding to the cell in the primary
%       sensor.
%
%   Written by Luke Allen (lrallen3@ncsu.edu), last updated 2023/06/13
%
potentialmatch = cell(height(primary_event_table), 1);
if isempty(secondary_event_table)
    return
end
for cur_event = 1:height(primary_event_table)
    for candidate_event = 1:height(secondary_event_table)
        % CONDITIONS:
        % 1. secondary event center period within primary event period
        % range
        % 2. primary event center period within secondary event period
        % range
        % 3. the earlier event's end time should be no more than maxdelay hours
        % before the later event's start time
        meets_conditions = ...
            secondary_event_table.CenterPeriod(candidate_event) >= ...
            primary_event_table.MinPeriod(cur_event) & ...
            secondary_event_table.CenterPeriod(candidate_event) <= ...
            primary_event_table.MaxPeriod(cur_event) & ...
            primary_event_table.CenterPeriod(cur_event) >= ...
            secondary_event_table.MinPeriod(candidate_event) & ...
            primary_event_table.CenterPeriod(cur_event) <= ...
            secondary_event_table.MaxPeriod(candidate_event) & ...
            min([primary_event_table.EndTime(cur_event), ...
            secondary_event_table.EndTime(candidate_event)]) + hours(maxdelay) > ...
            max([primary_event_table.StartTime(cur_event), ...
            secondary_event_table.StartTime(candidate_event)]);

        if meets_conditions
            potentialmatch{cur_event} = ...
                [potentialmatch{cur_event}, candidate_event];
        end
    end
end

end

