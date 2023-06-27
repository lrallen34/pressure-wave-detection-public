function plotevents_fromcatalog_simplified(eventcatalog, struct_of_tables, ...
    wt_struct, threshloc, outdir_main, minsensors)
%PLOTEVENTS_FROMCATALOG_simplified Creates plots of the wave events in a
%catalog and outputs .png images of the plots to a specified directory
%
%   General form: plotevents_fromcatalog_simplified(eventcatalog, struct_of_tables, ...
%    wt_struct, threshloc, outdir_main, minsensors)
%
%   INPUTS:
%       -eventcatalog (struct): catalog of wave events produced by
%       create_eventcatalog.m
%       -struct_of_tables (struct): structure array with tables of wave
%       events which was used to create the event catalog
%       -wt_struct (struct): structure array of the wavelet transforms of
%       pressure traces for each sensor in the network being processed
%       -threshloc (string): path to the text file containing mean wavelet
%       power values by wave period, which are then used to normalize the
%       wavelet transform (produced by calcthresh.m)
%       -outdir_main (string): path to the image output directory
%       -minsensors (numeric scalar): the minimum number of sensors that
%       must have captured an event for it to be plotted. Events which were
%       captured by fewer sensors will be skipped. Set to 1 to plot every
%       event.
%
%   OUTPUTS:
%       none (.png images will be output to subdirectories within
%       outdir_main). Three plots will be produced for each event: one with
%       extracted event and full pressure traces, one with the wavelet
%       power as a function of time and wave period, and one with the
%       normalized wavelet power as a function of time and wave period.
%
%   Written by Luke Allen (lrallen3@ncsu.edu), Dec 2022

% determine the 'central' sensor in each event to decide the order of
% sensors in each plot
[arrayofevents, ~, bestcenter] = determine_centralsensor(struct_of_tables);
fullnetwork = {struct_of_tables.sensor_id};

% read the mean wavelet power file
meanmat = readmatrix(threshloc);
period_vec = meanmat(:, 1);
meanvec = meanmat(:, 2);

% loop through each event in the catalog
for curevent = 1:length(eventcatalog)
    curID = eventcatalog(curevent).event_id;
    % if not enough sensors captured the event or it wasn't possible to                                                                                   -
    % calculate a wave speed (this can happen for events at the                                                                                    -
    % beginning/end of the processing period), skip it 
    if length(arrayofevents{curID}) >= minsensors && eventcatalog(curevent).wave_speed ~= 0
        disp(['plotting event ' num2str(curevent) '/' num2str(length(eventcatalog))])
        primary_sensor = fullnetwork{bestcenter(curevent)};
        all_sensors = fullnetwork(arrayofevents{curID});

        % shift the sensor array so the primary sensor is first
        primsensor_position = find(strcmp(all_sensors, primary_sensor));
        all_sensors = circshift(all_sensors, -primsensor_position+1);

        % find where each event is located in individual sensor event
        % tables and pull the needed data to plot
        dt_vec = cell(length(all_sensors), 1);
        event_series = cell(length(all_sensors), 1);
        full_trace = cell(length(all_sensors), 1);
        event_bounds = cell(length(all_sensors), 1);
        wavelet_transform = cell(length(all_sensors), 1);
        wt_norm = cell(length(all_sensors), 1);
        event_start = NaT(length(all_sensors), 1);
        event_center = NaT(length(all_sensors), 1);
        event_end = NaT(length(all_sensors), 1);
        event_maxperiod = NaN(length(all_sensors), 1);
        event_delay = NaN(length(all_sensors), 1);
        event_correlation = NaN(length(all_sensors), 1);
        for cursensor = 1:length(all_sensors)
            index_cursensor = find(strcmp({struct_of_tables.sensor_id}, all_sensors{cursensor}));
            eventsincursensor = struct_of_tables(index_cursensor).eventIDs;

            dt_vec{cursensor} = struct_of_tables(index_cursensor) ...
                .datetime_vecs;
            event_series{cursensor} = struct_of_tables(index_cursensor) ...
                .reconstructed_series{eventsincursensor == curID};
            full_trace{cursensor} = struct_of_tables(index_cursensor) ...
                .fullTrace;
            event_bounds{cursensor} = struct_of_tables(index_cursensor) ...
                .event_bounds{eventsincursensor == curID};
            wavelet_transform{cursensor} = abs(wt_struct(index_cursensor).values);
            wt_norm{cursensor} = wavelet_transform{cursensor} ./ ...
                meanvec(1:size(wavelet_transform{cursensor}, 1));
            event_start(cursensor) = struct_of_tables(index_cursensor)...
                .events_tables.StartTime(eventsincursensor == curID);
            event_center(cursensor) = struct_of_tables(index_cursensor)...
                .events_tables.CenterTime(eventsincursensor == curID);
            event_end(cursensor) = struct_of_tables(index_cursensor)...
                .events_tables.EndTime(eventsincursensor == curID);
            event_maxperiod(cursensor) = minutes(struct_of_tables(index_cursensor)...
                .events_tables.MaxPeriod(eventsincursensor == curID));
            event_delay(cursensor) = struct_of_tables(index_cursensor)...
                .delaytable.(['delayto_' primary_sensor])(eventsincursensor == curID);
            event_correlation(cursensor) = struct_of_tables(index_cursensor)...
                .corrtable.(['correlationwith_' primary_sensor])(eventsincursensor == curID);
        end

        % Plotting
        % define the start and end times of the plots
        plotstart = min(event_start) - minutes(30);
        plotend = max(event_end) + minutes(30);
        period_plotlims = [1 120];

        outdir_cur = [outdir_main '/' primary_sensor '_' ...
            datestr(event_center(1), 'yyyymmdd')];
        makeneededdirs(outdir_cur)

        tracefig = figure('visible', 'off');
        wtfig = figure('visible', 'off');
        wtnormfig = figure('visible', 'off');
        shedfig = figure('visible', 'off');

        % get the number of panels for the figure and resize accordingly
        num_panels = length(all_sensors);
        wtfig.Position = wtfig.Position .* [1 1 2 1.3*num_panels];
        wtnormfig.Position = wtnormfig.Position .* [1 1 4 1.3*num_panels];
        shedfig.Position = shedfig.Position .* [1 1 4 1.3*num_panels];
        if num_panels <= 2
            num_cols = 1;
        elseif num_panels <= 6
            num_cols = 2;
        else
            num_cols = 3;
        end
        set(tracefig, 'defaultAxesColorOrder', [0 0 0; 0 0.447 0.741])
        tracefig.Position = tracefig.Position .* [1 1 4 1.3*ceil(num_panels/num_cols)];

        ii = 1;
        for plotsensor = 1:length(all_sensors)

            set(0, 'CurrentFigure', tracefig);
            curax =  subplot(ceil(num_panels/num_cols), num_cols, ii);
            sensor_id = all_sensors{plotsensor};

            trace_to_plot = event_series{plotsensor};
            fulltrace_to_plot = full_trace{plotsensor};
            sensor_times = dt_vec{plotsensor};

            sensor_wt = wavelet_transform{plotsensor};
            sensor_wtnorm = wt_norm{plotsensor};
            sensor_bounds = event_bounds{plotsensor};

            if plotstart < sensor_times(1)
                padneeded = length(plotstart:seconds(10):sensor_times(1)-seconds(10));
                sensor_times = [(plotstart:seconds(10):sensor_times(1)-seconds(10))'; sensor_times];
                trace_to_plot = [zeros(length(sensor_times) - length(trace_to_plot), 1)', trace_to_plot];
                fulltrace_to_plot = [nan(length(sensor_times) - length(fulltrace_to_plot), 1); fulltrace_to_plot];
                sensor_wt = [nan(size(sensor_wtnorm,1), ...
                    length(sensor_times) - size(sensor_wt, 2)), sensor_wt];
                sensor_wtnorm = [nan(size(sensor_wtnorm,1), ...
                    length(sensor_times) - size(sensor_wtnorm, 2)), sensor_wtnorm];
                sensor_bounds(:, 2) = sensor_bounds(:, 2) + padneeded;
            end
            if plotend > sensor_times(end)
                sensor_times = [sensor_times; (sensor_times(end)+seconds(10):seconds(10):plotend)'];
                trace_to_plot = [trace_to_plot, zeros(length(sensor_times) - length(trace_to_plot), 1)'];
                fulltrace_to_plot = [fulltrace_to_plot; nan(length(sensor_times) - length(fulltrace_to_plot), 1)];
                sensor_wt = [sensor_wt, nan(size(sensor_wtnorm,1), ...
                    length(sensor_times) - size(sensor_wt, 2))];
                sensor_wtnorm = [sensor_wtnorm, nan(size(sensor_wtnorm,1), ...
                    length(sensor_times) - size(sensor_wtnorm, 2))];
            end
            
            % trim the arrays down
            ind_1 = find(sensor_times == plotstart - minutes(10));
            ind_2 = find(sensor_times == plotend + minutes(10));
            if isempty(ind_1) || isempty(ind_2); continue; end
            sensor_bounds(:, 2) = sensor_bounds(:, 2) - ind_1 + 1;
            sensor_times = sensor_times(ind_1:ind_2);
            trace_to_plot = trace_to_plot(ind_1:ind_2);
            fulltrace_to_plot = fulltrace_to_plot(ind_1:ind_2);
            sensor_wt = sensor_wt(:, ind_1:ind_2);
            sensor_wtnorm = sensor_wtnorm(:, ind_1:ind_2);
            if isempty(sensor_times); continue; end
            
            % calculate watershed transform of normalized wavelet power
            negative_wtnorm = -sensor_wtnorm;
            negative_wtnorm(isnan(negative_wtnorm)) = 0;
            negative_wtnorm = imhmax(negative_wtnorm, 1.5);
            wtnorm_watersheds = watershed(negative_wtnorm);

            plot(sensor_times, trace_to_plot, 'LineWidth', 2)
            xlim([plotstart plotend])
            grid on
            ylabel('Extracted event (hPa)')
            t = title([sensor_id ' (R = ' ...
                num2str(event_correlation(plotsensor)) ', lag = ' ...
                char(duration(0, 0, event_delay(plotsensor))) ')']);
            t.Interpreter = 'none';
            hold on
            yyaxis right
            ylabel('Total pressure (hPa)')
            plot(sensor_times, fulltrace_to_plot, 'LineWidth', 2)
            curax.FontSize = 20;
            

%             curax = subplot(num_panels, 2, 2*ii);
%             plot(sensor_times, fulltrace_to_plot, 'LineWidth', 2)
%             xlim([plotstart plotend])
%             grid on
%             ylabel('hPa')
%             t = title(['Full pressure trace from ' sensor_id]);
%             t.Interpreter = 'none';
%             curax.FontSize = 16;

            set(0, 'CurrentFigure', wtfig);

            curax = subplot(num_panels, 1, ii);
            p = pcolor(1:length(sensor_times), period_vec(1:size(sensor_wt,1)), abs(sensor_wt));
            p.EdgeColor = 'none';
            title(['Wavelet energy and event bounds for ' sensor_id])
            c = colorbar;
            c.Label.String = 'Wavelet Amplitude (hPa)';
            caxis([0 0.75])
            ylabel('Wave Period (min.)')
            colormap(cmocean('rain'))
            ylim(period_plotlims)
            xlim([find(sensor_times == plotstart) find(sensor_times == plotend)])
            set(gca, 'YScale', 'log')
            yticks(unique(round(logspace(log10(period_plotlims(1)), log10(period_plotlims(2)), 5))));
            hold on
            contour(1:length(sensor_times), ...
                period_vec(1:size(sensor_wt,1)), sensor_wtnorm, [5 5], ...
                '--', 'LineWidth', 1.5)
            contour(1:length(sensor_times), ...
                period_vec(1:size(sensor_wt,1)), sensor_wtnorm, [10 10], ...
                'LineWidth', 1.5)
            curx = sensor_bounds(:, 2);
            cury = period_vec(sensor_bounds(:, 1));
            plot(curx, cury, 'm', 'LineWidth', 3)
            curax.FontSize = 16;
            % manually define x tick labels to make this work with older matlab
            % versions
            xt = xticks;
            labelstrings = datestr(sensor_times(xt), 'HH:MM:SS');
            xticklabels(labelstrings)
            hold off

            set(0, 'CurrentFigure', wtnormfig)
            ax1 = subplot(num_panels, 1, ii);
            p = pcolor(ax1, 1:length(sensor_times), period_vec(1:size(sensor_wt,1)), sensor_wtnorm);
            p.EdgeColor = 'none';
            title(['Normalized wavelet energy and event bounds for ' sensor_id])
            c = colorbar;
            c.Label.String = 'Normalized Wavelet Amplitude';
            caxis([0 15])
            ylabel('Wave Period (min.)')
            colormap(cmocean('rain'))
            ylim(period_plotlims)
            xlim([find(sensor_times == plotstart) find(sensor_times == plotend)])
            set(gca, 'YScale', 'log')
            yticks(unique(round(logspace(log10(period_plotlims(1)), log10(period_plotlims(2)), 5))));
            hold on
            contour(1:length(sensor_times), ...
                period_vec(1:size(sensor_wt,1)), sensor_wtnorm, [5 5], ...
                '--', 'LineWidth', 1.5)
            contour(1:length(sensor_times), ...
                period_vec(1:size(sensor_wt,1)), sensor_wtnorm, [10 10], ...
                'LineWidth', 1.5)
            curx = sensor_bounds(:, 2);
            cury = period_vec(sensor_bounds(:, 1));
            plot(curx, cury, 'm', 'LineWidth', 3)
%             for curshed = watershed_values'
%                 [shedy, shedx] = find(wtnorm_watersheds == curshed, 1);
%                 contour = bwtraceboundary(wtnorm_watersheds == curshed, [shedy, shedx], 'W');
%                 plot(contour(:,2), contour(:,1), 'Color', [0.8500 0.3250 0.0980])
%             end
            wshedmask = boundarymask(wtnorm_watersheds);
            ax2 = subplot(num_panels, 1, ii);
            hold on
            surf(ax2, 1:length(sensor_times), ...
                period_vec(1:size(sensor_wt,1)), ...
                zeros(size(wtnorm_watersheds)), ...
                'EdgeColor', 'none', 'AlphaData', wshedmask, ...
                'FaceAlpha', 'flat', 'FaceColor', [0.8500 0.3250 0.0980]);
            ax1.FontSize = 16;
            % manually define x tick labels to make this work with older matlab
            % versions
            xt = xticks;
            labelstrings = datestr(sensor_times(xt), 'HH:MM:SS');
            xticklabels(labelstrings)
            hold off

            ii = ii+1;

        end

        set(0, 'CurrentFigure', tracefig);
        st = sgtitle(['Events corresponding to event at ' ...
            char(event_center(1)) ' in ' primary_sensor]);
        st.Interpreter = 'none';
        print('-r300', tracefig, [outdir_cur '/' ...
            datestr(event_center(1), 'yyyymmdd_HHMMSS')], ...
            '-dpng')

        close(tracefig)

        set(0, 'CurrentFigure', wtfig);
        st = sgtitle(['Events corresponding to event at ' ...
            char(event_center(1)) ' in ' primary_sensor]);
        st.Interpreter = 'none';
        print('-r300', wtfig, [outdir_cur '/' ...
            datestr(event_center(1), 'yyyymmdd_HHMMSS') '_wavelettransforms'], ...
            '-dpng')

        close(wtfig)

        set(0, 'CurrentFigure', wtnormfig);
        st = sgtitle(['Events corresponding to event at ' ...
            char(event_center(1)) ' in ' primary_sensor]);
        st.Interpreter = 'none';
        print('-r300', wtnormfig, [outdir_cur '/' ...
            datestr(event_center(1), 'yyyymmdd_HHMMSS') '_norm_wts'], ...
            '-dpng')

        close(wtnormfig)
    end
end

end
