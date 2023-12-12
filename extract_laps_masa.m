function lap_times = extract_laps_masa(varargin)

if isempty(varargin{1})
    plot_option= 0; % do not plot
else
    plot_option= 1;
end

if ~isempty(varargin(2))
    Behaviour = varargin{2};
end

if ~isempty(varargin(3))
    position = varargin{3};
end

% load extracted_position
parameters=list_of_parameters;

if plot_option
    figure;
end

if ~isempty(Behaviour.start_time_all)

    for track_id=1:length(position.linear)
        start_indices= [];

        x = position.linear(track_id).linear;  % position data
        t= position.t;
        
        this_track_indices = find(Behaviour.track_ID_all == track_id);

        for nlap = 1:length(this_track_indices)
            start_indices(nlap) = find(t == Behaviour.start_time_all(this_track_indices(nlap)));
        end

        lap_times(track_id).completeLaps_start = [];
        lap_times(track_id).completeLaps_end = [];
        lap_times(track_id).abortedLaps_start = [];
        lap_times(track_id).abortedLaps_end = [];
        lap_times(track_id).completeLaps_id = [];
        lap_times(track_id).abortedLaps_id = [];

%         lap_count = 1;
        for nlap = 1:length(start_indices)

            if nlap < length(start_indices)
                current_lap_x = x(start_indices(nlap):start_indices(nlap+1));
                current_lap_t = t(start_indices(nlap):start_indices(nlap+1));
            else
                current_lap_x = x(start_indices(nlap):end);
                current_lap_t = t(start_indices(nlap):end);
            end


            if length(current_lap_x) > 1 && length(current_lap_x(~isnan(current_lap_x))) >1% Only if getting more than 1 datapoint (maybe noise)
                on_track_x = current_lap_x(~isnan(current_lap_x));
                on_track_t = current_lap_t(~isnan(current_lap_x));
                
                if sum(on_track_x==0)>0
                    start_position = find(on_track_x == 0);
                    on_track_x = on_track_x(start_position(1):end);
                    on_track_t = on_track_t(start_position(1):end);
                end

                if on_track_x(end) ~= on_track_x(end-1)
                   on_track_x(end) = on_track_x(end-1);
                   on_track_t(end) = on_track_t(end-1);
                end

                [last_position last_position_index] = max(on_track_x);
                
                lap_times(track_id).start(nlap) = t(start_indices(nlap));
                lap_times(track_id).end(nlap) = on_track_t(end);
                lap_times(track_id).lap_id(nlap) = nlap;

                if last_position >= 139 % sometimes last lap ends before 140cm

                    lap_times(track_id).completeLaps_start = [lap_times(track_id).completeLaps_start t(start_indices(nlap))];
                    lap_times(track_id).completeLaps_end = [lap_times(track_id).completeLaps_end lap_times(track_id).end(nlap)];
                    lap_times(track_id).completeLaps_id = [lap_times(track_id).completeLaps_id nlap];
                else

                    lap_times(track_id).abortedLaps_start = [lap_times(track_id).abortedLaps_start t(start_indices(nlap))];
                    lap_times(track_id).abortedLaps_end = [lap_times(track_id).abortedLaps_end lap_times(track_id).end(nlap)];
                    lap_times(track_id).abortedLaps_id = [lap_times(track_id).abortedLaps_id nlap];
                end


%                 lap_count = lap_count + 1;
            end


        end

        %     if length(lap_times(track_id).start) == length(lap_times(track_id).end)
        %         lap_times(track_id).duration = lap_times(track_id).end  - lap_times(track_id).start; % lap durations
        %     else % probably missing one end, treat the last position as the end (very rare this will happen)
        %         lap_times(track_id).end = [lap_times(track_id).end t(end)];
        %         lap_times(track_id).duration = lap_times(track_id).end  - lap_times(track_id).start; % lap durations
        %     end

        for i=1:length(lap_times(track_id).start)
            % store x and t for lap (includes end zone)
            lap_times(track_id).lap(i).x= x(t>=lap_times(track_id).start(i) & t<=lap_times(track_id).end(i));
            lap_times(track_id).lap(i).t= t(t>=lap_times(track_id).start(i) & t<=lap_times(track_id).end(i));
            lap_times(track_id).lap(i).v_cm= position.v_cm(t>=lap_times(track_id).start(i) & t<=lap_times(track_id).end(i));

            %         % store x and t only for end zone
            %         lap_times(track_id).end_zone(i).x= lap_times(track_id).lap(i).x(lap_times(track_id).lap(i).x < (min_pos+delta) | lap_times(track_id).lap(i).x > (max_pos-delta));
            %         lap_times(track_id).end_zone(i).t= lap_times(track_id).lap(i).t(lap_times(track_id).lap(i).x < (min_pos+delta) | lap_times(track_id).lap(i).x > (max_pos-delta));
            %         lap_times(track_id).end_zone(i).v_cm=  lap_times(track_id).lap(i).v_cm(lap_times(track_id).lap(i).x < (min_pos+delta) | lap_times(track_id).lap(i).x > (max_pos-delta));
        end


        lap_times(track_id).total_number_of_laps = length(lap_times(track_id).start);

        if plot_option
            subplot(1,length(position.linear),track_id);
            plot([lap_times(track_id).lap.t],[lap_times(track_id).lap.x],'k.');
            hold on;
            %         scatter(t(start_indices),x(start_indices),'r','filled')
            %         scatter(t(start_indices),x(start_indices),'b','filled')
            scatter(lap_times(track_id).start,zeros(1,length(lap_times(track_id).start)),'g','filled')
            %         hold on
            %         scatter(t(start_indices(30)),x(start_indices(30)),'b','filled')
            %         scatter(t(start_indices(40)),x(start_indices(40)),'b','filled')
            %         scatter(t(start_indices(60)),x(start_indices(60)),'b','filled')
            %         scatter(t(start_indices(90)),x(start_indices(90)),'b','filled')
            ylim([-10 150])
            %         ylim([-0.5 10])
        end
    end

else
    for track_id=1:length(position.linear)

        x = position.linear(track_id).linear;  % position data
        t= position.t;

        max_pos = max(position.linear(track_id).linear);   % find max and min position
        min_pos = min(position.linear(track_id).linear);
        delta   = 0.5;  % segment track length


        % track_edge_indices = find( x<(min_pos+delta) | x>(max_pos-delta)); % find indices when animal at top or bottom of track
        start_edge_indices = find(x<=(min_pos+delta)); % find indices when animal at start of the track
        %     end_edge_indices = find(x>(max_pos-delta)); % find indices when animnal at end of the track

        start_indices = [start_edge_indices(1) start_edge_indices(find(abs([0 diff((start_edge_indices))])>200)) start_edge_indices(end)]; % looks for jumps larger than 200
        %     end_indices = [end_edge_indices(find(abs(diff((end_edge_indices)))>100)) end_edge_indices(end)];

        %     lap_times(track_id).start = t(start_indices); % lap start time
        %     lap_times(track_id).end  = t(end_indices);   % lap end time

        for nlap = 1:length(start_indices)
            %         if x(start_indices(nlap)) > x(start_indices(nlap)-1)
            %             start_indices(nlap) = [];
            %         end
            % Remove if the position changes 100 frames before the putative
            % start position.
            if sum(diff(x(start_indices(nlap)-101:start_indices(nlap)-1))) > 0
                start_indices(nlap) = [];
            end
        end


        lap_times(track_id).completeLaps_start = [];
        lap_times(track_id).completeLaps_end = [];
        lap_times(track_id).abortedLaps_start = [];
        lap_times(track_id).abortedLaps_end = [];
        lap_times(track_id).completeLaps_id = [];
        lap_times(track_id).abortedLaps_id = [];

        lap_count = 1;
        for nlap = 1:length(start_indices)

            if nlap < length(start_indices)
                current_lap_x = x(start_indices(nlap):start_indices(nlap+1)-1);
                current_lap_t = t(start_indices(nlap):start_indices(nlap+1)-1);
            else
                current_lap_x = x(start_indices(nlap):end);
                current_lap_t = t(start_indices(nlap):end);
            end

            if x(start_indices(nlap)) > x(start_indices(nlap)-1)
                continue
            end
            if length(current_lap_x) > 1 % Only if getting more than 1 datapoint (maybe noise)
                on_track_x = current_lap_x(~isnan(current_lap_x));
                [last_position last_position_index] = max(on_track_x);
                on_track_t = current_lap_t(~isnan(current_lap_x));
                lap_times(track_id).start(lap_count) = t(start_indices(nlap));
                lap_times(track_id).end(lap_count) = on_track_t(end);
                lap_times(track_id).lap_id(lap_count) = lap_count;

                if last_position >= 138.5 % sometimes last lap ends before 140cm

                    lap_times(track_id).completeLaps_start = [lap_times(track_id).completeLaps_start t(start_indices(lap_count))];
                    lap_times(track_id).completeLaps_end = [lap_times(track_id).completeLaps_end lap_times(track_id).end(lap_count)];
                    lap_times(track_id).completeLaps_id = [lap_times(track_id).completeLaps_id lap_count];
                else

                    lap_times(track_id).abortedLaps_start = [lap_times(track_id).abortedLaps_start t(start_indices(lap_count))];
                    lap_times(track_id).abortedLaps_end = [lap_times(track_id).abortedLaps_end lap_times(track_id).end(lap_count)];
                    lap_times(track_id).abortedLaps_id = [lap_times(track_id).abortedLaps_id lap_count];
                end


                lap_count = lap_count + 1;
            end


        end

        %     if length(lap_times(track_id).start) == length(lap_times(track_id).end)
        %         lap_times(track_id).duration = lap_times(track_id).end  - lap_times(track_id).start; % lap durations
        %     else % probably missing one end, treat the last position as the end (very rare this will happen)
        %         lap_times(track_id).end = [lap_times(track_id).end t(end)];
        %         lap_times(track_id).duration = lap_times(track_id).end  - lap_times(track_id).start; % lap durations
        %     end

        for i=1:length(lap_times(track_id).start)
            % store x and t for lap (includes end zone)
            lap_times(track_id).lap(i).x= x(t>=lap_times(track_id).start(i) & t<=lap_times(track_id).end(i));
            lap_times(track_id).lap(i).t= t(t>=lap_times(track_id).start(i) & t<=lap_times(track_id).end(i));
            lap_times(track_id).lap(i).v_cm= position.v_cm(t>=lap_times(track_id).start(i) & t<=lap_times(track_id).end(i));

            %         % store x and t only for end zone
            %         lap_times(track_id).end_zone(i).x= lap_times(track_id).lap(i).x(lap_times(track_id).lap(i).x < (min_pos+delta) | lap_times(track_id).lap(i).x > (max_pos-delta));
            %         lap_times(track_id).end_zone(i).t= lap_times(track_id).lap(i).t(lap_times(track_id).lap(i).x < (min_pos+delta) | lap_times(track_id).lap(i).x > (max_pos-delta));
            %         lap_times(track_id).end_zone(i).v_cm=  lap_times(track_id).lap(i).v_cm(lap_times(track_id).lap(i).x < (min_pos+delta) | lap_times(track_id).lap(i).x > (max_pos-delta));
        end


        lap_times(track_id).total_number_of_laps = length(lap_times(track_id).start);

        if plot_option
            subplot(1,length(position.linear),track_id);
            plot([lap_times(track_id).lap.t],[lap_times(track_id).lap.x],'k.');
            hold on;
            %         scatter(t(start_indices),x(start_indices),'r','filled')
            %         scatter(t(start_indices),x(start_indices),'b','filled')
            scatter(lap_times(track_id).start,zeros(1,length(lap_times(track_id).start)),'g','filled')
            %         hold on
            %         scatter(t(start_indices(30)),x(start_indices(30)),'b','filled')
            %         scatter(t(start_indices(40)),x(start_indices(40)),'b','filled')
            %         scatter(t(start_indices(60)),x(start_indices(60)),'b','filled')
            %         scatter(t(start_indices(90)),x(start_indices(90)),'b','filled')
            ylim([-10 150])
            %         ylim([-0.5 10])
        end
    end
end


if plot_option
    for track_id = 1:2
        if ~isempty(lap_times(track_id).abortedLaps_id)
            figure
            for nlap = 1:length(lap_times(track_id).abortedLaps_id)
                subplot(ceil(length(lap_times(track_id).abortedLaps_id)/5),5,nlap)
                plot(lap_times(track_id).lap(lap_times(track_id).abortedLaps_id(nlap)).t,lap_times(track_id).lap(lap_times(track_id).abortedLaps_id(nlap)).x)
                title(lap_times(track_id).abortedLaps_id(nlap))
            end
        end

        sgtitle(sprintf('Aborted laps track %i',track_id))
    end

end


% save extracted_laps lap_times

end