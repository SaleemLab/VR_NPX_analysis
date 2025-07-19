function session_behaviour_summary = plot_session_behaviour_summary(behaviour,varargin)
% SETTINGS
p = inputParser;
addParameter(p, 'task_info', []);
addParameter(p, 'session_name', []);
% addParameter(p, 'output', []);
% addParameter(p, 'plot_option', 0);
% addParameter(p, 'subject_id', 0);
% addParameter(p, 'time_window', 0.1);

parse(p, varargin{:});

task_info = p.Results.task_info;
session_name = p.Results.session_name;



if ~isfield(behaviour,'position')
    tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
    speed = [0; diff(behaviour.wheel_raw_input*tick_to_cm_conversion)];
    speed(speed<-100) = 0;% big change in speed often due to teleportation or wheel tick resetting
    speed(speed>100) = 0;
    behaviour.speed = speed./mean(diff(behaviour.wheel_time))';%
    behaviour.speed = smoothdata(behaviour.speed, 'movmean', [10, 0]);
    % behaviour.speed = speed';

    %% Virtual position and virtual speed
    % Convert raw wheel wheel position into virtual position (0 - 140cm for both tracks)
    % Position for grey screen period outside of running blocks will be nan
    % Position at the end of each track (also looks grey screen) is still
    % logged

    % Add track 1 position
    behaviour.position = nan(1,length(behaviour.wheel_position));
    behaviour.track_ID = zeros(1,length(behaviour.wheel_position));
    behaviour.position(find(behaviour.wheel_position <= -990)) = abs(behaviour.wheel_position(find(behaviour.wheel_position <= -990))+1140);
    behaviour.track_ID((find(behaviour.wheel_position <= -990))) = 1;

    % Add track 2 position
    behaviour.position(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10)) ...
        = abs(behaviour.wheel_position(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10))+140);
    behaviour.track_ID(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10)) = 2;

    behaviour.position(behaviour.position>140) = 140;


    behaviour.sglxTime = behaviour.wheel_time';
end


if ~isempty(task_info)
    fields = fieldnames(task_info);
    for i = 1:numel(fields)
        behaviour.(fields{i}) = task_info.(fields{i});
    end
end


%%%%%%%% Get per lap positona and speed

behaviour.complete_laps_id = [];
behaviour.aborted_laps_id = [];
start_indices = [];

x = behaviour.position;    % position data
t = behaviour.sglxTime;    % time

for nlap = 1:length(behaviour.track_ID_all)
    [~, start_indices(nlap)] = min(abs(t - behaviour.start_time_all(nlap)));
end

for nlap = 1:length(start_indices)

    if nlap < length(start_indices)
        current_lap_x = x(start_indices(nlap):start_indices(nlap+1));
        current_lap_t = t(start_indices(nlap):start_indices(nlap+1));
    else
        current_lap_x = x(start_indices(nlap):end);
        current_lap_t = t(start_indices(nlap):end);
    end

    if length(current_lap_x) > 1 && length(current_lap_x(~isnan(current_lap_x))) > 1
        on_track_x = current_lap_x(~isnan(current_lap_x));
        on_track_t = current_lap_t(~isnan(current_lap_x));

        if length(on_track_x) > 30
            end_frame = 30;

            if ~isempty(find(diff(on_track_x(1:end_frame)) > 5))
                jump_index = find(diff(on_track_x(1:end_frame)) > 5);
                jump_index = jump_index(end);
                on_track_t(1:jump_index) = [];
                on_track_x(1:jump_index) = [];
            end

            if length(on_track_x) > 1
                if ~isempty(find(diff(on_track_x(1:end_frame)) < -5))
                    jump_index = find(diff(on_track_x(1:end_frame)) < -5);
                    jump_index = jump_index(end);
                    on_track_t(1:jump_index) = [];
                    on_track_x(1:jump_index) = [];
                end
            end
        end

        if sum(on_track_x == 0) > 0
            start_position = find(on_track_x == 0);
            if start_position < length(on_track_x) - 10
                on_track_x = on_track_x(start_position(1):end);
                on_track_t = on_track_t(start_position(1):end);
            end
        end

        if on_track_x(end) ~= on_track_x(end-1)
            on_track_x(end) = on_track_x(end-1);
            on_track_t(end) = on_track_t(end-1);
        end

        [last_position, last_position_index] = max(on_track_x);

        if last_position_index * mean(diff(on_track_t)) < 0.1
            on_track_x(1:last_position_index) = [];
            on_track_t(1:last_position_index) = [];
            [last_position, last_position_index] = max(on_track_x);
        end

        if isfield(behaviour, 'reward_delivery_time')
            if ~isempty(ismember(behaviour.reward_delivery_time, on_track_t'))
                behaviour.rewarded_lap_id(ismember(behaviour.reward_delivery_time, on_track_t')') = nlap;
            end
        end

        if last_position >= 139
            behaviour.end_time_all(nlap) = on_track_t(last_position_index);
            behaviour.complete_laps_id = [behaviour.complete_laps_id, nlap];
        else
            if last_position_index * mean(diff(on_track_t)) < 0.1
                on_track_x(1:last_position_index) = [];
                on_track_t(1:last_position_index) = [];
                [last_position, last_position_index] = max(on_track_x);
            end

            behaviour.end_time_all(nlap) = on_track_t(end);
            behaviour.aborted_laps_id = [behaviour.aborted_laps_id, nlap];
        end
    end
end





clear lap_times
lap_times(1).lap(1).track1_lick_latency=[];
lap_times(1).lap(1).track2_lick_latency=[];
lap_times(2).lap(1).track1_lick_latency=[];
lap_times(2).lap(1).track2_lick_latency=[];

for nlap = 1:length(behaviour.lap_ID_all)
    current_track = behaviour.track_ID_all(nlap);
    current_lap = behaviour.lap_ID_all(nlap);
    if isfield(behaviour,'reward_track1_count')
        reward_lap = find(behaviour.(sprintf('reward_track%i_count',current_track)) == current_lap);
    else
        reward_lap = find(behaviour.(sprintf('track%i_count',current_track)) == current_lap);
    end

    if current_lap == 0
        continue
    end

    if current_track == 1
        correct_position = 1140;
    elseif current_track == 2
        correct_position = 140;
    end

    if isempty(reward_lap) | isempty(behaviour.reward_position)
        lap_times(current_track).lap(current_lap).reward_position = nan;
        lap_times(current_track).lap(current_lap).reward_type = 0;

        reward_position(nlap,1) = nan;
        reward_type(nlap,1) = 0;
    else
        reward_lap = reward_lap(1); % Only the first time, it turns to this lap count
        reward_position(nlap,1) = behaviour.reward_position(reward_lap,1)+correct_position;
        reward_type(nlap,1) = behaviour.reward_type(reward_lap,1);

        lap_times(current_track).lap(current_lap).reward_position = reward_position(nlap,1);
        lap_times(current_track).lap(current_lap).reward_type = reward_type(nlap,1);
    end
    
    % lap_times(current_track).lap(current_lap).x = 

    % lap_times(current_track).lap(current_lap).trial_type = behaviour.trial_type(nlap); % if 1 means active only trial, 2 means hybrid
    lap_times(current_track).lap(current_lap).lap_ID_all = behaviour.lap_ID_all(nlap); % if 1 means active only trial, 2 means hybrid

    current_lap_licks = find(behaviour.(sprintf('lick_track%i_count',current_track)) == current_lap & ...
        behaviour.lick_track_ID == current_track);

    track1_licks = find(behaviour.lick_state(current_lap_licks) == 1); % track1 means left lick
    track2_licks = behaviour.lick_state(current_lap_licks) == 2; % track2 means right lick

    track1_lick_position{nlap} = behaviour.lick_position(current_lap_licks(track1_licks))+correct_position;
    track2_lick_position{nlap} = behaviour.lick_position(current_lap_licks(track2_licks))+correct_position;


    if max(current_lap_licks) <= length(behaviour.lick_time)
        track1_lick_time{nlap} = behaviour.lick_time(current_lap_licks(track1_licks));
        track2_lick_time{nlap} = behaviour.lick_time(current_lap_licks(track2_licks));
    elseif isempty(current_lap_licks)
        track1_lick_time{nlap} = [nan];
        track2_lick_time{nlap} = [nan];
    elseif current_lap_licks(end-1) < length(behaviour.lick_time)
        track1_lick_time{nlap} = [behaviour.lick_time(current_lap_licks(track1_licks(1:end-1)));nan];
        track2_lick_time{nlap} = [behaviour.lick_time(current_lap_licks(track2_licks(1:end-1)));nan];

    else
        track1_lick_time{nlap} = [nan];
        track2_lick_time{nlap} = [nan];
  
        % track1_lick_time{nlap} = [nan];
        % track2_lick_time{nlap} = [nan];
    end
    
    lap_times(current_track).lap(current_lap).track1_lick_time = track1_lick_time{nlap};
    lap_times(current_track).lap(current_lap).track2_lick_time = track2_lick_time{nlap};



    if reward_lap <= length(reward_position)
        if ~isnan(reward_position(reward_lap,1))
            if ~isempty(track1_lick_position{nlap})
                if ~isempty(find(track1_lick_position{nlap}(1) > reward_position(reward_lap,1)))
                    lap_times(current_track).lap(current_lap).track1_lick_latency = track1_lick_time{nlap}(find(track1_lick_position{nlap}(find(track1_lick_position{nlap}(1) > reward_position(reward_lap,1))) > reward_position(reward_lap,1))) - behaviour.reward_delivery_time(reward_lap);
                end
            end


            if ~isempty(track2_lick_position{nlap})
                if ~isempty(find(track2_lick_position{nlap}(1) > reward_position(reward_lap,1)))
                    lap_times(current_track).lap(current_lap).track2_lick_latency = track2_lick_time{nlap}(find(track2_lick_position{nlap}(find(track2_lick_position{nlap}(1) > reward_position(reward_lap,1))) > reward_position(reward_lap,1))) - behaviour.reward_delivery_time(reward_lap);
                end
            end
        end
    end




    if isfield(behaviour,'manual_reward_time')
        this_lap =   find(behaviour.manual_reward_time>behaviour.start_time_all(nlap) & behaviour.manual_reward_time<behaviour.end_time_all(nlap));
        if isempty(this_lap)
            session_behaviour_summary.manual_trial{current_track}(current_lap) =0;
            session_behaviour_summary.manual_position{current_track}{current_lap} =[];
        else
            session_behaviour_summary.manual_trial{current_track}(current_lap) = 1;
            session_behaviour_summary.manual_position{current_track}{current_lap} = behaviour.manual_reward_position(this_lap)+correct_position;
        end

    end

    session_behaviour_summary.L_lick{current_track}{current_lap} =track1_lick_position{nlap};
    session_behaviour_summary.R_lick{current_track}{current_lap} =track2_lick_position{nlap};
end



% session_behaviour_summary.manual_trial


%% plot all trials

colorlines = [ ...
    0.2, 0.4, 0.8;    % blue
    0.85, 0.2, 0.2;   % red
];

landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
landmarkVertex(4,:,1) = [86 94 94 86]; %4th
landmarkVertex(5,:,1) = [106 114 114 106]; %5th
landmarkVertex(1,:,2) = [0 0 length(behaviour.lap_ID_all) length(behaviour.lap_ID_all)];
landmarkVertex(2,:,2) = [0 0 length(behaviour.lap_ID_all) length(behaviour.lap_ID_all)];
landmarkVertex(3,:,2) = [0 0 length(behaviour.lap_ID_all) length(behaviour.lap_ID_all)];
landmarkVertex(4,:,2) = [0 0 length(behaviour.lap_ID_all) length(behaviour.lap_ID_all)];
landmarkVertex(5,:,2) = [0 0 length(behaviour.lap_ID_all) length(behaviour.lap_ID_all)];

% colorlines reference:
%   Blue = [0.2, 0.4, 0.8]
%   Red  = [0.85, 0.2, 0.2]

% === Blue track (track 1) ===
patch_color{1} = [0.8, 0.9, 1.0;    % light blue
                  0.4, 0.6, 0.9;    % medium-light blue
                  0.2, 0.4, 0.8;    % fixed deep blue (colorlines)
                  0.4, 0.6, 0.9;    % repeat
                  0.2, 0.4, 0.8];   % repeat fixed color

% === Red track (track 2) ===
patch_color{2} = [1.0, 0.8, 0.8;    % light red
                  0.95, 0.4, 0.4;   % medium-light red
                  0.85, 0.2, 0.2;   % fixed deep red (colorlines)
                  0.95, 0.4, 0.4;   % repeat
                  0.85, 0.2, 0.2];  % repeat fixed color



fig = figure;
fig.Position = [580 150 1100 830];
fig.Name = sprintf('Lick behaviour summary %s',session_name);

rewardColor = [77/256,175/256,74/256;
    228/256,26/256,28/256];

xline(100,'r','LineWidth',3)
for nlap=1:length(behaviour.lap_ID_all)

    if reward_type(nlap) > 2
        rewardType_tmp = 2;
    else
        rewardType_tmp = 1;
    end

    if rewardType_tmp == 2
        scatter(reward_position(nlap),nlap,30-0.3,'o','MarkerFaceColor',patch_color{behaviour.track_ID_all(nlap)}(5,:),'MarkerEdgeColor','k');
    else
        scatter(reward_position(nlap),nlap,30-0.7,'o','MarkerFaceColor',rewardColor(1,:),'MarkerEdgeColor','k');
    end
    hold on
    if ~isempty(track1_lick_position{nlap})
        scatter(track1_lick_position{nlap},nlap-0.3, 10,'*','MarkerEdgeColor',patch_color{1}(5,:))
    end
    if ~isempty(track2_lick_position{nlap})
        scatter(track2_lick_position{nlap},nlap-0.7, 10,'*','MarkerEdgeColor',patch_color{2}(5,:))
    end
    % set(gca,'fontsize',14)
end


track_ids = behaviour.track_ID_all;

n_laps = length(track_ids);
if isempty(behaviour.track_ID_all)
    session_behaviour_summary = [];
end
% Find boundaries where track ID changes
transitions = [1; find(diff(track_ids) ~= 0) + 1; n_laps+1];

hold on
for b = 1:(length(transitions)-1)
    block_start = transitions(b);
    block_end   = transitions(b+1) - 1;
    track_id    = track_ids(block_start);

    % Shared Y values for all patches in block
    y_vals = [block_start-0.9, block_start-0.9, block_end+0.1, block_end+0.1];

    % Plot 5 column patches
    for iPatch = 1:5
        x_vals = landmarkVertex(iPatch,:,1);  % constant X coords
        patch(x_vals, y_vals, patch_color{track_id}(iPatch,:), ...
            'FaceAlpha', 0.6, 'EdgeAlpha', 0);  % Reduced alpha
    end
end
xline([100],'k','LineWidth',1)
ylim([0 length(behaviour.lap_ID_all)])
xlim([0,140])
xticks([0 30 50 70 90 110])
set(gca,'TickDir','out','Box','off','FontSize',14)
xlabel('Position (cm)')
ylabel('Lap')


%% track 1 and track 2 lick histogram

for track_id = 1:2
    this_track_laps = find(behaviour.track_ID_all == track_id);
    plot_count = 1;

    lickBin = 0:2:140;

    fig = figure;
    fig.Position = [580 150 1100 830*3/4];
    fig.Name = sprintf('Track %i lick bias %s',track_id,session_name);
    Trial_types = {'All','Passive reward','Active reward','Active only'};

    for type = 1:3

        if type == 1
            this_track_laps = find(behaviour.track_ID_all == track_id);
        elseif type == 2
            this_track_laps = intersect(find(behaviour.track_ID_all == track_id), find(reward_type <= 2));

            session_behaviour_summary.passive_laps{track_id} = behaviour.lap_ID_all(this_track_laps);
        elseif type == 3
            this_track_laps = intersect(find(behaviour.track_ID_all == track_id), find(reward_type > 2));
            session_behaviour_summary.active_laps{track_id} = behaviour.lap_ID_all(this_track_laps);
            % elseif type == 4
            %     this_track_laps = intersect(find(behaviour.track_ID_all == track_id), find(behaviour.trial_type == 1));
        end

        if isempty(this_track_laps)
            continue
        end

        ymin = [];
        ymax = [];

        subplot(3,2,plot_count)
        if ~isempty(cat(1,track1_lick_position{this_track_laps}))
            licks_histcount = histcounts(cat(1,track1_lick_position{this_track_laps}), lickBin) / length(this_track_laps);
            plot(lickBin(2:end), licks_histcount, 'Color', patch_color{1}(5,:))
            ymax = max(licks_histcount);
        else
            plot(lickBin(2:end), zeros(1,length(lickBin(2:end))), 'Color', patch_color{1}(5,:));
            ymax = 0.2;
        end
        hold on

        if ~isempty(cat(1,track2_lick_position{this_track_laps}))
            licks_histcount = histcounts(cat(1,track2_lick_position{this_track_laps}), lickBin) / length(this_track_laps);
            plot(lickBin(2:end), licks_histcount, 'Color', patch_color{2}(5,:))
            ymax = max([ymax max(licks_histcount)]);
        else
            plot(lickBin(2:end), zeros(1,length(lickBin(2:end))), 'Color', patch_color{2}(5,:));
        end

        if ymax==0
            ymax = 1;
            
        end

        xline(100, 'Color', patch_color{track_id}(5,:), 'LineWidth', 5)

        hold on
        landmarkVertex = zeros(5,4,2);
        landmarkVertex(1,:,1) = [26 34 34 26];
        landmarkVertex(2,:,1) = [46 54 54 46];
        landmarkVertex(3,:,1) = [66 74 74 66];
        landmarkVertex(4,:,1) = [86 94 94 86];
        landmarkVertex(5,:,1) = [106 114 114 106];
        landmarkVertex(:,:,2) = repmat([0 0 1.5*ymax 1.5*ymax], [5 1]);

        for iPatch = 1:5
            patch(landmarkVertex(iPatch,:,1), landmarkVertex(iPatch,:,2), patch_color{track_id}(iPatch,:), ...
                'FaceAlpha', 0.3, 'EdgeAlpha', 0)
        end

        xlim([0,140])
        ylim([0 1.5*ymax])
        xlabel('Position')
        ylabel('Mean lick number')
        xticks([0 30 50 70 90 110])
        title(sprintf('all track %i licks (%s)', track_id, Trial_types{type}))
        set(gca,'fontsize',14)
        plot_count = plot_count + 1;

        %% subplot for first licks
        subplot(3,2,plot_count)
        track1_first_click = nan(1,length(this_track_laps));
        track2_first_click = nan(1,length(this_track_laps));
        first_lick = nan(1,length(this_track_laps));

        for nlap = 1:length(this_track_laps)
            if ~isempty(track1_lick_position{this_track_laps(nlap)})
                track1_first_click(nlap) = track1_lick_position{this_track_laps(nlap)}(1);
            end
            if ~isempty(track2_lick_position{this_track_laps(nlap)})
                track2_first_click(nlap) = track2_lick_position{this_track_laps(nlap)}(1);
            end

            if isnan(track1_first_click(nlap)) && ~isnan(track2_first_click(nlap))
                first_lick(nlap) = track2_first_click(nlap);
            elseif ~isnan(track1_first_click(nlap)) && isnan(track2_first_click(nlap))
                first_lick(nlap) = track1_first_click(nlap);
            elseif ~isnan(track1_first_click(nlap)) && ~isnan(track2_first_click(nlap))
                first_lick(nlap) = min(track1_first_click(nlap), track2_first_click(nlap));
            end
        end

        ymax = [];
        if sum(~isnan(track1_first_click)) > 0
            licks_histcount = histcounts(track1_first_click, lickBin) / length(this_track_laps);
            plot(lickBin(2:end), licks_histcount, 'Color', patch_color{1}(5,:))
            ymax = max(licks_histcount);
        else
            plot(lickBin(2:end), zeros(1,length(lickBin(2:end))), 'Color', patch_color{1}(5,:));
            ymax = 0.2;
        end
        hold on

        if sum(~isnan(track2_first_click)) > 0
            licks_histcount = histcounts(track2_first_click, lickBin) / length(this_track_laps);
            plot(lickBin(2:end), licks_histcount, 'Color', patch_color{2}(5,:))
            ymax = max([ymax max(licks_histcount)]);
        else
            plot(lickBin(2:end), zeros(1,length(lickBin(2:end))), 'Color', patch_color{2}(5,:));
        end

        if ymax==0
            ymax = 1;
            
        end
        
        xline(100, 'Color', patch_color{track_id}(5,:), 'LineWidth', 5)

        hold on
        landmarkVertex(:,:,2) = repmat([0 0 1.5*ymax 1.5*ymax], [5 1]);
        for iPatch = 1:5
            patch(landmarkVertex(iPatch,:,1), landmarkVertex(iPatch,:,2), patch_color{track_id}(iPatch,:), ...
                'FaceAlpha', 0.3, 'EdgeAlpha', 0)
        end

        xlim([0,140])
        ylim([0 1.5*ymax])
        xlabel('Position')
        ylabel('Proportion of first licks')
        title(sprintf('track %i first licks (%s)', track_id, Trial_types{type}))
        set(gca,'fontsize',14)
        plot_count = plot_count + 1;


        if type == 1
        
            session_behaviour_summary.L_first_lick{track_id} = track1_first_click;
            session_behaviour_summary.R_first_lick{track_id} = track2_first_click;
        end
    end

    legend('Left licks','Right licks','Color','none')
    sgtitle(sprintf('Track %i lick histogram', track_id))

end




fig = figure
fig.Position = [300 150 600 600];
fig.Name = sprintf('Lap running speed vs position %s',session_name);

x_bin_width = 2;
x_bin_edges = 0:x_bin_width:140; % forces x_bins to be from 0 to 140cm
x_bin_centres = [(x_bin_edges(2)-x_bin_width/2):x_bin_width:(x_bin_edges(end-1)+x_bin_width/2)];
x = 1:x_bin_width:140;
ymax = [];

for track_id = 1:2
    speed_position = [];
    lap_mean_speed = [];
    subplot(2,1,track_id)
    hold on

    for nlap = 1:length(lap_times(track_id).lap)
        
        start_time = behaviour.start_time_all(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id);
        end_time =  behaviour.end_time_all(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id);
        if end_time < start_time
            % if find(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id)+1 <= length(behaviour.lap_ID_all)
            %     end_time = behaviour.start_time_all(find(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id)+1);
            % end
             lap_mean_speed(nlap) = nan;
             speed_position(nlap,:) = nan;
        else
            [~, start_idx] = min(abs(behaviour.sglxTime - start_time));
            [~, end_idx]   = min(abs(behaviour.sglxTime - end_time));

            lap_mean_speed(nlap) = median(behaviour.speed(start_idx:end_idx));
            [N,edges,x_bins] = histcounts(behaviour.position(start_idx:end_idx),x_bin_edges);
            this_lap_index = start_idx:end_idx;

            for nbin = 1:max(x_bins)
                speed_position(nlap,nbin) = median(behaviour.speed(this_lap_index(x_bins == nbin)));
            end
        end
    end

    session_behaviour_summary.lap_speed{track_id} = speed_position;
    session_behaviour_summary.x_bins = x_bin_centres;
    session_behaviour_summary.x_bin_edges= x_bin_edges;

    
    %             [~,sorted_lap_id{track_id}] = sort(lap_mean_speed);
    low_speed_laps = find(lap_mean_speed <= median(lap_mean_speed,'omitnan'));
    high_speed_laps = find(lap_mean_speed > median(lap_mean_speed,'omitnan'));
    
    speed_LSE = nanmean(speed_position(low_speed_laps,:))- nanstd(speed_position(low_speed_laps,:))/sqrt(length(lap_times(track_id).lap));
    speed_USE = nanmean(speed_position(low_speed_laps,:))+ nanstd(speed_position(low_speed_laps,:))/sqrt(length(lap_times(track_id).lap));

    %             plot(x, speed_LSE, 'k--', 'LineWidth', 1);hold on;
    %             plot(x, speed_USE, 'k--', 'LineWidth', 1);
    hold on
    plot(x, nanmean(speed_position(low_speed_laps,:)), 'Color',patch_color{track_id}(4,:), 'LineWidth', 1);hold on
    x2 = [x, fliplr(x)];
    inBetween = [speed_LSE, fliplr(speed_USE)];
    h(1) = fill(x2, inBetween, patch_color{track_id}(1,:),'FaceAlpha',0.7,'EdgeColor','none');


    speed_LSE = nanmean(speed_position(high_speed_laps,:))- nanstd(speed_position(high_speed_laps,:))./sqrt(length(lap_times(track_id).lap));
    speed_USE = nanmean(speed_position(high_speed_laps,:))+ nanstd(speed_position(high_speed_laps,:))./sqrt(length(lap_times(track_id).lap));
    ymax = 65;

    %             plot(x, speed_LSE, 'k--', 'LineWidth', 1);hold on;
    plot(x, nanmean(speed_position(high_speed_laps,:)), 'Color',patch_color{track_id}(5,:), 'LineWidth', 1);
    x2 = [x, fliplr(x)];
    inBetween = [speed_LSE, fliplr(speed_USE)];
    h(2) = fill(x2, inBetween, patch_color{track_id}(5,:),'FaceAlpha',0.2,'EdgeColor','none');
    ylim([0 ymax])
    xlabel('Position (cm)')
    ylabel('Speed cm/s')
    hold on

    xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',2)
    hold on
    landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
    landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
    landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
    landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
    landmarkVertex(4,:,1) = [86 94 94 86]; %4th
    landmarkVertex(5,:,1) = [106 114 114 106]; %5th
    landmarkVertex(:,:,2) = repmat([-1.25*ymax -1.25*ymax 1.25*ymax 1.25*ymax],[5 1]); %landmark Y coord

    for iPatch = 1:5
        patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',0.3,'EdgeAlpha',0)
    end
    % legend('Left licks','Right licks')
    xlim([0,140])
    xticks([0 30 50 70 90 110])

    set(gca,"TickDir","out",'box', 'off','Color','none')
    legend([h(1),h(2)],["Low speed","High speed"],'Color','none');
    set(gca,'fontsize',14)
end





fig = figure
fig.Position = [300 150 600 600];
fig.Name = sprintf('Lap running speed vs position (passive vs active) %s',session_name);

x_bin_width = 2;
x_bin_edges = 0:x_bin_width:140; % forces x_bins to be from 0 to 140cm
x_bin_centres = [(x_bin_edges(2)-x_bin_width/2):x_bin_width:(x_bin_edges(end-1)+x_bin_width/2)];
x = 1:x_bin_width:140;
ymax = [];

for track_id = 1:2
    speed_position = [];
    lap_mean_speed = [];
    subplot(2,1,track_id)
    hold on

    for nlap = 1:length(lap_times(track_id).lap)

        start_time = behaviour.start_time_all(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id);
        end_time =  behaviour.end_time_all(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id);
        if end_time < start_time
            % if find(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id)+1 <= length(behaviour.lap_ID_all)
            %     end_time = behaviour.start_time_all(find(behaviour.lap_ID_all==lap_times(track_id).lap(nlap).lap_ID_all & behaviour.track_ID_all==track_id)+1);
            % end
            lap_mean_speed(nlap) = nan;
            speed_position(nlap,:) = nan;
        else
            [~, start_idx] = min(abs(behaviour.sglxTime - start_time));
            [~, end_idx]   = min(abs(behaviour.sglxTime - end_time));

            lap_mean_speed(nlap) = median(behaviour.speed(start_idx:end_idx));
            [N,edges,x_bins] = histcounts(behaviour.position(start_idx:end_idx),x_bin_edges);
            this_lap_index = start_idx:end_idx;

            for nbin = 1:max(x_bins)
                speed_position(nlap,nbin) = median(behaviour.speed(this_lap_index(x_bins == nbin)));
            end

        end

    end
    %             [~,sorted_lap_id{track_id}] = sort(lap_mean_speed);
    low_speed_laps = session_behaviour_summary.passive_laps{track_id};
    high_speed_laps = session_behaviour_summary.active_laps{track_id};

    if length(high_speed_laps)<2
        continue
    end

    speed_LSE = nanmean(speed_position(low_speed_laps,:))- nanstd(speed_position(low_speed_laps,:))/sqrt(length(lap_times(track_id).lap));
    speed_USE = nanmean(speed_position(low_speed_laps,:))+ nanstd(speed_position(low_speed_laps,:))/sqrt(length(lap_times(track_id).lap));

    %             plot(x, speed_LSE, 'k--', 'LineWidth', 1);hold on;
    %             plot(x, speed_USE, 'k--', 'LineWidth', 1);
    hold on
    plot(x, nanmean(speed_position(low_speed_laps,:)), 'Color',patch_color{track_id}(4,:), 'LineWidth', 1);hold on
    x2 = [x, fliplr(x)];
    inBetween = [speed_LSE, fliplr(speed_USE)];
    h(1) = fill(x2, inBetween, patch_color{track_id}(1,:),'FaceAlpha',0.7,'EdgeColor','none');


    speed_LSE = nanmean(speed_position(high_speed_laps,:))- nanstd(speed_position(high_speed_laps,:))./sqrt(length(lap_times(track_id).lap));
    speed_USE = nanmean(speed_position(high_speed_laps,:))+ nanstd(speed_position(high_speed_laps,:))./sqrt(length(lap_times(track_id).lap));
    ymax = 65;

    %             plot(x, speed_LSE, 'k--', 'LineWidth', 1);hold on;
    plot(x, nanmean(speed_position(high_speed_laps,:)), 'Color',patch_color{track_id}(5,:), 'LineWidth', 1);
    x2 = [x, fliplr(x)];
    inBetween = [speed_LSE, fliplr(speed_USE)];
    h(2) = fill(x2, inBetween, patch_color{track_id}(5,:),'FaceAlpha',0.3,'EdgeColor','none');
    ylim([0 ymax])
    xlabel('Position (cm)')
    ylabel('Speed cm/s')
    hold on

    xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',2)
    hold on
    landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
    landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
    landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
    landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
    landmarkVertex(4,:,1) = [86 94 94 86]; %4th
    landmarkVertex(5,:,1) = [106 114 114 106]; %5th
    landmarkVertex(:,:,2) = repmat([-1.25*ymax -1.25*ymax 1.25*ymax 1.25*ymax],[5 1]); %landmark Y coord

    for iPatch = 1:5
        patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',0.3,'EdgeAlpha',0)
    end
    % legend('Left licks','Right licks')
    xlim([0,140])
    xticks([0 30 50 70 90 110])

    set(gca,"TickDir","out",'box', 'off','Color','none')
    legend([h(1),h(2)],["Passive","Active"],'Color','none');
    set(gca,'fontsize',14)
end






