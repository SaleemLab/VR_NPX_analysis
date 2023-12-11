clear all
%Plot training data
mouse = 'M23028';
% mouse = 'M23032';
Year = '23';
Month = '06';
Day = '02';
% date = '2023-05-25';
date = ['20',Year,'-',Month,'-',Day];

% path = ['\\live.rd.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\',mouse,'\training\',Day,Month,Year,'\'];
path = ['Z:\ibn-vision\DATA\SUBJECTS\',mouse,'\training\',Day,Month,Year,'\'];
% behaviour = readTrainingData(path,date);
Behaviour = read_training_data(path,date);

[Behaviour] = read_behaviour_data_without_alignment(path,date)
cd(path)
save bonsai_behaviour Behaviour

%% Extract laps
cd(path)
% data.linear(1).linear = nan(1,length(Behaviour.wheel_position));
% data.linear(1).linear(find(Behaviour.wheel_position <= -990)) = abs(Behaviour.wheel_position(find(Behaviour.wheel_position <= -990))+1140);
% data.linear(1).linear(data.linear(1).linear>140) = 140;
linear = Behaviour.wheel_position; linear (linear < 300) = nan; linear(linear > 140) = 140;
data.linear(1).linear = linear;
Behaviour.track1_position = data.linear(1).linear;
data.x = nan(1,length(Behaviour.wheel_position));
data.x(~isnan(data.linear(1).linear)) = data.linear(1).linear(~isnan(data.linear(1).linear));

if sum(find(Behaviour.wheel_position >= -141 & Behaviour.wheel_position < 0)) > 0 % If there is track 2
    data.linear(2).linear = nan(1,length(Behaviour.wheel_position));
    data.linear(2).linear(find(Behaviour.wheel_position >= -141 & Behaviour.wheel_position < 10)) ...
        = abs(Behaviour.wheel_position(find(Behaviour.wheel_position >= -141 & Behaviour.wheel_position < 10))+140);
    data.linear(2).linear(data.linear(2).linear>140) = 140;
    Behaviour.track2_position = data.linear(2).linear;
    data.x(~isnan(data.linear(2).linear)) = data.linear(2).linear(~isnan(data.linear(2).linear));
end

data.t = Behaviour.wheel_time';
data.v = zeros(1,length(Behaviour.wheel_time));
data.v(2:end) = diff(Behaviour.wheel_raw_input'*0.0613)/mean(diff(data.t));% instaneous speed
% data.v(2:end) = diff(data.x);
data.v_cm = data.v;
position = data;

save bonsai_behaviour Behaviour
save extracted_position position -v7.3

lap_times = extract_laps_masa(1,Behaviour,position)
save extracted_laps lap_times
close all
if isempty(Behaviour.start_time_all)
    all_laps_start_time = [lap_times(1).start lap_times(2).start];
    [~,correct_order] = sort(all_laps_start_time);
    Behaviour.start_time_all = all_laps_start_time(correct_order);

    track_ID_all = [ones(1,length(lap_times(1).start)) 2*ones(1,length(lap_times(2).start))];
    Behaviour.track_ID_all = track_ID_all(correct_order);
    lap_ID_all = [lap_times(1).lap_id lap_times(2).lap_id];
    Behaviour.lap_ID_all= lap_ID_all(correct_order);
end

%% Saving lap info


for nlap = 1:length(Behaviour.lap_ID_all)
    current_track = Behaviour.track_ID_all(nlap);
    current_lap = Behaviour.lap_ID_all(nlap);
    if current_lap == 0
        continue
    end
    reward_lap = find(Behaviour.(sprintf('track%i_count',current_track)) == current_lap);

    if current_track == 1
        correct_position = 1140;
    elseif current_track == 2
        correct_position = 140;
    end

    if isempty(reward_lap)
        lap_times(current_track).lap(current_lap).reward_position = nan;
        lap_times(current_track).lap(current_lap).reward_type = 0;

        reward_position(nlap,1) = nan;
        reward_type(nlap,1) = 0;
    else
        reward_lap = reward_lap(1); % Only the first time, it turns to this lap count
        reward_position(nlap,1) = Behaviour.reward_position(reward_lap,1)+correct_position;
        reward_type(nlap,1) = Behaviour.reward_type(reward_lap,1);

        lap_times(current_track).lap(current_lap).reward_position = reward_position(nlap,1);
        lap_times(current_track).lap(current_lap).reward_type = reward_type(nlap,1);
    end

    current_lap_licks = find(Behaviour.(sprintf('lick_track%i_count',current_track)) == current_lap & ...
        Behaviour.lick_track_ID == current_track);

    track1_licks = find(Behaviour.lick_state(current_lap_licks) == 1);
    track2_licks = Behaviour.lick_state(current_lap_licks) == 2;

    track1_lick_time{nlap} = Behaviour.lick_time(current_lap_licks(track1_licks));
    track2_lick_time{nlap} = Behaviour.lick_time(current_lap_licks(track2_licks));

    lap_times(current_track).lap(current_lap).track1_lick_time = track1_lick_time{nlap};
    lap_times(current_track).lap(current_lap).track2_lick_time = track2_lick_time{nlap};

    track1_lick_position{nlap} = Behaviour.lick_position(current_lap_licks(track1_licks))+correct_position;
    track2_lick_position{nlap} = Behaviour.lick_position(current_lap_licks(track2_licks))+correct_position;

    lap_times(current_track).lap(current_lap).track1_lick_position =  track1_lick_position{nlap};
    lap_times(current_track).lap(current_lap).track2_lick_position =  track2_lick_position{nlap};
end

%%


lick_ratio = [];
lick_ratio_all = [];
lickBin = 0:5:140;
for track_id = 1:2
    for type = 1

        if type == 1
            % All trials
            this_track_laps = find(Behaviour.track_ID_all == track_id);
        elseif type == 2
            % Find all passive trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type <= 2));

        elseif type == 3
            % Find all actively licked trials during hybrid trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type > 2));

        elseif type == 4
            % Find active only trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(Behaviour.trial_type == 1));
        end

        %         if track_id == 1
        %             this_track_licks = track1_lick_position;
        %         elseif track_id == 2
        %             this_track_licks = track2_lick_position;
        %         end

        for nlap = 1:length(this_track_laps)


            left_licks = track1_lick_position{this_track_laps(nlap)};
            right_licks = track2_lick_position{this_track_laps(nlap)};
            left_reward_lick = [];
            right_reward_lick = [];
            lick_ratio{track_id}(nlap) = 0;
            lick_ratio_all{track_id}(nlap) = 0;

            if sum(left_licks>=95 & left_licks<105) > 0
                left_reward_lick = left_licks(left_licks >= 95&left_licks <= 105);
                if ~isempty(left_reward_lick)
                    left_reward_lick = left_reward_lick(1);
                end
            end

            if sum(right_licks>=95 & right_licks<105) > 0
                right_reward_lick =  right_licks(right_licks >= 95&right_licks <= 105);
                if ~isempty(right_reward_lick)
                    right_reward_lick = right_reward_lick(1);
                end
            end


            if track_id == 1
                if sum(left_licks>95 & left_licks<105) > 0
                    licks_in_reward_zone{track_id} = this_track_laps(nlap);
                end

                if ~isempty([left_licks; right_licks])
                    lick_ratio_all{track_id}(nlap) = (length(left_licks)-length(right_licks))/(length(left_licks)+length(right_licks));
                else
                    lick_ratio{track_id}(nlap) = 0;
                    lick_ratio_all{track_id}(nlap) = 0;
                    lick_ratio_reward_zone{track_id}(nlap) = 0;
                    continue
                end

                left_licks = [left_licks(left_licks<95)];
                right_licks = [right_licks(right_licks<95)];

                if ~isempty([left_licks; right_licks])
                    lick_ratio{track_id}(nlap) = (length(left_licks)-length(right_licks))/(length(left_licks)+length(right_licks));
                    lick_ratio_reward_zone{track_id}(nlap) = (length(left_reward_lick)-length(right_reward_lick))/(length(left_reward_lick)+length(right_reward_lick));
                else
                    lick_ratio{track_id}(nlap) = 0;
                end

                if ~isempty([right_reward_lick; left_reward_lick])
                    lick_ratio_reward_zone{track_id}(nlap) = (length(left_reward_lick)-length(right_reward_lick))/(length(left_reward_lick)+length(right_reward_lick));
                else
                    lick_ratio_reward_zone{track_id}(nlap) = 0;
                end


            elseif track_id == 2
                if sum(right_licks>95 & right_licks<105) > 0
                    licks_in_reward_zone{track_id} = this_track_laps(nlap);
                end

                if ~isempty([left_licks; right_licks])
                    lick_ratio_all{track_id}(nlap) = (length(right_licks)-length(left_licks))/(length(left_licks)+length(right_licks));
                else
                    lick_ratio{track_id}(nlap) = 0;
                    lick_ratio_all{track_id}(nlap) = 0;
                    lick_ratio_reward_zone{track_id}(nlap) = 0;
                    continue
                end

                left_licks = [left_licks(left_licks<95)];
                right_licks = [right_licks(right_licks<95)];

                if ~isempty([left_licks; right_licks])
                    lick_ratio{track_id}(nlap) = (length(right_licks)-length(left_licks))/(length(left_licks)+length(right_licks));

                else
                    lick_ratio{track_id}(nlap) = 0;
                end


                if ~isempty([right_reward_lick; left_reward_lick])
                    lick_ratio_reward_zone{track_id}(nlap) = (length(right_reward_lick)-length(left_reward_lick))/(length(right_reward_lick)+length(left_reward_lick));
                else
                    lick_ratio_reward_zone{track_id}(nlap) = 0;
                end

                if type == 1
                    lap_times(track_id).lap(nlap).lick_ratio =  lick_ratio{track_id}(nlap);
                    lap_times(track_id).lap(nlap).lick_ratio_all =  lick_ratio_all{track_id}(nlap);
                end
            end

        end
    end

end

save('extracted_laps.mat','lap_times')

%% plot all trials

landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
landmarkVertex(4,:,1) = [86 94 94 86]; %4th
landmarkVertex(5,:,1) = [106 114 114 106]; %5th
landmarkVertex(1,:,2) = [0 0 length(Behaviour.lap_ID_all) length(Behaviour.lap_ID_all)]; %1st landmark Y coord
landmarkVertex(2,:,2) = [0 0 length(Behaviour.lap_ID_all) length(Behaviour.lap_ID_all)]; %2nd
landmarkVertex(3,:,2) = [0 0 length(Behaviour.lap_ID_all) length(Behaviour.lap_ID_all)]; %3rd
landmarkVertex(4,:,2) = [0 0 length(Behaviour.lap_ID_all) length(Behaviour.lap_ID_all)]; %4th
landmarkVertex(5,:,2) = [0 0 length(Behaviour.lap_ID_all) length(Behaviour.lap_ID_all)]; %5th

% patchColor = [251/256,180/256,174/256;
%               179/256,205/256,227/256;
%               204/256,235/256,197/256;
%               179/256,205/256,227/256;
%               204/256,235/256,197/256];
patch_color{2} = [161,218,180;
    44,127,184;
    37,52,148;
    44,127,184;
    37,52,148]/256;

patch_color{1} = [254,204,92;
    240,59,32;
    189,0,38;
    240,59,32;
    189,0,38]/256;

fig = figure;
fig.Position = [580 150 1100 830];
fig.Name = sprintf('Lick behaviour summary %s %s',mouse,date);
for iPatch = 1:5
    patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{1}(iPatch,:),'FaceAlpha',.3,'EdgeAlpha',0)
end
hold on

rewardColor = [77/256,175/256,74/256;
    228/256,26/256,28/256];

% nfig.position =

xline(100,'r','LineWidth',3)
for nlap=1:length(Behaviour.lap_ID_all)

    if reward_type(nlap) > 2
        rewardType_tmp = 2;
    else
        rewardType_tmp = 1;
    end

    if rewardType_tmp == 2
        scatter(reward_position(nlap),nlap,30-0.3,'o','MarkerFaceColor',patch_color{Behaviour.track_ID_all(nlap)}(5,:),'MarkerEdgeColor','k');
    else
        scatter(reward_position(nlap),nlap,30-0.7,'o','MarkerFaceColor',rewardColor(1,:),'MarkerEdgeColor','k');
    end
    hold on
    if ~isempty(track1_lick_position{nlap})
        scatter(track1_lick_position{nlap},nlap-0.3, 10,'*','MarkerEdgeColor',patch_color{1}(5,:))
    end
    if~isempty(track2_lick_position{nlap})
        scatter(track2_lick_position{nlap},nlap-0.7, 10,'*','MarkerEdgeColor',patch_color{2}(5,:))
    end

    landmarkVertex(:,:,2) = repmat([nlap-0.9 nlap-0.9 nlap-0.1 nlap-0.1],5,1);
    for iPatch = 1:5
        patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{Behaviour.track_ID_all(nlap)}(iPatch,:),'FaceAlpha',.5,'EdgeAlpha',0)
    end
    set(gca,'fontsize',14)
end
ylim([0 length(Behaviour.lap_ID_all)])
xlim([0,140])
xticks([0 30 50 70 90 110])


%% track 1 and track 2 lick histogram

for track_id = 1:2
    this_track_laps = find(Behaviour.track_ID_all == track_id);
    % track_2_laps = find(Behaviour.track_ID_all == 2);

    plot_count = 1;

    lickBin = 0:5:140;

    fig = figure;
    fig.Position = [580 150 1100 830];
    fig.Name = sprintf('Track %i lick bias %s %s',track_id,mouse,date);
    Trial_types = {'All','Passive reward','Active reward','Active only'};

    for type = 1:3

        if type == 1
            % All trials
            this_track_laps = find(Behaviour.track_ID_all == track_id);
        elseif type == 2
            % Find all passive trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type <= 2));

        elseif type == 3
            % Find all active trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type > 2));

        elseif type == 4
            % Find active only trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(Behaviour.trial_type == 1));
        end

        if isempty(this_track_laps)
            continue
        end
        ymin = [];
        ymax = [];

        subplot(4,2,plot_count)
        if ~isempty(cat(1,track1_lick_position{this_track_laps}))
            licks_histcount = histcounts(cat(1,track1_lick_position{this_track_laps}),lickBin)/length(this_track_laps);
            %     histogram(cat(1,track1_lick_position{track_1_laps}),lickBin,'FaceAlpha',0.5,'FaceColor',patch_color{1}(5,:),'Normalization','probability')
            %             bar(lickBin(2:end),licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{1}(5,:))
            plot(lickBin(2:end),licks_histcount,'Color',patch_color{1}(5,:))
            ymax= max(licks_histcount);
        else
            plot(lickBin(2:end),zeros(1,length(lickBin(2:end))),'Color',patch_color{1}(5,:));
            ymax = 0.2;
        end
        hold on

        if ~isempty(cat(1,track2_lick_position{this_track_laps}))
            licks_histcount = histcounts(cat(1,track2_lick_position{this_track_laps}),lickBin)/length(this_track_laps);
            %             bar(lickBin(2:end),-licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{2}(5,:))
            plot(lickBin(2:end),licks_histcount,'Color',patch_color{2}(5,:))
            %     histogram(cat(1,track2_lick_position{track_1_laps}),lickBin,'FaceAlpha',0.5,'FaceColor',patch_color{2}(5,:),'Normalization','probability')


            ymax= max([ymax max(licks_histcount)]);
        else % just plot zeors
            plot(lickBin(2:end),zeros(1,length(lickBin(2:end))),'Color',patch_color{2}(5,:));
            %     histogram(cat(1,track2_lick_position{track_1_laps}),lickBin,'FaceAlpha',0.5,'FaceColor',patch_color{2}(5,:),'Normalization','probability')

        end

        xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',5)

        hold on
        landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
        landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
        landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
        landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
        landmarkVertex(4,:,1) = [86 94 94 86]; %4th
        landmarkVertex(5,:,1) = [106 114 114 106]; %5th
        landmarkVertex(:,:,2) = repmat([0 0 1.5*ymax 1.5*ymax],[5 1]); %landmark Y coord

        for iPatch = 1:5
            patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',0.3,'EdgeAlpha',0)
        end
        % legend('Left licks','Right licks')
        xlim([0,140])
        ylim([0 1.5*ymax])
        xlabel('Position')
        ylabel('Mean lick number')
        xticks([0 30 50 70 90 110])
        title(sprintf('all track %i licks (%s)',track_id,Trial_types{type}))
        set(gca,'fontsize',14)
        plot_count = plot_count + 1;


        subplot(4,2,plot_count)
        track1_first_click = [];
        track2_first_click = [];
        for nlap = 1:length(this_track_laps)
            if isempty(track1_lick_position{this_track_laps(nlap)})
                track1_first_click(nlap) = nan;
            else
                track1_first_click(nlap) = track1_lick_position{this_track_laps(nlap)}(1);
            end

            if isempty(track2_lick_position{this_track_laps(nlap)})
                track2_first_click(nlap) = nan;
            else
                track2_first_click(nlap) = track2_lick_position{this_track_laps(nlap)}(1);
            end

            if isnan(track1_first_click(nlap)) & isnan(track2_first_click(nlap))
                first_lick(nlap) = nan;
            elseif ~isnan(track1_first_click(nlap)) & isnan(track2_first_click(nlap))
                first_lick(nlap) = track1_first_click(nlap);
            elseif isnan(track1_first_click(nlap)) & ~isnan(track2_first_click(nlap))
                first_lick(nlap) = track2_first_click(nlap);
            elseif track1_first_click(nlap) > track2_first_click(nlap)
                first_lick(nlap) = track2_first_click(nlap);
            elseif track1_first_click(nlap) < track2_first_click(nlap)
                first_lick(nlap) = track1_first_click(nlap);
            else
                first_lick(nlap) = nan;
            end
        end

        ymax = [];
        if sum(~isnan(track1_first_click)) > 0
            %     histogram(track1_first_click,lickBin,'FaceAlpha',0.3,'FaceColor',patch_color{1}(5,:),'Normalization','probability')
            licks_histcount = histcounts(track1_first_click,lickBin)/length(this_track_laps);
            %             bar(lickBin(2:end),licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{1}(5,:))
            plot(lickBin(2:end),licks_histcount,'Color',patch_color{1}(5,:))
            ymax= max(licks_histcount);
        else
            plot(lickBin(2:end),zeros(1,length(lickBin(2:end))),'Color',patch_color{1}(5,:));
            ymax = 0.2;
        end

        hold on
        if sum(~isnan(track2_first_click)) > 0
            %     histogram(track2_first_click,lickBin,'FaceAlpha',0.3,'FaceColor',patch_color{2}(5,:),'Normalization','probability')
            licks_histcount = histcounts(track2_first_click,lickBin)/length(this_track_laps);
            %             bar(lickBin(2:end),-licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{2}(5,:))
            plot(lickBin(2:end),licks_histcount,'Color',patch_color{2}(5,:))
            %                     ymax= max([ymax max(licks_histcount)]);

            ymax= max([ymax max(licks_histcount)]);
        else % just plot zeors
            plot(lickBin(2:end),zeros(1,length(lickBin(2:end))),'Color',patch_color{2}(5,:));
            %     histogram(cat(1,track2_lick_position{track_1_laps}),lickBin,'FaceAlpha',0.5,'FaceColor',patch_color{2}(5,:),'Normalization','probability')

        end

        xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',5)

        hold on
        landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
        landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
        landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
        landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
        landmarkVertex(4,:,1) = [86 94 94 86]; %4th
        landmarkVertex(5,:,1) = [106 114 114 106]; %5th
        landmarkVertex(:,:,2) = repmat([0 0 1.5*ymax 1.5*ymax],[5 1]); %landmark Y coord

        for iPatch = 1:5
            patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',0.3,'EdgeAlpha',0)
        end
        % legend('Left licks','Right licks')
        xlim([0,140])
        xticks([0 30 50 70 90 110])
        ylim([0 1.5*ymax])
        %         ylabel('Average lick number per lap')
        %         xlabel('Position (cm)')
        title(sprintf('track %i first licks (%s)',track_id,Trial_types{type}))
        xlabel('Position')
        ylabel('Proportion of first licks')
        plot_count = plot_count + 1;
        set(gca,'fontsize',14)
    end
    legend('Left licks','Right licks','Color','none')
    sgtitle(sprintf('Track %i lick histogram',track_id))
end




%% Percentage trials with active reward
fig = figure
fig.Position = [300 150 480 360];
fig.Name = sprintf('Proportion of active reward trials %s %s',mouse,date);


prop_active(1) = length(intersect(find(Behaviour.track_ID_all == 1),find(reward_type > 2)))...
    /length(find(Behaviour.track_ID_all == 1));

prop_active(2) = length(intersect(find(Behaviour.track_ID_all == 2),find(reward_type > 2)))...
    /length(find(Behaviour.track_ID_all == 2));


bar(1,prop_active(1),'FaceColor','r','FaceAlpha',0.5); hold on
bar(2,prop_active(2),'FaceColor','b','FaceAlpha',0.5)

ylim([0 1])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Track 1','Track 2'})
ylabel('proportion of trials with active reward')
hold on;
set(gca,"TickDir","out",'box', 'off','Color','none')
set(gca,'fontsize',14)


%% Lick bias
fig = figure
fig.Position = [300 150 1250 860];
fig.Name = sprintf('Lick bias ratio across two tracks %s %s',mouse,date);
%         fig.Position = [580 150 1100 830];
%         fig.Name = sprintf('Track %i lick bias %s %s',track_id,session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);
Trial_types = {'All','Passive','Active','Active only'};
count = 1;
for type = 1:3

    if type == 1
        % All trials
        this_type_laps = 1:length(Behaviour.track_ID_all);
    elseif type == 2
        % Find all passive trials
        this_type_laps = find(reward_type <= 2);

    elseif type == 3
        % Find all active trials
        this_type_laps = find(reward_type > 2);

    elseif type == 4
        % Find active only trials
        this_type_laps = find(Behaviour.trial_type == 1);
    end

%     this_type_laps = this_type_laps(this_type_laps < 100);
    track1_all_laps = find(Behaviour.track_ID_all == 1);
    track2_all_laps = find(Behaviour.track_ID_all == 2);
    [~,track1_laps,~]= intersect(track1_all_laps,this_type_laps);
    [~,track2_laps,~]= intersect(track2_all_laps,this_type_laps);

%     subplot(4,2,count)
%     bar(1,mean(lick_ratio{1}(track1_laps)))
%     hold on
%     bar(2,mean(lick_ratio{2}(track2_laps)))
%     errorbar([mean(lick_ratio{1}(track1_laps)),mean(lick_ratio{2}(track2_laps))],...
%         [std(lick_ratio{1}(track1_laps))/sqrt(length(track1_laps)),std(lick_ratio{2}(track2_laps))/sqrt(length(track2_laps))],...
%         '.','Color','k')
%     ylim([0 1.25*max([mean(lick_ratio{1}(track1_laps)) mean(lick_ratio{2}(track2_laps))])])
%     xticks([1 2])
%     xticklabels({'T1 lick bias','T2 lick bias'})
%     set(gca,"TickDir","out",'box', 'off','Color','none')
%     set(gca,'fontsize',14)
%     title(sprintf('%s Trials (anticipatory licks)',Trial_types{type}))
%     count = count + 1;
% 
% 
%     subplot(4,2,count)
%     bar(1,mean(lick_ratio_all{1}(track1_laps)))
%     hold on
%     bar(2,mean(lick_ratio_all{2}(track2_laps)))
%     errorbar([mean(lick_ratio_all{1}(track1_laps)),mean(lick_ratio_all{2}(track2_laps))],...
%         [std(lick_ratio_all{1}(track1_laps))/sqrt(length(track1_laps)),std(lick_ratio_all{2}(track2_laps))/sqrt(length(track2_laps))],...
%         '.','Color','k')
%     ylim([0 1.25*max([mean(lick_ratio_all{1}(track1_laps)) mean(lick_ratio_all{2}(track2_laps))])])
%     xticks([1 2])
%     xticklabels({'T1 bias','T2 bias'})
%     set(gca,"TickDir","out",'box', 'off','Color','none')
%     set(gca,'fontsize',14)
%     title(sprintf('%s Trials (all licks)',Trial_types{type}))
%     count = count + 1;


    subplot(4,3,count)
    histogram(lick_ratio{1}(track1_laps),-1:0.1:1,'FaceColor',patch_color{1}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
    hold on
    histogram(lick_ratio{2}(track2_laps),-1:0.1:1,'FaceColor',patch_color{2}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
    xlabel('Normalized lick bias')
    ylabel('Proportion of trials')
    ylim([0 1])
    set(gca,"TickDir","out",'box', 'off','Color','none')
    set(gca,'fontsize',12)
    title(sprintf('%s Trials (anticipatory licks)',Trial_types{type}))
    count = count + 1;

    subplot(4,3,count)
    histogram(lick_ratio_reward_zone{1}(track1_laps),-1:0.1:1,'FaceColor',patch_color{1}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
    hold on
    histogram(lick_ratio_reward_zone{2}(track2_laps),-1:0.1:1,'FaceColor',patch_color{2}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
    xlabel('Normalized lick bias')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none')
    set(gca,'fontsize',12)
    ylim([0 1])
    title(sprintf('%s Trials (reward first licks)',Trial_types{type}))
    count = count + 1;

    subplot(4,3,count)
    histogram(lick_ratio_all{1}(track1_laps),-1:0.1:1,'FaceColor',patch_color{1}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
    hold on
    histogram(lick_ratio_all{2}(track2_laps),-1:0.1:1,'FaceColor',patch_color{2}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
    xlabel('Normalized lick bias')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none')
    set(gca,'fontsize',12)
    ylim([0 1])
    title(sprintf('%s Trials (all licks)',Trial_types{type}))
    count = count + 1;


    %             beeswarm(ones(1,length(lick_ratio{1}(track1_laps)))',lick_ratio{1}(track1_laps)','MarkerFaceColor','r');
    %             hold on
    %             beeswarm(2*ones(1,length(lick_ratio{2}(track2_laps)))',lick_ratio{2}(track2_laps)','MarkerFaceColor','b');
    %             xlim([0.5 2.5])
    %             ylim([0 1.2])


    %             beeswarm(ones(1,length(lick_ratio_all{1}(track1_laps)))',lick_ratio_all{1}(track1_laps)','MarkerFaceColor','r');
    %             hold on
    %             beeswarm(2*ones(1,length(lick_ratio_all{2}(track2_laps)))',lick_ratio_all{2}(track2_laps)','MarkerFaceColor','b');
    %             xlim([0.5 2.5])

end

save_all_figures(path,[])







%% Every 10 laps


Trial_types = {'All','Passive','Active'};
for track_id = 1:2

    lickBin = 0:5:140;

    for type = 1:1

        if type == 1
            % All trials
            this_track_laps = find(Behaviour.track_ID_all == track_id);
        elseif type == 2
            % Find all passive trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type <= 2));

        elseif type == 3
            % Find all active trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type > 2));
        end

        %         this_track_laps = find(Behaviour.track_ID_all == track_id);
        % track_2_laps = find(Behaviour.track_ID_all == 2);
        lap_change = [1; diff(this_track_laps)];
        trial_blocks = [[this_track_laps(1); this_track_laps(find(lap_change>1))]';...
            [this_track_laps(find(lap_change>1)-1); this_track_laps(end)]']';

        figure;
        plot_count = 1;
        for nblock = 1:size(trial_blocks,1)

            this_block_laps = this_track_laps(this_track_laps>=trial_blocks(nblock,1)...
                &this_track_laps<=trial_blocks(nblock,2));

            if isempty(this_block_laps)
                continue
            end

            subplot(4,ceil(size(trial_blocks,1)/4),plot_count)
            if ~isempty(cat(1,track1_lick_position{this_block_laps}))
                licks_histcount = histcounts(cat(1,track1_lick_position{this_block_laps}),lickBin);
                %     histogram(cat(1,track1_lick_position{track_1_laps}),lickBin,'FaceAlpha',0.5,'FaceColor',patch_color{1}(5,:),'Normalization','probability')
                bar(lickBin(2:end),licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{1}(5,:))
                ymax= max(licks_histcount);
            end
            hold on

            if ~isempty(cat(1,track2_lick_position{this_block_laps}))
                licks_histcount = histcounts(cat(1,track2_lick_position{this_block_laps}),lickBin);
                bar(lickBin(2:end),-licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{2}(5,:))
                %     histogram(cat(1,track2_lick_position{track_1_laps}),lickBin,'FaceAlpha',0.5,'FaceColor',patch_color{2}(5,:),'Normalization','probability')
                hold on

                ymax= max([ymax max(licks_histcount)]);
            end

            landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
            landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
            landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
            landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
            landmarkVertex(4,:,1) = [86 94 94 86]; %4th
            landmarkVertex(5,:,1) = [106 114 114 106]; %5th
            landmarkVertex(:,:,2) = repmat([-1.5*ymax -1.5*ymax 1.5*ymax 1.5*ymax],[5 1]); %landmark Y coord

            xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',5)
            for iPatch = 1:5
                patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',.3,'EdgeAlpha',0)
            end
            %             legend('Left licks','Right licks')
            xlim([0,140])
            title(sprintf('Lap %i to %i Track %i licks (%s)',trial_blocks(nblock,1),trial_blocks(nblock,2),track_id,Trial_types{type}))
            plot_count = plot_count + 1;
        end
        sgtitle(sprintf('Track %i licks (%s)',track_id,Trial_types{type}))


        track1_first_click = [];
        track2_first_click = [];
        for nlap = 1:length(this_track_laps)
            if isempty(track1_lick_position{this_track_laps(nlap)})
                track1_first_click(nlap) = nan;
            else
                track1_first_click(nlap) = track1_lick_position{this_track_laps(nlap)}(1);
            end

            if isempty(track2_lick_position{this_track_laps(nlap)})
                track2_first_click(nlap) = nan;
            else
                track2_first_click(nlap) = track2_lick_position{this_track_laps(nlap)}(1);
            end
            if isnan(track1_first_click(nlap)) & isnan(track2_first_click(nlap))
                first_lick(nlap) = nan;
            elseif ~isnan(track1_first_click(nlap)) & isnan(track2_first_click(nlap))
                first_lick(nlap) = track1_first_click(nlap);
            elseif isnan(track1_first_click(nlap)) & ~isnan(track2_first_click(nlap))
                first_lick(nlap) = track2_first_click(nlap);
            elseif track1_first_click(nlap) > track2_first_click(nlap)
                first_lick(nlap) = track2_first_click(nlap);
            elseif track1_first_click(nlap) < track2_first_click(nlap)
                first_lick(nlap) = track1_first_click(nlap);
            else
                first_lick(nlap) = nan;
            end
        end

        figure;
        plot_count = 1;
        for nblock = 1:size(trial_blocks,1)

            this_block_laps = this_track_laps>=trial_blocks(nblock,1)...
                &this_track_laps<=trial_blocks(nblock,2);

            subplot(4,ceil(size(trial_blocks,1)/4),plot_count)
            if sum(~isnan(track1_first_click(this_block_laps))) > 0
                %     histogram(track1_first_click,lickBin,'FaceAlpha',0.3,'FaceColor',patch_color{1}(5,:),'Normalization','probability')
                licks_histcount = histcounts(track1_first_click(this_block_laps),lickBin);
                bar(lickBin(2:end),licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{1}(5,:))
                ymax= max(licks_histcount);
            end

            hold on
            if sum(~isnan(track2_first_click(this_block_laps))) > 0
                %     histogram(track2_first_click,lickBin,'FaceAlpha',0.3,'FaceColor',patch_color{2}(5,:),'Normalization','probability')
                licks_histcount = histcounts(track2_first_click(this_block_laps),lickBin);
                bar(lickBin(2:end),-licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{2}(5,:))
                ymax= max([ymax max(licks_histcount)]);
            end

            xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',5)
            hold on
            landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
            landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
            landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
            landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
            landmarkVertex(4,:,1) = [86 94 94 86]; %4th
            landmarkVertex(5,:,1) = [106 114 114 106]; %5th
            landmarkVertex(:,:,2) = repmat([-1.5*ymax -1.5*ymax 1.5*ymax 1.5*ymax],[5 1]); %landmark Y coord

            for iPatch = 1:5
                patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',.5,'EdgeAlpha',0)
            end
            % legend('Left licks','Right licks')
            xlim([0,140])
            title(sprintf('Lap %i to %i Track %i first licks (%s)',trial_blocks(nblock,1),trial_blocks(nblock,2),track_id,Trial_types{type}))
            plot_count = plot_count + 1;
        end
        sgtitle(sprintf('Track %i first licks (%s)',track_id,Trial_types{type}))
    end
end

%% Speed profile




%% Every 10 laps


Trial_types = {'All','Passive','Active'};
for track_id = 1:2

    lickBin = 0:5:140;

    for type = 1:1

        if type == 1
            % All trials
            this_track_laps = find(Behaviour.track_ID_all == track_id);
        elseif type == 2
            % Find all passive trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type <= 2));

        elseif type == 3
            % Find all active trials
            this_track_laps = intersect(find(Behaviour.track_ID_all == track_id),find(reward_type > 2));
        end

%         this_track_laps = find(Behaviour.track_ID_all == track_id);
        % track_2_laps = find(Behaviour.track_ID_all == 2);
        lap_change = [1; diff(this_track_laps)];
        trial_blocks = [[this_track_laps(1); this_track_laps(find(lap_change>1))]';...
            [this_track_laps(find(lap_change>1)-1); this_track_laps(end)]']';
        
        figure;
        plot_count = 1;
        for nblock = 1:size(trial_blocks,1)
            
            this_block_laps = this_track_laps(this_track_laps>=trial_blocks(nblock,1)...
                &this_track_laps<=trial_blocks(nblock,2));

            if isempty(this_block_laps)
                continue
            end
            
            lap_times
            subplot(4,ceil(size(trial_blocks,1)/4),plot_count)
            plot()
                ymax= max([ymax max(licks_histcount)]);


            landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
            landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
            landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
            landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
            landmarkVertex(4,:,1) = [86 94 94 86]; %4th
            landmarkVertex(5,:,1) = [106 114 114 106]; %5th
            landmarkVertex(:,:,2) = repmat([-1.5*ymax -1.5*ymax 1.5*ymax 1.5*ymax],[5 1]); %landmark Y coord

            xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',5)
            for iPatch = 1:5
                patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',.3,'EdgeAlpha',0)
            end
%             legend('Left licks','Right licks')
            xlim([0,140])
            title(sprintf('Lap %i to %i Track %i licks (%s)',trial_blocks(nblock,1),trial_blocks(nblock,2),track_id,Trial_types{type}))
            plot_count = plot_count + 1;
        end
         sgtitle(sprintf('Track %i licks (%s)',track_id,Trial_types{type}))


        track1_first_click = [];
        track2_first_click = [];
        for nlap = 1:length(this_track_laps)
            if isempty(track1_lick_position{this_track_laps(nlap)})
                track1_first_click(nlap) = nan;
            else
                track1_first_click(nlap) = track1_lick_position{this_track_laps(nlap)}(1);
            end

            if isempty(track2_lick_position{this_track_laps(nlap)})
                track2_first_click(nlap) = nan;
            else
                track2_first_click(nlap) = track2_lick_position{this_track_laps(nlap)}(1);
            end
            if isnan(track1_first_click(nlap)) & isnan(track2_first_click(nlap))
                first_lick(nlap) = nan;
            elseif ~isnan(track1_first_click(nlap)) & isnan(track2_first_click(nlap))
                first_lick(nlap) = track1_first_click(nlap);
            elseif isnan(track1_first_click(nlap)) & ~isnan(track2_first_click(nlap))
                first_lick(nlap) = track2_first_click(nlap);
            elseif track1_first_click(nlap) > track2_first_click(nlap)
                first_lick(nlap) = track2_first_click(nlap);
            elseif track1_first_click(nlap) < track2_first_click(nlap)
                first_lick(nlap) = track1_first_click(nlap);
            else
                first_lick(nlap) = nan;
            end
        end
        
        figure;
        plot_count = 1;
        for nblock = 1:size(trial_blocks,1)

            this_block_laps = this_track_laps>=trial_blocks(nblock,1)...
                &this_track_laps<=trial_blocks(nblock,2);

        subplot(4,ceil(size(trial_blocks,1)/4),plot_count)
            if sum(~isnan(track1_first_click(this_block_laps))) > 0
                %     histogram(track1_first_click,lickBin,'FaceAlpha',0.3,'FaceColor',patch_color{1}(5,:),'Normalization','probability')
                licks_histcount = histcounts(track1_first_click(this_block_laps),lickBin);
                bar(lickBin(2:end),licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{1}(5,:))
                ymax= max(licks_histcount);
            end

            hold on
            if sum(~isnan(track2_first_click(this_block_laps))) > 0
                %     histogram(track2_first_click,lickBin,'FaceAlpha',0.3,'FaceColor',patch_color{2}(5,:),'Normalization','probability')
                licks_histcount = histcounts(track2_first_click(this_block_laps),lickBin);
                bar(lickBin(2:end),-licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{2}(5,:))
                ymax= max([ymax max(licks_histcount)]);
            end

            xline(100,'Color',patch_color{track_id}(5,:),'LineWidth',5)
            hold on
            landmarkVertex=zeros(5,4,2); %vertex coordinate for patch of landmarks on the plot, third dimension is x and y.
            landmarkVertex(1,:,1) = [26 34 34 26]; %1st landmark X coord
            landmarkVertex(2,:,1) = [46 54 54 46]; %2nd
            landmarkVertex(3,:,1) = [66 74 74 66]; %3rd
            landmarkVertex(4,:,1) = [86 94 94 86]; %4th
            landmarkVertex(5,:,1) = [106 114 114 106]; %5th
            landmarkVertex(:,:,2) = repmat([-1.5*ymax -1.5*ymax 1.5*ymax 1.5*ymax],[5 1]); %landmark Y coord

            for iPatch = 1:5
                patch(landmarkVertex(iPatch,:,1),landmarkVertex(iPatch,:,2),patch_color{track_id}(iPatch,:),'FaceAlpha',.5,'EdgeAlpha',0)
            end
            % legend('Left licks','Right licks')
            xlim([0,140])
            title(sprintf('Lap %i to %i Track %i first licks (%s)',trial_blocks(nblock,1),trial_blocks(nblock,2),track_id,Trial_types{type}))
            plot_count = plot_count + 1;
        end
        sgtitle(sprintf('Track %i first licks (%s)',track_id,Trial_types{type}))
    end
end
%% 
for track_id = 1:2
    subplot(2,2,track_id)
    for nlap = 1:length(lap_times(track_id).lap)
        this_lap_speed = lap_times(track_id).lap(nlap).v_cm;
        this_lap_speed(this_lap_speed > 100 | this_lap_speed < -100) = 0;
        if lap_times(track_id).lap(nlap).reward_type>= 10 % active trials
            plot(lap_times(track_id).lap(nlap).x,this_lap_speed,'r')
        else
            plot(lap_times(track_id).lap(nlap).x,this_lap_speed,'k')
        end
        hold on
    end
    position_bin = 0:10:140;
    this_track_all_positions = [lap_times(track_id).lap(nlap).x];
    = histcounts(this_track_all_positions,position_bin)

    this_track_all_speeds = [lap_times(track_id).lap(:).v_cm];
    this_track_all_speeds(this_track_all_speeds > 100 | this_track_all_speeds < -100) = 0;

    licks_histcount = histcounts(track2_first_click(this_block_laps),lickBin);
    bar(lickBin(2:end),-licks_histcount,'FaceAlpha',0.9,'FaceColor',patch_color{2}(5,:))
    ymax= max([ymax max(licks_histcount)]);

    title(sprintf('Track %i',track_id))
end
