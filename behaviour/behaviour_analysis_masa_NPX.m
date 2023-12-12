
clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';



for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(session_info)
        continue
    end

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)

        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load extracted_laps
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        column = 1;
        %% Saving lap info


        for nlap = 1:length(Behaviour.lap_ID_all)
            current_track = Behaviour.track_ID_all(nlap);
            current_lap = Behaviour.lap_ID_all(nlap);
            reward_lap = find(Behaviour.(sprintf('track%i_count',current_track)) == current_lap);

            if current_lap == 0
                continue
            end

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
            lap_times(current_track).lap(current_lap).trial_type = Behaviour.trial_type(nlap); % if 1 means active only trial, 2 means hybrid
            lap_times(current_track).lap(current_lap).lap_ID_all = Behaviour.lap_ID_all(nlap); % if 1 means active only trial, 2 means hybrid

            current_lap_licks = find(Behaviour.(sprintf('lick_track%i_count',current_track)) == current_lap & ...
                Behaviour.lick_track_ID == current_track);

            track1_licks = find(Behaviour.lick_state(current_lap_licks) == 1); % track1 means left lick
            track2_licks = Behaviour.lick_state(current_lap_licks) == 2; % track2 means right lick

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

                    end

                    if type == 1
                        lap_times(track_id).lap(nlap).lick_ratio =  lick_ratio{track_id}(nlap);
                        lap_times(track_id).lap(nlap).lick_ratio_all =  lick_ratio_all{track_id}(nlap);
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
        fig.Name = sprintf('Lick behaviour summary %s %s',session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);

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
            fig.Name = sprintf('Track %i lick bias %s %s',track_id,session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);
            Trial_types = {'All','Passive reward','Active reward','Active only'};

            for type = 1:4

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



        fig = figure
        fig.Position = [300 150 600 600];
        fig.Name = sprintf('Lap running speed vs position %s %s',session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);
        
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
                lap_mean_speed(nlap) = median(lap_times(track_id).lap(nlap).v_cm((lap_times(track_id).lap(nlap).x > 0)));
                [N,edges,x_bins] = histcounts(lap_times(track_id).lap(nlap).x,x_bin_edges);
                for nbin = 1:max(x_bins)
                    speed_position(nlap,nbin) = median(lap_times(track_id).lap(nlap).v_cm(x_bins == nbin));
                end
            end
            %             [~,sorted_lap_id{track_id}] = sort(lap_mean_speed);
            low_speed_laps = find(lap_mean_speed <= median(lap_mean_speed));
            high_speed_laps = find(lap_mean_speed > median(lap_mean_speed));

            speed_LSE = nanmean(speed_position(low_speed_laps,:))- nanstd(speed_position(low_speed_laps,:))/sqrt(length(lap_times(track_id).lap));
            speed_USE = nanmean(speed_position(low_speed_laps,:))+ nanstd(speed_position(low_speed_laps,:))/sqrt(length(lap_times(track_id).lap));

            %             plot(x, speed_LSE, 'k--', 'LineWidth', 1);hold on;
            %             plot(x, speed_USE, 'k--', 'LineWidth', 1);
                        hold on
            plot(x, nanmean(speed_position(low_speed_laps,:)), 'Color',patch_color{track_id}(4,:), 'LineWidth', 1);hold on
            x2 = [x, fliplr(x)];
            inBetween = [speed_LSE, fliplr(speed_USE)];
            h(1) = fill(x2, inBetween, patch_color{track_id}(4,:),'FaceAlpha',0.2,'EdgeColor','none');


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
        fig.Name = sprintf('Mean Lick bias ratio %s %s',session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);
%         fig.Position = [580 150 1100 830];
%         fig.Name = sprintf('Track %i lick bias %s %s',track_id,session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);
        Trial_types = {'All','Passive','Active','Active only'};

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

            track1_all_laps = find(Behaviour.track_ID_all == 1);
            track2_all_laps = find(Behaviour.track_ID_all == 2);
            [~,track1_laps,~]= intersect(track1_all_laps,this_type_laps);
            [~,track2_laps,~]= intersect(track2_all_laps,this_type_laps);

            subplot(2,2,type)
            bar(1,mean(lick_ratio{1}(track1_laps)))
            hold on
            bar(2,mean(lick_ratio{2}(track2_laps)))
            errorbar([mean(lick_ratio{1}(track1_laps)),mean(lick_ratio{2}(track2_laps))],...
                [std(lick_ratio{1}(track1_laps))/sqrt(length(track1_laps)),std(lick_ratio{2}(track2_laps))/sqrt(length(track2_laps))],...
                '.','Color','k')
            if 1.25*max([mean(lick_ratio{1}(track1_laps)) mean(lick_ratio{2}(track2_laps))]) > 0
                ylim([0 1.25*max([mean(lick_ratio{1}(track1_laps)) mean(lick_ratio{2}(track2_laps))])])
            end


            xticks([1 2])
            xticklabels({'T1 lick bias','T2 lick bias'})
            set(gca,"TickDir","out",'box', 'off','Color','none')
            set(gca,'fontsize',14)
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

        %% Percentage trials with active reward
        fig = figure
        fig.Position = [300 150 480 360];
        fig.Name = sprintf('Proportion of active reward trials %s %s',session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);


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
        fig.Position = [300 150 1100 830];
        %         fig.Name = sprintf('Lick bias ratio across two tracks %s %s',mouse,date);
        %         fig.Position = [580 150 1100 830];
        fig.Name = sprintf('lick bias ratio across two tracks %s %s',session_info(n).probe(1).SUBJECT,session_info(n).probe(1).SESSION);
        Trial_types = {'All','Passive','Active','Active only'};
        count = 1;
        for type = 1:4

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
            ylabel('Proportion of trials')
            ylim([0 1])
            set(gca,"TickDir","out",'box', 'off','Color','none')
            set(gca,'fontsize',12)
            title(sprintf('%s Trials (reward first licks)',Trial_types{type}))
            count = count + 1;

            subplot(4,3,count)
            histogram(lick_ratio_all{1}(track1_laps),-1:0.1:1,'FaceColor',patch_color{1}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
            hold on
            histogram(lick_ratio_all{2}(track2_laps),-1:0.1:1,'FaceColor',patch_color{2}(5,:),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
            xlabel('Normalized lick bias')
            ylabel('Proportion of trials')
            ylim([0 1])
            set(gca,"TickDir","out",'box', 'off','Color','none')
            set(gca,'fontsize',12)
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
        %         save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project\behaviour',[])
        mkdir('behaviour')
        save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis','behaviour'),[])
       
        
    end

end