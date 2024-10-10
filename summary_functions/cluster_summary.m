function cluster_summary(clusters_all,cluster_ids,method,spatial_response,spatial_response_extended,varargin)
%currently just feed the place_fields containing all the clusters from all
%sessions

% purpose of this is just to input cluster_ids that you want to make a
% summary from and it should easily find them in place_fields to do it
p = inputParser;
addParameter(p, 'plot_behaviour', false);
addParameter(p, 'break_loop', false);
parse(p, varargin{:});


plot_behaviour = p.Results.plot_behaviour;
break_loop = p.Results.break_loop;
% first get the index of the cluster_ids
ids_index = ismember(clusters_all.cluster_id,cluster_ids);
session_count = clusters_all.session_count(ids_index);
if nargin < 4
    unique_session_count = unique(session_count);
    delay = 0;
    spatial_response = cell(sum(ids_index),2);
    spatial_response_extended = cell(sum(ids_index),2);
    cluster_counter = 0;
    for iS = unique_session_count
        clusters_this_session = clusters_all.session_count == iS & ids_index;
        no_clusters_this_session = sum(clusters_this_session);

        tvec = clusters_all.tvec{iS};
        start_time_all = clusters_all.start_time_all{iS};
        end_time_all = clusters_all.end_time_all{iS};
        track_ID_all = clusters_all.track_ID_all{iS};
        position = clusters_all.position{iS};
        speed = clusters_all.speed{iS};
        spike_times_this_session = clusters_all.spike_times(clusters_this_session);
        parfor iC = 1:no_clusters_this_session
            spike_times = spike_times_this_session{iC};
            spatial_response(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
                position,speed,track_ID_all,start_time_all,end_time_all,'within',delay);
            spatial_response_extended(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
                position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay);
        end
        cluster_counter = cluster_counter+no_clusters_this_session;
    end
end


no_clusters = length(cluster_ids);
iC = 1;
% Get the screen size
screenSize = get(0, 'ScreenSize');

% Calculate the position and size of the figure
figureWidth = screenSize(3) * 0.9;
figureHeight = screenSize(4) * 0.9;
figureX = (screenSize(3) - figureWidth) / 2;
figureY = (screenSize(4) - figureHeight) / 2;
%% Create the figure with the specified position and size
fig1 = figure('Position', [figureX figureY figureWidth figureHeight]);
if plot_behaviour
    fig2 = figure('Position', [figureX figureY figureWidth figureHeight]);
end
prev_session_count = [];
while true

    %this is where you get all the info for that cluster at a time
    % so basically this is like a object-based programming
    
    %task information of this session that this cluster is from
    current_cluster = cluster_spatial_info(clusters_all, cluster_ids(iC),...
        spatial_response(iC,:), spatial_response_extended(iC,:));

    if iC == 1
        figure(fig1);
        update_figure1(current_cluster,method);
        if plot_behaviour
            figure(fig2);
            update_figure2(current_cluster);
        end
    end

    figure(fig1);
    update_figure1(current_cluster,method);
    if plot_behaviour
        % Update the second figure if the session count has changed
        if isempty(prev_session_count) || current_cluster.session_count ~= prev_session_count
            figure(fig2);
            update_figure2(current_cluster);  % replace with your actual code for updating the figure
        end
    end
    if break_loop
        break;
    end
    % Update the previous session count
    prev_session_count = current_cluster.session_count;
    % Wait for a key press
    waitforbuttonpress;
    key = get(gcf, 'CurrentCharacter');

    % If right arrow key is pressed
    if double(key) == 29 && iC < no_clusters
        iC = iC + 1;
        disp(iC)
        % If left arrow key is pressed
    elseif double(key) == 28 && iC > 1
        iC = iC - 1;
        disp(iC)
        % If 'q' is pressed
    elseif lower(key) == 'q'
        break;
    end




end
    function update_figure1(c,method)
        %% Plot subplots with inbuilt plotting functions - all start with cluster_xxx()
        clf;
        %title gives information about the cluster properties
        sgtitle(['Cluster ID: ' num2str(c.id), 'Depth ', num2str(c.peak_depth),...
            'Region',c.region], 'FontSize', 8);

        %   T1: Context(X) - A - B - A - B
        %   T2: Context(Y) - C - B - C - B - Context(E)
        % Colors (RGB):
        %       A: blue series: [146,197,222]/255 [67,147,195]/255
        %       B: red series: [214,96,77]/255 [244,165,130]/255 (4 of them so backup [178,24,43]/255 [253,219,199]/255
        %       C: green series: [166,219,160]/255 [90,174,97]/255
        %       X: purple series: [194,165,207]/255 [153,112,171]/255
        %       Y: orange series: [255,165,0]/255 [255,140,0]/255

        switch method
            case 'within'
                ax1 = subplot(8,6,[1:2,7:8]);
                cluster_plot_lap(c.spatial_response{1}(c.rewarded_laps_t1,:)); %use normalized spatial response in lap plots
                colormap(ax1, flip(gray(256)));
                title('t1')
                yl = ylim;
                track_label('t1',yl(1));
                ax2 = subplot(8,6,[13:14,19:20]);
                if ~isempty(c.t2_index)
                    cluster_plot_lap(c.spatial_response{2}(c.rewarded_laps_t2,:));
                    colormap(ax2, flip(gray(256)));
                    title('t2')
                    yl = ylim;
                    track_label('t2',yl(1));
                    subplot(8,6,[25:26,31:32]);
                    cluster_plot_both_psth(c.spatial_response{1}(c.rewarded_laps_t1,:),c.spatial_response{2}(c.rewarded_laps_t2,:));%use smoothed response in psth to represent the firing rate of this cluster
                    yl = ylim;
                    track_label('both',yl(2));
                end

            case 'extended'
                subplot(4,3,1)
                cluster_plot_lap_extended(c.spatial_response_extended{1},c.track_ID_all(c.rewarded_laps_t1))
                subplot(4,3,4)
                cluster_plot_lap_extended(c.spatial_response_extended{2},c.track_ID_all(c.rewarded_laps_t2))
                subplot(4,3,7)
                cluster_plot_both_extended_psth(c.spatial_response_extended{1},c.spatial_response_extended{2})
        end
        ax_cl1 = subplot(8,6,[3,9]);
        landmark1_bin = [30:70; 70:110];
        h = comapre_landmark_response(c.spatial_response{1}(c.rewarded_laps_t1,:),landmark1_bin);
        title(['SMI_t1_l1=',num2str(c.SMI{1}(1,1))])
        set(ax_cl1.XAxis, 'visible', 'off');
        xline([0],'LineWidth',1.5,'Color',[146,197,222]/255,'HandleVisibility','off')
        yl = ylim;
        text(0,yl(1),'A','HorizontalAlignment','center','VerticalAlignment','bottom')
        lgd = legend([h(1:4)],{'Odd r1','Even r1','Odd r2','Even r2'}, 'Location', 'best');
        lgd.FontSize = 8;
        ax_cl2 = subplot(8,6,[4,10]);
        landmark2_bin = [50:90; 90:130];
        comapre_landmark_response(c.spatial_response{1}(c.rewarded_laps_t1,:),landmark2_bin);
        title(['SMI_t1_l2=',num2str(c.SMI{1}(1,2))])
        set(ax_cl2.XAxis, 'visible', 'off');
        xline([0],'LineWidth',1.5,'Color',[214,96,77]/255,'HandleVisibility','off')
        yl = ylim;
        text(0,yl(2),'B','HorizontalAlignment','center','VerticalAlignment','bottom')
        ax_cl3 = subplot(8,6,[15,21]);
        if ~isempty(c.t2_index)
            landmark1_bin = [30:70; 70:110];
            comapre_landmark_response(c.spatial_response{2}(c.rewarded_laps_t2,:),landmark1_bin);
            title(['SMI_t2_l1=',num2str(c.SMI{2}(1,1))])
            set(ax_cl3.XAxis, 'visible', 'off');
            xline([0],'LineWidth',1.5,'Color',[166,219,160]/255,'HandleVisibility','off')
            yl = ylim;
            text(0,yl(2),'C','HorizontalAlignment','center','VerticalAlignment','bottom')
            ax_cl4 = subplot(8,6,[16,22]);
            landmark2_bin = [50:90; 90:130];
            comapre_landmark_response(c.spatial_response{2}(c.rewarded_laps_t2,:),landmark2_bin);
            title(['SMI_t2_l2=',num2str(c.SMI{2}(1,2))])
            set(ax_cl4.XAxis, 'visible', 'off');
            xline([0],'LineWidth',1.5,'Color',[214,96,77]/255,'HandleVisibility','off')
            yl = ylim;
            text(0,yl(2),'B','HorizontalAlignment','center','VerticalAlignment','bottom')

            subplot(8,6,[27,33])
            landmark_contextual_bin = [10:50; 110:150];
            spatial_response_extended_t2_contextual = zeros(length(c.spatial_response_extended{2}),160);
            for iLap = 1:length(c.spatial_response_extended{2})
                if length(c.spatial_response_extended{2}{iLap}) >= 160
                    spatial_response_extended_t2_contextual(iLap,:) = c.spatial_response_extended{2}{iLap}(1:160);
                else
                    spatial_response_extended_t2_contextual(iLap,1:length(c.spatial_response_extended{2}{iLap})) = c.spatial_response_extended{2}{iLap}(1:end);
                end
            end
            comapre_landmark_response(spatial_response_extended_t2_contextual(c.rewarded_laps_t2,:),landmark_contextual_bin);
            title(['SMI_t2_context=',num2str(c.SMI{2}(1,3))])
            xline([0],'LineWidth',1.5,'Color',[255,165,0]/255,'HandleVisibility','off')
            yl = ylim;
            text(0,yl(2),'Y','HorizontalAlignment','center','VerticalAlignment','bottom')
            subplot(8,6,[28,34])
            compare_common_landmark_between_t1_and_t2(c.spatial_response{1}(c.rewarded_laps_t1,:),c.spatial_response{2}(c.rewarded_laps_t2,:));
            xline([0],'LineWidth',1.5,'Color',[214,96,77]/255,'HandleVisibility','off');
            yl = ylim;
            text(0,yl(2),'B','HorizontalAlignment','center','VerticalAlignment','bottom')
        end
        %subplot(4,3,5)
        %cluster_even_odd_psth(c.spatial_response{2})
        %title('T2')
        speed_profile = behaviour_profile_calculator(c.tvec,c.speed,c.position,c.start_time_all,c.end_time_all);
        t1_index = c.rewarded_laps_index & c.track_ID_all == 1;
        t2_index = c.rewarded_laps_index & c.track_ID_all == 2;
        [~,t1_speed_sort] = sort(mean(speed_profile(t1_index,:),2,'omitnan'));
        ax5 = subplot(8,6,[5:6,11:12]);
        %active_reward_laps_sort_t1 = [find(c.active_reward_laps_t1);find(~c.active_reward_laps_t1)];
        cluster_plot_lap(c.spatial_response{1}(c.rewarded_laps_t1,:),t1_speed_sort);
        colormap(ax5, flip(gray(256)));
        title('t1')
        yl = ylim;
        track_label('t1',yl(1));
        %yline(sum(c.active_reward_laps_t1),'Color','r','LineWidth',1.5)
        %text(sum(c.active_reward_laps_t1),0,'A','HorizontalAlignment','left','VerticalAlignment','middle')
        %cluster_block_psth(c.spatial_response{1},c.block_index_t1) % psth of each block
        ax6 = subplot(8,6,[17:18,23:24]);
        if ~isempty(c.t2_index)
            %active_reward_laps_sort_t2 = [find(c.active_reward_laps_t2);find(~c.active_reward_laps_t2)];
            [~,t2_speed_sort] = sort(mean(speed_profile(t2_index,:),2,'omitnan'));
            cluster_plot_lap(c.spatial_response{2}(c.rewarded_laps_t2,:),t2_speed_sort);
            colormap(ax6, flip(gray(256)));
            title('t2')
            yl = ylim;
            track_label('t2',yl(1));
            %yline(sum(c.active_reward_laps_t2),'Color','r','LineWidth',1.5)
            %text(sum(c.active_reward_laps_t2),0,'C','HorizontalAlignment','left','VerticalAlignment','middle')
        end
        %cluster_block_psth(c.spatial_response{2},c.block_index_t2)

        subplot(8,6,[29:30,35:36])
        h = cluster_plot_psth(c.spatial_response{1}(c.rewarded_laps_t1,:),c.active_reward_laps_t1);
        title('t1')
        yl = ylim;
        track_label('t1',yl(2));
        legend([h(1:2)],{['active reward',num2str(sum(c.active_reward_laps_t1))],['passive',num2str(sum(~c.active_reward_laps_t1))]}, 'Location', 'northwest');
        %cluster_block_dynamic(c.spatial_response{1},c.block_index_t1) % first, second, third, etc laps of each block
        subplot(8,6,[41:42,47:48])
        if ~isempty(c.t2_index)
            h = cluster_plot_psth(c.spatial_response{2}(c.rewarded_laps_t2,:),c.active_reward_laps_t2);
            title('t2')
            yl = ylim;
            track_label('t2',yl(2));
            legend([h(1:2)],{['active reward',num2str(sum(c.active_reward_laps_t2))],['passive',num2str(sum(~c.active_reward_laps_t2))]}, 'Location', 'northwest');
        end
        %cluster_block_dynamic(c.spatial_response_t2,c.block_index_t2)

        subplot(8,6,[37:38,43:44])
        plot(c.waveform,'LineWidth',3);
        hold off;
        title('peak channel waveform')

        subplot(8,6,[39:40,45:46])
        plot_speed_response(c.speed_response);
        title('speed response')
    end

    function update_figure2(c)
        clf;
        speed_profile = behaviour_profile_calculator(c.tvec,c.speed,c.position,c.start_time_all,c.end_time_all);
        if ~isempty(c.face_motion_SVD)
        face_svd1_profile = behaviour_profile_calculator(c.tvec,c.face_motion_SVD(:,1),c.position,c.start_time_all,c.end_time_all);
        face_svd2_profile = behaviour_profile_calculator(c.tvec,c.face_motion_SVD(:,2),c.position,c.start_time_all,c.end_time_all);
        face_svd3_profile = behaviour_profile_calculator(c.tvec,c.face_motion_SVD(:,3),c.position,c.start_time_all,c.end_time_all);
        face_svd4_profile = behaviour_profile_calculator(c.tvec,c.face_motion_SVD(:,4),c.position,c.start_time_all,c.end_time_all);
        face_svd5_profile = behaviour_profile_calculator(c.tvec,c.face_motion_SVD(:,5),c.position,c.start_time_all,c.end_time_all);
        pupil_sze_profile = behaviour_profile_calculator(c.tvec,c.pupil_size,c.position,c.start_time_all,c.end_time_all);
        pupil_angle_profile = behaviour_profile_calculator(c.tvec,c.pupil_angle,c.position,c.start_time_all,c.end_time_all);
        end
        reward_position = c.reward_position;
        [lap_lick_position,lap_reward_position] = lick_reward_position_calculator(c.lick_position,reward_position,c.lick_time,c.track_ID_all,c.start_time_all,c.end_time_all);
        t1_index = c.rewarded_laps_index & c.track_ID_all == 1;
        t2_index = c.rewarded_laps_index & c.track_ID_all == 2;
        subplot(8,6,[1:2,7:8])
        session_behaviour_plot(speed_profile(t1_index,:),lap_lick_position(t1_index),c.reward_type(t1_index),lap_reward_position(t1_index),c.block_index_t1);
        xlim([0,140])
        landmark_label(sum(t1_index),'t1')
        yl = ylim; track_label('t1',yl(2));
        reward_label('t1')

        subplot(8,6,[13:14,19:20])
        lick_psth(lap_lick_position(t1_index),c.active_reward_laps_t1)
        %active vs passive?
        yl = ylim;
        landmark_label(yl(2),'t1')
        reward_label('t1')

        subplot(8,6,[25:26,31:32])
        speed_profile_t1 = speed_profile(t1_index,:);
        speed_psth(speed_profile_t1(~c.active_reward_laps_t1,:),[[67,147,195]/255,0.1])
        yl = ylim; landmark_label(yl(2),'t1'); reward_label('t1');
        subplot(8,6,[37:38,43:44])
        if sum(c.active_reward_laps_t1) > 0
            speed_psth(speed_profile_t1(c.active_reward_laps_t1,:),[[244,165,130]/255,0.1])
            yl = ylim; landmark_label(yl(2),'t1'); reward_label('t1');
        end

        if ~isempty(c.t2_index)
            subplot(8,6,[3:4,9:10])
            session_behaviour_plot(speed_profile(t2_index,:),lap_lick_position(t2_index),c.reward_type(t2_index),lap_reward_position(t2_index),c.block_index_t2);
            xlim([0,140])
            landmark_label(sum(t2_index),'t2')
            yl = ylim; track_label('t2',yl(2));
            reward_label('t2')
            subplot(8,6,[15:16,21:22])
            lick_psth(lap_lick_position(t2_index),c.active_reward_laps_t2)
            yl = ylim;
            reward_label('t2')
            landmark_label(yl(2),'t2')
            subplot(8,6,[27:28,33:34])
            speed_profile_t2 = speed_profile(t2_index,:);
            speed_psth(speed_profile_t2(~c.active_reward_laps_t2,:),[[67,147,195]/255,0.1])
            yl = ylim; landmark_label(yl(2),'t2'); reward_label('t2');
            subplot(8,6,[39:40,45:46])
            if sum(c.active_reward_laps_t2) > 0
                speed_psth(speed_profile_t2(c.active_reward_laps_t2,:),[[253,219,199]/255,0.1])
                yl = ylim; landmark_label(yl(2),'t2'); reward_label('t2');
            end
        end

        if ~isempty(c.face_motion_SVD)
        subplot(8,6,[5:6,11:12])
        hold on;
        t1_face1_motion = face_svd1_profile(t1_index,:);
        t1_face2_motion = face_svd2_profile(t1_index,:);
        t1_face3_motion = face_svd3_profile(t1_index,:);
        t1_face4_motion = face_svd4_profile(t1_index,:);
        t1_face5_motion = face_svd5_profile(t1_index,:);
        color = ([254,217,142;
            161,218,180;
            65,182,196;
            44,127,184;
            37,52,148])./255;
        h(1) = shade_psth(1:140,mean(t1_face1_motion,1),std(t1_face1_motion,[],1)./sqrt(sum(t1_index)),color(1,:));
        h(2) = shade_psth(1:140,mean(t1_face2_motion,1),std(t1_face2_motion,[],1)./sqrt(sum(t1_index)),color(2,:));
        h(3) = shade_psth(1:140,mean(t1_face3_motion,1),std(t1_face3_motion,[],1)./sqrt(sum(t1_index)),color(3,:));
        h(4) = shade_psth(1:140,mean(t1_face4_motion,1),std(t1_face4_motion,[],1)./sqrt(sum(t1_index)),color(4,:));
        h(5) = shade_psth(1:140,mean(t1_face5_motion,1),std(t1_face5_motion,[],1)./sqrt(sum(t1_index)),color(5,:));
        yl = ylim; landmark_label(yl(2),'t1'); reward_label('t1');
        legend([h(1:5)],{'1st','2nd','3rd','4th','5th'}, 'Location', 'northwest');
        title('t1 face motion SVD')


        subplot(8,6,[17:18,23:24])
        if ~isempty(c.t2_index)
            hold on;
            t2_face1_motion = face_svd1_profile(t2_index,:);
            t2_face2_motion = face_svd2_profile(t2_index,:);
            t2_face3_motion = face_svd3_profile(t2_index,:);
            t2_face4_motion = face_svd4_profile(t2_index,:);
            t2_face5_motion = face_svd5_profile(t2_index,:);
            color = ([254,217,142;
                161,218,180;
                65,182,196;
                44,127,184;
                37,52,148])./255;
            h(1) = shade_psth(1:140,mean(t2_face1_motion,1),std(t2_face1_motion,[],1)./sqrt(sum(t2_index)),color(1,:));
            h(2) = shade_psth(1:140,mean(t2_face2_motion,1),std(t2_face2_motion,[],1)./sqrt(sum(t2_index)),color(2,:));
            h(3) = shade_psth(1:140,mean(t2_face3_motion,1),std(t2_face3_motion,[],1)./sqrt(sum(t2_index)),color(3,:));
            h(4) = shade_psth(1:140,mean(t2_face4_motion,1),std(t2_face4_motion,[],1)./sqrt(sum(t2_index)),color(4,:));
            h(5) = shade_psth(1:140,mean(t2_face5_motion,1),std(t2_face5_motion,[],1)./sqrt(sum(t2_index)),color(5,:));
            yl = ylim; landmark_label(yl(2),'t2'); reward_label('t2');
            legend([h(1:5)],{'1st','2nd','3rd','4th','5th'}, 'Location', 'northwest');
            title('t2 face motion SVD')

        end

        subplot(8,6,[29:30,35:36])
        hold on;
        t1_pupil_size = pupil_sze_profile(t1_index,:);
        t1_pupil_angle = pupil_angle_profile(t1_index,:);
        h(1) = shade_psth(1:140,mean(t1_pupil_size,1),std(t1_pupil_size,[],1)./sqrt(sum(t1_index)),'b');
        %h(2) = shade_psth(1:140,mean(t1_pupil_angle,1),std(t1_pupil_angle,[],1)./sqrt(sum(t1_index)),'r');

        title('t1 pupil size')
        yl = ylim; landmark_label(yl(2),'t1'); reward_label('t1');

        subplot(8,6,[41:42,47:48])
        if ~isempty(c.t2_index)
            hold on;
            t2_pupil_size = pupil_sze_profile(t2_index,:);
            t2_pupil_angle = pupil_angle_profile(t2_index,:);
            h(1) = shade_psth(1:140,mean(t2_pupil_size,1),std(t2_pupil_size,[],1)./sqrt(sum(t2_index)),'b');
            %h(2) = shade_psth(1:140,mean(t2_pupil_angle,1),std(t2_pupil_angle,[],1)./sqrt(sum(t2_index)),'r');
            title('t2 pupil size')
            yl = ylim; landmark_label(yl(2),'t2'); reward_label('t2');

        end
        %
        end
    end

%% functions for plotting
    function block_index = find_block_index(track_ID_all,track_ID)
        track_number = str2double(track_ID(2));
        lap_index = track_ID_all == track_number;
        if lap_index(1) == 1
            block_indices_start = [1; find(diff(lap_index) == 1)+1];
            block_indices_end = find(diff(lap_index) == -1);
        else
            block_indices_start = [find(diff(lap_index) == 1)+1];
            block_indices_end = find(diff(lap_index) == -1);
        end
        if lap_index(end) == 1
            block_indices_end = [block_indices_end; length(lap_index)];
        end
        block_length = [1;(block_indices_end - block_indices_start)+1];
        block_index = cumsum(block_length(1:end-1));
    end
    function cluster_plot_lap(spatial_response,lap_index)
        % this function plots the spatial response of the clusters in the specified track (track_ID)
        if nargin < 2
            lap_index = 1:size(spatial_response,1);
        end
        spatial_response = spatial_response(lap_index,:);
        imagesc(spatial_response);
        xlabel('position')
        ylabel('lap')
        try
            clim([min(mean(spatial_response,'omitnan')),max(mean(spatial_response,'omitnan'))]);
        catch
        end
        colorbar;

        % for i = 2:length(block_index)
        %     yline(block_index(i),'Color','r','LineWidth',1.5)
        % end
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        hold off
    end

    function cluster_plot_both_psth(spatial_response_t1,spatial_response_t2)
        % this function plots the spatial response of the clusters in the specified track (track_ID)
        no_bin_t1 = size(spatial_response_t1,2);
        no_bin_t2 = size(spatial_response_t2,2);
        spatial_response_t1_psth = mean(spatial_response_t1,1,'omitnan');
        spatial_response_t2_psth = mean(spatial_response_t2,1,'omitnan');
        h(1) = shade_psth(1:no_bin_t1,spatial_response_t1_psth,std(spatial_response_t1,[],1,'omitnan')./sqrt(size(spatial_response_t1,1)),'b');
        hold on;
        h(2) = shade_psth(1:no_bin_t2,spatial_response_t2_psth,std(spatial_response_t2,[],1,'omitnan')./sqrt(size(spatial_response_t2,1)),'r');
        ylim([min([min(spatial_response_t2_psth),min(spatial_response_t1_psth)])-2,max([max(spatial_response_t2_psth),max(spatial_response_t1_psth)])+2])
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('firing rate')
        xline([30 50 70 90 110 130],'LineWidth',1.5,'HandleVisibility','off')
        legend([h(1:2)],{'t1','t2'}, 'Location', 'northwest');
    end

    function h =  cluster_plot_psth(spatial_response,lap_index)
        spatial_response_psth_lap = mean(spatial_response(lap_index,:), 1, 'omitnan');
        spatial_response_psth_nonlap = mean(spatial_response(~lap_index,:), 1, 'omitnan');
        spatial_response_std_lap = std(spatial_response(lap_index,:), [], 1, 'omitnan')./sqrt(sum(lap_index));
        spatial_response_std_nonlap = std(spatial_response(~lap_index,:), [], 1, 'omitnan')./sqrt(sum(~lap_index));
        h(1) = shade_psth(1:size(spatial_response, 2), spatial_response_psth_lap,spatial_response_std_lap , 'b');
        hold on;
        h(2) = shade_psth(1:size(spatial_response, 2), spatial_response_psth_nonlap,spatial_response_std_nonlap, 'r');
        set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);
        xlabel('position');
        ylabel('firing rate');
    end

    function h = shade_psth(x,y,ysd,color)
        patch([x fliplr(x)], [(y+ysd) fliplr((y-ysd))],color,'FaceAlpha',0.2,'EdgeColor','none');
        hold on;
        h = plot(x,y,'Color', color, 'LineWidth', 2);
    end
    function cluster_plot_lap_extended(spatial_response,block_index)
        % this function plots the spatial response of the clusters in the specified track (track_ID)
        no_lap = length(spatial_response);
        for iL = 1:no_lap
            no_bin = length(spatial_response{iL});
            plot(1:no_bin,spatial_response{iL} + iL*ones(size(spatial_response{iL})))
            hold on
        end
        xlabel('position')
        ylabel('lap')
        title(['track 1 normalized'])
        for i = 2:length(block_index)
            yline(block_index(i),'Color','r','LineWidth',1.5)
        end
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
    end

    function cluster_plot_both_extended_psth(spatial_response_t1,spatial_response_t2)

        no_lap_t1 = length(spatial_response_t1);
        %for psth, only use the first 200cm - its misleading if only a few laps have long distance
        spatial_response_t1_mat = zeros(no_lap_t1,200);
        for iL = 1:no_lap
            spatial_response_t1_mat(iL,:) = spatial_response_t1{iL}(1:200);
        end
        no_lap_t2 = length(spatial_response_t2);
        spatial_response_t2_mat = zeros(no_lap,200);
        for iL = 1:no_lap
            spatial_response_t2_mat(iL,:) = spatial_response_t2{iL}(1:200);
        end
        spatial_response_t1_mat_psth = mean(spatial_response_t1_mat,1,'omitnan');
        spatial_response_t2_mat_psth = mean(spatial_response_t2_mat,1,'omitnan');
        shade_psth(1:200,spatial_response_t1_mat_psth,std(spatial_response_t1_mat,[],1,'omitnan'),'b');
        hold on;
        shade_psth(1:200,spatial_response_t2_mat_psth,std(spatial_response_t2_mat,[],1,'omitnan'),'r');
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('firing rate')
        ylim([min([min(spatial_response_t2_mat_psth),min(spatial_response_t1_mat_psth)])-2,max([max(spatial_response_t2_mat_psth),max(spatial_response_t1_mat_psth)])+2])
    end

    function h = comapre_landmark_response(spatial_response,landmark_bin)
        no_lap = size(spatial_response,1);
        all_laps = 1:no_lap;
        even_laps = all_laps(mod(all_laps, 2) == 0);
        odd_laps = all_laps(mod(all_laps, 2) == 1);
        spatial_response = squeeze(spatial_response);
        spatial_response_psth_odd = mean(spatial_response(odd_laps,:),1,'omitnan');
        odd_peaks = islocalmax(smooth(spatial_response_psth_odd),'MinSeparation',15);
        spatial_response_psth_even = mean(spatial_response(even_laps,:),1,'omitnan');
        even_peaks = islocalmax(smooth(spatial_response_psth_even),'MinSeparation',15);
        bins = -20:20;
        %spatial_response_std_odd = std(spatial_response(odd_laps,:),[],1,'omitnan')./sqrt(length(odd_laps));
        %spatial_response_std_even = std(spatial_response(even_laps,:),[],1,'omitnan')./sqrt(length(even_laps));
        h(1) = plot(bins,spatial_response_psth_odd(landmark_bin(1,:)),'Color',[214,96,77]/255,'LineWidth',2);
        %h(1) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_odd(landmark_bin(1,:)),spatial_response_std_odd(landmark_bin(1,:)),[214,96,77]/255);
        hold on;
        h(2) = plot(bins,spatial_response_psth_even(landmark_bin(1,:)),'Color',[244,165,130]/255,'LineWidth',2);
        h(3) = plot(bins,spatial_response_psth_odd(landmark_bin(2,:)),'Color',[67,147,195]/255,'LineWidth',2);
        h(4) = plot(bins,spatial_response_psth_even(landmark_bin(2,:)),'Color',[146,197,222]/255,'LineWidth',2);
        scatter(bins(odd_peaks(landmark_bin(1,:))),spatial_response_psth_odd(landmark_bin(1,odd_peaks(landmark_bin(1,:)))),20,'k','filled');
        scatter(bins(even_peaks(landmark_bin(1,:))),spatial_response_psth_even(landmark_bin(1,even_peaks(landmark_bin(1,:)))),20,'k','filled');
        scatter(bins(odd_peaks(landmark_bin(2,:))),spatial_response_psth_odd(landmark_bin(2,odd_peaks(landmark_bin(2,:)))),20,'k','filled');
        scatter(bins(even_peaks(landmark_bin(2,:))),spatial_response_psth_even(landmark_bin(2,even_peaks(landmark_bin(2,:)))),20,'k','filled');
        %h(2) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_even(landmark_bin(1,:)),spatial_response_std_even(landmark_bin(1,:)),[244,165,130]/255);
        %h(3) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_odd(landmark_bin(2,:)),spatial_response_std_odd(landmark_bin(2,:)),[67,147,195]/255);
        %h(4) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_even(landmark_bin(2,:)),spatial_response_std_even(landmark_bin(2,:)),[146,197,222]/255);
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xline([-10,10],'k','LineWidth',1.5,'LineStyle','--','HandleVisibility','off')
        ylabel('firing rate')

    end

    function h = compare_common_landmark_between_t1_and_t2(spatial_response_t1,spatial_response_t2)
        landmark_bin = [50:90; 90:130];
        spatial_response_psth_t1 = mean(spatial_response_t1,1,'omitnan');
        %spatial_response_std_t1 = std(spatial_response_t1,[],1,'omitnan')./sqrt(size(spatial_response_t1,1));
        spatial_response_psth_t2 = mean(spatial_response_t2,1,'omitnan');
        %spatial_response_std_t2 = std(spatial_response_t2,[],1,'omitnan')./sqrt(size(spatial_response_t2,1));
        %h(1) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_t1(landmark_bin(1,:)),spatial_response_std_t1(landmark_bin(1,:)),[214,96,77]/255);
        h(1) = plot(-20:20,spatial_response_psth_t1(landmark_bin(1,:)),'Color',[67,147,195]/255,'LineWidth',2);
        hold on;
        h(2) = plot(-20:20,spatial_response_psth_t1(landmark_bin(2,:)),'Color',[146,197,222]/255,'LineWidth',2);
        h(3) = plot(-20:20,spatial_response_psth_t2(landmark_bin(1,:)),'Color',[214,96,77]/255,'LineWidth',2);
        h(4) = plot(-20:20,spatial_response_psth_t2(landmark_bin(2,:)),'Color',[244,165,130]/255,'LineWidth',2);
        %{
 h(2) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_t1(landmark_bin(2,:)),spatial_response_std_t1(landmark_bin(2,:)),[244,165,130]/255);
    h(3) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_t2(landmark_bin(1,:)),spatial_response_std_t2(landmark_bin(1,:)),[67,147,195]/255);
    h(4) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_t2(landmark_bin(2,:)),spatial_response_std_t2(landmark_bin(2,:)),[146,197,222]/255);
     
        %}
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xline([-10,10],'k','LineWidth',1.5,'LineStyle','--','HandleVisibility','off')

        ylabel('firing rate')
        lgd = legend([h(1:4)],{'B1','B2','B3','B4'}, 'Location', 'best');
        lgd.FontSize = 8;


    end

    function cluster_even_odd_psth(spatial_response)
        % this function plots the spatial response of the clusters in the specified track (track_ID)
        no_lap = size(spatial_response,1);
        no_bin = size(spatial_response,2);
        all_laps = 1:no_lap;
        even_laps = all_laps(mod(all_laps, 2) == 0);
        odd_laps = all_laps(mod(all_laps, 2) == 1);
        h(1) = shade_psth(1:no_bin,mean(spatial_response(even_laps,:),1,'omitnan'),std(spatial_response(even_laps,:),[],1,'omitnan')./sqrt(length(even_laps)),'b');
        hold on;
        h(2) = shade_psth(1:no_bin,mean(spatial_response(odd_laps,:),1,'omitnan'),std(spatial_response(odd_laps,:),[],1,'omitnan')./sqrt(length(odd_laps)),'r');
        legend([h(1:2)],{'even','odd'})
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('firing rate')
        ylim([min(mean(spatial_response,'omitnan'))-5,max(mean(spatial_response,'omitnan'))+5])

    end

    function cluster_block_psth(spatial_response,block_index)
        no_block = length(block_index);
        block_psth = zeros(no_block,size(spatial_response,2));
        block_std = zeros(no_block,size(spatial_response,2));
        color = sky(no_block);
        for iB = 1:no_block-1
            block_psth(iB,:) = mean(spatial_response(block_index(iB):block_index(iB+1),:),1,'omitnan');
            block_std(iB,:) = std(spatial_response(block_index(iB):block_index(iB+1),:),[],1,'omitnan');
        end
        block_psth(no_block,:) = mean(spatial_response(block_index(end):end,:),1,'omitnan');
        block_psth(no_block,:) = std(spatial_response(block_index(end):end,:),1,'omitnan');
        for iB = 1:no_block-1
            hold on;
            %h(iB) = shade_psth(1:size(spatial_response,2),block_psth(iB,:),block_std(iB,:)./sqrt(block_index(iB+1)-block_index(iB)),color(iB,:));
            legend_names{iB} = ['block ',num2str(iB)];
            h(iB) = plot(1:size(spatial_response,2),block_psth(iB,:),'Color',color(iB,:));

        end
        hold on;
        %h(no_block) = shade_psth(1:size(spatial_response,2),block_psth(no_block,:),block_std(no_block,:)./sqrt(size(spatial_response,1)-block_index(end)),color(no_block,:));
        h(no_block) = plot(1:size(spatial_response,2),block_psth(no_block,:),'Color',color(no_block,:));
        legend_names{no_block} = ['block ',num2str(no_block)];
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('firing rate')
        legend([h(1:no_block)],legend_names);
        ylim([min(mean(spatial_response,'omitnan'))-5,max(mean(spatial_response,'omitnan'))+5])
        title('psth of each block')
    end

    function cluster_block_dynamic(spatial_response,block_index)
        no_block = length(block_index);
        max_lap = 6;
        lap_resp = zeros(max_lap,no_block,size(spatial_response,2));
        color = sky(max_lap);
        for iB = 1:no_block
            for iL = 1:max_lap
                if block_index(iB)+iL-1 <= size(spatial_response,1)
                    lap_resp(iL,iB,:) = spatial_response(block_index(iB)+iL-1,:);
                end
            end
        end
        lap_psth = squeeze(mean(lap_resp,2,'omitnan'));
        %lap_std =squeeze(std(lap_resp,[],2,'omitnan'));


        for iL = 1:max_lap
            h(iL) = plot(1:size(spatial_response,2),lap_psth(iL,:),'Color',color(iL,:));
            hold on;
            lgend_names{iL} = ['lap ',num2str(iL)];
        end
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('firing rate')
        title('lap psth')
        ylim([min(mean(spatial_response,'omitnan'))-5,max(mean(spatial_response,'omitnan'))+5])
        legend([h(1:max_lap)],lgend_names);
    end

    function plot_speed_response(speed_response)
        plot(1:2:49,speed_response{1},'LineWidth',2);
        hold on;
        plot(1:2:49,speed_response{2},'LineWidth',2);
        plot(1:2:49,speed_response{3},'LineWidth',2);
        legend('t1','t2','non-track')
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('speed')
        ylabel('firing rate')

    end

    function SMI = cluster_SMI(spatial_response_t1,spatial_response_t2)
        % caluculate the spatial modulation index of the cluster by even vs odd averages
        SMI = cell(1,2); % cell 1 for track 1 and cell 2 for track 2
        SMI_t1 = nan([1,2]);
        SMI_t2 = nan([1,3]);
        no_lap = size(spatial_response_t1,1);
        all_laps = 1:no_lap;
        even_laps = all_laps(mod(all_laps, 2) == 0);
        odd_laps = all_laps(mod(all_laps, 2) == 1);
        odd_avg = mean(spatial_response_t1(odd_laps,:),1,'omitnan');
        even_avg = mean(spatial_response_t1(even_laps,:),1,'omitnan');
        landmark_position = [41:60, 81:100; 61:80, 101:120];
        [~,odd_peak_location_l1] = max(odd_avg(landmark_position(1,:)),[],'omitnan');
        [~,odd_peak_location_l2] = max(odd_avg(landmark_position(2,:)),[],'omitnan');
        non_prefer_location_l1(odd_peak_location_l1 > 20) = odd_peak_location_l1 - 20;
        non_prefer_location_l1(odd_peak_location_l1 <= 20) = odd_peak_location_l1 + 20;
        non_prefer_location_l2(odd_peak_location_l2 > 20) = odd_peak_location_l2 - 20;
        non_prefer_location_l2(odd_peak_location_l2 <= 20) = odd_peak_location_l2 + 20;
        SMI_t1(1,1) = SMI_calculator(even_avg(landmark_position(1,odd_peak_location_l1)),even_avg(landmark_position(1,non_prefer_location_l1)));
        SMI_t1(1,2) = SMI_calculator(even_avg(landmark_position(2,odd_peak_location_l2)),even_avg(landmark_position(2,non_prefer_location_l2)));
        SMI{1,1} = SMI_t1;
        no_lap = size(spatial_response_t2,1);
        all_laps = 1:no_lap;
        even_laps = all_laps(mod(all_laps, 2) == 0);
        odd_laps = all_laps(mod(all_laps, 2) == 1);
        odd_avg = mean(spatial_response_t2(odd_laps,:),1,'omitnan');
        even_avg = mean(spatial_response_t2(even_laps,:),1,'omitnan');
        landmark_position = [41:60, 81:100; 61:80, 101:120; 21:40, 121:140];
        [~,odd_peak_location_l1] = max(odd_avg(landmark_position(1,:)),[],'omitnan');
        [~,odd_peak_location_l2] = max(odd_avg(landmark_position(2,:)),[],'omitnan');
        [~,odd_peak_location_l3] = max(odd_avg(landmark_position(3,:)),[],'omitnan');
        non_prefer_location_l1(odd_peak_location_l1 > 20) = odd_peak_location_l1 - 20;
        non_prefer_location_l1(odd_peak_location_l1 <= 20) = odd_peak_location_l1 + 20;
        non_prefer_location_l2(odd_peak_location_l2 > 20) = odd_peak_location_l2 - 20;
        non_prefer_location_l2(odd_peak_location_l2 <= 20) = odd_peak_location_l2 + 20;
        non_prefer_location_l3(odd_peak_location_l3 > 20) = odd_peak_location_l3 - 20;
        non_prefer_location_l3(odd_peak_location_l3 <= 20) = odd_peak_location_l3 + 20;
        SMI_t2(1,1) = SMI_calculator(even_avg(landmark_position(1,odd_peak_location_l1)),even_avg(landmark_position(1,non_prefer_location_l1)));
        SMI_t2(1,2) = SMI_calculator(even_avg(landmark_position(2,odd_peak_location_l2)),even_avg(landmark_position(2,non_prefer_location_l2)));
        SMI_t2(1,3) = SMI_calculator(even_avg(landmark_position(3,odd_peak_location_l3)),even_avg(landmark_position(3,non_prefer_location_l3)));
        SMI{1,2} = SMI_t2;

        function SMI = SMI_calculator(even_prefer,even_non_prefer)
            SMI = (even_prefer - even_non_prefer)/(even_prefer + even_non_prefer);
        end
    end

    function track_label(track_ID,y_lim)
        switch track_ID
            case 't1'
                xline(30,'LineWidth',1.5,'HandleVisibility','off')
                text(30, y_lim, 'X', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline([50 90],'LineWidth',1.5,'Color',landmark_color('A',1),'HandleVisibility','off')
                text(50, y_lim, 'A1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                text(90, y_lim, 'A2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline([70 110],'LineWidth',1.5,'Color',landmark_color('B',1),'HandleVisibility','off')
                text(70, y_lim, 'B1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                text(110, y_lim, 'B2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            case 't2'
                xline(30,'LineWidth',1.5,'HandleVisibility','off')
                text(30, y_lim, 'Y1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline([50 90],'LineWidth',1.5,'Color',landmark_color('C',1),'HandleVisibility','off')
                text(50, y_lim, 'C1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                text(90, y_lim, 'C2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline([70 110],'LineWidth',1.5,'Color',landmark_color('B',1),'HandleVisibility','off')
                text(70, y_lim, 'B3', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                text(110, y_lim, 'B4', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline(130,'LineWidth',1.5,'HandleVisibility','off')
                text(130, y_lim, 'Y2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            case 'both'
                xline(30,'LineWidth',1.5,'HandleVisibility','off')
                text(30, y_lim, 'X-Y', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline([50 90],'LineWidth',1.5,'Color',landmark_color('A',1),'HandleVisibility','off')
                text(50, y_lim, 'A-C1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                text(90, y_lim, 'A-C2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline([70 110],'LineWidth',1.5,'Color',landmark_color('B',1),'HandleVisibility','off')
                text(70, y_lim, 'B1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                text(110, y_lim, 'B2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                xline(130,'LineWidth',1.5,'LineStyle','--','HandleVisibility','off')
                text(130, y_lim, 'Y2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

        end
    end

    function landmark_label(no_lap,track_ID)

        switch track_ID
            case 't1'
                landmark_vertex =zeros(5,4,2);
                landmark_vertex(1,:,1) =[26 34 34 26];
                landmark_vertex(2,:,1) =[46 54 54 46];
                landmark_vertex(3,:,1) =[66 74 74 66];
                landmark_vertex(4,:,1) =[86 94 94 86];
                landmark_vertex(5,:,1) =[106 114 114 106];
                for i = 1:5
                    landmark_vertex(i,:,2) = [0 0 no_lap no_lap];
                end
                patch_color = [landmark_color('X',1);
                    landmark_color('A',1);
                    landmark_color('B',1);
                    landmark_color('A',1);
                    landmark_color('B',1)];
                hold on;
                for i = 1:5
                    patch(landmark_vertex(i,:,1),landmark_vertex(i,:,2),patch_color(i,:),'EdgeColor','none','FaceAlpha',0.1)
                end
                text(30,0,'X','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(50,0,'A1','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(70,0,'B1','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(90,0,'A2','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(110,0,'B2','HorizontalAlignment','center','VerticalAlignment','bottom')
            case 't2'
                landmark_vertex =zeros(6,4,2);
                landmark_vertex(1,:,1) =[26 34 34 26];
                landmark_vertex(2,:,1) =[46 54 54 46];
                landmark_vertex(3,:,1) =[66 74 74 66];
                landmark_vertex(4,:,1) =[86 94 94 86];
                landmark_vertex(5,:,1) =[106 114 114 106];
                landmark_vertex(6,:,1) =[126 134 134 126];
                for i = 1:6
                    landmark_vertex(i,:,2) = [0 0 no_lap no_lap];
                end
                patch_color = [landmark_color('Y',1);
                    landmark_color('C',1);
                    landmark_color('B',1);
                    landmark_color('C',1);
                    landmark_color('B',1);
                    landmark_color('Y',1)];
                hold on;
                for i = 1:6
                    patch(landmark_vertex(i,:,1),landmark_vertex(i,:,2),patch_color(i,:),'EdgeColor','none','FaceAlpha',0.1)
                end
                text(30,0,'Y1','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(50,0,'C1','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(70,0,'B3','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(90,0,'C2','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(110,0,'B4','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(130,0,'Y1','HorizontalAlignment','center','VerticalAlignment','bottom')

        end

    end
    function color_output = landmark_color(landmark,color_choice)
        colors = zeros(5,3,2);
        colors(:,:,1) = [146,197,222;
            214,96,77;
            166,219,160;
            194,165,207;
            255,165,0]/255;
        colors(:,:,2) = [67,147,195;
            244,165,130;
            90,174,97;
            153,112,171;
            255,140,0]/255;

        switch landmark
            case 'A'
                color_output = colors(1,:,color_choice);
            case 'B'
                color_output = colors(2,:,color_choice);
            case 'C'
                color_output = colors(3,:,color_choice);
            case 'X'
                color_output = colors(4,:,color_choice);
            case 'Y'
                color_output = colors(5,:,color_choice);
        end
        color_output = squeeze(color_output);
    end

    function behaviour_profile = behaviour_profile_calculator(tvec,behaviour,position,start_time,end_time)
        w = gausswin(15);
        w =     w / sum(w);
        behaviour(isnan(behaviour)) = 0;
        behaviour= filtfilt(w,1,behaviour')';
        bin_edges = 0:140;bin_centres = 0.5:1:139.5;
        behaviour_profile = zeros(length(start_time),length(bin_centres));
        for iL = 1:length(start_time)
            lap_start = start_time(iL);
            lap_end = end_time(iL);
            lap_tvec_index = tvec >= lap_start & tvec <= lap_end;
            lap_behaviour = behaviour(lap_tvec_index);
            lap_position = position(lap_tvec_index);
            position_index = discretize(lap_position,bin_edges);
            for iB = 1:length(bin_centres)
                bin_index = position_index == iB;
                behaviour_profile(iL,iB) = mean(lap_behaviour(bin_index),'omitnan');
            end
        end
        behaviour_profile = fillmissing(behaviour_profile,'linear',1);
    end
    function [lap_lick_position,lap_reward_position] = lick_reward_position_calculator(lick_position,reward_position,lick_time,track_ID_all,start_time_all,end_time_all)
        lap_lick_position = cell(length(start_time_all),1);
        lap_reward_position = zeros(length(start_time_all),1);
        for iL = 1:length(start_time_all)
            lap_start = start_time_all(iL);
            lap_end = end_time_all(iL);
            lap_lick_index = lick_time >= lap_start & lick_time <= lap_end;
            %turn lick position into 0-140cm
            if track_ID_all(iL) == 1
                lap_lick_position_tmp = lick_position(lap_lick_index)+1140;
                reward_zone_lick = lap_lick_position_tmp >= 115;
                if sum(reward_zone_lick) > 1
                    reward_zone_lick_index = find(reward_zone_lick);
                    lap_lick_position_tmp(reward_zone_lick_index(2:end)) = NaN;
                end
                lap_lick_position{iL} = lap_lick_position_tmp;
                lap_reward_position(iL) = reward_position(iL)+1140;
            elseif track_ID_all(iL) == 2
                lap_lick_position_tmp = lick_position(lap_lick_index)+140;
                reward_zone_lick = lap_lick_position_tmp >= 115;
                if sum(reward_zone_lick) > 1
                    reward_zone_lick_index = find(reward_zone_lick);
                    lap_lick_position_tmp(reward_zone_lick_index(2:end)) = NaN;
                end
                lap_lick_position{iL} = lap_lick_position_tmp;
                lap_reward_position(iL) = reward_position(iL)+140;
            end
        end
    end
    function session_behaviour_plot(speed_profile,lap_lick_position,reward_type,reward_position,block_index)
        no_lap = size(speed_profile,1);
        active_laps = reward_type == 10;
        for iL = 1:no_lap
            hold on;
            scatter(lap_lick_position{iL},ones(1,length(lap_lick_position{iL}))*(iL-1)+0.5,10,'k','*');
        end
        laps = 0:no_lap-1;
        scatter(reward_position(active_laps),0.5+laps(active_laps),30,'r','o','filled');
        scatter(reward_position(~active_laps),0.5+laps(~active_laps),30,'b','o','filled');
        %plot(1:140,speed_profile./20 + (laps)','LineWidth',2,'Color','k');
        for i = 2:length(block_index)
            yline(block_index(i),'Color','r','LineWidth',1.5)
            yline(block_index(i),'Color','r','LineWidth',1.5)
            yline(block_index(i),'Color','r','LineWidth',1.5)
        end
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('lap')
        ylim([0,no_lap])
    end
    function lick_psth(lap_lick_position,lap_index)
        position_edges = 0:140; bins = 0.5:1:139.5;
        lick_position = zeros(length(lap_lick_position),140);
        for iL = 1:length(lap_lick_position)
            lick_position(iL,:) = histcounts(lap_lick_position{iL},position_edges);
            lick_position(iL,:) = smooth(lick_position(iL,:));
        end

        shade_psth(bins,mean(lick_position(lap_index,:),1,'omitnan'),std(lick_position(lap_index,:),[],1)./sqrt(sum(lap_index)),'r');
        hold on;
        shade_psth(bins,mean(lick_position(~lap_index,:),1,'omitnan'),std(lick_position(~lap_index,:),[],1)./sqrt(sum(~lap_index)),'b');
        set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        xlabel('position')
        ylabel('lick rate')
    end
    function reward_label(track_ID)
        yl = ylim;
        switch track_ID
            case 't1'
                xline(120,'Color','r','LineWidth',3,'HandleVisibility','off')
                text(120,yl(2),'reward zone','HorizontalAlignment','center','VerticalAlignment','bottom','Color','r')
                x1 = [115, 120, 120, 115]; x2 = [120, 125, 125, 120];
                y = [yl(1),yl(1),yl(2),yl(2)];
                %patch(x1, y, 'r', 'FaceAlpha', 0.05,'EdgeColor','none');
                %patch(x2, y, 'b', 'FaceAlpha', 0.05,'EdgeColor','none');
            case 't2'
                xline(135,'Color','r','LineWidth',3,'HandleVisibility','off')
                text(135,yl(2),'reward zone','HorizontalAlignment','center','VerticalAlignment','bottom','Color','r')
                x1 = [130, 135,135,130]; x2 = [135,140, 140, 135];
                y = [yl(1),yl(1),yl(2),yl(2)];
        end

    end
    function speed_psth(speed_profile,color)
        mean_speed = mean(speed_profile,1,'omitnan');
        plot(1:140,speed_profile,'Color',color,'LineWidth',0.5);
        hold on;
        shade_psth(1:140,mean_speed,std(speed_profile,[],1,'omitnan')./sqrt(size(speed_profile,1)),'k');
        ylim([min(mean_speed)-2*max(std(speed_profile,[],1,'omitnan')),max(mean_speed)+2*max(std(speed_profile,[],1,'omitnan'))])
    end

end
