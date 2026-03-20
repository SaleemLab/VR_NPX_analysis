classdef cluster_spatial_info
    properties
        session_count
        tvec
        speed
        face_motion_SVD
        pupil_size
        pupil_angle
        position
        start_time_all
        end_time_all
        rewarded_laps
        track_ID_all
        rewarded_laps_index
        t1_index
        t2_index
        rewarded_laps_t1
        rewarded_laps_t2
        reward_type
        active_reward_laps_t1
        active_reward_laps_t2
        id
        raw_t1
        raw_t2
        spike_times
        spatial_response
        spatial_response_extended
        waveform
        speed_response
        reward_position
        lick_position
        lick_time
        block_index_t1
        block_index_t2
        SMI
        peak_depth
        region
    end

    methods
        function obj = cluster_spatial_info(clusters_all, cluster_id, spatial_response, spatial_response_extended)
            % Get the field names of clusters_all
            fields = fieldnames(clusters_all);

            % % Loop over the field names
            % for i = 1:numel(fields)
            %     % Add a dynamic property with the current field name
            %     obj.addprop(fields{i});
            % end
            obj.id = cluster_id;
            ids_index = ismember(clusters_all.cluster_id,cluster_id);
            obj.session_count = clusters_all.session_count(ids_index);
            obj.tvec = clusters_all.tvec{obj.session_count};
            obj.speed = clusters_all.speed{obj.session_count};
%             obj.face_motion_SVD = clusters_all.face_motion_SVD{obj.session_count};
%             obj.pupil_size = clusters_all.pupil_size{obj.session_count};
%             obj.pupil_angle = clusters_all.pupil_movement_angle{obj.session_count};
            obj.position = clusters_all.position{obj.session_count};
            obj.start_time_all = clusters_all.start_time_all{obj.session_count};
            obj.end_time_all = clusters_all.end_time_all{obj.session_count};
            obj.rewarded_laps = clusters_all.rewarded_lap_id{obj.session_count};
            obj.track_ID_all = clusters_all.track_ID_all{obj.session_count};
            obj.rewarded_laps_index = zeros(length(obj.track_ID_all),1);
            obj.rewarded_laps_index(obj.rewarded_laps) = 1;
            obj.rewarded_laps_index = logical(obj.rewarded_laps_index);
            obj.t1_index = ones(sum(obj.track_ID_all == 1),1);
            obj.t2_index = ones(sum(obj.track_ID_all == 2),1);
            obj.rewarded_laps_t1 = obj.t1_index & obj.rewarded_laps_index(obj.track_ID_all == 1);
            obj.rewarded_laps_t2 = obj.t2_index & obj.rewarded_laps_index(obj.track_ID_all == 2);
            obj.reward_type = nan(size(obj.rewarded_laps_index));
            obj.reward_type(obj.rewarded_laps_index) = clusters_all.reward_type{obj.session_count};
            obj.active_reward_laps_t1 = obj.reward_type(obj.track_ID_all == 1 & obj.rewarded_laps_index) == 10 ;
            obj.active_reward_laps_t2 = obj.reward_type(obj.track_ID_all == 2 & obj.rewarded_laps_index) == 10 ;
            obj.lick_time = clusters_all.lick_time{obj.session_count};
            %responses
            obj.raw_t1 = clusters_all.raw{ids_index,1};
            obj.raw_t2 = clusters_all.raw{ids_index,2};
            obj.spike_times = clusters_all.spike_times{ids_index};
            obj.spatial_response = spatial_response;
%             obj.spatial_response{1} = reshape(normalize(reshape(obj.spatial_response{1},1,[]),'range'),size(obj.spatial_response{1}));
%             obj.spatial_response{2} = reshape(normalize(reshape(obj.spatial_response{2},1,[]),'range'),size(obj.spatial_response{2}));
            obj.spatial_response_extended = spatial_response_extended;
            obj.waveform = clusters_all.peak_channel_waveforms(ids_index,:);
            obj.speed_response = cluster_speed_tuning(obj.spike_times,obj.tvec,obj.speed,...
                obj.track_ID_all,obj.start_time_all,obj.end_time_all);
            obj.reward_position = nan(size(obj.rewarded_laps_index));
            obj.reward_position(obj.rewarded_laps_index) = clusters_all.reward_position{obj.session_count};
            obj.lick_position = clusters_all.lick_position{obj.session_count};
            obj.peak_depth = clusters_all.peak_depth(ids_index,1);
            obj.region = clusters_all.region(ids_index,1);
            obj.block_index_t1 = find_block_index(obj,'t1');
            obj.block_index_t2 = find_block_index(obj,'t2');
            obj.SMI = cluster_SMI(obj);
        end

        function block_index = find_block_index(obj,track_ID)
            track_number = str2double(track_ID(2));
            track_ID_all_r = obj.track_ID_all(obj.rewarded_laps);
            lap_index = track_ID_all_r == track_number;
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
        function SMI = cluster_SMI(obj)
            % caluculate the spatial modulation index of the cluster by even vs odd averages
            % Normalize the responses by range (0,1) for each track individually across laps.
            % Calculate the average response for odd and even laps in each track.
            % Find the local peaks using islocalmax() function on odd average responses.
            % Find the peak location using max() of local peak values in the identical landmark regions.
            % Find the non-prefered peak of the other region.
            % Calculate the SMI using the formula (prefer - non-prefered)/(prefer + non-prefered) using even laps average responses.
            SMI = cell(1,2); % cell 1 for track 1 and cell 2 for track 2
            SMI_t1 = nan([1,2]);
            SMI_t2 = nan([1,3]);
            spatial_response_t1_norm = obj.spatial_response{1}(obj.rewarded_laps_t1,:);
            spatial_response_t2_norm = obj.spatial_response{2}(obj.rewarded_laps_t2,:);
            %spatial_response_t1_norm = reshape(normalize(reshape(obj.spatial_response{1}(obj.rewarded_laps_t1,:),1,[]),'range'),size(obj.spatial_response{1}(obj.rewarded_laps_t1,:)));
            %spatial_response_t2_norm = reshape(normalize(reshape(obj.spatial_response{2}(obj.rewarded_laps_t2,:),1,[]),'range'),size(obj.spatial_response{2}(obj.rewarded_laps_t2,:)));
            landmark_position = [41:60, 81:100; 61:80, 101:120];
            SMI_t1(1,1) = SMI_landmark(spatial_response_t1_norm,landmark_position(1,:));
            SMI_t1(1,2) = SMI_landmark(spatial_response_t1_norm,landmark_position(2,:));
            SMI{1,1} = SMI_t1;
            if size(spatial_response_t2_norm,1) > 0
                landmark_position = [41:60, 81:100; 61:80, 101:120; 21:40,121:140];
                SMI_t2(1,1) = SMI_landmark(spatial_response_t2_norm,landmark_position(1,:));
                SMI_t2(1,2) = SMI_landmark(spatial_response_t2_norm,landmark_position(2,:));
                SMI_t2(1,3) = SMI_landmark(spatial_response_t2_norm,landmark_position(3,:));
                SMI{1,2} = SMI_t2;
            else
                SMI{1,2} = SMI_t2;
            end

            function SMI = SMI_calculator(even_prefer,even_non_prefer)
                SMI = (even_prefer - even_non_prefer)/(even_prefer + even_non_prefer);
            end
            function SMI = SMI_landmark(spatial_response,landmark_position)
                no_lap = size(spatial_response,1);
                all_laps = 1:no_lap;
                even_laps = all_laps(mod(all_laps, 2) == 0);
                odd_laps = all_laps(mod(all_laps, 2) == 1);
                odd_avg = normalize(mean(spatial_response(odd_laps,:),1,'omitnan'),'range');
                even_avg = normalize(mean(spatial_response(even_laps,:),1,'omitnan'),'range');
                odd_peaks = islocalmax(smooth(odd_avg),'MinSeparation',15);
                potential_peak_location = landmark_position(:,odd_peaks(landmark_position));
                [~,odd_peak_location_index] = max(odd_avg(potential_peak_location),[],'omitnan');
                odd_peak_location = potential_peak_location(odd_peak_location_index);
                if odd_avg(odd_peak_location) >= 0.05
                    odd_peak_location_index = find(landmark_position(1,:) == odd_peak_location);
                    if odd_peak_location_index > 20
                        non_prefer_location= landmark_position(1,odd_peak_location_index - 20);
                    else
                        non_prefer_location = landmark_position(1,odd_peak_location_index + 20);
                    end
                    SMI= SMI_calculator(even_avg(odd_peak_location),even_avg(non_prefer_location));
                else
                    SMI = NaN;
                end

            end
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
            spatial_response_psth_even = mean(spatial_response(even_laps,:),1,'omitnan');
            spatial_response_std_odd = std(spatial_response(odd_laps,:),[],1,'omitnan')./sqrt(length(odd_laps));
            spatial_response_std_even = std(spatial_response(even_laps,:),[],1,'omitnan')./sqrt(length(even_laps));
            h(1) = plot(-20:20,spatial_response_psth_odd(landmark_bin(1,:)),'Color',[214,96,77]/255,'LineWidth',2);
            %h(1) = shade_psth(1:size(landmark_bin,2),spatial_response_psth_odd(landmark_bin(1,:)),spatial_response_std_odd(landmark_bin(1,:)),[214,96,77]/255);
            hold on;
            h(2) = plot(-20:20,spatial_response_psth_even(landmark_bin(1,:)),'Color',[244,165,130]/255,'LineWidth',2);
            h(3) = plot(-20:20,spatial_response_psth_odd(landmark_bin(2,:)),'Color',[67,147,195]/255,'LineWidth',2);
            h(4) = plot(-20:20,spatial_response_psth_even(landmark_bin(2,:)),'Color',[146,197,222]/255,'LineWidth',2);
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
            spatial_response_std_t1 = std(spatial_response_t1,[],1,'omitnan')./sqrt(size(spatial_response_t1,1));
            spatial_response_psth_t2 = mean(spatial_response_t2,1,'omitnan');
            spatial_response_std_t2 = std(spatial_response_t2,[],1,'omitnan')./sqrt(size(spatial_response_t2,1));
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
            lgd = legend([h(1:4)],{'t1 r1','t1 r2','t2 r1','t2 r2'}, 'Location', 'best');
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

        function speed_profile = speed_profile_calculator(tvec,speed,position,start_time,end_time)
            w = gausswin(15);
            w =     w / sum(w);
            speed(isnan(speed)) = 0;
            speed= filtfilt(w,1,speed')';
            bin_edges = 0:140;bin_centres = 0.5:1:139.5;
            speed_profile = zeros(length(start_time),length(bin_centres));
            for iL = 1:length(start_time)
                lap_start = start_time(iL);
                lap_end = end_time(iL);
                lap_tvec_index = tvec >= lap_start & tvec <= lap_end;
                lap_speed = speed(lap_tvec_index);
                lap_position = position(lap_tvec_index);
                position_index = discretize(lap_position,bin_edges);
                for iB = 1:length(bin_centres)
                    bin_index = position_index == iB;
                    speed_profile(iL,iB) = mean(lap_speed(bin_index),'omitnan');
                end
            end
            speed_profile = fillmissing(speed_profile,'linear',1);
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

    end
end