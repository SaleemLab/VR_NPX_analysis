classdef cluster_spatial_info
    properties
        session_count
        tvec
        speed
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
            obj.spatial_response_extended = spatial_response_extended;
            obj.waveform = clusters_all.peak_channel_waveforms(ids_index,:);
            obj.speed_response = cluster_speed_tuning(obj.spike_times,obj.tvec,obj.speed,...
                obj.track_ID_all,obj.start_time_all,obj.end_time_all);

            obj.reward_position = clusters_all.reward_position{obj.session_count};
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
    
        function SMI = cluster_SMI(obj)
            % caluculate the spatial modulation index of the cluster by even vs odd averages
            SMI = cell(1,2); % cell 1 for track 1 and cell 2 for track 2
            SMI_t1 = nan([1,2]);
            SMI_t2 = nan([1,3]);
            spatial_response_t1 = obj.spatial_response{1};
            spatial_response_t2 = obj.spatial_response{2};
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
                    xline([50 90],'LineWidth',1.5,'Color',[146,197,222]/255,'HandleVisibility','off')
                    text(50, y_lim, 'A', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    text(90, y_lim, 'A', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline([70 110],'LineWidth',1.5,'Color',[214,96,77]/255,'HandleVisibility','off')
                    text(70, y_lim, 'B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    text(110, y_lim, 'B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                case 't2'
                    xline(30,'LineWidth',1.5,'HandleVisibility','off')
                   text(30, y_lim, 'Y', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline([50 90],'LineWidth',1.5,'Color',[166,219,160]/255,'HandleVisibility','off')
                    text(50, y_lim, 'C', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    text(90, y_lim, 'C', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline([70 110],'LineWidth',1.5,'Color',[214,96,77]/255,'HandleVisibility','off')
                    text(70, y_lim, 'B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    text(110, y_lim, 'B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline(130,'LineWidth',1.5,'HandleVisibility','off')
                    text(130, y_lim, 'Y', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                case 'both'
                    xline(30,'LineWidth',1.5,'HandleVisibility','off')
                    text(30, y_lim, 'X-Y', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline([50 90],'LineWidth',1.5,'Color',[146,197,222]/255,'HandleVisibility','off')
                    text(50, y_lim, 'A-C', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    text(90, y_lim, 'A-C', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline([70 110],'LineWidth',1.5,'Color',[214,96,77]/255,'HandleVisibility','off')
                    text(70, y_lim, 'B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    text(110, y_lim, 'B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                    xline(130,'LineWidth',1.5,'LineStyle','--','HandleVisibility','off')
                    text(130, y_lim, 'Y', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    
            end
        end
    end
end