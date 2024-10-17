%% Spatial raster plot and spatial tuning curves & spatial modulation analysis

clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'RUN2';

for nsession = 1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);

        if contains(stimulus_name{n},'Masa2tracks')
%             load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'))
            session_clusters_RUN1 = session_clusters;
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'))
            session_clusters_RUN2 = session_clusters;
%             load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_original%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        end

        % reconstruct cluster structure and place field structure from
        % session clusters

        clusters_combined1= session_clusters_RUN1;
        clusters_combined1.spike_id=vertcat(session_clusters_RUN1.spike_id{:});
        clusters_combined1.spike_times=vertcat(session_clusters_RUN1.spike_times{:});
        [clusters_combined1.spike_times,index] =sort(clusters_combined1.spike_times);
        clusters_combined1.spike_id=clusters_combined1.spike_id(index);
        clusters_combined=clusters_combined1;

        clusters_combined2= session_clusters_RUN2;
        clusters_combined2.spike_id=vertcat(session_clusters_RUN2.spike_id{:});
        clusters_combined2.spike_times=vertcat(session_clusters_RUN2.spike_times{:});
        [clusters_combined2.spike_times,index] =sort(clusters_combined2.spike_times);
        clusters_combined2.spike_id=clusters_combined2.spike_id(index);
        

        % Cell with spatial tuning
        % ib = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
        %     | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));
        % ia = find(contains(clusters_combined.region,'HPC'))
        %         [C,ia,ic] = unique(clusters_combined.cluster_id);
        ia = find(clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95);
%         ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
%             | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'HPC'));
        C = clusters_combined1.cluster_id(ia);


        ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'V1_L'));

        speed = clusters_combined.speed{1};
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
        lick_speed = interp1(clusters_combined.sglxTime_uncorrected{1},speed,clusters_combined.lick_time{1},'nearest');

        x_window = [0 140];
        x_bin_size =2;
        place_fields_all = calculate_spatial_cells(clusters_combined,clusters_combined.tvec{1},...
            clusters_combined.position{1},clusters_combined.speed{1},clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_window,x_bin_size);

        ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'V1_L'));
        plot_place_cell_map_correlation(place_fields_all,ia,[],[],options)

       
        ratemap_matrix=[];
        average_maps=[];
        lap_correlation_pv=[];
        within_track_corr=[];
        lap_correlation=[];
        for track_id = 1:2
            %             ratemap_matrix{track_id} = [conv(place_fields_all(track_id).raw{good_cell_index}gaussianWindow,'same')];
            ratemap_matrix{track_id} = [place_fields_all(track_id).raw{ia}];

            within_track_corr{track_id} = [place_fields_all(track_id).within_track_corr{ia}];
            %             ratemaps_track1(nevent,:) = conv(ratemaps_track1(nevent,:),gaussianWindow,'same');
            within_track_corr{track_id} = reshape( within_track_corr{track_id},size(place_fields_all(track_id).within_track_corr{1},1),size(place_fields_all(track_id).within_track_corr{1},1),[]);%laps X position bins X cells
            lap_correlation{track_id}= squeeze(nanmean(within_track_corr{track_id},3));

%             for nlap = 1:size(place_fields_all(track_id).raw{1},1)
%                 for mlap = 1:size(place_fields_all(track_id).raw{1},1)
%                     lap_correlation_pv(track_id,nlap,mlap) = nanmean(nanmean(corr(squeeze(ratemap_matrix{track_id}(nlap,:,:)),squeeze(ratemap_matrix{track_id}(mlap,:,:)))));
%                 end
%             end
% 
% 
% 
%             for ncell = 1:length(ia)
%                 for nlap = 1:size(place_fields_all(track_id).raw{1},1)
%                     all_cells_corr(ncell,nlap) = corr(zscore(place_fields_all(track_id).raw{ia(ncell)}(nlap,:))',average_maps{track_id}(:,ncell));
%                 end
%             end
        end
        
        figure
        subplot(2,2,1)
        imagesc(lap_correlation{1}')
        title('Track 1')
        colorbar
        subplot(2,2,2)
        imagesc(lap_correlation{2}')
        colorbar
        title('Track 2')
        
        mean([place_fields_all(track_id).with_track_corr{ia}])
        %         ratemap_matrix = [place_fields_all(track_id).raw{place_fields(1).all_good_cells_LIBERAL}];
        %         ratemap_matrix = reshape(ratemap_matrix,size(place_fields_all(track_id).raw{1},1),[],length(place_fields(1).all_good_cells_LIBERAL));%laps X position bins X cells
        temp_map = normalize([average_maps{1}; average_maps{2}],'range');

%         figure
%         for iCell =1:49
%             subplot(7,7,iCell)
%             imagesc(clusters_combined1.within_track_corr{ia(iCell),1})
%             title(sprintf('%i',clusters_combined1.cluster_id(ia(iCell))))
%         end

        figure
        for iCell =1:49
            subplot(7,7,iCell)
            imagesc(clusters_combined1.spatial_response{ia(iCell),1})
            title(sprintf('%i',clusters_combined1.cluster_id(ia(iCell))))
        end
        sgtitle('RUN1 Track L')

        figure
        for iCell =1:49
            subplot(7,7,iCell)
            imagesc(clusters_combined2.spatial_response{ia(iCell),1})
            title(sprintf('%i',clusters_combined2.cluster_id(ia(iCell))))
        end
       sgtitle('RUN2 Track L')

        figure
        for iCell =1:49
            subplot(7,7,iCell)
            imagesc(clusters_combined1.spatial_response{ia(iCell),2})
            title(sprintf('%i',clusters_combined1.cluster_id(ia(iCell))))
        end
        sgtitle('RUN1 Track R')

        figure
        for iCell =1:49
            subplot(7,7,iCell)
            imagesc(clusters_combined2.spatial_response{ia(iCell),2})
            title(sprintf('%i',clusters_combined2.cluster_id(ia(iCell))))
        end
       sgtitle('RUN2 Track R')
%         figure
%         for iCell =1:49
%             subplot(7,7,iCell)
%             imagesc(clusters_combined2.within_track_corr{ia(iCell),1})
%             title(sprintf('%i',clusters_combined2.cluster_id(ia(iCell))))
%         end


        clear place_fields
        x_bin_size =2;
        spatial_response = calculate_raw_spatial_response(clusters_combined.spike_id,clusters_combined.cluster_id,clusters_combined.spike_times,clusters_combined.tvec{1},...
            clusters_combined.position{1},clusters_combined.speed{1},clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_bin_size);
        for track_id = 1:max(session_clusters.track_ID_all{1})
            place_fields(track_id).x_bin_edges = 0:x_bin_size:140;
            place_fields(track_id).x_bin_centres = x_bin_size/2:x_bin_size:140-x_bin_size/2;
            place_fields(track_id).raw = spatial_response(:,track_id);
        end


        Task_info.start_time_all = clusters_combined.start_time_all{1};
        Task_info.end_time_all = clusters_combined.end_time_all{1};
        Task_info.track_ID_all = clusters_combined.track_ID_all{1};

        Behaviour.tvec = clusters_combined.tvec{1};
        Behaviour.position = clusters_combined.position{1};
        Behaviour.speed = clusters_combined.speed{1};
%         place_fields=struct();
%         Behaviour.tvec =Behaviour.sglxTime_uncorrected ;% check if photodiode correction is causing this issue...
        plot_raster_both_track(clusters_combined.spike_times,clusters_combined.spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
            'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C);
%    plot_raster_both_track(clusters_combined.spike_times,clusters_combined.spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
%             'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'place_fields',place_fields);
        if  contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH',sprintf(erase(stimulus_name{n},'Masa2tracks_'))))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH',sprintf(erase(stimulus_name{n},'Masa2tracks_'))),[])
        else
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','Spatial PSTH'),[])
        end


        %
        %         % Spatial modulation
        %         x_bin_size = mean(diff(place_fields_V1_L(1).x_bin_centres));
        %         SMI = calculate_spatial_modulation_index(V1_clusters_L,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_V1_L,'subplot_xy',[3 1],'plot_option',1)
        %         SMI = calculate_spatial_modulation_index(V1_clusters_R,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_V1_R,'subplot_xy',[3 1],'plot_option',1)
        %
        %
        %         metric_param.region = @(x) contains(x,'HPC_L');
        %         [HPC_clusters_L,cluster_id] = select_clusters(clusters(1),metric_param);
        %
        %         SMI = calculate_spatial_modulation_index(HPC_clusters_L,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_HPC_L,'subplot_xy',[3 1],'plot_option',1)
        %         SMI = calculate_spatial_modulation_index(HPC_clusters_R,Task_info,Behaviour,[0 140],x_bin_size,'place_fields',place_fields_HPC_R,'subplot_xy',[3 1],'plot_option',1)
        %         %         calculate_spatial_modulation_index(place_fields_V1_L);

        spatial_cell_id = find((clusters_combined.peak_percentile(:,1)>0.95&clusters_combined.odd_even_stability(:,1)>0.95) ...
            | (clusters_combined.peak_percentile(:,2)>0.95&clusters_combined.odd_even_stability(:,2)>0.95));

        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);

            % plot populational map and PV correlation
            %%%%%% HPC
            options.region = 'HPC';
            if options.probe_hemisphere==1
                cluster_id=intersect(spatial_cell_id,find(contains(clusters_combined.region,'HPC_L')));
            else
                cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'HPC_R')));
            end

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if options.probe_hemisphere == 1
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_HPC_L%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_L.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            elseif options.probe_hemisphere == 2
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_HPC_R%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_R.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            end

            % plot populational map and PV correlation
            %%%%%% V1
            options.region = 'V1';

            if options.probe_hemisphere==1
                cluster_id=intersect(spatial_cell_id,find(contains(clusters_combined.region,'V1_L')));
            else
                cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'V1_R')));
            end
            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if options.probe_hemisphere == 1
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_V1_L%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_L.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            elseif options.probe_hemisphere == 2
                if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_V1_R%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                else
                    save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_R.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
                end
            end

        end

        if length(session_info(n).probe) > 1
            % plot populational map and PV correlation
            %%%%%% HPC
            options.probe_combined = 1;

            options.region = 'HPC';
            cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'HPC')))

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_HPC_combined%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_HPC_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            end

            % plot populational map and PV correlation
            %%%%%% V1
            options.region = 'V1';
            cluster_id= intersect(spatial_cell_id,find(contains(clusters_combined.region,'V1')))

            [~,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = ...
                plot_place_cell_map_correlation(place_fields,cluster_id,Task_info,Behaviour,options); % Roughly 6-7 mins for shuffle and plotting

            if contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
                save(fullfile(options.ANALYSIS_DATAPATH,sprintf('population_vector_corr_V1_combined%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'population_vector_corr_V1_combined.mat'),'PPvector','shuffled_globalRemap_PPvector','shuffled_rateRemap_PPvector');
            end
            options = rmfield(options,'probe_combined');
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'))== 0
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'))
        end
        
        if  contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map',sprintf(erase(stimulus_name{n},'Masa2tracks_'))))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map',sprintf(erase(stimulus_name{n},'Masa2tracks_'))),[])
        else
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','populational_map'),[])
        end

        close all

        %         %GLM analysis
        %         spatial_modulation_GLM_analysis(V1_clusters_L,[],Behaviour,Task_info);
        %
        %         spatial_modulation_GLM_analysis(V1_clusters_L,place_fields_V1_L(1).cluster_id(place_fields_V1_L(1).all_good_cells_LIBERAL),Behaviour,Task_info);
    end
end