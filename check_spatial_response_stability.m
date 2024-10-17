%% Check spatial response stability

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
            clusters_combined.position{1},speed,clusters_combined.track_ID_all{1},clusters_combined.start_time_all{1},clusters_combined.end_time_all{1},x_window,x_bin_size);
        %
        %         ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
        %             | clusters_combined1.odd_even_stability(:,2)>0.95)&contains(clusters_combined1.region,'V1_L'));
        %         plot_place_cell_map_correlation(place_fields_all,ia,[],[],options)

        ia = find((clusters_combined1.odd_even_stability(:,1)>0.95 ...
            | clusters_combined1.odd_even_stability(:,2)>0.95));
        ratemap_matrix=[];
        average_maps=[];
        lap_correlation_pv=[];
        within_track_corr=[];
        lap_correlation=[];
        pv_corr=[];
        for track_id = 1:2
            %             ratemap_matrix{track_id} = [conv(place_fields_all(track_id).raw{good_cell_index}gaussianWindow,'same')];
            ratemap_matrix{track_id} = [place_fields_all(track_id).raw{ia}];
            ratemap_matrix{track_id} = reshape( ratemap_matrix{track_id},size(place_fields_all(track_id).within_track_corr{1},1),[],length(ia));%laps X position bins X cells


            within_track_corr{track_id} = [place_fields_all(track_id).within_track_corr{ia}];
            %             ratemaps_track1(nevent,:) = conv(ratemaps_track1(nevent,:),gaussianWindow,'same');
            within_track_corr{track_id} = reshape( within_track_corr{track_id},size(place_fields_all(track_id).within_track_corr{1},1),size(place_fields_all(track_id).within_track_corr{1},1),[]);
            lap_correlation{track_id}= squeeze(nanmean(within_track_corr{track_id},3));

            for nlap = 1:size(place_fields_all(track_id).raw{1},1)
                for mlap = 1:size(place_fields_all(track_id).raw{1},1)
                    map1 = zscore(squeeze(ratemap_matrix{track_id}(nlap,:,:)));
                    map2 = zscore(squeeze(ratemap_matrix{track_id}(mlap,:,:)));
                    temp = corr(map1',map2');
                    x_bins = 1:size(map1,1);
                   
                    lap_correlation_pv{track_id}(nlap,mlap) = nanmean(diag(temp));
%                     lap_correlation_pv{track_id}(nlap,mlap) = nanmean(nanmean(corr(map1,map2)));
                end
            end

            for nlap = 1:size(place_fields_all(track_id).raw{1},1)
                map1 = zscore(squeeze(ratemap_matrix{track_id}(nlap,:,:)));
                map2 = zscore(squeeze(nanmean(ratemap_matrix{track_id}(:,:,:))));
                temp = corr(map1',map2');
                x_bins = 1:size(map1,1);
                
                pv_corr(track_id,nlap) = nanmean(diag(temp));
            end
            %
            %
            %
            %             for ncell = 1:length(ia)
            %                 for nlap = 1:size(place_fields_all(track_id).raw{1},1)
            %                     all_cells_corr(ncell,nlap) = corr(zscore(place_fields_all(track_id).raw{ia(ncell)}(nlap,:))',average_maps{track_id}(:,ncell));
            %                 end
            %             end
        end

        fig = figure
        fig.Position = [500 70 830 600];
        sgtitle(sprintf('%s %s %s cell map remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s %s cell map remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});

        subplot(3,2,1)
        percentile_95 = prctile(reshape(lap_correlation{1},1,[]),95);
        imagesc(lap_correlation{1}')
        clim([0 percentile_95])
        title('Track 1 mean lap by lap single cell corr')
        colorbar

        subplot(3,2,2)
        percentile_95 = prctile(reshape(lap_correlation{2},1,[]),95);
        imagesc(lap_correlation{2}')
        clim([0 percentile_95])
        colorbar
        title('Track 2 mean lap by lap single cell corr')

        
        subplot(3,2,3)
        percentile_95 = prctile(reshape(lap_correlation_pv{1},1,[]),95);
        imagesc(lap_correlation_pv{1}')
        clim([0 percentile_95])
        title('Track 1 mean lap to lap PV corr')
        colorbar
        
        subplot(3,2,4)
        percentile_95 = prctile(reshape(lap_correlation_pv{2},1,[]),95);
        imagesc(lap_correlation_pv{2}')
        clim([0 percentile_95])
        colorbar
        title('Track 2 mean lap to lap PV corr')

        subplot(3,2,5)
%         percentile_95 = prctile(reshape(pv_corr{1},1,[]),95);
        plot(movmedian(pv_corr(1,:)',3))
%       ylim([-percentile_95 percentile_95])
        title('Track 1 mean lap to template PV corr')
        box off
%         colorbar

        subplot(3,2,6)
%         percentile_95 = prctile(reshape(pv_corr(2,:),1,[]),95);
        plot(movmedian(pv_corr(2,:)',3))
%         ylim([-percentile_95 percentile_95])
        title('Track 2 mean lap to template PV corr')
        box off
%         colorbar

        
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