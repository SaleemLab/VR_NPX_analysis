%% Main pipeline for viewing events
%% UP DOWN state and ripple and spindle analysis

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'Sleep';
experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
        if isempty(DIR)
            continue
        end

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks3;
        elseif contains(stimulus_name{n},'Sleep')
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
%             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
%             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            clusters=clusters_ks3;
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            clusters=clusters_ks3;
        end
        clear CA1_clusters V1_clusters 
        
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        DIR2 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'));
        
        DIR3 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters.mat'));
        
        session_clusters_RUN=[];
        session_clusters_RUN1=[];
        session_clusters_RUN2=[];
        
        if ~isempty(DIR)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            session_clusters_RUN=session_clusters;
        end

        if ~isempty(DIR1)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
            session_clusters_RUN1=session_clusters;
            session_clusters_RUN=session_clusters_RUN1;
        end

        if ~isempty(DIR2)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN2.mat'));
            session_clusters_RUN2=session_clusters;
        end
        
        session_clusters=[];
        if ~isempty(DIR3)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters.mat'));
        end

        % From cell structure back to spike times and spike id
        session_clusters_RUN.spike_id=vertcat(session_clusters_RUN.spike_id{:});
        session_clusters_RUN.spike_times=vertcat(session_clusters_RUN.spike_times{:});
        [session_clusters_RUN.spike_times,index] =sort(session_clusters_RUN.spike_times);
        session_clusters_RUN.spike_id=session_clusters_RUN.spike_id(index);

        params = create_cluster_selection_params('sorting_option','masa');
        clear selected_clusters
        for nprobe = 1:length(clusters)
            selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters
        end

        for nprobe = 1:length(clusters)
            % Convert to unique spike/cluster id
            selected_clusters(nprobe).spike_id = selected_clusters(nprobe).spike_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
            selected_clusters(nprobe).cluster_id = selected_clusters(nprobe).cluster_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
            if session_info(n).probe(nprobe).probe_hemisphere==1
                selected_clusters(nprobe).region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'L'));
                selected_clusters(nprobe).spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'L'),:);
                selected_clusters(nprobe).odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'L'),:);
                selected_clusters(nprobe).peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'L'),:);
            elseif session_info(n).probe(nprobe).probe_hemisphere==2
                selected_clusters(nprobe).region = session_clusters_RUN.region(contains(session_clusters_RUN.region,'R'));
                selected_clusters(nprobe).spatial_response = session_clusters_RUN.spatial_response(contains(session_clusters_RUN.region,'R'),:);
                selected_clusters(nprobe).odd_even_stability = session_clusters_RUN.odd_even_stability(contains(session_clusters_RUN.region,'R'),:);
                selected_clusters(nprobe).peak_percentile = session_clusters_RUN.peak_percentile(contains(session_clusters_RUN.region,'R'),:);
            end
        end

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(selected_clusters(1),selected_clusters(2));
        else
            clusters_combined = selected_clusters;
        end
        clear selected_clusters
%         spatial_cell_index = find((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
%             | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));
        spatial_cell_index = find(session_clusters_RUN.odd_even_stability(:,1)>0.95 ...
            | session_clusters_RUN.odd_even_stability(:,2)>0.95);

        metric_param =[];
        metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));
        [selected_clusters,cluster_id] = select_clusters(session_clusters_RUN,metric_param);
        
        x_window = [0 140];
        x_bin_width = 2;
        place_fields = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
            selected_clusters.position{1},selected_clusters.speed{1},selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);

        
        for nprobe = 1:2
            if ~isempty(behavioural_state(nprobe).SWS)
                [V1_reactivations(nprobe).SWS_offset,V1_reactivations(nprobe).SWS_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state(nprobe).SWS);
                V1_reactivations(nprobe).SWS_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).SWS_index)';
            end
        end

        [reactivations_combined.SWS_offset,reactivations_combined.SWS_index] = RestrictInts(reactivations_combined.offset',behavioural_state(1).SWS);
        reactivations_combined.SWS_onset = reactivations_combined.onset(reactivations_combined.SWS_index)';

        cortex_LFP=[];
        CA1_LFP=[];
        for nprobe = 1:length(session_clusters_RUN.probe_hemisphere)
            probe_no=session_clusters_RUN.probe_hemisphere(nprobe);

            if isfield(LFP(probe_no),'L5')                

                if ~isempty(LFP(probe_no).L5)
                    bad_channels=[];
                    all_shanks = 1:size(LFP(probe_no).L5_power,1);
                    for nshank = 1:size(LFP(probe_no).L5_power,1)
                        bad_channels(nshank,:) = LFP(probe_no).L5_power(nshank,:)>3*mean(LFP(probe_no).L5_power(all_shanks~=nshank,:));
                    end
                    bad_channels = sum(bad_channels,2)>4;
                    [~,best_channel]=max(LFP(probe_no).L5_power(~bad_channels,7));
                    good_channels = find(~bad_channels);
                    cortex_LFP{probe_no} = LFP(probe_no).L5(good_channels(best_channel),:);

                else
                    bad_channels=[];
                    all_shanks = 1:size(LFP(probe_no).L4_power,1);
                    for nshank = 1:size(LFP(probe_no).L4_power,1)
                        bad_channels(nshank,:) = LFP(probe_no).L4_power(nshank,:)>3*mean(LFP(probe_no).L4_power(all_shanks~=nshank,:));
                    end
                    bad_channels = sum(bad_channels,2)>4;
                    [~,best_channel]=max(LFP(probe_no).L4_power(~bad_channels,7));
                    good_channels = find(~bad_channels);
                    cortex_LFP{probe_no} = LFP(probe_no).L4(good_channels(best_channel),:);
                end

            elseif isfield(LFP(probe_no),'L4')
                if ~isempty(LFP(probe_no).L4)
                    bad_channels=[];
                    all_shanks = 1:size(LFP(probe_no).L4_power,1);
                    for nshank = 1:size(LFP(probe_no).L4_power,1)
                        bad_channels(nshank,:) = LFP(probe_no).L4_power(nshank,:)>3*mean(LFP(probe_no).L4_power(all_shanks~=nshank,:));
                    end
                    bad_channels = sum(bad_channels,2)>4;
                    [~,best_channel]=max(LFP(probe_no).L4_power(~bad_channels,7));
                    good_channels = find(~bad_channels);
                    cortex_LFP{probe_no} = LFP(probe_no).L4(good_channels(best_channel),:);
                else
                    cortex_LFP = [];
                    disp('cortex LFP is missing')
                end

            elseif isfield(LFP(probe_no),'MEC')

            end

            if isfield(LFP(probe_no),'CA1')
                bad_channels=[];
                all_shanks = 1:size(LFP(probe_no).CA1_power,1);
                for nshank = 1:size(LFP(probe_no).CA1_power,1)
                    bad_channels(nshank,:) = LFP(probe_no).CA1_power(nshank,:)>3*mean(LFP(probe_no).CA1_power(all_shanks~=nshank,:));
                end
                bad_channels = sum(bad_channels,2)>4;
                
                [~,best_channel]=max(LFP(probe_no).CA1_power(~bad_channels,6));
                good_channels = find(~bad_channels);

                CA1_LFP{probe_no} = LFP(probe_no).CA1(good_channels(best_channel),:);
            end
        end

        spatial_cell_index = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
            | clusters_combined.odd_even_stability(:,2)>0.95);

        metric_param =[];
        metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
        metric_param.region = @(x) contains(x,'V1_L');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        V1_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'V1_R');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        V1_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'HPC_L');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'HPC_R');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];

        metric_param.region = @(x) contains(x,'HPC');
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        HPC_spikes{3}=[selected_clusters.spike_id selected_clusters.spike_times];
        group_name_V1 = {'V1_L','V1_R'};
        group_name_HPC = {'HPC_L','HPC_R','HPC'};
        

        % Sleep
        mobility_thresholded = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))])>2000;%1 second movemean windows

        mob_index = find(mobility_thresholded==1);
        for nindex = 1:length(mob_index)
            if nindex<length(mob_index)
                if Behaviour.tvec(mob_index(nindex+1))-Behaviour.tvec(mob_index(nindex))<1 % if immobility less than 1 second, treats that as movement period
                    mobility_thresholded(mob_index(nindex):mob_index(nindex+1))=1;
                end
            end
        end
        %             plot(Behaviour.tvec,Behaviour.mobility);hold on;
        %             plot(Behaviour.tvec,mobility*500000);plot(Behaviour.tvec,abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]))
        %             plot(Behaviour.mobility_zscore)
        mobility_thresholded = interp1(Behaviour.tvec,double(mobility_thresholded),LFP(probe_no).tvec,'linear');
        Behaviour.mobility_zscore = abs([0 diff(movmean(Behaviour.mobility,1/mean(diff(Behaviour.tvec))))]);% Diff of pixel change
        Behaviour.mobility_zscore(isnan(Behaviour.mobility_zscore))=mean(Behaviour.mobility_zscore,'omitnan');
        Behaviour.mobility_zscore=zscore(Behaviour.mobility_zscore);

        Behaviour.mobility = mobility_thresholded;


        nprobe = 2
        figure
        LFP_SR = 1/mean(diff(LFP(nprobe).tvec));
        hold on;
        plot(LFP(nprobe).tvec,CA1_LFP{1},'r')

        figure
        [s,f,t] = stft(cortex_LFP{1},LFP_SR,'Window',2*hann(LFP_SR),'FrequencyRange','onesided');
        index = find(f<50);
        %         sss= abs(s).^2;
        S_mag_smoothed = imgaussfilt(abs(s(index,:)), 10); % Gaussian smoothing with sigma=2
        figure;
        imagesc(t, f(index), S_mag_smoothed);
        
%         set(gca,YScale="log",...
%             YDir="reverse",View=[0 300])
%         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed,1,[]),1) prctile(reshape(S_mag_smoothed,1,[]),99)])
        colorbar
        hold on;

        plot(LFP(nprobe).tvec,CA1_LFP{2}+1000,'b')

        plot(LFP(nprobe).tvec,cortex_LFP{1}+2000,'r')
        hold on;
        plot(LFP(nprobe).tvec,cortex_LFP{2}+3000,'b')
        
        hold on
        for nwin = 1:length(behavioural_state(nprobe).SWS)
            [~,index1]=min(abs(LFP(nprobe).tvec-behavioural_state(nprobe).SWS(nwin,1)));
            [~,index2]=min(abs(LFP(nprobe).tvec-behavioural_state(nprobe).SWS(nwin,2)));
            
            yline(index1,'r')
            yline(index2,'k')

        end


        
        ripples

        [s,f,t] = stft()
        
    end

end

zero_meaned_log_odds = [decoded_ripple_events(1).track(1).replay_events(:).z_log_odds] - mean([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds]);
histogram(zero_meaned_log_odds,100)


% histogram([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds],100)
find(zero_meaned_log_odds<-1)
find(zero_meaned_log_odds>1)

zero_meaned_log_odds = zscore([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds]);
T2_events = find(zero_meaned_log_odds<-1);
figure
count=1;
for event = 1:5:length(T2_events)

    subplot(5,5,count)
    T1_data = decoded_ripple_events(1).track(1).replay_events(T2_events(event));
    T2_data = decoded_ripple_events(1).track(2).replay_events(T2_events(event));

    timebins = T1_data.timebins_edges(1:end-1);
    imagesc(timebins,...
        [1:2*size(T2_data.replay,1)],...
        [T1_data.replay; T2_data.replay])
    hold on
    colormap(flipud(gray))
    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/5)
    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])


    [b,index]=min(abs(T1_data.onset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    [b,index]=min(abs(T1_data.offset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    count = count + 1;
end


T1_events = find(zero_meaned_log_odds>1);
figure
count = 1;
for event = 1:length(T1_events)
    subplot(5,5,count)
    T1_data = decoded_ripple_events(1).track(1).replay_events(T1_events(event));
    T2_data = decoded_ripple_events(1).track(2).replay_events(T1_events(event));

    timebins = T1_data.timebins_edges(1:end-1);
    imagesc(timebins,...
        [1:2*size(size(T2_data.replay,1))],...
        [T1_data.replay; T2_data.replay])
    colormap(flipud(gray))
    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])

    hold on
    [b,index]=min(abs(T1_data.onset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    [b,index]=min(abs(T1_data.offset-timebins));
    xline(timebins(index)-mean(diff(timebins)/2))
    count = count + 1;
end




