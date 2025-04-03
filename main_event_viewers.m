%% Main pipeline for viewing sleep events

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
            
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks4;
        elseif contains(stimulus_name{n},'Sleep')
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
%             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
%             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
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
        
        % DIR3 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_Sleep.mat'));

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
        if contains(stimulus_name{n},'Sleep')
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
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

        %         for nprobe = 1:2
        %             if ~isempty(behavioural_state(nprobe).SWS)
        %                 [V1_reactivations(nprobe).SWS_offset,V1_reactivations(nprobe).SWS_index] = RestrictInts(V1_reactivations(nprobe).offset',behavioural_state(nprobe).SWS);
        %                 V1_reactivations(nprobe).SWS_onset = V1_reactivations(nprobe).onset(V1_reactivations(nprobe).SWS_index)';
        %             end
        %         end
        %
        %         [reactivations_combined.SWS_offset,reactivations_combined.SWS_index] = RestrictInts(reactivations_combined.offset',behavioural_state(1).SWS);
        %         reactivations_combined.SWS_onset = reactivations_combined.onset(reactivations_combined.SWS_index)';

        cortex_LFP=[];
        CA1_LFP=[];
        cortex_LFP_filtered=[];
        CA1_LFP_ripple_filtered=[];
        CA1_LFP_filtered=[];
        LFP_SR = 1/mean(diff(LFP(nprobe).tvec));
        for nprobe = 1:length(session_clusters_RUN.probe_hemisphere)
            probe_no=session_clusters_RUN.probe_hemisphere(nprobe);

            if isfield(LFP(probe_no),'best_V1')                

                if ~isempty(LFP(probe_no).best_V1)
                    bad_channels=[];
                    all_shanks = 1:size(LFP(probe_no).best_V1_power,1);
                    for nshank = 1:size(LFP(probe_no).best_V1_power,1)
                        bad_channels(nshank,:) = LFP(probe_no).best_V1_power(nshank,:)>3*mean(LFP(probe_no).best_V1_power(all_shanks~=nshank,:));
                    end
                    bad_channels = sum(bad_channels,2)>4;
                    [~,best_channel]=max(LFP(probe_no).best_V1_power(~bad_channels,1));
                    good_channels = find(~bad_channels);
                    cortex_LFP{probe_no} = LFP(probe_no).best_V1(good_channels(best_channel),:);

                else
                    bad_channels=[];
                    all_shanks = 1:size(LFP(probe_no).L4_power,1);
                    for nshank = 1:size(LFP(probe_no).L4_power,1)
                        bad_channels(nshank,:) = LFP(probe_no).L4_power(nshank,:)>3*mean(LFP(probe_no).L4_power(all_shanks~=nshank,:));
                    end
                    bad_channels = sum(bad_channels,2)>4;
                    [~,best_channel]=max(LFP(probe_no).L4_power(~bad_channels,1));
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
                    [~,best_channel]=max(LFP(probe_no).L4_power(~bad_channels,1));
                    good_channels = find(~bad_channels);
                    cortex_LFP{probe_no} = LFP(probe_no).L4(good_channels(best_channel),:);
                else
                    cortex_LFP = [];
                    disp('cortex LFP is missing')
                end

            elseif isfield(LFP(probe_no),'MEC')

            end

            passband = [9 17];
            filter_type  = 'bandpass';
            filter_order = round(6*LFP_SR/(max(passband)-min(passband)));  % creates filter for ripple
            norm_freq_range = passband/(LFP_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
            b_ripple = fir1(filter_order, norm_freq_range,filter_type);
            cortex_LFP_filtered{probe_no} = filtfilt(b_ripple,1,cortex_LFP{probe_no});
            %         cortex_LFP{probe_no} = zscore(abs(hilbert(signal)));

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
                
                passband = [0.5 300];
                filter_type  = 'bandpass';
                filter_order = round(6*LFP_SR/(max(passband)-min(passband)));  % creates filter for ripple
                norm_freq_range = passband/(LFP_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                b_ripple = fir1(filter_order, norm_freq_range,filter_type);
                signal = filtfilt(b_ripple,1,CA1_LFP{probe_no});
                %         zscored_ripple = zscore(abs(hilbert(signal)));
                CA1_LFP_filtered{probe_no}=signal;

                passband = [150 300];
                filter_type  = 'bandpass';
                filter_order = round(6*LFP_SR/(max(passband)-min(passband)));  % creates filter for ripple
                norm_freq_range = passband/(LFP_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                b_ripple = fir1(filter_order, norm_freq_range,filter_type);
                signal = filtfilt(b_ripple,1,CA1_LFP{probe_no});
                %         zscored_ripple = zscore(abs(hilbert(signal)));
                CA1_LFP_ripple_filtered{probe_no}=signal;
            end
        end

       
        spatial_cell_index = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
            | clusters_combined.odd_even_stability(:,2)>0.95);

        metric_param =[];
%         metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
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
        

        %%%% Sleep
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

        %%% Get spike counts

        tvec = LFP(1).tvec;
        tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];
        CA1_spike_counts=[];
        V1_spike_counts=[];
        w = gausswin(0.02*1/mean(diff(tvec)));
        w = w / sum(w);

        % [NREM_status,interval,~] = InIntervals(tvec,behavioural_state_merged.SWS);
        % spike_times_sleep = NREM_status;

        for nprobe = 1:length(V1_spikes)
            spike_times = V1_spikes{nprobe}(:,2);
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            spike_times_sleep = spike_times;
            V1_spike_counts{nprobe} = filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';

            spike_times = HPC_spikes{nprobe}(:,2);
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            spike_times_sleep = spike_times;
            CA1_spike_counts{nprobe} = filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';
        end
        

        %%%%%% TF plot
        nprobe = 2
        LFP_SR = 1/mean(diff(LFP(nprobe).tvec));

        [s,f,t] = stft(cortex_LFP{1},LFP_SR,'Window',hann(LFP_SR),'FrequencyRange','onesided','FFTLength',10*round(LFP_SR));
        %         index = find(f<20);
        index = find(f>1 & f<20);
        %         sss= abs(s).^2;
        S_mag_smoothed_V1 = imgaussfilt(abs(s(index,:)), 10); % Gaussian smoothing with sigma=2

        [s,f,t] = stft(CA1_LFP{1},LFP_SR,'Window',hann(LFP_SR),'FrequencyRange','onesided','FFTLength',10*round(LFP_SR));
        index = find(f>1 & f<20);
        %         sss= abs(s).^2;
        S_mag_smoothed_HPC = imgaussfilt(abs(s(index,:)), 10); % Gaussian smoothing with sigma=2




        tindex = find(t>behavioural_state(1).REM(1,1)-100 & t<behavioural_state(1).REM(1,1)+200);

        figure
        subplot(3,1,1)
        imagesc(t(tindex), f(index), S_mag_smoothed_V1(:,tindex));

        %         set(gca,YScale="log",...
        %             YDir="reverse",View=[0 300])
        %         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed,1,[]),1) prctile(reshape(S_mag_smoothed,1,[]),99)])
        colorbar
        colormap(flipud(bone))
        hold on

        yline(f(min(find(f>1&f<3))),'r')
        yline(f(max(find(f>1&f<3))),'r')
        title('V1 LFP Time frequency plot')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(3,1,2)
        imagesc(t(tindex), f(index), S_mag_smoothed_HPC(:,tindex));

        %         set(gca,YScale="log",...
        %             YDir="reverse",View=[0 300])
        %         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed_HPC,1,[]),1) prctile(reshape(S_mag_smoothed_HPC,1,[]),99)])
        colorbar
        colormap(flipud(bone))
        hold on

        yline(f(min(find(f>4&f<10))),'r')
        yline(f(max(find(f>4&f<10))),'r')

        xline(min(t(t>behavioural_state(1).REM(1,1))),'k')
        xline(max(t(t<behavioural_state(1).REM(1,2))),'k')
        title('HPC LFP Time frequency plot')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(3,1,3)
        interp1_mobility = interp1(Behaviour.tvec,Behaviour.mobility_zscore,t(tindex));
        plot(t(tindex),interp1_mobility)
        xlabel('Time(s)')
        ylabel('Mobility zscored')
        for nevent =1:size(behavioural_state(1).SWS,1)
            xline(min(t(t>behavioural_state(1).SWS(nevent,1))),'r')
            xline(max(t(t<behavioural_state(1).SWS(nevent,2))),'b')
        end
        xlim([min(t(tindex)) max(t(tindex))])
        colorbar
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        %%%%%% Sleep TF plot
        
        tindex = 1:length(t);

        figure
        subplot(3,1,1)
        imagesc(t(tindex), f(index), S_mag_smoothed_V1(:,tindex));

        %         set(gca,YScale="log",...
        %             YDir="reverse",View=[0 300])
        %         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed,1,[]),1) prctile(reshape(S_mag_smoothed,1,[]),99)])
        colorbar
        colormap(flipud(bone))
        hold on

        yline(f(min(find(f>1&f<3))),'r')
        yline(f(max(find(f>1&f<3))),'r')
        xlabel('Time(s)')
        ylabel('Frqeuncy')
        title('V1 LFP Time frequency plot')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(3,1,2)
        imagesc(t(tindex), f(index), S_mag_smoothed_HPC(:,tindex));

        %         set(gca,YScale="log",...
        %             YDir="reverse",View=[0 300])
        %         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed_HPC,1,[]),1) prctile(reshape(S_mag_smoothed_HPC,1,[]),99)])
        colorbar
        colormap(flipud(bone))
        hold on

        yline(f(min(find(f>4&f<10))),'r')
        yline(f(max(find(f>4&f<10))),'r')
%         for nevent =1:size(behavioural_state(1).REM,1)
%             xline(min(t(t>behavioural_state(1).REM(nevent,1))),'r')
%             xline(max(t(t<behavioural_state(1).REM(nevent,2))),'b')
%         end
        xlabel('Time(s)')
        ylabel('Frqeuncy')
        title('HPC LFP Time frequency plot')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(3,1,3)
        interp1_mobility = interp1(Behaviour.tvec,Behaviour.mobility_zscore,t(tindex));
        plot(t(tindex),interp1_mobility)
        xlabel('Time(s)')
        ylabel('Mobility zscored')
        for nevent =1:size(behavioural_state(1).SWS,1)
            xline(min(t(t>behavioural_state(1).SWS(nevent,1))),'r')
            xline(max(t(t<behavioural_state(1).SWS(nevent,2))),'b')
        end
        xlim([min(t(tindex)) max(t(tindex))])
        colorbar
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


      

        %%%%%% peri ripple LFP and spiking
        for nprobe = 1:2
            [ripples(nprobe).DOWN_UP_transition_offset,ripples(nprobe).DOWN_UP_transition_index] = RestrictInts(ripples(nprobe).offset,[slow_waves(1).UP_ints(:,1) slow_waves(1).UP_ints(:,1)+0.2]);
            ripples(nprobe).DOWN_UP_transition_onset = ripples(nprobe).onset(ripples(nprobe).DOWN_UP_transition_index)';
        end

        for nprobe = 1:2
            [ripples(nprobe).UP_DOWN_transition_offset,ripples(nprobe).UP_DOWN_transition_index] = RestrictInts(ripples(nprobe).offset,[slow_waves(1).UP_ints(:,2)-0.4 slow_waves(1).UP_ints(:,2)+0.1]);
            ripples(nprobe).UP_DOWN_transition_onset = ripples(nprobe).onset(ripples(nprobe).UP_DOWN_transition_index)';
        end

% slow_waves = [];

        event_times = [spindles(1).onset(spindles(1).peak_zscore>4) spindles(1).offset(spindles(1).peak_zscore>4)];
        %         tindex = find(tvec>ripples(1).UP_DOWN_transition_onset(min(nevent))-0.5 & tvec<ripples(1).UP_DOWN_transition_onset(max(nevent))+0.5);
        tindex = find(tvec>event_times(min(nevent),1)-0.1 & tvec<event_times(max(nevent),2)+0.1);

        count = 1;
        nevent =1:48;
        for count = 1:48
            subplot(7,7,count)
            tindex = find(tvec>event_times(nevent(count),1)-0.2 & tvec<event_times(nevent(count),2)+0.2);
            plot(tvec(tindex),cortex_LFP_filtered{1}(tindex),'k','LineWidth',1.5)
            ylim([min(cortex_LFP_filtered{1}(tindex))-50 max(cortex_LFP_filtered{1}(tindex))+50])
        end




        nevent =200
        event_times = [ripples(1).UP_DOWN_transition_onset' ripples(1).UP_DOWN_transition_offset];
        %         tindex = find(tvec>ripples(1).UP_DOWN_transition_onset(min(nevent))-0.5 & tvec<ripples(1).UP_DOWN_transition_onset(max(nevent))+0.5);
        tindex = find(tvec>event_times(min(nevent),1)-0.1 & tvec<event_times(max(nevent),2)+0.1);


        count = 1;
        nevent =200:225;
        for count = 1:25
            subplot(5,5,count)
            tindex = find(tvec>event_times(nevent(count),1)-0.05 & tvec<event_times(nevent(count),2)+0.05);
            plot(tvec(tindex),LFP(1).best_HPC(1,tindex),'k','LineWidth',1.5)
            ylim([min(LFP(1).best_HPC(1,tindex))-50 max(LFP(1).best_HPC(1,tindex))+50])
        end

        figure
        subplot(3,1,1)
        plot(tvec(tindex),2*cortex_LFP_filtered{1}(tindex)+1500,'r')
        hold on
        plot(tvec(tindex),2*cortex_LFP_filtered{2}(tindex)+2500,'b')
        hold on;
        plot(tvec(tindex),LFP(1).best_HPC(1,tindex)+800,'k')
        % plot(tvec(tindex),CA1_LFP_filtered{1}(tindex)+500,'k')
        plot(tvec(tindex),CA1_LFP_ripple_filtered{1}(tindex)*2,'k')
  
        for n = nevent
            %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
            %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
            %             xline(min(tvec(tvec>event_times(n,2))),'b')
            %             xline(max(tvec(tvec<event_times(n,1))),'r')
            rectangle('Position',[min(tvec(tvec>event_times(n,1))) -0.5 ...
                min(tvec(tvec>event_times(n,2)))-max(tvec(tvec<event_times(n,1))),...
                2500],...
                'FaceColor',[1 0 0 0.2],'EdgeColor','none')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(3,1,2)
        %         CA1_spike_counts{3} = zscore(CA1_spike_counts{1}+CA1_spike_counts{2});
        %         V1_spike_counts{3} = zscore(V1_spike_counts{1}+V1_spike_counts{2});

        plot(tvec(tindex),2*V1_spike_counts{1}(tindex)+2,'r')
        hold on
        plot(tvec(tindex),2*V1_spike_counts{2}(tindex),'b')
        title('V1 spiking')
        for n = nevent
            %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
            %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
            %             xline(min(tvec(tvec>event_times(n,2))),'b')
            %             xline(max(tvec(tvec<event_times(n,1))),'r')
            rectangle('Position',[min(tvec(tvec>event_times(n,1))) -0.5 ...
                min(tvec(tvec>event_times(n,2)))-max(tvec(tvec<event_times(n,1))),...
                2+max(V1_spike_counts{1}(tindex))],...
                'FaceColor',[1 0 0 0.2],'EdgeColor','none')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(3,1,3)
        plot(tvec(tindex),CA1_spike_counts{1}(tindex)+3,'r')
        hold on
        plot(tvec(tindex),CA1_spike_counts{2}(tindex),'b')
        title('HPC spiking')
        for n = nevent
            %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
            %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
            rectangle('Position',[min(tvec(tvec>event_times(n,1))) -0.5 ...
                min(tvec(tvec>event_times(n,2)))-max(tvec(tvec<event_times(n,1))),...
                3+max(CA1_spike_counts{1}(tindex))],...
                'FaceColor',[1 0 0 0.2],'EdgeColor','none')

            %             xline(min(tvec(tvec>event_times(n,2))),'b')
            %             xline(max(tvec(tvec<event_times(n,1))),'r')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



        
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

%% 
figure
nexttile
%         plot_perievent_event_histogram(ripples(2).SWS_peaktimes,ripples(1).SWS_peaktimes,'twin',[-1 1],'event_name','Right ripple')
time_wondows = [-2 2];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes,ripples(2).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
title('Probability of Left ripple during Right HPC')
ylabel('probability')
xlabel('Time(s)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
%         plot_perievent_event_histogram(slow_waves(2).ints.UP(:,1),slow_waves(1).ints.UP(:,1),'twin',[-1 1],'event_name','Right UP events')
time_wondows = [-2 2];
time_bin = 0.02;
probabilities = calculate_event_probability(slow_waves(1).ints.UP(:,1),slow_waves(2).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
title('Probability of UP state from Left V1 during Right V1 Upstates')
title('Up state from Left V1 relative to Up state from Right V1')
ylabel('probability')
xlabel('Time(s)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
%         plot_perievent_event_histogram(spindles(2).SWS_peaktimes,spindles(1).SWS_peaktimes,'twin',[-1 1],'event_name','Right spindles events')
time_wondows = [-2 2];
time_bin = 0.05;
probabilities = calculate_event_probability(spindles(1).SWS_peaktimes,spindles(2).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
title('Probability of left spindles during right spindles')
title('Spindles from Left V1 relative to Spindles from Right V1')
ylabel('probability')
xlabel('Time(s)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
%         plot_perievent_event_histogram(V1_reactivations(2).SWS_onset,V1_reactivations(1).SWS_onset,'twin',[-1 1],'event_name','Right V1 populational burtsting events')
time_wondows = [-2 2];
time_bin = 0.02;
probabilities = calculate_event_probability(V1_reactivations(1).SWS_onset,V1_reactivations(2).SWS_onset, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
title('Probability of left V1 bursting during right V1 bursting')
ylabel('probability')
xlabel('Time(s)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
%         plot_perievent_event_histogram(reactivations(2).SWS_onset,reactivations(1).SWS_onset,'twin',[-1 1],'event_name','Right HPC populational burtsting events')
time_wondows = [-2 2];
time_bin = 0.02;
probabilities = calculate_event_probability(reactivations(2).SWS_onset,reactivations(1).SWS_onset, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,3))
title('Probability of left HPC bursting during right HPC bursting')
ylabel('probability')
xlabel('Time(s)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
sgtitle('Left and Right events Synchronisation')


figure
nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(1).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Left ripples during UP states')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(1).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Left ripples during DOWN states')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(1).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Right ripples during UP states')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(1).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Right ripples relative to DOWN states')
sgtitle('Left V1')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
sgtitle('Left V1 slow waves and HPC interaction')

figure
nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(2).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Left ripples relative to UP states')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, slow_waves(2).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Left ripples relative to DOWN states')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(2).ints.UP(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Right ripples relative to UP states')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, slow_waves(2).ints.DOWN(:,1), [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Right ripples relative to DOWN states')
sgtitle('Right V1')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

sgtitle('Right V1 slow waves and HPC interaction')


figure
nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, spindles(1).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Left ripples during Left spindles')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, spindles(1).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Right ripples during Left spindles')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
probabilities = calculate_event_probability(ripples(1).SWS_peaktimes, spindles(2).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Left ripples during Right spindles')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
time_wondows = [-1 1];
time_bin = 0.02;
probabilities = calculate_event_probability(ripples(2).SWS_peaktimes, spindles(2).SWS_peaktimes, [time_wondows(1):time_bin:time_wondows(2)])
plot([time_wondows(1)+time_bin/2:time_bin:time_wondows(2)-time_bin/2],smooth(probabilities,5))
title('probability of Right ripples during Right spindles')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

sgtitle('V1 spindles and HPC interaction')

% Spike relative to events

spatial_cell_index = find(clusters_combined.odd_even_stability(:,1)>0.95 ...
    | clusters_combined.odd_even_stability(:,2)>0.95);

metric_param =[];
metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
metric_param.region = @(x) contains(x,'V1_L');
[selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
all_spikes{1}=[selected_clusters.spike_id selected_clusters.spike_times];

metric_param.region = @(x) contains(x,'V1_R');
[selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
all_spikes{2}=[selected_clusters.spike_id selected_clusters.spike_times];

metric_param.region = @(x) contains(x,'HPC');
[selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
all_spikes{3}=[selected_clusters.spike_id selected_clusters.spike_times];
group_name = {'V1_L','V1_R','HPC'};

% UP
plot_perievent_spiketime_histogram(all_spikes,slow_waves(1).ints.UP(:,1),'group','by cell zscore','group_name',group_name,'event_name','left UP onset','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,slow_waves(1).ints.UP(:,1),'group','by region','group_name',group_name,'event_name','left UP onset','twin',[-1 1])

plot_perievent_spiketime_histogram(all_spikes,slow_waves(2).ints.UP(:,1),'group','by cell zscore','group_name',group_name,'event_name','right UP onset','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,slow_waves(2).ints.UP(:,1),'group','by region','group_name',group_name,'event_name','right UP onset','twin',[-1 1])

% Ripples
plot_perievent_spiketime_histogram(all_spikes,ripples(1).SWS_peaktimes,'group','by cell zscore','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,ripples(1).SWS_peaktimes,'group','by region','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,ripples(2).SWS_peaktimes,'group','by cell zscore','group_name',group_name,'event_name','Right Ripples','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,ripples(2).SWS_peaktimes,'group','by region','group_name',group_name,'event_name','Right Ripples','twin',[-1 1])


plot_perievent_spiketime_histogram(all_spikes,ripples(1).awake_peaktimes,'group','by cell zscore','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,ripples(1).awake_peaktimes,'group','by region','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])

% Ripples
metric_param =[];
metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
metric_param.region = @(x) contains(x,'V1');
[selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);

ia = find(cluster_id);
C = clusters_combined.cluster_id(ia)
event_id = [ones(1,length(ripples(1).SWS_peaktimes))];
event_times = [ripples(1).SWS_peaktimes];
plot_perievent_spiketimes(clusters_combined.spike_times,clusters_combined.spike_id,[],[],[5 1],[-1 1],0.02,...
    'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times,...
    'event_id',event_id,'event_label','ripple','place_fields',place_fields,'plot_option','by time');

%

for nprobe = 1:2
    [ripples(nprobe).DOWN_UP_transition_offset,ripples(nprobe).DOWN_UP_transition_index] = RestrictInts(ripples(nprobe).offset,[slow_waves(1).ints.DOWN(:,1) slow_waves(1).ints.DOWN(:,2)]);
    ripples(nprobe).DOWN_UP_transition_onset = ripples(nprobe).onset(ripples(nprobe).DOWN_UP_transition_index)';
end

for nprobe = 1:2
    [ripples(nprobe).UP_DOWN_transition_offset,ripples(nprobe).UP_DOWN_transition_index] = RestrictInts(ripples(nprobe).offset,[slow_waves(1).ints.UP(:,1) slow_waves(1).ints.UP(:,1)+0.15]);
    ripples(nprobe).UP_DOWN_transition_onset = ripples(nprobe).onset(ripples(nprobe).UP_DOWN_transition_index)';
end

T1_events = ripples(nprobe).DOWN_UP_transition_onset'; 
T2_events = ripples(nprobe).UP_DOWN_transition_onset';

slow_waves(1).ints.UP(:,1)

find(V1_reactivations(2).ripple_peak(V1_reactivations(2).SWS_index)>=3)

zero_meaned_log_odds = zscore([decoded_ripple_events(1).track(1).replay_events(:).z_log_odds]);

T2_events = ripples(2).peaktimes(find(zero_meaned_log_odds<-1));
[T2_events,~] = RestrictInts(T2_events,behavioural_state(1).SWS);

T1_events = ripples(2).peaktimes(find(zero_meaned_log_odds>1));
[T1_events,~] = RestrictInts(T1_events,behavioural_state(1).SWS);




T2_events = ripples(2).peaktimes(...
    find(reactivation_strength(2).track(1).strength_percentile>0.95 & reactivation_strength(2).track(2).strength_percentile<0.95));
[T2_events,~] = RestrictInts(T2_events,behavioural_state(1).SWS);

T1_events = ripples(2).peaktimes(...
    find(reactivation_strength(2).track(1).strength_percentile<0.95 & reactivation_strength(2).track(2).strength_percentile>0.95));
[T1_events,~] = RestrictInts(T1_events,behavioural_state(1).SWS);

%         T1_events = ripples(1).peaktimes(find(zero_meaned_log_odds>0.5));
plot_perievent_spiketime_histogram(all_spikes,T2_events,'group','by cell zscore','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])
plot_perievent_spiketime_histogram(all_spikes,T2_events,'group','by region','group_name',group_name,'event_name','Left Ripples','twin',[-1 1])


metric_param =[];
metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
metric_param.region = @(x) contains(x,'V1');
[selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
ia = find(cluster_id);
%         ia = ia((41:60));
C = clusters_combined.cluster_id(ia);

event_id = [ones(1,length(T1_events)) 2*ones(1,length(T2_events))];
event_times = [T1_events; T2_events];
[event_times,index] = sort(event_times);
event_id=event_id(index);

plot_perievent_spiketimes(clusters_combined.spike_times,clusters_combined.spike_id,[],[],[5 1],[-1 1],0.02,...
    'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C,'event_times',event_times,...
    'event_id',event_id,'event_label','ripple T1','place_fields',place_fields,'plot_option','by time');


