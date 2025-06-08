%% Main pipeline for viewing sleep events

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
Stimulus_type = 'Sleep';
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]

for nsession =6
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
        cortex_LFP_SO_filtered=[];
        CA1_LFP_ripple_filtered=[];
        CA1_LFP_filtered=[];
        LFP_SR = 1/mean(diff(LFP(nprobe).tvec));

        for nprobe = 1:length(session_clusters_RUN.probe_hemisphere)

            probe_no=session_clusters_RUN.probe_hemisphere(nprobe);



            cortex_LFP{probe_no} = LFP(probe_no).best_V1(LFP(probe_no).best_V1_channel ==  slow_waves(probe_no).best_channel,:);

            passband = [0.5 17];
            filter_type  = 'bandpass';
            filter_order = round(6*LFP_SR/(max(passband)-min(passband)));  % creates filter for ripple
            norm_freq_range = passband/(LFP_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
            b_ripple = fir1(filter_order, norm_freq_range,filter_type);
            cortex_LFP_filtered{probe_no} = filtfilt(b_ripple,1,cortex_LFP{probe_no});
            %         cortex_LFP{probe_no} = zscore(abs(hilbert(signal)));

            passband = [0.5 4];
            filter_type  = 'bandpass';
            filter_order = round(6*LFP_SR/(max(passband)-min(passband)));  % creates filter for ripple
            norm_freq_range = passband/(LFP_SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
            b_ripple = fir1(filter_order, norm_freq_range,filter_type);
            cortex_LFP_SO_filtered{probe_no} = filtfilt(b_ripple,1,cortex_LFP{probe_no});
            %         cortex_LFP{probe_no} = zscore(abs(hilbert(signal)));



            CA1_LFP{probe_no} = LFP(probe_no).best_HPC(ripples(probe_no).best_channel,:);

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



        %
        %         %%%%%  Sleep detection
        %
        %         speedTreshold = 1;
        %         speed = double(session_clusters.mobility_thresholded{1});
        %         speed = interp1(Behaviour.tvec,speed,tvec,'linear');
        %
        %         [freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
        %             [tvec' LFP(probe_no).best_V1(LFP(probe_no).best_V1_channel ==  slow_waves(probe_no).best_channel,:)'],[tvec' LFP(probe_no).best_HPC(ripples(probe_no).best_channel,:)'],...
        %             [tvec' speed'],speedTreshold);


        %%% Get spike counts

        tvec = LFP(1).tvec;
        tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];
        CA1_spike_counts=[];
        V1_spike_counts=[];
        w = gausswin(0.05*1/mean(diff(tvec)));
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

        [s,f,t] = stft(cortex_LFP{2},LFP_SR,'Window',hann(LFP_SR),'FrequencyRange','onesided','FFTLength',10*round(LFP_SR));
        %         index = find(f<20);
        index = find(f>1 & f<30);
        %         sss= abs(s).^2;
%         S_mag_smoothed_V1 = imgaussfilt(abs(s(index,:)), 20); % Gaussian smoothing with sigma=2
        S_mag_smoothed_V1 = abs(s(index,:));
        S_mag_smoothed_V1 = movmedian(S_mag_smoothed_V1,20,2);
        S_mag_smoothed_V1 = movmedian(S_mag_smoothed_V1,20,1);

        [s,f,t] = stft(CA1_LFP{2},LFP_SR,'Window',hann(LFP_SR),'FrequencyRange','onesided','FFTLength',10*round(LFP_SR));
        index = find(f>1 & f<30);
        %         sss= abs(s).^2;
%         S_mag_smoothed_HPC = imgaussfilt(abs(s(index,:)),20); % Gaussian smoothing with sigma=2
        S_mag_smoothed_HPC = abs(s(index,:));
        S_mag_smoothed_HPC = movmedian(S_mag_smoothed_HPC,20,2);
        S_mag_smoothed_HPC = movmedian(S_mag_smoothed_HPC,20,1);
% 
%         % Interpolate to log-frequency scale
%         f_log = logspace(log10(min(f(index))), log10(max(f(index))), 200);  % log scale
%         S_log_V1 = interp1(f(index), S_mag_smoothed_V1, f_log, 'linear', 'extrap');
% 
%         S_log_HPC = interp1(f(index), S_mag_smoothed_HPC, f_log, 'linear', 'extrap');
% 



        freq_orig = f(index);         % size: [N_freq]
        time_orig = t;                % size: [1 x N_time]
        S_orig = S_mag_smoothed_V1;   % size: [N_freq x N_time]

        % Define new target grids
        %         f_log = logspace(log10(min(freq_orig)), log10(max(freq_orig)), 200);  % log-freq axis
        f_log = 2 .^ linspace(log2(min(freq_orig)), log2(max(freq_orig)), 200);
        S_log_V1 = interp1(freq_orig, S_mag_smoothed_V1, f_log, 'linear', 'extrap');
        S_log_HPC = interp1(freq_orig, S_mag_smoothed_HPC, f_log, 'linear', 'extrap');



        tindex = find(t>behavioural_state_merged.REM(2,1)-290 & t<behavioural_state_merged.REM(2,1)+300);
        tindex = find(t>3800 & t<5500)
% tindex = 1:length(t);

        theta_delta_ratio = mean(S_mag_smoothed_HPC(f>4 & f<12,:),1)./mean(S_mag_smoothed_HPC(f>0.5 & f<4,:),1);



        fig = figure;
        fig.Name = 'NREM and REM sleep scoring based on V1 and HPC LFP';
        fig.Position = [880 110 950 840];
        subplot(4,1,1)
        imagesc(t(tindex), freq_orig, S_mag_smoothed_V1(:,tindex));
        set(gca, 'YDir', 'normal');               % Ensure Y-axis increases upward
        %         set(gca,YScale="log",...
        %             YDir="reverse",View=[0 300])
        %         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed_V1,1,[]),1) prctile(reshape(S_mag_smoothed_V1,1,[]),90)])
        colorbar
        colormap(flipud(bone))
        hold on

        yl = ylim;
        hold on
        for nevent = 1:size(behavioural_state_merged.SWS,1)
            x_start = behavioural_state_merged.SWS(nevent,1);
            x_end = behavioural_state_merged.SWS(nevent,2);
            fill([x_start x_end x_end x_start], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

%         yline(f(min(find(f>0.5&f<4))),'r')
%         yline(f(max(find(f>0.5&f<4))),'r')
        ylabel('Frequency (Hz)')
        title('V1 LFP Time frequency plot')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(4,1,2)
        imagesc(t(tindex), freq_orig, S_mag_smoothed_HPC(:,tindex));
        set(gca, 'YDir', 'normal');               % Ensure Y-axis increases upward
        %         set(gca,YScale="log",...
        %             YDir="reverse",View=[0 300])
        %         clim([min(reshape(abs(s).^2,1,[])) max(reshape(abs(s).^2,1,[]))])
        clim([prctile(reshape(S_mag_smoothed_HPC,1,[]),1) prctile(reshape(S_mag_smoothed_HPC,1,[]),90)])
        colorbar
        colormap(flipud(bone))
        hold on

%         yline(f(min(find(f>4&f<12))),'b')
%         yline(f(max(find(f>4&f<12))),'b')

        yl = ylim;
        for nevent = 1:size(behavioural_state_merged.REM,1)
            x_start = behavioural_state_merged.REM(nevent,1);
            x_end = behavioural_state_merged.REM(nevent,2);
            fill([x_start x_end x_end x_start], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        ylabel('Frequency (Hz)')
        title('HPC LFP Time frequency plot')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        subplot(4,1,3)
        plot(t(tindex),theta_delta_ratio(tindex),'k')
        yl = ylim;
        hold on
        for nevent = 1:size(behavioural_state_merged.REM,1)
            x_start = behavioural_state_merged.REM(nevent,1);
            x_end = behavioural_state_merged.REM(nevent,2);
            fill([x_start x_end x_end x_start], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        xlim([t(tindex(1)) t(tindex(end))])
        colorbar
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        ylabel('Theta/Delta ratio')

        subplot(4,1,4)
        interp1_mobility = interp1(Behaviour.tvec,Behaviour.mobility_zscore,t(tindex));
        plot(t(tindex),interp1_mobility,'k')
        xlabel('Time(s)')
        ylabel('Mobility zscored')
        yl = ylim;
        hold on
        for nevent = 1:size(behavioural_state_merged.SWS,1)
            x_start = behavioural_state_merged.SWS(nevent,1);
            x_end = behavioural_state_merged.SWS(nevent,2);
            fill([x_start x_end x_end x_start], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        hold on
        for nevent = 1:size(behavioural_state_merged.REM,1)
            x_start = behavioural_state_merged.REM(nevent,1);
            x_end = behavioural_state_merged.REM(nevent,2);
            fill([x_start x_end x_end x_start], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        xlim([min(t(tindex)) max(t(tindex))])
        colorbar
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','figure'),[],'ContentType','vector')


      

        %%%%%% peri ripple LFP and spiking
        for nprobe = 1:2
            [ripples(nprobe).DOWN_UP_transition_offset,ripples(nprobe).DOWN_UP_transition_index] = RestrictInts(ripples(nprobe).offset,[slow_waves(1).UP_ints(:,1) slow_waves(1).UP_ints(:,1)+0.3]);
            ripples(nprobe).DOWN_UP_transition_onset = ripples(nprobe).onset(ripples(nprobe).DOWN_UP_transition_index)';
        end

        for nprobe = 1:2
            [ripples(nprobe).UP_DOWN_transition_offset,ripples(nprobe).UP_DOWN_transition_index] = RestrictInts(ripples(nprobe).offset,[slow_waves(1).UP_ints(:,2)-0.4 slow_waves(1).UP_ints(:,2)+0.3]);
            ripples(nprobe).UP_DOWN_transition_onset = ripples(nprobe).onset(ripples(nprobe).UP_DOWN_transition_index)';
        end

        % slow_waves = [];

        %         event_times = [spindles(1).onset(spindles(1).peak_zscore>4) spindles(1).offset(spindles(1).peak_zscore>4)];
        %         %         tindex = find(tvec>ripples(1).UP_DOWN_transition_onset(min(nevent))-0.5 & tvec<ripples(1).UP_DOWN_transition_onset(max(nevent))+0.5);
        %         tindex = find(tvec>event_times(min(nevent),1)-0.1 & tvec<event_times(max(nevent),2)+0.1);
        %
        %         count = 1;
        %         nevent =1:48;
        %         for count = 1:48
        %             subplot(7,7,count)
        %             tindex = find(tvec>event_times(nevent(count),1)-0.2 & tvec<event_times(nevent(count),2)+0.2);
        %             plot(tvec(tindex),cortex_LFP_filtered{1}(tindex),'k','LineWidth',1.5)
        %             ylim([min(cortex_LFP_filtered{1}(tindex))-50 max(cortex_LFP_filtered{1}(tindex))+50])
        %         end
        %
        %
        %
        %
        %         nevent =200
        %         event_times = [ripples(1).UP_DOWN_transition_onset' ripples(1).UP_DOWN_transition_offset];
        %         %         tindex = find(tvec>ripples(1).UP_DOWN_transition_onset(min(nevent))-0.5 & tvec<ripples(1).UP_DOWN_transition_onset(max(nevent))+0.5);
        %         tindex = find(tvec>event_times(min(nevent),1)-0.1 & tvec<event_times(max(nevent),2)+0.1);
        %
        %
        %         count = 1;
        %         nevent =200:225;
        %         for count = 1:25
        %             subplot(5,5,count)
        %             tindex = find(tvec>event_times(nevent(count),1)-0.05 & tvec<event_times(nevent(count),2)+0.05);
        %             plot(tvec(tindex),LFP(1).best_HPC(1,tindex),'k','LineWidth',1.5)
        %             ylim([min(LFP(1).best_HPC(1,tindex))-50 max(LFP(1).best_HPC(1,tindex))+50])
        %         end
%         event_times = [ripples(1).onset(ripples(1).UP_DOWN_transition_index ==1 & ripples(1).peak_zscore>10) ripples(1).offset(ripples(1).UP_DOWN_transition_index ==1 & ripples(1).peak_zscore>10)];

        % Get relevant ripple onset times
        down_up_onsets = ripples(1).onset(ripples(1).DOWN_UP_transition_index == 1 & ripples(1).peak_zscore > 5);
        up_down_onsets = ripples(1).onset(ripples(1).UP_DOWN_transition_index == 1 & ripples(1).peak_zscore > 5);

        threshold = 0.5;  % 300 ms in seconds
        event_times = [];

        for i = 1:length(down_up_onsets)
            % Look for UP_DOWN events that come after this DOWN_UP within 300ms
            dt = up_down_onsets - down_up_onsets(i);
            idx = find(dt > 0 & dt <= threshold);

            for j = 1:length(idx)
                event_times = [event_times; down_up_onsets(i), up_down_onsets(idx(j))];
            end
        end


        
        ripple_times = [ripples(1).SWS_onset ripples(1).SWS_offset];
        ripple_times = [ripples(1).SWS_onset ripples(1).SWS_offset];
        UP_times = slow_waves(1).UP_ints;
        DOWN_times = slow_waves(1).DOWN_ints;
        DOWN_times2 = slow_waves(2).DOWN_ints;

        nevent = 3;
        %         nevent = 101:104;
        colour_lines = [0,90,50;74,20,134]/256;

        spike_data = {V1_spikes{1},V1_spikes{2},HPC_spikes{1},HPC_spikes{2}};
        region_names = {'Left V1','Right V1','Left HPC','Right HPC'};


        for nfig = 5
            tindex = find(tvec>event_times(min(nfig),1)-0.8 & tvec<event_times(max(nfig),2)+0.8);
            fig = figure;
            fig.Position = [740 130 440 850]
            fig.Name = sprintf('raster raw race during UP DOWN ripples %i',nfig)





            %         colour_lines = [0,90,50;74,20,134]/256;
            colour_lines = [...
                0.5000, 0.6768, 0.5984;  % lighter green
                0.6445, 0.5391, 0.7617   % lighter purple
                ];
            colour_lines = [...
                0.7000, 0.7886, 0.7492;  % even lighter green
                0.7571, 0.7373, 0.8809   % even lighter purple
                ];

            for nregion = 1:4
                % Map each neuron ID to a row index
                start_time = tvec(tindex(1));
                end_time = tvec(tindex(end));
%                 if nregion <3
%                     bin_size = 0.01; % 1 ms resolution
%                 else
                    bin_size = 0.005; % 1 ms resolution
%                 end
                time_bins =  tvec(tindex(1)):bin_size: tvec(tindex(end));
                num_bins = length(time_bins);
                

                subplot(4,1,nregion)
                neuron_ids = unique(spike_data{nregion}(:,1));
                num_neurons = length(neuron_ids);
                neuron_map = containers.Map(neuron_ids, 1:num_neurons);

                % Initialize raster matrix
                raster_matrix = zeros(num_neurons, num_bins);

                % Fill raster matrix
                for i = 1:size(spike_data{nregion}, 1)
                    neuron = spike_data{nregion}(i,1);
                    time = spike_data{nregion}(i,2);

                    if time >= start_time && time <= end_time
                        row = neuron_map(neuron);
                        col = floor((time - start_time)/bin_size) + 1;
                        raster_matrix(row, col) = 1;
                    end
                end
                % Remove neurons (rows) that did not fire at all
                firing_neurons_mask = any(raster_matrix, 2);  % True for firing neurons
                raster_matrix = raster_matrix(firing_neurons_mask, :);
                num_active = size(raster_matrix, 1);

                % Pad with blank rows at random positions to reach 60 rows
                desired_rows = 65;
                num_pad = desired_rows - num_active;
                pad_matrix = zeros(num_pad, num_bins);

                % Randomly choose positions to insert blank rows
                insert_idx = sort(randperm(desired_rows, num_pad));
                new_raster = zeros(desired_rows, num_bins);

                active_idx = setdiff(1:desired_rows, insert_idx);

                % Fill new raster
                new_raster(active_idx, :) = raster_matrix;
                % Blank rows already zero
                raster_matrix = new_raster;

                imagesc(time_bins, 1:size(raster_matrix,1), raster_matrix);
                % Custom RGB color for spikes (e.g., red)
                bg_rgb = [1 1 1];     % white background
                if nregion == 1|nregion == 3
                    spike_rgb = colour_lines(1,:)  % green
                else
                    spike_rgb = colour_lines(2,:) % purple
                end
                colormap(gca, [bg_rgb; spike_rgb]);

                colorbar off;
                hold on
                yl = ylim;
                for iEvent = 1:size(DOWN_times,1)
                    x_start = DOWN_times(iEvent,1);
                    x_end = DOWN_times(iEvent,2);
                    x_rect = [x_start x_end x_end x_start];
                    y_rect = [yl(1) yl(1) yl(2) yl(2)];
                    patch('XData', x_rect, ...
                        'YData', y_rect, ...
                        'FaceColor','b', ...
                        'FaceAlpha', 0.2, ...
                        'EdgeColor', 'none');
                end
            
                hold on
                for iEvent = 1:size(ripple_times,1)
                    x_start = ripple_times(iEvent,1);
                    x_end = ripple_times(iEvent,2);
                    fill([x_start x_end x_end x_start], ...
                        [yl(1) yl(1) yl(2) yl(2)], ...
                        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end


                xlabel('Time (ms)');
                ylabel('Neuron');
%                 yticks(1:num_neurons);
%                 yticklabels(neuron_ids);
                title(region_names{nregion})


                xlim([tvec(tindex(1)) tvec(tindex(end))])
                set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
            end


            colour_lines = [0,90,50;74,20,134]/256;

            fig = figure;
            fig.Position = [740 130 440 850]
            fig.Name = sprintf('LFP and MUA raw race during UP DOWN ripples %i',nfig)

            subplot(5,1,[1 2 3])
            plot(tvec(tindex),4*cortex_LFP_SO_filtered{1}(tindex)+4000,'Color',colour_lines(1,:))
            hold on
            plot(tvec(tindex),6*cortex_LFP_SO_filtered{2}(tindex)+3000,'Color',colour_lines(2,:))
            hold on;
            plot(tvec(tindex),CA1_LFP_filtered{1}(tindex)+2000,'Color',colour_lines(1,:))
            % plot(tvec(tindex),CA1_LFP_filtered{1}(tindex)+500,'k')
            plot(tvec(tindex),CA1_LFP_ripple_filtered{1}(tindex)*2+1500,'Color',colour_lines(1,:))

            plot(tvec(tindex),2*CA1_LFP_filtered{2}(tindex)+500,'Color',colour_lines(2,:))
            % plot(tvec(tindex),CA1_LFP_filtered{1}(tindex)+500,'k')
            plot(tvec(tindex),CA1_LFP_ripple_filtered{2}(tindex)*3,'Color',colour_lines(2,:))


            yl = ylim;
            this_colour = colour_lines(1,:);
            for iEvent = 1:size(DOWN_times,1)
                x_start = DOWN_times(iEvent,1);
                x_end = DOWN_times(iEvent,2);
                x_rect = [x_start x_end x_end x_start];
                y_rect = [yl(1) yl(1) yl(2) yl(2)];
                patch('XData', x_rect, ...
                    'YData', y_rect, ...
                    'FaceColor','b', ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
            end

%             for iEvent = 1:size(UP_times,1)
%                 x_start = UP_times(iEvent,1);
%                 x_end = UP_times(iEvent,2);
%                 x_rect = [x_start x_end x_end x_start];
%                 y_rect = [yl(1) yl(1) yl(2) yl(2)];
%                 patch('XData', x_rect, ...
%                     'YData', y_rect, ...
%                     'FaceColor','r', ...
%                     'FaceAlpha', 0.2, ...
%                     'EdgeColor', 'none');
%             end

            yl = ylim;
            this_colour = colour_lines(1,:);
            for iEvent = 1:size(DOWN_times2,1)
                x_start = DOWN_times2(iEvent,1);
                x_end = DOWN_times2(iEvent,2);
                x_rect = [x_start x_end x_end x_start];
                y_rect = [yl(1) yl(1) yl(2) yl(2)];
                patch('XData', x_rect, ...
                    'YData', y_rect, ...
                    'FaceColor', colour_lines(2,:), ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
            end

            for iEvent = 1:size(ripple_times,1)
                x_start = ripple_times(iEvent,1);
                x_end = ripple_times(iEvent,2);
                fill([x_start x_end x_end x_start], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
            %             for iEvent = 1:size(UP_times,1)
            %                 x_start = UP_times(iEvent,1);
            %                 x_end = UP_times(iEvent,2);
            %                 fill([x_start x_end x_end x_start], ...
            %                     [yl(1) yl(1) yl(2) yl(2)], ...
            %                     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            %             end

            xlim([tvec(tindex(1)) tvec(tindex(end))])
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            subplot(5,1,4)
            %         CA1_spike_counts{3} = zscore(CA1_spike_counts{1}+CA1_spike_counts{2});
            %         V1_spike_counts{3} = zscore(V1_spike_counts{1}+V1_spike_counts{2});

            plot(tvec(tindex),2*V1_spike_counts{1}(tindex)+2,'Color',colour_lines(1,:))
            hold on
            plot(tvec(tindex),2*V1_spike_counts{2}(tindex),'Color',colour_lines(2,:))
            title('V1 spiking')
            yl = ylim;

            for iEvent = 1:size(DOWN_times,1)
                x_start = DOWN_times(iEvent,1);
                x_end = DOWN_times(iEvent,2);
                x_rect = [x_start x_end x_end x_start];
                y_rect = [yl(1) yl(1) yl(2) yl(2)];
                patch('XData', x_rect, ...
                    'YData', y_rect, ...
                      'FaceColor','b', ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
            end

            yl = ylim;
            this_colour = colour_lines(1,:);
            for iEvent = 1:size(DOWN_times2,1)
                x_start = DOWN_times2(iEvent,1);
                x_end = DOWN_times2(iEvent,2);
                x_rect = [x_start x_end x_end x_start];
                y_rect = [yl(1) yl(1) yl(2) yl(2)];
                patch('XData', x_rect, ...
                    'YData', y_rect, ...
                    'FaceColor', colour_lines(2,:), ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
            end

            for iEvent = 1:size(ripple_times,1)
                x_start = ripple_times(iEvent,1);
                x_end = ripple_times(iEvent,2);
                fill([x_start x_end x_end x_start], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end


            xlim([tvec(tindex(1)) tvec(tindex(end))])
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

            subplot(5,1,5)
            plot(tvec(tindex),CA1_spike_counts{1}(tindex)+3,'Color',colour_lines(1,:))
            hold on
            plot(tvec(tindex),CA1_spike_counts{2}(tindex),'Color',colour_lines(2,:))
            title('HPC spiking')
            yl = ylim;


            for iEvent = 1:size(DOWN_times,1)
                x_start = DOWN_times(iEvent,1);
                x_end = DOWN_times(iEvent,2);
                x_rect = [x_start x_end x_end x_start];
                y_rect = [yl(1) yl(1) yl(2) yl(2)];
                patch('XData', x_rect, ...
                    'YData', y_rect, ...
                      'FaceColor','b', ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
            end

            yl = ylim;
            this_colour = colour_lines(1,:);
            for iEvent = 1:size(DOWN_times2,1)
                x_start = DOWN_times2(iEvent,1);
                x_end = DOWN_times2(iEvent,2);
                x_rect = [x_start x_end x_end x_start];
                y_rect = [yl(1) yl(1) yl(2) yl(2)];
                patch('XData', x_rect, ...
                    'YData', y_rect, ...
                    'FaceColor', colour_lines(2,:), ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
            end

            for iEvent = 1:size(ripple_times,1)
                x_start = ripple_times(iEvent,1);
                x_end = ripple_times(iEvent,2);
                fill([x_start x_end x_end x_start], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
            xlim([tvec(tindex(1)) tvec(tindex(end))])
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        end

        save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','Raw trace M24016'),[],'ContentType','vector')

    end

end

