%%%%%%% Main Time frequency analysis

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]
session_count = 0;

slow_waves_TF_all = struct();
ripples_TF_all = struct();
spindles_TF_all = struct();

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


    if length(stimulus_name)>1 % Based on if POST or PRE
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
        else

            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));

            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
            else
                n =1;
            end
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
    options = session_info(n).probe(1);

    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
    if isempty(DIR)
        continue
    end

    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
    DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

    if ~isempty(DIR)
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        session_clusters_RUN=session_clusters;
        clear session_clusters
    end

    if ~isempty(DIR1)
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        session_clusters_RUN=session_clusters;
        clear session_clusters
    end

    if contains(stimulus_name{n},'Masa2tracks')
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        % %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         clusters=clusters_ks4;
    elseif contains(stimulus_name{n},'Sleep')
        %         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_amplitude.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_phase.mat'));

        
    else
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
        %         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    end

    
    timevec = [TF_phase_HPC(:).timebin];
    freqs = [TF_amp_V1(:).freq];

    event_types = {'ripples','spindles','UP_ints','DOWN_ints'};

    baseline_win = [-2 -1.5];  % baseline time window
    baseline_idx = timevec >= baseline_win(1) & timevec <= baseline_win(2);


    UP_index1 = slow_waves(1).UP_ints(:,2) - slow_waves(1).UP_ints(:,1) <= 2;
    UP_index2 =slow_waves(2).UP_ints(:,2) - slow_waves(2).UP_ints(:,1) <= 2;

    [C,DOWN_index1,ib] = intersect(slow_waves(1).DOWN_ints(:,1), slow_waves(1).UP_ints(UP_index1,2));
    [C,DOWN_index2,ib] = intersect(slow_waves(2).DOWN_ints(:,1), slow_waves(2).UP_ints(UP_index1,2));


    for e = 1:length(event_types)
        if contains(event_types{e},'ripples')
            index1 = ripples(1).SWS_index;
            index2 = ripples(2).SWS_index;
        elseif contains(event_types{e},'spindles')
            index1 = ripples(1).SWS_index;
            index2 = ripples(2).SWS_index;

        elseif contains(event_types{e},'UP')
            index1 = UP_index1;
            index2 = UP_index2;
        elseif contains(event_types{e},'DOWN')
            index1 = DOWN_index1;
            index2 = DOWN_index2;
        end

        TF_amp_V1_ipsi = squeeze([cat(4,TF_amp_V1(1).(event_types{e})(1,:,:,index1), TF_amp_V1(2).(event_types{e})(2,:,:,index2))]);
        TF_amp_V1_ipsi = normalize_tf(TF_amp_V1_ipsi, baseline_idx, 'dB');
        mean_TF_amp_V1_ipsi = mean(TF_amp_V1_ipsi,3,'omitnan');
        [~, sig_mask_V1_ipsi] = permutation_TF_test(TF_amp_V1_ipsi, 1000, 0.05);


        TF_amp_V1_contra = squeeze([cat(4,TF_amp_V1(1).(event_types{e})(2,:,:,index1), TF_amp_V1(2).(event_types{e})(1,:,:,index2))]);
        TF_amp_V1_contra = normalize_tf(TF_amp_V1_contra, baseline_idx, 'dB');
        mean_TF_amp_V1_contra = mean(TF_amp_V1_contra,3,'omitnan');
        [~, sig_mask_V1_contra] = permutation_TF_test(TF_amp_V1_contra, 1000, 0.05);


        TF_amp_HPC_ipsi = squeeze([cat(4,TF_amp_HPC(1).(event_types{e})(1,:,:,index1), TF_amp_HPC(2).(event_types{e})(2,:,:,index2))]);
        TF_amp_HPC_ipsi = normalize_tf(TF_amp_HPC_ipsi, baseline_idx, 'dB');
        mean_TF_amp_HPC_ipsi = mean(TF_amp_HPC_ipsi,3,'omitnan');
        [~, sig_mask_HPC_ipsi] = permutation_TF_test(TF_amp_HPC_ipsi, 1000, 0.05);


        TF_amp_HPC_contra = squeeze([cat(4,TF_amp_HPC(1).(event_types{e})(2,:,:,index1), TF_amp_HPC(2).(event_types{e})(1,:,:,index2))]);
        TF_amp_HPC_contra = normalize_tf(TF_amp_HPC_contra, baseline_idx, 'dB');
        mean_TF_amp_HPC_contra = mean(TF_amp_HPC_contra,3,'omitnan');
        [~, sig_mask_HPC_contra] = permutation_TF_test(TF_amp_HPC_contra, 1000, 0.05);




        % TF_amp_V1_ipsi = mean(squeeze([cat(4,TF_amp_V1(1).(event_types{e})(1,:,:,:), TF_amp_V1(2).(event_types{e})(2,:,:,:))]),3,'omitnan');
        % TF_amp_V1_contra = mean(squeeze([cat(4,TF_amp_V1(1).(event_types{e})(2,:,:,:), TF_amp_V1(2).(event_types{e})(1,:,:,:))]),3,'omitnan');
        % TF_amp_HPC_ipsi = mean(squeeze([cat(4,TF_amp_HPC(1).(event_types{e})(1,:,:,:), TF_amp_HPC(2).(event_types{e})(2,:,:,:))]),3,'omitnan');
        % TF_amp_HPC_contra = mean(squeeze([cat(4,TF_amp_HPC(1).(event_types{e})(2,:,:,:), TF_amp_HPC(2).(event_types{e})(1,:,:,:))]),3,'omitnan');

        % === V1 Phase ===
        TF_phase_V1_ipsi = squeeze(cat(4, ...
            TF_phase_V1(1).(event_types{e})(1,:,:,index1), ... % V1 L to L events
            TF_phase_V1(2).(event_types{e})(2,:,:,index2)  ... % V1 R to R events
            ));  % size: [freq x time x events]

        TF_phase_V1_contra = squeeze(cat(4, ...
            TF_phase_V1(1).(event_types{e})(2,:,:,index1), ... % V1 L to R events
            TF_phase_V1(2).(event_types{e})(1,:,:,index2)  ... % V1 R to L events
            ));  % size: [freq x time x events]

        % === HPC Phase ===
        TF_phase_HPC_ipsi = squeeze(cat(4, ...
            TF_phase_HPC(1).(event_types{e})(1,:,:,index1), ...
            TF_phase_HPC(2).(event_types{e})(2,:,:,index2)  ...
            ));

        TF_phase_HPC_contra = squeeze(cat(4, ...
            TF_phase_HPC(1).(event_types{e})(2,:,:,index1), ...
            TF_phase_HPC(2).(event_types{e})(1,:,:,index2)  ...
            ));

        % % Compute average complex phase vector
        % avg_phase_V1_ipsi = angle(mean(exp(1i * TF_phase_V1_ipsi), 3, 'omitnan'));     % [freq x time]
        % avg_phase_V1_contra = angle(mean(exp(1i * TF_phase_V1_contra), 3, 'omitnan'));
        % 
        % avg_phase_HPC_ipsi = angle(mean(exp(1i * TF_phase_HPC_ipsi), 3, 'omitnan'));
        % avg_phase_HPC_contra = angle(mean(exp(1i * TF_phase_HPC_contra), 3, 'omitnan'));



        % timevec = linspace(-2, 2, size(TF_amp_V1_ipsi, 2));  % time axis
        % freq_band = [0.5 300];
        % freqs = logspace(log10(freq_band(1)), log10(freq_band(2)), 100);

        figure('Name',['TF amplitude - ' event_types{e}],'Position',[100 100 1200 500])

        subplot(2,2,1)
        % contourf(timevec, log10(freqs), TF_amp_V1_ipsi,40,'linecolor','none')
        contourf(timevec, log2(freqs), mean_TF_amp_V1_ipsi,40,'linecolor','none')
        axis xy; title(['Ipsi V1 - ' event_types{e}])
        xlabel('Time (s)'), ylabel('log_{2}(Freq Hz)')
        % hold on; contour(timevec, log10(freqs), sig_mask_ipsi, [1 1], 'r--')
        % set(gca,'ylim',[0 300]);
        set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
            'YTickLabel', {'1','2','4','8','16','32','64','128','256'});
        hold on;
        
        colorbar
        xlim([-2 2])
        ylim([min(log2(freqs)) max(log2(freqs))])
        % clim([-1 3])
        yline(log2([9 15]))
        yline(log2([1 4]))
        yline(log2([125 300]))

        subplot(2,2,2)
        contourf(timevec, log2(freqs), mean_TF_amp_V1_contra,40,'linecolor','none')
        axis xy; title(['Contra V1 - ' event_types{e}])
        xlabel('Time (s)'), ylabel('log_{2}(Freq Hz)')
        % hold on; contour(timevec, log10(freqs), sig_mask_ipsi, [1 1], 'r--')
        % set(gca,'ylim',[0 300]);
        set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
            'YTickLabel', {'1','2','4','8','16','32','64','128','256'});
        colorbar
        xlim([-2 2])
        ylim([min(log2(freqs)) max(log2(freqs))])
        % clim([-1 3])
        yline(log2([9 15]))
        yline(log2([1 4]))
        yline(log2([125 300]))


        subplot(2,2,3)
        contourf(timevec, log2(freqs), mean_TF_amp_HPC_ipsi,40,'linecolor','none')
        axis xy; title(['Ipsi HPC - ' event_types{e}])
        xlabel('Time (s)'), ylabel('log_{2}(Freq Hz)')
        % hold on; contour(timevec, log10(freqs), sig_mask_ipsi, [1 1], 'r--')
        % set(gca,'ylim',[0 300]);
        set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
            'YTickLabel', {'1','2','4','8','16','32','64','128','256'});
        colorbar
        xlim([-2 2])
        ylim([min(log2(freqs)) max(log2(freqs))])
        % clim([-1 3])
        yline(log2([9 15]))
        yline(log2([1 4]))
        yline(log2([125 300]))

        subplot(2,2,4)
        contourf(timevec, log2(freqs), mean_TF_amp_HPC_contra,40,'linecolor','none')
        axis xy; title(['Ipsi HPC - ' event_types{e}])
        xlabel('Time (s)'), ylabel('log_{2}(Freq Hz)')
        % hold on; contour(timevec, log10(freqs), sig_mask_ipsi, [1 1], 'r--')
        % set(gca,'ylim',[0 300]);
        set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
            'YTickLabel', {'1','2','4','8','16','32','64','128','256'});
        colorbar
         xlim([-2 2])
        ylim([min(log2(freqs)) max(log2(freqs))])
        % clim([-1 3])
        yline(log2([9 15]))
        yline(log2([1 4]))
        yline(log2([125 300]))

        % 
        % figure('Name',['TF phase - ' event_types{e}],'Position',[100 100 1200 500])
        % 
        % subplot(2,2,1)
        % imagesc(timevec, log10(freqs), avg_phase_V1_ipsi)
        % axis xy; title(['IPSI V1 - ' event_types{e}])
        % xlabel('Time (s)'), ylabel('log_{10}(Freq Hz)')
        % % hold on; contour(timevec, log10(freqs), sig_mask_ipsi, [1 1], 'r--')
        % colorbar
        % 
        % subplot(2,2,2)
        % imagesc(timevec, log10(freqs), avg_phase_V1_contra)
        % axis xy; title(['CONTRA V1 - ' event_types{e}])
        % xlabel('Time (s)'), ylabel('log_{10}(Freq Hz)')
        % % hold on; contour(timevec, log10(freqs), sig_mask_contra, [1 1], 'r--')
        % colorbar
        % 
        % 
        % subplot(2,2,3)
        % imagesc(timevec, log10(freqs), avg_phase_HPC_ipsi)
        % axis xy; title(['IPSI HPC - ' event_types{e}])
        % xlabel('Time (s)'), ylabel('log_{10}(Freq Hz)')
        % % hold on; contour(timevec, log10(freqs), sig_mask_ipsi, [1 1], 'r--')
        % colorbar
        % 
        % subplot(2,2,4)
        % imagesc(timevec, log10(freqs), avg_phase_HPC_contra)
        % axis xy; title(['CONTRA HPC - ' event_types{e}])
        % xlabel('Time (s)'), ylabel('log_{10}(Freq Hz)')
        % % hold on; contour(timevec, log10(freqs), sig_mask_contra, [1 1], 'r--')
        % colorbar

    end


    TF_amp_V1
    TF_amp_HPC
    TF_phase_V1
    TF_phase_HPC
end

