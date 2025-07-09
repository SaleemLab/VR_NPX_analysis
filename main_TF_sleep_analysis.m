%%%%%%% Main Time frequency analysis (Average activities + extract -200-0 and 0-200ms phase coupling and average SO/Spindle/ripples power)

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

ripples_TF_stats = struct();
UP_TF_stats = struct();
DOWN_TF_stats = struct();
spindles_TF_stats = struct();



for nsession =1:length(experiment_info)
    tic
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

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_amplitude_ft.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_TF_phase_ft.mat'));

        
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

    event_types = {'ripples','spindles','UP','DOWN'};

    baseline_win = [-2 -1.5];  % baseline time window
    baseline_idx = timevec >= baseline_win(1) & timevec <= baseline_win(2);


    UP_index1 = slow_waves(1).UP_ints(:,2) - slow_waves(1).UP_ints(:,1) <= 2;
    UP_index2 =slow_waves(2).UP_ints(:,2) - slow_waves(2).UP_ints(:,1) <= 2;

    [C,DOWN_index1,ib] = intersect(slow_waves(1).DOWN_ints(:,1), slow_waves(1).UP_ints(UP_index1,2));
    [C,DOWN_index2,ib] = intersect(slow_waves(2).DOWN_ints(:,1), slow_waves(2).UP_ints(UP_index2,2));
   
    
    for e = 1:length(event_types)
        if contains(event_types{e},'ripples')
            index1 = ripples(1).SWS_index;
            index2 = ripples(2).SWS_index;

            [status_DOWN,interval,index] = InIntervals([ripples(1).SWS_onset; ripples(2).SWS_onset],[[slow_waves(1).DOWN_ints(DOWN_index1,1)-0.2 slow_waves(1).DOWN_ints(DOWN_index1,1)];...
                [slow_waves(2).DOWN_ints(DOWN_index2,1)-0.2 slow_waves(2).DOWN_ints(DOWN_index2,1)]]);

            [status_UP,interval,index] = InIntervals([ripples(1).SWS_onset; ripples(2).SWS_onset],[[slow_waves(1).UP_ints(UP_index1,1) slow_waves(1).UP_ints(UP_index1,1)+0.2];...
                [slow_waves(2).UP_ints(UP_index2,1) slow_waves(2).UP_ints(UP_index2,1)+0.2]]);

        elseif contains(event_types{e},'spindles')
            index1 = spindles(1).SWS_index;
            index2 = spindles(2).SWS_index;
            [status_DOWN,interval,index] = InIntervals([slow_waves(1).DOWN_ints(DOWN_index1,1);slow_waves(2).DOWN_ints(DOWN_index2,1)],[ripples(1).onset ripples(1).onset+0.2]);

        elseif contains(event_types{e},'UP')
            index1 = UP_index1;
            index2 = UP_index2;
            [status_UP,interval,index] = InIntervals([slow_waves(1).UP_ints(UP_index1,1);slow_waves(2).UP_ints(UP_index2,1)],[ripples(1).onset-0.2 ripples(1).onset]);


        elseif contains(event_types{e},'DOWN')
            index1 = DOWN_index1;
            index2 = DOWN_index2;
        end

        if isempty(index1)
            TF_amp_V1_ipsi = squeeze(TF_amp_V1(2).(event_types{e})(2,:,:,index2));
            TF_amp_V1_ipsi = normalize_tf(TF_amp_V1_ipsi, baseline_idx, 'dB');
            % mean_TF_amp_V1_ipsi = mean(TF_amp_V1_ipsi,3,'omitnan');
            % [~, sig_mask_V1_ipsi] = permutation_TF_test(TF_amp_V1_ipsi, 1000);


            TF_amp_V1_contra = squeeze(TF_amp_V1(2).(event_types{e})(1,:,:,index2));
            TF_amp_V1_contra = normalize_tf(TF_amp_V1_contra, baseline_idx, 'dB');
            % mean_TF_amp_V1_contra = mean(TF_amp_V1_contra,3,'omitnan');
            % [~, sig_mask_V1_contra] = permutation_TF_test(TF_amp_V1_contra, 1000);


            TF_amp_HPC_ipsi = squeeze(TF_amp_HPC(2).(event_types{e})(2,:,:,index2));
            TF_amp_HPC_ipsi = normalize_tf(TF_amp_HPC_ipsi, baseline_idx, 'dB');
            % mean_TF_amp_V1_ipsi = mean(TF_amp_V1_ipsi,3,'omitnan');
            % [~, sig_mask_V1_ipsi] = permutation_TF_test(TF_amp_V1_ipsi, 1000);


            TF_amp_HPC_contra = squeeze(TF_amp_HPC(2).(event_types{e})(1,:,:,index2));
            TF_amp_HPC_contra = normalize_tf(TF_amp_HPC_contra, baseline_idx, 'dB');


            % === V1 Phase ===
            TF_phase_V1_ipsi = squeeze(cat(4, ...
                TF_phase_V1(2).(event_types{e})(2,:,:,index2)  ... % V1 R to R events
                ));  % size: [freq x time x events]

            TF_phase_V1_contra = squeeze(cat(4, ...
                TF_phase_V1(2).(event_types{e})(1,:,:,index2)  ... % V1 R to L events
                ));  % size: [freq x time x events]

            % === HPC Phase ===
            TF_phase_HPC_ipsi = squeeze(cat(4, ...
                TF_phase_HPC(2).(event_types{e})(2,:,:,index2)  ...
                ));

            TF_phase_HPC_contra = squeeze(cat(4, ...
                TF_phase_HPC(2).(event_types{e})(1,:,:,index2)  ...
                ));

        else
            TF_amp_V1_ipsi = squeeze([cat(4,TF_amp_V1(1).(event_types{e})(1,:,:,index1), TF_amp_V1(2).(event_types{e})(2,:,:,index2))]);
            TF_amp_V1_ipsi = normalize_tf(TF_amp_V1_ipsi, baseline_idx, 'dB');
            % mean_TF_amp_V1_ipsi = mean(TF_amp_V1_ipsi,3,'omitnan');
            % [~, sig_mask_V1_ipsi] = permutation_TF_test(TF_amp_V1_ipsi, 1000);


            TF_amp_V1_contra = squeeze([cat(4,TF_amp_V1(1).(event_types{e})(2,:,:,index1), TF_amp_V1(2).(event_types{e})(1,:,:,index2))]);
            TF_amp_V1_contra = normalize_tf(TF_amp_V1_contra, baseline_idx, 'dB');
            % mean_TF_amp_V1_contra = mean(TF_amp_V1_contra,3,'omitnan');
            % [~, sig_mask_V1_contra] = permutation_TF_test(TF_amp_V1_contra, 1000);


            TF_amp_HPC_ipsi = squeeze([cat(4,TF_amp_HPC(1).(event_types{e})(1,:,:,index1), TF_amp_HPC(2).(event_types{e})(2,:,:,index2))]);
            TF_amp_HPC_ipsi = normalize_tf(TF_amp_HPC_ipsi, baseline_idx, 'dB');
            % mean_TF_amp_HPC_ipsi = mean(TF_amp_HPC_ipsi,3,'omitnan');
            % [~, sig_mask_HPC_ipsi] = permutation_TF_test(TF_amp_HPC_ipsi, 1000);


            TF_amp_HPC_contra = squeeze([cat(4,TF_amp_HPC(1).(event_types{e})(2,:,:,index1), TF_amp_HPC(2).(event_types{e})(1,:,:,index2))]);
            TF_amp_HPC_contra = normalize_tf(TF_amp_HPC_contra, baseline_idx, 'dB');
            % mean_TF_amp_HPC_contra = mean(TF_amp_HPC_contra,3,'omitnan');
            % [~, sig_mask_HPC_contra] = permutation_TF_test(TF_amp_HPC_contra, 1000);




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
        end



        if contains(event_types{e},'ripples')

            Types = {'All','UP','DOWN'};
            num_figs = 3;
            % figure('Name',['TF amplitude - ' event_types{e}],'Position',[100 100 1200 500])
        elseif contains(event_types{e},'spindles')
            Types = {'All'};
            num_figs = 1;
            % figure('Name',['TF amplitude - ' event_types{e}],'Position',[100 100 1200 500])
        else
            Types = {'All','ripples'};
            num_figs = 2;
        end


        for n = 1:num_figs

            % Determine suffix and event_index
            if contains(Types{n}, 'UP') 
                suffix = '_UP';
                event_index = find(status_UP);
            elseif contains(Types{n}, 'DOWN') 
                suffix = '_DOWN';
                event_index = find(status_DOWN);
            elseif contains(event_types{e},'DOWN')&contains(Types{n}, 'ripples') 
                suffix = '_ripples';
                event_index = find(status_DOWN);
            elseif contains(event_types{e},'UP')&contains(Types{n}, 'ripples')
                suffix = '_ripples';
                event_index = find(status_UP);
            else
                suffix = '_all';
                event_index = 1:size(TF_amp_V1_ipsi,3);
            end

            % Dynamically get current event TF stats structure
            curr_TF_stats = eval([event_types{e}, '_TF_stats']);

            % === Compute mean maps ===
            mean_TF_amp_V1_ipsi    = mean(TF_amp_V1_ipsi(:,:,event_index), 3, 'omitnan');
            mean_TF_amp_V1_contra  = mean(TF_amp_V1_contra(:,:,event_index), 3, 'omitnan');
            mean_TF_amp_HPC_ipsi   = mean(TF_amp_HPC_ipsi(:,:,event_index), 3, 'omitnan');
            mean_TF_amp_HPC_contra = mean(TF_amp_HPC_contra(:,:,event_index), 3, 'omitnan');

            PLV_V1 = abs(mean(exp(1i * (TF_phase_V1_ipsi(:,:,event_index) - TF_phase_V1_contra(:,:,event_index))), 3, 'omitnan'));
            PLV_HPC = abs(mean(exp(1i * (TF_phase_HPC_ipsi(:,:,event_index) - TF_phase_HPC_contra(:,:,event_index))), 3, 'omitnan'));

            if contains(event_types{e}, 'ripples')
                PLV_V1_HPC_ipsi  = abs(mean(exp(1i * (TF_phase_HPC_ipsi(:,:,event_index) - TF_phase_V1_ipsi(:,:,event_index))), 3, 'omitnan'));
                PLV_V1_HPC_contra = abs(mean(exp(1i * (TF_phase_HPC_ipsi(:,:,event_index) - TF_phase_V1_contra(:,:,event_index))), 3, 'omitnan'));
            else
                PLV_V1_HPC_ipsi  = abs(mean(exp(1i * (TF_phase_V1_ipsi(:,:,event_index) - TF_phase_HPC_ipsi(:,:,event_index))), 3, 'omitnan'));
                PLV_V1_HPC_contra = abs(mean(exp(1i * (TF_phase_V1_ipsi(:,:,event_index) - TF_phase_HPC_contra(:,:,event_index))), 3, 'omitnan'));
            end


            %%%%%% amp corr
            %%%% V1 ipsi and contra
            nFreq = size(TF_amp_V1_ipsi,1);
            nTime = size(TF_amp_V1_ipsi,2);

            amp_corr_V1 = nan(nFreq, nTime);
            amp_ipsi = TF_amp_V1_ipsi(:,:,event_index);amp_contra = TF_amp_V1_contra(:,:,event_index);
            for f = 1:nFreq
                for t = 1:nTime
                    x = squeeze(amp_ipsi(f,t,:));     % [nEvent x 1]
                    y = squeeze(amp_contra(f,t,:));   % [nEvent x 1]

                    valid = ~isnan(x) & ~isnan(y);
                    if sum(valid) >= 10
                        amp_corr_V1(f,t) = corr(x(valid), y(valid), 'Type', 'Spearman');
                    end
                end
            end

            %%%% HPC ipsi and contra
            amp_corr_HPC = nan(nFreq, nTime);
            amp_ipsi = TF_amp_HPC_ipsi(:,:,event_index);amp_contra = TF_amp_HPC_contra(:,:,event_index);
            for f = 1:nFreq
                for t = 1:nTime
                    x = squeeze(amp_ipsi(f,t,:));     % [nEvent x 1]
                    y = squeeze(amp_contra(f,t,:));   % [nEvent x 1]

                    valid = ~isnan(x) & ~isnan(y);
                    if sum(valid) >= 10
                        amp_corr_HPC(f,t) = corr(x(valid), y(valid), 'Type', 'Spearman');
                    end
                end
            end

            %%%% V1-HPC ipsi
            amp_corr_V1_HPC = nan(2,nFreq, nTime);
            if contains(event_types{e}, 'ripples')
                amp_ipsi = TF_amp_HPC_ipsi(:,:,event_index);amp_contra = TF_amp_V1_ipsi(:,:,event_index);
            else
                amp_ipsi = TF_amp_V1_ipsi(:,:,event_index);amp_contra = TF_amp_HPC_ipsi(:,:,event_index);
            end

            for f = 1:nFreq
                for t = 1:nTime
                    x = squeeze(amp_ipsi(f,t,:));     % [nEvent x 1]
                    y = squeeze(amp_contra(f,t,:));   % [nEvent x 1]

                    valid = ~isnan(x) & ~isnan(y);
                    if sum(valid) >= 10
                        amp_corr_V1_HPC(1,f,t) = corr(x(valid), y(valid), 'Type', 'Spearman');
                    end
                end
            end


            %%%% V1-HPC contra
%             amp_corr_V1_HPC = nan(2,nFreq, nTime);
            if contains(event_types{e}, 'ripples')
                amp_ipsi = TF_amp_HPC_ipsi(:,:,event_index);amp_contra = TF_amp_V1_contra(:,:,event_index);
            else
                amp_ipsi = TF_amp_V1_ipsi(:,:,event_index);amp_contra = TF_amp_HPC_contra(:,:,event_index);
            end

            for f = 1:nFreq
                for t = 1:nTime
                    x = squeeze(amp_ipsi(f,t,:));     % [nEvent x 1]
                    y = squeeze(amp_contra(f,t,:));   % [nEvent x 1]

                    valid = ~isnan(x) & ~isnan(y);
                    if sum(valid) >= 10
                        amp_corr_V1_HPC(2,f,t) = corr(x(valid), y(valid), 'Type', 'Spearman');
                    end
                end
            end


            % === Store mean maps ===
            curr_TF_stats.(['mean_V1_amp', suffix]){nsession}(1,:,:) = mean_TF_amp_V1_ipsi;
            curr_TF_stats.(['mean_V1_amp', suffix]){nsession}(2,:,:) = mean_TF_amp_V1_contra;
            curr_TF_stats.(['mean_HPC_amp', suffix]){nsession}(1,:,:) = mean_TF_amp_HPC_ipsi;
            curr_TF_stats.(['mean_HPC_amp', suffix]){nsession}(2,:,:) = mean_TF_amp_HPC_contra;

            curr_TF_stats.(['mean_PLV_V1', suffix]){nsession} = PLV_V1;
            curr_TF_stats.(['mean_PLV_HPC', suffix]){nsession} = PLV_HPC;
            curr_TF_stats.(['mean_PLV_V1_HPC', suffix]){nsession}(1,:,:) = PLV_V1_HPC_ipsi;
            curr_TF_stats.(['mean_PLV_V1_HPC', suffix]){nsession}(2,:,:) = PLV_V1_HPC_contra;

            curr_TF_stats.(['mean_amp_corr_V1', suffix]){nsession} = amp_corr_V1;
            curr_TF_stats.(['mean_amp_corr_HPC', suffix]){nsession} = amp_corr_HPC;
            curr_TF_stats.(['mean_amp_corr_V1_HPC', suffix]){nsession} = amp_corr_V1_HPC;
%             curr_TF_stats.(['mean_amp_corr_V1_HPC', suffix]){nsession}(2,:,:) = PLV_V1_HPC_contra;


            if strcmp(suffix, '_all')
                bin_edges = linspace(-1, 1, 11);  % 10 bins
                nBins = length(bin_edges) - 1;
                nEvents = length(event_index);

                % Frequency bands
                filterparms = [];
                filterparms.deltafilter      = [1 4];
                filterparms.thetafilter      = [4 12];
                filterparms.spindlesfilter   = [9 17];
                filterparms.lowgammafilter   = [30 60];
                filterparms.highgammafilter  = [60 100];
                filterparms.ripplesfilter    = [125 300];
                band_names = fieldnames(filterparms);
                nBands = numel(band_names);

                % Preallocate: [hemi x band x bin x event]
                V1_amp_all = nan(2, nBands, nBins, nEvents);
                HPC_amp_all = nan(2, nBands, nBins, nEvents);
                V1_HPC_PLV_all = nan(2, nBands, nBins, nEvents);
                V1_PLV_all = nan(nBands, nBins, nEvents);
                HPC_PLV_all = nan(nBands, nBins, nEvents);
                V1_HPC_amp_corr_all = nan(2, nBands, nBins, nEvents);
                V1_amp_corr_all = nan(nBands, nBins, nEvents);
                HPC_amp_corr_all = nan(nBands, nBins, nEvents);

                V1_phase_median_all = nan(2, nBands, nBins, nEvents);  % [hemi x band x bin x event]
                HPC_phase_median_all = nan(2, nBands, nBins, nEvents);

                for f = 1:nBands
                    band = filterparms.(band_names{f});
                    freq_idx = find(freqs >= band(1) & freqs <= band(2));

                    for b = 1:nBins
                        tIdx = find(timevec >= bin_edges(b) & timevec < bin_edges(b+1));

                        for ev = 1:nEvents
                            ev_idx = event_index(ev);

                            % Amplitude
                            V1_amp_all(1,f,b,ev) = median(median(TF_amp_V1_ipsi(freq_idx, tIdx, ev_idx), 2, 'omitnan'), 'omitnan');
                            V1_amp_all(2,f,b,ev) = median(median(TF_amp_V1_contra(freq_idx, tIdx, ev_idx), 2, 'omitnan'), 'omitnan');

                            HPC_amp_all(1,f,b,ev) = median(median(TF_amp_HPC_ipsi(freq_idx, tIdx, ev_idx), 2, 'omitnan'), 'omitnan');
                            HPC_amp_all(2,f,b,ev) = median(median(TF_amp_HPC_contra(freq_idx, tIdx, ev_idx), 2, 'omitnan'), 'omitnan');

                            % Phase Median (circular)
                            % V1
                            V1_phase_ipsi = TF_phase_V1_ipsi(freq_idx, tIdx, ev_idx);
                            V1_phase_contra = TF_phase_V1_contra(freq_idx, tIdx, ev_idx);

                            % Take circular median over all [freq x time] points
                            V1_phase_median_ipsi = angle(median(exp(1i * V1_phase_ipsi(:)), 'omitnan'));
                            V1_phase_median_contra = angle(median(exp(1i * V1_phase_contra(:)), 'omitnan'));

                            % HPC
                            HPC_phase_ipsi = TF_phase_HPC_ipsi(freq_idx, tIdx, ev_idx);
                            HPC_phase_contra = TF_phase_HPC_contra(freq_idx, tIdx, ev_idx);

                            HPC_phase_median_ipsi = angle(median(exp(1i * HPC_phase_ipsi(:)), 'omitnan'));
                            HPC_phase_median_contra = angle(median(exp(1i * HPC_phase_contra(:)), 'omitnan'));

                            % Store
                            V1_phase_median_all(1,f,b,ev) = V1_phase_median_ipsi;
                            V1_phase_median_all(2,f,b,ev) = V1_phase_median_contra;

                            HPC_phase_median_all(1,f,b,ev) = HPC_phase_median_ipsi;
                            HPC_phase_median_all(2,f,b,ev) = HPC_phase_median_contra;

                            % PLV
                            if contains(event_types{e}, 'ripples')
                                phase_diff_ipsi = exp(1i * (TF_phase_HPC_ipsi(freq_idx, tIdx, ev_idx) - TF_phase_V1_ipsi(freq_idx, tIdx, ev_idx)));
                                phase_diff_contra = exp(1i * (TF_phase_HPC_ipsi(freq_idx, tIdx, ev_idx) - TF_phase_V1_contra(freq_idx, tIdx, ev_idx)));
                            else
                                phase_diff_ipsi = exp(1i * (TF_phase_V1_ipsi(freq_idx, tIdx, ev_idx) - TF_phase_HPC_ipsi(freq_idx, tIdx, ev_idx)));
                                phase_diff_contra = exp(1i * (TF_phase_V1_ipsi(freq_idx, tIdx, ev_idx) - TF_phase_HPC_contra(freq_idx, tIdx, ev_idx)));
                            end

                            V1_HPC_PLV_all(1,f,b,ev) = abs(mean(phase_diff_ipsi(:), 'omitnan'));
                            V1_HPC_PLV_all(2,f,b,ev) = abs(mean(phase_diff_contra(:), 'omitnan'));

                            % V1 PLV
                            ph_diff = exp(1i * (TF_phase_V1_ipsi(freq_idx, tIdx, ev_idx) - TF_phase_V1_contra(freq_idx, tIdx, ev_idx)));
                            V1_PLV_all(f,b,ev) = abs(mean(ph_diff(:), 'omitnan'));

                            % HPC PLV
                            ph_diff = exp(1i * (TF_phase_HPC_ipsi(freq_idx, tIdx, ev_idx) - TF_phase_HPC_contra(freq_idx, tIdx, ev_idx)));
                            HPC_PLV_all(f,b,ev) = abs(mean(ph_diff(:), 'omitnan'));



                            % amp corr
                            if contains(event_types{e}, 'ripples')
                                amp_corr_ipsi = corr(double(mean(TF_amp_HPC_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))',double(mean(TF_amp_V1_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))', 'Type', 'Spearman');
                                amp_corr_contra = corr(double(mean(TF_amp_HPC_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))',double(mean(TF_amp_V1_contra(freq_idx, tIdx, ev_idx), 1, 'omitnan'))', 'Type', 'Spearman');
                            else
                                amp_corr_ipsi = corr(double(mean(TF_amp_V1_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))',double(mean(TF_amp_HPC_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))', 'Type', 'Spearman');
                                amp_corr_contra = corr(double(mean(TF_amp_V1_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))',double(mean(TF_amp_HPC_contra(freq_idx, tIdx, ev_idx), 1, 'omitnan'))', 'Type', 'Spearman');
                            end

                            V1_HPC_amp_corr_all(1,f,b,ev) = amp_corr_ipsi;
                            V1_HPC_amp_corr_all(2,f,b,ev) =amp_corr_contra;

                            % V1 PLV
                            amp_corr_all = corr(double(mean(TF_amp_HPC_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))',double(mean(TF_amp_HPC_contra(freq_idx, tIdx, ev_idx), 1, 'omitnan'))', 'Type', 'Spearman');
                            V1_amp_corr_all(f,b,ev) = amp_corr_all;

                            % HPC PLV
                            amp_corr_all = corr(double(mean(TF_amp_V1_ipsi(freq_idx, tIdx, ev_idx), 1, 'omitnan'))',double(mean(TF_amp_V1_contra(freq_idx, tIdx, ev_idx), 1, 'omitnan'))', 'Type', 'Spearman');
                            HPC_amp_corr_all(f,b,ev) = amp_corr_all;
                        end
                    end
                end

                % Save into structure
                curr_TF_stats.V1_amp{nsession} = V1_amp_all;          % [2 x band x bin x event]
                curr_TF_stats.HPC_amp{nsession} = HPC_amp_all;

                curr_TF_stats.V1_HPC_PLV{nsession} = V1_HPC_PLV_all;
                curr_TF_stats.V1_PLV{nsession} = V1_PLV_all;
                curr_TF_stats.HPC_PLV{nsession} = HPC_PLV_all;

                curr_TF_stats.V1_HPC_amp_corr{nsession} = V1_HPC_amp_corr_all;
                curr_TF_stats.V1_amp_corr{nsession} = V1_amp_corr_all;
                curr_TF_stats.HPC_amp_corr{nsession} = HPC_amp_corr_all;

                curr_TF_stats.V1_phase_median{nsession} = V1_phase_median_all;
                curr_TF_stats.HPC_phase_median{nsession} = HPC_phase_median_all;


            end

            curr_TF_stats.freq{nsession} =    freqs;
            curr_TF_stats.timebin{nsession} =    timevec;
          
            % === Write updated structure back ===
            eval([event_types{e}, '_TF_stats = curr_TF_stats;']);

            fig = figure('Name', [options.SUBJECT, ' ', options.SESSION, ' ',options.StimulusName, ...
                ' TF amplitude - ', event_types{e}, ' ', Types{n}], ...
                'Position', [100 100 1200 500]);

            % --- V1 ipsi
            subplot(2,2,1)
            %             contourf(timevec, log2(freqs), mean_TF_amp_V1_ipsi, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), mean_TF_amp_V1_ipsi)
            axis xy; title(['Ipsi V1 - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            clim([-5 3])
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            % Sig mask
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_V1_amp_sig_mask', suffix]){nsession}(1,:,:);
            %             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)

            % --- V1 contra
            subplot(2,2,2)
            %             contourf(timevec, log2(freqs), mean_TF_amp_V1_contra, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), mean_TF_amp_V1_contra)
            axis xy; title(['Contra V1 - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([-2 2])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             hold on
%             sig_mask = curr_TF_stats.(['mean_V1_amp_sig_mask', suffix]){nsession}(2,:,:);
%             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)

            % --- HPC ipsi
            subplot(2,2,3)
%             contourf(timevec, log2(freqs), mean_TF_amp_HPC_ipsi, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), mean_TF_amp_HPC_ipsi)
            axis xy; title(['Ipsi HPC - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             hold on
%             sig_mask = curr_TF_stats.(['mean_HPC_amp_sig_mask', suffix]){nsession}(1,:,:);
%             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)

            % --- HPC contra
            subplot(2,2,4)
            imagesc(timevec, log2(freqs), mean_TF_amp_HPC_contra)
%             contourf(timevec, log2(freqs), mean_TF_amp_HPC_contra, 40, 'linecolor', 'none')
            axis xy; title(['Contra HPC - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([-2 5])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_HPC_amp_sig_mask', suffix]){nsession}(2,:,:);
            %             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)




            fig = figure('Name', [options.SUBJECT, ' ', options.SESSION, ' ',options.StimulusName, ...
                ' TF PLV - ', event_types{e}, ' ', Types{n}], ...
                'Position', [100 100 1200 500]);

            % --- V1 PLV
            subplot(2,2,1)
            %             contourf(timevec, log2(freqs), PLV_V1, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), PLV_V1)
            axis xy; title(['V1 PLV - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
            %             contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

            % --- HPC PLV
            subplot(2,2,2)
            %             contourf(timevec, log2(freqs), PLV_HPC, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), PLV_HPC)
            axis xy; title(['HPC PLV - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_HPC', suffix]){nsession};
            %             contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

            % --- V1-HPC Ipsi
            subplot(2,2,3)
            %             contourf(timevec, log2(freqs), PLV_V1_HPC_ipsi, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), PLV_V1_HPC_ipsi)

            axis xy; title(['Ipsi V1-HPC PLV - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1_HPC', suffix]){nsession}(1,:,:);
            %             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)

            % --- V1-HPC Contra
            subplot(2,2,4)
            %             contourf(timevec, log2(freqs), PLV_V1_HPC_contra, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), PLV_V1_HPC_contra)
            axis xy; title(['Contra V1-HPC PLV - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1_HPC', suffix]){nsession}(2,:,:);
            %             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)



            fig = figure('Name', [options.SUBJECT, ' ', options.SESSION, ' ',options.StimulusName, ...
                ' TF amp corr - ', event_types{e}, ' ', Types{n}], ...
                'Position', [100 100 1200 500]);

            % --- V1 amp
            subplot(2,2,1)
            %             contourf(timevec, log2(freqs), PLV_V1, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), amp_corr_V1)
            axis xy; title(['V1 amp corr - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
            %             contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

            % --- HPC amp
            subplot(2,2,2)
            %             contourf(timevec, log2(freqs), PLV_HPC, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), amp_corr_HPC)
            axis xy; title(['HPC amp corr - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_HPC', suffix]){nsession};
            %             contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

            % --- V1-HPC Ipsi amp
            subplot(2,2,3)
            %             contourf(timevec, log2(freqs), PLV_V1_HPC_ipsi, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), squeeze(amp_corr_V1_HPC(1,:,:)))

            axis xy; title(['Ipsi V1-HPC amp corr - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1_HPC', suffix]){nsession}(1,:,:);
            %             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)

            % --- V1-HPC Contra amp
            subplot(2,2,4)
            %             contourf(timevec, log2(freqs), PLV_V1_HPC_contra, 40, 'linecolor', 'none')
            imagesc(timevec, log2(freqs), squeeze(amp_corr_V1_HPC(2,:,:)))
            axis xy; title(['Contra V1-HPC amp corr - ' event_types{e}])
            xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
            set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
                'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
            colorbar
            xline(0,'r--')
            xlim([-1.5 1.5])
            ylim([min(log2(freqs)) max(log2(freqs))])
            clim([0 0.6])
            yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
            set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
            %             hold on
            %             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1_HPC', suffix]){nsession}(2,:,:);
            %             contour(timevec, log2(freqs), squeeze(sig_mask), [1 1], 'r', 'LineWidth', 1.2)

        end



        if  contains(options.StimulusName,'RUN1')|contains(options.StimulusName,'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(options.StimulusName,'Masa2tracks_'))))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(options.StimulusName,'Masa2tracks_'))),[])
        else
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',options.StimulusName)))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',options.StimulusName)),[])
        end

%         eval([event_types{e}, '_TF_stats = curr_TF_stats;']);
    end
   toc
end
% 
% ripples_TF_stats = struct();
% UP_TF_stats = struct();
% DOWN_TF_stats = struct();
% spindles_TF_stats = struct();

% ripples_all = rmfield(ripples_all, 'detectorinfo');
% spindles_all = rmfield(spindles_all, 'detectorinfo');
if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

if contains(Stimulus_type,'Sleep') & ~contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'ripples_TF_stats_POST.mat'),'ripples_TF_stats','-v7.3')
    save(fullfile(analysis_folder,'UP_TF_stats_POST.mat'),'UP_TF_stats','-v7.3')
    save(fullfile(analysis_folder,'DOWN_TF_stats_POST.mat'),'DOWN_TF_stats','-v7.3')
    save(fullfile(analysis_folder,'spindles_TF_stats_POST.mat'),'spindles_TF_stats','-v7.3')
elseif contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'ripples_TF_stats_PRE.mat'),'ripples_TF_stats','-v7.3')
    save(fullfile(analysis_folder,'UP_TF_stats_PRE.mat'),'UP_TF_stats','-v7.3')
    save(fullfile(analysis_folder,'DOWN_TF_stats_PRE.mat'),'DOWN_TF_stats','-v7.3')
    save(fullfile(analysis_folder,'spindles_TF_stats_PRE.mat'),'spindles_TF_stats','-v7.3')
else
%     save(fullfile(analysis_folder,'slow_waves_all.mat'),'slow_waves_all','-v7.3')
%     save(fullfile(analysis_folder,'ripples_all.mat'),'ripples_all','-v7.3')
%     save(fullfile(analysis_folder,'spindles_all.mat'),'spindles_all','-v7.3')
%     save(fullfile(analysis_folder,'behavioural_state_merged_all.mat'),'behavioural_state_merged_all','-v7.3')
end

% 