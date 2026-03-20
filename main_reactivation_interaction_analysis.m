%% Extract and Process reactivation data


%% Extract Log odds bayesian bias

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

session_count = 0;

bayesian_reactivation_all= struct();
bayesian_reactivation_V1_all= struct();

bin_size = 0.02;
time_edges = -1:bin_size:1;
time_centres = time_edges(1)+bin_size/2:bin_size:time_edges(end)-bin_size/2;

for nsession =1:length(experiment_info)

    tic
    disp(sprintf('session %i',nsession))
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    if length(stimulus_name)>1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2'));
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

    if contains(stimulus_name{n},'Sleep')
        load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'),'decoded_ripple_events');
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_shuffled.mat'),'decoded_ripple_events_shuffled');
        load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_V1.mat'),'decoded_ripple_events_V1');
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_V1_shuffled.mat'),'decoded_ripple_events_V1_shuffled');
    end

    for nprobe = 1:length(decoded_ripple_events_V1)

        % Store bin quality and threshold
        bayesian_reactivation_V1_all(nprobe).good_bins{nsession} = decoded_ripple_events_V1(nprobe).good_bins;
        % bayesian_reactivation_V1_all(nprobe).PV_threshold{nsession} = PV_threshold;

        % -- Compute bias and log-odds distributions for z-scoring
        summed_probability_distribution = [
            decoded_ripple_events_V1(nprobe).track(1).replay_events(:).summed_probability;
            decoded_ripple_events_V1(nprobe).track(2).replay_events(:).summed_probability
            ];
        bias_distribution = summed_probability_distribution(1,:) ./ sum(summed_probability_distribution, 1, 'omitnan');
        bias_mean = mean(bias_distribution, 'omitnan');
        bias_std = std(bias_distribution, 0, 'omitnan');

        log_odds_distribution = log(summed_probability_distribution(1,:) ./ summed_probability_distribution(2,:));

        Capping = max([prctile(log_odds_distribution,99.5) abs(prctile(log_odds_distribution,0.5))]);

        log_odds_distribution(log_odds_distribution>Capping)=Capping;
        log_odds_distribution(log_odds_distribution<-Capping)=-Capping;

        log_odds_mean = mean(log_odds_distribution, 'omitnan');
        log_odds_std = std(log_odds_distribution, 0, 'omitnan');
      
        % -- Initialize accumulation
        if session_count == 1 || ~isfield(bayesian_reactivation_V1_all(nprobe), 'decoded_matrix')
            bayesian_reactivation_V1_all(nprobe).decoded_matrix = [];
            bayesian_reactivation_V1_all(nprobe).summed_probability = [];
            % bayesian_reactivation_V1_all(nprobe).bias = [];
            % bayesian_reactivation_V1_all(nprobe).z_bias = [];
            bayesian_reactivation_V1_all(nprobe).log_odds = [];
            bayesian_reactivation_V1_all(nprobe).z_log_odds = [];
            bayesian_reactivation_V1_all(nprobe).z_log_odds_shuffled = [];
            % bayesian_reactivation_V1_all(nprobe).log_odds_percentile = [];
        end

        for nevent = 1:length(decoded_ripple_events_V1(nprobe).track(1).replay_events)
            % Decode info
            replay1 = decoded_ripple_events_V1(nprobe).track(1).replay_events(nevent).replay;
            replay2 = decoded_ripple_events_V1(nprobe).track(2).replay_events(nevent).replay;
            sumprob1 = decoded_ripple_events_V1(nprobe).track(1).replay_events(nevent).summed_probability;
            sumprob2 = decoded_ripple_events_V1(nprobe).track(2).replay_events(nevent).summed_probability;

            % Bin alignment
            time_edges = decoded_ripple_events_V1(nprobe).track(1).replay_events(nevent).timebins_edges;
            onset_time = decoded_ripple_events_V1(nprobe).track(2).replay_events(nevent).onset;
            [~, onset_bin] = min(abs(time_edges - onset_time));

            target_bins = 100;
            center_bin  = 51;

            [npos, nbins] = size(replay1);

            % Map original bins -> new 1:100 bins so that onset_bin -> 51
            src_bins = 1:nbins;                         % original bin indices
            shift    = center_bin - onset_bin;          % e.g. 51 - onset_bin
            dst_bins = src_bins + shift;                % where each bin will land

            % Only keep bins that land inside 1:100
            valid_idx = dst_bins >= 1 & dst_bins <= target_bins;

            % ---- Replays (npos x nbins -> npos x 100) ----
            mat1 = nan(npos, target_bins);
            mat2 = nan(npos, target_bins);

            mat1(:, dst_bins(valid_idx)) = replay1(:, valid_idx);
            mat2(:, dst_bins(valid_idx)) = replay2(:, valid_idx);

            combined_matrix = cat(1, mat1, mat2);

            bayesian_reactivation_V1_all(nprobe).decoded_matrix(:,:,end+1) = combined_matrix;
            % bin_range = [];
            sp_nan = nan(2, 100);
            sp_nan(1, dst_bins(valid_idx)) = sumprob1(valid_idx);
            sp_nan(2, dst_bins(valid_idx)) = sumprob2(valid_idx);
            bayesian_reactivation_V1_all(nprobe).summed_probability(:,:,end+1) = sp_nan;

            % % Bias
            % total_prob = sum(sp_nan, 1, 'omitnan');
            % bias = sp_nan(1,:) ./ total_prob;
            % z_bias = (bias - bias_mean) ./ bias_std;
            % bayesian_reactivation_V1_all(nprobe).bias(:,end+1) = bias;
            % bayesian_reactivation_V1_all(nprobe).z_bias(:,end+1) = z_bias;

            % Log odds and shuffled z-score
            % shuffled1 = decoded_ripple_events_V1_shuffled(nprobe).track(1).replay_events(nevent).summed_probability;
            % shuffled2 = decoded_ripple_events_V1_shuffled(nprobe).track(2).replay_events(nevent).summed_probability;
            % sp_shuf = log(shuffled1 ./ shuffled2);
            %
            % shuf_nan = nan(size(sp_shuf,1), 50);
            % for s = 1:size(sp_shuf,1)
            %     shuf_nan(s, bin_range(valid_idx)) = sp_shuf(s, valid_idx);
            % end
            %
            data = log(sp_nan(1,:) ./ sp_nan(2,:));
            data(data>Capping)=Capping;
            data(Capping<-Capping)=-Capping;

            bayesian_reactivation_V1_all(nprobe).z_log_odds(:,end+1) = ...
                (data - log_odds_mean) ./ log_odds_std;

            z_shuf_full = nan(1,target_bins);
            z_shuf_full(dst_bins(valid_idx)) = decoded_ripple_events_V1(nprobe).track(1).replay_events(nevent).z_log_odds(valid_idx);

            bayesian_reactivation_V1_all(nprobe).z_log_odds_shuffled(:,end+1) = z_shuf_full;

        end
    end


    %     save(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events_V1.mat'),'decoded_ripple_events_V1');


    %%%%%%%% HPC
    for nprobe = 1:length(decoded_ripple_events)

        % Store bin quality and threshold
        bayesian_reactivation_all(nprobe).good_bins{nsession} = decoded_ripple_events(nprobe).good_bins;
        % bayesian_reactivation_all(nprobe).PV_threshold{nsession} = PV_threshold;

        % -- Compute bias and log-odds distributions for z-scoring
        summed_probability_distribution = [
            decoded_ripple_events(nprobe).track(1).replay_events(:).summed_probability;
            decoded_ripple_events(nprobe).track(2).replay_events(:).summed_probability
            ];
        bias_distribution = summed_probability_distribution(1,:) ./ ...
            sum(summed_probability_distribution, 1, 'omitnan');
        bias_mean = mean(bias_distribution, 'omitnan');
        bias_std  = std(bias_distribution, 0, 'omitnan');


        log_odds_distribution = log( summed_probability_distribution(1,:) ./ ...
            summed_probability_distribution(2,:) );
        Capping = max([prctile(log_odds_distribution,99.5) abs(prctile(log_odds_distribution,0.5))]);

        log_odds_distribution(log_odds_distribution>Capping)=Capping;
        log_odds_distribution(log_odds_distribution<-Capping)=-Capping;

        log_odds_mean = mean(log_odds_distribution, 'omitnan');
        log_odds_std = std(log_odds_distribution, 0, 'omitnan');

        % -- Initialize accumulation
        if session_count == 1 || ~isfield(bayesian_reactivation_all(nprobe), 'decoded_matrix')
            bayesian_reactivation_all(nprobe).decoded_matrix = [];
            bayesian_reactivation_all(nprobe).summed_probability = [];
            % bayesian_reactivation_all(nprobe).bias = [];
            % bayesian_reactivation_all(nprobe).z_bias = [];
            % bayesian_reactivation_all(nprobe).log_odds = [];
            bayesian_reactivation_all(nprobe).z_log_odds = [];
            bayesian_reactivation_all(nprobe).z_log_odds_shuffled = [];
            % bayesian_reactivation_all(nprobe).log_odds_percentile = [];
        end

        for nevent = 1:length(decoded_ripple_events(nprobe).track(1).replay_events)
            % Decode info
            replay1  = decoded_ripple_events(nprobe).track(1).replay_events(nevent).replay;
            replay2  = decoded_ripple_events(nprobe).track(2).replay_events(nevent).replay;
            sumprob1 = decoded_ripple_events(nprobe).track(1).replay_events(nevent).summed_probability;
            sumprob2 = decoded_ripple_events(nprobe).track(2).replay_events(nevent).summed_probability;

            % Bin alignment
            time_edges = decoded_ripple_events(nprobe).track(1).replay_events(nevent).timebins_edges;
            onset_time = decoded_ripple_events(nprobe).track(2).replay_events(nevent).onset;
            [~, onset_bin] = min(abs(time_edges - onset_time));

            target_bins = 100;
            center_bin  = 51;

            [npos, nbins] = size(replay1);

            % Map original bins -> new 1:100 bins so that onset_bin -> 51
            src_bins = 1:nbins;                % original bin indices
            shift    = center_bin - onset_bin; % e.g. 51 - onset_bin
            dst_bins = src_bins + shift;       % where each bin will land

            % Only keep bins that land inside 1:100
            valid_idx = dst_bins >= 1 & dst_bins <= target_bins;

            % ---- Replays (npos x nbins -> npos x 100) ----
            mat1 = nan(npos, target_bins);
            mat2 = nan(npos, target_bins);

            mat1(:, dst_bins(valid_idx)) = replay1(:, valid_idx);
            mat2(:, dst_bins(valid_idx)) = replay2(:, valid_idx);

            combined_matrix = cat(1, mat1, mat2);
            bayesian_reactivation_all(nprobe).decoded_matrix(:,:,end+1) = combined_matrix;

            sp_nan = nan(2, target_bins);
            sp_nan(1, dst_bins(valid_idx)) = sumprob1(valid_idx);
            sp_nan(2, dst_bins(valid_idx)) = sumprob2(valid_idx);
            bayesian_reactivation_all(nprobe).summed_probability(:,:,end+1) = sp_nan;

            % % Bias
            % total_prob = sum(sp_nan, 1, 'omitnan');
            % bias  = sp_nan(1,:) ./ total_prob;
            % z_bias = (bias - bias_mean) ./ bias_std;
            % 
            % bayesian_reactivation_all(nprobe).bias(:,end+1)   = bias;
            % bayesian_reactivation_all(nprobe).z_bias(:,end+1) = z_bias;

            % Log odds and shuffled z-score
            data = log(sp_nan(1,:) ./ sp_nan(2,:));
            data(data>Capping)=Capping;
            data(Capping<-Capping)=-Capping;

            bayesian_reactivation_all(nprobe).z_log_odds(:,end+1) = ...
                (data - log_odds_mean) ./ log_odds_std;

            % Shuffled z-log-odds, padded/shifted into 1x100
            z_shuf_full = nan(1, target_bins);
            z_shuf_src  = decoded_ripple_events(nprobe).track(1).replay_events(nevent).z_log_odds;
            z_shuf_src  = z_shuf_src(:).';  % ensure row

            z_shuf_full(1, dst_bins(valid_idx)) = z_shuf_src(1, valid_idx);

            bayesian_reactivation_all(nprobe).z_log_odds_shuffled(:,end+1) = z_shuf_full;

        end
    end

    %     save(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'),'decoded_ripple_events');
end

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

save(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'),'bayesian_reactivation_all')
save(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'),'bayesian_reactivation_V1_all')





%% Grab PLS KDE RUN vaidation summary
clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'Sleep';
% [1 2 3 4 9 10 12 14]
% Stimulus_types_all = {'RUN'};
% Stimulus_types_all = {'RUN','POST'};
KDE_RUN_validations_all = struct();
% KDE_RUN_validations_POST_all = struc();

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'KDE_RUN_validations*.mat'));
        % DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'KDE_RUN_validations_all_bins*.mat'));
        
        if isempty(DIR)
            continue
        end

        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_RUN_validations_V1.mat'),'KDE_RUN_validations_V1');
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_RUN_validations.mat'),'KDE_RUN_validations');
        % load(fullfile(options.ANALYSIS_DATAPATH,'KDE_RUN_validations_V1_all_bins.mat'),'KDE_RUN_validations_V1');
        % load(fullfile(options.ANALYSIS_DATAPATH,'KDE_RUN_validations_all_bins.mat'),'KDE_RUN_validations');
        % 20ms ROC AUC decoding performance
        index = 2;
        % index = 1;
        KDE_RUN_validations_all.FPR(nsession,:) = (KDE_RUN_validations{index}.FPR);

        KDE_RUN_validations_all.HPC_AUC(1,nsession) = mean(KDE_RUN_validations{index}.AUC_real);
        KDE_RUN_validations_all.HPC_AUC_shuffled(1,nsession) = prctile(KDE_RUN_validations{index}.AUC_shuffle,97.5);
        % KDE_RUN_validations_all.HPC_AUC_shuffled(1,nsession) = mean(KDE_RUN_validations{index}.AUC_shuffle);

        KDE_RUN_validations_all.HPC_TPR(1,nsession,:) = mean(KDE_RUN_validations{index}.TPR_real);
        % KDE_RUN_validations_all.HPC_TPR_shuffled(1,nsession,:) = mean(KDE_RUN_validations{index}.TPR_shuffle);
        KDE_RUN_validations_all.HPC_TPR_shuffled(1,nsession,:) = prctile(KDE_RUN_validations{index}.TPR_shuffle,97.5);
        
        KDE_RUN_validations_all.V1_AUC(1,nsession) = mean(KDE_RUN_validations_V1{index}.AUC_real);
        % KDE_RUN_validations_all.V1_AUC_shuffled(1,nsession) = mean(KDE_RUN_validations_V1{index}.AUC_shuffle);
        KDE_RUN_validations_all.V1_AUC_shuffled(1,nsession) = prctile(KDE_RUN_validations_V1{index}.AUC_shuffle,97.5);

        KDE_RUN_validations_all.V1_TPR(1,nsession,:) = mean(KDE_RUN_validations_V1{index}.TPR_real);
        % KDE_RUN_validations_all.V1_TPR_shuffled(1,nsession,:) = mean(KDE_RUN_validations_V1{index}.TPR_shuffle);
        KDE_RUN_validations_all.V1_TPR_shuffled(1,nsession,:) = prctile(KDE_RUN_validations_V1{index}.TPR_shuffle,97.5);
        
        % 100ms ROC AUC decoding performance
        KDE_RUN_validations_all.HPC_AUC(2,nsession) = mean(KDE_RUN_validations{end}.AUC_real);
        % KDE_RUN_validations_all.HPC_AUC_shuffled(2,nsession) = mean(KDE_RUN_validations{end}.AUC_shuffle);
        KDE_RUN_validations_all.HPC_AUC_shuffled(2,nsession) = prctile(KDE_RUN_validations{end}.AUC_shuffle,97.5);

        KDE_RUN_validations_all.HPC_TPR(2,nsession,:) = mean(KDE_RUN_validations{end}.TPR_real);
        % KDE_RUN_validations_all.HPC_TPR_shuffled(2,nsession,:) = mean(KDE_RUN_validations{end}.TPR_shuffle);
        KDE_RUN_validations_all.HPC_TPR_shuffled(2,nsession,:) = prctile(KDE_RUN_validations{end}.TPR_shuffle,97.5);


        KDE_RUN_validations_all.V1_AUC(2,nsession) = mean(KDE_RUN_validations_V1{end}.AUC_real);
        % KDE_RUN_validations_all.V1_AUC_shuffled(2,nsession) = mean(KDE_RUN_validations_V1{end}.AUC_shuffle);
        KDE_RUN_validations_all.V1_AUC_shuffled(2,nsession) = prctile(KDE_RUN_validations_V1{end}.AUC_shuffle,97.5);

        KDE_RUN_validations_all.V1_TPR(2,nsession,:) = mean(KDE_RUN_validations_V1{end}.TPR_real);
        % KDE_RUN_validations_all.V1_TPR_shuffled(2,nsession,:) = mean(KDE_RUN_validations_V1{end}.TPR_shuffle);
        KDE_RUN_validations_all.V1_TPR_shuffled(2,nsession,:) = prctile(KDE_RUN_validations_V1{end}.TPR_shuffle,97.5);
    end
end

% spindles_all = rmfield(spindles_all, 'detectorinfo');
if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
% save(fullfile(analysis_folder,'KDE_RUN_validations_all_bins_all_95percentile.mat'),'KDE_RUN_validations_all','-v7.3')
save(fullfile(analysis_folder,'KDE_RUN_validations_all_95percentile.mat'),'KDE_RUN_validations_all','-v7.3')


% save(fullfile(analysis_folder,'KDE_RUN_validations_cell_id_all.mat'),'KDE_RUN_validations_all','-v7.3')
save(fullfile(analysis_folder,'KDE_RUN_validations_all_bins_all.mat'),'KDE_RUN_validations_all','-v7.3')
% save(fullfile(analysis_folder,'KDE_RUN_validations_all.mat'),'KDE_RUN_validations_all')


%% Plotting ROC AUC of RUN 2 track discrimination 
% spindles_all = rmfield(spindles_all, 'detectorinfo');
if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
% load(fullfile(analysis_folder,'KDE_RUN_validations_all.mat'),'KDE_RUN_validations_all')
load(fullfile(analysis_folder,'KDE_RUN_validations_all_bins_all.mat'),'KDE_RUN_validations_all')

load(fullfile(analysis_folder,'KDE_RUN_validations_all_bins_all_95percentile.mat'),'KDE_RUN_validations_all')
% load(fullfile(analysis_folder,'KDE_RUN_validations_all_95percentile.mat'),'KDE_RUN_validations_all')


timebins = [20,100];
% title_text = 'PLS contextual discrimination HPC RUN1';
% title_text = 'PLS contextual discrimination HPC RUN1 scatter';
% title_text = 'PLS contextual discrimination HPC RUN1 scatter (all bins)';
% title_text = 'PLS contextual discrimination HPC RUN1 scatter (95th percentile)';
title_text = 'PLS contextual discrimination HPC RUN1 scatter (all bins 95th percentile)';

fig = figure('Name', title_text, 'Position', [200 100 640 580]); hold on;

% % title_text = 'PLS contextual discrimination RUN1';
% for n = 1:2
%     nexttile
%     fpr = KDE_RUN_validations_all.FPR(1,:);
%     TPR_real = squeeze(mean(KDE_RUN_validations_all.HPC_TPR(n,:,:),2))';
%     TPR_shuf = squeeze(mean(KDE_RUN_validations_all.HPC_TPR_shuffled(n,:,:),2))';
%     CI_real = (std(squeeze(KDE_RUN_validations_all.HPC_TPR(n,:,:)),'omitnan'))./sqrt(size(KDE_RUN_validations_all.HPC_TPR,2));
%     CI_shuf = (std(squeeze(KDE_RUN_validations_all.HPC_TPR_shuffled(n,:,:)),'omitnan'))./sqrt(size(KDE_RUN_validations_all.HPC_TPR_shuffled,2));
%     hold on
%     plot([0 1],[0 1],'--k')
%     PLOT(1) = fill([fpr fliplr(fpr)], [TPR_real+CI_real fliplr(TPR_real-CI_real)], ...
%         [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     PLOT(2) = fill([fpr fliplr(fpr)], [TPR_shuf+CI_shuf fliplr(TPR_shuf-CI_shuf)], ...
%         'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     plot(fpr, TPR_real, 'Color', [231,41,138]/256);
%     plot(fpr, TPR_shuf, 'k', 'LineWidth', 1.5);
%     xlabel('True positive rate')
%     ylabel('False Positive rate')
%     %     xline(0.5, '--k'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
%     title(sprintf('ROC curve HPC %i ms bins',timebins(n))); legend(PLOT(1:2),{ 'Real', 'Shuffle'},'box', 'off');
%     set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% end

line_w = 0.2;

for n = 1:2
    nexttile
    % Extract base data
    fpr = KDE_RUN_validations_all.FPR(1,:);
    % These should be (sessions x fpr_points)
    data_real = squeeze(KDE_RUN_validations_all.HPC_TPR(n,:,:));
    data_shuf = squeeze(KDE_RUN_validations_all.HPC_TPR_shuffled(n,:,:));
    
    % Calculate Means
    TPR_real_mean = mean(data_real, 1, 'omitnan');
    TPR_shuf_mean = mean(data_shuf, 1, 'omitnan');
    
    hold on
    % 1. Plot the diagonal chance line
    plot([0 1],[0 1],'--k', 'HandleVisibility', 'off')
    
    % 2. Plot individual session lines (Low Alpha)
    % Color definitions
    real_col = [231,41,138]/255;
    shuf_col = [0, 0, 0]/255;
    
    % Plot all individual Real lines
    p_indiv_real = plot(fpr, data_real', 'Color', [real_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % Plot all individual Shuffled lines
    p_indiv_shuf = plot(fpr, data_shuf', 'Color', [shuf_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % 3. Plot Mean lines (Thick)
    PLOT(1) = plot(fpr, TPR_real_mean, 'Color', real_col, 'LineWidth', 2.5);
    PLOT(2) = plot(fpr, TPR_shuf_mean, 'Color', shuf_col, 'LineWidth', 2.5);
    
    % Formatting
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    title(sprintf('ROC curve HPC %i ms bins', timebins(n)));
    
    legend(PLOT(1:2), {'Real Mean', 'Shuffle Mean'}, 'Location', 'southeast', 'Box', 'off');
    
    xlim([0 1]); ylim([0 1]);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end



for n = 1:2
    nexttile
    % Extract data for the current condition
    real_data = KDE_RUN_validations_all.HPC_AUC(n,:);
    shuffled_data = KDE_RUN_validations_all.HPC_AUC_shuffled(n,:);
    
    % Combine for easier mean calculation
    % AUC = [real_data; shuffled_data];
    mean_real = mean(real_data, 'omitnan');
    mean_shuf = mean(shuffled_data, 'omitnan');
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos = [1 2];

    hold on

    % 1. Plot the Bars (Mean only, no error bars)
    plot([x_pos(1)-line_w, x_pos(1)+line_w], [mean_real, mean_real], 'Color', bar_colors(1,:), 'LineWidth', 4);
    plot([x_pos(2)-line_w, x_pos(2)+line_w], [mean_shuf, mean_shuf], 'Color', bar_colors(2,:), 'LineWidth', 4);

    % 2. Plot Raw Data Points and Connections with Jitter
    rng(1); % Seed for consistent jitter appearance
    jitter_strength = 0.2; % Adjusted for better visibility of distributions
    
    for s = 1:length(real_data)
        if ~isnan(real_data(s)) && ~isnan(shuffled_data(s))
            % Calculate jittered x-positions
            % This spreads the dots randomly around the x_pos center
            x1 = x_pos(1) + (rand - 0.5) * jitter_strength;
            x2 = x_pos(2) + (rand - 0.5) * jitter_strength;
            
            % Plot the connection line first (so it's behind the dots)
            plot([x1, x2], [real_data(s), shuffled_data(s)], 'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.5);
            
            % Plot the individual dots
            scatter(x1, real_data(s), 50, bar_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.3);
            scatter(x2, shuffled_data(s), 50, bar_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.3);
        end
    end

    % Statistics and Formatting
    [p, h] = signrank(real_data, shuffled_data, 'tail', 'right');
    
    xlim([-0.5 3.5])
    xticks([1 2]);
    xticklabels({'Real', 'Shuffled'});
    ylabel('AUC');
    ylim([0 1.1]);
    yticks([0:0.25:1])
    title(['Condition ' num2str(n) ' (p=' num2str(p, '%.4e') ')']);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end
% for n = 1:2
%     nexttile
% 
%     AUC = [KDE_RUN_validations_all.HPC_AUC(n,:);KDE_RUN_validations_all.HPC_AUC_shuffled(n,:)];
%     bar_colors = [231,41,138; 0, 0, 0]/255;
%     x_pos=[1 2];
%     for i = 1:2
% 
%         mean_AUC = mean(AUC(i,:),'omitnan');
%         auc_errors = [std(AUC(i,:),'omitnan')./sqrt(size(AUC(i,:),2)),...  % lower error
%             std(AUC(i,:),'omitnan')./sqrt(size(AUC(i,:),2))];   % upper error
%         hold on
%         bar(x_pos(i), mean_AUC, 0.4, 'FaceAlpha',0.5, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
%         errorbar(x_pos(i), mean_AUC, auc_errors(1), auc_errors(2), 'k', 'linestyle', 'none', 'linewidth', 1);
%     end
% 
%     [p,h,stats] = signrank(AUC(1,:),AUC(2,:),'tail','right');
% 
% 
%     xticks([1 2]);
%     xticklabels({'Real', 'Shuffled'});
%     ylabel('AUC');
%     ylim([0 1]);
%     title('Mean AUC – Real vs Shuffled');
%     set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% end


%%%%%%%%% V1
timebins = [20,100];
% title_text = 'PLS contextual discrimination V1 RUN1';
%  title_text = 'PLS contextual discrimination V1 RUN1 scatter';
 % title_text = 'PLS contextual discrimination V1 RUN1 scatter (all bins)';
 % title_text = 'PLS contextual discrimination V1 RUN1 scatter (95th percentile)';
title_text = 'PLS contextual discrimination V1 RUN1 scatter (all bins 95th percentile)';
fig = figure('Name', title_text, 'Position', [200 100 640 580]); hold on;

% % title_text = 'PLS contextual discrimination RUN1';
% for n = 1:2
%     nexttile
%     fpr = KDE_RUN_validations_all.FPR(1,:);
%     TPR_real = squeeze(mean(KDE_RUN_validations_all.V1_TPR(n,:,:),2))';
%     TPR_shuf = squeeze(mean(KDE_RUN_validations_all.V1_TPR_shuffled(n,:,:),2))';
%     CI_real = (std(squeeze(KDE_RUN_validations_all.V1_TPR(n,:,:)),'omitnan'))./sqrt(size(KDE_RUN_validations_all.V1_TPR,2));
%     CI_shuf = (std(squeeze(KDE_RUN_validations_all.V1_TPR_shuffled(n,:,:)),'omitnan'))./sqrt(size(KDE_RUN_validations_all.V1_TPR_shuffled,2));
%     hold on
%     plot([0 1],[0 1],'--k')
%     PLOT(1) = fill([fpr fliplr(fpr)], [TPR_real+CI_real fliplr(TPR_real-CI_real)], ...
%         [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     PLOT(2) = fill([fpr fliplr(fpr)], [TPR_shuf+CI_shuf fliplr(TPR_shuf-CI_shuf)], ...
%         'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     plot(fpr, TPR_real, 'Color', [231,41,138]/256);
%     plot(fpr, TPR_shuf, 'k', 'LineWidth', 1.5);
%     xlabel('True positive rate')
%     ylabel('False Positive rate')
%     %     xline(0.5, '--k'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
%     title(sprintf('ROC curve V1 %i ms bins',timebins(n))); legend(PLOT(1:2),{ 'Real', 'Shuffle'},'box', 'off');
%     set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% end


for n = 1:2
    nexttile
    % Extract base data
    fpr = KDE_RUN_validations_all.FPR(1,:);
    % These should be (sessions x fpr_points)
    data_real = squeeze(KDE_RUN_validations_all.V1_TPR(n,:,:));
    data_shuf = squeeze(KDE_RUN_validations_all.V1_TPR_shuffled(n,:,:));
    
    % Calculate Means
    TPR_real_mean = mean(data_real, 1, 'omitnan');
    TPR_shuf_mean = mean(data_shuf, 1, 'omitnan');
    
    hold on
    % 1. Plot the diagonal chance line
    plot([0 1],[0 1],'--k', 'HandleVisibility', 'off')
    
    % 2. Plot individual session lines (Low Alpha)
    % Color definitions
    real_col = [231,41,138]/255;
    shuf_col = [0, 0, 0]/255;
    
    % Plot all individual Real lines
    p_indiv_real = plot(fpr, data_real', 'Color', [real_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % Plot all individual Shuffled lines
    p_indiv_shuf = plot(fpr, data_shuf', 'Color', [shuf_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % 3. Plot Mean lines (Thick)
    PLOT(1) = plot(fpr, TPR_real_mean, 'Color', real_col, 'LineWidth', 2.5);
    PLOT(2) = plot(fpr, TPR_shuf_mean, 'Color', shuf_col, 'LineWidth', 2.5);
    
    % Formatting
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    title(sprintf('ROC curve V1 %i ms bins', timebins(n)));
    
    legend(PLOT(1:2), {'Real Mean', 'Shuffle Mean'}, 'Location', 'southeast', 'Box', 'off');
    
    xlim([0 1]); ylim([0 1]);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end



for n = 1:2
    nexttile
    % Extract data for the current condition
    real_data = KDE_RUN_validations_all.V1_AUC(n,:);
    shuffled_data = KDE_RUN_validations_all.V1_AUC_shuffled(n,:);
    
    % Combine for easier mean calculation
    % AUC = [real_data; shuffled_data];
    mean_real = mean(real_data, 'omitnan');
    mean_shuf = mean(shuffled_data, 'omitnan');
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos = [1 2];
    
    hold on

    % 1. Plot the Bars (Mean only, no error bars)
    plot([x_pos(1)-line_w, x_pos(1)+line_w], [mean_real, mean_real], 'Color', bar_colors(1,:), 'LineWidth', 4);
    plot([x_pos(2)-line_w, x_pos(2)+line_w], [mean_shuf, mean_shuf], 'Color', bar_colors(2,:), 'LineWidth', 4);

    % 2. Plot Raw Data Points and Connections with Jitter
    rng(1); % Seed for consistent jitter appearance
    jitter_strength = 0.2; % Adjusted for better visibility of distributions
    
    for s = 1:length(real_data)
        if ~isnan(real_data(s)) && ~isnan(shuffled_data(s))
            % Calculate jittered x-positions
            % This spreads the dots randomly around the x_pos center
            x1 = x_pos(1) + (rand - 0.5) * jitter_strength;
            x2 = x_pos(2) + (rand - 0.5) * jitter_strength;
            
            % Plot the connection line first (so it's behind the dots)
            plot([x1, x2], [real_data(s), shuffled_data(s)], 'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.5);
            
            % Plot the individual dots
            scatter(x1, real_data(s), 50, bar_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.3);
            scatter(x2, shuffled_data(s), 50, bar_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.3);
        end
    end

    % Statistics and Formatting
    [p, h] = signrank(real_data, shuffled_data, 'tail', 'right');
    
    xlim([-0.5 3.5])
    xticks([1 2]);
    xticklabels({'Real', 'Shuffled'});
    ylabel('AUC');
    ylim([0 1.1]);
    yticks([0:0.25:1])
    title(['Condition ' num2str(n) ' (p=' num2str(p, '%.4e') ')']);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end

% for n = 1:2
%     nexttile
% 
%     AUC = [KDE_RUN_validations_all.V1_AUC(n,:);KDE_RUN_validations_all.V1_AUC_shuffled(n,:)];
%     bar_colors = [231,41,138; 0, 0, 0]/255;
%     x_pos=[1 2];
%     for i = 1:2
% 
%         mean_AUC = mean(AUC(i,:),'omitnan');
%         auc_errors = [std(AUC(i,:),'omitnan')./sqrt(size(AUC(i,:),2)),...  % lower error
%             std(AUC(i,:),'omitnan')./sqrt(size(AUC(i,:),2))];   % upper error
%         hold on
%         bar(x_pos(i), mean_AUC, 0.4, 'FaceAlpha',0.5, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
%         errorbar(x_pos(i), mean_AUC, auc_errors(1), auc_errors(2), 'k', 'linestyle', 'none', 'linewidth', 1);
%     end
% 
% 
%     [p,h,stats] = signrank(AUC(1,:),AUC(2,:),'tail','right');
% 
%     p
% 
%     xticks([1 2]);
%     xticklabels({'Real', 'Shuffled'});
%     ylabel('AUC');
%     ylim([0 1]);
%     title('Mean AUC – Real vs Shuffled');
%     set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% end


save_all_figures(fullfile(analysis_folder,'V1-HPC spatial representation'),[],'ContentType','vector')



%% Add on PLS KDE regression with log_odds and z-scored values
% -- Initialization and setup
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS = {'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS, option);
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
% experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

session_count = 0;

% KDE_reactivation_V1_all = struct();
KDE_reactivation_V1_UP_all = struct();
KDE_reactivation_V1_DOWN_all = struct();

KDE_reactivation_UP_all = struct();
KDE_reactivation_DOWN_all = struct();
% KDE_reactivation_all = struct();

for nsession = 1:length(experiment_info)
    tic
    fprintf('session %i\n', nsession);
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName, Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName, Stimulus_type));
    if isempty(stimulus_name), continue; end

    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT}, option);
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));

    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH, '..', 'best_channels.mat'));

    if length(stimulus_name) > 1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            n = length(stimulus_name) > 1 && find(contains(stimulus_name,'_2')) || 1;
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
    options = session_info(n).probe(1);

    if isempty(dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'))), continue; end

    if exist(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'), 'file')
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        session_clusters_RUN = session_clusters;
        clear session_clusters
    end
    if exist(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'), 'file')
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        session_clusters_RUN = session_clusters;
        clear session_clusters
    end

    % --- V1 data ---
    if contains(stimulus_name{n},'Sleep')
        % load(fullfile(options.ANALYSIS_DATAPATH,'PLS_KDE_reactivation_V1.mat'),'KDE_reactivation_V1');
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_V1_DOWN.mat'),'KDE_reactivation_V1_DOWN');
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_V1_UP.mat'),'KDE_reactivation_V1_UP');
    end

    for nprobe = 1:length(KDE_reactivation_V1)
        structures = {
            % 'KDE_reactivation_V1','KDE_reactivation_V1_all';
            'KDE_reactivation_V1_UP','KDE_reactivation_V1_UP_all';
            'KDE_reactivation_V1_DOWN','KDE_reactivation_V1_DOWN_all'
        };

        all_log_odds = [
            % log(KDE_reactivation_V1(nprobe).event_T1_probability ./ KDE_reactivation_V1(nprobe).event_T2_probability);
            log(KDE_reactivation_V1_UP(nprobe).event_T1_probability ./ KDE_reactivation_V1_UP(nprobe).event_T2_probability);
            log(KDE_reactivation_V1_DOWN(nprobe).event_T1_probability ./ KDE_reactivation_V1_DOWN(nprobe).event_T2_probability)
        ];
        all_log_odds = all_log_odds((isfinite(all_log_odds)));

        all_bias = [
            % KDE_reactivation_V1(nprobe).event_bias;
            KDE_reactivation_V1_UP(nprobe).event_bias;
            KDE_reactivation_V1_DOWN(nprobe).event_bias
        ];

        for s = 1:size(structures,1)
            src = eval(structures{s,1});
            dst = structures{s,2};

            real_bias = src(nprobe).event_bias(:)';
            % shuffled = src(nprobe).event_T1_probability_shuffled ./ (src(nprobe).event_T1_probability_shuffled + src(nprobe).event_T2_probability_shuffled);
            % nbin = size(shuffled, 2);
            % zscored_bias_shuffled = nan(1, nbin);
            % percentile = nan(1, nbin);
            % for i = 1:nbin
            %     valid = shuffled(:, i);
            %     valid = valid(~isnan(valid));
            %     if ~isempty(valid) && ~isnan(real_bias(i))
            %         percentile(i) = sum(valid < real_bias(i), 'omitnan') / sum(~isnan(valid)) * 100;
            %         zscored_bias_shuffled(i) = (real_bias(i) - mean(valid, 'omitnan')) / std(valid, 'omitnan');
            %     end
            % end

            log_odds = log(src(nprobe).event_T1_probability ./ src(nprobe).event_T2_probability);
            % log_odds_shuffled = log(src(nprobe).event_T1_probability_shuffled ./ src(nprobe).event_T2_probability_shuffled);
            % zscored_log_odds_shuffled = nan(1, size(log_odds_shuffled, 2));
            % log_odds_percentile = nan(1, size(log_odds_shuffled, 2));
            % for i = 1:size(log_odds_shuffled,2)
            %     valid = log_odds_shuffled(:, i);
            %     valid = valid(~isnan(valid));
            %     if ~isempty(valid) && ~isnan(log_odds(i))
            %         zscored_log_odds_shuffled(i) = (log_odds(i) - mean(valid, 'omitnan')) / std(valid, 'omitnan');
            %         log_odds_percentile(i) = sum(valid < log_odds(i), 'omitnan') / sum(~isnan(valid)) * 100;
            %     end
            % end

            eval([dst '(nprobe).bias{nsession} = real_bias;']);
            eval([dst '(nprobe).zscored_bias{nsession} = (real_bias - mean(all_bias,''omitnan'')) ./ std(all_bias,''omitnan'');']);
            % eval([dst '(nprobe).zscored_bias_shuffled{nsession} = zscored_bias_shuffled;']);
            % eval([dst '(nprobe).percentile{nsession} = percentile;']);

            eval([dst '(nprobe).log_odds{nsession} = log_odds;']);
            eval([dst '(nprobe).zscored_log_odds{nsession} = (log_odds - mean(all_log_odds,''omitnan'')) ./ std(all_log_odds,''omitnan'');']);
            % eval([dst '(nprobe).zscored_log_odds_shuffled{nsession} = zscored_log_odds_shuffled;']);
            % eval([dst '(nprobe).log_odds_percentile{nsession} = log_odds_percentile;']);

            % New: Save event metadata
            eval([dst '(nprobe).event_bins{nsession} = src(nprobe).event_bins;']);
            eval([dst '(nprobe).event_id{nsession} = src(nprobe).event_id;']);
        end
    end
    clear KDE_reactivation_V1 KDE_reactivation_V1_DOWN KDE_reactivation_V1_UP

    % --- General (non-V1) data ---
    % load(fullfile(options.ANALYSIS_DATAPATH,'PLS_KDE_reactivation.mat'),'KDE_reactivation');
    load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_DOWN.mat'),'KDE_reactivation_DOWN');
    load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_UP.mat'),'KDE_reactivation_UP');

    for nprobe = 1:length(KDE_reactivation)
        structures = {
            % 'KDE_reactivation','KDE_reactivation_all';
            'KDE_reactivation_UP','KDE_reactivation_UP_all';
            'KDE_reactivation_DOWN','KDE_reactivation_DOWN_all'
        };

        all_log_odds = [
            % log(KDE_reactivation(nprobe).event_T1_probability ./ KDE_reactivation(nprobe).event_T2_probability);
            log(KDE_reactivation_UP(nprobe).event_T1_probability ./ KDE_reactivation_UP(nprobe).event_T2_probability);
            log(KDE_reactivation_DOWN(nprobe).event_T1_probability ./ KDE_reactivation_DOWN(nprobe).event_T2_probability)
        ];
        all_log_odds = all_log_odds((isfinite(all_log_odds)));

        all_bias = [
            % KDE_reactivation(nprobe).event_bias;
            KDE_reactivation_UP(nprobe).event_bias;
            KDE_reactivation_DOWN(nprobe).event_bias
        ];

        for s = 1:size(structures,1)
            src = eval(structures{s,1});
            dst = structures{s,2};

            real_bias = src(nprobe).event_bias(:)';
            % shuffled = src(nprobe).event_T1_probability_shuffled ./ (src(nprobe).event_T1_probability_shuffled + src(nprobe).event_T2_probability_shuffled);
            % nbin = size(shuffled, 2);
            % zscored_bias_shuffled = nan(1, nbin);
            % percentile = nan(1, nbin);
            % for i = 1:nbin
            %     valid = shuffled(:, i);
            %     valid = valid(~isnan(valid));
            %     if ~isempty(valid) && ~isnan(real_bias(i))
            %         percentile(i) = sum(valid < real_bias(i), 'omitnan') / sum(~isnan(valid)) * 100;
            %         zscored_bias_shuffled(i) = (real_bias(i) - mean(valid, 'omitnan')) / std(valid, 'omitnan');
            %     end
            % end

            log_odds = log(src(nprobe).event_T1_probability ./ src(nprobe).event_T2_probability);
            % log_odds_shuffled = log(src(nprobe).event_T1_probability_shuffled ./ src(nprobe).event_T2_probability_shuffled);
            % zscored_log_odds_shuffled = nan(1, size(log_odds_shuffled, 2));
            % log_odds_percentile = nan(1, size(log_odds_shuffled, 2));
            % for i = 1:size(log_odds_shuffled,2)
            %     valid = log_odds_shuffled(:, i);
            %     valid = valid(~isnan(valid));
            %     if ~isempty(valid) && ~isnan(log_odds(i))
            %         zscored_log_odds_shuffled(i) = (log_odds(i) - mean(valid, 'omitnan')) / std(valid, 'omitnan');
            %         log_odds_percentile(i) = sum(valid < log_odds(i), 'omitnan') / sum(~isnan(valid)) * 100;
            %     end
            % end

            eval([dst '(nprobe).bias{nsession} = real_bias;']);
            eval([dst '(nprobe).zscored_bias{nsession} = (real_bias - mean(all_bias,''omitnan'')) ./ std(all_bias,''omitnan'');']);
            % eval([dst '(nprobe).zscored_bias_shuffled{nsession} = zscored_bias_shuffled;']);
            % eval([dst '(nprobe).percentile{nsession} = percentile;']);

            eval([dst '(nprobe).log_odds{nsession} = log_odds;']);
            eval([dst '(nprobe).zscored_log_odds{nsession} = (log_odds - mean(all_log_odds,''omitnan'')) ./ std(all_log_odds,''omitnan'');']);
            % eval([dst '(nprobe).zscored_log_odds_shuffled{nsession} = zscored_log_odds_shuffled;']);
            % eval([dst '(nprobe).log_odds_percentile{nsession} = log_odds_percentile;']);

            % New: Save event time metadata
            eval([dst '(nprobe).event_bins{nsession} = src(nprobe).event_bins;']);
            eval([dst '(nprobe).event_id{nsession} = src(nprobe).event_id;']);
        end
    end
    clear KDE_reactivation KDE_reactivation_DOWN KDE_reactivation_UP
    toc
end

% Save all results
if exist('D:\\corticohippocampal_replay','dir')
    analysis_folder = 'D:\\corticohippocampal_replay';
elseif exist('P:\\corticohippocampal_replay','dir')
    analysis_folder = 'P:\\corticohippocampal_replay';
end

save(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
save(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
% save(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
save(fullfile(analysis_folder,'KDE_reactivation_V1_DOWN_all_POST.mat'),'KDE_reactivation_V1_DOWN_all')
save(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
% save(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')

%% Extract all UP state log odds

SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';
load(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;

KDE_log_odds_UP = struct();
for nprobe = 1:2
    for nsession =1:length(experiment_info)
        event_id = find(slow_waves_all(nprobe).UP_session_count == nsession);
        temp_event_id = KDE_reactivation_UP_all(nprobe).event_id{nsession};
        temp_log_odds = KDE_reactivation_UP_all(nprobe).zscored_log_odds{nsession};
        temp1 = temp_log_odds(isfinite(temp_log_odds));
        temp_log_odds(temp_log_odds>=inf) = prctile(temp1,99.5);
        temp_log_odds(temp_log_odds<=-inf) = prctile(temp1,0.5);

        temp_log_odds_V1 = KDE_reactivation_V1_UP_all(nprobe).zscored_log_odds{nsession};
        temp1 = temp_log_odds_V1(isfinite(temp_log_odds_V1));
        temp_log_odds_V1(temp_log_odds_V1>=inf) = prctile(temp1,99.5);
        temp_log_odds_V1(temp_log_odds_V1<=-inf) = prctile(temp1,0.5);

        for nevent = 1:max(KDE_reactivation_UP_all(nprobe).event_id{nsession})
%             event_duration = diff(slow_waves_all(1).UP_ints(event_id(nevent),:));
%             temp = temp_log_odds(temp_event_id == nevent);
%             nbins = floor(event_duration/0.01); 
            KDE_log_odds_UP(nprobe).HPC{event_id(nevent)} = temp_log_odds(temp_event_id == nevent);
            KDE_log_odds_UP(nprobe).V1{event_id(nevent)} = temp_log_odds_V1(temp_event_id == nevent);
        end
    end
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_log_odds_UP.mat'),'KDE_log_odds_UP');

%% Add on PLS KDE regression Ripples

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

session_count = 0;

KDE_reactivation_V1_ripples_all= struct();
KDE_reactivation_ripples_all = struct();

for nsession =1:length(experiment_info)

    tic
    disp(sprintf('session %i',nsession))
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    if length(stimulus_name)>1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2'));
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

    if contains(stimulus_name{n},'Sleep')
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_V1_ripples.mat'),'KDE_reactivation_V1_ripples');
        load(fullfile(options.ANALYSIS_DATAPATH,'KDE_reactivation_ripples.mat'),'KDE_reactivation_ripples');
    end

    for nprobe = 1:length(KDE_reactivation_V1_ripples)
        real_bias = KDE_reactivation_V1_ripples(nprobe).event_bias(:)';
        bias_distribution = [real_bias];

        log_odds = log(KDE_reactivation_V1_ripples(nprobe).event_T1_probability ./ KDE_reactivation_V1_ripples(nprobe).event_T2_probability)';
        log_odds_distribution = log_odds((isfinite(log_odds)));

        KDE_reactivation_V1_ripples_all(nprobe).event_bins{nsession} = KDE_reactivation_V1_ripples(nprobe).event_bins;
        KDE_reactivation_V1_ripples_all(nprobe).event_id{nsession} = KDE_reactivation_V1_ripples(nprobe).event_id;

        KDE_reactivation_V1_ripples_all(nprobe).bias{nsession} = real_bias;
        KDE_reactivation_V1_ripples_all(nprobe).zscored_bias{nsession} = (real_bias - mean(bias_distribution,'omitnan')) ./ std(bias_distribution,'omitnan');

        KDE_reactivation_V1_ripples_all(nprobe).log_odds{nsession} = log_odds;
        KDE_reactivation_V1_ripples_all(nprobe).zscored_log_odds{nsession} = (log_odds - mean(log_odds_distribution,'omitnan')) ./ std(log_odds_distribution,'omitnan');
    end

    for nprobe = 1:length(KDE_reactivation_ripples)
        real_bias = KDE_reactivation_ripples(nprobe).event_bias(:)';
        bias_distribution = [real_bias];

        log_odds = log(KDE_reactivation_ripples(nprobe).event_T1_probability ./ KDE_reactivation_ripples(nprobe).event_T2_probability)';
        log_odds_distribution = log_odds((isfinite(log_odds)));

        KDE_reactivation_ripples_all(nprobe).event_bins{nsession} = KDE_reactivation_ripples(nprobe).event_bins;
        KDE_reactivation_ripples_all(nprobe).event_id{nsession} = KDE_reactivation_ripples(nprobe).event_id;

        KDE_reactivation_ripples_all(nprobe).bias{nsession} = real_bias;
        KDE_reactivation_ripples_all(nprobe).zscored_bias{nsession} = (real_bias - mean(bias_distribution,'omitnan')) ./ std(bias_distribution,'omitnan');

        KDE_reactivation_ripples_all(nprobe).log_odds{nsession} = log_odds;
        KDE_reactivation_ripples_all(nprobe).zscored_log_odds{nsession} = (log_odds - mean(log_odds_distribution,'omitnan')) ./ std(log_odds_distribution,'omitnan');
    end
    clear KDE_reactivation_ripples KDE_reactivation_V1_ripples
    toc
end

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

save(fullfile(analysis_folder,'KDE_reactivation_ripples_all_POST.mat'),'KDE_reactivation_ripples_all')
save(fullfile(analysis_folder,'KDE_reactivation_V1_ripples_all_POST.mat'),'KDE_reactivation_V1_ripples_all')


%%
%%%%%
%%%%%
%%%%%
%%%%%
%% Process and extract KDE reactivation PSTH

clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
% load(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_DOWN_all_POST.mat'),'KDE_reactivation_V1_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
% load(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')


%%%% Combined V1 and HPC PSTH extraction for both bias and log odds

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

% Parameters
psth_window = [-1 1];
psth_step = 0.01;
nbins = round(diff(psth_window) / psth_step);
half_bins = nbins / 2;

psth_window = [-1 1];
psth_step = 0.01;
% bins = psth_window(1)+psth_step/2 : psth_step: psth_window(end)-psth_step/2; 
% Initialize PSTH structure
KDE_reactivation_PSTH = struct();

for region = ["V1", "HPC"]
    for nprobe = 1:2
        up_bias_all = [];
        down_bias_all = [];
        up_bias_all_z = [];
        down_bias_all_z = [];
        % up_bias_all_zshuff = [];
        % down_bias_all_zshuff = [];

        up_logodds_all = [];
        down_logodds_all = [];
        up_logodds_all_z = [];
        down_logodds_all_z = [];
        % up_logodds_all_zshuff = [];
        % down_logodds_all_zshuff = [];

        if region == "V1"
            up_struct = KDE_reactivation_V1_UP_all;
            down_struct = KDE_reactivation_V1_DOWN_all;
        else
            up_struct = KDE_reactivation_UP_all;
            down_struct = KDE_reactivation_DOWN_all;
        end

        for nsession = 1:length(up_struct(nprobe).bias)
            bias_up = up_struct(nprobe).bias{nsession};
            bias_down = down_struct(nprobe).bias{nsession};
            z_up = up_struct(nprobe).zscored_bias{nsession};
            z_down = down_struct(nprobe).zscored_bias{nsession};
            % zshuff_up = up_struct(nprobe).zscored_bias_shuffled{nsession};
            % zshuff_down = down_struct(nprobe).zscored_bias_shuffled{nsession};

            logodds_up = up_struct(nprobe).log_odds{nsession};
            logodds_down = down_struct(nprobe).log_odds{nsession};
            zlog_up = up_struct(nprobe).zscored_log_odds{nsession};
            zlog_down = down_struct(nprobe).zscored_log_odds{nsession};
            % zshuff_log_up = up_struct(nprobe).zscored_log_odds_shuffled{nsession};
            % zshuff_log_down = down_struct(nprobe).zscored_log_odds_shuffled{nsession};

            time_up = up_struct(nprobe).event_bins{nsession}(:,1)';
            time_down = down_struct(nprobe).event_bins{nsession}(:,1)';
            event_id_up = up_struct(nprobe).event_id{nsession};
            event_id_down = down_struct(nprobe).event_id{nsession};

            up_all = find(slow_waves_all(nprobe).UP_session_count == nsession);
            %             UP_event_index = intersect(up_all, probability(nprobe).UP_all_index);
            UP_event_index = up_all;

            temp_index = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
            [~, ia, ib] = intersect(slow_waves_all(nprobe).DOWN_ints(temp_index, 2), ...
                slow_waves_all(nprobe).UP_ints(UP_event_index, 1));

            up_ints = slow_waves_all(nprobe).UP_ints(up_all, :);
            previous_DOWN_event_index = ia;
            UP_event_index = find(ismember(up_all, UP_event_index));


            down_all = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
            %             DOWN_event_index = intersect(down_all, probability(nprobe).DOWN_all_index);
            %             DOWN_event_index = find(ismember(down_all, DOWN_event_index));
            DOWN_event_index = down_all;
            DOWN_event_index = find(ismember(down_all, DOWN_event_index));
            down_ints = slow_waves_all(nprobe).DOWN_ints(down_all, :);

            for i = 1:length(UP_event_index)
                up_bias = nan(1, nbins);
                up_z = nan(1, nbins);
                % up_zshuff = nan(1, nbins);

                up_log = nan(1, nbins);
                up_log_z = nan(1, nbins);
                % up_log_zshuff = nan(1, nbins);

                if i <= length(previous_DOWN_event_index)
                    eid = previous_DOWN_event_index(i);
                    event_end = down_ints(eid, 2);
                    rel_time = time_down - event_end;
                    idx = find(event_id_down == eid & rel_time >= -1 & rel_time < 0);
                    rel_times = rel_time(idx);
                    bin_idx = round((rel_times + 1) / psth_step) + 1;
                    valid = bin_idx > 0 & bin_idx <= half_bins;

                    up_bias(bin_idx(valid)) = bias_down(idx(valid));
                    up_z(bin_idx(valid)) = z_down(idx(valid));
                    % up_zshuff(bin_idx(valid)) = zshuff_down(idx(valid));

                    up_log(bin_idx(valid)) = logodds_down(idx(valid));
                    up_log_z(bin_idx(valid)) = zlog_down(idx(valid));
                    % up_log_zshuff(bin_idx(valid)) = zshuff_log_down(idx(valid));
                end

                eid = UP_event_index(i);
                event_start = up_ints(eid, 1);
                rel_time = time_up - event_start;
                idx = find(event_id_up == eid & rel_time >= 0 & rel_time < 1);
                rel_times = rel_time(idx);
                bin_idx = round(rel_times / psth_step) + half_bins;
                valid = bin_idx > half_bins & bin_idx <= nbins;

                up_bias(bin_idx(valid)) = bias_up(idx(valid));
                up_z(bin_idx(valid)) = z_up(idx(valid));
                % up_zshuff(bin_idx(valid)) = zshuff_up(idx(valid));

                up_log(bin_idx(valid)) = logodds_up(idx(valid));
                up_log_z(bin_idx(valid)) = zlog_up(idx(valid));
                % up_log_zshuff(bin_idx(valid)) = zshuff_log_up(idx(valid));

                up_bias_all = [up_bias_all; up_bias];
                up_bias_all_z = [up_bias_all_z; up_z];
                % up_bias_all_zshuff = [up_bias_all_zshuff; up_zshuff];

                up_logodds_all = [up_logodds_all; up_log];
                up_logodds_all_z = [up_logodds_all_z; up_log_z];
                % up_logodds_all_zshuff = [up_logodds_all_zshuff; up_log_zshuff];
            end

            for i = 1:length(DOWN_event_index)
                down_bias = nan(1, nbins);
                down_z = nan(1, nbins);
                down_zshuff = nan(1, nbins);

                down_log = nan(1, nbins);
                down_log_z = nan(1, nbins);
                % down_log_zshuff = nan(1, nbins);

                if i <= length(UP_event_index)
                    eid = UP_event_index(i);
                    event_end = up_ints(eid, 2);
                    rel_time = time_up - event_end;
                    idx = find(event_id_up == eid & rel_time >= -1 & rel_time < 0);
                    rel_times = rel_time(idx);
                    bin_idx = round((rel_times + 1) / psth_step) + 1;
                    valid = bin_idx > 0 & bin_idx <= half_bins;

                    down_bias(bin_idx(valid)) = bias_up(idx(valid));
                    down_z(bin_idx(valid)) = z_up(idx(valid));
                    % down_zshuff(bin_idx(valid)) = zshuff_up(idx(valid));

                    down_log(bin_idx(valid)) = logodds_up(idx(valid));
                    down_log_z(bin_idx(valid)) = zlog_up(idx(valid));
                    % down_log_zshuff(bin_idx(valid)) = zshuff_log_up(idx(valid));
                end

                eid = DOWN_event_index(i);
                event_start = down_ints(eid, 1);
                rel_time = time_down - event_start;
                idx = find(event_id_down == eid & rel_time >= 0 & rel_time < 1);
                rel_times = rel_time(idx);
                bin_idx = round(rel_times / psth_step) + half_bins;
                valid = bin_idx > half_bins & bin_idx <= nbins;

                down_bias(bin_idx(valid)) = bias_down(idx(valid));
                down_z(bin_idx(valid)) = z_down(idx(valid));
                % down_zshuff(bin_idx(valid)) = zshuff_down(idx(valid));

                down_log(bin_idx(valid)) = logodds_down(idx(valid));
                down_log_z(bin_idx(valid)) = zlog_down(idx(valid));
                % down_log_zshuff(bin_idx(valid)) = zshuff_log_down(idx(valid));

                down_bias_all = [down_bias_all; down_bias];
                down_bias_all_z = [down_bias_all_z; down_z];
                % down_bias_all_zshuff = [down_bias_all_zshuff; down_zshuff];

                down_logodds_all = [down_logodds_all; down_log];
                down_logodds_all_z = [down_logodds_all_z; down_log_z];
                % down_logodds_all_zshuff = [down_logodds_all_zshuff; down_log_zshuff];
            end
        end

        prefix = sprintf('%s_', region);
        KDE_reactivation_PSTH(nprobe).([prefix 'UP']) = up_bias_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN']) = down_bias_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'UP_z']) = up_bias_all_z;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_z']) = down_bias_all_z;
        % KDE_reactivation_PSTH(nprobe).([prefix 'UP_zshuff']) = up_bias_all_zshuff;
        % KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_zshuff']) = down_bias_all_zshuff;

        KDE_reactivation_PSTH(nprobe).([prefix 'UP_log_odds']) = up_logodds_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_log_odds']) = down_logodds_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'UP_log_odds_z']) = up_logodds_all_z;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_log_odds_z']) = down_logodds_all_z;
        % KDE_reactivation_PSTH(nprobe).([prefix 'UP_log_odds_zshuff']) = up_logodds_all_zshuff;
        % KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_log_odds_zshuff']) = down_logodds_all_zshuff;
    end
end
for nprobe = 1:2
    KDE_reactivation_PSTH(nprobe).tvec = psth_window(1)+psth_step/2 : psth_step: psth_window(end)-psth_step/2;
end
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_PSTH.mat'),'KDE_reactivation_PSTH','-v7.3');

% 
% %%%%% HPC ripples
% % Parameters
% ripple_window = [0, 0.2];  % seconds
% psth_step = 0.01;          % 10 ms
% nbins = round(diff(ripple_window) / psth_step);
% 
% % Metrics to extract
% % metric_names = {'bias', 'zscored_bias', 'zscored_bias_shuffled'};
% metric_names = {'bias', 'zscored_bias', 'log_odds', 'zscored_log_odds','zscored_log_odds_shuffled','log_odds_percentile'};
% % Output containers
% ripple_bias_masked_SWS = struct();
% ripple_bias_masked_nonSWS = struct();
% for region = ["HPC", "V1"]
%     for metric = metric_names
%         ripple_bias_masked_SWS.(region).(metric{1})     = cell(1, 2);
%         ripple_bias_masked_nonSWS.(region).(metric{1})  = cell(1, 2);
%     end
% end
% 
% % Main loop
% for isSWS = [1, 0]
%     for nprobe = 1:2
%         for region = ["HPC", "V1"]
% 
%             if contains(region,'V1')
%                 temp = KDE_reactivation_V1_all;
%             else
%                 temp = KDE_reactivation_all;
%             end
% 
%             for metric = metric_names
%                 metric_str = metric{1};
%                 masked_all = [];
% 
%                 for nsession = 1:length(temp(nprobe).(metric_str))
%                     % Get all ripple indices for this session
%                     session_mask = ripples_all(nprobe).session_count == nsession;
%                     ripple_idx_all = find(session_mask);
% 
%                     % Get SWS or non-SWS ripple subset
%                     state_mask = ripples_all(nprobe).SWS_index == isSWS;
%                     [ripple_idx_all,ripple_idx_session,~] = intersect(ripple_idx_all,find(state_mask));
% 
%                     if isempty(ripple_idx_session), continue; end
% 
%                     ripple_times = [ ...
%                         ripples_all(nprobe).onset(ripple_idx_all), ...
%                         ripples_all(nprobe).offset(ripple_idx_all)];
% 
%                     % Get data
%                     bias_vals = temp(nprobe).(metric_str){nsession};
%                     time_bins = temp(nprobe).event_bins{nsession}(:,1)';
%                     event_ids = temp(nprobe).event_id{nsession};
% 
%                     for i = 1:length(ripple_idx_session)
%                         ripple_id = ripple_idx_session(i);
%                         t_start = ripple_times(i, 1);
% 
%                         rel_time = time_bins - t_start;
%                         idx = find(event_ids == ripple_id & rel_time >= 0 & rel_time <= 0.2);
% 
%                         ripple_bias = nan(1, nbins);
% 
%                         ripple_bias(1:length(idx)) = bias_vals(idx);
%                         % Mask if another ripple starts within this window
%                         all_starts = ripples_all(nprobe).onset(ripple_idx_all);
%                         next_ripple = all_starts > t_start & all_starts < (t_start + 0.2);
%                         if any(next_ripple)
%                             t_next = min(all_starts(next_ripple));
%                             mask_start = round((t_next - t_start) / psth_step) + 1;
%                             ripple_bias(mask_start:end) = nan;
%                         end
% 
%                         masked_all = [masked_all; ripple_bias];
%                     end
%                 end
% 
%                 % Save result
%                 if isSWS
%                     ripple_bias_masked_SWS.(region).(metric_str){nprobe} = masked_all;
%                 else
%                     ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = masked_all;
%                 end
%             end
%         end
%     end
% end
% 
% 
% % SWS ripples
% % SWS ripples
% % bias
% KDE_reactivation_content.HPC_ripples                  = [ripple_bias_masked_SWS.HPC.bias{1}; ripple_bias_masked_SWS.HPC.bias{2}];
% KDE_reactivation_content.V1_ripples                   = [ripple_bias_masked_SWS.V1.bias{1}; ripple_bias_masked_SWS.V1.bias{2}];
% KDE_reactivation_content.HPC_awake_ripples            = [ripple_bias_masked_nonSWS.HPC.bias{1}; ripple_bias_masked_nonSWS.HPC.bias{2}];
% KDE_reactivation_content.V1_awake_ripples             = [ripple_bias_masked_nonSWS.V1.bias{1}; ripple_bias_masked_nonSWS.V1.bias{2}];
% 
% % zscored_bias
% KDE_reactivation_content.HPC_z_ripples                = [ripple_bias_masked_SWS.HPC.zscored_bias{1}; ripple_bias_masked_SWS.HPC.zscored_bias{2}];
% KDE_reactivation_content.V1_z_ripples                 = [ripple_bias_masked_SWS.V1.zscored_bias{1}; ripple_bias_masked_SWS.V1.zscored_bias{2}];
% KDE_reactivation_content.HPC_z_awake_ripples          = [ripple_bias_masked_nonSWS.HPC.zscored_bias{1}; ripple_bias_masked_nonSWS.HPC.zscored_bias{2}];
% KDE_reactivation_content.V1_z_awake_ripples           = [ripple_bias_masked_nonSWS.V1.zscored_bias{1}; ripple_bias_masked_nonSWS.V1.zscored_bias{2}];
% 
% % log_odds
% KDE_reactivation_content.HPC_logodds_ripples         = [ripple_bias_masked_SWS.HPC.log_odds{1}; ripple_bias_masked_SWS.HPC.log_odds{2}];
% KDE_reactivation_content.V1_logodds_ripples          = [ripple_bias_masked_SWS.V1.log_odds{1}; ripple_bias_masked_SWS.V1.log_odds{2}];
% KDE_reactivation_content.HPC_logodds_awake_ripples   = [ripple_bias_masked_nonSWS.HPC.log_odds{1}; ripple_bias_masked_nonSWS.HPC.log_odds{2}];
% KDE_reactivation_content.V1_logodds_awake_ripples    = [ripple_bias_masked_nonSWS.V1.log_odds{1}; ripple_bias_masked_nonSWS.V1.log_odds{2}];
% 
% % zscored_log_odds
% KDE_reactivation_content.HPC_z_logodds_ripples       = [ripple_bias_masked_SWS.HPC.zscored_log_odds{1}; ripple_bias_masked_SWS.HPC.zscored_log_odds{2}];
% KDE_reactivation_content.V1_z_logodds_ripples        = [ripple_bias_masked_SWS.V1.zscored_log_odds{1}; ripple_bias_masked_SWS.V1.zscored_log_odds{2}];
% KDE_reactivation_content.HPC_z_logodds_awake_ripples = [ripple_bias_masked_nonSWS.HPC.zscored_log_odds{1}; ripple_bias_masked_nonSWS.HPC.zscored_log_odds{2}];
% KDE_reactivation_content.V1_z_logodds_awake_ripples  = [ripple_bias_masked_nonSWS.V1.zscored_log_odds{1}; ripple_bias_masked_nonSWS.V1.zscored_log_odds{2}];
% 
% % zscored_log_odds_shuffled
% KDE_reactivation_content.HPC_zshuff_logodds_ripples       = [ripple_bias_masked_SWS.HPC.zscored_log_odds_shuffled{1}; ripple_bias_masked_SWS.HPC.zscored_log_odds_shuffled{2}];
% KDE_reactivation_content.V1_zshuff_logodds_ripples        = [ripple_bias_masked_SWS.V1.zscored_log_odds_shuffled{1}; ripple_bias_masked_SWS.V1.zscored_log_odds_shuffled{2}];
% KDE_reactivation_content.HPC_zshuff_logodds_awake_ripples = [ripple_bias_masked_nonSWS.HPC.zscored_log_odds_shuffled{1}; ripple_bias_masked_nonSWS.HPC.zscored_log_odds_shuffled{2}];
% KDE_reactivation_content.V1_zshuff_logodds_awake_ripples  = [ripple_bias_masked_nonSWS.V1.zscored_log_odds_shuffled{1}; ripple_bias_masked_nonSWS.V1.zscored_log_odds_shuffled{2}];
% 
% % log_odds_percentile
% KDE_reactivation_content.HPC_logodds_percentile_ripples         = [ripple_bias_masked_SWS.HPC.log_odds_percentile{1}; ripple_bias_masked_SWS.HPC.log_odds_percentile{2}];
% KDE_reactivation_content.V1_logodds_percentile_ripples          = [ripple_bias_masked_SWS.V1.log_odds_percentile{1}; ripple_bias_masked_SWS.V1.log_odds_percentile{2}];
% KDE_reactivation_content.HPC_logodds_percentile_awake_ripples   = [ripple_bias_masked_nonSWS.HPC.log_odds_percentile{1}; ripple_bias_masked_nonSWS.HPC.log_odds_percentile{2}];
% KDE_reactivation_content.V1_logodds_percentile_awake_ripples    = [ripple_bias_masked_nonSWS.V1.log_odds_percentile{1}; ripple_bias_masked_nonSWS.V1.log_odds_percentile{2}];
% 
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'),'KDE_reactivation_content');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
load(fullfile(analysis_folder,'KDE_reactivation_ripples_all_POST.mat'),'KDE_reactivation_ripples_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_ripples_all_POST.mat'),'KDE_reactivation_V1_ripples_all')

ripple_window = [-1, 1];
psth_step = 0.01;
nbins = round(diff(ripple_window) / psth_step);
metric_names = {'bias', 'zscored_bias', 'log_odds', 'zscored_log_odds'};
regions = ["HPC", "V1"];

% Initialize
ripple_bias_masked_SWS = struct();
ripple_bias_masked_nonSWS = struct();
KDE_reactivation_ripples_PSTH = struct();
nan_mask = cell(1,2);          % for probe 1,2 SWS
nan_mask_awake = cell(1,2);    % for probe 1,2 awake

% Loop over states
for isSWS = [1, 0]
    for nprobe = 1:2
        for region = regions
            for metric = metric_names
                metric_str = metric{1};
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = [];
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = [];
                end
            end
        end

        nan_mask_probe = [];  % accumulate here

        for region = regions
            if region == "V1"
                temp = KDE_reactivation_V1_ripples_all;
            else
                temp = KDE_reactivation_ripples_all;
            end

            for metric = metric_names
                metric_str = metric{1};
                metric_all = [];

                for nsession = 1:length(temp(nprobe).(metric_str))
                    session_mask = ripples_all(nprobe).session_count == nsession;
                    ripple_idx_all = find(session_mask);
                    state_mask = ripples_all(nprobe).SWS_index == isSWS;
                    [ripple_idx_all, ripple_idx_session] = intersect(ripple_idx_all, find(state_mask));
                    if isempty(ripple_idx_session), continue; end

                    ripple_times = [ripples_all(nprobe).onset(ripple_idx_all), ...
                                    ripples_all(nprobe).offset(ripple_idx_all)];
                    metric_vals = temp(nprobe).(metric_str){nsession};
                    time_bins = temp(nprobe).event_bins{nsession}(:,1)';
                    event_ids = temp(nprobe).event_id{nsession};

                    for i = 1:length(ripple_idx_session)
                        ripple_id = ripple_idx_session(i);
                        t_start = ripple_times(i, 1);
                        rel_time = time_bins - t_start;
                        idx = find(event_ids == ripple_id & rel_time >= -1 & rel_time < 1);

                        ripple_vec = nan(1, nbins);
                        ripple_vec(1:length(idx)) = metric_vals(idx);
                        metric_all = [metric_all; ripple_vec];

                        % Only compute mask once (on HPC region, one metric)
                        if strcmp(region, 'HPC') && strcmp(metric_str, 'bias')
                            nan_row = zeros(1, nbins);  % 0 = keep, NaN = mask

                            all_onsets = ripple_times(:,1);
                            all_offsets = ripple_times(:,2);

                            % Mask next ripple
                            next_ripples = all_onsets > t_start & all_onsets < (t_start + 1);
                            if any(next_ripples)
                                t_next = min(all_onsets(next_ripples));
                                mask_start = round((t_next - t_start) / psth_step) + 101;
                                if mask_start <= nbins
                                    nan_row(mask_start:end) = NaN;
                                end
                            end

                            % Mask previous ripple
                            prev_ripples = all_offsets < t_start & all_offsets > (t_start - 1);
                            if any(prev_ripples)
                                t_prev = max(all_offsets(prev_ripples));
                                mask_end = round((t_prev - t_start) / psth_step) + 101;
                                if mask_end >= 1
                                    nan_row(1:mask_end) = NaN;
                                end
                            end

                            nan_mask_probe = [nan_mask_probe; nan_row];
                        end
                    end
                end

                % Save metric data per region/state
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = metric_all;
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = metric_all;
                end
            end
        end

        % Store mask for this probe
        if isSWS
            nan_mask{nprobe} = nan_mask_probe;
        else
            nan_mask_awake{nprobe} = nan_mask_probe;
        end
    end
end

% Combine and store
for region = regions
    region = char(region);
    KDE_reactivation_ripples_PSTH.([region '_ripples'])             = [ripple_bias_masked_SWS.(region).bias{1};              ripple_bias_masked_SWS.(region).bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_ripples'])           = [ripple_bias_masked_SWS.(region).zscored_bias{1};      ripple_bias_masked_SWS.(region).zscored_bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_logodds_ripples'])     = [ripple_bias_masked_SWS.(region).log_odds{1};          ripple_bias_masked_SWS.(region).log_odds{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_logodds_ripples'])   = [ripple_bias_masked_SWS.(region).zscored_log_odds{1};  ripple_bias_masked_SWS.(region).zscored_log_odds{2}];

    KDE_reactivation_ripples_PSTH.([region '_awake_ripples'])             = [ripple_bias_masked_nonSWS.(region).bias{1};              ripple_bias_masked_nonSWS.(region).bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_awake_ripples'])           = [ripple_bias_masked_nonSWS.(region).zscored_bias{1};      ripple_bias_masked_nonSWS.(region).zscored_bias{2}];
    KDE_reactivation_ripples_PSTH.([region '_logodds_awake_ripples'])     = [ripple_bias_masked_nonSWS.(region).log_odds{1};          ripple_bias_masked_nonSWS.(region).log_odds{2}];
    KDE_reactivation_ripples_PSTH.([region '_z_logodds_awake_ripples'])   = [ripple_bias_masked_nonSWS.(region).zscored_log_odds{1};  ripple_bias_masked_nonSWS.(region).zscored_log_odds{2}];
end

% Add nan masks
KDE_reactivation_ripples_PSTH.nan_mask        = [nan_mask{1}; nan_mask{2}];
KDE_reactivation_ripples_PSTH.nan_mask_awake  = [nan_mask_awake{1}; nan_mask_awake{2}];

% Save
save(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');


%% Add on RRR decoded log odds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -- Initialization and setup
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS = {'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS, option);
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
% experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

session_count = 0;

RRR_reactivation_V1_ripples_all = struct();
RRR_reactivation_V1_UP_all = struct();
RRR_reactivation_V1_DOWN_all = struct();

RRR_reactivation_UP_all = struct();
RRR_reactivation_DOWN_all = struct();
RRR_reactivation_ripples_all = struct();

for nsession = 1:length(experiment_info)
    tic
    fprintf('session %i\n', nsession);
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName, Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName, Stimulus_type));
    if isempty(stimulus_name), continue; end

    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT}, option);
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));

    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH, '..', 'best_channels.mat'));

    if length(stimulus_name) > 1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            n = length(stimulus_name) > 1 && find(contains(stimulus_name,'_2')) || 1;
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
    options = session_info(n).probe(1);

    if isempty(dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'))), continue; end

    if exist(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'), 'file')
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        session_clusters_RUN = session_clusters;
        clear session_clusters
    end
    if exist(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'), 'file')
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        session_clusters_RUN = session_clusters;
        clear session_clusters
    end

    % --- V1 data ---
    if contains(stimulus_name{n},'Sleep')
        load(fullfile(options.ANALYSIS_DATAPATH,'RRR_reactivation_V1_ripples.mat'),'RRR_reactivation_V1_ripples');
        load(fullfile(options.ANALYSIS_DATAPATH,'RRR_reactivation_V1_DOWN.mat'),'RRR_reactivation_V1_DOWN');
        load(fullfile(options.ANALYSIS_DATAPATH,'RRR_reactivation_V1_UP.mat'),'RRR_reactivation_V1_UP');
    end

    for nprobe = 1:length(RRR_reactivation_V1_ripples)
        structures = {
            'RRR_reactivation_V1_ripples','RRR_reactivation_V1_ripples_all';
            'RRR_reactivation_V1_UP','RRR_reactivation_V1_UP_all';
            'RRR_reactivation_V1_DOWN','RRR_reactivation_V1_DOWN_all'
        };

        all_log_odds = [
            % log(KDE_reactivation_V1(nprobe).event_T1_probability ./ KDE_reactivation_V1(nprobe).event_T2_probability);
            RRR_reactivation_V1_UP(nprobe).log_odds;
            RRR_reactivation_V1_DOWN(nprobe).log_odds
        ];
        all_log_odds = all_log_odds((isfinite(all_log_odds)));

        all_log_odds_ripples = [
            % log(KDE_reactivation(nprobe).event_T1_probability ./ KDE_reactivation(nprobe).event_T2_probability);
            RRR_reactivation_V1_ripples(nprobe).log_odds;
            ];
        all_log_odds_ripples = all_log_odds_ripples((isfinite(all_log_odds_ripples)));

        for s = 1:size(structures,1)
            src = eval(structures{s,1});
            dst = structures{s,2};

            log_odds = src(nprobe).log_odds;
            % log_odds_shuffled = log(src(nprobe).event_T1_probability_shuffled ./ src(nprobe).event_T2_probability_shuffled);
            % zscored_log_odds_shuffled = nan(1, size(log_odds_shuffled, 2));
            % log_odds_percentile = nan(1, size(log_odds_shuffled, 2));
            % for i = 1:size(log_odds_shuffled,2)
            %     valid = log_odds_shuffled(:, i);
            %     valid = valid(~isnan(valid));
            %     if ~isempty(valid) && ~isnan(log_odds(i))
            %         zscored_log_odds_shuffled(i) = (log_odds(i) - mean(valid, 'omitnan')) / std(valid, 'omitnan');
            %         log_odds_percentile(i) = sum(valid < log_odds(i), 'omitnan') / sum(~isnan(valid)) * 100;
            %     end
            % end

            eval([dst '(nprobe).log_odds{nsession} = log_odds;']);

            if contains(structures{s,1},'ripple')
                eval([dst '(nprobe).zscored_log_odds{nsession} = (log_odds - mean(all_log_odds_ripples,''omitnan'')) ./ std(all_log_odds_ripples,''omitnan'');']);
            else
                eval([dst '(nprobe).zscored_log_odds{nsession} = (log_odds - mean(all_log_odds,''omitnan'')) ./ std(all_log_odds,''omitnan'');']);
            end
            % eval([dst '(nprobe).zscored_log_odds_shuffled{nsession} = zscored_log_odds_shuffled;']);
            % eval([dst '(nprobe).log_odds_percentile{nsession} = log_odds_percentile;']);

            % New: Save event metadata
            eval([dst '(nprobe).event_bins{nsession} = src(nprobe).event_bins;']);
            eval([dst '(nprobe).event_id{nsession} = src(nprobe).event_id;']);
        end
    end
    clear KDE_reactivation_V1 KDE_reactivation_V1_DOWN KDE_reactivation_V1_UP

    % --- General (non-V1) data ---
    load(fullfile(options.ANALYSIS_DATAPATH,'RRR_reactivation_ripples.mat'),'RRR_reactivation_ripples');
    load(fullfile(options.ANALYSIS_DATAPATH,'RRR_reactivation_DOWN.mat'),'RRR_reactivation_DOWN');
    load(fullfile(options.ANALYSIS_DATAPATH,'RRR_reactivation_UP.mat'),'RRR_reactivation_UP');

    for nprobe = 1:length(RRR_reactivation_ripples)
        structures = {
            'RRR_reactivation_ripples','RRR_reactivation_ripples_all';
            'RRR_reactivation_UP','RRR_reactivation_UP_all';
            'RRR_reactivation_DOWN','RRR_reactivation_DOWN_all'
        };

        all_log_odds = [
            % log(KDE_reactivation(nprobe).event_T1_probability ./ KDE_reactivation(nprobe).event_T2_probability);
            RRR_reactivation_UP(nprobe).log_odds;
            RRR_reactivation_DOWN(nprobe).log_odds
        ];

        all_log_odds_ripples = [
            % log(KDE_reactivation(nprobe).event_T1_probability ./ KDE_reactivation(nprobe).event_T2_probability);
            RRR_reactivation_ripples(nprobe).log_odds;
            ];

        all_log_odds = all_log_odds((isfinite(all_log_odds)));
        all_log_odds_ripples = all_log_odds_ripples((isfinite(all_log_odds_ripples)));

        for s = 1:size(structures,1)
            src = eval(structures{s,1});
            dst = structures{s,2};

            log_odds = src(nprobe).log_odds;
            % log_odds_shuffled = log(src(nprobe).event_T1_probability_shuffled ./ src(nprobe).event_T2_probability_shuffled);
            % zscored_log_odds_shuffled = nan(1, size(log_odds_shuffled, 2));
            % log_odds_percentile = nan(1, size(log_odds_shuffled, 2));
            % for i = 1:size(log_odds_shuffled,2)
            %     valid = log_odds_shuffled(:, i);
            %     valid = valid(~isnan(valid));
            %     if ~isempty(valid) && ~isnan(log_odds(i))
            %         zscored_log_odds_shuffled(i) = (log_odds(i) - mean(valid, 'omitnan')) / std(valid, 'omitnan');
            %         log_odds_percentile(i) = sum(valid < log_odds(i), 'omitnan') / sum(~isnan(valid)) * 100;
            %     end
            % end

            eval([dst '(nprobe).log_odds{nsession} = log_odds;']);
            if contains(structures{s,1},'ripple')
                eval([dst '(nprobe).zscored_log_odds{nsession} = (log_odds - mean(all_log_odds_ripples,''omitnan'')) ./ std(all_log_odds_ripples,''omitnan'');']);
            else
                eval([dst '(nprobe).zscored_log_odds{nsession} = (log_odds - mean(all_log_odds,''omitnan'')) ./ std(all_log_odds,''omitnan'');']);
            end

            % eval([dst '(nprobe).zscored_log_odds_shuffled{nsession} = zscored_log_odds_shuffled;']);
            % eval([dst '(nprobe).log_odds_percentile{nsession} = log_odds_percentile;']);

            % New: Save event time metadata
            eval([dst '(nprobe).event_bins{nsession} = src(nprobe).event_bins;']);
            eval([dst '(nprobe).event_id{nsession} = src(nprobe).event_id;']);
        end
    end
    clear RRR_reactivation RRR_reactivation_DOWN RRR_reactivation_UP
    toc
end
% KDE_reactivation_V1
% Save all results
if exist('D:\\corticohippocampal_replay','dir')
    analysis_folder = 'D:\\corticohippocampal_replay';
elseif exist('P:\\corticohippocampal_replay','dir')
    analysis_folder = 'P:\\corticohippocampal_replay';
end

save(fullfile(analysis_folder,'RRR_reactivation_DOWN_all_POST.mat'),'RRR_reactivation_DOWN_all')
save(fullfile(analysis_folder,'RRR_reactivation_UP_all_POST.mat'),'RRR_reactivation_UP_all')
% save(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
save(fullfile(analysis_folder,'RRR_reactivation_V1_DOWN_all_POST.mat'),'RRR_reactivation_V1_DOWN_all')
save(fullfile(analysis_folder,'RRR_reactivation_V1_UP_all_POST.mat'),'RRR_reactivation_V1_UP_all')
% save(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')
save(fullfile(analysis_folder,'RRR_reactivation_ripples_all_POST.mat'),'RRR_reactivation_ripples_all')
save(fullfile(analysis_folder,'RRR_reactivation_V1_ripples_all_POST.mat'),'RRR_reactivation_V1_ripples_all')

%%%%%
%% Process and extract RRR reactivation PSTH

clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'RRR_reactivation_DOWN_all_POST.mat'),'RRR_reactivation_DOWN_all')
load(fullfile(analysis_folder,'RRR_reactivation_UP_all_POST.mat'),'RRR_reactivation_UP_all')
% load(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
load(fullfile(analysis_folder,'RRR_reactivation_V1_DOWN_all_POST.mat'),'RRR_reactivation_V1_DOWN_all')
load(fullfile(analysis_folder,'RRR_reactivation_V1_UP_all_POST.mat'),'RRR_reactivation_V1_UP_all')
% load(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')

% KDE_reactivation_V1
%%%% Combined V1 and HPC PSTH extraction for both bias and log odds

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

% Parameters
psth_window = [-1 1];
psth_step = 0.01;
nbins = round(diff(psth_window) / psth_step);
half_bins = nbins / 2;

% Initialize PSTH structure
KDE_reactivation_PSTH = struct();

for region = ["V1", "HPC"]
    for nprobe = 1:2
        up_bias_all = [];
        down_bias_all = [];
        up_bias_all_z = [];
        down_bias_all_z = [];
        % up_bias_all_zshuff = [];
        % down_bias_all_zshuff = [];

        up_logodds_all = [];
        down_logodds_all = [];
        up_logodds_all_z = [];
        down_logodds_all_z = [];
        % up_logodds_all_zshuff = [];
        % down_logodds_all_zshuff = [];

        if region == "V1"
            up_struct = KDE_reactivation_V1_UP_all;
            down_struct = KDE_reactivation_V1_DOWN_all;
        else
            up_struct = KDE_reactivation_UP_all;
            down_struct = KDE_reactivation_DOWN_all;
        end

        for nsession = 1:length(up_struct(nprobe).bias)
            bias_up = up_struct(nprobe).bias{nsession};
            bias_down = down_struct(nprobe).bias{nsession};
            z_up = up_struct(nprobe).zscored_bias{nsession};
            z_down = down_struct(nprobe).zscored_bias{nsession};
            % zshuff_up = up_struct(nprobe).zscored_bias_shuffled{nsession};
            % zshuff_down = down_struct(nprobe).zscored_bias_shuffled{nsession};

            logodds_up = up_struct(nprobe).log_odds{nsession};
            logodds_down = down_struct(nprobe).log_odds{nsession};
            zlog_up = up_struct(nprobe).zscored_log_odds{nsession};
            zlog_down = down_struct(nprobe).zscored_log_odds{nsession};
            % zshuff_log_up = up_struct(nprobe).zscored_log_odds_shuffled{nsession};
            % zshuff_log_down = down_struct(nprobe).zscored_log_odds_shuffled{nsession};

            time_up = up_struct(nprobe).event_bins{nsession}(:,1)';
            time_down = down_struct(nprobe).event_bins{nsession}(:,1)';
            event_id_up = up_struct(nprobe).event_id{nsession};
            event_id_down = down_struct(nprobe).event_id{nsession};

            up_all = find(slow_waves_all(nprobe).UP_session_count == nsession);
            UP_event_index = intersect(up_all, probability(nprobe).UP_all_index);

            temp_index = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
            [~, ia, ib] = intersect(slow_waves_all(nprobe).DOWN_ints(temp_index, 2), ...
                                    slow_waves_all(nprobe).UP_ints(UP_event_index, 1));

            up_ints = slow_waves_all(nprobe).UP_ints(up_all, :);
            previous_DOWN_event_index = ia;
            UP_event_index = find(ismember(up_all, UP_event_index));

            down_all = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
            DOWN_event_index = intersect(down_all, probability(nprobe).DOWN_all_index);
            DOWN_event_index = find(ismember(down_all, DOWN_event_index));
            down_ints = slow_waves_all(nprobe).DOWN_ints(down_all, :);

            for i = 1:length(UP_event_index)
                up_bias = nan(1, nbins);
                up_z = nan(1, nbins);
                % up_zshuff = nan(1, nbins);

                up_log = nan(1, nbins);
                up_log_z = nan(1, nbins);
                % up_log_zshuff = nan(1, nbins);

                if i <= length(previous_DOWN_event_index)
                    eid = previous_DOWN_event_index(i);
                    event_end = down_ints(eid, 2);
                    rel_time = time_down - event_end;
                    idx = find(event_id_down == eid & rel_time >= -1 & rel_time < 0);
                    rel_times = rel_time(idx);
                    bin_idx = round((rel_times + 1) / psth_step) + 1;
                    valid = bin_idx > 0 & bin_idx <= half_bins;

                    up_bias(bin_idx(valid)) = bias_down(idx(valid));
                    up_z(bin_idx(valid)) = z_down(idx(valid));
                    % up_zshuff(bin_idx(valid)) = zshuff_down(idx(valid));

                    up_log(bin_idx(valid)) = logodds_down(idx(valid));
                    up_log_z(bin_idx(valid)) = zlog_down(idx(valid));
                    % up_log_zshuff(bin_idx(valid)) = zshuff_log_down(idx(valid));
                end

                eid = UP_event_index(i);
                event_start = up_ints(eid, 1);
                rel_time = time_up - event_start;
                idx = find(event_id_up == eid & rel_time >= 0 & rel_time < 1);
                rel_times = rel_time(idx);
                bin_idx = round(rel_times / psth_step) + half_bins;
                valid = bin_idx > half_bins & bin_idx <= nbins;

                up_bias(bin_idx(valid)) = bias_up(idx(valid));
                up_z(bin_idx(valid)) = z_up(idx(valid));
                % up_zshuff(bin_idx(valid)) = zshuff_up(idx(valid));

                up_log(bin_idx(valid)) = logodds_up(idx(valid));
                up_log_z(bin_idx(valid)) = zlog_up(idx(valid));
                % up_log_zshuff(bin_idx(valid)) = zshuff_log_up(idx(valid));

                up_bias_all = [up_bias_all; up_bias];
                up_bias_all_z = [up_bias_all_z; up_z];
                % up_bias_all_zshuff = [up_bias_all_zshuff; up_zshuff];

                up_logodds_all = [up_logodds_all; up_log];
                up_logodds_all_z = [up_logodds_all_z; up_log_z];
                % up_logodds_all_zshuff = [up_logodds_all_zshuff; up_log_zshuff];
            end

            for i = 1:length(DOWN_event_index)
                down_bias = nan(1, nbins);
                down_z = nan(1, nbins);
                down_zshuff = nan(1, nbins);

                down_log = nan(1, nbins);
                down_log_z = nan(1, nbins);
                % down_log_zshuff = nan(1, nbins);

                if i <= length(UP_event_index)
                    eid = UP_event_index(i);
                    event_end = up_ints(eid, 2);
                    rel_time = time_up - event_end;
                    idx = find(event_id_up == eid & rel_time >= -1 & rel_time < 0);
                    rel_times = rel_time(idx);
                    bin_idx = round((rel_times + 1) / psth_step) + 1;
                    valid = bin_idx > 0 & bin_idx <= half_bins;

                    down_bias(bin_idx(valid)) = bias_up(idx(valid));
                    down_z(bin_idx(valid)) = z_up(idx(valid));
                    % down_zshuff(bin_idx(valid)) = zshuff_up(idx(valid));

                    down_log(bin_idx(valid)) = logodds_up(idx(valid));
                    down_log_z(bin_idx(valid)) = zlog_up(idx(valid));
                    % down_log_zshuff(bin_idx(valid)) = zshuff_log_up(idx(valid));
                end

                eid = DOWN_event_index(i);
                event_start = down_ints(eid, 1);
                rel_time = time_down - event_start;
                idx = find(event_id_down == eid & rel_time >= 0 & rel_time < 1);
                rel_times = rel_time(idx);
                bin_idx = round(rel_times / psth_step) + half_bins;
                valid = bin_idx > half_bins & bin_idx <= nbins;

                down_bias(bin_idx(valid)) = bias_down(idx(valid));
                down_z(bin_idx(valid)) = z_down(idx(valid));
                % down_zshuff(bin_idx(valid)) = zshuff_down(idx(valid));

                down_log(bin_idx(valid)) = logodds_down(idx(valid));
                down_log_z(bin_idx(valid)) = zlog_down(idx(valid));
                % down_log_zshuff(bin_idx(valid)) = zshuff_log_down(idx(valid));

                down_bias_all = [down_bias_all; down_bias];
                down_bias_all_z = [down_bias_all_z; down_z];
                % down_bias_all_zshuff = [down_bias_all_zshuff; down_zshuff];

                down_logodds_all = [down_logodds_all; down_log];
                down_logodds_all_z = [down_logodds_all_z; down_log_z];
                % down_logodds_all_zshuff = [down_logodds_all_zshuff; down_log_zshuff];
            end
        end

        prefix = sprintf('%s_', region);
        KDE_reactivation_PSTH(nprobe).([prefix 'UP']) = up_bias_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN']) = down_bias_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'UP_z']) = up_bias_all_z;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_z']) = down_bias_all_z;
        % KDE_reactivation_PSTH(nprobe).([prefix 'UP_zshuff']) = up_bias_all_zshuff;
        % KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_zshuff']) = down_bias_all_zshuff;

        KDE_reactivation_PSTH(nprobe).([prefix 'UP_log_odds']) = up_logodds_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_log_odds']) = down_logodds_all;
        KDE_reactivation_PSTH(nprobe).([prefix 'UP_log_odds_z']) = up_logodds_all_z;
        KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_log_odds_z']) = down_logodds_all_z;
        % KDE_reactivation_PSTH(nprobe).([prefix 'UP_log_odds_zshuff']) = up_logodds_all_zshuff;
        % KDE_reactivation_PSTH(nprobe).([prefix 'DOWN_log_odds_zshuff']) = down_logodds_all_zshuff;
    end
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','RRR_reactivation_PSTH.mat'),'RRR_reactivation_PSTH','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
% load(fullfile(analysis_folder,'KDE_reactivation_ripples_all_POST.mat'),'KDE_reactivation_ripples_all')
load(fullfile(analysis_folder,'RRR_reactivation_ripples_all_POST.mat'),'RRR_reactivation_ripples_all')
load(fullfile(analysis_folder,'RRR_reactivation_V1_ripples_all_POST.mat'),'RRR_reactivation_V1_ripples_all')
load(fullfile(analysis_folder,'ripples_all_POST.mat'))

ripple_window = [-1, 1];
psth_step = 0.01;
nbins = round(diff(ripple_window) / psth_step);
metric_names = {'log_odds', 'zscored_log_odds'};
regions = ["HPC", "V1"];

% Initialize
ripple_bias_masked_SWS = struct();
ripple_bias_masked_nonSWS = struct();
KDE_reactivation_ripples_PSTH = struct();
nan_mask = cell(1,2);          % for probe 1,2 SWS
nan_mask_awake = cell(1,2);    % for probe 1,2 awake

% Loop over states
for isSWS = [1, 0]
    for nprobe = 1:2
        for region = regions
            for metric = metric_names
                metric_str = metric{1};
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = [];
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = [];
                end
            end
        end

        nan_mask_probe = [];  % accumulate here

        for region = regions
            if region == "V1"
                temp = RRR_reactivation_V1_ripples_all;
            else
                temp = RRR_reactivation_ripples_all;
            end

            for metric = metric_names
                metric_str = metric{1};
                metric_all = [];

                for nsession = 1:length(temp(nprobe).(metric_str))
                    session_mask = ripples_all(nprobe).session_count == nsession;
                    ripple_idx_all = find(session_mask);
                    state_mask = ripples_all(nprobe).SWS_index == isSWS;
                    [ripple_idx_all, ripple_idx_session] = intersect(ripple_idx_all, find(state_mask));
                    if isempty(ripple_idx_session), continue; end

                    ripple_times = [ripples_all(nprobe).onset(ripple_idx_all), ...
                                    ripples_all(nprobe).offset(ripple_idx_all)];
                    metric_vals = temp(nprobe).(metric_str){nsession};
                    time_bins = temp(nprobe).event_bins{nsession}(:,1)';
                    event_ids = temp(nprobe).event_id{nsession};

                    for i = 1:length(ripple_idx_session)
                        ripple_id = ripple_idx_session(i);
                        t_start = ripple_times(i, 1);
                        rel_time = time_bins - t_start;
                        idx = find(event_ids == ripple_id & rel_time >= -1 & rel_time < 1);

                        ripple_vec = nan(1, nbins);
                        ripple_vec(1:length(idx)) = metric_vals(idx);
                        metric_all = [metric_all; ripple_vec];

                        % Only compute mask once (on HPC region, one metric)
                        if strcmp(region, 'HPC') && strcmp(metric_str, 'log_odds')
                            nan_row = zeros(1, nbins);  % 0 = keep, NaN = mask

                            all_onsets = ripple_times(:,1);
                            all_offsets = ripple_times(:,2);

                            % Mask next ripple
                            next_ripples = all_onsets > t_start & all_onsets < (t_start + 1);
                            if any(next_ripples)
                                t_next = min(all_onsets(next_ripples));
                                mask_start = round((t_next - t_start) / psth_step) + 101;
                                if mask_start <= nbins
                                    nan_row(mask_start:end) = NaN;
                                end
                            end

                            % Mask previous ripple
                            prev_ripples = all_offsets < t_start & all_offsets > (t_start - 1);
                            if any(prev_ripples)
                                t_prev = max(all_offsets(prev_ripples));
                                mask_end = round((t_prev - t_start) / psth_step) + 101;
                                if mask_end >= 1
                                    nan_row(1:mask_end) = NaN;
                                end
                            end

                            nan_mask_probe = [nan_mask_probe; nan_row];
                        end
                    end
                end

                % Save metric data per region/state
                if isSWS
                    ripple_bias_masked_SWS.(region).(metric_str){nprobe} = metric_all;
                else
                    ripple_bias_masked_nonSWS.(region).(metric_str){nprobe} = metric_all;
                end
            end
        end

        % Store mask for this probe
        if isSWS
            nan_mask{nprobe} = nan_mask_probe;
        else
            nan_mask_awake{nprobe} = nan_mask_probe;
        end
    end
end

% Combine and store
for region = regions
    region = char(region);
%     RRR_reactivation_ripples_PSTH.([region '_ripples'])             = [ripple_bias_masked_SWS.(region).bias{1};              ripple_bias_masked_SWS.(region).bias{2}];
%     RRR_reactivation_ripples_PSTH.([region '_z_ripples'])           = [ripple_bias_masked_SWS.(region).zscored_bias{1};      ripple_bias_masked_SWS.(region).zscored_bias{2}];
    RRR_reactivation_ripples_PSTH.([region '_logodds_ripples'])     = [ripple_bias_masked_SWS.(region).log_odds{1};          ripple_bias_masked_SWS.(region).log_odds{2}];
    RRR_reactivation_ripples_PSTH.([region '_z_logodds_ripples'])   = [ripple_bias_masked_SWS.(region).zscored_log_odds{1};  ripple_bias_masked_SWS.(region).zscored_log_odds{2}];

%     KDE_reactivation_ripples_PSTH.([region '_awake_ripples'])             = [ripple_bias_masked_nonSWS.(region).bias{1};              ripple_bias_masked_nonSWS.(region).bias{2}];
%     KDE_reactivation_ripples_PSTH.([region '_z_awake_ripples'])           = [ripple_bias_masked_nonSWS.(region).zscored_bias{1};      ripple_bias_masked_nonSWS.(region).zscored_bias{2}];
    RRR_reactivation_ripples_PSTH.([region '_logodds_awake_ripples'])     = [ripple_bias_masked_nonSWS.(region).log_odds{1};          ripple_bias_masked_nonSWS.(region).log_odds{2}];
    RRR_reactivation_ripples_PSTH.([region '_z_logodds_awake_ripples'])   = [ripple_bias_masked_nonSWS.(region).zscored_log_odds{1};  ripple_bias_masked_nonSWS.(region).zscored_log_odds{2}];
end

% Add nan masks
RRR_reactivation_ripples_PSTH.nan_mask        = [nan_mask{1}; nan_mask{2}];
RRR_reactivation_ripples_PSTH.nan_mask_awake  = [nan_mask_awake{1}; nan_mask_awake{2}];

% Save
save(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'RRR_reactivation_ripples_PSTH.mat'), 'RRR_reactivation_ripples_PSTH');
% load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');

