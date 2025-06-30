%% Context-selective ripple modulation in V1 and HPC
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
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

% load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
% load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))
% load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
load(fullfile(analysis_folder,'session_clusters_all_POST.mat'))

%     load(fullfile(analysis_folder,'session_clusters_all_POST.mat'),'session_clusters_all','-v7.3')
sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);

%% Ripple modulation in V1 and HPC
ripple_modulation_PSTH_all = [];
psthBinSize = 0.01;
windows = [-1 1];
for nsession = 1:length(sessions_to_process)
    %     ripple_modulation_PSTH_all{nsession} = [];

    all_clusters = session_clusters_all.spatial_cell_id{nsession};
    %       plot(unique(session_clusters_all.spike_id{nsession}));hold on;plot(all_clusters)

    %         for nprobe = 1:length(ripples_all)
    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];
    tic
    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession},session_clusters_all.spike_id{nsession},windows,psthBinSize,...
        'unit_id',all_clusters,'event_times',event_times,'event_id',event_id,'saving_PSTH',0,'shuffle_option',1);
    ripple_modulation_PSTH_all{nsession} = ripple_modulation;
    toc
end

save(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all','-v7.3')
%  = struct();


%% context-selective ripple modulation + Spatial correlation and ripple correlation in V1 and HPC 
load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')
% load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
%  = struct();

%%%% KDE bias
timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;
bins_to_use = bin_centers>0 & bin_centers<0.1;


session_id = [ripples_all(1).session_count(ripples_all(1).SWS_index); ripples_all(2).session_count(ripples_all(2).SWS_index)];

z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);
clear z_bias1


psthBinSize = 0.01;
windows = [-1 1];

z_V1_population_ripple_PSTH{1} = [];
z_V1_population_ripple_PSTH{2} = [];

context_modulation_all = struct();
% context_corr_all = struct();

for nsession = 1:length(sessions_to_process)
    all_clusters = session_clusters_all.spatial_cell_id{nsession};

    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession},session_clusters_all.spike_id{nsession},windows,psthBinSize,...
        'unit_id',all_clusters,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);
    % bins = ripple_modulation.bins > 0 & ripple_modulation.bins < 0.5;
   
    % all_clusters = session_clusters_all.spatial_cell_id{nsession};
    all_regions = session_clusters_all.region{nsession};
    V1_id = find(contains(all_regions,'V1'));
    HPC_id = find(contains(all_regions,'HPC'));

    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));

    %%%%%%% Populational zscored firing rate difference
    T1_V1_cell =  find(session_clusters_all.mean_FR{nsession}(V1_id,1) >  session_clusters_all.mean_FR{nsession}(V1_id,2));
    T2_V1_cell =  find(session_clusters_all.mean_FR{nsession}(V1_id,1) <  session_clusters_all.mean_FR{nsession}(V1_id,2));
    % T1_HPC_cell =  session_clusters_all.mean_FR{nsession}(HPC_id,1) >  session_clusters_all.mean_FR{nsession}(HPC_id,2);
    % T2_HPC_cell =  session_clusters_all.mean_FR{nsession}(HPC_id,1) >  session_clusters_all.mean_FR{nsession}(HPC_id,2);
    z_V1_population_ripple_PSTH{1} = [z_V1_population_ripple_PSTH{1}; squeeze(mean(ripple_modulation.PSTH_zscored(T1_V1_cell,:,:),1,'omitnan'))];
    z_V1_population_ripple_PSTH{2} = [z_V1_population_ripple_PSTH{2}; squeeze(mean(ripple_modulation.PSTH_zscored(T2_V1_cell,:,:),1,'omitnan'))];
    % z_V1_population_ripple_PSTH{2} = [z_V1_population_ripple_PSTH{2}; squeeze(mean(ripple_modulation.PSTH_zscored(T2_V1_cell,:,:),1,'omitnan'))];
    % mean(z_V1_population_ripple_PSTH{1}(T1_index,:)-z_V1_population_ripple_PSTH{2}(T1_index,:),'omitnan') -  mean(z_V1_population_ripple_PSTH{1}(T2_index,:)-z_V1_population_ripple_PSTH{2}(T2_index,:),'omitnan')


    %%%%%%% Context selective ripple modulation
    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH{nsession}(2,nCell,:))';

        % Pre- and post-ripple firing rate modulation
        context_modulation_all.PRE_ripple_FR{nsession}(1,nCell) = ...
            mean(context_modulation_all.PSTH_diff{nsession}(nCell, ripple_modulation.bins > -0.2 & ripple_modulation.bins < 0));

        context_modulation_all.POST_ripple_FR{nsession}(1,nCell) = ...
            mean(context_modulation_all.PSTH_diff{nsession}(nCell, ripple_modulation.bins > 0 & ripple_modulation.bins < 0.2));

        % Firing rates on both tracks
        context_modulation_all.FR_track{nsession}(1,nCell) = session_clusters_all.mean_FR{nsession}(nCell,1);
        context_modulation_all.FR_track{nsession}(2,nCell) = session_clusters_all.mean_FR{nsession}(nCell,2);

        % Z-scored track FR (relative to combined distribution)
        FR_distribution = reshape([session_clusters_all.spatial_response_raw{nsession}{nCell,1}; ...
            session_clusters_all.spatial_response_raw{nsession}{nCell,2}], 1, []);
        context_modulation_all.z_FR_track{nsession}(1,nCell) = ...
            (session_clusters_all.mean_FR{nsession}(nCell,1) - mean(FR_distribution)) ./ std(FR_distribution);

        context_modulation_all.z_FR_track{nsession}(2,nCell) = ...
            (session_clusters_all.mean_FR{nsession}(nCell,2) - mean(FR_distribution)) ./ std(FR_distribution);

        context_modulation_all.region{nsession}(nCell) = all_regions(nCell);

        context_modulation_all.session_id{nsession}(nCell) = nsession;
    end

end



    % subplot(2,2,1)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.POST_ripple_FR(V1_id))
    % subplot(2,2,2)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
    % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial corr vs ripple corr
context_corr_all = struct();

for nsession = 1:length(sessions_to_process)
    all_clusters = session_clusters_all.spatial_cell_id{nsession};

    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    % ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession},session_clusters_all.spike_id{nsession},windows,psthBinSize,...
    %     'unit_id',all_clusters,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);
    % bins = ripple_modulation.bins > 0 & ripple_modulation.bins < 0.5;
   
    % all_clusters = session_clusters_all.spatial_cell_id{nsession};
    all_regions = session_clusters_all.region{nsession};
    V1_id = find(contains(all_regions,'V1'));
    HPC_id = find(contains(all_regions,'HPC'));

    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));
    %%%%%%%% Spatial Corr context diff and Ripple corr context diff
    V1_cell_spatial_psth=[];
    V1_cell_psth=[];

    HPC_cell_spatial_psth=[];
    HPC_cell_psth=[];

    % Ensure parallel pool is open
    if isempty(gcp('nocreate'))
        parpool;
    end


    tic
    [pair_corr, pair_pval,lags] = compute_ripple_xcorr_cell_pair(...
        session_clusters_all.spike_times{nsession}, ...
        session_clusters_all.spike_id{nsession}, ...
        all_clusters(V1_id), all_clusters(HPC_id), event_times,'lag_range',[0],'window',0.2,'step', 0.2,'shuffle_option',0);

    context_corr_all.ripple_corr_all{nsession} = pair_corr;
    context_corr_all.ripple_corr_pval_all{nsession} = pair_pval;

    [pair_corr_T1, pair_pval,lags] = compute_ripple_xcorr_cell_pair(...
        session_clusters_all.spike_times{nsession}, ...
        session_clusters_all.spike_id{nsession}, ...
        all_clusters(V1_id), all_clusters(HPC_id), event_times(T1_index),'lag_range',[0],'window',0.2,'step', 0.2,'shuffle_option',0);

    context_corr_all.ripple_corr_T1{nsession} = pair_corr_T1;
    context_corr_all.ripple_corr_pval_T1{nsession} = pair_pval;

    [pair_corr_T2,pair_pval,lags] = compute_ripple_xcorr_cell_pair(...
        session_clusters_all.spike_times{nsession}, ...
        session_clusters_all.spike_id{nsession}, ...
        all_clusters(V1_id), all_clusters(HPC_id), event_times(T2_index),'lag_range',[0],'window',0.2,'step', 0.2,'shuffle_option',0);

    context_corr_all.ripple_corr_T2{nsession} = pair_corr_T2;
    context_corr_all.ripple_corr_pval_T2{nsession} = pair_pval;
    toc

    pair_pval_T1 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
    pair_corr_T1 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
    pair_pval_T2 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
    pair_corr_T2 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)

    spatial_raw  = session_clusters_all.spatial_response_raw{nsession};

    HPC_spatial = cell(length(HPC_id), 2);
    for mCell = 1:length(HPC_id)
        HPC_spatial{mCell,1} = reshape(spatial_raw {HPC_id(mCell), 1}, 1, []);
        HPC_spatial{mCell,2} = reshape(spatial_raw {HPC_id(mCell), 2}, 1, []);
    end

    for nCell = 1:length(V1_id)
        V1_spatial{1} = reshape(spatial_raw {V1_id(nCell), 1}, 1, []);
        V1_spatial{2} = reshape(spatial_raw {V1_id(nCell), 2}, 1, []);

        % === Preload HPC cell responses for this loop ===


        parfor mCell = 1:length(HPC_id)
            [r1, p1] = corr(HPC_spatial{mCell,1}', V1_spatial{1}', 'type', 'Spearman', 'Rows', 'complete');
            [r2, p2] = corr(HPC_spatial{mCell,2}', V1_spatial{2}', 'type', 'Spearman', 'Rows', 'complete');

            pair_pval_T1(mCell, nCell) = p1;
            pair_corr_T1(mCell, nCell) = r1;
            pair_pval_T2(mCell, nCell) = p2;
            pair_corr_T2(mCell, nCell) = r2;
        end
    end


    % pair_pval = reshape(pair_pval, nHPC, nV1); % Reshape to 2D [nHPC x nV1]

    context_corr_all.spatial_corr_T1{nsession}= pair_corr_T1;
    context_corr_all.spatial_corr_pval_T1{nsession} = pair_pval_T1;

    context_corr_all.spatial_corr_T2{nsession} = pair_corr_T2;
    context_corr_all.spatial_corr_pval_T2{nsession} = pair_pval_T1;
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_corr_all.mat'),'context_corr_all')
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')

%% Plotting context selecitve ripple modulation







%% Plotting context selecitve spatial corr and ripple corr






%% 