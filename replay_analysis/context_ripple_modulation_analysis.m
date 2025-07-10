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


%% context-selective ripple modulation
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
windows = [-1.5 1.5];

z_V1_population_ripple_PSTH{1} = [];
z_V1_population_ripple_PSTH{2} = [];

context_modulation_all = struct();
% context_corr_all = struct();
tic
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

    % % Only grab ripple-modulated cells
    % V1_id = intersect(V1_id,find(ripple_modulation_PSTH_all{nsession}(1).ripple_modulation_percentile >= 0.95 | ripple_modulation_PSTH_all{nsession}(2).ripple_modulation_percentile >= 0.95...
    %     | ripple_modulation_PSTH_all{nsession}(1).modulation_percentile_PRE >=0.95 |ripple_modulation_PSTH_all{nsession}(2).modulation_percentile_PRE >=0.95));



    HPC_id = find(contains(all_regions,'HPC'));

    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));

    %%%%%%% Populational zscored firing rate difference


    V1_id_R = find(contains(all_regions,'V1_R'));
    T1_V1_cell =  find(session_clusters_all.mean_FR{nsession}(V1_id_R,1) >  session_clusters_all.mean_FR{nsession}(V1_id_R,2));
    V1_id_L = find(contains(all_regions,'V1_L'));
    T2_V1_cell =  find(session_clusters_all.mean_FR{nsession}(V1_id_L,1) <  session_clusters_all.mean_FR{nsession}(V1_id_L,2));
    % T1_HPC_cell =  session_clusters_all.mean_FR{nsession}(HPC_id,1) >  session_clusters_all.mean_FR{nsession}(HPC_id,2);
    % T2_HPC_cell =  session_clusters_all.mean_FR{nsession}(HPC_id,1) >  session_clusters_all.mean_FR{nsession}(HPC_id,2);
    % z_V1_population_ripple_PSTH{1} = [z_V1_population_ripple_PSTH{1}; squeeze(mean(ripple_modulation.PSTH_zscored(V1_id_R(T1_V1_cell),:,ripple_modulation.bins>-1 & ripple_modulation.bins<1),1,'omitnan'))];
    % z_V1_population_ripple_PSTH{2} = [z_V1_population_ripple_PSTH{2}; squeeze(mean(ripple_modulation.PSTH_zscored(V1_id_L(T2_V1_cell),:,ripple_modulation.bins>-1 & ripple_modulation.bins<1),1,'omitnan'))];

    %
    % T1_V1_cell =  find(session_clusters_all.mean_FR{nsession}(V1_id,1) >  session_clusters_all.mean_FR{nsession}(V1_id,2));
    % T2_V1_cell =  find(session_clusters_all.mean_FR{nsession}(V1_id,1) <  session_clusters_all.mean_FR{nsession}(V1_id,2));
    % % T1_HPC_cell =  session_clusters_all.mean_FR{nsession}(HPC_id,1) >  session_clusters_all.mean_FR{nsession}(HPC_id,2);
    % % T2_HPC_cell =  session_clusters_all.mean_FR{nsession}(HPC_id,1) >  session_clusters_all.mean_FR{nsession}(HPC_id,2);
    % z_V1_population_ripple_PSTH{1} = [z_V1_population_ripple_PSTH{1}; squeeze(mean(ripple_modulation.PSTH_zscored(V1_id(T1_V1_cell),:,ripple_modulation.bins>-1 & ripple_modulation.bins<1),1,'omitnan'))];
    % z_V1_population_ripple_PSTH{2} = [z_V1_population_ripple_PSTH{2}; squeeze(mean(ripple_modulation.PSTH_zscored(V1_id(T2_V1_cell),:,ripple_modulation.bins>-1 & ripple_modulation.bins<1),1,'omitnan'))];
    % 
    
    temp1 = squeeze(mean(ripple_modulation.PSTH(V1_id_R(T1_V1_cell),:,ripple_modulation.bins>-1 & ripple_modulation.bins<1),1,'omitnan'));
    temp1 = temp1-mean(temp1,"all",'omitnan')./std(temp1,0,'all','omitnan');

    temp2 = squeeze(mean(ripple_modulation.PSTH(V1_id_L(T2_V1_cell),:,ripple_modulation.bins>-1 & ripple_modulation.bins<1),1,'omitnan'));
    temp2 = temp2-mean(temp2,"all",'omitnan')./std(temp2,0,'all','omitnan');

    z_V1_population_ripple_PSTH{1} = [z_V1_population_ripple_PSTH{1}; temp1];
    z_V1_population_ripple_PSTH{2} = [z_V1_population_ripple_PSTH{2}; temp2];
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


        % Firing rates on both tracks
        context_modulation_all.FR_track{nsession}(1,nCell) = session_clusters_all.mean_FR{nsession}(nCell,1);
        context_modulation_all.FR_track{nsession}(2,nCell) = session_clusters_all.mean_FR{nsession}(nCell,2);

        % Z-scored track FR (relative to combined distribution)
        FR_distribution = reshape([session_clusters_all.spatial_response_raw{nsession}{nCell,1}; ...
            session_clusters_all.spatial_response_raw{nsession}{nCell,2}], 1, []);
        context_modulation_all.z_FR_track{nsession}(1,nCell) = ...
            (mean(mean(session_clusters_all.spatial_response_raw{nsession}{nCell,1})) - mean(FR_distribution)) ./ std(FR_distribution);

        context_modulation_all.z_FR_track{nsession}(2,nCell) = ...
            (session_clusters_all.mean_FR{nsession}(nCell,2) - mean(FR_distribution)) ./ std(FR_distribution);

        context_modulation_all.region{nsession}(nCell) = all_regions(nCell);

        context_modulation_all.session_id{nsession}(nCell) = nsession;
    end

    %%%%%%%%%%%%%%%%% Low ripples
    ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)]';
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');

    power_threshold = prctile(ripple_powers,[25 75]);
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & ripple_powers <= power_threshold(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & ripple_powers <= power_threshold(1));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_low_ripple{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_low_ripple{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_low_ripple{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_low_ripple{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_low_ripple{nsession}(2,nCell,:))';
    end

    %%%%%%%%%% High Ripples
    ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)]';
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');

    power_threshold = prctile(ripple_powers,[25 75]);
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & ripple_powers >= power_threshold(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & ripple_powers >= power_threshold(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_ripple{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_ripple{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_ripple{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_ripple{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_ripple{nsession}(2,nCell,:))';
    end


    context_modulation_all.ripple_modulation_percentile{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).ripple_modulation_percentile; ripple_modulation_PSTH_all{nsession}(2).ripple_modulation_percentile]);
    context_modulation_all.modulation_percentile_PRE{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).modulation_percentile_PRE; ripple_modulation_PSTH_all{nsession}(2).modulation_percentile_PRE]);
end
context_modulation_all.timebin = ripple_modulation.bins;

toc
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')

    % subplot(2,2,1)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.POST_ripple_FR(V1_id))
    % subplot(2,2,2)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
%     % 
%%  Spatial correlation and ripple correlation in V1 and HPC 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial corr vs ripple corr
% context_corr_all = struct();
% 
% for nsession = 1:length(sessions_to_process)
%     tic
%     all_clusters = session_clusters_all.spatial_cell_id{nsession};
% 
%     event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
%     event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];
% 
%     % ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession},session_clusters_all.spike_id{nsession},windows,psthBinSize,...
%     %     'unit_id',all_clusters,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);
%     % bins = ripple_modulation.bins > 0 & ripple_modulation.bins < 0.5;
% 
%     % all_clusters = session_clusters_all.spatial_cell_id{nsession};
%     all_regions = session_clusters_all.region{nsession};
%     V1_id = find(contains(all_regions,'V1'));
%     HPC_id = find(contains(all_regions,'HPC'));
% 
%     % Get mean bias for each ripple event in HPC and get within session
%     % Track 1 and Track 2 biased ripple events based on HPC bias
%     mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
%     log_odds_threshold = prctile(mean_bias,[20 80]);
%     T1_index = find(mean_bias > log_odds_threshold(2));
%     T2_index = find(mean_bias < log_odds_threshold(1));
%     %%%%%%%% Spatial Corr context diff and Ripple corr context diff
%     V1_cell_spatial_psth=[];
%     V1_cell_psth=[];
% 
%     HPC_cell_spatial_psth=[];
%     HPC_cell_psth=[];
% 
%     % Ensure parallel pool is open
%     if isempty(gcp('nocreate'))
%         parpool;
%     end
% 
% 
% 
%     [pair_corr, pair_pval,lags] = compute_ripple_xcorr_cell_pair(...
%         session_clusters_all.spike_times{nsession}, ...
%         session_clusters_all.spike_id{nsession}, ...
%         all_clusters(V1_id), all_clusters(HPC_id), event_times,'lag_range',[0],'window',0.2,'step', 0.2,'shuffle_option',0);
% 
%     context_corr_all.ripple_corr_all{nsession} = pair_corr;
%     context_corr_all.ripple_corr_pval_all{nsession} = pair_pval;
% 
%     [pair_corr_T1, pair_pval,lags] = compute_ripple_xcorr_cell_pair(...
%         session_clusters_all.spike_times{nsession}, ...
%         session_clusters_all.spike_id{nsession}, ...
%         all_clusters(V1_id), all_clusters(HPC_id), event_times(T1_index),'lag_range',[0],'window',0.2,'step', 0.2,'shuffle_option',0);
% 
%     context_corr_all.ripple_corr_T1{nsession} = pair_corr_T1;
%     context_corr_all.ripple_corr_pval_T1{nsession} = pair_pval;
% 
%     [pair_corr_T2,pair_pval,lags] = compute_ripple_xcorr_cell_pair(...
%         session_clusters_all.spike_times{nsession}, ...
%         session_clusters_all.spike_id{nsession}, ...
%         all_clusters(V1_id), all_clusters(HPC_id), event_times(T2_index),'lag_range',[0],'window',0.2,'step', 0.2,'shuffle_option',0);
% 
%     context_corr_all.ripple_corr_T2{nsession} = pair_corr_T2;
%     context_corr_all.ripple_corr_pval_T2{nsession} = pair_pval;
% 
% 
%     pair_pval_T1 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
%     pair_corr_T1 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
%     pair_pval_T2 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
%     pair_corr_T2 = nan(length(HPC_id),length(V1_id));           % 1D array (linear index)
% 
%     spatial_raw  = session_clusters_all.spatial_response_raw{nsession};
% 
%     HPC_spatial = cell(length(HPC_id), 2);
%     for mCell = 1:length(HPC_id)
%         HPC_spatial{mCell,1} = reshape(spatial_raw {HPC_id(mCell), 1}, 1, []);
%         HPC_spatial{mCell,2} = reshape(spatial_raw {HPC_id(mCell), 2}, 1, []);
%     end
% 
%     for nCell = 1:length(V1_id)
%         V1_spatial{1} = reshape(spatial_raw {V1_id(nCell), 1}, 1, []);
%         V1_spatial{2} = reshape(spatial_raw {V1_id(nCell), 2}, 1, []);
% 
%         % === Preload HPC cell responses for this loop ===
% 
% 
%         parfor mCell = 1:length(HPC_id)
%             [r1, p1] = corr(HPC_spatial{mCell,1}', V1_spatial{1}', 'type', 'Spearman', 'Rows', 'complete');
%             [r2, p2] = corr(HPC_spatial{mCell,2}', V1_spatial{2}', 'type', 'Spearman', 'Rows', 'complete');
% 
%             pair_pval_T1(mCell, nCell) = p1;
%             pair_corr_T1(mCell, nCell) = r1;
%             pair_pval_T2(mCell, nCell) = p2;
%             pair_corr_T2(mCell, nCell) = r2;
%         end
%     end
% 
% 
%     % pair_pval = reshape(pair_pval, nHPC, nV1); % Reshape to 2D [nHPC x nV1]
% 
%     context_corr_all.spatial_corr_T1{nsession}= pair_corr_T1;
%     context_corr_all.spatial_corr_pval_T1{nsession} = pair_pval_T1;
% 
%     context_corr_all.spatial_corr_T2{nsession} = pair_corr_T2;
%     context_corr_all.spatial_corr_pval_T2{nsession} = pair_pval_T1;
%     toc
% end
% 
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_corr_all.mat'),'context_corr_all')
% % save(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')


%% Plotting context selecitve ripple modulation
% scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')
load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')


ripple_modulation_percentile = [context_modulation_all.ripple_modulation_percentile{:}];
PRE_ripple_modulation_percentile = [context_modulation_all.modulation_percentile_PRE{:}];
ripple_modulation_percentile = max([ripple_modulation_percentile; PRE_ripple_modulation_percentile])>0.95;
ripple_modulation_id = ripple_modulation_percentile>0.95;

session_id_all = [context_modulation_all.session_id{:}];
regions_all = [context_modulation_all.region{:}];
regions_all(ripple_modulation_id==0)=nan;

FR_track = [context_modulation_all.FR_track{:}];
z_FR_track = [context_modulation_all.z_FR_track{:}];
z_FR_track_diff = z_FR_track(1,:)-z_FR_track(2,:);

FR_track = [context_modulation_all.FR_track{:}];
FR_track_diff = FR_track(1,:)-FR_track(2,:);


ripple_types = {'ALL','LOW','HIGH'};
PSTH_fields = {'PSTH_diff','PSTH_diff_low_ripple','PSTH_diff_high_ripple'};

% Preallocate
all_vars_V1 = struct();
all_vars_HPC = struct();

for iType = 1:numel(ripple_types)
    tag = ripple_types{iType};
    PSTH_data = vertcat(context_modulation_all.(PSTH_fields{iType}){:});

    % Compute FR diffs
    POST = mean(PSTH_data(:, context_modulation_all.timebin > 0 & context_modulation_all.timebin < 0.2), 2, 'omitnan');
    PRE = mean(PSTH_data(:, context_modulation_all.timebin > -0.2 & context_modulation_all.timebin < 0), 2, 'omitnan');
    SHIFT = mean(PSTH_data(:, context_modulation_all.timebin > -1 & context_modulation_all.timebin < -0.8), 2, 'omitnan');

    % Index
    isV1 = contains(regions_all,'V1');
    isHPC = contains(regions_all,'HPC');

    % Store for V1
    all_vars_V1.(['POST_' tag]) = normalize(double(POST(isV1)));
    all_vars_V1.(['PRE_' tag]) = normalize(double(PRE(isV1)));
    all_vars_V1.(['SHIFT_' tag]) = normalize(double(SHIFT(isV1)));

    % Store for HPC
    all_vars_HPC.(['POST_' tag]) = normalize(double(POST(isHPC)));
    all_vars_HPC.(['PRE_' tag]) = normalize(double(PRE(isHPC)));
    all_vars_HPC.(['SHIFT_' tag]) = normalize(double(SHIFT(isHPC)));
end

% Add z_FR_track_diff and subjectID
z_FR_V1 = normalize(double(z_FR_track_diff(contains(regions_all,'V1'))'));
z_FR_HPC = normalize(double(z_FR_track_diff(contains(regions_all,'HPC'))'));
subjectID_V1 = categorical(session_id_all(contains(regions_all,'V1'))');
subjectID_HPC = categorical(session_id_all(contains(regions_all,'HPC'))');

% Construct tables
tbl_V1 = table(z_FR_V1, ...
    all_vars_V1.POST_ALL, all_vars_V1.POST_LOW, all_vars_V1.POST_HIGH, ...
    all_vars_V1.PRE_ALL, all_vars_V1.PRE_LOW, all_vars_V1.PRE_HIGH, ...
    all_vars_V1.SHIFT_ALL, all_vars_V1.SHIFT_LOW, all_vars_V1.SHIFT_HIGH, ...
    subjectID_V1, ...
    'VariableNames', {'z_FR_track_diff', ...
    'POST_ripple_FR_diff_ALL','POST_ripple_FR_diff_LOW','POST_ripple_FR_diff_HIGH', ...
    'PRE_ripple_FR_diff_ALL','PRE_ripple_FR_diff_LOW','PRE_ripple_FR_diff_HIGH', ...
    'SHIFT_ripple_FR_diff_ALL','SHIFT_ripple_FR_diff_LOW','SHIFT_ripple_FR_diff_HIGH', ...
    'subjectID'});

tbl_HPC = table(z_FR_HPC, ...
    all_vars_HPC.POST_ALL, all_vars_HPC.POST_LOW, all_vars_HPC.POST_HIGH, ...
    all_vars_HPC.PRE_ALL, all_vars_HPC.PRE_LOW, all_vars_HPC.PRE_HIGH, ...
    all_vars_HPC.SHIFT_ALL, all_vars_HPC.SHIFT_LOW, all_vars_HPC.SHIFT_HIGH, ...
    subjectID_HPC, ...
    'VariableNames', {'z_FR_track_diff', ...
    'POST_ripple_FR_diff_ALL','POST_ripple_FR_diff_LOW','POST_ripple_FR_diff_HIGH', ...
    'PRE_ripple_FR_diff_ALL','PRE_ripple_FR_diff_LOW','PRE_ripple_FR_diff_HIGH', ...
    'SHIFT_ripple_FR_diff_ALL','SHIFT_ripple_FR_diff_LOW','SHIFT_ripple_FR_diff_HIGH', ...
    'subjectID'});

ModelList = {
    'POST_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)';
    'POST_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)';
    'POST_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)';

    'PRE_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)';
    'PRE_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)';
    'PRE_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)';

    'SHIFT_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)';
    'SHIFT_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)';
    'SHIFT_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)';

    'z_FR_track_diff ~ POST_ripple_FR_diff_ALL + PRE_ripple_FR_diff_ALL + SHIFT_ripple_FR_diff_ALL + (1|subjectID)';
    'z_FR_track_diff ~ POST_ripple_FR_diff_LOW + PRE_ripple_FR_diff_LOW + SHIFT_ripple_FR_diff_LOW + (1|subjectID)';
    'z_FR_track_diff ~ POST_ripple_FR_diff_HIGH + PRE_ripple_FR_diff_HIGH + SHIFT_ripple_FR_diff_HIGH + (1|subjectID)';

    'z_FR_track_diff ~ POST_ripple_FR_diff_LOW + POST_ripple_FR_diff_HIGH + (1|subjectID)';
    'z_FR_track_diff ~ PRE_ripple_FR_diff_LOW + PRE_ripple_FR_diff_HIGH + (1|subjectID)';
    'z_FR_track_diff ~ SHIFT_ripple_FR_diff_LOW + SHIFT_ripple_FR_diff_HIGH + (1|subjectID)';
    };

clear output
tbl_all = {tbl_V1,tbl_HPC};
for nregion = 1:2
    parfor iBoot = 1:1000
        local_b = cell(length(ModelList),1);
        local_tstat = cell(length(ModelList),1);
        local_pval = cell(length(ModelList),1);
        local_variable = cell(length(ModelList),1);
        local_R2 = zeros(length(ModelList),1);

        tic


        for m = 1:length(ModelList)
            tbl = tbl_all{nregion};
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            index = randsample(s, 1:height(tbl), height(tbl), true);
            tbl = tbl(index,:);

            glme = fitlme(tbl, ModelList{m});
            local_b{m} = glme.Coefficients.Estimate(2:end);
            local_tstat{m} = glme.Coefficients.tStat(2:end);
            local_pval{m} = glme.Coefficients.pValue(2:end);
            local_variable{m} = glme.CoefficientNames(2:end);

            if isprop(glme, 'Rsquared') && isfield(glme.Rsquared, 'Adjusted')
                local_R2(m) = glme.Rsquared.Adjusted;
            else
                local_R2(m) = NaN;  % fallback
            end
        end

        output(iBoot).b = local_b;
        output(iBoot).R2 = local_R2;
        output(iBoot).t_stat = local_tstat;
        output(iBoot).pval = local_pval;
        output(iBoot).variable = local_variable;
        output(iBoot).model = ModelList;
        %     output(iBoot).type = modelNames;
        toc
    end
    if nregion == 1
        output_V1 = output;
    else
        output_HPC = output;
    end
   
end
 save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_ripple_modulation_glme.mat'),'output_V1','output_HPC');
% z_FR_track_diff(contains(regions_all,'V1'))



colorlines = [ ...
    0.85, 0.2, 0.2;  % red
    0.2, 0.4, 0.8    % blue
];

nfig = figure;
nfig.Name = 'Context selective FR difference in V1 and HPC';
nfig.Position = [623          88        1215         846];

nexttile
% Assume z_FR_track is loaded and size is [2 x N]
a = FR_track(1,:);  % Track L firing rate
b = FR_track(2,:);  % Track R firing rate
% a = z_FR_track(1,:);  % Track L firing rate
% b = z_FR_track(2,:);  % Track R firing rate
diff_ab = a - b;

% Define plot limits
% lims = [min([a b]) max([a b])];
% lims = lims + [-1 1]*0.1*range(lims); % padding

% Plot scatter
is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
hold on
% scatter(a(contains(regions_all,'V1_R')), b(contains(regions_all,'V1_R')), 'b','filled', 'MarkerFaceAlpha', 0.2)
% scatter(a(contains(regions_all,'V1_L')), b(contains(regions_all,'V1_L')), 'r','filled', 'MarkerFaceAlpha', 0.2)
scatter(a, b, 15, scolors, 'filled', 'MarkerFaceAlpha', 0.2)

% legend('V1 R','V1 L','box','off')
title('Track L - Track R FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

plot([0 prctile([a b],99)], [0 prctile([a b],99)], 'k--', 'LineWidth', 1.5)  % identity line
plot([0 50], [50 0], 'k--', 'LineWidth', 1.5)  % identity line
scatter([10 20 30 40 50], [50 40 30 20 10], 20, 'k', 'filled', 'MarkerFaceAlpha', 1)
axis equal
xlim([0 60])
ylim([0 60])

% xlim(lims); ylim(lims)
xlabel('Track L FR (Hz)')
ylabel('Track R FR (Hz)')
title('Track L vs R with (L - R) projection')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal

nexttile
FR_diff = a-b;
histogram(FR_diff(contains(regions_all,'V1_R')),[-20:0.5:20],'Normalization','probability','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
histogram(FR_diff(contains(regions_all,'V1_L')),[-20:0.5:20],'Normalization','probability','FaceColor',colorlines(2,:),'EdgeColor','none');
[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'HPC_R')),FR_diff(contains(regions_all,'HPC_L')));
text(prctile(a,90),0.1,sprintf('p = %.3e',pval))
ylim([0 0.12])
xline(0,'k--')
xlabel('Track L - Track R FR diff (Hz)')
ylabel('proportion of cells')
legend('V1 R','V1 L','box','off')
title('Track L - Track R FR diff')
set(gca,'TickDir','out','Box','off','FontSize',12)

% axis equal
nexttile
hold on;
FR_diff = a-b;
histogram(FR_diff(contains(regions_all,'V1_L')),[-20:0.5:20],'Normalization','cdf','FaceColor',colorlines(2,:),'EdgeColor','none');
histogram(FR_diff(contains(regions_all,'V1_R')),[-20:0.5:20],'Normalization','cdf','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;

ylim([0 1])
xline(0,'k--')
xlabel('Track L - Track R FR diff (Hz)')
ylabel('proportion of cells')
legend('V1 L','V1 R','box','off')
title('Track L - Track R FR diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal


nexttile
a = z_FR_track(1,:);  % Track L firing rate
b = z_FR_track(2,:);  % Track R firing rate
FR_diff = a-b;
histogram(FR_diff(contains(regions_all,'V1_R')),[-2:0.1:2],'Normalization','probability','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
histogram(FR_diff(contains(regions_all,'V1_L')),[-2:0.1:2],'Normalization','probability','FaceColor',colorlines(2,:),'EdgeColor','none');

[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'V1_R')),FR_diff(contains(regions_all,'V1_L')));
text(prctile(a,90),0.1,sprintf('p = %.3e',pval))
ylim([0 0.12])
xline(0,'k--')
xlabel('Track L - Track R FR diff (z)')
ylabel('cum prop of cells')
legend('V1 R','V1 L','box','off')
title('Track L - Track R FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal

nexttile
a = z_FR_track(1,:);  % Track L firing rate
b = z_FR_track(2,:);  % Track R firing rate
FR_diff = a-b;
hold on;
histogram(FR_diff(contains(regions_all,'V1_L')),[-2:0.1:2],'Normalization','cdf','FaceColor',colorlines(2,:),'EdgeColor','none');
histogram(FR_diff(contains(regions_all,'V1_R')),[-2:0.1:2],'Normalization','cdf','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
ylim([0 1])
xline(0,'k--')
xlabel('Track L - Track R FR diff (z)')
ylabel('cum prop of cells')
legend('V1 L','V1 R','box','off')
title('Track L - Track R FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal


nexttile
% Plot scatter
% Assume z_FR_track is loaded and size is [2 x N]
a = FR_track(1,:);  % Track L firing rate
b = FR_track(2,:);  % Track R firing rate
% a = z_FR_track(1,:);  % Track L firing rate
% b = z_FR_track(2,:);  % Track R firing rate
diff_ab = a - b;
is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
hold on
% scatter(a(contains(regions_all,'V1_R')), b(contains(regions_all,'V1_R')), 'b','filled', 'MarkerFaceAlpha', 0.2)
% scatter(a(contains(regions_all,'V1_L')), b(contains(regions_all,'V1_L')), 'r','filled', 'MarkerFaceAlpha', 0.2)
scatter(a, b, 15, scolors, 'filled', 'MarkerFaceAlpha', 0.2)

% legend('V1 R','V1 L','box','off')
title('Track L - Track R FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

plot([0 prctile([a b],99)], [0 prctile([a b],99)], 'k--', 'LineWidth', 1.5)  % identity line
plot([0 50], [50 0], 'k--', 'LineWidth', 1.5)  % identity line
scatter([10 20 30 40 50], [50 40 30 20 10], 20, 'k', 'filled', 'MarkerFaceAlpha', 1)
axis equal
xlim([0 60])
ylim([0 60])


% xlim(lims); ylim(lims)
xlabel('Track L FR (Hz)')
ylabel('Track R FR (Hz)')
title('Track L vs R with (L - R) projection')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal

nexttile
FR_diff = a-b;
histogram(FR_diff(contains(regions_all,'HPC_R')),[-20:0.5:20],'Normalization','probability','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
histogram(FR_diff(contains(regions_all,'HPC_L')),[-20:0.5:20],'Normalization','probability','FaceColor',colorlines(2,:),'EdgeColor','none');

[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'HPC_R')),FR_diff(contains(regions_all,'HPC_L')));
text(prctile(a,90),0.1,sprintf('p = %.3e',pval))
ylim([0 0.2])
xline(0,'k--')
xlabel('Track L - Track R FR diff (Hz)')
ylabel('proportion of cells')
legend('HPC R','HPC L','box','off')
title('Track L - Track R FR diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal


nexttile
FR_diff = a-b;
hold on;histogram(FR_diff(contains(regions_all,'HPC_L')),[-20:0.5:20],'Normalization','cdf','FaceColor',colorlines(2,:),'EdgeColor','none');
histogram(FR_diff(contains(regions_all,'HPC_R')),[-20:0.5:20],'Normalization','cdf','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
ylim([0 1])
xline(0,'k--')
xlabel('Track L - Track R FR diff (Hz)')
ylabel('cum prop of cells')
legend('HPC L','HPC R','box','off')
title('Track L - Track R FR diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal


nexttile
a = z_FR_track(1,:);  % Track L firing rate
b = z_FR_track(2,:);  % Track R firing rate
FR_diff = a-b;
histogram(FR_diff(contains(regions_all,'HPC_R')),[-2:0.1:2],'Normalization','probability','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
histogram(FR_diff(contains(regions_all,'HPC_L')),[-2:0.1:2],'Normalization','probability','FaceColor',colorlines(2,:),'EdgeColor','none');
[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'HPC_R')),FR_diff(contains(regions_all,'HPC_L')));
text(prctile(a,90),0.1,sprintf('p = %.3e',pval))
ylim([0 0.2])
xline(0,'k--')
xlabel('Track L - Track R FR diff (z)')
ylabel('proportion of cells')
legend('HPC R','HPC L','box','off')
title('Track L - Track R FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal


nexttile
a = z_FR_track(1,:);  % Track L firing rate
b = z_FR_track(2,:);  % Track R firing rate
FR_diff = a-b;

hold on;histogram(FR_diff(contains(regions_all,'HPC_L')),[-2:0.1:2],'Normalization','cdf','FaceColor',colorlines(2,:),'EdgeColor','none');
histogram(FR_diff(contains(regions_all,'HPC_R')),[-2:0.1:2],'Normalization','cdf','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
ylim([0 1])
xline(0,'k--')
xlabel('Track L - Track R FR diff (z)')
ylabel('cum prop of cells')
legend('HPC L','HPC R','box','off')
title('Track L - Track R FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal






%%%%%%%%%%% ALL ripples

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression';
nfig.Position = [   842   345   954   578];




subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
xline(0,'k--')
yline(0,'k--')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.45 0.45])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'V1')));

scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-1 1])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.7 0.7])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



%%%%%%%%%%%%% Low ripple


all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression {low ripples}';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
xline(0,'k--')
yline(0,'k--')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.45 0.45])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'V1')));

scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-1 1])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.7 0.7])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)






%%%%%%%%%%%%% High ripple

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression {high ripples}';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
xline(0,'k--')
yline(0,'k--')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.45 0.45])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'V1')));

scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-1 1])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.7 0.7])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



%%%% Loop though bins


%%%%%%%%%%% kstest2
ks2stat_all = nan(1000,10);

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 (high ripples bins) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');

% timewindows = -1:0.2:1;
timewindows = -1:0.2:1;
for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
    text(0.4,0.05,sprintf('p = %.3e',pval))
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ylim([0 1])
    xline(0,'k--','LineWidth',2)
    xlabel('Track L - Track R Ripple FR diff (z)')
    ylabel('cum prop of cells')
    legend('V1 Track R prefering','V1 Track L prefering','box','off')
    title(sprintf('%.1f to %.1fs relative to ripples',timewindows(nbin),timewindows(nbin+1)))
    set(gca,'TickDir','out','Box','off','FontSize',12)
end



%%%%%%%%%%% kstest2
ks2stat_all = nan(1000,10);

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 (low ripples bins) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');

timewindows = -1:0.2:1;

for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
    text(0.4,0.05,sprintf('p = %.3e',pval))
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ylim([0 1])
    xline(0,'k--','LineWidth',2)
    xlabel('Track L - Track R Ripple FR diff (z)')
    ylabel('cum prop of cells')
    legend('V1 Track R prefering','V1 Track L prefering','box','off')
    title(sprintf('%.1f to %.1fs relative to ripples',timewindows(nbin),timewindows(nbin+1)))
    set(gca,'TickDir','out','Box','off','FontSize',12)
   

end



%%%%%%%%%%% KS max difference bootstrap
bin_width = 0.1;    % 50 ms
step_size = 0.01;   % 10 ms
t_start = -0.5;
t_end = 0.5;

% Bin centers
bin_centers = t_start:step_size:t_end;

% Bin edges (each bin spans 100 ms centered at bin_centers)
timewindows = [bin_centers - bin_width/2; bin_centers + bin_width/2];
ks2stat_all = nan(1000,length(timewindows)-1);

%%% Low ripples
all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
for nbin = 1:length(timewindows)-1
    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(1,nbin)&context_modulation_all.timebin<timewindows(2,nbin)),2,'omitnan');

    T1_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0);
    T2_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0);
    num_neurons = min([length(T1_FR_dff) length(T2_FR_dff)]);


    parfor iBoot = 1:1000
        temp1 = [];temp2 = [];
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index1 = randsample(s, 1:length(T1_FR_dff), num_neurons, true);
        temp1 = T1_FR_dff(index1);
        index2 = randsample(s, 1:length(T2_FR_dff), num_neurons, true);
        temp2 = T2_FR_dff(index2);
        [h,pval,ks2stat_boot] = kstest2(temp1,temp2);
        ks2stat_all(iBoot,nbin) = ks2stat_boot;
    end
end
ks2stat_low_ripples = ks2stat_all;

%%%% High ripples
all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
ks2stat_all = nan(1000,length(timewindows)-1);

for nbin = 1:length(timewindows)-1
    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(1,nbin)&context_modulation_all.timebin<timewindows(2,nbin)),2,'omitnan');

    T1_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0);
    T2_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0);
    num_neurons = min([length(T1_FR_dff) length(T2_FR_dff)]);


    parfor iBoot = 1:1000
        temp1 = [];temp2 = [];
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index1 = randsample(s, 1:length(T1_FR_dff), num_neurons, true);
        temp1 = T1_FR_dff(index1);
        index2 = randsample(s, 1:length(T2_FR_dff), num_neurons, true);
        temp2 = T2_FR_dff(index2);
        [h,pval,ks2stat_boot] = kstest2(temp1,temp2);
        ks2stat_all(iBoot,nbin) = ks2stat_boot;
    end
end
ks2stat_high_ripples = ks2stat_all;



nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 KS difference time series (high vs low ripple)';
nfig.Position = [ 1150         322         363         282];
timewindows = [bin_centers - bin_width/2; bin_centers + bin_width/2];
% x = timewindows(1:end-1)+mean(diff(timewindows))/2;
x = mean(timewindows);
x = x(1:end-1)+mean(diff(x))/2;
colour_lines = [158,202,225;33,113,181]/256;% two blue 

nexttile
% x = timewindows(1:end-1)+mean(diff(timewindows))/2;
plot(x,mean(ks2stat_low_ripples),'Color',colour_lines(1,:));hold on;
ci_low  = prctile(ks2stat_low_ripples, 2.5, 1);
ci_high = prctile(ks2stat_low_ripples, 97.5, 1);

% Fill 95% CI
F(1) = fill([x, fliplr(x)], [ci_low, fliplr(ci_high)], ...
    colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

plot(x,mean(ks2stat_high_ripples),'Color',colour_lines(end,:));hold on;
ci_low  = prctile(ks2stat_high_ripples, 2.5, 1);
ci_high = prctile(ks2stat_high_ripples, 97.5, 1);

% Fill 95% CI
F(2)=fill([x, fliplr(x)], [ci_low, fliplr(ci_high)], ...
    colour_lines(end,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
set(gca,'TickDir','out','Box','off','FontSize',12)
xticks([-0.5 -0.25 0 0.25 0.5])
xline(0,'r')

legend(F(1:2),{'low ripples','high ripples'},'box','off')
xlabel('Time relative to ripple onset (s)')
ylabel('Maximum empirical cumulative distribution difference');

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])

%% Plotting context selecitve ripple modulation (Mixed effect regression)

%%%%%%%%%%%%% 

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_ripple_modulation_glme.mat'),'output_V1','output_HPC');
% summary_table = UP_DOWN_ripples_glme_summary_table(output_V1,[1:15],'context_ripple_modulation_glme_summary.csv');
summary_table = UP_DOWN_ripples_glme_summary_table(output_V1,[1:15]);

%%%%%%%%% V1 all vs low ripple vs high ripple
model_index == 11

%%%%%%%%% V1 all vs low ripple vs high ripple
% b_all = []; t_all = []; R2_all = []; pval_all = [];
b_low = []; t_low = []; R2_low = []; pval_low = [];
b_high = []; t_high = []; R2_high = []; pval_high = [];

for nBoot = 1:1000
    % b_all(nBoot) = output_all(nBoot).b{7};
    % t_all(nBoot) = output_all(nBoot).t_stat{7};
    % R2_all(nBoot) = output_all(nBoot).R2(7);
    % pval_all(nBoot) = output_all(nBoot).pval{7};

    b_low(nBoot) = output_V1(nBoot).b{2};
    t_low(nBoot) = output_V1(nBoot).t_stat{2};
    R2_low(nBoot) = output_V1(nBoot).R2(2);
    pval_low(nBoot) = output_V1(nBoot).pval{2};


    b_high(nBoot) = output_V1(nBoot).b{3};
    t_high(nBoot) = output_V1(nBoot).t_stat{3};
    R2_high(nBoot) = output_V1(nBoot).R2(3);
    pval_high(nBoot) = output_V1(nBoot).pval{3};
end

R2_diff = R2_high - R2_low;
ci = prctile(R2_diff, [2.5 97.5]);
barData = mean(R2_diff);


ci = prctile(t_V1delta, [2.5 97.5]);

nfig = figure;
nfig.Name = 'Context selective ripple modulation glme regression';
nfig.Position = [71 280  1815  662];
bar([1 2 3], barData, 0.4, 'FaceColor', customColors(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
errorbar([1 2 3], barData, barData - ci(1), barData - ci(2), 'k', 'LineWidth', 1.5)
xlim([0.5 1.5]); xticks(1); xticklabels({'V1 delta power'})
ylabel('Bootstrapped t-statistic')
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off





%%%%%%%%% V1 
b_all = []; t_all = []; R2_all = []; pval_all = [];
b_low = []; t_low = []; R2_low = []; pval_low = [];
b_high = []; t_high = []; R2_high = []; pval_high = [];

for nBoot = 1:1000
    b_all(nBoot) = output_all(nBoot).b{1};
    t_all(nBoot) = output_all(nBoot).t_stat{1};
    R2_all(nBoot) = output_all(nBoot).R2(1);
    pval_all(nBoot) = output_all(nBoot).pval{1};

    b_low(nBoot) = output_low(nBoot).b{1};
    t_low(nBoot) = output_low(nBoot).t_stat{1};
    R2_low(nBoot) = output_low(nBoot).R2(1);
    pval_low(nBoot) = output_low(nBoot).pval{1};


    b_high(nBoot) = output_high(nBoot).b{1};
    t_high(nBoot) = output_high(nBoot).t_stat{1};
    R2_high(nBoot) = output_high(nBoot).R2(1);
    pval_high(nBoot) = output_high(nBoot).pval{1};
end






%%%%%%%%% Baseline
b_all = []; t_all = []; R2_all = []; pval_all = [];
b_low = []; t_low = []; R2_low = []; pval_low = [];
b_high = []; t_high = []; R2_high = []; pval_high = [];

for nBoot = 1:1000
    b_all(nBoot) = output_all(nBoot).b{1};
    t_all(nBoot) = output_all(nBoot).t_stat{1};
    R2_all(nBoot) = output_all(nBoot).R2(1);
    pval_all(nBoot) = output_all(nBoot).pval{1};

    b_low(nBoot) = output_low(nBoot).b{1};
    t_low(nBoot) = output_low(nBoot).t_stat{1};
    R2_low(nBoot) = output_low(nBoot).R2(1);
    pval_low(nBoot) = output_low(nBoot).pval{1};


    b_high(nBoot) = output_high(nBoot).b{1};
    t_high(nBoot) = output_high(nBoot).t_stat{1};
    R2_high(nBoot) = output_high(nBoot).R2(1);
    pval_high(nBoot) = output_high(nBoot).pval{1};
end
%% Plotting context selecitve spatial corr and ripple corr




%% 