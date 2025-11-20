%% How reactivation coherence modulates V1 activity
%%%%%
%%%%% Does Neurons active before 
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
% load(fullfile(analysis_folder,'ripples_TF_stats_POST.mat'))
%     load(fullfile(analysis_folder,'session_clusters_all_POST.mat'),'session_clusters_all','-v7.3')
sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);


%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end

%% V1 activity before and after coherent reactivations
load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')
% load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','RRR_reactivation_ripples_PSTH.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripple_normalised_duration.mat'),'DOWN_normalised_duration','UP_normalised_duration','SO_normalised_duration');

%  = struct();

%%%% KDE bias
timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;
bins_to_use = bin_centers>0 & bin_centers<0.1;
% bins_to_use

session_id = [ripples_all(1).session_count(ripples_all(1).SWS_index); ripples_all(2).session_count(ripples_all(2).SWS_index)];
%%% PLS KDE
z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);
clear z_bias1

z_bias_KDE = z_bias;
z_bias_V1_KDE = z_bias_V1;

%%% RRR 
z_bias = RRR_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = RRR_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);
clear z_bias1
z_bias_RRR = z_bias;
z_bias_V1_RRR = z_bias_V1;

psthBinSize = 0.01;
windows = [-1.5 1.5];

timebins = [5 6 10]; % Timebin of the LFP metric relative to ripples where 1 is -1 to -0.8s and 10 is 0.8 to 1s


z_V1_population_ripple_PSTH{1} = [];
z_V1_population_ripple_PSTH{2} = [];

V1_reactivation_modulation_all = struct();
% context_corr_all = struct();
tic
for nsession = 1:length(sessions_to_process)

    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    %%%%%%%%%% HPC
    cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
    spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

    spike_id = session_clusters_all.spike_id{nsession}(spike_index);
    spike_id(spike_id>0)=1;

    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
        'unit_id',1,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);

    ripple_MUA_PSTH_all{1}{nsession}(hemi,:,:) = zscore(squeeze(ripple_modulation(1).PSTH),0,2);
    ripple_MUA_PSTH_all{2}{nsession}(hemi,:,:) = zscore(squeeze(ripple_modulation(2).PSTH),0,2);


    spike_id = session_clusters_all.spike_id{nsession}(spike_index);
    %         spike_id(spike_id>0)=1;

    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
        'unit_id',cell_id,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);


    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias_KDE(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[25 75]);
    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));

    mean_bias = mean(z_bias_V1_KDE(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[25 75]);
    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));

    bin_centers>0 & bin_centers<0.1;

    z_bias_KDE

    z_bias_V1_KDE

end
