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
load(fullfile(analysis_folder,'ripples_TF_stats_POST.mat'))
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


%% context-selective ripple modulation based on V1 bias
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
bins_to_use = bin_centers>0 & bin_centers<0.2;


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

timebins = [5 6 10]; % Timebin of the LFP metric relative to ripples where 1 is -1 to -0.8s and 10 is 0.8 to 1s


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

    % Get mean bias for each ripple event in V1 and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias_V1(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
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
    amplitudes = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)]';
    % mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');

    power_threshold = prctile(amplitudes,[25 75]);
    log_odds_threshold = prctile(mean_bias,[20 80]);

    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes <= power_threshold(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes <= power_threshold(1));

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
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes >= power_threshold(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes >= power_threshold(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_ripple{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_ripple{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_ripple{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_ripple{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_ripple{nsession}(2,nCell,:))';
    end


    %%%%%%%%%%%% spindles amplitude
    amplitudes = [];
    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        amplitudes = [amplitudes ripples_all(nprobe).spindle_amplitude_ripple_onset{nsession}(cortex_ref_shank(nsession,:),ia)];
    end

    power_threshold1 = prctile(amplitudes(1,:),[25 75]);
    power_threshold2 = prctile(amplitudes(2,:),[25 75]);

    %%%% High
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) >= power_threshold2(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) >= power_threshold1(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_spindle{nsession}(2,nCell,:))';
    end

    %%%% Low
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) <= power_threshold2(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) <= power_threshold1(1));


    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_low_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_low_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_low_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_low_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_low_spindle{nsession}(2,nCell,:))';
    end


    %%%%%%%%%%%% SO power
    amplitudes=[];

    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        amplitudes = [amplitudes ripples_all(nprobe).SO_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),ia)];
    end

    power_threshold1 = prctile(amplitudes(1,:),[25 75]);
    power_threshold2 = prctile(amplitudes(2,:),[25 75]);

    %%%% High
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) >= power_threshold2(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) >= power_threshold1(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_SO{nsession}(2,nCell,:))';
    end

    %%%% Low
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) <= power_threshold2(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) <= power_threshold1(1));


    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_low_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_low_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_low_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_low_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_low_SO{nsession}(2,nCell,:))';
    end



    %%%%%%%%%%%% SO phase
    event_phase=[];

    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        event_phase = [event_phase ripples_all(nprobe).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),ia)];
    end
    
    is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
    is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_peak_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_peak_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_peak_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_peak_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_peak_SO{nsession}(2,nCell,:))';
    end



    is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
    is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_trough_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_trough_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_trough_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_trough_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_trough_SO{nsession}(2,nCell,:))';
    end



    %%%%%%%%%%%% spindle phase
    event_phase=[];

    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        event_phase = [event_phase ripples_all(nprobe).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),ia)];
    end
    
    is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
    is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_peak_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_peak_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_peak_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_peak_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_peak_spindle{nsession}(2,nCell,:))';
    end



    is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
    is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_trough_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_trough_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_trough_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_trough_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_trough_spindle{nsession}(2,nCell,:))';
    end




    hemi_id = [ones(1,sum((ripples_all(1).SWS_index==1 & ripples_all(1).session_count==nsession)>0)) 2*ones(1,sum((ripples_all(2).SWS_index==1 & ripples_all(2).session_count==nsession)>0))];

    %%%%%%%%%%%% SO power (different timing)
    log_odds_threshold = prctile(mean_bias,[20 80]);

    for n = 1:3
        temp = [squeeze(ripples_TF_stats.V1_amp{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_amp{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_amp{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_amp{nsession}(1,1,timebins(n),hemi_id == 2))'];
        temp_percentile = tiedrank(temp.').' / size(temp, 2);

        %%%%% LOW
        T1_index = find(mean_bias > log_odds_threshold(2) & temp_percentile(2,:) <= 0.25);
        T2_index = find(mean_bias < log_odds_threshold(1) & temp_percentile(1,:) <= 0.25);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_low_SO_time{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_low_SO_time{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_low_SO_time{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_low_SO_time{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_low_SO_time{nsession}(2,nCell,:))';
        end

        %%%%% High
        T1_index = find(mean_bias > log_odds_threshold(2) & temp_percentile(2,:) >= 0.75);
        T2_index = find(mean_bias < log_odds_threshold(1) & temp_percentile(1,:) >= 0.75);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_high_SO_time{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_high_SO_time{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_high_SO_time{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_high_SO_time{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_high_SO_time{nsession}(2,nCell,:))';
        end

    end


    %%%%%%%%%%%% SO phase with different times
    event_phase = [];
    for n = 1:3
        event_phase = [squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 2))'];
        % temp_percentile = tiedrank(temp.').' / size(temp, 2);

        is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
        is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_peak_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_peak_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_peak_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_peak_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_peak_SO_time{n}{nsession}(2,nCell,:))';
        end

   

        is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
        is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_trough_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_trough_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_trough_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_trough_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_trough_SO_time{n}{nsession}(2,nCell,:))';
        end
    end

    %%%%%%%%%%%% SO phase (Sync)

    event_phase = [];
    for n = 1:3
        event_phase = [squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 2))'];
        % temp_percentile = tiedrank(temp.').' / size(temp, 2);

        is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
        is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;

        T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1 & is_peak_phase_1 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1 & is_peak_phase_2 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_sync_peak_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(2,nCell,:))';
        end

   

        is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
        is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1 & is_trough_phase_1 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1 & is_trough_phase_2 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_sync_trough_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(2,nCell,:))';
        end
    end


    %%%%%%%%%%%% SO phase (Anti phase)
    event_phase = [];
    for n = 1:3
        event_phase = [squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 2))'];
        % temp_percentile = tiedrank(temp.').' / size(temp, 2);

        is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
        is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;

        T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1 & is_peak_phase_1 == 0);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1 & is_peak_phase_2 == 0);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_dsync_peak_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(2,nCell,:))';
        end



        is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
        is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1 & is_trough_phase_1 == 0);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1 & is_trough_phase_2 == 0);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_dsync_trough_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(2,nCell,:))';
        end
    end

    %%%%%%%%%%%% SO V1 HPC PLV
    temp = [];
    for n = 1:3
        temp = [squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,1,timebins(n),hemi_id == 2))'];
        temp_percentile = tiedrank(temp.').' / size(temp, 2);

        SO_index1 = temp_percentile(1,:) < 0.25;
        SO_index2 = temp_percentile(2,:) < 0.25;

        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_low_V1_HPC_SO_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:))';
        end

        SO_index1 = temp_percentile(1,:) > 0.75;
        SO_index2 = temp_percentile(2,:) > 0.75;


        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_high_V1_HPC_SO_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:))';
        end
    end

    %%%%%%%%%%%% Spindle V1 HPC PLV

    temp = [];
    for n = 1:3
        temp = [squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,3,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,3,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,3,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,3,timebins(n),hemi_id == 2))'];
        temp_percentile = tiedrank(temp.').' / size(temp, 2);

        SO_index1 = temp_percentile(1,:) < 0.25;
        SO_index2 = temp_percentile(2,:) < 0.25;

        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_low_V1_HPC_spindle_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:))';
        end

        SO_index1 = temp_percentile(1,:) > 0.75;
        SO_index2 = temp_percentile(2,:) > 0.75;


        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_high_V1_HPC_spindle_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:))';
        end
    end



    context_modulation_all.ripple_modulation_percentile{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).ripple_modulation_percentile; ripple_modulation_PSTH_all{nsession}(2).ripple_modulation_percentile]);
    context_modulation_all.modulation_percentile_PRE{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).modulation_percentile_PRE; ripple_modulation_PSTH_all{nsession}(2).modulation_percentile_PRE]);
end
context_modulation_all.timebin = ripple_modulation.bins;

toc

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all_V1_bias.mat'),'context_modulation_all')







%% context-selective ripple modulation based on V1 bias
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
bins_to_use = bin_centers>-0.2 & bin_centers<0;


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

timebins = [5 6 10]; % Timebin of the LFP metric relative to ripples where 1 is -1 to -0.8s and 10 is 0.8 to 1s


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

    % Get mean bias for each ripple event in V1 and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias_V1(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
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
    amplitudes = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)]';
    % mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');

    power_threshold = prctile(amplitudes,[25 75]);
    log_odds_threshold = prctile(mean_bias,[20 80]);

    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes <= power_threshold(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes <= power_threshold(1));

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
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes >= power_threshold(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes >= power_threshold(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_ripple{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_ripple{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_ripple{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_ripple{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_ripple{nsession}(2,nCell,:))';
    end


    %%%%%%%%%%%% spindles amplitude
    amplitudes = [];
    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        amplitudes = [amplitudes ripples_all(nprobe).spindle_amplitude_ripple_onset{nsession}(cortex_ref_shank(nsession,:),ia)];
    end

    power_threshold1 = prctile(amplitudes(1,:),[25 75]);
    power_threshold2 = prctile(amplitudes(2,:),[25 75]);

    %%%% High
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) >= power_threshold2(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) >= power_threshold1(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_spindle{nsession}(2,nCell,:))';
    end

    %%%% Low
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) <= power_threshold2(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) <= power_threshold1(1));


    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_low_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_low_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_low_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_low_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_low_spindle{nsession}(2,nCell,:))';
    end


    %%%%%%%%%%%% SO power
    amplitudes=[];

    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        amplitudes = [amplitudes ripples_all(nprobe).SO_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),ia)];
    end

    power_threshold1 = prctile(amplitudes(1,:),[25 75]);
    power_threshold2 = prctile(amplitudes(2,:),[25 75]);

    %%%% High
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) >= power_threshold2(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) >= power_threshold1(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_high_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_high_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_high_SO{nsession}(2,nCell,:))';
    end

    %%%% Low
    log_odds_threshold = prctile(mean_bias,[20 80]);
    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes(2,:) <= power_threshold2(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes(1,:) <= power_threshold1(1));


    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_low_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_low_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_low_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_low_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_low_SO{nsession}(2,nCell,:))';
    end



    %%%%%%%%%%%% SO phase
    event_phase=[];

    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        event_phase = [event_phase ripples_all(nprobe).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),ia)];
    end
    
    is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
    is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_peak_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_peak_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_peak_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_peak_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_peak_SO{nsession}(2,nCell,:))';
    end



    is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
    is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_trough_SO{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_trough_SO{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_trough_SO{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_trough_SO{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_trough_SO{nsession}(2,nCell,:))';
    end


    %%%%%%%%%%%% spindle phase
    event_phase=[];

    for nprobe = 1:2
        session_event_index = find(ripples_all(nprobe).session_count == nsession);
        [C,ia,ib] = intersect(session_event_index,find(ripples_all(nprobe).session_count == nsession&ripples_all(nprobe).SWS_index==1));
        event_phase = [event_phase ripples_all(nprobe).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),ia)];
    end
    
    is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
    is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_peak_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_peak_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_peak_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_peak_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_peak_spindle{nsession}(2,nCell,:))';
    end



    is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
    is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
    T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1);
    T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1);

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_trough_spindle{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
        context_modulation_all.PSTH_trough_spindle{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

        % Difference PSTH (stim1 - stim2)
        context_modulation_all.PSTH_diff_trough_spindle{nsession}(nCell,:) = ...
            squeeze(context_modulation_all.PSTH_trough_spindle{nsession}(1,nCell,:))' - ...
            squeeze(context_modulation_all.PSTH_trough_spindle{nsession}(2,nCell,:))';
    end


    hemi_id = [ones(1,sum((ripples_all(1).SWS_index==1 & ripples_all(1).session_count==nsession)>0)) 2*ones(1,sum((ripples_all(2).SWS_index==1 & ripples_all(2).session_count==nsession)>0))];

    %%%%%%%%%%%% SO power (different timing)
    log_odds_threshold = prctile(mean_bias,[20 80]);

    for n = 1:3
        temp = [squeeze(ripples_TF_stats.V1_amp{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_amp{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_amp{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_amp{nsession}(1,1,timebins(n),hemi_id == 2))'];
        temp_percentile = tiedrank(temp.').' / size(temp, 2);

        %%%%% LOW
        T1_index = find(mean_bias > log_odds_threshold(2) & temp_percentile(2,:) <= 0.25);
        T2_index = find(mean_bias < log_odds_threshold(1) & temp_percentile(1,:) <= 0.25);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_low_SO_time{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_low_SO_time{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_low_SO_time{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_low_SO_time{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_low_SO_time{nsession}(2,nCell,:))';
        end

        %%%%% High
        T1_index = find(mean_bias > log_odds_threshold(2) & temp_percentile(2,:) >= 0.75);
        T2_index = find(mean_bias < log_odds_threshold(1) & temp_percentile(1,:) >= 0.75);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_high_SO_time{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_high_SO_time{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_high_SO_time{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_high_SO_time{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_high_SO_time{nsession}(2,nCell,:))';
        end

    end


    %%%%%%%%%%%% SO phase with different times
    event_phase = [];
    for n = 1:3
        event_phase = [squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 2))'];
        % temp_percentile = tiedrank(temp.').' / size(temp, 2);

        is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
        is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_peak_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_peak_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_peak_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_peak_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_peak_SO_time{n}{nsession}(2,nCell,:))';
        end

   

        is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
        is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_trough_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_trough_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_trough_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_trough_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_trough_SO_time{n}{nsession}(2,nCell,:))';
        end
    end

    %%%%%%%%%%%% SO phase (Sync)

    event_phase = [];
    for n = 1:3
        event_phase = [squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 2))'];
        % temp_percentile = tiedrank(temp.').' / size(temp, 2);

        is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
        is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;

        T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1 & is_peak_phase_1 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1 & is_peak_phase_2 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_sync_peak_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_sync_peak_SO_time{n}{nsession}(2,nCell,:))';
        end

   

        is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
        is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1 & is_trough_phase_1 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1 & is_trough_phase_2 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_sync_trough_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_sync_trough_SO_time{n}{nsession}(2,nCell,:))';
        end
    end


    %%%%%%%%%%%% SO phase (Anti phase)
    event_phase = [];
    for n = 1:3
        event_phase = [squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_phase_median{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_phase_median{nsession}(1,1,timebins(n),hemi_id == 2))'];
        % temp_percentile = tiedrank(temp.').' / size(temp, 2);

        is_peak_phase_1 = event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2;
        is_peak_phase_2 =  event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2;

        T1_index = find(mean_bias > log_odds_threshold(2) & is_peak_phase_2 == 1 & is_peak_phase_1 == 0);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_peak_phase_1 == 1 & is_peak_phase_2 == 0);
        
        % mainly session 4
        if length(T1_index)==1
            T1_index = [T1_index T1_index];
        end

        if length(T2_index)==1
            T2_index = [T2_index T2_index];
        end

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_dsync_peak_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_dsync_peak_SO_time{n}{nsession}(2,nCell,:))';
        end



        is_trough_phase_1 = event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi;
        is_trough_phase_2 =  event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi;
        T1_index = find(mean_bias > log_odds_threshold(2) & is_trough_phase_2 == 1 & is_trough_phase_1 == 0);
        T2_index = find(mean_bias < log_odds_threshold(1) & is_trough_phase_1 == 1 & is_trough_phase_2 == 0);
        
        % mainly session 4
        if length(T1_index)==1
            T1_index = [T1_index T1_index];
        end

        if length(T2_index)==1
            T2_index = [T2_index T2_index];
        end

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_dsync_trough_SO_time{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_dsync_trough_SO_time{n}{nsession}(2,nCell,:))';
        end
    end

    %%%%%%%%%%%% SO V1 HPC PLV
    temp = [];
    for n = 1:3
        temp = [squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,1,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,1,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,1,timebins(n),hemi_id == 2))'];
        temp_percentile = tiedrank(temp.').' / size(temp, 2);

        SO_index1 = temp_percentile(1,:) < 0.25;
        SO_index2 = temp_percentile(2,:) < 0.25;

        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_low_V1_HPC_SO_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:))';
        end

        SO_index1 = temp_percentile(1,:) > 0.75;
        SO_index2 = temp_percentile(2,:) > 0.75;


        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_high_V1_HPC_SO_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_SO_PLV{n}{nsession}(2,nCell,:))';
        end
    end

    %%%%%%%%%%%% Spindle V1 HPC PLV

    temp = [];
    for n = 1:3
        temp = [squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,3,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,3,timebins(n),hemi_id == 2))'; ...
            squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(2,3,timebins(n),hemi_id == 1))' squeeze(ripples_TF_stats.V1_HPC_PLV{nsession}(1,3,timebins(n),hemi_id == 2))'];
        temp_percentile = tiedrank(temp.').' / size(temp, 2);

        SO_index1 = temp_percentile(1,:) < 0.25;
        SO_index2 = temp_percentile(2,:) < 0.25;

        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_low_V1_HPC_spindle_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_low_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:))';
        end

        SO_index1 = temp_percentile(1,:) > 0.75;
        SO_index2 = temp_percentile(2,:) > 0.75;


        T1_index = find(mean_bias > log_odds_threshold(2) & SO_index2 == 1);
        T2_index = find(mean_bias < log_odds_threshold(1) & SO_index1 == 1);

        for nCell = 1:length(all_clusters)
            % Ripple PSTH
            context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T1_index,:)));
            context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH_zscored(nCell,T2_index,:)));

            % Difference PSTH (stim1 - stim2)
            context_modulation_all.PSTH_diff_high_V1_HPC_spindle_PLV{n}{nsession}(nCell,:) = ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(1,nCell,:))' - ...
                squeeze(context_modulation_all.PSTH_high_V1_HPC_spindle_PLV{n}{nsession}(2,nCell,:))';
        end
    end



    context_modulation_all.ripple_modulation_percentile{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).ripple_modulation_percentile; ripple_modulation_PSTH_all{nsession}(2).ripple_modulation_percentile]);
    context_modulation_all.modulation_percentile_PRE{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).modulation_percentile_PRE; ripple_modulation_PSTH_all{nsession}(2).modulation_percentile_PRE]);
end
context_modulation_all.timebin = ripple_modulation.bins;

toc

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all_V1_bias_PRE_ripple.mat'),'context_modulation_all')

    % subplot(2,2,1)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.POST_ripple_FR(V1_id))
    % subplot(2,2,2)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
%     % 



%% Plotting context selecitve ripple modulation based on V1 bias
% scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all_V1_bias.mat'),'context_modulation_all')

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
FR_track_ratio = (FR_track(1,:)-FR_track(2,:))./( FR_track(1,:)+FR_track(2,:));






load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all_V1_bias.mat'),'context_modulation_all')

%%%%%%%%%%%%%%% Spindle power
ripple_types = {'LOW','HIGH'};
PSTH_fields = {'PSTH_diff_low_spindle','PSTH_diff_high_spindle'};
ripple_name = 'spindle power';
plot_context_selective_ripple_modulation(context_modulation_all,PSTH_fields,ripple_types,ripple_name,'bias_option', 'V1')
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])

%%%%%%%%%%%%%%% Spindle power
ripple_types = {'peak','trough'};
PSTH_fields = {'PSTH_diff_low_spindle','PSTH_diff_high_spindle'};
ripple_name = 'spindle power';
plot_context_selective_ripple_modulation(context_modulation_all,PSTH_fields,ripple_types,ripple_name,'bias_option', 'V1')
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])




load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all_V1_bias_PRE_ripple.mat'),'context_modulation_all')

%%%%%%%%%%%%%%% Spindle power
ripple_types = {'LOW','HIGH'};
PSTH_fields = {'PSTH_diff_low_spindle','PSTH_diff_high_spindle'};
ripple_name = 'spindle power';
plot_context_selective_ripple_modulation(context_modulation_all,PSTH_fields,ripple_types,ripple_name,'bias_option', 'V1')
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])




% 
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')
% 
% %%%%%%%%%%%%%%% Spindle power
% ripple_types = {'LOW','HIGH'};
% PSTH_fields = {'PSTH_diff_low_ripple','PSTH_diff_high_ripple'};
% ripple_name = 'ripple power';
% plot_context_selective_ripple_modulation(context_modulation_all,PSTH_fields,ripple_types,ripple_name,'bias_option', 'HPC')


%% 