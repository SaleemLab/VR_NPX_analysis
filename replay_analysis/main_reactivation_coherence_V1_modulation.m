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

%% high ripple V1 modulation
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

% timebins = [5 6 10]; % Timebin of the LFP metric relative to ripples where 1 is -1 to -0.8s and 10 is 0.8 to 1s

V1_reactivation_modulation_all = struct();
% context_corr_all = struct();
tic
hemispheres = {'L','R'};


%%%% high ripple
for hemi = 1:2
    for track_id = 1:2
        V1_SUA_reactivation_PSTH{hemi}{track_id}=[];
        V1_SUA_reactivation_PSTH_odd{hemi}{track_id}=[];
        V1_SUA_reactivation_PSTH_even{hemi}{track_id}=[];
        V1_SUA_track_difference{hemi}=[];
        V1_SUA_z_track_difference{hemi}=[];
        cell_session_id{hemi} = [];
    end
end

% 
for nsession = 1:length(sessions_to_process)

    for hemi = 1:2
        track_difference=[];
        z_track_difference=[];
        % cell_index

        cell_index = find((contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi})));

        % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
        %     session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1'));
        % 
        % cell_index = find((session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
        %     session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') &...
        %     contains(session_clusters_all.region{nsession},hemispheres{hemi})));

        cell_session_id{hemi} = [cell_session_id{hemi} nsession*ones(1,length(cell_index))];

        for ncell = 1:length(cell_index)
            FR_distribution = reshape([session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),1}; ...
                session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),2}], 1, []);

            z_track_difference(ncell) = (mean(mean(session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),1})) - mean(FR_distribution)) ./ std(FR_distribution) ...
                - (mean(mean(session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),2})) - mean(FR_distribution)) ./ std(FR_distribution);
        end

        V1_SUA_z_track_difference{hemi} = [V1_SUA_z_track_difference{hemi};  z_track_difference'];

        track_difference{hemi} = session_clusters_all.mean_FR{nsession}(cell_index,abs(hemi-3)) - session_clusters_all.mean_FR{nsession}(cell_index,hemi);
        V1_SUA_track_difference{hemi} = [V1_SUA_track_difference{hemi};  track_difference];
    end
end


for nsession = 1:length(sessions_to_process)

    ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    threshold = prctile(ripple_powers,25);
    % threshold = 0;
    [event_ids_first,event_ids_second] = merge_bilateral_ripple_events(event_id,event_times,0.05);
    high_ripple_index = find( ripple_powers < threshold);
    event_times = event_times(high_ripple_index);
    event_id = event_id(high_ripple_index);

    % event_ids_first = intersect(high_ripple_index,event_ids_first);
    % event_times = event_times(event_ids_first);
    % event_id = event_id(event_ids_first);

    %%%%%%%%%% V1
    V1_zPSTH=cell(1,2);
    V1_PSTH=cell(1,2);
    track_difference=cell(1,2);
    for hemi = 1:2
        % cell_index

        cell_index = (contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi}));

        % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
        %     contains(session_clusters_all.region{nsession},'V1'));

        cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
            contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi}));

        % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
        %     session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1'));

        % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
        %     session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') &...
        %     contains(session_clusters_all.region{nsession},hemispheres{hemi}));

        cell_id = session_clusters_all.spatial_cell_id{nsession}(cell_index);

        cell_session_id{hemi} = [cell_session_id{hemi} nsession*ones(1,length(cell_index))];

        for ncell = 1:length(cell_index)
            FR_distribution = reshape([session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),1}; ...
                session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),2}], 1, []);

            z_track_difference(ncell) = (mean(mean(session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),1})) - mean(FR_distribution)) ./ std(FR_distribution) ...
                - (mean(mean(session_clusters_all.spatial_response_raw{nsession}{cell_index(ncell),2})) - mean(FR_distribution)) ./ std(FR_distribution);
        end

        V1_SUA_z_track_difference{hemi} = [V1_SUA_z_track_difference{hemi};  z_track_difference'];

        track_difference{hemi} = session_clusters_all.mean_FR{nsession}(cell_index,abs(hemi-3)) - session_clusters_all.mean_FR{nsession}(cell_index,hemi);
        V1_SUA_track_difference{hemi} = [V1_SUA_track_difference{hemi};  track_difference];


        spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

        spike_id = session_clusters_all.spike_id{nsession}(spike_index);
        % spike_id(spike_id>0)=1;

        ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
            'unit_id',cell_id,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);

        % V1_PSTH{hemi} = [zscore(squeeze(ripple_modulation(1).PSTH),0,2) zscore(squeeze(ripple_modulation(2).PSTH),0,2)];
        V1_PSTH{hemi} = [ripple_modulation(1).PSTH_zscored ripple_modulation(2).PSTH_zscored];
        % V1_PSTH{hemi} = ;
        
    end

    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias
    mean_bias = mean(z_bias_KDE(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[25 75]);

    mean_bias = mean_bias(event_ids_first);

    log_odds_threshold = prctile(mean_bias,[25 75]);
    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));

    for hemi = 1:2
        % for track_id = 1:2
        V1_SUA_reactivation_PSTH{hemi}{1}=[V1_SUA_reactivation_PSTH{hemi}{1}; squeeze(mean(V1_PSTH{hemi}(:,T1_index,:),2))];
        V1_SUA_reactivation_PSTH{hemi}{2}=[V1_SUA_reactivation_PSTH{hemi}{2}; squeeze(mean(V1_PSTH{hemi}(:,T2_index,:),2))];

        V1_SUA_reactivation_PSTH_odd{hemi}{1}=[V1_SUA_reactivation_PSTH_odd{hemi}{1}; squeeze(mean(V1_PSTH{hemi}(:,T1_index(1:2:end),:),2))];
        V1_SUA_reactivation_PSTH_odd{hemi}{2}=[V1_SUA_reactivation_PSTH_odd{hemi}{2}; squeeze(mean(V1_PSTH{hemi}(:,T2_index(1:2:end),:),2))];

        V1_SUA_reactivation_PSTH_even{hemi}{1}=[V1_SUA_reactivation_PSTH_even{hemi}{1}; squeeze(mean(V1_PSTH{hemi}(:,T1_index(2:2:end),:),2))];
        V1_SUA_reactivation_PSTH_even{hemi}{2}=[V1_SUA_reactivation_PSTH_even{hemi}{2}; squeeze(mean(V1_PSTH{hemi}(:,T2_index(2:2:end),:),2))];
    end

    
    % mean_bias = mean(z_bias_V1_KDE(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    % log_odds_threshold = prctile(mean_bias,[25 75]);
    % T1_index_V1 = find(mean_bias > log_odds_threshold(2));
    % T2_index_V1 = find(mean_bias < log_odds_threshold(1));
    % 
    % T1_index = intersect(T1_index,T1_index_V1);
    % T2_index = intersect(T2_index,T2_index_V1);

    % bin_centers>0 & bin_centers<0.1;

    % z_bias_KDE

    % z_bias_V1_KDE
    % for hemi = 1:2
    % 
    %     squeeze(mean(V1_PSTH{hemi}(:,T2_index,:),2));
    % end
end
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all_low.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all_high.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% 
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');


% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_track.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_track_high.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_track_low.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');


% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_high.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_low.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),'V1_SUA_reactivation_PSTH','V1_SUA_reactivation_PSTH');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_high.mat'),'V1_SUA_reactivation_PSTH');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_low.mat'),'V1_SUA_reactivation_PSTH');

% ripple_modulation(1).bins(100)

% 
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_high.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_low.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_low.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% V1_SUA_reactivation_PSTH_low = V1_SUA_reactivation_PSTH;
% V1_SUA_reactivation_PSTH_high = V1_SUA_reactivation_PSTH;


% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_track.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');



load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all_low.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all_high.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');



%%% Colormap
% track1_color = [0 0 1]; % blue
% track2_color = [1 0 0]; % red
track1_color = [145,191,219]/256; % blue
track2_color = [252,141,89]/256; % red


n_half = 80;     % Number of steps from red to white and white to blue
n_white = 40;     % Number of extra white steps to make the white region broader


% Extra white section
white = ones(n_white, 3);

% Red to white
r2w = [ ...
    linspace(track2_color(1), white(1), n_half)', ...
    linspace(track2_color(2), white(2), n_half)', ...
    linspace(track2_color(3), white(3), n_half)' ];

% White to blue
w2b =[ ...
    linspace(white(1), track1_color(1), n_half)', ...
    linspace(white(2), track1_color(2), n_half)', ...
    linspace(white(3), track1_color(3), n_half)' ];


% Combine
red_white_blue = [r2w; white; w2b];
red_white_blue = flip(red_white_blue);
% 
%%%%%%%%%%%%%%%%%%%%%
% set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
% xlim([-0.5 0.5])
% colorbar;
% colormap((red_white_blue))
% set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
%%%%%%%%%%%%%%%%%%%%%

Fig = figure;
Fig.Name = 'All V1 response to track-selective reactivations (even by odd high)';
% Fig.Name = 'All V1 response to track-selective reactivations (even by odd)';


% Fig.Name = 'Track-selective V1 response to track-selective reactivations (even by odd)';
% Fig.Name = 'Track-selective V1 response to track-selective reactivations (0-0.2s even by odd)';
% Fig.Name = 'Track-selective V1 response to track-selective reactivations (-0.2-0s even by odd)';

% Fig.Name = 'V1 response to track-selective reactivations (even by odd)';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s even by odd)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s even by odd)';

% Fig.Name = 'V1 response to track-selective reactivations (even by odd high)';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s even by odd high)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s even by odd high)';

% Fig.Name = 'V1 response to track-selective reactivations';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
% Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]

psthBinSize = 0.01;
windows = [-1.5 1.5];
x_bins = windows(1)+ timebin/2 : timebin:windows(end)-timebin/2;
ripple_modulation(1).bins = x_bins;

psthBinSize = 0.01;
windows = [-1.5 1.5];
x_bins = windows(1)+ timebin/2 : timebin:windows(end)-timebin/2;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins >-0.2 & x_bins < 0.2;
% figure

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{1}{2}(:,bins_selected),2));

% [~,index] = max((V1_SUA_reactivation_PSTH_odd{1}{2}(:,bins_selected))');
% [~,index] = sort(index);
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{1}{1}(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{1}{2}(:,bins_selected),2));

% [~,index] = max((V1_SUA_reactivation_PSTH_odd{1}{2}(:,bins_selected))');
% [~,index] = sort(index);
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{1}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,3)
% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{2}{1}(:,bins_selected),2));

% [~,index] = max((V1_SUA_reactivation_PSTH_odd{2}{1}(:,bins_selected))');
% [~,index] = sort(index);
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{2}{1}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,4)
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{2}{1}(:,bins_selected),2));

% [~,index] = max((V1_SUA_reactivation_PSTH_odd{2}{1}(:,bins_selected))');
% [~,index] = sort(index);
% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{2}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



Fig = figure;
% Fig.Name = 'All V1 response to track-selective reactivations';
Fig.Name = 'All V1 response to track-selective reactivations (low)';
% Fig.Name = 'All V1 response to track-selective reactivations (high)';
% Fig.Name = 'Track-selective V1 response to track-selective reactivations';
% Fig.Name = 'V1 response to track-selective reactivations';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > -0.2 & x_bins < 0.2;
% figure


subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{1}(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,3)
% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2));
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{1}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,4)
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2));
% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])




Fig = figure;
% Fig.Name = 'All V1 response to track-selective reactivations';
Fig.Name = 'All V1 response to track-selective reactivations (sorted by diff)';
% Fig.Name = 'All V1 response to track-selective reactivations (high sorted by diff)';

% Fig.Name = 'Track-selective V1 response to track-selective reactivations';
% Fig.Name = 'V1 response to track-selective reactivations';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > -0.2 & x_bins < 0.2;
% figure


subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{1}(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,3)
% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{2}{2}(:,bins_selected),2));

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{1}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,4)
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{2}{2}(:,bins_selected),2));

% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



Fig = figure;
Fig.Name = 'All V1 response diff to track-selective reactivations';

% Fig.Name = 'All V1 response diff to track-selective reactivations (high)';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > 0 & x_bins < 0.2;
% figure

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{2}(index,:) - V1_SUA_reactivation_PSTH{1}{1}(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track R - Track L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2));
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{2}{2}(:,bins_selected),2));

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{1}(index,:)-V1_SUA_reactivation_PSTH{2}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track L - Track R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])

%% Track Preferring neuron ripple modulation
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')

psthBinSize = 0.01;
windows = [-1.5 1.5];
x_bins = windows(1)+ timebin/2 : timebin:windows(end)-timebin/2;
ripple_modulation(1).bins = x_bins;

PRE_windows = x_bins > -0.2 & x_bins < 0;
all_windows = x_bins > -0.2 & x_bins < 0.2;
POST_windows = x_bins > 0 & x_bins < 0.2;

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
V1_SUA_reactivation_PSTH_all = V1_SUA_reactivation_PSTH;
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all_low.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
V1_SUA_reactivation_PSTH_low = V1_SUA_reactivation_PSTH;
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_all_high.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
V1_SUA_reactivation_PSTH_high = V1_SUA_reactivation_PSTH;
% 
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% V1_SUA_reactivation_PSTH_all = V1_SUA_reactivation_PSTH;
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_low.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% V1_SUA_reactivation_PSTH_low = V1_SUA_reactivation_PSTH;
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_high.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% V1_SUA_reactivation_PSTH_high = V1_SUA_reactivation_PSTH;



%%%%%%% Mixed effect model ripple modulation
windows_all = {all_windows,PRE_windows,POST_windows};

%%% All ripples
track_difference = [V1_SUA_z_track_difference{1}; V1_SUA_z_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH_all{1}{1}(:,:); V1_SUA_reactivation_PSTH_all{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH_all{1}{2}(:,:); V1_SUA_reactivation_PSTH_all{2}{2}(:,:)];
psth_diff = psth1-psth2;

session_count = [cell_session_id{1} cell_session_id{2}];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

ripple_modulation_lme = struct();
for nwin = 1:3
    y = double(mean(psth_diff(:,windows_all{nwin}),2,'omitnan'));
    x = track_difference;

    tbl = table(zscore(x), zscore(y,[],'omitnan'),subject_id, 'VariableName',{'track_difference','ripple_difference','animal_ids'});
    lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');
    ripple_modulation_lme(nwin).window = windows_all{nwin};
    % ripple_modulation_lme(nwin).window = lme.Coefficients.Estimate(2);
    ripple_modulation_lme(nwin).b = lme.Coefficients.Estimate(2);
    ripple_modulation_lme(nwin).t_stat = lme.Coefficients.tStat(2);
    ripple_modulation_lme(nwin).p = lme.Coefficients.pValue(2);
    % ripple_modulation_lme(nwin).b_CI(1) = lme.Coefficients.Lower(2);
    % ripple_modulation_lme(nwin).b_CI(2) = lme.Coefficients.Upper(2);
    % ripple_modulation_lme(nwin).p(2) = lme.Coefficients.pValue(2);

    % Mixed effect model with shuffling
    beta_shuffled = nan(1,1000);
    t_shuffled = nan(1,1000);
    beta_boot = nan(1,1000);
    t_boot = nan(1,1000);

    parfor iBoot = 1:1000
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        % idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);
        idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);

        %Run Linear Mixed-Effects Model (bootstrapped)
        tbl = table(zscore(x(idx)), zscore(y(idx),[],'omitnan'),subject_id(idx), 'VariableName',{'track_difference','ripple_difference','animal_ids'});
        lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');
        beta_boot(iBoot) = lme.Coefficients.Estimate(2);
        t_boot(iBoot) = lme.Coefficients.tStat(2);

        %Run Linear Mixed-Effects Model (shuffled)
        idx = datasample(1:length(subject_id), length(subject_id), 'Replace', false);
        tbl = table(zscore(x), zscore(y(idx),[],'omitnan'),subject_id(idx), 'VariableName',{'track_difference','ripple_difference','animal_ids'});
        lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');

        beta_shuffled(iBoot) = lme.Coefficients.Estimate(2);
        t_shuffled(iBoot) = lme.Coefficients.tStat(2);
        toc
    end
    ripple_modulation_lme(nwin).b_shuffle = beta_shuffled;
    ripple_modulation_lme(nwin).t_stat_shuffle = t_shuffled;

    ripple_modulation_lme(nwin).b_boot = beta_boot;
    ripple_modulation_lme(nwin).t_stat_boot = t_boot;
end



%%% low ripples
track_difference = [V1_SUA_z_track_difference{1}; V1_SUA_z_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH_low{1}{1}(:,:); V1_SUA_reactivation_PSTH_low{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH_low{1}{2}(:,:); V1_SUA_reactivation_PSTH_low{2}{2}(:,:)];
psth_diff = psth1-psth2;

session_count = [cell_session_id{1} cell_session_id{2}];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

ripple_modulation_lme_low = struct();
for nwin = 1:3
    y = double(mean(psth_diff(:,windows_all{nwin}),2,'omitnan'));
    x = track_difference;

    tbl = table(zscore(x), zscore(y,[],'omitnan'),subject_id, 'VariableName',{'track_difference','ripple_difference','animal_ids'});
    lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');
    ripple_modulation_lme_low(nwin).window = windows_all{nwin};
    % ripple_modulation_lme_low(nwin).window = lme.Coefficients.Estimate(2);
    ripple_modulation_lme_low(nwin).b = lme.Coefficients.Estimate(2);
    ripple_modulation_lme_low(nwin).t_stat = lme.Coefficients.tStat(2);
    ripple_modulation_lme_low(nwin).p = lme.Coefficients.pValue(2);
    % ripple_modulation_lme(nwin).b_CI(1) = lme.Coefficients.Lower(2);
    % ripple_modulation_lme(nwin).b_CI(2) = lme.Coefficients.Upper(2);
    % ripple_modulation_lme(nwin).p(2) = lme.Coefficients.pValue(2);
% end
    % Mixed effect model with shuffling
    beta_shuffled = nan(1,1000);
    t_shuffled = nan(1,1000);
    beta_boot = nan(1,1000);
    t_boot = nan(1,1000);

    parfor iBoot = 1:1000
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        % idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);
        idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);

        %Run Linear Mixed-Effects Model (bootstrapped)
        tbl = table(zscore(x(idx)), zscore(y(idx),[],'omitnan'),subject_id(idx), 'VariableName',{'track_difference','ripple_difference','animal_ids'});
        lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');
        beta_boot(iBoot) = lme.Coefficients.Estimate(2);
        t_boot(iBoot) = lme.Coefficients.tStat(2);

        %Run Linear Mixed-Effects Model (shuffled)
        idx = datasample(1:length(subject_id), length(subject_id), 'Replace', false);
        tbl = table(zscore(x), zscore(y(idx),[],'omitnan'),subject_id(idx), 'VariableName',{'track_difference','ripple_difference','animal_ids'});
        lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');

        beta_shuffled(iBoot) = lme.Coefficients.Estimate(2);
        t_shuffled(iBoot) = lme.Coefficients.tStat(2);
        toc
    end
    ripple_modulation_lme_low(nwin).b_shuffle = beta_shuffled;
    ripple_modulation_lme_low(nwin).t_stat_shuffle = t_shuffled;

    ripple_modulation_lme_low(nwin).b_boot = beta_boot;
    ripple_modulation_lme_low(nwin).t_stat_boot = t_boot;
end



%%% high ripples
track_difference = [V1_SUA_z_track_difference{1}; V1_SUA_z_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH_high{1}{1}(:,:); V1_SUA_reactivation_PSTH_high{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH_high{1}{2}(:,:); V1_SUA_reactivation_PSTH_high{2}{2}(:,:)];
psth_diff = psth1-psth2;

session_count = [cell_session_id{1} cell_session_id{2}];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

ripple_modulation_lme_high = struct();
for nwin = 1:3
    y = double(mean(psth_diff(:,windows_all{nwin}),2,'omitnan'));
    x = track_difference;

    tbl = table(zscore(x), zscore(y,[],'omitnan'),subject_id, 'VariableName',{'track_difference','ripple_difference','animal_ids'});
    lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');
    ripple_modulation_lme_high(nwin).window = windows_all{nwin};
    % ripple_modulation_lme_high(nwin).window = lme.Coefficients.Estimate(2);
    ripple_modulation_lme_high(nwin).b = lme.Coefficients.Estimate(2);
    ripple_modulation_lme_high(nwin).t_stat = lme.Coefficients.tStat(2);
    ripple_modulation_lme_high(nwin).p = lme.Coefficients.pValue(2);
% end
    % ripple_modulation_lme(nwin).b_CI(1) = lme.Coefficients.Lower(2);
    % ripple_modulation_lme(nwin).b_CI(2) = lme.Coefficients.Upper(2);
    % ripple_modulation_lme(nwin).p(2) = lme.Coefficients.pValue(2);

    % Mixed effect model with shuffling
    beta_shuffled = nan(1,1000);
    t_shuffled = nan(1,1000);
    beta_boot = nan(1,1000);
    t_boot = nan(1,1000);

    parfor iBoot = 1:1000
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        % idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);
        idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);

        %Run Linear Mixed-Effects Model (bootstrapped)
        tbl = table(zscore(x(idx)), zscore(y(idx),[],'omitnan'),subject_id(idx), 'VariableName',{'track_difference','ripple_difference','animal_ids'});
        lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');
        beta_boot(iBoot) = lme.Coefficients.Estimate(2);
        t_boot(iBoot) = lme.Coefficients.tStat(2);

        %Run Linear Mixed-Effects Model (shuffled)
        idx = datasample(1:length(subject_id), length(subject_id), 'Replace', false);
        tbl = table(zscore(x), zscore(y(idx),[],'omitnan'),subject_id(idx), 'VariableName',{'track_difference','ripple_difference','animal_ids'});
        lme = fitlme(tbl, 'ripple_difference ~ track_difference + (1|animal_ids)');

        beta_shuffled(iBoot) = lme.Coefficients.Estimate(2);
        t_shuffled(iBoot) = lme.Coefficients.tStat(2);
        toc
    end
    ripple_modulation_lme_high(nwin).b_shuffle = beta_shuffled;
    ripple_modulation_lme_high(nwin).t_stat_shuffle = t_shuffled;

    ripple_modulation_lme_high(nwin).b_boot = beta_boot;
    ripple_modulation_lme_high(nwin).t_stat_boot = t_boot;
end



save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_ripple_modulation_all.mat'),'ripple_modulation_lme',...
    'ripple_modulation_lme_low','ripple_modulation_lme_high');
% 
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_ripple_modulation.mat'),'ripple_modulation_lme',...
%     'ripple_modulation_lme_low','ripple_modulation_lme_high');

% 
% for nwin = 1:3
%     ripple_modulation_lme(nwin).window = find(windows_all{nwin});
%     ripple_modulation_lme_low(nwin).window =find(windows_all{nwin});
%     ripple_modulation_lme_high(nwin).window = find(windows_all{nwin});
% end

%%%%  ripples
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;


% nwin = 1;
Fig = figure;
Fig.Name = 'context-selective ripple modulation glme summary';

Fig.Position = [60 60 530 930];
counter = 0;
twin_text = {'All','PRE','POST'};

for nwin = 1:3
    counter = counter+1;
    subplot(3,2,counter)

    % 1. Calculate Means
    p  = [ripple_modulation_lme(nwin).p; ripple_modulation_lme_low(nwin).p; ripple_modulation_lme_high(nwin).p]';
    beta_boot = [ripple_modulation_lme(nwin).b_boot; ripple_modulation_lme_low(nwin).b_boot; ripple_modulation_lme_high(nwin).b_boot]';
    beta_shuffled = [ripple_modulation_lme(nwin).b_shuffle; ripple_modulation_lme_low(nwin).b_shuffle; ripple_modulation_lme_high(nwin).b_shuffle]';

    mean_boot = [ripple_modulation_lme(nwin).b ripple_modulation_lme_low(nwin).b ripple_modulation_lme_high(nwin).b];
    mean_shuf = mean(beta_shuffled, 'omitnan');

    % 2. Calculate 95% Confidence Intervals (Percentile Method)
    % This calculates the distance from the mean to the 2.5th and 97.5th percentiles
    ci_boot = [mean_boot - prctile(beta_boot, 2.5); prctile(beta_boot, 97.5) - mean_boot];
    ci_shuf = [mean_shuf - prctile(beta_shuffled, 2.5); prctile(beta_shuffled, 97.5) - mean_shuf];

    % 3. Plotting
    hold on;
    x = [1, 2, 3];
    % colors = [231, 41, 138; 256/2 256/2 256/2]/256; % Blue for Boot, Gray for Shuffled
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.15;    % Distance from the center integer (half the gap between bars)

    for i = 1:3
        % Plot Bar
        bar(x(i) + group_offset, mean_boot(i),bar_width, 'FaceColor', colour_lines(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        % Plot 95% CI Error Bar
        errorbar(x(i)+ group_offset, mean_boot(i), ci_boot(1,i), ci_boot(2,i), 'k', 'LineWidth', 1.5, 'CapSize', 10);

        % Plot Bar
        bar(x(i) - group_offset, mean_shuf(i),bar_width, 'FaceColor', colour_lines(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        % Plot 95% CI Error Bar
        errorbar(x(i)- group_offset, mean_shuf(i), ci_shuf(1,i), ci_shuf(2,i), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(x(i),0.07,sprintf('p = %.3e',p(i)))
    end

    % Formatting
    xticks([1 2 3]);
    xticklabels({'All', 'low','high'});
    ylabel('Standardized Beta');
    title('Mean Beta with 95% CI');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



    counter = counter+1;
    subplot(3,2,counter)

    % 1. Calculate Means
    p  = [ripple_modulation_lme(nwin).p; ripple_modulation_lme_low(nwin).p; ripple_modulation_lme_high(nwin).p]';
    beta_boot = [ripple_modulation_lme(nwin).t_stat_boot; ripple_modulation_lme_low(nwin).t_stat_boot; ripple_modulation_lme_high(nwin).t_stat_boot]';
    beta_shuffled = [ripple_modulation_lme(nwin).t_stat_shuffle; ripple_modulation_lme_low(nwin).t_stat_shuffle; ripple_modulation_lme_high(nwin).t_stat_shuffle]';

    mean_boot = [ripple_modulation_lme(nwin).t_stat ripple_modulation_lme_low(nwin).t_stat ripple_modulation_lme_high(nwin).t_stat];
    mean_shuf = mean(beta_shuffled, 'omitnan');

    % 2. Calculate 95% Confidence Intervals (Percentile Method)
    % This calculates the distance from the mean to the 2.5th and 97.5th percentiles
    ci_boot = [mean_boot - prctile(beta_boot, 2.5); prctile(beta_boot, 97.5) - mean_boot];
    ci_shuf = [mean_shuf - prctile(beta_shuffled, 2.5); prctile(beta_shuffled, 97.5) - mean_shuf];

    % 3. Plotting
    hold on;
    x = [1, 2, 3];
    % colors = [231, 41, 138; 256/2 256/2 256/2]/256; % Blue for Boot, Gray for Shuffled
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.15;    % Distance from the center integer (half the gap between bars)

    for i = 1:3
        % Plot Bar
        bar(x(i) + group_offset, mean_boot(i),bar_width, 'FaceColor', colour_lines(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        % Plot 95% CI Error Bar
        errorbar(x(i)+ group_offset, mean_boot(i), ci_boot(1,i), ci_boot(2,i), 'k', 'LineWidth', 1.5, 'CapSize', 10);

        % Plot Bar
        bar(x(i) - group_offset, mean_shuf(i),bar_width, 'FaceColor', colour_lines(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        % Plot 95% CI Error Bar
        errorbar(x(i)- group_offset, mean_shuf(i), ci_shuf(1,i), ci_shuf(2,i), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(x(i),0.07,sprintf('p = %.3e',p(i)))
    end

    % Formatting
    xticks([1 2 3]);
    xticklabels({'All', 'low','high'});
    ylabel('t stat');
    % title('Mean t stat with 95% CI');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])


%%%%% Linear regression

% nwin = 1;
Fig = figure;
Fig.Name = 'context-selective ripple modulation glme scatter';

Fig.Position = [60 60 1270 930];
counter = 0;
twin_text = {'All','PRE','POST'};
%%% All ripples
track_difference = [V1_SUA_z_track_difference{1}; V1_SUA_z_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH_all{1}{1}(:,:); V1_SUA_reactivation_PSTH_all{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH_all{1}{2}(:,:); V1_SUA_reactivation_PSTH_all{2}{2}(:,:)];
psth_diff = psth1-psth2;

session_count = [cell_session_id{1} cell_session_id{2}];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

for nwin = 1:3
    p  = ripple_modulation_lme(nwin).p;
    nexttile
    y = double(mean(psth_diff(:,windows_all{nwin}),2,'omitnan'));

    x = track_difference(~isnan(y));
    y(isnan(y))=[];

    fitted = polyfit(x, y, 1);

    % 3. Calculate fitted values
    y_fit = polyval(fitted, x);

    scatter(x,y,20,'k',...
        'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
    hold on;plot(x,y_fit,'r-','LineWidth',2)
    ylim([-0.2 0.2]);
    xlim([-2.5 2.5])
    xlabel('Track difference (z)')
    ylabel('Ripple difference (z)')
    text(-2,0.1,sprintf('p = %.3e',p))
    title([twin_text{nwin},' all'])
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end


%%% low ripples
track_difference = [V1_SUA_z_track_difference{1}; V1_SUA_z_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH_low{1}{1}(:,:); V1_SUA_reactivation_PSTH_low{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH_low{1}{2}(:,:); V1_SUA_reactivation_PSTH_low{2}{2}(:,:)];
psth_diff = psth1-psth2;

session_count = [cell_session_id{1} cell_session_id{2}];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

for nwin = 1:3
    p  = ripple_modulation_lme_low(nwin).p;
    nexttile
    y = double(mean(psth_diff(:,windows_all{nwin}),2,'omitnan'));

    x = track_difference(~isnan(y));
    y(isnan(y))=[];

    fitted = polyfit(x, y, 1);

    % 3. Calculate fitted values
    y_fit = polyval(fitted, x);

    scatter(x,y,20,'k',...
        'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
    hold on;plot(x,y_fit,'r-','LineWidth',2)
    ylim([-0.2 0.2]);
    xlim([-2.5 2.5])
    xlabel('Track difference (z)')
    ylabel('Ripple difference (z)')
    text(-2,0.1,sprintf('p = %.3e',p))
    title([twin_text{nwin},' low'])
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end


%%% high ripples
track_difference = [V1_SUA_z_track_difference{1}; V1_SUA_z_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH_high{1}{1}(:,:); V1_SUA_reactivation_PSTH_high{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH_high{1}{2}(:,:); V1_SUA_reactivation_PSTH_high{2}{2}(:,:)];
psth_diff = psth1-psth2;

session_count = [cell_session_id{1} cell_session_id{2}];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

for nwin = 1:3
    p  = ripple_modulation_lme_high(nwin).p;
    nexttile
    y = double(mean(psth_diff(:,windows_all{nwin}),2,'omitnan'));

    x = track_difference(~isnan(y));
    y(isnan(y))=[];

    fitted = polyfit(x, y, 1);

    % 3. Calculate fitted values
    y_fit = polyval(fitted, x);

    scatter(x,y,20,'k',...
        'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
    hold on;plot(x,y_fit,'r-','LineWidth',2)
    ylim([-0.2 0.2]);
    xlim([-2.5 2.5])
    xlabel('Track difference (z)')
    ylabel('Ripple difference (z)')
    text(-2,0.1,sprintf('p = %.3e',p))
    title([twin_text{nwin},' high'])
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



%%%%%%%%%%%%%%%%%%%%%% PSTH visualisation ripple
Fig = figure;
% Fig.Name = 'All track preferring V1 response to track-selective reactivations';
% Fig.Name = 'All track preferring V1 response to track-selective reactivations (high)';
Fig.Name = 'All track preferring V1 response to track-selective reactivations (low)';

% Fig.Name = 'All V1 response diff to track-selective reactivations (high)';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > 0 & x_bins < 0.2;
% figure

track_difference = [-V1_SUA_track_difference{1}; V1_SUA_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH{1}{1}(:,:); V1_SUA_reactivation_PSTH{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH{1}{2}(:,:); V1_SUA_reactivation_PSTH{2}{2}(:,:)];

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));
% [~,index] = sort(track_difference);
index = find(track_difference>0);

[~,index1] = sort(nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth1(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference>0);

[~,index1] = sort(nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth2(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(1,4,3)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

[~,index1] = sort(nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth1(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



subplot(1,4,4)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

[~,index1] = sort(nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth2(index,:));


% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])





Fig = figure;
% Fig.Name = 'All track preferring V1 response to track-selective reactivations (sorted by diff)';
% Fig.Name = 'All track preferring V1 response to track-selective reactivations (sorted by diff low)';
Fig.Name = 'All track preferring V1 response to track-selective reactivations (sorted by diff high)';
% Fig.Name = 'All V1 response diff to track-selective reactivations (high)';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > 0 & x_bins < 0.2;
% figure

% cell_session_id

session_id = [cell_session_id{1} cell_session_id{2}];
track_difference = [-V1_SUA_track_difference{1}; V1_SUA_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH{1}{1}(:,:); V1_SUA_reactivation_PSTH{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH{1}{2}(:,:); V1_SUA_reactivation_PSTH{2}{2}(:,:)];

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));
% [~,index] = sort(track_difference);
index = find(track_difference>0);

[~,index1] = sort(nanmean(psth1(index,bins_selected),2)-nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth1(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference>0);

[~,index1] = sort(nanmean(psth1(index,bins_selected),2)-nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth2(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(1,4,3)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

[~,index1] = sort(nanmean(psth2(index,bins_selected),2)-nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth1(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



subplot(1,4,4)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

[~,index1] = sort(nanmean(psth2(index,bins_selected),2)-nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth2(index,:));


% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



Fig = figure;
% Fig.Name = 'All track preferring V1 diff response to track-selective reactivations';
Fig.Name = 'All track preferring V1 diff response to track-selective reactivations (low)';
% Fig.Name = 'All track preferring V1 diff response to track-selective reactivations (high)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > 0 & x_bins < 0.2;
% figure

track_difference = [-V1_SUA_track_difference{1}; V1_SUA_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH{1}{1}(:,:); V1_SUA_reactivation_PSTH{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH{1}{2}(:,:); V1_SUA_reactivation_PSTH{2}{2}(:,:)];

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));
% [~,index] = sort(track_difference);
index = find(track_difference>0);

[~,index1] = sort(nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth1(index,:)-psth2(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track L-R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

[~,index1] = sort(nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth2(index,:)-psth1(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track R-L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])


Fig = figure;
Fig.Name = 'All track preferring V1 diff response to track-selective reactivations (sorted by diff)';
Fig.Name = 'All track preferring V1 diff response to track-selective reactivations (sorted by diff low)';
Fig.Name = 'All track preferring V1 diff response to track-selective reactivations (sorted by diff high)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > 0 & x_bins < 0.2;
% figure

track_difference = [-V1_SUA_track_difference{1}; V1_SUA_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH{1}{1}(:,:); V1_SUA_reactivation_PSTH{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH{1}{2}(:,:); V1_SUA_reactivation_PSTH{2}{2}(:,:)];

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));
% [~,index] = sort(track_difference);
index = find(track_difference>0);

[~,index1] = sort(nanmean(psth1(index,bins_selected),2)-nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth1(index,:)-psth2(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track L-R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

[~,index1] = sort(nanmean(psth2(index,bins_selected),2)-nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,psth2(index,:)-psth1(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.2 0.2])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track R-L ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Colormap
% track1_color = [0 0 1]; % blue
% track2_color = [1 0 0]; % red
track1_color = [145,191,219]/256; % blue
track2_color = [252,141,89]/256; % red
white = ones(n_white, 3);

n_half = 80;     % Number of steps from red to white and white to blue
n_white = 40;     % Number of extra white steps to make the white region broader

% Red to white
r2w = [ ...
    linspace(track2_color(1), white(1), n_half)', ...
    linspace(track2_color(2), white(2), n_half)', ...
    linspace(track2_color(3), white(3), n_half)' ];

% White to blue
w2b =[ ...
    linspace(white(1), track1_color(1), n_half)', ...
    linspace(white(2), track1_color(2), n_half)', ...
    linspace(white(3), track1_color(3), n_half)' ];

% Combine
red_white_blue = [r2w; white; w2b];
red_white_blue = flip(red_white_blue);

Fig = figure;
Fig.Name = 'All track preferring V1 diff response to L-R track-selective reactivations (high smoothed)';
Fig.Position = [60 60 1280 406];
% Fig.Position = [60 60 1876 1023];
Fig.Position = [60 60 1296 953];


% [595 645 1296 406]
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins > 0 & x_bins < 0.2;
% figure

track_difference = [-V1_SUA_track_difference{1}; V1_SUA_track_difference{2}];% Track L - Track R;
psth1 = [V1_SUA_reactivation_PSTH{1}{1}(:,:); V1_SUA_reactivation_PSTH{2}{1}(:,:)];
psth2 = [V1_SUA_reactivation_PSTH{1}{2}(:,:); V1_SUA_reactivation_PSTH{2}{2}(:,:)];

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2)-nanmean(V1_SUA_reactivation_PSTH{1}{1}(:,bins_selected),2));
% [~,index] = sort(track_difference);
index = find(track_difference>0);

% [~,index1] = sort(nanmean(psth1(index,bins_selected),2)-nanmean(psth2(index,bins_selected),2));
[~,index1] = sort(nanmean(psth1(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,movmedian(psth1(index,:)-psth2(index,:),10,2,'omitnan'));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.15 0.15])
xlim([-0.5 0.5])
xline(0,'--')
title('All L prefering V1 responding to Track L-R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
index = find(track_difference<0);

% [~,index1] = sort(nanmean(psth2(index,bins_selected),2)-nanmean(psth1(index,bins_selected),2));
[~,index1] = sort(nanmean(psth2(index,bins_selected),2));
index = index(index1);

x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
% h = imagesc(x_data,y_data,psth2(index,:)-psth1(index,:));
h = imagesc(x_data,y_data,movmedian(psth1(index,:)-psth2(index,:),10,2,'omitnan'));
% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.15 0.15])
xlim([-0.5 0.5])
xline(0,'--')
title('All R prefering V1 responding to Track L-R ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])








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

% timebins = [5 6 10]; % Timebin of the LFP metric relative to ripples where 1 is -1 to -0.8s and 10 is 0.8 to 1s

% V1_reactivation_modulation_all = struct();
% context_corr_all = struct();
tic
hemispheres = {'L','R'};

clear V1_reactivation_coherence_modulation
V1_SUA_track_difference=cell(1,2);
%%%%  ripple
for track_id = 1:2
    V1_reactivation_coherence_modulation(track_id).event_type = [];

    V1_reactivation_coherence_modulation(track_id).event_type_PRE = [];

    V1_reactivation_coherence_modulation(track_id).session_id= [];

    V1_reactivation_coherence_modulation(track_id).V1_bias = [];

    V1_reactivation_coherence_modulation(track_id).V1_bias_PRE = [];

    V1_reactivation_coherence_modulation(track_id).neurons_active = [];

    V1_reactivation_coherence_modulation(track_id).neurons_active_PRE = [];

    V1_reactivation_coherence_modulation(track_id).PSTH1 = [];
    V1_reactivation_coherence_modulation(track_id).PSTH2 = [];
    V1_reactivation_coherence_modulation(track_id).HPC_PSTH = [];

    V1_reactivation_coherence_modulation(track_id).PSTH1_zscored = [];
    V1_reactivation_coherence_modulation(track_id).PSTH2_zscored = [];
    V1_reactivation_coherence_modulation(track_id).HPC_PSTH_zscored = [];

    % V1_reactivation_coherence_modulation(track_id).PSTH2_PRE = [];

    V1_reactivation_coherence_modulation(track_id).event_id = [];

end


counter = 1;
for nsession = 1:length(sessions_to_process)
    tic
    ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    threshold = prctile(ripple_powers,75);
    [event_ids_first,event_ids_second] = merge_bilateral_ripple_events(event_id,event_times,0.05);
    % high_ripple_index = find( ripple_powers > threshold);
    % event_ids_first = intersect(high_ripple_index,event_ids_first);


    event_times = event_times(event_ids_first);
    event_id = event_id(event_ids_first);

    %%%%%%%%%% V1
    V1_zPSTH=cell(1,2);
    V1_PSTH=cell(1,2);
    HPC_PSTH=cell(1,1);
    HPC_zPSTH=cell(1,1);
    track_difference=cell(1,2);


    for hemi = 1:2
        % cell_index

        % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
        %     contains(session_clusters_all.region{nsession},'V1') &...
        %     contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        cell_index = (contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi}));

        cell_id = session_clusters_all.spatial_cell_id{nsession}(cell_index);

        track_difference{hemi} = session_clusters_all.mean_FR{nsession}(cell_index,abs(hemi-3)) - session_clusters_all.mean_FR{nsession}(cell_index,hemi);
        V1_SUA_track_difference{hemi} = [V1_SUA_track_difference{hemi};  track_difference{hemi}];

        spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

        spike_id = session_clusters_all.spike_id{nsession}(spike_index);
        % spike_id(spike_id>0)=1;

        ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
            'unit_id',cell_id,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0,'smooth_option',0);

        % V1_PSTH{hemi} = [zscore(squeeze(ripple_modulation(1).PSTH),0,2) zscore(squeeze(ripple_modulation(2).PSTH),0,2)];
        V1_PSTH{hemi} = [squeeze(ripple_modulation(1).PSTH).*psthBinSize squeeze(ripple_modulation(2).PSTH).*psthBinSize];
        V1_zPSTH{hemi} = [squeeze(ripple_modulation(1).PSTH_zscored) squeeze(ripple_modulation(2).PSTH_zscored)];

    end


    cell_index = (contains(session_clusters_all.region{nsession},'HPC'));
    % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') &...
    %     contains(session_clusters_all.region{nsession},hemispheres{hemi}));
    cell_id = session_clusters_all.spatial_cell_id{nsession}(cell_index);

    spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);
    spike_id = session_clusters_all.spike_id{nsession}(spike_index);
    % spike_id(spike_id>0)=1;

    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
        'unit_id',cell_id,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0,'smooth_option',0);

    % V1_PSTH{hemi} = [zscore(squeeze(ripple_modulation(1).PSTH),0,2) zscore(squeeze(ripple_modulation(2).PSTH),0,2)];
    HPC_PSTH = [squeeze(ripple_modulation(1).PSTH).*psthBinSize squeeze(ripple_modulation(2).PSTH).*psthBinSize];
    HPC_zPSTH = [squeeze(ripple_modulation(1).PSTH_zscored) squeeze(ripple_modulation(2).PSTH_zscored)];


    xbins = ripple_modulation(1).bins;
    % tvec = ripple_modulation(1).;
    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias

    HPC_bias_this_session = z_bias_KDE(:,session_id == sessions_to_process(nsession));
    HPC_bias_this_session = HPC_bias_this_session(:,event_ids_first);

    mean_bias = mean(z_bias_KDE(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    % log_odds_threshold = prctile(mean_bias,[25 75]);

    mean_bias = mean_bias(event_ids_first);

    log_odds_threshold = prctile(mean_bias,[25 75]);
    event_index = [];
    event_index{1} = find(mean_bias > log_odds_threshold(2));
    event_index{2} = find(mean_bias < log_odds_threshold(1));

    % log_odds_threshold = prctile(mean_bias,[25 50 75]);
    % event_index = [];
    % event_index{1} = find(mean_bias < log_odds_threshold(3)&mean_bias > log_odds_threshold(2));
    % event_index{2} = find(mean_bias > log_odds_threshold(1) & mean_bias < log_odds_threshold(2));


    V1_bias_this_session = z_bias_V1_KDE(:,session_id == sessions_to_process(nsession));
    V1_bias_this_session = V1_bias_this_session(:,event_ids_first);

    mean_V1_bias = mean(z_bias_V1_KDE(bin_centers>0 & bin_centers<0.2,session_id == sessions_to_process(nsession)),1,'omitnan');
    mean_V1_bias = mean_V1_bias(event_ids_first);

    mean_V1_bias_PRE = mean(z_bias_V1_KDE(bin_centers>-0.2 & bin_centers<0,session_id == sessions_to_process(nsession)),1,'omitnan');
    mean_V1_bias_PRE = mean_V1_bias_PRE(event_ids_first);

    log_odds_threshold = prctile(abs(mean_V1_bias),0:20:99.5);


    event_id = find(session_id == sessions_to_process(nsession));
    event_id = event_id(event_ids_first);
    event_index_all=[];
    event_index_all{1} = event_id(event_index{1})';
    event_index_all{2} = event_id(event_index{2})';

    % T1_index = intersect(T1_index,T1_index_V1);
    % T2_index = intersect(T2_index,T2_index_V1);

    for track_id = 1:2
        neurons_active_PRE=[];
        neurons_active=[];
        neurons_FR_PRE=[];
        neurons_FR=[];

        event_type = [];
        event_type_PRE = [];
        event_type_V1 = [];

        bias_PRE = [];
        bias = [];

        PSTH_all = [];
        PSTH_z_all = [];
        PSTH1_all = [];
        PSTH1_z_all = [];
        PSTH2_all = [];
        PSTH2_z_all = [];

        for nevent = 1:length(event_index{track_id})
            bias_PRE(nevent) = mean(V1_bias_this_session(bin_centers<0&bin_centers>-0.2,event_index{track_id}(nevent)),'omitnan');
            bias(nevent) = mean(V1_bias_this_session(bin_centers>0&bin_centers<0.2,event_index{track_id}(nevent)),'omitnan');

            PSTH1 = squeeze(V1_PSTH{1}(:,event_index{track_id}(nevent),xbins>0&xbins<0.2));
            PSTH_PRE1 = squeeze(V1_PSTH{1}(:,event_index{track_id}(nevent),xbins>-0.2&xbins<0));

            PSTH2 = squeeze(V1_PSTH{2}(:,event_index{track_id}(nevent),xbins>0&xbins<0.2));
            PSTH_PRE2 = squeeze(V1_PSTH{2}(:,event_index{track_id}(nevent),xbins>-0.2&xbins<0));

            neurons_active(1,nevent) = sum(sum(PSTH1,2)>0)./size(PSTH1,1);
            neurons_active_PRE(1,nevent) = sum(sum(PSTH_PRE1,2)>0)./size(PSTH_PRE1,1);
            neurons_active(2,nevent) = sum(sum(PSTH2,2)>0)./size(PSTH2,1);
            neurons_active_PRE(2,nevent) = sum(sum(PSTH_PRE2,2)>0)./size(PSTH_PRE2,1);


            PSTH1 = squeeze(V1_PSTH{1}(:,event_index{track_id}(nevent),xbins>=-1&xbins<=1));
            % PSTH_PRE1 = squeeze(V1_PSTH{1}(:,event_index{track_id}(nevent),xbins>-0.2&xbins<0));

            PSTH2 = squeeze(V1_PSTH{2}(:,event_index{track_id}(nevent),xbins>=-1&xbins<=1));
            % PSTH_PRE2 = squeeze(V1_PSTH{2}(:,event_index{track_id}(nevent),xbins>-0.2&xbins<0));
            PSTH = squeeze(HPC_PSTH(:,event_index{track_id}(nevent),xbins>=-1&xbins<=1));

            PSTH1_all = [PSTH1_all {PSTH1}];
            PSTH2_all = [PSTH2_all {PSTH2}];
            PSTH_all = [PSTH_all {PSTH}];


            PSTH1 = squeeze(V1_zPSTH{1}(:,event_index{track_id}(nevent),xbins>=-1&xbins<=1));
            PSTH2 = squeeze(V1_zPSTH{2}(:,event_index{track_id}(nevent),xbins>=-1&xbins<=1));
            PSTH = squeeze(HPC_zPSTH(:,event_index{track_id}(nevent),xbins>=-1&xbins<=1));

            PSTH1_z_all = [PSTH1_z_all {PSTH1}];
            PSTH2_z_all = [PSTH2_z_all {PSTH2}];
            PSTH_z_all = [PSTH_z_all {PSTH}];
        end


        for nthreshold = 1:length(log_odds_threshold)-1
            if track_id == 1
                V1_index = find(mean_V1_bias > log_odds_threshold(nthreshold));
            else
                V1_index = find(mean_V1_bias < -log_odds_threshold(nthreshold+1));
            end

            if track_id == 1
                V1_index_PRE = find(mean_V1_bias_PRE > log_odds_threshold(nthreshold));
            else
                V1_index_PRE = find(mean_V1_bias_PRE < -log_odds_threshold(nthreshold+1));
            end


            for nevent = 1:length(event_index{track_id})
                if sum(ismember(V1_index,event_index{track_id}(nevent)))==1 % Coherent events
                    event_type(nthreshold,nevent) = 1;
                else
                    event_type(nthreshold,nevent) = 0;
                end


                if sum(ismember(V1_index_PRE,event_index{track_id}(nevent)))==1 % Coherent events
                    event_type_PRE(nthreshold,nevent) = 1;
                else
                    event_type_PRE(nthreshold,nevent) = 0;
                end



                % if sum(ismember(V1_index_PRE,event_index{track_id}(nevent)))==1 % Coherent events
                %     event_type_PRE(nthreshold,nevent) = 1;
                % else
                %     event_type_PRE(nthreshold,nevent) = 0;
                % end
            end
        end

        % nevent = nevent+1;sum(ismember(V1_index,event_index{track_id}(nevent)))==1

        % figure
        % subplot(3,1,1)
        % imagesc(xbins(xbins>-0.3&xbins<0.3),1:size(V1_PSTH{2},1),squeeze(V1_PSTH{2}(:,event_index{track_id}(nevent),xbins>-0.3&xbins<0.3)));
        % xlim([-0.3 0.3])
        % subplot(3,1,2)
        % imagesc(xbins(xbins>-0.3&xbins<0.3),1:size(V1_PSTH{1},1),squeeze(V1_PSTH{1}(:,event_index{track_id}(nevent),xbins>-0.3&xbins<0.3)));
        % xlim([-0.3 0.3])
        % subplot(3,1,3)
        % plot(bin_centers(bin_centers>-0.3&bin_centers<0.3),HPC_bias_this_session(bin_centers>-0.3&bin_centers<0.3,event_index{track_id}(nevent)));
        % hold on;plot(bin_centers(bin_centers>-0.3&bin_centers<0.3),V1_bias_this_session(bin_centers>-0.3&bin_centers<0.3,event_index{track_id}(nevent)))
        % xlim([-0.3 0.3])

        % V1_reactivation_coherence_modulation(track_id).event_type_V1 = [V1_reactivation_coherence_modulation(track_id).event_type_V1 event_type_V1];

        % V1_reactivation_coherence_modulation(track_id).mean_V1_bias_PRE = [V1_reactivation_coherence_modulation(track_id).mean_V1_bias_PRE mean_V1_bias_PRE];

        V1_reactivation_coherence_modulation(track_id).event_type = [V1_reactivation_coherence_modulation(track_id).event_type event_type];

        V1_reactivation_coherence_modulation(track_id).event_type_PRE = [V1_reactivation_coherence_modulation(track_id).event_type_PRE event_type_PRE];

        V1_reactivation_coherence_modulation(track_id).session_id = [V1_reactivation_coherence_modulation(track_id).session_id nsession*ones(size(bias))];

        V1_reactivation_coherence_modulation(track_id).V1_bias = [V1_reactivation_coherence_modulation(track_id).V1_bias bias];

        V1_reactivation_coherence_modulation(track_id).V1_bias_PRE = [V1_reactivation_coherence_modulation(track_id).V1_bias_PRE bias_PRE];

        V1_reactivation_coherence_modulation(track_id).neurons_active = [V1_reactivation_coherence_modulation(track_id).neurons_active neurons_active];

        V1_reactivation_coherence_modulation(track_id).neurons_active_PRE = [V1_reactivation_coherence_modulation(track_id).neurons_active_PRE neurons_active_PRE];


        V1_reactivation_coherence_modulation(track_id).PSTH1 = [V1_reactivation_coherence_modulation(track_id).PSTH1 PSTH1_all];
        V1_reactivation_coherence_modulation(track_id).PSTH2 = [V1_reactivation_coherence_modulation(track_id).PSTH2 PSTH2_all];
        V1_reactivation_coherence_modulation(track_id).HPC_PSTH = [V1_reactivation_coherence_modulation(track_id).HPC_PSTH PSTH_all];

        V1_reactivation_coherence_modulation(track_id).PSTH1_zscored = [V1_reactivation_coherence_modulation(track_id).PSTH1_zscored PSTH1_z_all];
        V1_reactivation_coherence_modulation(track_id).PSTH2_zscored = [V1_reactivation_coherence_modulation(track_id).PSTH2_zscored PSTH2_z_all];
        V1_reactivation_coherence_modulation(track_id).HPC_PSTH_zscored = [V1_reactivation_coherence_modulation(track_id).HPC_PSTH_zscored PSTH_z_all];

        % V1_reactivation_coherence_modulation(track_id).PSTH2_PRE = [V1_reactivation_coherence_modulation(track_id).PSTH2_PRE PSTH2_PRE_all];

        V1_reactivation_coherence_modulation(track_id).event_id = [V1_reactivation_coherence_modulation(track_id).event_id  event_index_all{track_id}];
    end
    % V1_reactivation_coherence_modulation(track_id).event_id = [V1_reactivation_coherence_modulation(track_id).event_id  event_index_all{track_id}];

    % for hemi = 1:2
    %     % for track_id = 1:2
    %     V1_SUA_reactivation_PSTH{hemi}{1}=[V1_SUA_reactivation_PSTH{hemi}{1}; squeeze(mean(V1_PSTH{hemi}(:,T1_index,:),2))];
    %     V1_SUA_reactivation_PSTH{hemi}{2}=[V1_SUA_reactivation_PSTH{hemi}{2}; squeeze(mean(V1_PSTH{hemi}(:,T2_index,:),2))];
    % end
    %
    %
    toc
end

V1_reactivation_coherence_modulation(1).time_bins =  xbins(xbins>=-1&xbins<=1);
V1_reactivation_coherence_modulation(2).time_bins =  xbins(xbins>=-1&xbins<=1);

% V1_reactivation_coherence_modulation_backup = V1_reactivation_coherence_modulation;
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'),'V1_reactivation_coherence_modulation','-v7.3');

% V1_reactivation_coherence_modulation_low = V1_reactivation_coherence_modulation;
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation_low.mat'),'V1_reactivation_coherence_modulation_low','-v7.3');

%% Neurons active vs neurons active before
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'),'V1_reactivation_coherence_modulation');
% V1_reactivation_coherence_modulation1 = V1_reactivation_coherence_modulation;

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation_low.mat'),'V1_reactivation_coherence_modulation_low');
V1_reactivation_coherence_modulation = V1_reactivation_coherence_modulation_low;

tvec = V1_reactivation_coherence_modulation(1).time_bins;
hemispheres = {'L','R'};

% V1_reactivation_modulation_all

% Clear the struct first to ensure a clean slate
% clear V1_reactivation_coherence_modulation

for hemi = 1:2

    for track_id = 1:2
        
        % 1. Initialize Neuron Containers
        V1_reactivation_coherence_modulation(track_id).POST_neuron{hemi} = [];
        V1_reactivation_coherence_modulation(track_id).PRE_neuron{hemi} = [];
        
        % 2. Initialize Correlation & Similarity Metrics
        V1_reactivation_coherence_modulation(track_id).PRE_POST_corr{hemi} = [];
        V1_reactivation_coherence_modulation(track_id).PRE_POST_cos{hemi} = [];
        
        % 3. Initialize Set Metrics (Jaccard, etc.)
        V1_reactivation_coherence_modulation(track_id).jaccard_scores{hemi} = [];
        V1_reactivation_coherence_modulation(track_id).persistence_scores{hemi} = [];
        V1_reactivation_coherence_modulation(track_id).recruitment_scores{hemi} = [];
        
        % 4. Initialize Neuron session ID
        % Note: I put this inside the loop so it is stored per track/hemi
        V1_reactivation_coherence_modulation(track_id).neuron_session_ID{hemi} = [];
    end

end


for nsession = 1:max(V1_reactivation_coherence_modulation(1).session_id)
    tic
    track_difference=[];
    % for track_id = 1:2
    % CCA(track_id).V1{nsession} = [];
    % CCA(track_id).HPC{nsession} = [];
    % CCA(track_id).corr_coef{nsession} = [];

    for hemi = 1:2
        % hemi = 1;
        cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
            contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        cell_id = session_clusters_all.spatial_cell_id{nsession}(cell_index);
        track_difference{hemi} = session_clusters_all.mean_FR{nsession}(cell_index,abs(hemi-3)) - session_clusters_all.mean_FR{nsession}(cell_index,hemi);
    end

    %%%% L V1
    for track_id= 1:2
        selected_event_id = find(V1_reactivation_coherence_modulation(track_id).session_id == nsession);

        % Prepare V1 PSTH (Concatenate Hemispheres)
        psth1_cells = {V1_reactivation_coherence_modulation(track_id).PSTH1{selected_event_id}};
        psth2_cells = {V1_reactivation_coherence_modulation(track_id).PSTH2{selected_event_id}};
        PSTH = cellfun(@(a, b) [a; b], psth1_cells, psth2_cells, 'UniformOutput', false);
        PSTH = cat(3, PSTH{:});
        PSTH = permute(PSTH, [3, 1, 2]); % Result: [Events x Neurons x Time]

        % Prepare HPC PSTH
        HPC_PSTH = {V1_reactivation_coherence_modulation(track_id).HPC_PSTH{selected_event_id}};
        HPC_PSTH = cat(3, HPC_PSTH{:});
        HPC_PSTH = permute(HPC_PSTH, [3, 1, 2]);

        hemi_id = [ones(1,size(V1_reactivation_coherence_modulation(track_id).PSTH1{selected_event_id(1)},1))...
            2*ones(1,size(V1_reactivation_coherence_modulation(track_id).PSTH2{selected_event_id(1)},1))];
        
        pre_PSTH = PSTH(:,:,tvec>-0.2 & tvec<0);
        post_PSTH = PSTH(:,:,tvec>-0 & tvec<0.2);

        active_cells = (sum(sum(pre_PSTH,3)) + sum(sum(post_PSTH,3)) > 10);

        for hemi = 1:2
            % Save Session IDs for the neurons in this hemisphere
            V1_reactivation_coherence_modulation(track_id).neuron_session_ID{hemi} = ...
                [V1_reactivation_coherence_modulation(track_id).neuron_session_ID{hemi}, nsession*ones(1,sum(hemi_id==hemi))];

            % Save POST spike counts
            V1_reactivation_coherence_modulation(track_id).POST_neuron{hemi} = ...
                [V1_reactivation_coherence_modulation(track_id).POST_neuron{hemi}, sum(sum(post_PSTH(:,hemi_id==hemi,:),3))];

            % Save PRE spike counts
            V1_reactivation_coherence_modulation(track_id).PRE_neuron{hemi} = ...
                [V1_reactivation_coherence_modulation(track_id).PRE_neuron{hemi}, sum(sum(pre_PSTH(:,hemi_id==hemi,:),3))];
        end



        for nevent = 1:size(PSTH,1)
            pre_PSTH = sum(squeeze(PSTH(nevent,:,tvec>-0.2 & tvec<0)),2);
            post_PSTH = sum(squeeze(PSTH(nevent,:,tvec>-0 & tvec<0.2)),2);
            % pre_PSTH = sum(squeeze(PSTH(nevent,:,tvec>-1 & tvec<-0.8)),2);
            % post_PSTH = sum(squeeze(PSTH(nevent,:,tvec>-0.8 & tvec<-0.6)),2);

            % active_cells = find(sum(pre_PSTH,2)+sum(post_PSTH,2)>0);
            for hemi = 1:2
                POST_active = find(sum(post_PSTH(hemi_id==hemi,:),2));
                PRE_active = find(sum(pre_PSTH(hemi_id==hemi,:),2));
                PRE_POST_active = intersect(PRE_active,POST_active);
                PRE_only_active = PRE_active(~ismember(PRE_active,POST_active));
                POST_only_active = POST_active(~ismember(POST_active,PRE_active));

                active_cells = find(sum(squeeze(PSTH(nevent,hemi_id==hemi,:)),2)>0);% Cells that fired during a 2s windows.

                %%%% Jaccard Index (The Jaccard Index measures the overlap between the two sets)
                denom = length(PRE_POST_active) + length(PRE_only_active) + length(POST_only_active);
                if denom > 0
                    temp = length(PRE_POST_active) / denom;
                else
                    temp = NaN;
                end
                V1_reactivation_coherence_modulation(track_id).jaccard_scores{hemi} = ...
                    [V1_reactivation_coherence_modulation(track_id).jaccard_scores{hemi}, temp];

                %%%% Persistence Score (proportion of neuron active PRE and POST were active during PRE)
                if length(PRE_active) > 0
                    temp = length(PRE_POST_active) / length(PRE_active);
                else
                    temp = NaN;
                end
                V1_reactivation_coherence_modulation(track_id).persistence_scores{hemi} = ...
                    [V1_reactivation_coherence_modulation(track_id).persistence_scores{hemi}, temp];

                %%%% Recruitment Score (proportion of neuron active POST were new recuited)
                if length(POST_active) > 0
                    temp = length(POST_only_active) / (length(active_cells)-length(PRE_only_active));
                    % temp = length(POST_only_active) / length(POST_active);
                else
                    temp = NaN;
                end
                V1_reactivation_coherence_modulation(track_id).recruitment_scores{hemi} = ...
                    [V1_reactivation_coherence_modulation(track_id).recruitment_scores{hemi}, temp];

                %%%% Correlation
                temp = corr(pre_PSTH(hemi_id==hemi), post_PSTH(hemi_id==hemi));
                V1_reactivation_coherence_modulation(track_id).PRE_POST_corr{hemi} = ...
                    [V1_reactivation_coherence_modulation(track_id).PRE_POST_corr{hemi}, temp];

                %%%% Cosine Similarity
                vec_pre = pre_PSTH(hemi_id==hemi);
                vec_post = post_PSTH(hemi_id==hemi);

                % Check for zero-vectors to avoid divide-by-zero errors
                if norm(vec_pre) > 0 && norm(vec_post) > 0
                    temp = dot(vec_pre', vec_post') / (norm(vec_pre') * norm(vec_post'));
                else
                    temp = NaN; % or 0, depending on preference
                end

                V1_reactivation_coherence_modulation(track_id).PRE_POST_cos{hemi} = ...
                    [V1_reactivation_coherence_modulation(track_id).PRE_POST_cos{hemi}, temp];
            end

            % figure;subplot(2,2,1);imagesc(PRE_corr);colormap(flipud(gray));
            % subplot(2,2,2);imagesc(POST_corr);colormap(flipud(gray));
            % subplot(2,2,3);imagesc(POST_corr-PRE_corr);colormap(flipud(gray));
            % subplot(2,2,4);imagesc(PRE_POST_corr);colormap(flipud(gray));
            %
            % figure;subplot(2,2,1);imagesc(tvec,1:size(PSTH,2),squeeze(PSTH(nevent,:,:)));hold on;
            % xline(-0.2,'r--');xline(0.2,'r--');xline(0,'r');colormap(flipud(gray));
            % yline(sum(hemi_id==1),'k','LineWidth',1)
            % subplot(2,2,3);plot(bin_centers,z_bias_V1_KDE(:,V1_reactivation_coherence_modulation(track_id).event_id(selected_event_id(nevent))));
            % hold on;plot(bin_centers,z_bias_KDE(:,V1_reactivation_coherence_modulation(track_id).event_id(selected_event_id(nevent))));
            % xline(-0.2,'r--');xline(0.2,'r--');xline(0,'r');xlim([-1 1])
            % 
            % subplot(2,2,2);imagesc(tvec,1:size(HPC_PSTH,2),squeeze(HPC_PSTH(nevent,:,:)));hold on;
            % xline(-0.2,'r--');xline(0.2,'r--');xline(0,'r');colormap(flipud(gray));
            % 
            % subplot(2,2,4);plot(bin_centers,z_bias_V1_KDE(:,V1_reactivation_coherence_modulation(track_id).event_id(selected_event_id(nevent))));
            % hold on;plot(bin_centers,z_bias_KDE(:,V1_reactivation_coherence_modulation(track_id).event_id(selected_event_id(nevent))));
            % xline(-0.2,'r--');xline(0.2,'r--');xline(0,'r');xlim([-1 1])

        end

    end
end


V1_reactivation_coherence_modulation_low = V1_reactivation_coherence_modulation;
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation_low.mat'),'V1_reactivation_coherence_modulation_low','-v7.3');


% V1_reactivation_coherence_modulation_low = V1_reactivation_coherence_modulation;
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'),'V1_reactivation_coherence_modulation','-v7.3');

%% Compare persistence, recuirtment score (low vs high bias ripple events)
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation_backup.mat'),'V1_reactivation_coherence_modulation');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'),'V1_reactivation_coherence_modulation');
V1_reactivation_coherence_modulation1 = V1_reactivation_coherence_modulation;

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation_low.mat'),'V1_reactivation_coherence_modulation_low');
V1_reactivation_coherence_modulation = V1_reactivation_coherence_modulation_low;
V1_reactivation_coherence_modulation = V1_reactivation_coherence_modulation1;

% V1_reactivation_coherence_modulation1(1).event_type(nbin,:)
%%%%%%%%%%%
nbin = 1;
T1_coherent = V1_reactivation_coherence_modulation(1).event_type(nbin,:)==1;
T1_coherent_PRE = V1_reactivation_coherence_modulation(1).event_type_PRE(nbin,:)==1;

T2_coherent = V1_reactivation_coherence_modulation(2).event_type(nbin,:)==1;
T2_coherent_PRE = V1_reactivation_coherence_modulation(2).event_type_PRE(nbin,:)==1;

%%%%%%%%%%%%
nbin = 1;
T1_coherent_low = V1_reactivation_coherence_modulation_low(1).event_type(nbin,:)==1;
T1_coherent_PRE_low = V1_reactivation_coherence_modulation_low(1).event_type_PRE(nbin,:)==1;

T2_coherent_low = V1_reactivation_coherence_modulation_low(2).event_type(nbin,:)==1;
T2_coherent_PRE_low = V1_reactivation_coherence_modulation_low(2).event_type_PRE(nbin,:)==1;

%%%%%%%%%%%%%
index1 = T1_coherent_PRE==1&T1_coherent==1;index2 = T2_coherent_PRE==1&T2_coherent==1;
index1_low = T1_coherent_PRE_low==1&T1_coherent_low==1;index2_low = T2_coherent_PRE_low==1&T2_coherent_low==1;

index1 = T1_coherent_PRE==1;index2 = T2_coherent_PRE==1;
index1_low = T1_coherent_PRE_low==1;index2_low = T2_coherent_PRE_low==1;

index1 = T1_coherent_PRE>=0;index2 = T2_coherent_PRE>=0;
index1_low = T1_coherent_PRE_low>=0;index2_low = T2_coherent_PRE_low>=0;

clear ripple_powers_all ripple_powers_all_low;
ripple_threshold = prctile(ripple_info.ripple_power,[25 75]);
ripple_powers_all{1} = ripple_info.ripple_power(V1_reactivation_coherence_modulation(1).event_id);
ripple_powers_all{2} = ripple_info.ripple_power(V1_reactivation_coherence_modulation(2).event_id);
ripple_powers_all_low{1} = ripple_info.ripple_power(V1_reactivation_coherence_modulation_low(1).event_id);
ripple_powers_all_low{2} = ripple_info.ripple_power(V1_reactivation_coherence_modulation_low(2).event_id);
index1 = ripple_powers_all{1}>ripple_threshold(end);index2 = ripple_powers_all{2}>ripple_threshold(end);
index1_low = ripple_powers_all_low{1} < ripple_threshold(1);index2_low = ripple_powers_all_low{2}<ripple_threshold(1);


p(1) =ranksum(V1_reactivation_coherence_modulation(1).persistence_scores{1}(index1),V1_reactivation_coherence_modulation_low(1).persistence_scores{1}(index1_low),'tail','left')
p(2) =ranksum(V1_reactivation_coherence_modulation(1).persistence_scores{2}(index1),V1_reactivation_coherence_modulation_low(1).persistence_scores{2}(index1_low),'tail','right')
% [~,p(2)] =kstest2(V1_reactivation_coherence_modulation(1).persistence_scores{2}(index1),V1_reactivation_coherence_modulation_low(1).persistence_scores{2}(index1_low),'tail','smaller')
p(3) =ranksum(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).persistence_scores{1}(index2_low),'tail','right')
% [~,p(3)] =kstest2(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).persistence_scores{1}(index2_low),'tail','smaller')
p(4) =ranksum(V1_reactivation_coherence_modulation(2).persistence_scores{2}(index2),V1_reactivation_coherence_modulation_low(2).persistence_scores{2}(index2_low),'tail','left')

p(5) =ranksum(V1_reactivation_coherence_modulation(1).recruitment_scores{1}(index1),V1_reactivation_coherence_modulation_low(1).recruitment_scores{1}(index1_low),'tail','left')
p(6) =ranksum(V1_reactivation_coherence_modulation(1).recruitment_scores{2}(index1),V1_reactivation_coherence_modulation_low(1).recruitment_scores{2}(index1_low),'tail','right')
p(7) =ranksum(V1_reactivation_coherence_modulation(2).recruitment_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).recruitment_scores{1}(index2_low),'tail','right')
p(8) =ranksum(V1_reactivation_coherence_modulation(2).recruitment_scores{2}(index2),V1_reactivation_coherence_modulation_low(2).recruitment_scores{2}(index2_low),'tail','left')


p(9) =ranksum(V1_reactivation_coherence_modulation(1).jaccard_scores{1}(index1),V1_reactivation_coherence_modulation_low(1).jaccard_scores{1}(index1_low),'tail','left')
p(10) =ranksum(V1_reactivation_coherence_modulation(1).jaccard_scores{2}(index1),V1_reactivation_coherence_modulation_low(1).jaccard_scores{2}(index1_low),'tail','right')
p(11) =ranksum(V1_reactivation_coherence_modulation(2).jaccard_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).jaccard_scores{1}(index2_low),'tail','right')
p(12) =ranksum(V1_reactivation_coherence_modulation(2).jaccard_scores{2}(index2),V1_reactivation_coherence_modulation_low(2).jaccard_scores{2}(index2_low),'tail','left')

p(13) =ranksum(V1_reactivation_coherence_modulation(1).PRE_POST_cos{1}(index1),V1_reactivation_coherence_modulation_low(1).PRE_POST_cos{1}(index1_low),'tail','left')
p(14) =ranksum(V1_reactivation_coherence_modulation(1).PRE_POST_cos{2}(index1),V1_reactivation_coherence_modulation_low(1).PRE_POST_cos{2}(index1_low),'tail','right')
p(15) =ranksum(V1_reactivation_coherence_modulation(2).PRE_POST_cos{1}(index2),V1_reactivation_coherence_modulation_low(2).PRE_POST_cos{1}(index2_low),'tail','right')
p(16) =ranksum(V1_reactivation_coherence_modulation(2).PRE_POST_cos{2}(index2),V1_reactivation_coherence_modulation_low(2).PRE_POST_cos{2}(index2_low),'tail','left')


% dist1 = V1_reactivation_coherence_modulation(1).V1_bias(index1)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(index1);
% dist2 = V1_reactivation_coherence_modulation_low(1).V1_bias(index1_low)-V1_reactivation_coherence_modulation_low(1).V1_bias_PRE(index1_low);
% p(17) =ranksum(dist1,dist2,'tail','right');
% 
% dist1 = V1_reactivation_coherence_modulation(2).V1_bias(index2)-V1_reactivation_coherence_modulation(2).V1_bias_PRE(index2);
% dist2 = V1_reactivation_coherence_modulation_low(2).V1_bias(index2_low)-V1_reactivation_coherence_modulation_low(2).V1_bias_PRE(index2_low);
% p(18) =ranksum(dist1,dist2,'tail','left');
dist1 = V1_reactivation_coherence_modulation(1).V1_bias(index1);
dist2 = V1_reactivation_coherence_modulation_low(1).V1_bias(index1_low);
p(17) =ranksum(dist1,dist2,'tail','right');

dist1 = V1_reactivation_coherence_modulation(2).V1_bias(index2);
dist2 = V1_reactivation_coherence_modulation_low(2).V1_bias(index2_low);
p(18) =ranksum(dist1,dist2,'tail','left');


% 
% subplot(2,2,1)
% histogram(V1_reactivation_coherence_modulation(1).persistence_scores{2}(index1)-V1_reactivation_coherence_modulation(1).persistence_scores{1}(index1));
% hold on;histogram(V1_reactivation_coherence_modulation_low(1).persistence_scores{2}(index1_low)-V1_reactivation_coherence_modulation_low(1).persistence_scores{1}(index1_low));
% 
% subplot(2,2,2)
% histogram(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2)-V1_reactivation_coherence_modulation(2).persistence_scores{2}(index2));
% hold on;histogram(V1_reactivation_coherence_modulation_low(2).persistence_scores{1}(index2_low)-V1_reactivation_coherence_modulation_low(2).persistence_scores{2}(index2_low));

% fig = figure('Name','V1 neuron persistence vs recruitment following ripple low vs high bias (all)');
% fig = figure('Name','V1 neuron persistence vs recruitment following ripple low vs high bias (PRE coherent)');
fig = figure('Name','V1 neuron persistence vs recruitment following ripple low vs high bias (PRE and POST coherent)');
fig.Position = [500 400 1130 510];
subplot(2,4,1)
histogram(V1_reactivation_coherence_modulation(1).persistence_scores{2}(index1),0:0.05:1,'Normalization','probability');hold on;histogram(V1_reactivation_coherence_modulation_low(1).persistence_scores{2}(index1_low),0:0.05:1,'Normalization','probability');
title('Track L ripple persistence score')
text(0.1,0.1,sprintf('p = %.3e',p(2)))
legend({'V1 R high','V1 R low'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,2)
histogram(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2),0:0.05:1,'Normalization','probability');hold on;histogram(V1_reactivation_coherence_modulation_low(2).persistence_scores{1}(index2_low),0:0.05:1,'Normalization','probability');
title('Track R ripple persistence score')
text(0.1,0.1,sprintf('p = %.3e',p(3)))
legend({'V1 L high','V1 L low'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(2,4,3)
histogram(V1_reactivation_coherence_modulation(1).recruitment_scores{2}(index1),0:0.05:1,'Normalization','probability');hold on;histogram(V1_reactivation_coherence_modulation_low(1).recruitment_scores{2}(index1_low),0:0.05:1,'Normalization','probability');
text(0.1,0.1,sprintf('p = %.3e',p(6)))
title('Track L ripple recruitment score')
legend({'V1 R high','V1 R low'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(2,4,4)
histogram(V1_reactivation_coherence_modulation(2).recruitment_scores{1}(index2),0:0.05:1,'Normalization','probability');hold on;histogram(V1_reactivation_coherence_modulation_low(2).recruitment_scores{1}(index2_low),0:0.05:1,'Normalization','probability');
text(0.1,0.1,sprintf('p = %.3e',p(7)))
title('Track R ripple recruitment score')
legend({'V1 L high','V1 L low'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,5)
histogram(V1_reactivation_coherence_modulation(1).PRE_POST_cos{2}(index1),0:0.05:1,'Normalization','probability');hold on;histogram(V1_reactivation_coherence_modulation_low(1).PRE_POST_cos{2}(index1_low),0:0.05:1,'Normalization','probability');
text(0.1,0.1,sprintf('p = %.3e',p(14)))
title('Track L ripple cos sim')
legend({'V1 R high','V1 R low'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(2,4,6)
histogram(V1_reactivation_coherence_modulation(2).PRE_POST_cos{1}(index2),0:0.05:1,'Normalization','probability');hold on;histogram(V1_reactivation_coherence_modulation_low(2).PRE_POST_cos{1}(index2_low),0:0.05:1,'Normalization','probability');
text(0.1,0.1,sprintf('p = %.3e',p(15)))
title('Track R ripple cos sim')
legend({'V1 L high','V1 L low'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



% 
% 
% %% V1 log odds PRE vs POST (POST coherent or incoherent)
% 
% T1_V1_L_biased = V1_reactivation_coherence_modulation(1).V1_bias>0;
% T1_V1_R_biased = V1_reactivation_coherence_modulation(1).V1_bias<0;
% T2_V1_L_biased = V1_reactivation_coherence_modulation(2).V1_bias>0;
% T2_V1_R_biased = V1_reactivation_coherence_modulation(2).V1_bias<0;
% 
% 
% T1_V1_L_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE>0;
% T1_V1_R_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE<0;
% T2_V1_L_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE>0;
% T2_V1_R_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE<0;
% 
% fig = figure('Name','V1 log odds PRE vs POST (POST coherent or incoherent)');
% fig.Position = [500 400 730 510];
% 
% subplot(2,2,1)
% scatter(V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased),V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased),...
%     'filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1)
% hold on
% scatter(V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased),V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased),...
%     'filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1)
% 
% % histogram()
% ylabel('PRE V1 log odds')
% xlabel('POST V1 log odds')
% xlim([-0.1 3])
% ylim([-3 3])
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('Track L ripple')
% 
% subplot(2,2,2)
% hold on
% scatter(V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased),V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased),...
%     'filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1)
% scatter(V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased),V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased),...
%     'filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1)
% % histogram()
% % ylabel('PRE V1 bias')
% % xlabel('POST V1 bias')
% ylabel('PRE V1 log odds')
% xlabel('POST V1 log odds')
% xlim([-3 0.1])
% ylim([-3 3])
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('Track R ripple')
% 
% 
% subplot(2,2,3)
% % dist1 = V1_reactivation_coherence_modulation(1).V1_bias(index1);
% % dist2 = V1_reactivation_coherence_modulation_low(1).V1_bias(index1_low);
% % dist1 = V1_reactivation_coherence_modulation(1).V1_bias(V1_reactivation_coherence_modulation(1).V1_bias_PRE>0)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(V1_reactivation_coherence_modulation(1).V1_bias_PRE>0);
% % dist2 = V1_reactivation_coherence_modulation(2).V1_bias(V1_reactivation_coherence_modulation(2).V1_bias_PRE>0)-V1_reactivation_coherence_modulation(2).V1_bias_PRE(V1_reactivation_coherence_modulation(2).V1_bias_PRE>0);
% dist1 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased);
% dist2 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased);
% 
% p(17) =ranksum(dist1,dist2,'tail','right');
% 
% histogram(dist1,-3:0.05:3,'Normalization','probability');hold on;histogram(dist2,-3:0.05:3,'Normalization','probability');
% text(1,0.05,sprintf('p = %.3e',p(17)))
% title('Track L ripple POST V1 bias')
% legend({'coherent','incoherent'},'box','off')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% xlabel('V1 Log odds (Track L -> +)')
% 
% 
% subplot(2,2,4)
% % dist1 = V1_reactivation_coherence_modulation(2).V1_bias(index2);
% % dist2 = V1_reactivation_coherence_modulation_low(2).V1_bias(index2_low);
% % dist2 = V1_reactivation_coherence_modulation(1).V1_bias(V1_reactivation_coherence_modulation(1).V1_bias_PRE<0)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(V1_reactivation_coherence_modulation(1).V1_bias_PRE<0);
% % dist1 = V1_reactivation_coherence_modulation(2).V1_bias(V1_reactivation_coherence_modulation(2).V1_bias_PRE<0)-V1_reactivation_coherence_modulation(2).V1_bias_PRE(V1_reactivation_coherence_modulation(2).V1_bias_PRE<0);
% 
% dist2 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased);
% dist1 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased);
% p(18) =ranksum(dist1,dist2,'tail','left');
% 
% histogram(dist1,-3:0.05:3,'Normalization','probability');hold on;histogram(dist2,-3:0.05:3,'Normalization','probability');
% text(-1,0.05,sprintf('p = %.3e',p(18)))
% title('Track R ripple POST V1 bias')
% legend({'coherent','incoherent'},'box','off')
% xlabel('V1 Log odds (- <-Track R)')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 


%% V1 log odds PRE vs POST (PRE coherent or incoherent)
nbin = 1;
T1_V1_L_biased = V1_reactivation_coherence_modulation(1).V1_bias>0;
T1_V1_R_biased = V1_reactivation_coherence_modulation(1).V1_bias<0;
T2_V1_L_biased = V1_reactivation_coherence_modulation(2).V1_bias>0;
T2_V1_R_biased = V1_reactivation_coherence_modulation(2).V1_bias<0;


T1_V1_L_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE>0;
T1_V1_R_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE<0;
T2_V1_L_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE>0;
T2_V1_R_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE<0;

% T1_V1_L_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE>0.5;
% T1_V1_R_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE<-0.5;
% T2_V1_L_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE>0.5;
% T2_V1_R_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE<-0.5;

mean_bias = mean(z_bias_KDE(bin_centers>0 & bin_centers<0.1,:),'omitnan');
T1_V1_L_biased_PRE_HPC_bias = mean_bias(V1_reactivation_coherence_modulation(1).event_id(V1_reactivation_coherence_modulation(1).V1_bias_PRE>0));
T1_V1_R_biased_PRE_HPC_bias =  mean_bias(V1_reactivation_coherence_modulation(1).event_id(V1_reactivation_coherence_modulation(1).V1_bias_PRE<0));
T2_V1_L_biased_PRE_HPC_bias =  mean_bias(V1_reactivation_coherence_modulation(2).event_id(V1_reactivation_coherence_modulation(2).V1_bias_PRE>0));
T2_V1_R_biased_PRE_HPC_bias =  mean_bias(V1_reactivation_coherence_modulation(2).event_id(V1_reactivation_coherence_modulation(2).V1_bias_PRE<0));


fig = figure('Name','V1 log odds PRE vs POST (PRE coherent or incoherent)');
fig.Position = [500 400 730 510];

subplot(2,2,1)
scatter(V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE),V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased_PRE),10,...
    'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on
scatter(V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE),V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE),10,...
    'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)

% histogram()
ylabel('PRE V1 log odds')
xlabel('POST V1 log odds')
ylim([-0.1 3])
xlim([-2.5 2.5])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Track L ripple')

subplot(2,2,2)
hold on
scatter(V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased_PRE),V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased_PRE),10,...
    'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
scatter(V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased_PRE),V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased_PRE),10,...
    'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
% histogram()
% ylabel('PRE V1 bias')
% xlabel('POST V1 bias')
ylabel('PRE V1 log odds')
xlabel('POST V1 log odds')
ylim([-3 0.1])
xlim([-2.5 2.5])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Track R ripple')


subplot(2,2,3)
% dist1 = V1_reactivation_coherence_modulation(1).V1_bias(index1);
% dist2 = V1_reactivation_coherence_modulation_low(1).V1_bias(index1_low);
% dist1 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased_PRE);
% dist2 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE)-V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE);
dist1 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE);
dist2 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE);

p(17) =ranksum(dist1,dist2,'tail','right');
histogram(dist1,-2.5:0.05:2.5,'Normalization','probability');hold on;histogram(dist2,-2:0.05:2,'Normalization','probability');
xlim([-2.5 2.5])
text(1,0.05,sprintf('p = %.3e',p(17)))
title('Track L ripple POST V1 bias')
legend({'coherent','incoherent'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('V1 Log odds (Track L -> +)')

subplot(2,2,4)
% dist1 = V1_reactivation_coherence_modulation(2).V1_bias(index2);
% dist2 = V1_reactivation_coherence_modulation_low(2).V1_bias(index2_low);
% dist2 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased_PRE)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased_PRE);
% dist1 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased_PRE)-V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased_PRE);

dist2 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased_PRE);
dist1 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased_PRE);
p(18) =ranksum(dist1,dist2,'tail','left');

histogram(dist1,-2.5:0.05:2.5,'Normalization','probability');hold on;histogram(dist2,-2:0.05:2,'Normalization','probability');
xlim([-2.5 2.5])
text(1,0.05,sprintf('p = %.3e',p(18)))
title('Track R ripple POST V1 bias')
legend({'coherent','incoherent'},'box','off')
xlabel('V1 Log odds (- <-Track R)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



%% HPC ripple modulation of V1 PRE vs POST log odds 

nbin = 1;
T1_V1_L_biased = V1_reactivation_coherence_modulation(1).V1_bias>0;
T1_V1_R_biased = V1_reactivation_coherence_modulation(1).V1_bias<0;
T2_V1_L_biased = V1_reactivation_coherence_modulation(2).V1_bias>0;
T2_V1_R_biased = V1_reactivation_coherence_modulation(2).V1_bias<0;


T1_V1_L_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE>0;
T1_V1_R_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE<0;
T2_V1_L_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE>0;
T2_V1_R_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE<0;

% T1_V1_L_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE>0.5;
% T1_V1_R_biased_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE<-0.5;
% T2_V1_L_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE>0.5;
% T2_V1_R_biased_PRE = V1_reactivation_coherence_modulation(2).V1_bias_PRE<-0.5;

mean_bias = mean(z_bias_KDE(bin_centers>0 & bin_centers<0.1,:),'omitnan');
T1_V1_L_biased_PRE_HPC_bias = mean_bias(V1_reactivation_coherence_modulation(1).event_id(V1_reactivation_coherence_modulation(1).V1_bias_PRE>0));
T1_V1_R_biased_PRE_HPC_bias =  mean_bias(V1_reactivation_coherence_modulation(1).event_id(V1_reactivation_coherence_modulation(1).V1_bias_PRE<0));
T2_V1_L_biased_PRE_HPC_bias =  mean_bias(V1_reactivation_coherence_modulation(2).event_id(V1_reactivation_coherence_modulation(2).V1_bias_PRE>0));
T2_V1_R_biased_PRE_HPC_bias =  mean_bias(V1_reactivation_coherence_modulation(2).event_id(V1_reactivation_coherence_modulation(2).V1_bias_PRE<0));


x1 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE);% Coherent
x2 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE);% Incoherent
y1 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased_PRE);% Coherent
y2 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE);% Incoherent

log_odds_thresholds = prctile([y1 y2],0:10:100);
log_odds_thresholds1 = log_odds_thresholds;

x1 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased_PRE);% Coherent
x2 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased_PRE);% Incoherent
y1 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased_PRE);% Coherent
y2 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased_PRE);% Incoherent

log_odds_thresholds = prctile([y1 y2],0:10:100);
log_odds_thresholds2 = log_odds_thresholds;


x1 = [V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased_PRE) V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased_PRE)];% Coherent
x2 =  [V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased_PRE) V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE)];% incoherent
y1 = [V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE) V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased_PRE)];% Coherent
y2 =  [V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased_PRE) V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE)];% incoherent

fig = figure('Name','V1 log odds PRE vs POST');
fig.Position = [500 400 730 510];
subplot(2,2,1)
p = polyfit([x1 x2],[y1 y2], 1);
% 2. Calculate the fitted values (the line)
b_fit = polyval(p, [x1 x2]);

hold on
scatter([x1 x2],[y1 y2],10,'k',...
    'filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)
plot([x1 x2], b_fit, 'r-', 'LineWidth', 2);            % Regression line

ylim([-2 2])
xlim([-2 2])
xlabel('PRE V1 log odds')
ylabel('POST V1 log odds')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,2)
hold on
scatter(x1,y1,10,'m',...
    'filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)
scatter(x2,y2,10,'k',...
    'filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)

% scatter(V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE),V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE),10,...
%     'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)

ylim([-2 2])
xlim([-2 2])
xlabel('PRE V1 log odds')
ylabel('POST V1 log odds')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


fig = figure('Name','V1 log odds PRE distribution');
fig.Position = [500 400 730/2 510/3];

dist1= [V1_reactivation_coherence_modulation(1).V1_bias_PRE V1_reactivation_coherence_modulation(2).V1_bias_PRE];

histogram(dist1,-2.5:0.05:2.5,'Normalization','probability','EdgeAlpha',0)
hold on
xline(log_odds_thresholds2,'r')
hold on
xline(log_odds_thresholds1,'b')
xlim([-2 2])
xlabel('PRE V1 Log odds')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])


%%%%%%%%%%%%%% Mixed effect model
%%%%%%%%%%%%%%
ripple_track = [2*ones(1,length(V1_reactivation_coherence_modulation(1).V1_bias)) ones(1,length(V1_reactivation_coherence_modulation(2).V1_bias))]';
v1_pre = [V1_reactivation_coherence_modulation(1).V1_bias V1_reactivation_coherence_modulation(2).V1_bias]';
v1_post = [V1_reactivation_coherence_modulation(1).V1_bias_PRE V1_reactivation_coherence_modulation(2).V1_bias_PRE]';
session_count = [V1_reactivation_coherence_modulation(1).session_id V1_reactivation_coherence_modulation(2).session_id]';
% subject_id = session_count;
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);



% Mixed effect model with shuffling
beta_shuffled = nan(1,1000);
t_shuffled = nan(1,1000);
beta_boot = nan(1,1000);
t_boot = nan(1,1000);

parfor iBoot = 1:1000
    tic
    s = RandStream('philox4x32_10', 'Seed', iBoot);
    % idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);
    idx = datasample(1:length(subject_id), length(subject_id), 'Replace', false);

    tbl = table(subject_id, zscore(v1_pre), ripple_track(idx), zscore(v1_post),'VariableNames',{'subject_id','v1_pre','ripple_track','v1_post'});
    tbl.subject_id = categorical(tbl.subject_id);
    tbl.ripple_track = categorical(tbl.ripple_track); % Important for the model

    %Run Linear Mixed-Effects Model
    % Formula: Post ~ Pre + Ripple + Interaction + (Random Intercept per Animal)
    lme = fitlme(tbl, 'v1_post ~ v1_pre + ripple_track + (1|subject_id)');
    beta_shuffled(iBoot) = lme.Coefficients.Estimate(3);
    t_shuffled(iBoot) = lme.Coefficients.tStat(3);

    idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);

    tbl = table(subject_id(idx), zscore(v1_pre(idx)), ripple_track(idx), zscore(v1_post(idx)),'VariableNames',{'subject_id','v1_pre','ripple_track','v1_post'});
    tbl.subject_id = categorical(tbl.subject_id);
    tbl.ripple_track = categorical(tbl.ripple_track); % Important for the model

    %Run Linear Mixed-Effects Model
    % Formula: Post ~ Pre + Ripple + Interaction + (Random Intercept per Animal)
    lme = fitlme(tbl, 'v1_post ~ v1_pre + ripple_track + (1|subject_id)');
    beta_boot(iBoot) = lme.Coefficients.Estimate(3);
    t_boot(iBoot) = lme.Coefficients.tStat(3);
    toc
end


% Mixed effect model 
tbl = table(subject_id(idx), v1_pre((idx)), ripple_track((idx)), v1_post((idx)),'VariableNames',{'subject_id','v1_pre','ripple_track','v1_post'});
tbl.subject_id = categorical(tbl.subject_id);
tbl.ripple_track = categorical(tbl.ripple_track); % Important for the model
lme = fitlme(tbl, 'v1_post ~ v1_pre * ripple_track + (1|subject_id)');
p = lme.Coefficients.pValue(3);


fig = figure('Name','Ripple modulation of POST V1 log odds Beta coefficient and t-stat');
fig.Position = [500 400 330 600];
subplot(2,2,1)
% 1. Calculate Means
mean_boot = mean(beta_boot, 'omitnan');
mean_shuf = mean(beta_shuffled, 'omitnan');

% 2. Calculate 95% Confidence Intervals (Percentile Method)
% This calculates the distance from the mean to the 2.5th and 97.5th percentiles
ci_boot = [mean_boot - prctile(beta_boot, 2.5), prctile(beta_boot, 97.5) - mean_boot];
ci_shuf = [mean_shuf - prctile(beta_shuffled, 2.5), prctile(beta_shuffled, 97.5) - mean_shuf];

% Prepare data for plotting
means = [mean_boot, mean_shuf];
lower_err = [ci_boot(1), ci_shuf(1)];
upper_err = [ci_boot(2), ci_shuf(2)];

% 3. Plotting
hold on;
x = [1, 2];
colors = [231, 41, 138; 256/2 256/2 256/2]/256; % Blue for Boot, Gray for Shuffled

for i = 1:2
    % Plot Bar
    bar(x(i), means(i), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    % Plot 95% CI Error Bar
    errorbar(x(i), means(i), lower_err(i), upper_err(i), 'k', 'LineWidth', 1.5, 'CapSize', 10);
end

% Formatting
xticks([1 2]);
xticklabels({'Beta Boot', 'Beta Shuffled'});
ylabel('Beta Value');
title('Mean Beta with 95% CI');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
text(1,0.07,sprintf('p = %.3e',p))


subplot(2,2,2)
% 1. Calculate Means
mean_boot = mean(t_boot, 'omitnan');
mean_shuf = mean(t_shuffled, 'omitnan');

% 2. Calculate 95% Confidence Intervals (Percentile Method)
% This calculates the distance from the mean to the 2.5th and 97.5th percentiles
ci_boot = [mean_boot - prctile(t_boot, 2.5), prctile(t_boot, 97.5) - mean_boot];
ci_shuf = [mean_shuf - prctile(t_shuffled, 2.5), prctile(t_shuffled, 97.5) - mean_shuf];

% Prepare data for plotting
means = [mean_boot, mean_shuf];
lower_err = [ci_boot(1), ci_shuf(1)];
upper_err = [ci_boot(2), ci_shuf(2)];

% 3. Plotting
hold on;
x = [1, 2];
% colors = [0.2 0.4 0.8; 0.5 0.5 0.5]; % Blue for Boot, Gray for Shuffled

for i = 1:2
    % Plot Bar
    bar(x(i), means(i), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    % Plot 95% CI Error Bar
    errorbar(x(i), means(i), lower_err(i), upper_err(i), 'k', 'LineWidth', 1.5, 'CapSize', 10);
end

% Formatting
xticks([1 2]);
xticklabels({'t stat Boot', 't stat Shuffled'});
ylabel('t-stat Value');
title('Mean t-stat with 95% CI');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
text(1,2,sprintf('p = %.3e',p))

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])


%%%%%%% ranksum POST ripple V1 log odds difference modulated by ripple HPC log
%%%%%%% odds 
nBins = 10;
clear PRE_POST_hist PRE_POST_p PRE_POST_mean

bin_edges = -2:0.1:2;
bin_centers = bin_edges(1)+mean(diff(bin_edges))/2:mean(diff(bin_edges)):bin_edges(end)-mean(diff(bin_edges))/2;

x1 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE);% Coherent
x2 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE);% Incoherent
y1 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased_PRE);% Coherent
y2 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE);% Incoherent

log_odds_thresholds = prctile([y1 y2],0:100/nBins:100);
log_odds_thresholds1 = log_odds_thresholds;
for nbin = 1:length(log_odds_thresholds)-1
        dist1 = x1(y1 >= log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1));
        dist2 = x2(y2 >= log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1));

        PRE_POST_p(1,nbin) = ranksum(dist1,dist2,'tail','right');
        % [~,PRE_POST_p(1,nbin)] = kstest2(dist1,dist2,'Tail','smaller');
        % [~,PRE_POST_p(1,nbin)] = kstest([dist1-dist2],'Tail','smaller');
        % PRE_POST_mean(1,2,nbin) = x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1);

        PRE_POST_mean(1,1,nbin) = mean(x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1)),'omitnan');
        PRE_POST_mean(1,2,nbin) = mean(x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1)),'omitnan');

        PRE_POST_hist{1}(1,nbin,:) = movmean(histcounts(x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1)),-2:0.1:2,'Normalization', 'probability'),4);
        PRE_POST_hist{1}(2,nbin,:) = movmean(histcounts(x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1)),-2:0.1:2,'Normalization', 'probability'),4);
        % 
        % nexttile
        % histogram(dist1,-2:0.1:2,'Normalization', 'probability');hold on;histogram(dist2,-2:0.1:2,'Normalization', 'probability')
end

x1 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_R_biased_PRE);% Coherent
x2 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_R_biased_PRE);% Incoherent
y1 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_R_biased_PRE);% Coherent
y2 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_R_biased_PRE);% Incoherent

log_odds_thresholds = prctile([y1 y2],0:100/nBins:100);
log_odds_thresholds2 = log_odds_thresholds;
for nbin = 1:length(log_odds_thresholds)-1
        dist1 = x1(y1 >= log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1));
        dist2 = x2(y2 >= log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1));

        PRE_POST_p(2,nbin) = ranksum(dist1,dist2,'tail','left');
        % [~,PRE_POST_p(2,nbin)] = kstest2(dist1,dist2,'Tail','larger');
        % PRE_POST_mean(1,2,nbin) = x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1);

        PRE_POST_mean(2,1,nbin) = mean(x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1)),'omitnan');
        PRE_POST_mean(2,2,nbin) = mean(x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1)),'omitnan');

        PRE_POST_hist{2}(1,nbin,:) = movmean(histcounts(x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1)),-2:0.1:2,'Normalization', 'probability'),4);
        PRE_POST_hist{2}(2,nbin,:) = movmean(histcounts(x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1)),-2:0.1:2,'Normalization', 'probability'),4);
end

% bin_edges = -3+0.05/2:0.05:3-0.05/2;
bin_edges = -2:0.1:2;
bin_centers = bin_edges(1)+mean(diff(bin_edges))/2:mean(diff(bin_edges)):bin_edges(end)-mean(diff(bin_edges))/2;



% fig = figure('Name','V1 log odds POST ripple modulation (PRE coherent or incoherent)');
fig = figure('Name','V1 log odds POST ripple modulation (percentile) (PRE coherent or incoherent)');
fig.Position = [500 400 1800 510];
subplot(2,4,1)
imagesc(1:nBins,bin_centers,squeeze(PRE_POST_hist{1}(1,:,:))')
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds1([1 5 10]))
colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([0 0.1])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
title('L Ripple with pre-ripple L bias in V1')

subplot(2,4,2)
imagesc(1:nBins,bin_centers,squeeze(PRE_POST_hist{1}(2,:,:))')
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds1([1 5 10]))
colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([0 0.1])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
title('R Ripple with pre-ripple L bias in V1')

subplot(2,4,3)
imagesc(1:nBins,bin_centers,squeeze(PRE_POST_hist{2}(1,:,:))')
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds2([2 6 11]))

colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([0 0.1])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
title('R Ripple with pre-ripple R bias in V1')

subplot(2,4,4)
imagesc(1:nBins,bin_centers,squeeze(PRE_POST_hist{2}(2,:,:))')
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds2([2 6 11]))
colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([0 0.1])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
title('L Ripple with pre-ripple R bias in V1')


subplot(2,4,5)
imagesc(1:nBins,bin_centers,squeeze(PRE_POST_hist{1}(1,:,:))'-squeeze(PRE_POST_hist{1}(2,:,:))')
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% imagesc(bin_centers,1:10,squeeze(PRE_POST_hist{1}(1,:,:))-squeeze(PRE_POST_hist{1}(2,:,:)))
% yticks([1 5 10])
% yticklabels(log_odds_thresholds1([1 5 10]))
colormap(flipud(gray));yline(0,'r--');colorbar;clim([-0.02 0.02])
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ripple with pre-ripple L bias in V1 (L-R)')

subplot(2,4,6)
imagesc(1:nBins,bin_centers,squeeze(PRE_POST_hist{2}(1,:,:))'-squeeze(PRE_POST_hist{2}(2,:,:))')
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds2([2 6 11]))
colormap(flipud(gray));yline(0,'r--');yline(0,'r--');colorbar;clim([-0.02 0.02])
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ripple with pre-ripple R bias in V1 (R-L)')

subplot(2,4,7)
% bar(log_odds_thresholds1(1:end-1),-log10((PRE_POST_ranksum(1,:))));hold on;
% bar(log_odds_thresholds2(2:end),-log10(PRE_POST_ranksum(2,:)))
bar(nBins+1:2*nBins,-log10((PRE_POST_p(1,:))));hold on;
bar(1:nBins,-log10(PRE_POST_p(2,:)))
yline(-log10(0.05),'--r','p < 0.05')
xlabel('PRE ripple V1 log odds')
ylabel('-log10(p)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ripple with different pre-ripple in V1')

subplot(2,4,8)
bar(nBins+1:2*nBins,squeeze(squeeze(PRE_POST_mean(1,1,:)))- squeeze(squeeze(PRE_POST_mean(1,2,:))));hold on;
% bar(log_odds_thresholds1(1:end-1),squeeze(PRE_POST_mean(1,2,:)))
bar(1:nBins,squeeze(squeeze(PRE_POST_mean(2,1,:)))- squeeze(squeeze(PRE_POST_mean(2,2,:))));hold on;
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds diff')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Mean POST ripple V1 log odds diff')

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])


% fig = figure('Name','V1 log odds POST ripple modulation (PRE coherent or incoherent)');
fig = figure('Name','V1 log odds POST ripple modulation (percentile plot merged) (PRE coherent or incoherent)');
fig.Position = [500 400 1800*2/3 510];
subplot(2,4,1)
imagesc([-nBins:-1 1:nBins],bin_centers,[squeeze(PRE_POST_hist{2}(1,:,:))' squeeze(PRE_POST_hist{1}(1,:,:))']);
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds1([1 5 10]))
% ylim([-1.5 1.5])
colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([0 0.1])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
title('Ripple matching PRE V1 bias')

subplot(2,4,2)
imagesc([-nBins:-1 1:nBins],bin_centers,[squeeze(PRE_POST_hist{2}(2,:,:))' squeeze(PRE_POST_hist{1}(2,:,:))'])
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
% yticks([1 5 10])
% yticklabels(log_odds_thresholds1([1 5 10]))
% ylim([-1.5 1.5])
colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([0 0.1])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
title('Ripple not matching PRE V1 bias')

subplot(2,4,5)
dist1 = [squeeze(PRE_POST_hist{2}(1,:,:))' squeeze(PRE_POST_hist{1}(1,:,:))'];
dist2 = [squeeze(PRE_POST_hist{2}(2,:,:))' squeeze(PRE_POST_hist{1}(2,:,:))'];
imagesc([-nBins:-1 1:nBins],bin_centers,[dist1-dist2])
% imagesc(1:10,bin_centers,squeeze(PRE_POST_hist{1}(1,:,:))'-squeeze(PRE_POST_hist{1}(2,:,:))')
% imagesc(bin_centers,1:10,squeeze(PRE_POST_hist{1}(1,:,:))-squeeze(PRE_POST_hist{1}(2,:,:)))
% yticks([1 5 10])
% yticklabels(log_odds_thresholds1([1 5 10]))
set(gca, 'YDir', 'normal'); % Forces Y-axis to increase upward
colormap(flipud(gray));yline(0,'r--');xline(0,'r--');colorbar;clim([-0.015 0.015])
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ripple with pre-ripple L bias in V1 (L-R)')

% subplot(2,4,6)

subplot(2,4,7)
% bar(log_odds_thresholds1(1:end-1),-log10((PRE_POST_ranksum(1,:))));hold on;
% bar(log_odds_thresholds2(2:end),-log10(PRE_POST_ranksum(2,:)))
bar(nBins+1:2*nBins,-log10((PRE_POST_p(1,:))));hold on;
bar(1:nBins,-log10(PRE_POST_p(2,:)))
yline(-log10(0.05),'--r','p < 0.05')
xlabel('PRE ripple V1 log odds')
ylabel('-log10(p)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ripple with different pre-ripple in V1')

subplot(2,4,8)
bar(nBins+1:2*nBins,squeeze(squeeze(PRE_POST_mean(1,1,:)))- squeeze(squeeze(PRE_POST_mean(1,2,:))));hold on;
% bar(log_odds_thresholds1(1:end-1),squeeze(PRE_POST_mean(1,2,:)))
bar(1:nBins,squeeze(squeeze(PRE_POST_mean(2,1,:)))- squeeze(squeeze(PRE_POST_mean(2,2,:))));hold on;
ylim([-0.18 0.18])
xlabel('PRE ripple V1 log odds')
ylabel('POST ripple V1 log odds diff')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Mean POST ripple V1 log odds diff')



save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% combined

bin_edges = -2:0.1:2;
bin_centers = bin_edges(1)+mean(diff(bin_edges))/2:mean(diff(bin_edges)):bin_edges(end)-mean(diff(bin_edges))/2;

x1 = V1_reactivation_coherence_modulation(1).V1_bias(T1_V1_L_biased_PRE);% Coherent
x2 = V1_reactivation_coherence_modulation(2).V1_bias(T2_V1_L_biased_PRE);% Incoherent
y1 = V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_V1_L_biased_PRE);% Coherent
y2 = V1_reactivation_coherence_modulation(2).V1_bias_PRE(T2_V1_L_biased_PRE);% Incoherent

log_odds_thresholds = prctile([y1 y2],0:10:100);
log_odds_thresholds1 = log_odds_thresholds;
for nbin = 1:length(log_odds_thresholds)-1
        dist1 = x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1));
        dist2 = x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1));

        % PRE_POST_ranksum(1,nbin) = ranksum(dist1,dist2,'tail','right');
        [~,PRE_POST_p(1,nbin)] = kstest2(dist1,dist2,'Tail','smaller');
        % PRE_POST_mean(1,2,nbin) = x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1);

        PRE_POST_mean(1,1,nbin) = mean(x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1)),'omitnan');
        PRE_POST_mean(1,2,nbin) = mean(x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1)),'omitnan');

        PRE_POST_hist{1}(1,nbin,:) = movmean(histcounts(x1(y1 > log_odds_thresholds(nbin) & y1 < log_odds_thresholds(nbin+1)),-2:0.1:2,'Normalization', 'probability'),4);
        PRE_POST_hist{1}(2,nbin,:) = movmean(histcounts(x2(y2 > log_odds_thresholds(nbin) & y2 < log_odds_thresholds(nbin+1)),-2:0.1:2,'Normalization', 'probability'),4);

        % nexttile
        % histogram(dist1,-2:0.1:2,'Normalization', 'probability');hold on;histogram(dist2,-2:0.1:2,'Normalization', 'probability')
end

%%
dist1 = V1_reactivation_coherence_modulation(1).PRE_POST_cos{2}(V1_reactivation_coherence_modulation(1).V1_bias_PRE>0); % coherent
dist2 = V1_reactivation_coherence_modulation(2).PRE_POST_cos{2}(V1_reactivation_coherence_modulation(2).V1_bias_PRE>0); % non-coherent
% dist1 = V1_reactivation_coherence_modulation_low(1).jaccard_scores{2}(V1_reactivation_coherence_modulation_low(1).V1_bias>0);
% dist2 = V1_reactivation_coherence_modulation_low(2).jaccard_scores{2}(V1_reactivation_coherence_modulation_low(2).V1_bias>0);
subplot(2,4,1)
p(19) =ranksum(dist1,dist2,'tail','right');
histogram(dist1,0:0.05:1,'Normalization','probability');hold on;histogram(dist2,0:0.05:1,'Normalization','probability');
text(0.1,0.1,sprintf('p = %.3e',p(19)))
title('V1 R jaccard scores')
legend({'Track L ripple','Track R ripple'},'box','off')
xlabel('Jaccard score')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



dist1 = V1_reactivation_coherence_modulation(1).PRE_POST_cos{1}(V1_reactivation_coherence_modulation(1).V1_bias_PRE<0);% non coherent
dist2 = V1_reactivation_coherence_modulation(2).PRE_POST_cos{1}(V1_reactivation_coherence_modulation(2).V1_bias_PRE<0);% coherent
% dist1 = V1_reactivation_coherence_modulation_low(1).jaccard_scores{1}(V1_reactivation_coherence_modulation_low(1).V1_bias<0);% non coherent
% dist2 = V1_reactivation_coherence_modulation_low(2).jaccard_scores{1}(V1_reactivation_coherence_modulation_low(2).V1_bias<0);% coherent
subplot(2,4,2)
p(20) =ranksum(dist1,dist2,'tail','left');
histogram(dist1,0:0.05:1,'Normalization','probability');hold on;histogram(dist2,0:0.05:1,'Normalization','probability');
text(0.1,0.1,sprintf('p = %.3e',p(20)))
title('Track R ripple POST V1 bias')
legend({'Track L ripple','Track R ripple'},'box','off')
title('V1 L jaccard scores')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])

%%

%% Compare persistence, recuirtment score
p =ranksum(V1_reactivation_coherence_modulation(1).PRE_POST_corr{1},V1_reactivation_coherence_modulation(2).PRE_POST_corr{1},'tail','left')
p =ranksum(V1_reactivation_coherence_modulation(1).PRE_POST_corr{2},V1_reactivation_coherence_modulation(2).PRE_POST_corr{2},'tail','right')


p =ranksum(V1_reactivation_coherence_modulation(1).PRE_POST_cos{1},V1_reactivation_coherence_modulation(2).PRE_POST_cos{1},'tail','left')
p =ranksum(V1_reactivation_coherence_modulation(1).PRE_POST_cos{2},V1_reactivation_coherence_modulation(2).PRE_POST_cos{2},'tail','right')


p =ranksum(V1_reactivation_coherence_modulation(1).persistence_scores{1},V1_reactivation_coherence_modulation(2).persistence_scores{1},'tail','left')
p =ranksum(V1_reactivation_coherence_modulation(1).persistence_scores{2},V1_reactivation_coherence_modulation(2).persistence_scores{2},'tail','right')


p =ranksum(V1_reactivation_coherence_modulation(1).recruitment_scores{1},V1_reactivation_coherence_modulation(2).recruitment_scores{1},'tail','left')
p =ranksum(V1_reactivation_coherence_modulation(1).recruitment_scores{2},V1_reactivation_coherence_modulation(2).recruitment_scores{2},'tail','right')


p =ranksum(V1_reactivation_coherence_modulation(1).jaccard_scores{1},V1_reactivation_coherence_modulation(2).jaccard_scores{1},'tail','left')
p =ranksum(V1_reactivation_coherence_modulation(1).jaccard_scores{2},V1_reactivation_coherence_modulation(2).jaccard_scores{2},'tail','right')


T1_index = ismember(ripple_info.event_id,V1_reactivation_coherence_modulation(1).event_id);
T2_index = ismember(ripple_info.event_id,V1_reactivation_coherence_modulation(2).event_id);

clear ripple_powers_all;
ripple_threshold = prctile(ripple_info.ripple_power,75);
ripple_powers_all{1} = ripple_info.ripple_power(V1_reactivation_coherence_modulation(1).event_id);
ripple_powers_all{2} = ripple_info.ripple_power(V1_reactivation_coherence_modulation(2).event_id);


clear SO_phase_all;
% ripple_threshold = prctile(ripple_info.ripple_power,75);
SO_phase_all{1} = ripple_info.SO_phase(V1_reactivation_coherence_modulation(1).event_id,2);
SO_phase_all{2} = ripple_info.SO_phase(V1_reactivation_coherence_modulation(2).event_id,1);


% fig = figure('Name','V1 neuron persistence vs recruitment following ripple (POST incoherent)');
% fig = figure('Name','V1 neuron persistence vs recruitment following ripple (PRE and POST incoherent)');
% fig = figure('Name','V1 neuron persistence vs recruitment following ripple (PRE incoherent)');

% fig = figure('Name','V1 neuron persistence vs recruitment following ripple (POST coherent)');
% fig = figure('Name','V1 neuron persistence vs recruitment following ripple (PRE and POST coherent)');
fig = figure('Name','V1 neuron persistence vs recruitment following ripple (PRE coherent)');

% fig = figure('Name','V1 neuron persistence vs recruitment following low bias ripple (PRE coherent)');
fig.Position = [500 400 1130 510];

index1 = T1_coherent_PRE==1;index2 = T2_coherent_PRE==1;
% index1 = T1_coherent==1;index2 = T2_coherent==1;
% index1 = T1_coherent_PRE==1&T1_coherent==1;index2 = T2_coherent_PRE==1&T2_coherent==1;
% index1 = V1_reactivation_coherence_modulation(1).neurons_active(1,:)>0&V1_reactivation_coherence_modulation(1).neurons_active(2,:)>0&T1_coherent_PRE==1&T1_coherent==1;
% index2 = V1_reactivation_coherence_modulation(2).neurons_active(1,:)>0&V1_reactivation_coherence_modulation(2).neurons_active(2,:)>0&T2_coherent_PRE==1&T2_coherent==1;
% index1 = V1_reactivation_coherence_modulation(1).neurons_active(1,:)>0&V1_reactivation_coherence_modulation(1).neurons_active(2,:)>0;
% index2 = V1_reactivation_coherence_modulation(2).neurons_active(1,:)>0&V1_reactivation_coherence_modulation(2).neurons_active(2,:)>0;
% index1 = T1_coherent_PRE==1&T1_coherent==1;index2 = T2_coherent_PRE==1&T2_coherent==1;

% index1 = T1_coherent_PRE==0;index2 = T2_coherent_PRE==0;
% index1 = T1_coherent==0;index2 = T2_coherent==0;
% index1 = ripple_powers_all{1}>ripple_threshold(end);index2 =ripple_powers_all{2}>ripple_threshold(end);
% 
% % Trough
% index1 = SO_phase_all{1}>pi/2 | SO_phase_all{1}<-pi/2; index2 =SO_phase_all{2}>pi/2 | SO_phase_all{2}<-pi/2;
% 
% % Peak
% index1 = SO_phase_all{1}>-pi/2 & SO_phase_all{1}<pi/2; index2 =SO_phase_all{2}>-pi/2 & SO_phase_all{2}<pi/2;




% index1 = T1_coherent>=0;index2 = T2_coherent>=0;
p=[];
p(1) =signrank(V1_reactivation_coherence_modulation(1).persistence_scores{1}(index1),V1_reactivation_coherence_modulation(1).persistence_scores{2}(index1),'tail','left')
p(2) =signrank(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2),V1_reactivation_coherence_modulation(2).persistence_scores{2}(index2),'tail','right')
p(3) =signrank(V1_reactivation_coherence_modulation(1).recruitment_scores{1}(index1),V1_reactivation_coherence_modulation(1).recruitment_scores{2}(index1),'tail','left')
p(4) =signrank(V1_reactivation_coherence_modulation(2).recruitment_scores{1}(index2),V1_reactivation_coherence_modulation(2).recruitment_scores{2}(index2),'tail','right')
p(5) =signrank(V1_reactivation_coherence_modulation(1).PRE_POST_cos{1}(index1),V1_reactivation_coherence_modulation(1).PRE_POST_cos{2}(index1),'tail','left')
p(6) =signrank(V1_reactivation_coherence_modulation(2).PRE_POST_cos{1}(index2),V1_reactivation_coherence_modulation(2).PRE_POST_cos{2}(index2),'tail','right')
p(7) =signrank(V1_reactivation_coherence_modulation(1).jaccard_scores{1}(index1),V1_reactivation_coherence_modulation(1).jaccard_scores{2}(index1),'tail','left')
p(8) =signrank(V1_reactivation_coherence_modulation(2).jaccard_scores{1}(index2),V1_reactivation_coherence_modulation(2).jaccard_scores{2}(index2),'tail','right')



% 
% p(2) =ranksum(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).persistence_scores{2}(index2_low),'tail','right')
% p(3) =signrank(V1_reactivation_coherence_modulation(1).recruitment_scores{1}(index1),V1_reactivation_coherence_modulation_low(1).recruitment_scores{2}(index1),'tail','left')
% p(4) =signrank(V1_reactivation_coherence_modulation(2).recruitment_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).recruitment_scores{2}(index2),'tail','right')
% p(5) =signrank(V1_reactivation_coherence_modulation(1).PRE_POST_cos{1}(index1),V1_reactivation_coherence_modulation_low(1).PRE_POST_cos{2}(index1),'tail','left')
% p(6) =signrank(V1_reactivation_coherence_modulation(2).PRE_POST_cos{1}(index2),V1_reactivation_coherence_modulation_low(2).PRE_POST_cos{2}(index2),'tail','right')
% p(7) =signrank(V1_reactivation_coherence_modulation(1).jaccard_scores{1}(index1),V1_reactivation_coherence_modulation_low(1).jaccard_scores{2}(index1),'tail','left')
% p(8) =signrank(V1_reactivation_coherence_modulation(2).jaccard_scores{1}(index2),V1_reactivation_coherence_modulation_low(2).jaccard_scores{2}(index2),'tail','right')


subplot(2,4,1)
histogram(V1_reactivation_coherence_modulation(1).persistence_scores{1}(index1),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(1).persistence_scores{2}(index1),0:0.05:1,'Normalization','probability')
title('Track L ripple persistence score')
text(0.1,0.1,sprintf('p = %.3e',p(1)))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,2)
hold on;
histogram(V1_reactivation_coherence_modulation(2).persistence_scores{1}(index2),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(2).persistence_scores{2}(index2),0:0.05:1,'Normalization','probability')
title('Track R ripple persistence score')
text(0.1,0.1,sprintf('p = %.3e',p(2)))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,3)
histogram(V1_reactivation_coherence_modulation(1).recruitment_scores{1}(index1),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(1).recruitment_scores{2}(index1),0:0.05:1,'Normalization','probability')
title('Track L ripple recuirtment score')
text(0.1,0.1,sprintf('p = %.3e',p(3)))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,4)
hold on;
histogram(V1_reactivation_coherence_modulation(2).recruitment_scores{1}(index2),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(2).recruitment_scores{2}(index2),0:0.05:1,'Normalization','probability')
title('Track R ripple recuirtment score')
text(0.1,0.1,sprintf('p = %.3e',p(4)))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,5)
histogram(V1_reactivation_coherence_modulation(1).PRE_POST_cos{1}(index1),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(1).PRE_POST_cos{2}(index1),0:0.05:1,'Normalization','probability')
title('Track L ripple cos sim')
text(0.1,0.1,sprintf('p = %.3e',p(5)))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,6)
hold on;
histogram(V1_reactivation_coherence_modulation(2).PRE_POST_cos{1}(index2),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(2).PRE_POST_cos{2}(index2),0:0.05:1,'Normalization','probability')
title('Track R ripple cos sim')
text(0.1,0.1,sprintf('p = %.3e',p(6)))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,7)
histogram(V1_reactivation_coherence_modulation(1).jaccard_scores{1}(index1),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(1).jaccard_scores{2}(index1),0:0.05:1,'Normalization','probability')
text(0.1,0.1,sprintf('p = %.3e',p(7)))
title('Track R ripple jaccard scores')
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,4,8)
histogram(V1_reactivation_coherence_modulation(2).jaccard_scores{1}(index2),0:0.05:1,'Normalization','probability')
hold on;
histogram(V1_reactivation_coherence_modulation(2).jaccard_scores{2}(index2),0:0.05:1,'Normalization','probability')
text(0.1,0.1,sprintf('p = %.3e',p(8)))
title('Track R ripple jaccard scores')
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])



%% Distribution of neurons active before vs after ripple

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),'V1_reactivation_coherence_modulation.mat');

% T1_coherent_PRE = V1_reactivation_coherence_modulation(1).V1_bias_PRE>0;
nbin = 1;
T1_coherent = V1_reactivation_coherence_modulation(1).event_type(nbin,:)==1;
T1_coherent_PRE = V1_reactivation_coherence_modulation(1).event_type_PRE(nbin,:)==1;

T2_coherent = V1_reactivation_coherence_modulation(2).event_type(nbin,:)==1;
T2_coherent_PRE = V1_reactivation_coherence_modulation(2).event_type_PRE(nbin,:)==1;


figure
subplot(2,2,1)
[N_PRE,x,y]= histcounts2(V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE),...
    V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE)...
    ,0:0.05:1,-1:0.05:1);
% max_N = prctile(reshape(N_PRE,1,[]),99);
% N_PRE(N >=max_N) = max_N;
imagesc(x,y,N_PRE);
set(gca, 'YDir', 'normal');
xlabel('PRE proportion of active neurons')
% ylim([0 1])


subplot(2,2,2)
[N,x,y]= histcounts2(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE),...
    V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE)...
    ,0:0.05:1,-1:0.05:1);
% max_N = prctile(reshape(N,1,[]),99);
% N(N >=max_N) = max_N;

imagesc(x,y,N);
set(gca, 'YDir', 'normal');
xlabel('POST proportion of active neurons')



subplot(2,2,3)
[N,x,y]= histcounts2(V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE),...
    V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE)...
    ,0:0.05:1,0:0.05:1);
% max_N = prctile(reshape(N,1,[]),99);
% N(N >=max_N) = max_N;

imagesc(x,y,N);
set(gca, 'YDir', 'normal');
xlabel('PRE proportion of active neurons')
ylabel('POST proportion of active neurons')

subplot(2,2,3)
imagesc(x,y,N - N_PRE );
set(gca, 'YDir', 'normal');




scatter(V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE),...
    V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE),'filled','MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0')

scatter(V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE),...
    V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE),'filled','MarkerFaceAlpha','0.01','MarkerEdgeAlpha','0')

histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE),'Normalization','probability')
hold on;histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,:)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,:),'Normalization','probability')
hold on;histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent),'Normalization','probability')

hold on;histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,:)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,:),'Normalization','probability')
hold on;histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,:)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,:),'Normalization','probability')

hold on;histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,:)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,:),'Normalization','probability')
hold on;histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,:)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,:),'Normalization','probability')


%%%%%%%%%%%%

fig = figure
fig.Name = 'Increase in proportion of context-selective V1 neurons active following coherent reactivation';
fig.Position = [800 100 1050 720];

subplot(2,2,1)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent&T1_coherent_PRE);
dist2 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent&T1_coherent_PRE);
p = signrank(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent),-1:0.04:1,'Normalization','probability')
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('Coherent Track L reactivation')


subplot(2,2,2)
dist1 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent);
p = signrank(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent),-1:0.04:1,'Normalization','probability')
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('Coherent Track R reactivation')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



subplot(2,2,3)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent);
p = ranksum(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('V1 L neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent);
p = ranksum(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('V1 R neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])





fig = figure
fig.Name = 'Increase in proportion of context-selective V1 neurons active following coherent reactivation';
fig.Position = [800 100 1050 720];

subplot(2,2,1)
index = V1_reactivation_coherence_modulation(1).neurons_active(1,:)>0 & V1_reactivation_coherence_modulation(1).neurons_active(2,:)>0;
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,index & T1_coherent_PRE & T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,index&T1_coherent_PRE & T1_coherent);
% index = V1_reactivation_coherence_modulation(1).neurons_active(1,:)>0 & V1_reactivation_coherence_modulation(1).neurons_active(2,:)>0;
dist2 = V1_reactivation_coherence_modulation(1).neurons_active(2,index &T1_coherent_PRE & T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,index &T1_coherent_PRE & T1_coherent);
p = signrank(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent_PRE & T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent_PRE & T1_coherent),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent_PRE & T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent_PRE & T1_coherent),-1:0.04:1,'Normalization','probability')
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('Coherent Track L reactivation')


subplot(2,2,2)
dist1 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent_PRE & T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent_PRE& T2_coherent);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent_PRE & T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent_PRE& T2_coherent);
p = signrank(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent_PRE& T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent_PRE& T2_coherent),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent_PRE& T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent_PRE& T2_coherent),-1:0.04:1,'Normalization','probability')
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('Coherent Track R reactivation')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



subplot(2,2,3)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent_PRE&T1_coherent);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent_PRE &T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent_PRE&T2_coherent);
p = ranksum(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent&T1_coherent_PRE),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent_PRE&T2_coherent)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent_PRE&T2_coherent),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('V1 L neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent&T1_coherent_PRE);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent&T2_coherent_PRE)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent&T2_coherent_PRE);
p = ranksum(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent&T1_coherent_PRE),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent&T2_coherent_PRE)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent&T2_coherent_PRE),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('V1 R neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)








fig = figure
fig.Name = 'Increase in proportion of context-selective V1 neurons active following incoherent reactivation';
fig.Position = [800 100 1050 720];

subplot(2,2,1)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent==0);
dist2 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent==0);
p = signrank(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent==0),-1:0.04:1,'Normalization','probability')
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('Coherent Track L reactivation')


subplot(2,2,2)
dist1 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent==0);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent==0);
p = signrank(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent==0),-1:0.04:1,'Normalization','probability')
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('Coherent Track R reactivation')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



subplot(2,2,3)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent==0);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent==0);
p = ranksum(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('V1 L neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent==0);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent==0);
p = ranksum(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('V1 R neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])




T1_coherent = V1_reactivation_coherence_modulation(1).event_type(1,:)==1;
T1_coherent_PRE = V1_reactivation_coherence_modulation(1).event_type_PRE(1,:)==1;

T2_coherent = V1_reactivation_coherence_modulation(2).event_type(1,:)==1;
T2_coherent_PRE = V1_reactivation_coherence_modulation(2).event_type_PRE(1,:)==1;

fig = figure
fig.Name = 'Proportion of context-selective V1 neurons active following all reactivations';
fig.Position = [800 100 1050 720];

subplot(2,2,1)
index = V1_reactivation_coherence_modulation(1).neurons_active(1,:)>0 & V1_reactivation_coherence_modulation(1).neurons_active(2,:)>0;
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,index);
dist2 = V1_reactivation_coherence_modulation(1).neurons_active(2,index);
p = signrank(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,index),0:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,index),0:0.04:1,'Normalization','probability')
% xlim([-0.8 0.8])
text(0.4,0.05,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('Track L reactivation')


subplot(2,2,2)
index = V1_reactivation_coherence_modulation(2).neurons_active(1,:)>0 & V1_reactivation_coherence_modulation(2).neurons_active(2,:)>0;
dist1 = V1_reactivation_coherence_modulation(2).neurons_active(1,index);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,index);
p = signrank(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,index),0:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,index),0:0.04:1,'Normalization','probability')
xlim([0 1])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'V1 L','V1 R'},'box','off')
xlabel('proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('Track R reactivation')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



subplot(2,2,3)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent==0);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent==0);
p = ranksum(dist1,dist2,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(1,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(1,T2_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')
title('V1 L neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent==0);
dist2 = V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent==0);
p = ranksum(dist2,dist1,'tail','left');

histogram(V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent==0)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
histogram(V1_reactivation_coherence_modulation(2).neurons_active(2,T2_coherent==0)-V1_reactivation_coherence_modulation(2).neurons_active_PRE(2,T2_coherent==0),-1:0.04:1,'Normalization','probability');hold on;
xlim([-0.8 0.8])
text(0.4,0.1,sprintf('p = %.3e',p))
legend({'Coherent Track L event','Coherent Track R event'},'box','off')
xlabel('Difference in proportion of neurons active following ripple')
ylabel('proportion of ripple events')

title('V1 R neurons')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])






dist1 = V1_reactivation_coherence_modulation(1).V1_bias(:);
dist2 = V1_reactivation_coherence_modulation(2).V1_bias(:);
p = signrank(dist2,dist1,'tail','left');


V1_reactivation_coherence_modulation(1).V1_bias(T1_coherent&T1_coherent_PRE)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_coherent&T1_coherent_PRE)


V1_reactivation_coherence_modulation(1).V1_bias(T1_coherent)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(T1_coherent)
V1_reactivation_coherence_modulation(1).V1_bias(:)-V1_reactivation_coherence_modulation(1).V1_bias_PRE(:)

histogram(V1_reactivation_coherence_modulation(1).V1_bias(:))
hold on
histogram(V1_reactivation_coherence_modulation(2).V1_bias(:))





%%% Coherent events but coherent vs not coherent before ripple
figure;track_id = 1;histogram(V1_reactivation_coherence_modulation(track_id).V1_bias(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==1 & V1_reactivation_coherence_modulation(track_id).event_type_PRE(1,:)==0),-2:0.1:2)
hold on;histogram(V1_reactivation_coherence_modulation(track_id).V1_bias(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==1 & V1_reactivation_coherence_modulation(track_id).event_type_PRE(1,:)==1),-2:0.1:2)

%%% Not Coherent events but coherent vs not coherent before ripple
figure;track_id = 1;histogram(V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type_PRE(1,:)==0 & V1_reactivation_coherence_modulation(track_id).event_type_PRE(1,:)==1),-2:0.1:2)
hold on;histogram(V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==0 & V1_reactivation_coherence_modulation(track_id).event_type_PRE(1,:)==1),-2:0.1:2)

%%% Coherent events -> increase?
figure;track_id = 1;
histogram(V1_reactivation_coherence_modulation(track_id).V1_bias(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==1) -...
    V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==1),-2:0.1:2)

% figure;track_id = 1;
hold on
histogram(V1_reactivation_coherence_modulation(track_id).V1_bias(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==0) -...
    V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==0),-2:0.1:2)


figure;track_id = 1;
hold on
histogram(V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==0),-2:0.1:2)
histogram(V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==1),-2:0.1:2)

% figure;track_id = 1;

histogram(V1_reactivation_coherence_modulation(track_id).V1_bias(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==0) -...
    V1_reactivation_coherence_modulation(track_id).V1_bias_PRE(V1_reactivation_coherence_modulation(track_id).event_type(1,:)==0),-2:0.1:2)











%% CCA analysis
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'));
nbin = 1;
T1_coherent = V1_reactivation_coherence_modulation(1).event_type(nbin,:)==1;
T1_coherent_PRE = V1_reactivation_coherence_modulation(1).event_type_PRE(nbin,:)==1;

T2_coherent = V1_reactivation_coherence_modulation(2).event_type(nbin,:)==1;
T2_coherent_PRE = V1_reactivation_coherence_modulation(2).event_type_PRE(nbin,:)==1;


ripple_powers_all = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];

clear CCA
% warning('off', 'stats:canoncorr:NotFullRank'); % Place it here
for nsession = 1:max(V1_reactivation_coherence_modulation(1).session_id)
    tic

    for track_id = 1:2
        CCA(track_id).V1{nsession} = [];
        CCA(track_id).HPC{nsession} = [];
        CCA(track_id).corr_coef{nsession} = [];


        selected_event_id = find(V1_reactivation_coherence_modulation(track_id).session_id == nsession);
        % ripple_powers = ripple_powers_all(V1_reactivation_coherence_modulation(track_id).event_id(selected_event_id));
        % ripple_powers >

        psth1_cells = {V1_reactivation_coherence_modulation(track_id).PSTH1_zscored{selected_event_id}};
        psth2_cells = {V1_reactivation_coherence_modulation(track_id).PSTH2_zscored{selected_event_id}};
        combined_cells = cellfun(@(a, b) [a; b], psth1_cells, psth2_cells, 'UniformOutput', false);
        temp_3d = cat(3, combined_cells{:});
        psth_v1 = permute(temp_3d, [3, 1, 2]);

        psth_tmp =  {V1_reactivation_coherence_modulation(track_id).HPC_PSTH_zscored{selected_event_id(:)}};
        temp_3d = cat(3, psth_tmp{:});
        psth_hpc = permute(temp_3d, [3, 1, 2]);


        % Time information
        bin_size = bin_centers(2) - bin_centers(1); % Calculate bin size (e.g., 0.01s)

        % Define the time lags to test (e.g., -0.2s to 0.2s in steps of bin_size)
        lag_times = -0.2:bin_size:0.2;
        lag_bins = round(lag_times / bin_size); % Convert time to bin indices
        n_lags = length(lag_bins);
        n_time = size(psth_v1,3);

        % Pre-allocate storage for results
        % Store the maximum canonical correlation (r) for each lag
        max_r_per_lag = nan(size(psth_v1,3),n_lags);

        %%% Time-Resolved Lagged CCA Loop
        cca_map = nan(n_time, n_lags); % Pre-allocate heatmap
        CCA_V1 = nan(n_time, n_lags,size(psth_v1,2),8); % Pre-allocate heatmap
        CCA_HPC = nan(n_time, n_lags,size(psth_hpc,2),8); % Pre-allocate heatmap
        % v1_ncell= size(psth_v1,2);
        % hpc_ncell= size(psth_hpc,2);
        % if v1_ncell>hpc_ncell
        %     CCA_V1 = nan(n_time, n_lags,size(psth_v1,2),10); % Pre-allocate heatmap
        %     CCA_HPC = nan(n_time, n_lags,size(psth_hpc,2),10); % Pre-allocate heatmap
        % else
        %     CCA_V1 = nan(n_time, n_lags,size(psth_v1,2),10); % Pre-allocate heatmap
        %     CCA_HPC = nan(n_time, n_lags,size(psth_hpc,2),10); % Pre-allocate heatmap
        % end

        t_win = 0.05;

        for t = 1:size(psth_v1,3)
            tidx = find(bin_centers<bin_centers(t)+t_win & bin_centers>=bin_centers(t));

            for i_lag = 1:n_lags
                lag = lag_bins(i_lag);
                % Determine the time indices to use for each population based on the lag
                % The lag shifts the HPC data relative to the V1 data.

                % if V1 timedows are within the time windows
                if tidx(1)+lag <=0 | tidx(end)+lag >length(bin_centers)
                    continue
                end

                % Find lagged V1 time indices
                v1_time_idx =  tidx(:)+lag;

                % Find HPC time indices
                hpc_time_idx = tidx(:);

                % 3. Extract and Reshape Data
                % The goal is to reshape the [nEvent x nNeuron x nTime] data into a
                % 2D matrix of [ (nEvent * nTime) x nNeuron ]

                % Extract the relevant time slices for V1 and HPC
                v1_lagged_data = sum(psth_v1(:, :, v1_time_idx),3);
                hpc_lagged_data = sum(psth_hpc(:, :, hpc_time_idx),3);

                % 4. Perform CCA
                % canoncorr requires the number of observations (rows) to be greater
                % than the number of variables (columns). This is typically true here
                % since nEvent * nTime >> nNeuron.
                warnState = warning('off', 'stats:canoncorr:NotFullRank');
                [A, B, r, ~, ~] = canoncorr(v1_lagged_data,hpc_lagged_data);

                % Store the maximum canonical correlation (r(1))
                max_r_per_lag(t,i_lag) = r(1);
                CCA_V1(t,i_lag,:,:) = A(:,1:8);
                CCA_HPC(t,i_lag,:,:) = B(:,1:8);

            end
        end

        CCA(track_id).V1{nsession} = CCA_V1;
        CCA(track_id).HPC{nsession} = CCA_HPC;
        CCA(track_id).corr_coef{nsession} = max_r_per_lag;
        % 
        % for nevent = 1:length(selected_event_id)
        % 
        %     event_id = V1_reactivation_coherence_modulation(track_id).event_id(selected_event_id(nevent));
        % 
        %     PSTH1= V1_reactivation_coherence_modulation(track_id).PSTH1{(selected_event_id(nevent))};
        %     PSTH2 = V1_reactivation_coherence_modulation(track_id).PSTH2{(selected_event_id(nevent))};
        % 
        %     subplot(4,1,1)
        %     imagesc(V1_reactivation_coherence_modulation(track_id).time_bins,1:size(PSTH2,1)+size(PSTH2,1),[PSTH2;PSTH1]);hold on;
        %     yline(size(PSTH2,1)+0.5,'r')
        %     xlim([bin_centers(1) bin_centers(end)])
        %     xlim([-1 1])
        %     xline(0,'r--')
        % 
        %     subplot(4,1,2)
        %     plot(bin_centers,z_bias_V1_KDE(:,event_id));hold on
        %     plot(bin_centers,z_bias_KDE(:,event_id));
        %     xlim([-1 1])
        %     xline(0,'r--')
        %     % plot(V1_reactivation_coherence_modulation(track_id).PSTH1{100})
        % end
    end
    toc
end
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_CCA.mat'),'CCA','-v7.3');

for track_id = 1:2
    figure
    for nsession = 1:22
        CCA_matrix = CCA(track_id).corr_coef{nsession};
        % V1_leading = mean(CCA_matrix(:,1:end/2),2);
        % V1_lagging = mean(CCA_matrix(:,end/2+1:end),2);
        V1_leading = normalize(mean(CCA_matrix(:,1:end/2),2),'range');
        V1_lagging = normalize(mean(CCA_matrix(:,end/2+1:end),2),'range');
        V1_before = normalize(mean(CCA_matrix(bin_centers>-0.2 & bin_centers<0,:)),'range');
        V1_after = normalize(mean(CCA_matrix(bin_centers>0 & bin_centers<0.2,:)),'range');

        % V1_leading(bin_centers>)
        subplot(5,5,nsession)
        % plot(bin_centers,V1_leading);hold on; plot(bin_centers,V1_lagging)
        % xlim([-0.5 0.5])

        % plot(lag_times,V1_before);hold on; plot(lag_times,V1_after)
        % xlim([-0.2 0.2])
        imagesc(CCA(track_id).corr_coef{nsession});
        dist = reshape(CCA(track_id).corr_coef{nsession},1,[]);

        if prctile(dist,20)<1
            clim([prctile(dist,10) prctile(dist,90)]);
            % colorbar
        end
        % 
        colorbar
    end
end

