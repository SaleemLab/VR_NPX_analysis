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
    end
end


for nsession = 1:length(sessions_to_process)

    ripple_powers = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    threshold = prctile(ripple_powers,25);
    % threshold = 0;
    [event_ids_first,event_ids_second] = merge_bilateral_ripple_events(event_id,event_times,0.05);
%     high_ripple_index = find( ripple_powers > threshold);
%     event_ids_first = intersect(high_ripple_index,event_ids_first);


    event_times = event_times(event_ids_first);
    event_id = event_id(event_ids_first);

    %%%%%%%%%% V1
    V1_PSTH=cell(1,2);
    track_difference=cell(1,2);
    for hemi = 1:2
        % cell_index

        cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
            session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi}));

        cell_id = session_clusters_all.spatial_cell_id{nsession}(cell_index);

        track_difference{hemi} = session_clusters_all.mean_FR{nsession}(cell_index,abs(hemi-3)) - session_clusters_all.mean_FR{nsession}(cell_index,hemi);
        V1_SUA_track_difference{hemi} = [V1_SUA_track_difference{hemi};  track_difference{hemi}];

        spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

        spike_id = session_clusters_all.spike_id{nsession}(spike_index);
        % spike_id(spike_id>0)=1;

        ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
            'unit_id',cell_id,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);

        V1_PSTH{hemi} = [zscore(squeeze(ripple_modulation(1).PSTH),0,2) zscore(squeeze(ripple_modulation(2).PSTH),0,2)];
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
    for hemi = 1:2
        
        squeeze(mean(V1_PSTH{hemi}(:,T2_index,:),2));
    end
end

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


load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH.mat'),'V1_SUA_reactivation_PSTH',...
    'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_high.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_SUA_reactivation_PSTH_low.mat'),'V1_SUA_reactivation_PSTH',...
%     'V1_SUA_reactivation_PSTH_odd','V1_SUA_reactivation_PSTH_even');
% V1_SUA_reactivation_PSTH_low = V1_SUA_reactivation_PSTH;
V1_SUA_reactivation_PSTH_high = V1_SUA_reactivation_PSTH;
%%% Colormap
track1_color = [0 0 1]; % blue
track2_color = [1 0 0]; % red

n_half = 80;     % Number of steps from red to white and white to blue
n_white = 60;     % Number of extra white steps to make the white region broader

% Red to white
r2w = [linspace(1,1,n_half)', linspace(0,1,n_half)', linspace(0,1,n_half)'];

% White to blue
w2b = [linspace(1,0,n_half)', linspace(1,0,n_half)', linspace(1,1,n_half)'];

% Extra white section
white = ones(n_white, 3);

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
Fig.Name = 'V1 response to track-selective reactivations (even by odd)';
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
x_bins = ripple_modulation(1).bins;
% bins_selected = x_bins > -0.5 & x_bins < 0.5;
bins_selected = x_bins >-0.2 & x_bins < 0.2;
% figure

subplot(1,4,1)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{1}{2}(:,bins_selected),2));
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
h = imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{1}{1}(index,:));

% xticks
% xticklabels([])
colorbar;colormap((red_white_blue));clim([-0.3 0.3])
xlim([-0.5 0.5])
xline(0,'--')
title('L V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,2)
% [~,index] = sort(V1_SUA_track_difference{1});
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{1}{2}(:,bins_selected),2));
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{1}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.3 0.3])
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
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{2}{1}(index,:));
colorbar;colormap((red_white_blue));clim([-0.3 0.3])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 1 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(1,4,4)
[~,index] = sort(nanmean(V1_SUA_reactivation_PSTH_odd{2}{1}(:,bins_selected),2));
% [~,index] = sort(V1_SUA_track_difference{2});
% [~,index] = sort(V1_SUA_track_difference{1});
x_data = ripple_modulation(1).bins;
y_data = 1:length(index);
imagesc(x_data,y_data,V1_SUA_reactivation_PSTH_even{2}{2}(index,:));
colorbar;colormap((red_white_blue));clim([-0.25 0.25])
% xlim([100 200])
xlim([-0.5 0.5])
xline(0,'--')
title('R V1 Track 2 ripples')
xlabel('Time relative to ripple onset (s)')
ylabel('Cell ID')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])

% 
% Fig = figure;
% % Fig.Name = 'V1 response to track-selective reactivations';
% Fig.Name = 'V1 response to track-selective reactivations (0-0.2s)';
% % Fig.Name = 'V1 response to track-selective reactivations (-0.2-0s)';
% % Fig.Name = 'V1 response to track-selective reactivations (high ripples)';
% % Fig.Name = 'V1 response to track-selective reactivations (low ripples)';
% Fig.Position = [60 60 1280 406];
% % Fig.Position = [60 60 1876 1023];
% Fig.Position = [60 60 1296 953];
% 
% 
% % [595 645 1296 406]
% x_bins = ripple_modulation(1).bins;
% % bins_selected = x_bins > -0.5 & x_bins < 0.5;
% bins_selected = x_bins > 0 & x_bins < 0.2;
% % figure
% 
% subplot(1,4,1)
% % [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% x_data = ripple_modulation(1).bins;
% y_data = 1:length(index);
% h = imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{1}(index,:));
% 
% % xticks
% % xticklabels([])
% colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% xlim([-0.5 0.5])
% xline(0,'--')
% title('L V1 Track 1 ripples')
% xlabel('Time relative to ripple onset (s)')
% ylabel('Cell ID')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% subplot(1,4,2)
% % [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{1}{2}(:,bins_selected),2));
% x_data = ripple_modulation(1).bins;
% y_data = 1:length(index);
% imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{1}{2}(index,:));
% colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% % xlim([100 200])
% xlim([-0.5 0.5])
% xline(0,'--')
% title('L V1 Track 2 ripples')
% xlabel('Time relative to ripple onset (s)')
% ylabel('Cell ID')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% subplot(1,4,3)
% % [~,index] = sort(V1_SUA_track_difference{2});
% % [~,index] = sort(V1_SUA_track_difference{1});
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2));
% x_data = ripple_modulation(1).bins;
% y_data = 1:length(index);
% imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{1}(index,:));
% colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% % xlim([100 200])
% xlim([-0.5 0.5])
% xline(0,'--')
% title('R V1 Track 1 ripples')
% xlabel('Time relative to ripple onset (s)')
% ylabel('Cell ID')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% subplot(1,4,4)
% [~,index] = sort(nanmean(V1_SUA_reactivation_PSTH{2}{1}(:,bins_selected),2));
% % [~,index] = sort(V1_SUA_track_difference{2});
% % [~,index] = sort(V1_SUA_track_difference{1});
% x_data = ripple_modulation(1).bins;
% y_data = 1:length(index);
% imagesc(x_data,y_data,V1_SUA_reactivation_PSTH{2}{2}(index,:));
% colorbar;colormap((red_white_blue));clim([-0.2 0.2])
% % xlim([100 200])
% xlim([-0.5 0.5])
% xline(0,'--')
% title('R V1 Track 2 ripples')
% xlabel('Time relative to ripple onset (s)')
% ylabel('Cell ID')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA'),[])

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

V1_reactivation_modulation_all = struct();
% context_corr_all = struct();
tic
hemispheres = {'L','R'};

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

    % V1_reactivation_coherence_modulation(track_id).PSTH1_PRE = [];

    V1_reactivation_coherence_modulation(track_id).PSTH2 = [];
    V1_reactivation_coherence_modulation(track_id).HPC_PSTH = [];

    % V1_reactivation_coherence_modulation(track_id).PSTH2_PRE = [];

    V1_reactivation_coherence_modulation(track_id).event_id = [];

end


counter = 1;
for nsession = 1:length(sessions_to_process)

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
    V1_PSTH=cell(1,2);
    HPC_PSTH=cell(1,1);
    track_difference=cell(1,2);


    for hemi = 1:2
        % cell_index

        cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3)) > session_clusters_all.mean_FR{nsession}(:,hemi) &...
            session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') &...
            contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        % cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') &...
        %     contains(session_clusters_all.region{nsession},hemispheres{hemi}));

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

    end


    cell_index = (session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'HPC'));
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


    xbins = ripple_modulation(1).bins;
    % tvec = ripple_modulation(1).;
    % Get mean bias for each ripple event in HPC and get within session
    % Track 1 and Track 2 biased ripple events based on HPC bias

    HPC_bias_this_session = z_bias_KDE(:,session_id == sessions_to_process(nsession));
    HPC_bias_this_session = HPC_bias_this_session(:,event_ids_first);
    
    mean_bias = mean(z_bias_KDE(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');
    log_odds_threshold = prctile(mean_bias,[25 75]);

    mean_bias = mean_bias(event_ids_first);

    log_odds_threshold = prctile(mean_bias,[25 75]);
    event_index = [];
    event_index{1} = find(mean_bias > log_odds_threshold(2));
    event_index{2} = find(mean_bias < log_odds_threshold(1));
    

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
        PSTH1_all = [];
        PSTH1_PRE_all = [];
        PSTH2_all = [];
        PSTH2_PRE_all = [];

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
            % PSTH1_PRE_all = [PSTH1_PRE_all {PSTH_PRE1}];
            PSTH2_all = [PSTH2_all {PSTH2}];
            % PSTH2_PRE_all = [PSTH2_PRE_all {PSTH_PRE2}];
            PSTH_all = [PSTH_all {PSTH}];
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

        V1_reactivation_coherence_modulation(track_id).event_type = [V1_reactivation_coherence_modulation(track_id).event_type event_type];

        V1_reactivation_coherence_modulation(track_id).event_type_PRE = [V1_reactivation_coherence_modulation(track_id).event_type_PRE event_type_PRE];

        V1_reactivation_coherence_modulation(track_id).session_id = [V1_reactivation_coherence_modulation(track_id).session_id nsession*ones(size(bias))];

        V1_reactivation_coherence_modulation(track_id).V1_bias = [V1_reactivation_coherence_modulation(track_id).V1_bias bias];

        V1_reactivation_coherence_modulation(track_id).V1_bias_PRE = [V1_reactivation_coherence_modulation(track_id).V1_bias_PRE bias_PRE];

        V1_reactivation_coherence_modulation(track_id).neurons_active = [V1_reactivation_coherence_modulation(track_id).neurons_active neurons_active];

        V1_reactivation_coherence_modulation(track_id).neurons_active_PRE = [V1_reactivation_coherence_modulation(track_id).neurons_active_PRE neurons_active_PRE];


        V1_reactivation_coherence_modulation(track_id).PSTH1 = [V1_reactivation_coherence_modulation(track_id).PSTH1 PSTH1_all];

        % V1_reactivation_coherence_modulation(track_id).PSTH1_PRE = [V1_reactivation_coherence_modulation(track_id).PSTH1_PRE PSTH1_PRE_all];

        V1_reactivation_coherence_modulation(track_id).PSTH2 = [V1_reactivation_coherence_modulation(track_id).PSTH2 PSTH2_all];
        V1_reactivation_coherence_modulation(track_id).HPC_PSTH = [V1_reactivation_coherence_modulation(track_id).HPC_PSTH PSTH_all];

        % V1_reactivation_coherence_modulation(track_id).PSTH2_PRE = [V1_reactivation_coherence_modulation(track_id).PSTH2_PRE PSTH2_PRE_all];

        V1_reactivation_coherence_modulation(track_id).event_id = [V1_reactivation_coherence_modulation(track_id).event_id  event_index_all{track_id}];

    end

    % for hemi = 1:2
    %     % for track_id = 1:2
    %     V1_SUA_reactivation_PSTH{hemi}{1}=[V1_SUA_reactivation_PSTH{hemi}{1}; squeeze(mean(V1_PSTH{hemi}(:,T1_index,:),2))];
    %     V1_SUA_reactivation_PSTH{hemi}{2}=[V1_SUA_reactivation_PSTH{hemi}{2}; squeeze(mean(V1_PSTH{hemi}(:,T2_index,:),2))];
    % end
    % 
    % 
end

V1_reactivation_coherence_modulation(1).time_bins =  xbins(xbins>=-1&xbins<=1);
V1_reactivation_coherence_modulation(2).time_bins =  xbins(xbins>=-1&xbins<=1);

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'),'V1_reactivation_coherence_modulation','-v7.3');



%% CCA analysis
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation coherence SUA','V1_reactivation_coherence_modulation.mat'));
nbin = 1;
T1_coherent = V1_reactivation_coherence_modulation(1).event_type(nbin,:)==1;
T1_coherent_PRE = V1_reactivation_coherence_modulation(1).event_type_PRE(nbin,:)==1;

T2_coherent = V1_reactivation_coherence_modulation(2).event_type(nbin,:)==1;
T2_coherent_PRE = V1_reactivation_coherence_modulation(2).event_type_PRE(nbin,:)==1;


clear CCA
% warning('off', 'stats:canoncorr:NotFullRank'); % Place it here
for nsession = 1:max(V1_reactivation_coherence_modulation(1).session_id)
    tic

    for track_id = 1:2
        CCA(track_id).V1{nsession} = [];
        CCA(track_id).HPC{nsession} = [];
        CCA(track_id).corr_coef{nsession} = [];


        selected_event_id = find(V1_reactivation_coherence_modulation(track_id).session_id == nsession);

        psth1_cells = {V1_reactivation_coherence_modulation(1).PSTH1{selected_event_id}};
        psth2_cells = {V1_reactivation_coherence_modulation(1).PSTH2{selected_event_id}};
        combined_cells = cellfun(@(a, b) [a; b], psth1_cells, psth2_cells, 'UniformOutput', false);
        temp_3d = cat(3, combined_cells{:});
        psth_v1 = permute(temp_3d, [3, 1, 2]);

        psth_tmp =  {V1_reactivation_coherence_modulation(track_id).HPC_PSTH{selected_event_id(:)}};
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
        % V1_leading = normalize(mean(CCA_matrix(:,1:end/2),2),'range');
        % V1_lagging = normalize(mean(CCA_matrix(:,end/2+1:end),2),'range');
        V1_before = normalize(mean(CCA_matrix(bin_centers>-0.2 & bin_centers<0,:)),'range');
        V1_after = normalize(mean(CCA_matrix(bin_centers>0 & bin_centers<0.2,:)),'range');

        % V1_leading(bin_centers>)
        subplot(5,5,nsession)
        % plot(bin_centers,V1_leading);hold on; plot(bin_centers,V1_lagging)
        % xlim([-0.2 0.2])

        plot(lag_times,V1_before);hold on; plot(lag_times,V1_after)
        xlim([-0.2 0.2])
        % imagesc(CCA(track_id).corr_coef{nsession});
        % dist = reshape(CCA(track_id).corr_coef{nsession},1,[]);
        % 
        % if prctile(dist,20)<1
        %     clim([prctile(dist,20) prctile(dist,80)]);
        % end
        % 
        % colorbar
    end
end


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
dist1 = V1_reactivation_coherence_modulation(1).neurons_active(1,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(1,T1_coherent);
dist2 = V1_reactivation_coherence_modulation(1).neurons_active(2,T1_coherent)-V1_reactivation_coherence_modulation(1).neurons_active_PRE(2,T1_coherent);
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


