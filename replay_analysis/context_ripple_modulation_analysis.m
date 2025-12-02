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
% load(fullfile(analysis_folder,'session_clusters_all_POST.mat'))
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




%% Ripple PSTH MUA in V1 and HPC
sessions_to_process = 1:22;
ripple_MUA_PSTH_all = [];
ripple_MSUA_PSTH_all = [];

ripple_V1_MUA_PSTH_all = [];
ripple_V1_MSUA_PSTH_all = [];

% ripple_PSTH_MUA_all = [];
psthBinSize = 0.01;
windows = [-2 2];
hemispheres = {'L','R'};

for nsession = 1:length(sessions_to_process)
    %     ripple_modulation_PSTH_all{nsession} = [];

    %     all_clusters = session_clusters_all.spatial_cell_id{nsession};
    ripple_MUA_PSTH_all{1}{nsession} = [];
    ripple_MUA_PSTH_all{2}{nsession} = [];
    tic
    for hemi = 1:2

        %       plot(unique(session_clusters_all.spike_id{nsession}));hold on;plot(all_clusters)

        %         for nprobe = 1:length(ripples_all)
        event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
        event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];
        %         merge_bilateral_ripple_events(event_id,event_times,0.3);
        remove_id = find(diff(event_times)<0.2)+1;
        event_times(remove_id)=[];
        event_id(remove_id)=[];

        %%%%%%%%%% HPC
        cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'HPC') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
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

        clear cell_PSTH
        for ncell = 1:length(cell_id)
            cell_PSTH(1,ncell,:) = mean(zscore(squeeze(ripple_modulation(1).PSTH(ncell,:,:)),0,2),'omitnan');
            cell_PSTH(2,ncell,:) = mean(zscore(squeeze(ripple_modulation(2).PSTH(ncell,:,:)),0,2),'omitnan');
        end

        ripple_MSUA_PSTH_all{1}{nsession}(hemi,:) = mean(squeeze(cell_PSTH(1,:,:)),'omitnan');
        ripple_MSUA_PSTH_all{2}{nsession}(hemi,:) = mean(squeeze(cell_PSTH(2,:,:)),'omitnan');


        %%%%%%%%%%%% V1
        cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

        spike_id = session_clusters_all.spike_id{nsession}(spike_index);
        spike_id(spike_id>0)=1;

        ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
            'unit_id',1,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);

        ripple_V1_MUA_PSTH_all{1}{nsession}(hemi,:,:) = zscore(squeeze(ripple_modulation(1).PSTH),0,2);
        ripple_V1_MUA_PSTH_all{2}{nsession}(hemi,:,:) = zscore(squeeze(ripple_modulation(2).PSTH),0,2);


        spike_id = session_clusters_all.spike_id{nsession}(spike_index);
%         spike_id(spike_id>0)=1;

        ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
            'unit_id',cell_id,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);

        clear cell_PSTH
        for ncell = 1:length(cell_id)
            cell_PSTH(1,ncell,:) = mean(zscore(squeeze(ripple_modulation(1).PSTH(ncell,:,:)),0,2),'omitnan');
            cell_PSTH(2,ncell,:) = mean(zscore(squeeze(ripple_modulation(2).PSTH(ncell,:,:)),0,2),'omitnan');
        end

        ripple_V1_MSUA_PSTH_all{1}{nsession}(hemi,:) = mean(squeeze(cell_PSTH(1,:,:)),'omitnan');
        ripple_V1_MSUA_PSTH_all{2}{nsession}(hemi,:) = mean(squeeze(cell_PSTH(2,:,:)),'omitnan');

    end
    toc
end

save(fullfile(analysis_folder,'ripple_MUA_PSTH_all_POST.mat'),'ripple_MUA_PSTH_all','ripple_MSUA_PSTH_all','-v7.3')


temp_PSTH1 = [];temp_PSTH2 = [];temp_PSTH =[];
for hemi = 1:2
    temp_PSTH1{hemi}=[];
    temp_PSTH2{hemi}=[];
    for nsession = 1:22
        temp_PSTH1{hemi} = [temp_PSTH1{hemi}; squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(1,:,:))];
        temp_PSTH2{hemi} =  [temp_PSTH2{hemi}; squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(2,:,:))];
    end
end

x = ripple_modulation(1).bins;
temp_PSTH(1,:,:)= [temp_PSTH1{1}; temp_PSTH2{2}];
temp_PSTH(2,:,:) = [temp_PSTH2{1}; temp_PSTH1{2}];


ipsi_lower = prctile(temp_PSTH(1,:,:),25);
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

%%%%%%%%%%%%%%%%

temp_PSTH1 = [];temp_PSTH2 = [];temp_PSTH =[];
temp_HPC_PSTH1 = [];temp_HPC_PSTH2 = [];temp_HPC_PSTH =[];
for nsession = 1:22
    for hemi = 1:2
        temp_PSTH1(hemi,nsession,:) = mean(squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(1,:,:)));
        temp_PSTH2(hemi,nsession,:) =  mean(squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(2,:,:)));
    end

    temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    %     temp_PSTH2(hemi,nsession,:) =  mean(squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(2,:,:)));

    for hemi = 1:2
        temp_HPC_PSTH1(hemi,nsession,:) = mean(squeeze(ripple_MUA_PSTH_all{hemi}{nsession}(1,:,:)));
        temp_HPC_PSTH2(hemi,nsession,:) =  mean(squeeze(ripple_MUA_PSTH_all{hemi}{nsession}(2,:,:)));
    end

    temp_HPC_PSTH(1,nsession,:) = mean([squeeze(temp_HPC_PSTH1(1,nsession,:)) squeeze(temp_HPC_PSTH2(2,nsession,:))],2);
    temp_HPC_PSTH(2,nsession,:) = mean([squeeze(temp_HPC_PSTH1(2,nsession,:)) squeeze(temp_HPC_PSTH2(1,nsession,:))],2);
end

% temp_PSTH1 = [];temp_PSTH2 = [];temp_PSTH =[];
% 
% for nsession = 1:22
%     for hemi = 1:2
%         temp_PSTH1(hemi,nsession,:) = ripple_V1_MSUA_PSTH_all{hemi}{nsession}(1,:) ;
%         temp_PSTH2(hemi,nsession,:) = ripple_V1_MSUA_PSTH_all{hemi}{nsession}(2,:);
%     end
%     temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
%     temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
% end



fig = figure;
fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold)';
fig.Position = [988 440 820 530];


ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(2,2,1)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1 1])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('V1 MUA (z)');
title('V1');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'Ipsi', 'Contra'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-0.5 0.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(2,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1 1])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'Ipsi', 'Contra'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-0.5 0.5])

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])


%% Ripple modulation in V1 and HPC
ripple_modulation_PSTH_all = [];
ripple_PSTH_MUA_all = [];
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


    [event_ids_first,event_ids_second] = merge_bilateral_ripple_events(event_id,event_times,0.05);
    event_times = event_times(event_ids_first);
    event_id = event_id(event_ids_first);

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
    log_odds_threshold = prctile(mean_bias,[25 75]);
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
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');

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
        amplitudes = [amplitudes ripples_all(nprobe).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),ia)];
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



    context_modulation_all.ripple_modulation_percentile{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).ripple_modulation_percentile; ripple_modulation_PSTH_all{nsession}(2).ripple_modulation_percentile]);
    context_modulation_all.modulation_percentile_PRE{nsession} = max([ripple_modulation_PSTH_all{nsession}(1).modulation_percentile_PRE; ripple_modulation_PSTH_all{nsession}(2).modulation_percentile_PRE]);
end
context_modulation_all.timebin = ripple_modulation.bins;

toc

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')

    % subplot(2,2,1)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.POST_ripple_FR(V1_id))
    % subplot(2,2,2)
    % scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
%     % 

%% Plotting context selecitve ripple modulation (low vs high ripples)
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
FR_track_ratio = (FR_track(1,:)-FR_track(2,:))./( FR_track(1,:)+FR_track(2,:));


colorlines = [ ...
    0.2, 0.4, 0.8;    % blue
    0.85, 0.2, 0.2;
];

    % colorlines = [215,25,28;44,123,182]/256;

% colorlines = [ ...
%     0.85, 0.2, 0.2;  % red
%     0.2, 0.4, 0.8    % blue
% ];

nfig = figure;
nfig.Name = 'Context selective FR difference in V1 and HPC';
nfig.Position = [623          88        950         900];
tiledlayout(nfig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
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
[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'V1_R')),FR_diff(contains(regions_all,'V1_L')),'Tail','smaller');
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

[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'V1_R')),FR_diff(contains(regions_all,'V1_L')),'Tail','smaller');
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
% a = FR_track_ratio(1,:);  % Track L firing rate
% b = FR_track_ratio(2,:);  % Track R firing rate
FR_diff = FR_track_ratio;
histogram(FR_diff(contains(regions_all,'V1_R')),[-1:0.05:1],'Normalization','probability','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
histogram(FR_diff(contains(regions_all,'V1_L')),[-1:0.05:1],'Normalization','probability','FaceColor',colorlines(2,:),'EdgeColor','none');

[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'V1_R')),FR_diff(contains(regions_all,'V1_L')),'Tail','smaller');
text(prctile(a,90),0.1,sprintf('p = %.3e',pval))
% ylim([0 0.12])
xline(0,'k--')
xlabel('Track L- Track R FR ratio diff')
ylabel('cum prop of cells')
legend('V1 R','V1 L','box','off')
title('Track L - Track R FR ratio diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal

nexttile
% a = z_FR_track(1,:);  % Track L firing rate
% b = z_FR_track(2,:);  % Track R firing rate
FR_diff = FR_track_ratio;
hold on;
histogram(FR_diff(contains(regions_all,'V1_L')),[-1:0.05:1],'Normalization','cdf','FaceColor',colorlines(2,:),'EdgeColor','none');
histogram(FR_diff(contains(regions_all,'V1_R')),[-1:0.05:1],'Normalization','cdf','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
ylim([0 1])
xline(0,'k--')
xlabel('Track L- Track R FR ratio diff')
ylabel('cum prop of cells')
legend('V1 L','V1 R','box','off')
title('Track L - Track R FR ratio diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal


%%%%%%%%%%% HPC
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

[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'HPC_R')),FR_diff(contains(regions_all,'HPC_L')),'Tail','smaller');
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
[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'HPC_R')),FR_diff(contains(regions_all,'HPC_L')),'Tail','smaller');
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


nexttile
% a = FR_track_ratio(1,:);  % Track L firing rate
% b = FR_track_ratio(2,:);  % Track R firing rate
FR_diff = FR_track_ratio;
histogram(FR_diff(contains(regions_all,'HPC_R')),[-2:0.1:2],'Normalization','probability','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
histogram(FR_diff(contains(regions_all,'HPC_L')),[-2:0.1:2],'Normalization','probability','FaceColor',colorlines(2,:),'EdgeColor','none');

[h,pval,ks2stat] = kstest2(FR_diff(contains(regions_all,'HPC_R')),FR_diff(contains(regions_all,'HPC_L')),'Tail','smaller');
text(prctile(a,90),0.1,sprintf('p = %.3e',pval))
% ylim([0 0.12])
xline(0,'k--')
xlabel('Track L- Track R FR ratio diff')
ylabel('cum prop of cells')
legend('HPC R','HPC L','box','off')
title('Track L - Track R FR ratio diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal

nexttile
% a = z_FR_track(1,:);  % Track L firing rate
% b = z_FR_track(2,:);  % Track R firing rate
FR_diff = FR_track_ratio;
hold on;
histogram(FR_diff(contains(regions_all,'HPC_L')),[-1:0.05:1],'Normalization','cdf','FaceColor',colorlines(2,:),'EdgeColor','none');
histogram(FR_diff(contains(regions_all,'HPC_R')),[-1:0.05:1],'Normalization','cdf','FaceColor',colorlines(1,:),'EdgeColor','none');hold on;
ylim([0 1])
xline(0,'k--')
xlabel('Track L- Track R FR ratio diff')
ylabel('cum prop of cells')
legend('HPC L','HPC R','box','off')
title('Track L - Track R FR ratio diff')
set(gca,'TickDir','out','Box','off','FontSize',12)
% axis equal



%%%%%%%%%%%%%%%%%% context selective ripple modulation based on different
%%%%%%%%%%%%%%%%%% conditions (ripple power or spindle power or SO phase and etc)

ripple_types = {'ALL','LOW','HIGH'};
PSTH_fields = {'PSTH_diff','PSTH_diff_low_ripple','PSTH_diff_high_ripple'};
% conditions = {'Ripple power',''}
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
FR_V1_ratio = normalize(double(FR_track_ratio(contains(regions_all,'V1'))'));
FR_HPC_ratio = normalize(double(FR_track_ratio(contains(regions_all,'HPC'))'));

subjectID_V1 = categorical(session_id_all(contains(regions_all,'V1'))');
subjectID_HPC = categorical(session_id_all(contains(regions_all,'HPC'))');

% Construct tables
tbl_V1 = table(z_FR_V1,FR_V1_ratio,...
    all_vars_V1.POST_ALL, all_vars_V1.POST_LOW, all_vars_V1.POST_HIGH, ...
    all_vars_V1.PRE_ALL, all_vars_V1.PRE_LOW, all_vars_V1.PRE_HIGH, ...
    all_vars_V1.SHIFT_ALL, all_vars_V1.SHIFT_LOW, all_vars_V1.SHIFT_HIGH, ...
    subjectID_V1, ...
    'VariableNames', {'z_FR_track_diff','FR_track_ratio' ...
    'POST_ripple_FR_diff_ALL','POST_ripple_FR_diff_LOW','POST_ripple_FR_diff_HIGH', ...
    'PRE_ripple_FR_diff_ALL','PRE_ripple_FR_diff_LOW','PRE_ripple_FR_diff_HIGH', ...
    'SHIFT_ripple_FR_diff_ALL','SHIFT_ripple_FR_diff_LOW','SHIFT_ripple_FR_diff_HIGH', ...
    'subjectID'});

tbl_HPC = table(z_FR_HPC,FR_HPC_ratio,...
    all_vars_HPC.POST_ALL, all_vars_HPC.POST_LOW, all_vars_HPC.POST_HIGH, ...
    all_vars_HPC.PRE_ALL, all_vars_HPC.PRE_LOW, all_vars_HPC.PRE_HIGH, ...
    all_vars_HPC.SHIFT_ALL, all_vars_HPC.SHIFT_LOW, all_vars_HPC.SHIFT_HIGH, ...
    subjectID_HPC, ...
    'VariableNames', {'z_FR_track_diff','FR_track_ratio', ...
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


%%%%%%%%%%% ALL ripples

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression ratio difference';
nfig.Position = [   842   345   954   578];




subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
% is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0;
% is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_ALL ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.45 0.45])
xlabel('Track FR ratio diff')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_ALL ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('Track FR ratio diff')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_ALL ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('Track FR ratio diff')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');
% is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
% is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_ALL ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-1 1])
xlabel('Track FR ratio diff')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_ALL ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.7 0.7])
xlabel('Track FR ratio diff')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_ALL ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('Track FR ratio diff')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)




%%%%%%%%%%% ALL ripples

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (Track prefering)';
nfig.Position = [   842   345   954   578];




subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0;

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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
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




%%%%%%%%%%%%% low ripple ratio difference

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (low ripples) ratio difference';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
% is_V1R = contains(regions_all, 'V1')&FR_track_ratio>0;
% is_V1L = contains(regions_all, 'V1')&FR_track_ratio<0;


scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';
X = double(FR_track_ratio(contains(regions_all,'V1')))';
% FR_track_ratio = FR_track(1,:)-FR_track(2,:)./( FR_track(1,:)+FR_track(2,:));
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
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_LOW ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.45 0.45])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_LOW ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_LOW ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');

% is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
% is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_LOW ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-1 1])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_LOW ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.7 0.7])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_LOW ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)




%%%%%%%%%%%%% Low ripple


all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (low ripples Track prefering)';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');

is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0;

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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;

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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
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
% 
% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
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
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples) ratio difference';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
% is_V1R = contains(regions_all, 'V1')&FR_track_ratio>0;
% is_V1L = contains(regions_all, 'V1')&FR_track_ratio<0;


scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';
X = double(FR_track_ratio(contains(regions_all,'V1')))';
% FR_track_ratio = FR_track(1,:)-FR_track(2,:)./( FR_track(1,:)+FR_track(2,:));
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
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_HIGH ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.45 0.45])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_HIGH ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(FR_track_ratio(contains(regions_all,'V1')))';
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
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_HIGH ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

is_V1R = contains(regions_all, 'HPC_R');
is_V1L = contains(regions_all, 'HPC_L');

% is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
% is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_HIGH ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-1 1])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_HIGH ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-1 1])
ylim([-0.7 0.7])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(FR_track_ratio(contains(regions_all,'HPC')))';
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
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_HIGH ~ FR_track_ratio + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-1 1])
ylim([-0.4 0.4])
xlabel('FR_track_ratio')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)




%%%%%%%%%%%%% High ripple

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples)';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

is_V1R = contains(regions_all, 'V1_R');
is_V1L = contains(regions_all, 'V1_L');
% is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0;
% is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0;


scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
% X = double(FR_track_ratio(contains(regions_all,'V1')))';
% FR_track_ratio = FR_track(1,:)-FR_track(2,:)./( FR_track(1,:)+FR_track(2,:));
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
% is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;

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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
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






%%%%%%%%%%%%% High ripple

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples Track prefering)';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0;


scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
% X = double(FR_track_ratio(contains(regions_all,'V1')))';
% FR_track_ratio = FR_track(1,:)-FR_track(2,:)./( FR_track(1,:)+FR_track(2,:));
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
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
% 
% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;

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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
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
    xlim([-0.3 0.3])
    text(0.3,0.05,sprintf('p = %.3e',pval))
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ylim([0 1])
    xline(0,'k--','LineWidth',2)
    xlabel('Track L - Track R Ripple FR diff (z)')
    ylabel('prop of cells')
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
    xlim([-0.3 0.3])
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



%%%%%%% Plot others
ripple_types = {'ALL','LOW','HIGH'};
PSTH_fields = {'PSTH_diff_low_spindle','PSTH_diff_high_spindle'};
ripple_name = 'spindle power';
plot_context_selective_ripple_modulations(context_modulation_all,PSTH_fields,ripple_types,ripple_name)




%%
%%%%%%%%% Cell types (putative low FR and high FR neurons)


%%%%%%%%%%%%% High ripple

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples pyramidal cell)';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';
% 
% is_V1R = contains(regions_all, 'V1')&all_cell_types'~=2&all_mean_FR'<=10;
% is_V1L = contains(regions_all, 'V1')&all_cell_types'~=2&all_mean_FR'<=10;
% is_V1 = contains(regions_all, 'V1')&all_cell_types'~=2&all_mean_FR'<=5;
% is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_cell_types'~=2&all_mean_FR'<=5;
% is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_cell_types'~=2&all_mean_FR'<=5;
is_V1 = contains(regions_all, 'V1')&all_mean_FR'<=5;
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_mean_FR'<=5;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_mean_FR'<=5;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(is_V1))';
% X = double(FR_track_ratio(contains(regions_all,'V1')))';
% FR_track_ratio = FR_track(1,:)-FR_track(2,:)./( FR_track(1,:)+FR_track(2,:));
Y = double(POST_ripple_FR_diff(is_V1));
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

text(prctile(X,99), prctile(Y,99), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.45 0.45])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(z_FR_track_diff(is_V1))';
Y = double(PRE_ripple_FR_diff(is_V1));

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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(is_V1))';
Y = double(shifted_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

% is_V1 = contains(regions_all, 'HPC')&all_cell_types'~=2&all_mean_FR'<=5;
% is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff<0&all_cell_types'~=2&all_mean_FR'<=5;
% is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff>0&all_cell_types'~=2&all_mean_FR'<=5;
is_V1 = contains(regions_all, 'HPC')&all_mean_FR'<=5;
is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0&all_mean_FR'<=5;
is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0&all_mean_FR'<=5;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(is_V1))';
Y = double(POST_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(z_FR_track_diff(is_V1))';
Y = double(PRE_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(is_V1))';
Y = double(shifted_ripple_FR_diff(is_V1));
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






nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples interneruon cell)';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';
% 
% is_V1R = contains(regions_all, 'V1')&all_cell_types'~=2&all_mean_FR'<=10;
% is_V1L = contains(regions_all, 'V1')&all_cell_types'~=2&all_mean_FR'<=10;
% is_V1 = contains(regions_all, 'V1')&all_cell_types'>1&all_mean_FR'>5;
% is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_cell_types'>1&all_mean_FR'>5;
% is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_cell_types'>1&all_mean_FR'>5;
is_V1 = contains(regions_all, 'V1')&all_mean_FR'>5;
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_mean_FR'>5;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_mean_FR'>5;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(is_V1))';
% X = double(FR_track_ratio(contains(regions_all,'V1')))';
% FR_track_ratio = FR_track(1,:)-FR_track(2,:)./( FR_track(1,:)+FR_track(2,:));
Y = double(POST_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(z_FR_track_diff(is_V1))';
Y = double(PRE_ripple_FR_diff(is_V1));

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

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(is_V1))';
Y = double(shifted_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

is_V1 = contains(regions_all, 'HPC')&all_mean_FR'>5;
is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0&all_mean_FR'>5;
is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0&all_mean_FR'>5;
% is_V1 = contains(regions_all, 'HPC')&all_cell_types'>1&all_mean_FR'>5;
% is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0&all_cell_types'>1&all_mean_FR'>5;
% is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0&all_cell_types'>1&all_mean_FR'>5;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(is_V1))';
Y = double(POST_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(z_FR_track_diff(is_V1))';
Y = double(PRE_ripple_FR_diff(is_V1));
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

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(is_V1))';
Y = double(shifted_ripple_FR_diff(is_V1));
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








%%%% KS (pyramidal cells)


%%%%%%%%%%% kstest2
ks2stat_all = nan(1000,10);

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 (high ripples pyramidal) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
% POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
% PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% % shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');
is_V1 = contains(regions_all, 'V1')&all_cell_types'~=2&all_mean_FR'<=5;
% is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_cell_types'~=2&all_mean_FR'<=5;
% is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_cell_types'~=2&all_mean_FR'<=5;
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_mean_FR'<=5;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_mean_FR'<=5;
% timewindows = -1:0.2:1;
timewindows = -1:0.2:1;
for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(is_V1R),POST_ripple_FR_diff(is_V1L),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(is_V1L),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(is_V1R),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
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
nfig.Name = 'Context selective ripple modulation in V1 (low ripples pyramidal) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
% POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
% PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% % shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');

timewindows = -1:0.2:1;

for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(is_V1R),POST_ripple_FR_diff(is_V1L),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(is_V1L),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(is_V1R),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
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

% 
% 
% %%%%%%%%%%% KS max difference bootstrap
% bin_width = 0.1;    % 50 ms
% step_size = 0.01;   % 10 ms
% t_start = -0.5;
% t_end = 0.5;
% 
% % Bin centers
% bin_centers = t_start:step_size:t_end;
% 
% % Bin edges (each bin spans 100 ms centered at bin_centers)
% timewindows = [bin_centers - bin_width/2; bin_centers + bin_width/2];
% ks2stat_all = nan(1000,length(timewindows)-1);
% 
% %%% Low ripples
% all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
% for nbin = 1:length(timewindows)-1
%     POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(1,nbin)&context_modulation_all.timebin<timewindows(2,nbin)),2,'omitnan');
% 
%     T1_FR_dff = POST_ripple_FR_diff(is_V1R);
%     T2_FR_dff = POST_ripple_FR_diff(is_V1L);
%     num_neurons = min([length(T1_FR_dff) length(T2_FR_dff)]);
% 
% 
%     parfor iBoot = 1:1000
%         temp1 = [];temp2 = [];
%         s = RandStream('philox4x32_10', 'Seed', iBoot);
%         index1 = randsample(s, 1:length(T1_FR_dff), num_neurons, true);
%         temp1 = T1_FR_dff(index1);
%         index2 = randsample(s, 1:length(T2_FR_dff), num_neurons, true);
%         temp2 = T2_FR_dff(index2);
%         [h,pval,ks2stat_boot] = kstest2(temp1,temp2);
%         ks2stat_all(iBoot,nbin) = ks2stat_boot;
%     end
% end
% ks2stat_low_ripples = ks2stat_all;
% 
% %%%% High ripples
% all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
% ks2stat_all = nan(1000,length(timewindows)-1);
% 
% for nbin = 1:length(timewindows)-1
%     POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(1,nbin)&context_modulation_all.timebin<timewindows(2,nbin)),2,'omitnan');
% 
%     T1_FR_dff = POST_ripple_FR_diff(is_V1R);
%     T2_FR_dff = POST_ripple_FR_diff(is_V1L);
%     num_neurons = min([length(T1_FR_dff) length(T2_FR_dff)]);
% 
% 
%     parfor iBoot = 1:1000
%         temp1 = [];temp2 = [];
%         s = RandStream('philox4x32_10', 'Seed', iBoot);
%         index1 = randsample(s, 1:length(T1_FR_dff), num_neurons, true);
%         temp1 = T1_FR_dff(index1);
%         index2 = randsample(s, 1:length(T2_FR_dff), num_neurons, true);
%         temp2 = T2_FR_dff(index2);
%         [h,pval,ks2stat_boot] = kstest2(temp1,temp2);
%         ks2stat_all(iBoot,nbin) = ks2stat_boot;
%     end
% end
% ks2stat_high_ripples = ks2stat_all;
% 
% 
% 
% nfig = figure;
% nfig.Name = 'Context selective ripple modulation in V1 KS difference time series (high vs low ripple pyramidal)';
% nfig.Position = [ 1150         322         363         282];
% timewindows = [bin_centers - bin_width/2; bin_centers + bin_width/2];
% % x = timewindows(1:end-1)+mean(diff(timewindows))/2;
% x = mean(timewindows);
% x = x(1:end-1)+mean(diff(x))/2;
% colour_lines = [158,202,225;33,113,181]/256;% two blue 
% 
% nexttile
% % x = timewindows(1:end-1)+mean(diff(timewindows))/2;
% plot(x,mean(ks2stat_low_ripples),'Color',colour_lines(1,:));hold on;
% ci_low  = prctile(ks2stat_low_ripples, 2.5, 1);
% ci_high = prctile(ks2stat_low_ripples, 97.5, 1);
% 
% % Fill 95% CI
% F(1) = fill([x, fliplr(x)], [ci_low, fliplr(ci_high)], ...
%     colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% 
% plot(x,mean(ks2stat_high_ripples),'Color',colour_lines(end,:));hold on;
% ci_low  = prctile(ks2stat_high_ripples, 2.5, 1);
% ci_high = prctile(ks2stat_high_ripples, 97.5, 1);
% 
% % Fill 95% CI
% F(2)=fill([x, fliplr(x)], [ci_low, fliplr(ci_high)], ...
%     colour_lines(end,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% set(gca,'TickDir','out','Box','off','FontSize',12)
% xticks([-0.5 -0.25 0 0.25 0.5])
% xline(0,'r')
% 
% legend(F(1:2),{'low ripples','high ripples'},'box','off')
% xlabel('Time relative to ripple onset (s)')
% ylabel('Maximum empirical cumulative distribution difference');



%%%% Loop though bins (interneruon)


%%%%%%%%%%% kstest2
ks2stat_all = nan(1000,10);

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 (high ripples interneruon) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
% POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
% PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% % shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');
is_V1 = contains(regions_all, 'V1')&all_cell_types'>1&all_mean_FR'>5;
% is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_cell_types'>1&all_mean_FR'>5;
% is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_cell_types'>1&all_mean_FR'>5;
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0&all_mean_FR'>5;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0&all_mean_FR'>5;

% timewindows = -1:0.2:1;
timewindows = -1:0.2:1;
for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(is_V1R),POST_ripple_FR_diff(is_V1L),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(is_V1L),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(is_V1R),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
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
nfig.Name = 'Context selective ripple modulation in V1 (low ripples interneruon) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
% POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
% PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% % shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');

timewindows = -1:0.2:1;

for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(is_V1R),POST_ripple_FR_diff(is_V1L),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(is_V1L),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(is_V1R),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
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


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])


%% 