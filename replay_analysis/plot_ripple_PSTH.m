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
        [~,index] = sort(event_times);
        event_times = event_times(index);
        event_id = event_id(index);

        remove_id = find(diff(event_times)<0.3)+1;
        event_times(remove_id)=[];
        event_id(remove_id)=[];

        %%%%%%%%%% HPC
        %         cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'HPC') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        cell_id = session_clusters_all.spatial_cell_id{nsession}(contains(session_clusters_all.region{nsession},'HPC') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
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
        % clear cell_PSTH
        % for ncell = 1:length(cell_id)
        %     cell_PSTH(1,ncell,:) = squeeze(mean(ripple_modulation(1).PSTH_zscored(ncell,:,:),'omitnan'));
        %     cell_PSTH(2,ncell,:) = squeeze(mean(ripple_modulation(2).PSTH_zscored(ncell,:,:),'omitnan'));
        % end
        ripple_MSUA_PSTH_all{1}{nsession}(hemi,:) = mean(squeeze(cell_PSTH(1,:,:)),'omitnan');
        ripple_MSUA_PSTH_all{2}{nsession}(hemi,:) = mean(squeeze(cell_PSTH(2,:,:)),'omitnan');


        %%%%%%%%%%%% V1
        cell_id = session_clusters_all.spatial_cell_id{nsession}(contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
%         cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
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

% save(fullfile(analysis_folder,'ripple_MUA_PSTH_all_POST_100ms_interval.mat'),'ripple_MUA_PSTH_all','ripple_MSUA_PSTH_all','-v7.3')
% save(fullfile(analysis_folder,'ripple_V1_MUA_PSTH_all_POST_100ms_interval.mat'),'ripple_V1_MUA_PSTH_all','ripple_V1_MSUA_PSTH_all','-v7.3')
save(fullfile(analysis_folder,'ripple_MUA_PSTH_all_POST_300ms_interval.mat'),'ripple_MUA_PSTH_all','ripple_MSUA_PSTH_all','-v7.3')
save(fullfile(analysis_folder,'ripple_V1_MUA_PSTH_all_POST_300ms_interval.mat'),'ripple_V1_MUA_PSTH_all','ripple_V1_MSUA_PSTH_all','-v7.3')

% load(fullfile(analysis_folder,'ripple_MUA_PSTH_all_POST.mat'),'ripple_MUA_PSTH_all','ripple_MSUA_PSTH_all')
% load(fullfile(analysis_folder,'ripple_V1_MUA_PSTH_all_POST.mat'),'ripple_V1_MUA_PSTH_all','ripple_V1_MSUA_PSTH_all')

% ripple_V1_MSUA_PSTH_all
temp_PSTH1 = [];temp_PSTH2 = [];temp_PSTH =[];
for hemi = 1:2
    temp_PSTH1{hemi}=[];
    temp_PSTH2{hemi}=[];
    for nsession = 1:22
        temp_PSTH1{hemi} = [temp_PSTH1{hemi}; squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(1,:,:))];
        temp_PSTH2{hemi} =  [temp_PSTH2{hemi}; squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(2,:,:))];
    end
end

windows = [-2 2];
psthBinSize = 0.01;
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;
% x = ripple_modulation(1).bins;
temp_PSTH(1,:,:)= [temp_PSTH1{1}; temp_PSTH1{2}];
temp_PSTH(2,:,:) = [temp_PSTH2{1}; temp_PSTH2{2}];

% temp_PSTH(1,:,:)= [temp_PSTH1{1}; temp_PSTH2{2}];
% temp_PSTH(2,:,:) = [temp_PSTH2{1}; temp_PSTH1{2}];
% 

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

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH1(2,nsession,:))],2);
    temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH2(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    
    %     temp_PSTH2(hemi,nsession,:) =  mean(squeeze(ripple_V1_MUA_PSTH_all{hemi}{nsession}(2,:,:)));

    for hemi = 1:2
        temp_HPC_PSTH1(hemi,nsession,:) = mean(squeeze(ripple_MUA_PSTH_all{hemi}{nsession}(1,:,:)));
        temp_HPC_PSTH2(hemi,nsession,:) =  mean(squeeze(ripple_MUA_PSTH_all{hemi}{nsession}(2,:,:)));
    end

%     temp_HPC_PSTH(1,nsession,:) = mean([squeeze(temp_HPC_PSTH1(1,nsession,:)) squeeze(temp_HPC_PSTH2(2,nsession,:))],2);
%     temp_HPC_PSTH(2,nsession,:) = mean([squeeze(temp_HPC_PSTH1(2,nsession,:)) squeeze(temp_HPC_PSTH2(1,nsession,:))],2);
    temp_HPC_PSTH(1,nsession,:) = mean([squeeze(temp_HPC_PSTH1(1,nsession,:)) squeeze(temp_HPC_PSTH2(1,nsession,:))],2);
    temp_HPC_PSTH(2,nsession,:) = mean([squeeze(temp_HPC_PSTH1(2,nsession,:)) squeeze(temp_HPC_PSTH2(2,nsession,:))],2);
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
figure
hold on
plot(squeeze(mean(temp_PSTH(:,:,:),1))','b')
% plot(squeeze(temp_PSTH(1,:,:))','b')

plot(squeeze(mean(mean(temp_PSTH(:,:,:),1)))','r')


fig = figure;
fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold)';
fig.Position = [988 440 820 530];

% sessions_to_process() = []
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
xlim([-1.5 1.5])

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
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'Ipsi', 'Contra'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])


%% Ripple PSTH MUA in V1 and HPC (ripple features)


clear ripple_PSTH
ripple_PSTH.V1_MUA = [];
ripple_PSTH.HPC_MUA = [];

% ripple_PSTH_MUA_all = [];
psthBinSize = 0.01;
windows = [-2 2];
hemispheres = {'L','R'};
% 
% load(fullfile(analysis_folder,'ripples_all_best_V1_SO_POST.mat'))
% ripples_all_SO = ripples_all;

load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_best_POST.mat'))


% 
% 
% % ripples_all
% for nsession = 1:length(sessions_to_process)
%     %     ripple_modulation_PSTH_all{nsession} = [];
%     tic
% 
%     event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
%     event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];
% 
%     % [~,index] = sort(event_times);
%     % event_times = event_times(index);
%     % event_id = event_id(index);
%     % spindle_power = spindle_power(:,index);
%     % ripple_power = ripple_power(index);
%     % SO_phase = SO_phase(:,index);
%     %
%     % event_id(event_id>0) = 1;
%     % remove_id = find(diff(event_times)<0.3)+1;
%     % event_times(remove_id)=[];
%     % event_id(remove_id)=[];
%     % spindle_power(:,remove_id)=[];
%     % ripple_power(remove_id)=[];
%     % SO_phase(:,remove_id)=[];
% 
%     HPC_MUA = [];    V1_MUA = [];
%     ripple_PSTH(1).ripple_hemi{nsession} = event_id;
%     %%% HPC
%     cell_id = session_clusters_all.spatial_cell_id{nsession}(contains(session_clusters_all.region{nsession},'HPC'));
%     %         cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
%     spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);
% 
%     spike_id = session_clusters_all.spike_id{nsession}(spike_index);
%     spike_id(spike_id>0)=1;
% 
%     ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
%         'unit_id',1,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);
% 
%     ripple_PSTH(1).HPC_MUA{nsession} = zscore(squeeze(ripple_modulation(1).PSTH),0,2);
%     ripple_PSTH(2).HPC_MUA{nsession} = zscore(squeeze(ripple_modulation(2).PSTH),0,2);
% 
%     %%% V1
%     for hemi = 1:2
%         cell_id = session_clusters_all.spatial_cell_id{nsession}(contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
%         %         cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
%         spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);
% 
%         spike_id = session_clusters_all.spike_id{nsession}(spike_index);
%         spike_id(spike_id>0)=1;
% 
%         ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
%             'unit_id',1,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);
% 
%         ripple_PSTH(1).V1_MUA{nsession}(hemi,:,:) = zscore(squeeze(ripple_modulation(1).PSTH),0,2);
%         ripple_PSTH(2).V1_MUA{nsession}(hemi,:,:) = zscore(squeeze(ripple_modulation(2).PSTH),0,2);
%     end
% end
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_PSTH_POST.mat'),'ripple_PSTH')



sessions_to_process = 1:22;
V1_PSTH_spindle_power = [];
V1_PSTH_ripple_power = [];
V1_PSTH_SO = [];

contra_V1_PSTH_SO = [];
contra_V1_PSTH_spindle_power = [];

ipsi_contra_diff_V1_PSTH_SO = [];
ipsi_contra_diff_V1_PSTH_spindle_power = [];

HPC_PSTH_spindle_power = [];
HPC_PSTH_ripple_power = [];
HPC_PSTH_SO = [];

% contra_HPC_PSTH_SO = [];
% contra_HPC_PSTH_spindle_power = [];
% ripples_all
for nsession = 1:length(sessions_to_process)
    %     ripple_modulation_PSTH_all{nsession} = [];
    tic
    %     all_clusters = session_clusters_all.spatial_cell_id{nsession};


    %       plot(unique(session_clusters_all.spike_id{nsession}));hold on;plot(all_clusters)

    %         for nprobe = 1:length(ripples_all)
    temp = [];
    for n = 1:2
        session_event_index = find(ripples_all(n).session_count == nsession);
        [~,ia,~] = intersect(session_event_index,find(ripples_all(n).session_count == nsession&ripples_all(n).SWS_index==1));
        temp{n} = ripples_all(n).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),ia);
    end
    spindle_power = [temp{1} temp{2}];

        temp = [];
    for n = 1:2
        session_event_index = find(ripples_all(n).session_count == nsession);
        [~,ia,~] = intersect(session_event_index,find(ripples_all(n).session_count == nsession&ripples_all(n).SWS_index==1));
        temp{n} = ripples_all(n).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),ia);
    end
    SO_phase = [temp{1} temp{2}];

    ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];

    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); 2*ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    [~,index] = sort(event_times);
    event_times = event_times(index);
    event_id = event_id(index);
    spindle_power = spindle_power(:,index);
    ripple_power = ripple_power(index);
    SO_phase = SO_phase(:,index);

    event_id(event_id>0) = 1;
    remove_id = find(diff(event_times)<0.3)+1;
    event_times(remove_id)=[];
    event_id(remove_id)=[];
    spindle_power(:,remove_id)=[];
    ripple_power(remove_id)=[];
    SO_phase(:,remove_id)=[];

    HPC_MUA = [];    V1_MUA = [];
    ripple_PSTH.ripple_hemi{nsession} = event_id;
    %%% HPC
    cell_id = session_clusters_all.spatial_cell_id{nsession}(contains(session_clusters_all.region{nsession},'HPC'));
    %         cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
    spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

    spike_id = session_clusters_all.spike_id{nsession}(spike_index);
    spike_id(spike_id>0)=1;

    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
        'unit_id',1,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);


    HPC_MUA = [zscore(squeeze(ripple_modulation(1).PSTH),0,2) ];
    % HPC_MUA = squeeze(ripple_modulation(1).PSTH);

    %%% V1
    for hemi = 1:2
        cell_id = session_clusters_all.spatial_cell_id{nsession}(contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        %         cell_id = session_clusters_all.spatial_cell_id{nsession}(session_clusters_all.mean_FR{nsession}(:,abs(hemi-3))<=50 & contains(session_clusters_all.region{nsession},'V1') & contains(session_clusters_all.region{nsession},hemispheres{hemi}));
        spike_index = ismember(session_clusters_all.spike_id{nsession},cell_id);

        spike_id = session_clusters_all.spike_id{nsession}(spike_index);
        spike_id(spike_id>0)=1;

        ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession}(spike_index),spike_id,windows,psthBinSize,...
            'unit_id',1,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);
        V1_MUA{hemi} = zscore(squeeze(ripple_modulation(1).PSTH),0,2);
        % V1_MUA{hemi} = squeeze(ripple_modulation(1).PSTH);
    end

    temp =[];%ipsi
    %%% ipsilateral spindle power
    for hemi = 1:2
        % thresholds = prctile(ripple_power,[25 75]);
        % low_events = ripple_power < thresholds(2);
        % high_events = ripple_power > thresholds(2);
        %
        thresholds = prctile(spindle_power(hemi,:),[25 75]);
        low_events = spindle_power(hemi,:) < thresholds(1);
        high_events = spindle_power(hemi,:) > thresholds(2);

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:));
    end

    V1_PSTH_spindle_power(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    V1_PSTH_spindle_power(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);

    temp =[];%contralaterals
    %%% contralateral spindle power
    for hemi = 1:2
        thresholds = prctile(spindle_power(abs(hemi-3),:),[25 75]); % contra
        low_events = spindle_power(abs(hemi-3),:) < thresholds(1);
        high_events = spindle_power(abs(hemi-3),:) > thresholds(2);

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:));
    end

    contra_V1_PSTH_spindle_power(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    contra_V1_PSTH_spindle_power(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);


    temp =[];%ipsi-contra
    %%% ipsi-contra spindle power
    for hemi = 1:2
        thresholds = prctile(spindle_power(hemi,:),[25 75]);
        low_events = spindle_power(hemi,:) < thresholds(1);
        high_events = spindle_power(hemi,:) > thresholds(2);

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:)-V1_MUA{abs(hemi-3)}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:)-V1_MUA{abs(hemi-3)}(high_events,:));
    end
    ipsi_contra_diff_V1_PSTH_spindle_power(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    ipsi_contra_diff_V1_PSTH_spindle_power(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);

    temp =[];
    %%% ipsilateral SO phase
    for hemi = 1:2
        % thresholds = prctile(ripple_power,[25 75]);
        % low_events = ripple_power < thresholds(2);
        % high_events = ripple_power > thresholds(2);
        %
        % thresholds = prctile(SO_phase(hemi,:),[25 75]);
        low_events = SO_phase(hemi,:) > -pi/2 & SO_phase(hemi,:) < pi/2;
        high_events = SO_phase(hemi,:) > pi/2 | SO_phase(hemi,:) < -pi/2;

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:));
    end

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    V1_PSTH_SO(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    V1_PSTH_SO(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);

    temp=[];
    %%% contralateral SO phase
    for hemi = 1:2
        % thresholds = prctile(ripple_power,[25 75]);
        % low_events = ripple_power < thresholds(2);
        % high_events = ripple_power > thresholds(2);
        %
        % thresholds = prctile(SO_phase(hemi,:),[25 75]);
        low_events = SO_phase(abs(hemi-3),:) > -pi/2 & SO_phase(abs(hemi-3),:) < pi/2;
        high_events = SO_phase(abs(hemi-3),:) > pi/2 | SO_phase(abs(hemi-3),:) < -pi/2;

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:));
    end

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    contra_V1_PSTH_SO(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    contra_V1_PSTH_SO(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);


    temp =[];
    %%% ipsi-contra SO phase
    for hemi = 1:2
        % thresholds = prctile(ripple_power,[25 75]);
        % low_events = ripple_power < thresholds(2);
        % high_events = ripple_power > thresholds(2);
        %
        % thresholds = prctile(SO_phase(hemi,:),[25 75]);
        low_events = SO_phase(hemi,:) > -pi/2 & SO_phase(hemi,:) < pi/2;
        high_events = SO_phase(hemi,:) > pi/2 | SO_phase(hemi,:) < -pi/2;

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:)-V1_MUA{abs(hemi-3)}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:)-V1_MUA{abs(hemi-3)}(high_events,:));
    end

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    ipsi_contra_diff_V1_PSTH_SO(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    ipsi_contra_diff_V1_PSTH_SO(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);

    temp =[];
    %%% power
    for hemi = 1:2
        thresholds = prctile(ripple_power,[25 75]);
        low_events = ripple_power < thresholds(2);
        high_events = ripple_power > thresholds(2);

        temp{hemi}(1,:) = mean(V1_MUA{hemi}(low_events,:));
        temp{hemi}(2,:) = mean(V1_MUA{hemi}(high_events,:));
    end

    V1_PSTH_ripple_power(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    V1_PSTH_ripple_power(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);


    temp =[];
    %%% spindle power
    for hemi = 1:2
        % thresholds = prctile(ripple_power,[25 75]);
        % low_events = ripple_power < thresholds(2);
        % high_events = ripple_power > thresholds(2);
        %
        thresholds = prctile(spindle_power(hemi,:),[25 75]);
        low_events = spindle_power(hemi,:) < thresholds(1);
        high_events = spindle_power(hemi,:) > thresholds(2);

        temp{hemi}(1,:) = mean(HPC_MUA(low_events,:));
        temp{hemi}(2,:) = mean(HPC_MUA(high_events,:));
    end

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    HPC_PSTH_spindle_power(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    HPC_PSTH_spindle_power(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);

    temp =[];
    %%% SO phase
    for hemi = 1:2
        % thresholds = prctile(ripple_power,[25 75]);
        % low_events = ripple_power < thresholds(2);
        % high_events = ripple_power > thresholds(2);
        %
        low_events = SO_phase(hemi,:) > -pi/2 & SO_phase(hemi,:) < pi/2;
        high_events = SO_phase(hemi,:) > pi/2 | SO_phase(hemi,:) < -pi/2;

        temp{hemi}(1,:) = mean(HPC_MUA(low_events,:));
        temp{hemi}(2,:) = mean(HPC_MUA(high_events,:));
    end

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    HPC_PSTH_SO(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    HPC_PSTH_SO(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);


    temp =[];
    %%% ripple power
    for hemi = 1:2
        thresholds = prctile(ripple_power,[25 75]);
        low_events = ripple_power < thresholds(2);
        high_events = ripple_power > thresholds(2);

        temp{hemi}(1,:) = mean(HPC_MUA(low_events,:));
        temp{hemi}(2,:) = mean(HPC_MUA(high_events,:));
    end

    % temp_PSTH(1,nsession,:) = mean([squeeze(temp_PSTH1(1,nsession,:)) squeeze(temp_PSTH2(2,nsession,:))],2);
    % temp_PSTH(2,nsession,:) = mean([squeeze(temp_PSTH1(2,nsession,:)) squeeze(temp_PSTH2(1,nsession,:))],2);
    HPC_PSTH_ripple_power(1,nsession,:) = mean([temp{1}(1,:); temp{2}(1,:)]);
    HPC_PSTH_ripple_power(2,nsession,:) = mean([temp{1}(2,:); temp{2}(2,:)]);
    toc
end




% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_MUA_PSTH_all_POST_spindle_power_ripple_power_300ms_interval_raw.mat'),'V1_PSTH_spindle_power','V1_PSTH_ripple_power','V1_PSTH_SO',...
%     'HPC_PSTH_spindle_power','HPC_PSTH_ripple_power','HPC_PSTH_SO','-v7.3')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_MUA_PSTH_all_POST_spindle_power_ripple_power_300ms_interval_raw.mat'),'V1_PSTH_spindle_power','V1_PSTH_ripple_power','V1_PSTH_SO',...
%     'HPC_PSTH_spindle_power','HPC_PSTH_ripple_power','HPC_PSTH_SO')

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_MUA_PSTH_all_POST_spindle_power_ripple_power_300ms_interval.mat'),'V1_PSTH_spindle_power','V1_PSTH_ripple_power','V1_PSTH_SO',...
%     'HPC_PSTH_spindle_power','HPC_PSTH_ripple_power','HPC_PSTH_SO','-v7.3')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_MUA_PSTH_all_POST_spindle_power_ripple_power_300ms_interval.mat'))


save(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_MUA_PSTH_all_POST_ripple_spindle_SO_300ms_interval.mat'),...
    'V1_PSTH_spindle_power','V1_PSTH_ripple_power','V1_PSTH_SO','contra_V1_PSTH_spindle_power','contra_V1_PSTH_SO',...
    'ipsi_contra_diff_V1_PSTH_spindle_power','ipsi_contra_diff_V1_PSTH_SO',...
    'HPC_PSTH_spindle_power','HPC_PSTH_ripple_power','HPC_PSTH_SO','-v7.3')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_MUA_PSTH_all_POST_ripple_spindle_SO_300ms_interval.mat'))

% save(fullfile(analysis_folder,'ripple_MUA_PSTH_all_POST_spindle_power_ripple_power.mat'),'V1_PSTH_spindle_power','V1_PSTH_ripple_power',...
%     'HPC_PSTH_spindle_power','HPC_PSTH_ripple_power','-v7.3')
% load(fullfile(analysis_folder,'ripple_MUA_PSTH_all_POST_spindle_power_ripple_power.mat'),'V1_PSTH_spindle_power','V1_PSTH_ripple_power',...
%     'HPC_PSTH_spindle_power','HPC_PSTH_ripple_power')

% save(fullfile(analysis_folder,'ripple_V1_MUA_PSTH_all_POST_300ms_interval.mat'),'ripple_V1_MUA_PSTH_all','ripple_V1_MSUA_PSTH_all','-v7.3')



%%%%%% Spindle
% temp_PSTH = V1_PSTH_spindle_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_spindle_power.*psthBinSize;
temp_PSTH = V1_PSTH_spindle_power;
temp_HPC_PSTH = HPC_PSTH_spindle_power;
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold spindle power)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power) raw FR';
fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power)';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'low spindle', 'high spindle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'low spindle', 'high spindle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])


subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])





%%%%%% SO phase
% temp_PSTH = V1_PSTH_spindle_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_spindle_power.*psthBinSize;
temp_PSTH = V1_PSTH_SO;
temp_HPC_PSTH = HPC_PSTH_SO;
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold SO phase)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold SO phase) raw FR';
fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold SO phase)';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'peak', 'trough'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'peak', 'trough'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%%%% Ripple power
% temp_PSTH = V1_PSTH_ripple_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_ripple_power.*psthBinSize;
temp_PSTH = V1_PSTH_ripple_power;
temp_HPC_PSTH = HPC_PSTH_ripple_power
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold ripple power)';
fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold ripple power)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold ripple power) raw FR';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'low ripple', 'high ripple'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'low ripple', 'high ripple'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])

%%%%%% CONTRA Spindle
% temp_PSTH = V1_PSTH_spindle_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_spindle_power.*psthBinSize;
temp_PSTH = contra_V1_PSTH_spindle_power;
temp_HPC_PSTH = HPC_PSTH_spindle_power;
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold spindle power)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power) raw FR';
fig.Name = 'contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power)';
% fig.Name = 'ipsi-contra diff ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power)';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'low spindle', 'high spindle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'low spindle', 'high spindle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])


subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])




%%%%%% ipsi-contra diff Spindle
% temp_PSTH = V1_PSTH_spindle_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_spindle_power.*psthBinSize;
temp_PSTH = ipsi_contra_diff_V1_PSTH_spindle_power;
temp_HPC_PSTH = HPC_PSTH_spindle_power;
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold spindle power)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power) raw FR';
% fig.Name = 'contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power)';
fig.Name = 'ipsi-contra diff ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold spindle power)';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'low spindle', 'high spindle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'low spindle', 'high spindle'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])


subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])





%%%%%% SO phase
% temp_PSTH = V1_PSTH_spindle_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_spindle_power.*psthBinSize;
temp_PSTH = contra_V1_PSTH_SO;
temp_HPC_PSTH = HPC_PSTH_SO;
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold SO phase)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold SO phase) raw FR';
fig.Name = 'contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold SO phase)';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'peak', 'trough'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'peak', 'trough'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%%%% SO phase
% temp_PSTH = V1_PSTH_spindle_power.*psthBinSize;
% temp_HPC_PSTH = HPC_PSTH_spindle_power.*psthBinSize;
temp_PSTH = ipsi_contra_diff_V1_PSTH_SO;
temp_HPC_PSTH = HPC_PSTH_SO;
fig = figure;
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.1s inter-rippe threshold SO phase)';
% fig.Name = 'ipsi-contra ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold SO phase) raw FR';
fig.Name = 'ipsi contra diff ripple V1 and HPC MUA PSTH (0.3s inter-rippe threshold SO phase)';
fig.Position = [600 200 820 800];

% sessions_to_process() = []
ipsi_lower = mean(squeeze(temp_PSTH(1,:,:))) - std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_PSTH(1,:,:))) + std(squeeze(temp_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_PSTH(2,:,:))) - std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_PSTH(2,:,:))) + std(squeeze(temp_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_PSTH(2,:,:)));

psthBinSize = 0.01;
windows = [-2 2];
x = windows(1)+psthBinSize/2:psthBinSize:windows(end)-psthBinSize/2;

c_ipsi = [35,139,69]/256; 
c_contra = [106,81,163]/256;
% subplot(3,2,1)


subplot(3,2,1)
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
legend([h1 h2], {'peak', 'trough'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

ipsi_lower = mean(squeeze(temp_HPC_PSTH(1,:,:))) - std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_upper = mean(squeeze(temp_HPC_PSTH(1,:,:))) + std(squeeze(temp_HPC_PSTH(1,:,:)))/sqrt(length(sessions_to_process));
ipsi_mean = mean(squeeze(temp_HPC_PSTH(1,:,:)));

contra_lower = mean(squeeze(temp_HPC_PSTH(2,:,:))) - std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_upper = mean(squeeze(temp_HPC_PSTH(2,:,:))) + std(squeeze(temp_HPC_PSTH(2,:,:)))/sqrt(length(sessions_to_process));
contra_mean = mean(squeeze(temp_HPC_PSTH(2,:,:)));

subplot(3,2,2)
hold on;
% --- Plot Ipsi (Shading first, then Line) ---
fill([x fliplr(x)], [ipsi_lower fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h1 = plot(x, ipsi_mean, 'Color', c_ipsi, 'LineWidth', 2);

% --- Plot Contra (Shading first, then Line) ---
fill([x fliplr(x)], [contra_lower fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
h2 = plot(x, contra_mean, 'Color', c_contra, 'LineWidth', 2);
xlim([-1.5 1.5])
xline([0],'--')
xlabel('Time relative to ripple onset (s)');
ylabel('HPC MUA (z)');
title('HPC');
% Add a legend (using the plot handles h1 and h2)
legend([h1 h2], {'peak', 'trough'}, 'Location', 'best', 'Box', 'off');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-1.5 1.5])

subplot(3,2,3)
data1 = squeeze(mean(temp_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(3,2,4)
data1 = squeeze(mean(temp_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);

rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');

hold off;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,5)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>0 & x<0.2),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>0 & x<0.2),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'post low', 'post high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 1, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(3,2,6)
data1 = squeeze(mean(temp_HPC_PSTH(1,:,x>-0.2 & x<0),3))';
data2 = squeeze(mean(temp_HPC_PSTH(2,:,x>-0.2 & x<0),3))';
[p_val, h_stats] = signrank(data1, data2);
fprintf('Wilcoxon signed-rank test p-value: %.4f\n', p_val);
% 
rng(1); % Seed for consistent look
jitter_strength = 0.1; 
x1 = 1 + (rand(size(data1)) - 0.5) * jitter_strength;
x2 = 1.5 + (rand(size(data2)) - 0.5) * jitter_strength;

% Define X-axis positions for the two groups
x_pos = [x1, x2]';

% Plot the connecting lines first (so they are in the background)
% We concatenate data into a 2-column matrix and plot its transpose
plot(x_pos, [data1, data2]', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);

% Plot individual data points with transparency (Alpha)
hold on
s1 = scatter(x1, data1, 100, c_ipsi, 'filled', 'MarkerFaceAlpha', 0.5); % Green
s2 = scatter(x2, data2, 100, c_contra, 'filled', 'MarkerFaceAlpha', 0.5); % Pink

% 4. Formatting the Plot
set(gca, 'XTick', [1 1.5], 'XTickLabel', {'pre low', 'pre high'});
ylabel('Mean PSTH Activity');
xlim([0.5 2]);
grid off;
box off;
p_string = sprintf('p = %.4e', p_val);
text(1.5, 0.2, p_string, 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])
