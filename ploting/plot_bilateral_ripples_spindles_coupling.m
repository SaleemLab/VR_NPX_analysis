function plot_ipsi_contra_ripples_spindles_coupling(slow_waves_all,ripples_all,spindles_all)

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
% load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
% load(fullfile(analysis_folder,'ripples_all_POST.mat'))
% load(fullfile(analysis_folder,'spindles_all_POST.mat'))
% load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));



% Find reference channel/shank
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spindle_phase=[];
% for probe_no = 1:2
%     spindle_phase{probe_no}=[];
%     for nsession = 1:length(sessions_to_process)
%         spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,probe_no),:)];
%     end
% end
% spindle_phase = [spindle_phase{1}(ripples_all(1).SWS_index==1) spindle_phase{2}(ripples_all(2).SWS_index==1)];
% ripple_info.spindle_phase = spindle_phase';


%%%%%%%%%%%%%%%%%% Probability of ripples relative to spindles
load(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_whole_baseline_combined.mat'),'probability');
probability_psth_whole_baseline = probability;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_baseline.mat'),'probability');
% probability_psth_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_whole_combined.mat'),'probability');
probability_psth_whole = probability;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','spindles_ripples_probability_combined.mat'),'probability');
% probability_psth = probability;


time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


%%%%%%%%%%%%%%%%%%%%%
ipsi_probability = [probability_psth_whole(1).ripple; probability_psth_whole(2).ripple];
% contra_probability = [probability_psth_whole(1).R_ripple; probability_psth_whole(2).L_ripple];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).ripple; probability_psth_whole_baseline(2).ripple];
% contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripple; probability_psth_whole_baseline(2).L_ripple];

%%%%% calculate shuffled baseline
%%% Ipsi
binnedArray = ipsi_probability_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;



% clear probability_merged
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

probability_merged.x = x;
probability_merged.ripples_spindles = [];
% probability_merged.contra_ripples_spindles = [];
% probability_merged.ipsi_contra_diff_ripples_spindles = [];

event_idx = [];
peak_power = [spindles_all(1).peak_zscore(spindles_all(1).SWS_index); spindles_all(2).peak_zscore(spindles_all(2).SWS_index)];

event_idx{1} = {find(peak_power<prctile(peak_power,50)),find(peak_power>prctile(peak_power,50))};

event_idx{2} = {(1:size(ipsi_probability,1))'};

title_names = {'Ipsi-contra spindles ripples by spindle powers','Left-Right combined ipsi contra ripple distribution around spindle onset'};
group_name = [];
group_name{1} = {'low power','high power'};


for i = 1:length(event_idx)
    for ngroup = 1:length(event_idx{i})
        index =event_idx{i}{ngroup};

        binnedArray1 = ipsi_probability(index,:);
        % binnedArray2 = contra_probability(index,:);
        % binnedArray3 = binnedArray1-binnedArray2;

        temp1=[];
        temp2=[];
        temp3=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            % temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            % temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        probability_merged.ripples_spindles{i}{ngroup} = temp1;
        % probability_merged.contra_ripples_spindles{i}{ngroup}= temp2;
        % probability_merged.ipsi_contra_diff_ripples_spindles{i}{ngroup} = temp3;
    end
end

probability_merged.ripples_baseline_spindles = ipsi_baseline_bootstrap;
% probability_merged.contra_ripples_baseline_spindles = contra_baseline_bootstrap;
% probability_merged.ipsi_contra_diff_ripples_baseline_spindles = ipsi_contra_diff_baseline_bootstrap;
probability_merged.ripples_spindles_groups = [title_names];
probability_merged.ripples_spindles_index = [event_idx];



% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','spindles_ripples_probability_merged.mat'),'probability_merged')
save(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH','SO_ripples_spindles_probability_merged.mat'),'probability_merged')



time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

%%%%%%%%% Plot spindle onset
for ngroup = 1:length(event_idx)-1
    fig = figure('Color','w');
    fig.Position = [350 59 1650 465];
    fig.Name =title_names{ngroup};

    % if ngroup ==1
    %     colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    % elseif ngroup ==5
    %     colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    % else
    if ngroup ==1 
        colour_lines = [161,217,155;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;74,20,134]/256;% 5 purple for
    end


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ripples_spindles{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ripples_baseline_spindles;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([0 0.15])
    title('ipsi ripples')
    xlabel('Time relative to spindle onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    % 
    % if ngroup ==1 
    %     % colour_lines = [161,217,155;0,90,50]/256;% 5 green for
    %     colour_lines = [188,189,220;74,20,134]/256;% 5 purple for
    % end

    nexttile
  
    nexttile
   
end


%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All spindles rpples
event_averaging_scale = 10;
duration = [spindles_all(1).SWS_offset-spindles_all(1).SWS_onset; spindles_all(2).SWS_offset-spindles_all(2).SWS_onset];

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650 465];
fig.Name = 'Left-Right combined ipsi contra ripple distribution around spindle onset'

colour_lines = [0,90,50;74,20,134]/256; % Green Purple


nexttile
% duration = event_times(:,2) - event_times(:,1);
[~,sorted_index] = sort(duration);
imagesc(event_averaging_scale*movmean(ipsi_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)

hold on
% Convert duration (in seconds) to number of bins
duration_in_bins = duration(sorted_index) / 0.02;

% Add to center point (DOWN-UP at bin 101)
duration_bin_position = 51 + duration_in_bins;
% Plot yellow dashed line
plot(flip(duration_bin_position), flip(1:numel(duration)), 'r--', 'LineWidth', 1)


clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to spindle onset (s)')
ylabel('Event sorted by spindle duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('ipsi ripples')

nexttile
% imagesc(event_averaging_scale*movmean(contra_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([1.5 25.5 50.5 75.5 100.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([-1 -0.5 0 0.5 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Time relative to spindle onset (s)')
% ylabel('Event sorted by DOWN duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('contra ripples')

nexttile
clear ERROR_SHADE

binnedArray = probability_merged.ripples_spindles{end}{1};
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)


% baseline
binnedArray = probability_merged.ripples_baseline_spindles;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.14])
% title('ipsi ripples')
xlabel('Time relative to spindle onset (s)')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'real','shuffled'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH'),[],'ContentType','vector')


% 
% %%%%%%%%%%%%%%%%%%%%%
% ipsi_probability = [probability_psth(1).L_ripple; probability_psth(2).R_ripple];
% contra_probability = [probability_psth(1).R_ripple; probability_psth(2).L_ripple];
% 
% ipsi_probability_baseline = [probability_psth_baseline(1).L_ripple; probability_psth_baseline(2).R_ripple];
% contra_probability_baseline = [probability_psth_baseline(1).R_ripple; probability_psth_baseline(2).L_ripple];
% 
% %%%%% calculate shuffled baseline
% %%% Ipsi
% binnedArray = ipsi_probability_baseline;
% temp=[];
% parfor iBoot = 1:1000
%     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
%     event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
%     temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
%     % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
% end
% 
% ipsi_baseline_bootstrap = temp;
% 
% %%% Contra
% binnedArray = contra_probability_baseline;
% temp=[];
% parfor iBoot = 1:1000
%     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
%     event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
%     temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
%     % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
% end
% 
% contra_baseline_bootstrap = temp;
% 
% %%% Ipsi-Contra baseline
% binnedArray = ipsi_probability_baseline-contra_probability_baseline;
% temp=[];
% parfor iBoot = 1:1000
%     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
%     event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
%     temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
%     % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
% end
% 
% ipsi_contra_diff_baseline_bootstrap = temp;
% 
% 
% 
% 
% event_idx = [];
% peak_power = [spindles_all(1).peak_zscore(spindles_all(1).SWS_index); spindles_all(2).peak_zscore(spindles_all(2).SWS_index)];
% 
% event_idx{1} = {find(peak_power<prctile(peak_power,50)),find(peak_power>prctile(peak_power,50))};
% 
% event_idx{2} = {(1:size(ipsi_probability,1))'};
% 
% title_names = {'Ipsi-contra spindles ripple onset by spindle powers','Left-Right combined ipsi contra ripple onset around spindle onset'};
% group_name = [];
% group_name{1} = {'low power','high power'};
% 
% 
% % clear probability_merged
% time_wondows = [-1 1];
% time_bin = 0.02;
% x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
% 
% probability_merged.x = x;
% probability_merged.ipsi_ripple_onset_spindles = [];
% probability_merged.contra_ripple_onset_spindles = [];
% probability_merged.ipsi_contra_diff_ripple_onset_spindles = [];
% 
% for i = 1:length(event_idx)
%     for ngroup = 1:length(event_idx{i})
%         index =event_idx{i}{ngroup};
% 
%         binnedArray1 = ipsi_probability(index,:);
%         binnedArray2 = contra_probability(index,:);
%         binnedArray3 = binnedArray1-binnedArray2;
% 
%         temp1=[];
%         temp2=[];
%         temp3=[];
%         parfor iBoot = 1:1000
%             s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
%             event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
%             temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
%             temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
%             temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
%             % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
%         end
% 
%         probability_merged.ipsi_ripple_onset_spindles{i}{ngroup} = temp1;
%         probability_merged.contra_ripple_onset_spindles{i}{ngroup}= temp2;
%         probability_merged.ipsi_contra_diff_ripple_onset_spindles{i}{ngroup} = temp3;
%     end
% end
% 
% probability_merged.ipsi_ripple_onset_baseline_spindles = ipsi_baseline_bootstrap;
% probability_merged.contra_ripple_onset_baseline_spindles = contra_baseline_bootstrap;
% probability_merged.ipsi_contra_diff_ripple_onset_baseline_spindles = ipsi_contra_diff_baseline_bootstrap;
% probability_merged.ripple_onset_spindles_groups = [title_names];
% probability_merged.ripple_onset_spindles_index = [event_idx];
% 
% 
% 
% 
% time_wondows = [-1 1];
% time_bin = 0.02;
% x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
% 
% %%%%%%%%% Plot DU transition
% for ngroup = 1:length(event_idx)-1
%     fig = figure('Color','w');
%     fig.Position = [350 59 1650 465];
%     fig.Name =title_names{ngroup};
% 
%     % if ngroup ==1
%     %     colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
%     % elseif ngroup ==5
%     %     colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
%     % else
%     if ngroup ==1 
%         colour_lines = [161,217,155;0,90,50]/256;% 5 green for
%         % colour_lines = [188,189,220;74,20,134]/256;% 5 purple for
%     end
% 
% 
%     nexttile
%     clear ERROR_SHADE
%     for i = 1:length(event_idx{ngroup})
% 
%         binnedArray = probability_merged.ipsi_ripple_onset_spindles{ngroup}{i};
% 
%         y = mean(binnedArray,'omitnan');
%         LCI = prctile(binnedArray,2.5);
%         UCI = prctile(binnedArray,97.5);
% 
%         PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
%         ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
%     end
% 
%     % baseline
%     binnedArray = probability_merged.ipsi_ripples_baseline_spindles;
%     y = mean(binnedArray,'omitnan');
%     %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
%     LCI = prctile(binnedArray,2.5);
%     UCI = prctile(binnedArray,97.5);
% 
%     PLOT = plot(x,y,'k');hold on;
%     ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
% %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
% 
%     % xline(0,'r')
%     ylim([0 0.11])
%     title('ipsi ripples')
%     xlabel('Time relative to spindle onset (s)')
%     ylabel('Probability')
%     set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
%     if ngroup ==1 
%         % colour_lines = [161,217,155;0,90,50]/256;% 5 green for
%         colour_lines = [188,189,220;74,20,134]/256;% 5 purple for
%     end
% 
%     nexttile
%     clear ERROR_SHADE
%     for i = 1:length(event_idx{ngroup})
% 
%         binnedArray = probability_merged.contra_ripple_onset_spindles{ngroup}{i};
% 
%         y = mean(binnedArray,'omitnan');
%         LCI = prctile(binnedArray,2.5);
%         UCI = prctile(binnedArray,97.5);
% 
%         PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
%         ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
%     end
% 
%     % baseline
%     binnedArray = probability_merged.contra_ripples_baseline_spindles;
%     y = mean(binnedArray,'omitnan');
%     %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
%     LCI = prctile(binnedArray,2.5);
%     UCI = prctile(binnedArray,97.5);
%     PLOT = plot(x,y,'k');hold on;
%     ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
% %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
%     ylim([0 0.11])
% 
%     % xline(0,'r')
%     title('contra ripples')
%     xlabel('Time relative to spindle onset (s)')
%     ylabel('Probability')
%     set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
%     if ngroup ==1 
%         colour_lines = [161,217,155;0,90,50]/256;% 5 green for
%         % colour_lines = [188,189,220;74,20,134]/256;% 5 purple for
%     end
% 
%     nexttile
%     clear ERROR_SHADE
%     for i = 1:length(event_idx{ngroup})
% 
%         binnedArray = probability_merged.ipsi_contra_diff_ripple_onset_spindles{ngroup}{i};
% 
%         y = mean(binnedArray,'omitnan');
%         LCI = prctile(binnedArray,2.5);
%         UCI = prctile(binnedArray,97.5);
% 
%         PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
%         ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
%     end
% 
%     % baseline
%     binnedArray = probability_merged.ipsi_contra_diff_ripple_onset_baseline_spindles;
%     y = mean(binnedArray,'omitnan');
%     %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
%     LCI = prctile(binnedArray,2.5);
%     UCI = prctile(binnedArray,97.5);
% 
%     PLOT = plot(x,y,'k');hold on;
%     ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
%     ylim([-0.04 0.04])
% 
% 
%     % xline(0,'r')
%     title('Ipsi-contra diff')
%     xlabel('Time relative to spindle onset (s)')
%     ylabel('Probability')
%     set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% end
% 
% %%%%%%%%%%%%
% %%%%%%%%%%%%
% %%%%%%%%%%%% All spindles rpples
% event_averaging_scale = 15;
% duration = [spindles_all(1).SWS_offset-spindles_all(1).SWS_onset; spindles_all(2).SWS_offset-spindles_all(2).SWS_onset];
% 
% % for ngroup = 1:length(event_idx)
% ngroup = 2;
% fig = figure('Color','w');
% fig.Position = [350 59 1650 465];
% fig.Name = 'Left-Right combined ipsi contra ripple onset distribution around spindle onset'
% 
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple
% 
% 
% nexttile
% % duration = event_times(:,2) - event_times(:,1);
% [~,sorted_index] = sort(duration);
% imagesc(event_averaging_scale*movmean(ipsi_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([1.5 25.5 50.5 75.5 100.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([-1 -0.5 0 0.5 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Time relative to spindle onset (s)')
% ylabel('Event sorted by spindle duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('ipsi ripples')
% 
% nexttile
% imagesc(event_averaging_scale*movmean(contra_probability(sorted_index,:),event_averaging_scale,1,'omitnan'))
% % imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
% xticks([1.5 25.5 50.5 75.5 100.5])
% % xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
% xticklabels([-1 -0.5 0 0.5 1])
% xline(50.5,'r',LineWidth=1)
% clim([0 1])
% colorbar
% colormap(flipud(gray))
% xlabel('Time relative to UP-DOWN transition (s)')
% ylabel('Event sorted by DOWN duration')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% title('contra ripples')
% 
% nexttile
% clear ERROR_SHADE
% 
% binnedArray = probability_merged.ipsi_ripple_onset_spindles{end}{1};
% y = mean(binnedArray,'omitnan');
% LCI = prctile(binnedArray,2.5);
% UCI = prctile(binnedArray,97.5);
% 
% PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
% ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)
% 
% binnedArray = probability_merged.contra_ripple_onset_spindles{end}{1};
% y = mean(binnedArray,'omitnan');
% LCI = prctile(binnedArray,2.5);
% UCI = prctile(binnedArray,97.5);
% 
% PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
% ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)
% 
% % baseline
% binnedArray = probability_merged.ipsi_ripple_onset_baseline_spindles;
% y = mean(binnedArray,'omitnan');
% %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
% LCI = prctile(binnedArray,2.5);
% UCI = prctile(binnedArray,97.5);
% 
% PLOT = plot(x,y,'k');hold on;
% ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
% %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
% 
% % xline(0,'r')
% ylim([0 0.04])
% % title('ipsi ripples')
% xlabel('Time relative to spindle onset (s)')
% ylabel('Probability')
% legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')
% 
% 
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','spindles_ripples_probability_merged.mat'),'probability_merged')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% DOWN UP with and without spindles
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole.mat'),'probability');
spindles_probability = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline_combined.mat'),'probability');
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'),'probability');
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'),'probability');
probability_psth = probability;






% ipsi_ripples = [probability_psth(1).L_ripples_UP; probability_psth(2).R_ripples_UP];
% contra_ripples = [probability_psth(1).R_ripples_UP; probability_psth(2).L_ripples_UP];
% ripples_combined = ipsi_ripples+contra_ripples;

ripples_combined = [probability_psth_whole(1).ripples_UP;probability_psth_whole(2).ripples_UP];

ipsi_spindles = [spindles_probability(1).L_spindles_UP; spindles_probability(2).R_spindles_UP];
contra_spindles = [spindles_probability(1).R_spindles_UP; spindles_probability(2).L_spindles_UP];
combined_spindles = ipsi_spindles+contra_spindles;

[nEvents, nBins] = size(ipsi_ripples);

ipsi_eventIndices=[];
contra_eventIndices=[];
eventIndices=[];

for ngroup = 1:3
    firstRipple = NaN(nEvents,1);
    firstSpindle = NaN(nEvents,1);
    spindleBeforeRipple = NaN(nEvents,1);

    for i = 1:nEvents
        r = find(ripples_combined(i, :) >0 , 1, 'first');  % first ripple
        
        if ngroup == 1
            s = find(ipsi_spindles(i, :) >0 , 1, 'first'); % first spindle
        elseif ngroup == 2
            s = find(contra_spindles(i, :) >0 , 1, 'first'); % first spindle
        else
            s = find(combined_spindles(i, :) >0 , 1, 'first'); % first spindle
        end

        if ~isempty(s)
            firstSpindle(i) = s;
        end

        if ~isempty(r)
            firstRipple(i) = r;
            if ~isempty(s) & s <= r
                spindleBeforeRipple(i) = s;
            end
        end
    end

    if ngroup == 1
        % Get event indices where spindle happened before ripple

        ipsi_eventIndices{1} = find(~isnan(firstSpindle));
        ipsi_eventIndices{2} = find(~isnan(spindleBeforeRipple));
    elseif ngroup == 2
        contra_eventIndices{1} = find(~isnan(firstSpindle));
        contra_eventIndices{2} = find(~isnan(spindleBeforeRipple));
    else
        eventIndices{1}= find(~isnan(firstSpindle));
        eventIndices{2}= find(~isnan(spindleBeforeRipple));
    end
end

event_idx = {};  % Will hold 8 groups: 4 experimental, 4 control

all_event_ids = (1:nEvents)';

% Define experimental sets
groups = {
    ipsi_eventIndices{1};     % 1. Ipsi Spindle present
    ipsi_eventIndices{2};     % 2. Ipsi Spindle before ripple
    contra_eventIndices{1};   % 3. Contra Spindle present
    contra_eventIndices{2};   % 4. Contra Spindle before ripple
    };

group_name{1} = {'with spindles','without spindles'};
group_name{2} = {'with spindles','without spindles'};
group_name{3} = {'with spindles','without spindles'};
group_name{4} = {'with spindles','without spindles'};

title_names{1} = 'Ipsi-contra DOWN_UP ripples with ipsi spindles';
title_names{2} = 'Ipsi-contra DOWN_UP ripples with ipsi spindles before ripples';
title_names{3} = 'Ipsi-contra DOWN_UP ripples with contra spindles before ripples';
title_names{4} = 'Ipsi-contra DOWN_UP ripples with contra spindles before ripples';

for g = 1:4
    exp_ids = groups{g};
    ctrl_pool = setdiff(all_event_ids, exp_ids);

    % Save both with and without spindles
    event_idx{g}{1} = groups{g};   
    event_idx{g}{2} = setdiff(all_event_ids, groups{g});  
end




% ipsi_probability = [probability_psth(1).L_ripples_UP; probability_psth(2).R_ripples_UP];

ipsi_probability = [probability_psth_whole(1).ripples_UP; probability_psth_whole(2).ripples_UP];
% contra_probability = [probability_psth_whole(1).R_ripples_UP; probability_psth_whole(2).L_ripples_UP];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).ripples_UP; probability_psth_whole_baseline(2).ripples_UP];
% contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_UP; probability_psth_whole_baseline(2).L_ripples_UP];

%%%%% calculate shuffled baseline
%%% Ipsi
binnedArray = ipsi_probability_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),length(event_idx{1}{1}));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;
% 
% %%% Contra
% binnedArray = contra_probability_baseline;
% temp=[];
% parfor iBoot = 1:1000
%     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
%     event_id = datasample(s,1:size(binnedArray,1),length(event_idx{1}{1}));
%     temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
%     % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
% end
% 
% contra_baseline_bootstrap = temp;
% 
% %%% Ipsi-Contra baseline
% binnedArray = ipsi_probability_baseline-contra_probability_baseline;
% temp=[];
% parfor iBoot = 1:1000
%     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
%     event_id = datasample(s,1:size(binnedArray,1),length(event_idx{1}{1}));
%     temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
%     % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
% end
% 
% ipsi_contra_diff_baseline_bootstrap = temp;


% clear probability_merged
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

probability_merged.x = x;
probability_merged.ripples_spindles_UP = [];
% probability_merged.contra_ripples_UP = [];
% probability_merged.ipsi_contra_diff_ripples_UP = [];

for ngroup = 1:length(event_idx)

    group_index = [];
 
    group_index = event_idx{ngroup};


    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_probability(index,:);
        % binnedArray2 = contra_probability(index,:);
        % binnedArray3 = binnedArray1-binnedArray2;

        temp1=[];
        temp2=[];
        temp3=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),length(group_index{1}));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            % temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            % temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        probability_merged.ripples_spindles_UP{ngroup}{i} = temp1;
        % probability_merged.contra_ripples_UP{ngroup}{i} = temp2;
        % probability_merged.ipsi_contra_diff_ripples_UP{ngroup}{i} = temp3;
    end
end


probability_merged.ripples_spindles_baseline_UP = ipsi_baseline_bootstrap;
% probability_merged.contra_ripples_baseline_UP = contra_baseline_bootstrap;
% probability_merged.ipsi_contra_diff_ripples_baseline_UP = ipsi_contra_diff_baseline_bootstrap;
probability_merged.ripples_UP_groups = [title_names ];
probability_merged.ripples_UP_index = [event_idx];


time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 465];
    fig.Name =title_names{ngroup};

    % if ngroup ==1
    %     colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    % elseif ngroup ==5
    %     colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    % else
        colour_lines = [161,217,155;0,90,50]/256;% 5 green for
        % colour_lines = [188,189,220;74,20,134]/256;% 5 purple for

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = probability_merged.ripples_spindles_UP{ngroup}{i};

        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = probability_merged.ripples_spindles_baseline_UP;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([0 0.1])
    title('ipsi ripples')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % colour_lines = [161,217,155;0,90,50]/256;% 5 green for
    colour_lines = [188,189,220;74,20,134]/256;% 5 purple for

    nexttile
%     clear ERROR_SHADE
%     for i = 1:length(event_idx{ngroup})
% 
%         binnedArray = probability_merged.contra_ripples_UP{ngroup}{i};
% 
%         y = mean(binnedArray,'omitnan');
%         LCI = prctile(binnedArray,2.5);
%         UCI = prctile(binnedArray,97.5);
% 
%         PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
%         ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
%         xline(0,'r',LineWidth=1)
%     end
% 
%     % baseline
%     binnedArray = probability_merged.contra_ripples_baseline_UP;
%     y = mean(binnedArray,'omitnan');
%     %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
%     LCI = prctile(binnedArray,2.5);
%     UCI = prctile(binnedArray,97.5);
%     PLOT = plot(x,y,'k');hold on;
%     ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
% %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
%     ylim([0 0.08])
% 
%     % xline(0,'r')
%     title('contra ripples')
%     xlabel('Time relative to DOWN-UP transition (s)')
%     ylabel('Probability')
%     set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
%         colour_lines = [161,217,155;0,90,50]/256;% 5 green for
%         % colour_lines = [188,189,220;74,20,134]/256;% 5 purple for

    nexttile
    % clear ERROR_SHADE
    % for i = 1:length(event_idx{ngroup})
    % 
    %     binnedArray = probability_merged.ipsi_contra_diff_ripples_UP{ngroup}{i};
    % 
    %     y = mean(binnedArray,'omitnan');
    %     LCI = prctile(binnedArray,2.5);
    %     UCI = prctile(binnedArray,97.5);
    % 
    %     PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
    %     ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
    %     xline(0,'r',LineWidth=1)
    % end
    % 
    % % baseline
    % binnedArray = probability_merged.ipsi_contra_diff_ripples_baseline_UP;
    % y = mean(binnedArray,'omitnan');
    % %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(binnedArray,2.5);
    % UCI = prctile(binnedArray,97.5);
    % 
    % PLOT = plot(x,y,'k');hold on;
    % ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    % ylim([-0.04 0.04])
    % 
    % 
    % % xline(0,'r')
    % title('Ipsi-contra diff')
    % xlabel('Time relative to DOWN-UP transition (s)')
    % ylabel('Probability')
    % set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')
% 
% 
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','spindles_ripples_probability_merged.mat'),'probability_merged')


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH','SO_ripples_spindles_probability_merged.mat'),'probability_merged')
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH'),[],'ContentType','vector')


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole_baseline.mat'),'probability');
% probability_psth_whole_baseline = probability;
% 
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_baseline.mat'),'probability');
% probability_psth_baseline = probability;
% 
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'),'probability');
% probability_psth_whole = probability;
% 
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability.mat'),'probability');
% probability_psth = probability;