% save_all_figures('P:\corticohippocampal_replay\V1-HPC bilateral interaction\UP_DOWN_ripples_lme\mixed effect regression (full windows control)',[],'SVG_option',1)
% save_all_figures('P:\corticohippocampal_replay\V1-HPC bilateral interaction\UP_DOWN_ripples_lme\mixed effect regression (full windows)',[],'SVG_option',1)
% save_all_figures('P:\corticohippocampal_replay\V1-HPC bilateral interaction\UP_DOWN_ripples_lme\mixed effect regression (150ms lag windows control)',[],'SVG_option',1)
% save_all_figures('P:\corticohippocampal_replay\V1-HPC bilateral interaction\UP_DOWN_ripples_lme\mixed effect regression (150ms lag windows)',[],'SVG_option',1)

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline_combined.mat'));
probability_psth_whole_baseline = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;


% 
% [ripple_probability_UP,event_index,normalized_duration,temp] = calculate_relative_event_probability(merged_event_info.UP_ints,merged_event_info.ripples_ints(ripple_info.event_id,:),30,0,'Range',[0 1]);
% hold on;
% plot(0.025:0.05:0.975,mean(temp,1))
% 
%%%%%% Ripple PSTH
clear ripple_probability
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

% ripple_probability_UP = [probability_psth_whole(1).ripples_UP; probability_psth_whole(2).ripples_UP];
% ripple_probability_UP_baseline = [probability_psth_whole_baseline(1).ripples_UP; probability_psth_whole_baseline(2).ripples_UP];
% ripple_probability_DOWN = [probability_psth_whole(1).ripples_DOWN; probability_psth_whole(2).ripples_DOWN];
% ripple_probability_DOWN_baseline = [probability_psth_whole_baseline(1).ripples_DOWN; probability_psth_whole_baseline(2).ripples_DOWN];

ripple_probability.UP = [probability_psth_whole(1).ripples_UP; probability_psth_whole(2).ripples_UP];
ripple_probability.UP_baseline = [probability_psth_whole_baseline(1).ripples_UP; probability_psth_whole_baseline(2).ripples_UP];
ripple_probability.DOWN = [probability_psth_whole(1).ripples_DOWN; probability_psth_whole(2).ripples_DOWN];
ripple_probability.DOWN_baseline = [probability_psth_whole_baseline(1).ripples_DOWN; probability_psth_whole_baseline(2).ripples_DOWN];

ripple_probability_UP = ripple_probability.UP;
ripple_probability_UP_baseline = ripple_probability.UP_baseline;
ripple_probability_DOWN = ripple_probability.DOWN;
ripple_probability_DOWN_baseline = ripple_probability.DOWN_baseline;


% bootstrap distribution
tempUP = [];
tempDOWN = [];

parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_probability_UP,1),size(ripple_probability_UP,1));
    tempUP(iBoot,:) = sum(ripple_probability_UP(event_id,:),'omitnan')./sum(~isnan(ripple_probability_UP(event_id,:)));

    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_probability_DOWN,1),size(ripple_probability_DOWN,1));
    tempDOWN(iBoot,:) = sum(ripple_probability_DOWN(event_id,:),'omitnan')./sum(~isnan(ripple_probability_DOWN(event_id,:)));
end
ripple_probability.UP_boot = tempUP;
ripple_probability.DOWN_boot = tempDOWN;

% baseline bootstrap distribution
tempUP = [];
tempDOWN = [];

parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_probability_UP_baseline,1),size(ripple_probability_UP_baseline,1));
    tempUP(iBoot,:) = sum(ripple_probability_UP_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_probability_UP_baseline(event_id,:)));

    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_probability_DOWN_baseline,1),size(ripple_probability_DOWN_baseline,1));
    tempDOWN(iBoot,:) = sum(ripple_probability_DOWN_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_probability_DOWN_baseline(event_id,:)));
end
ripple_probability.UP_baseline_boot = tempUP;
ripple_probability.DOWN_baseline_boot = tempDOWN;

%%%%%% Normalised Probability
% num_bins = 25+10;
% Range=[-0.2 1.2];
% bin_edges = linspace(Range(1), Range(2), num_bins+1);
% bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;
num_bins = 50+20;
Range=[-0.2 1.2];
bin_edges = linspace(Range(1), Range(2), num_bins+1);
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;



[~,event_index,normalized_duration,ripple_count_UP] = calculate_relative_event_probability(merged_event_info.UP_ints,merged_event_info.ripples_ints(ripple_info.event_id,:),num_bins,0,'Range',Range);
[~,event_index,normalized_duration,ripple_count_DOWN] = calculate_relative_event_probability(merged_event_info.DOWN_ints,merged_event_info.ripples_ints(ripple_info.event_id,:),num_bins,0,'Range',Range);

[~,~,~,ripple_count_UP_baseline] = calculate_relative_event_probability(merged_event_info.UP_ints+3,merged_event_info.ripples_ints(ripple_info.event_id,:),num_bins,0,'Range',Range);
[~,~,~,ripple_count_DOWN_baseline] = calculate_relative_event_probability(merged_event_info.DOWN_ints+3,merged_event_info.ripples_ints(ripple_info.event_id,:),num_bins,0,'Range',Range);

ripple_probability.norm_UP = ripple_count_UP;
ripple_probability.norm_DOWN = ripple_count_DOWN;
ripple_probability.norm_UP_baseline = ripple_count_UP_baseline;
ripple_probability.norm_DOWN_baseline = ripple_count_DOWN_baseline;

% bootstrap distribution
tempUP = [];
tempDOWN = [];

parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_count_UP,1),size(ripple_count_UP,1));
    tempUP(iBoot,:) = sum(ripple_count_UP(event_id,:),'omitnan')./sum(~isnan(ripple_count_UP(event_id,:)));

    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_count_DOWN,1),size(ripple_count_DOWN,1));
    tempDOWN(iBoot,:) = sum(ripple_count_DOWN(event_id,:),'omitnan')./sum(~isnan(ripple_count_DOWN(event_id,:)));
end

ripple_probability.norm_UP_boot = tempUP;
ripple_probability.norm_DOWN_boot = tempDOWN;

% Baseline bootstrap distribution
tempUP = [];
tempDOWN = [];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_count_UP_baseline,1),size(ripple_count_UP_baseline,1));
    tempUP(iBoot,:) = sum(ripple_count_UP_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_count_UP_baseline(event_id,:)));

    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(ripple_count_DOWN_baseline,1),size(ripple_count_DOWN_baseline,1));
    tempDOWN(iBoot,:) = sum(ripple_count_DOWN_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_count_DOWN_baseline(event_id,:)));
end

ripple_probability.norm_UP_baseline_boot = tempUP;
ripple_probability.norm_DOWN_baseline_boot = tempDOWN;
ripple_probability.norm_baseline_boot = (tempUP+tempDOWN)./2;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripple_probability.mat'),'ripple_probability');

% 
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% probability_normalised_whole = probability;

figure
subplot(2,2,1)
plot(bin_centres,mean(ripple_count_UP,1))
hold on;
plot(bin_centres,mean([ripple_count_UP_baseline; ripple_count_DOWN_baseline],1))
xlim([-0.2 1.2])
xticks(-0.2:0.2:1.2)
ylim([0 0.11])
xline([0 1])

subplot(2,2,2)
plot(bin_centres,mean(ripple_count_DOWN,1))

hold on
plot(bin_centres,mean([ripple_count_UP_baseline; ripple_count_DOWN_baseline],1))
xlim([-0.2 1.2])
xticks(-0.2:0.2:1.2)
ylim([0 0.11])
xline([0 1])


% figure
% subplot(2,2,1)
% plot(bin_centres,mean(ripple_probability.norm_UP,1))
% hold on;
% plot(bin_centres,mean([ripple_probability.norm_UP_baseline; ripple_probability.norm_DOWN_baseline],1))
% xlim([-0.2 1.2])
% xticks(-0.2:0.2:1.2)
% ylim([0 0.11])
% xline([0 1])
% 
% subplot(2,2,2)
% plot(bin_centres,mean(ripple_probability.norm_DOWN,1))
% hold on
% plot(bin_centres,mean([ripple_probability.norm_UP_baseline; ripple_probability.norm_DOWN_baseline],1))
% xlim([-0.2 1.2])
% xticks(-0.2:0.2:1.2)
% ylim([0 0.11])
% xline([0 1])


%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All DOWN UP rpples
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
event_averaging_scale = 10;

fig = figure('Color','w');
fig.Position = [350 59 1650 465];
fig.Name = 'Left-Right combined ripple distribution around DOWN-UP transition (10 movemean)'
% fig.Name = 'Left-Right combined ripple distribution around DOWN-UP transition (15 movemean without pre nan)'
% fig.Name = 'Left-Right combined ripple distribution around DOWN-UP transition (by previous DOWN) (15 movemean)'
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
nexttile
%%%% By UP duration
event_times = merged_event_info.UP_ints;
duration = event_times(:,2) - event_times(:,1);
[~,sorted_index] = sort(duration);
temp_matrix = event_averaging_scale*movmean(ripple_probability.UP(sorted_index,:),event_averaging_scale,1,'omitnan');
% temp_matrix = ripple_probability.UP(sorted_index,:);
% Convert duration (in seconds) to number of bins
duration_in_bins = duration(sorted_index) / 0.02;
% Add to center point (DOWN-UP at bin 101)
duration_bin_position = 51 + duration_in_bins;


% for nevent = 1:length(temp_matrix)
%     temp_matrix(nevent,1:round(duration_bin_position(nevent))) = nan;
% end

h= imagesc(temp_matrix)
set(h, 'AlphaData', ~isnan(temp_matrix));
set(gca, 'Color', 'w');
% imagesc(ripple_probability(sorted_index,:))
hold on

% Plot yellow dashed line
plot(flip(duration_bin_position), flip(1:numel(duration)), 'r--', 'LineWidth', 1)


% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
xticks([0.5 13 25.5 38 50.5 62.5 75 87.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xline(50.5,'r',LineWidth=1)
xlim([35 100])

clim([0 1])
colorbar
colormap(flipud(gray))
% colormap(flipud(hot))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('ripples')

nexttile
clear ERROR_SHADE

binnedArray = ripple_probability.UP_boot;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
% 
% binnedArray = probability_merged.contra_ripples_UP{end}{1};
% y = mean(binnedArray,'omitnan');
% LCI = prctile(binnedArray,2.5);
% UCI = prctile(binnedArray,97.5);
% 
% PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
% ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

% baseline
binnedArray = ripple_probability.UP_baseline_boot;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.07])
% title('ipsi ripples')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'real','baseline'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xlim([-0.3 1])

nexttile

% ripple_probability_normalised_baseline = [probability_normalised_whole(1).ripples_UP_shuffled; probability_normalised_whole(2).ripples_UP_shuffled];
clear ERROR_SHADE

num_bins = 50+20;
Range=[-0.2 1.2];
bin_edges = linspace(Range(1), Range(2), num_bins+1);
x = bin_edges(1:end-1) + diff(bin_edges)/2;


binnedArray = ripple_probability.norm_UP_boot;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = ripple_probability.norm_baseline_boot;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.085])
% xlim([-0.025 1.025])
xlim([0 1])

xticks([0 0.25 0.5 0.75 1])
% title('ipsi ripples')
xlabel('Normalised UP Duration')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'real','baseline'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% xticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')





%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% All UP DOWN rpples


time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

event_averaging_scale = 10;

% for ngroup = 1:length(event_idx)
ngroup = 2;
fig = figure('Color','w');
fig.Position = [350 59 1650 465];
fig.Name = 'Left-Right combined ripple distribution around UP-DOWN transition';
fig.Name = 'Left-Right combined ripple distribution around UP-DOWN transition (15 movemean)';
fig.Name = 'Left-Right combined ripple distribution around UP-DOWN transition (10 movemean)'
colour_lines = [0,90,50;74,20,134]/256; % Green Purple


nexttile
event_times = merged_event_info.DOWN_ints;
duration = event_times(:,2) - event_times(:,1);
[~,sorted_index] = sort(duration);
temp_matrix = event_averaging_scale*movmean(ripple_probability.DOWN(sorted_index,:),event_averaging_scale,1,'omitnan');

% Convert duration (in seconds) to number of bins
duration_in_bins = duration(sorted_index) / 0.02;

% Add to center point (DOWN-UP at bin 101)
duration_bin_position = 51 + duration_in_bins;

% for nevent = 1:length(temp_matrix)
%     if round(duration_bin_position(nevent))>0
%     temp_matrix(nevent,round(duration_bin_position(nevent)):end) = nan;
%     end
% end

imagesc(temp_matrix)
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
hold on

% Plot yellow dashed line
plot(flip(duration_bin_position), flip(1:numel(duration)), 'r--', 'LineWidth', 1)

xticks([0.5 13 25.5 38 50.5 62.5 75 82.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('ipsi ripples')
xlim([35 100])

nexttile
clear ERROR_SHADE

binnedArray = ripple_probability.DOWN_boot;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = ripple_probability.DOWN_baseline_boot;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
xticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
% xline(0,'r')
ylim([0 0.11])
% title('ipsi ripples')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'real','baseline'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlim([-0.3 1])

nexttile
clear ERROR_SHADE

num_bins = 50+20;
Range=[-0.2 1.2];
bin_edges = linspace(Range(1), Range(2), num_bins+1);
x = bin_edges(1:end-1) + diff(bin_edges)/2;

binnedArray = ripple_probability.norm_DOWN_boot;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
% 
% binnedArray = probability_merged.contra_ripples_UP{end}{1};
% y = mean(binnedArray,'omitnan');
% LCI = prctile(binnedArray,2.5);
% UCI = prctile(binnedArray,97.5);
% 
% PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
% ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xline(0,'r',LineWidth=1)

% baseline
binnedArray = ripple_probability.norm_baseline_boot;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

% xline(0,'r')
ylim([0 0.085])
xlim([0 1])
xticks([0 0.25 0.5 0.75 1])
% title('ipsi ripples')
xlabel('Normalised DOWN Duration')
ylabel('Probability')
legend([ERROR_SHADE(1:end)],{'real','baseline'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH'),[],'ContentType','vector')


%% By ipsi and contra

load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

% Set up time windows and variables
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

% Extract absolute time probabilities
ipsi_ripples_UP = [probability_psth_whole(1).L_ripples_UP; probability_psth_whole(2).R_ripples_UP];
contra_ripples_UP = [probability_psth_whole(1).R_ripples_UP; probability_psth_whole(2).L_ripples_UP];
ipsi_ripples_DOWN = [probability_psth_whole(1).L_ripples_DOWN; probability_psth_whole(2).R_ripples_DOWN];
contra_ripples_DOWN = [probability_psth_whole(1).R_ripples_DOWN; probability_psth_whole(2).L_ripples_DOWN];

ipsi_ripples_UP_baseline = [probability_psth_whole_baseline(1).L_ripples_UP; probability_psth_whole_baseline(2).R_ripples_UP];
contra_ripples_UP_baseline = [probability_psth_whole_baseline(1).R_ripples_UP; probability_psth_whole_baseline(2).L_ripples_UP];
ipsi_ripples_DOWN_baseline = [probability_psth_whole_baseline(1).L_ripples_DOWN; probability_psth_whole_baseline(2).R_ripples_DOWN];
contra_ripples_DOWN_baseline = [probability_psth_whole_baseline(1).R_ripples_DOWN; probability_psth_whole_baseline(2).L_ripples_DOWN];

% Load ripple_info for relative probability calculation
if ~exist('ripple_info', 'var')
    try
        load(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_info.mat'),'ripple_info');
    catch
        warning('Could not load ripple_info.mat, using all ripples instead.');
        ripple_info = struct();
        ripple_info.event_id = 1:size(merged_event_info.ripples_ints, 1);
    end
end

% Set up bins and range for normalised probability
num_bins = 50+20;
Range = [-0.2 1.2];
bin_edges = linspace(Range(1), Range(2), num_bins+1);
bin_centres = bin_edges(1:end-1) + diff(bin_edges)/2;

% Filter events by hemisphere
L_UP_idx = (merged_event_info.UP_hemisphere_id == 1);
R_UP_idx = (merged_event_info.UP_hemisphere_id == 2);
L_DOWN_idx = (merged_event_info.DOWN_hemisphere_id == 1);
R_DOWN_idx = (merged_event_info.DOWN_hemisphere_id == 2);

L_UP_ints = merged_event_info.UP_ints(L_UP_idx, :);
R_UP_ints = merged_event_info.UP_ints(R_UP_idx, :);
L_DOWN_ints = merged_event_info.DOWN_ints(L_DOWN_idx, :);
R_DOWN_ints = merged_event_info.DOWN_ints(R_DOWN_idx, :);

% Filter ripples by hemisphere
subset_ripple_ids = ripple_info.event_id;
L_ripple_indices = merged_event_info.ripples_hemisphere_id == 1;
R_ripple_indices = merged_event_info.ripples_hemisphere_id == 2;
L_ripples_ints = merged_event_info.ripples_ints(L_ripple_indices, :);
R_ripples_ints = merged_event_info.ripples_ints(R_ripple_indices, :);

% Calculate relative event probabilities
[~,~,~,ripple_count_L_UP_L_ripple] = calculate_relative_event_probability(L_UP_ints, L_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_UP_R_ripple] = calculate_relative_event_probability(R_UP_ints, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_L_UP_R_ripple] = calculate_relative_event_probability(L_UP_ints, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_UP_L_ripple] = calculate_relative_event_probability(R_UP_ints, L_ripples_ints, num_bins, 0, 'Range', Range);

[~,~,~,ripple_count_L_DOWN_L_ripple] = calculate_relative_event_probability(L_DOWN_ints, L_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_DOWN_R_ripple] = calculate_relative_event_probability(R_DOWN_ints, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_L_DOWN_R_ripple] = calculate_relative_event_probability(L_DOWN_ints, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_DOWN_L_ripple] = calculate_relative_event_probability(R_DOWN_ints, L_ripples_ints, num_bins, 0, 'Range', Range);

% Baseline relative event probabilities (+3s shift)
[~,~,~,ripple_count_L_UP_baseline_L_ripple] = calculate_relative_event_probability(L_UP_ints+3, L_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_UP_baseline_R_ripple] = calculate_relative_event_probability(R_UP_ints+3, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_L_UP_baseline_R_ripple] = calculate_relative_event_probability(L_UP_ints+3, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_UP_baseline_L_ripple] = calculate_relative_event_probability(R_UP_ints+3, L_ripples_ints, num_bins, 0, 'Range', Range);

[~,~,~,ripple_count_L_DOWN_baseline_L_ripple] = calculate_relative_event_probability(L_DOWN_ints+3, L_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_DOWN_baseline_R_ripple] = calculate_relative_event_probability(R_DOWN_ints+3, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_L_DOWN_baseline_R_ripple] = calculate_relative_event_probability(L_DOWN_ints+3, R_ripples_ints, num_bins, 0, 'Range', Range);
[~,~,~,ripple_count_R_DOWN_baseline_L_ripple] = calculate_relative_event_probability(R_DOWN_ints+3, L_ripples_ints, num_bins, 0, 'Range', Range);

% Combine into ipsi vs contra for normalised time
ripple_count_ipsi_UP = [ripple_count_L_UP_L_ripple; ripple_count_R_UP_R_ripple];
ripple_count_contra_UP = [ripple_count_L_UP_R_ripple; ripple_count_R_UP_L_ripple];
ripple_count_ipsi_DOWN = [ripple_count_L_DOWN_L_ripple; ripple_count_R_DOWN_R_ripple];
ripple_count_contra_DOWN = [ripple_count_L_DOWN_R_ripple; ripple_count_R_DOWN_L_ripple];

ripple_count_ipsi_UP_baseline = [ripple_count_L_UP_baseline_L_ripple; ripple_count_R_UP_baseline_R_ripple];
ripple_count_contra_UP_baseline = [ripple_count_L_UP_baseline_R_ripple; ripple_count_R_UP_baseline_L_ripple];
ripple_count_ipsi_DOWN_baseline = [ripple_count_L_DOWN_baseline_L_ripple; ripple_count_R_DOWN_baseline_R_ripple];
ripple_count_contra_DOWN_baseline = [ripple_count_L_DOWN_baseline_R_ripple; ripple_count_R_DOWN_baseline_L_ripple];

% Bootstrapping Absolute Time Ripples
temp_ipsi_UP = []; temp_contra_UP = []; temp_ipsi_DOWN = []; temp_contra_DOWN = [];
temp_ipsi_UP_baseline = []; temp_contra_UP_baseline = []; temp_ipsi_DOWN_baseline = []; temp_contra_DOWN_baseline = [];

parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot);
    
    event_id = datasample(s,1:size(ipsi_ripples_UP,1),size(ipsi_ripples_UP,1));
    temp_ipsi_UP(iBoot,:) = sum(ipsi_ripples_UP(event_id,:),'omitnan')./sum(~isnan(ipsi_ripples_UP(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_ripples_UP,1),size(contra_ripples_UP,1));
    temp_contra_UP(iBoot,:) = sum(contra_ripples_UP(event_id,:),'omitnan')./sum(~isnan(contra_ripples_UP(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_ripples_DOWN,1),size(ipsi_ripples_DOWN,1));
    temp_ipsi_DOWN(iBoot,:) = sum(ipsi_ripples_DOWN(event_id,:),'omitnan')./sum(~isnan(ipsi_ripples_DOWN(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_ripples_DOWN,1),size(contra_ripples_DOWN,1));
    temp_contra_DOWN(iBoot,:) = sum(contra_ripples_DOWN(event_id,:),'omitnan')./sum(~isnan(contra_ripples_DOWN(event_id,:)));
    
    % Baselines
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_ripples_UP_baseline,1),size(ipsi_ripples_UP_baseline,1));
    temp_ipsi_UP_baseline(iBoot,:) = sum(ipsi_ripples_UP_baseline(event_id,:),'omitnan')./sum(~isnan(ipsi_ripples_UP_baseline(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_ripples_UP_baseline,1),size(contra_ripples_UP_baseline,1));
    temp_contra_UP_baseline(iBoot,:) = sum(contra_ripples_UP_baseline(event_id,:),'omitnan')./sum(~isnan(contra_ripples_UP_baseline(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_ripples_DOWN_baseline,1),size(ipsi_ripples_DOWN_baseline,1));
    temp_ipsi_DOWN_baseline(iBoot,:) = sum(ipsi_ripples_DOWN_baseline(event_id,:),'omitnan')./sum(~isnan(ipsi_ripples_DOWN_baseline(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_ripples_DOWN_baseline,1),size(contra_ripples_DOWN_baseline,1));
    temp_contra_DOWN_baseline(iBoot,:) = sum(contra_ripples_DOWN_baseline(event_id,:),'omitnan')./sum(~isnan(contra_ripples_DOWN_baseline(event_id,:)));
end

% Bootstrapping Relative (Normalised) Time Ripples
temp_ipsi_UP_norm = []; temp_contra_UP_norm = []; temp_ipsi_DOWN_norm = []; temp_contra_DOWN_norm = [];
temp_ipsi_UP_baseline_norm = []; temp_contra_UP_baseline_norm = []; temp_ipsi_DOWN_baseline_norm = []; temp_contra_DOWN_baseline_norm = [];

parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot);
    
    event_id = datasample(s,1:size(ripple_count_ipsi_UP,1),size(ripple_count_ipsi_UP,1));
    temp_ipsi_UP_norm(iBoot,:) = sum(ripple_count_ipsi_UP(event_id,:),'omitnan')./sum(~isnan(ripple_count_ipsi_UP(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_contra_UP,1),size(ripple_count_contra_UP,1));
    temp_contra_UP_norm(iBoot,:) = sum(ripple_count_contra_UP(event_id,:),'omitnan')./sum(~isnan(ripple_count_contra_UP(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_ipsi_DOWN,1),size(ripple_count_ipsi_DOWN,1));
    temp_ipsi_DOWN_norm(iBoot,:) = sum(ripple_count_ipsi_DOWN(event_id,:),'omitnan')./sum(~isnan(ripple_count_ipsi_DOWN(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_contra_DOWN,1),size(ripple_count_contra_DOWN,1));
    temp_contra_DOWN_norm(iBoot,:) = sum(ripple_count_contra_DOWN(event_id,:),'omitnan')./sum(~isnan(ripple_count_contra_DOWN(event_id,:)));
    
    % Baselines
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_ipsi_UP_baseline,1),size(ripple_count_ipsi_UP_baseline,1));
    temp_ipsi_UP_baseline_norm(iBoot,:) = sum(ripple_count_ipsi_UP_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_count_ipsi_UP_baseline(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_contra_UP_baseline,1),size(ripple_count_contra_UP_baseline,1));
    temp_contra_UP_baseline_norm(iBoot,:) = sum(ripple_count_contra_UP_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_count_contra_UP_baseline(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_ipsi_DOWN_baseline,1),size(ripple_count_ipsi_DOWN_baseline,1));
    temp_ipsi_DOWN_baseline_norm(iBoot,:) = sum(ripple_count_ipsi_DOWN_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_count_ipsi_DOWN_baseline(event_id,:)));
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ripple_count_contra_DOWN_baseline,1),size(ripple_count_contra_DOWN_baseline,1));
    temp_contra_DOWN_baseline_norm(iBoot,:) = sum(ripple_count_contra_DOWN_baseline(event_id,:),'omitnan')./sum(~isnan(ripple_count_contra_DOWN_baseline(event_id,:)));
end

% Save ripple probability structure
ipsi_contra_ripple_probability = struct();
ipsi_contra_ripple_probability.ipsi_UP = ipsi_ripples_UP;
ipsi_contra_ripple_probability.contra_UP = contra_ripples_UP;
ipsi_contra_ripple_probability.ipsi_DOWN = ipsi_ripples_DOWN;
ipsi_contra_ripple_probability.contra_DOWN = contra_ripples_DOWN;
ipsi_contra_ripple_probability.ipsi_UP_baseline = ipsi_ripples_UP_baseline;
ipsi_contra_ripple_probability.contra_UP_baseline = contra_ripples_UP_baseline;
ipsi_contra_ripple_probability.ipsi_DOWN_baseline = ipsi_ripples_DOWN_baseline;
ipsi_contra_ripple_probability.contra_DOWN_baseline = contra_ripples_DOWN_baseline;

ipsi_contra_ripple_probability.ipsi_UP_boot = temp_ipsi_UP;
ipsi_contra_ripple_probability.contra_UP_boot = temp_contra_UP;
ipsi_contra_ripple_probability.ipsi_DOWN_boot = temp_ipsi_DOWN;
ipsi_contra_ripple_probability.contra_DOWN_boot = temp_contra_DOWN;
ipsi_contra_ripple_probability.ipsi_UP_baseline_boot = temp_ipsi_UP_baseline;
ipsi_contra_ripple_probability.contra_UP_baseline_boot = temp_contra_UP_baseline;
ipsi_contra_ripple_probability.ipsi_DOWN_baseline_boot = temp_ipsi_DOWN_baseline;
ipsi_contra_ripple_probability.contra_DOWN_baseline_boot = temp_contra_DOWN_baseline;

ipsi_contra_ripple_probability.norm_ipsi_UP = ripple_count_ipsi_UP;
ipsi_contra_ripple_probability.norm_contra_UP = ripple_count_contra_UP;
ipsi_contra_ripple_probability.norm_ipsi_DOWN = ripple_count_ipsi_DOWN;
ipsi_contra_ripple_probability.norm_contra_DOWN = ripple_count_contra_DOWN;
ipsi_contra_ripple_probability.norm_ipsi_UP_baseline = ripple_count_ipsi_UP_baseline;
ipsi_contra_ripple_probability.norm_contra_UP_baseline = ripple_count_contra_UP_baseline;
ipsi_contra_ripple_probability.norm_ipsi_DOWN_baseline = ripple_count_ipsi_DOWN_baseline;
ipsi_contra_ripple_probability.norm_contra_DOWN_baseline = ripple_count_contra_DOWN_baseline;

ipsi_contra_ripple_probability.norm_ipsi_UP_boot = temp_ipsi_UP_norm;
ipsi_contra_ripple_probability.norm_contra_UP_boot = temp_contra_UP_norm;
ipsi_contra_ripple_probability.norm_ipsi_DOWN_boot = temp_ipsi_DOWN_norm;
ipsi_contra_ripple_probability.norm_contra_DOWN_boot = temp_contra_DOWN_norm;
ipsi_contra_ripple_probability.norm_ipsi_UP_baseline_boot = temp_ipsi_UP_baseline_norm;
ipsi_contra_ripple_probability.norm_contra_UP_baseline_boot = temp_contra_UP_baseline_norm;
ipsi_contra_ripple_probability.norm_ipsi_DOWN_baseline_boot = temp_ipsi_DOWN_baseline_norm;
ipsi_contra_ripple_probability.norm_contra_DOWN_baseline_boot = temp_contra_DOWN_baseline_norm;
ipsi_contra_ripple_probability.norm_ipsi_baseline_boot = (temp_ipsi_UP_baseline_norm + temp_ipsi_DOWN_baseline_norm)./2;

save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ipsi_contra_ripple_probability.mat'),'ipsi_contra_ripple_probability');

% Extract V1 and HPC MUA and group into ipsi vs contra
x_mua = PSTH_MUA(1).timebins;

ipsi_V1_MUA_UP = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA_UP = [PSTH_MUA(1).R_V1_UP; PSTH_MUA(2).L_V1_UP];
ipsi_HPC_MUA_UP = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA_UP = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];

ipsi_V1_MUA_UP_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
contra_V1_MUA_UP_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];
ipsi_HPC_MUA_UP_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_UP_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];

ipsi_V1_MUA_DOWN = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_V1_MUA_DOWN = [PSTH_MUA(1).R_V1_DOWN; PSTH_MUA(2).L_V1_DOWN];
ipsi_HPC_MUA_DOWN = [PSTH_MUA(1).L_HPC_DOWN; PSTH_MUA(2).R_HPC_DOWN];
contra_HPC_MUA_DOWN = [PSTH_MUA(1).R_HPC_DOWN; PSTH_MUA(2).L_HPC_DOWN];

ipsi_V1_MUA_DOWN_baseline = [PSTH_MUA_baseline(1).L_V1_DOWN; PSTH_MUA_baseline(2).R_V1_DOWN];
contra_V1_MUA_DOWN_baseline = [PSTH_MUA_baseline(1).R_V1_DOWN; PSTH_MUA_baseline(2).L_V1_DOWN];
ipsi_HPC_MUA_DOWN_baseline = [PSTH_MUA_baseline(1).L_HPC_DOWN; PSTH_MUA_baseline(2).R_HPC_DOWN];
contra_HPC_MUA_DOWN_baseline = [PSTH_MUA_baseline(1).R_HPC_DOWN; PSTH_MUA_baseline(2).L_HPC_DOWN];

% Bootstrap V1 & HPC MUA for UP and DOWN transitions
boot_ipsi_V1_MUA_UP = []; boot_contra_V1_MUA_UP = []; boot_ipsi_HPC_MUA_UP = []; boot_contra_HPC_MUA_UP = [];
boot_ipsi_V1_MUA_UP_baseline = []; boot_contra_V1_MUA_UP_baseline = []; boot_ipsi_HPC_MUA_UP_baseline = []; boot_contra_HPC_MUA_UP_baseline = [];

boot_ipsi_V1_MUA_DOWN = []; boot_contra_V1_MUA_DOWN = []; boot_ipsi_HPC_MUA_DOWN = []; boot_contra_HPC_MUA_DOWN = [];
boot_ipsi_V1_MUA_DOWN_baseline = []; boot_contra_V1_MUA_DOWN_baseline = []; boot_ipsi_HPC_MUA_DOWN_baseline = []; boot_contra_HPC_MUA_DOWN_baseline = [];

parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot);
    % UP states
    event_id = datasample(s,1:size(ipsi_V1_MUA_UP,1),size(ipsi_V1_MUA_UP,1));
    boot_ipsi_V1_MUA_UP(iBoot,:) = mean(ipsi_V1_MUA_UP(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_V1_MUA_UP,1),size(contra_V1_MUA_UP,1));
    boot_contra_V1_MUA_UP(iBoot,:) = mean(contra_V1_MUA_UP(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_HPC_MUA_UP,1),size(ipsi_HPC_MUA_UP,1));
    boot_ipsi_HPC_MUA_UP(iBoot,:) = mean(ipsi_HPC_MUA_UP(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_HPC_MUA_UP,1),size(contra_HPC_MUA_UP,1));
    boot_contra_HPC_MUA_UP(iBoot,:) = mean(contra_HPC_MUA_UP(event_id,:),1,'omitnan');
    
    % Baselines UP
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_V1_MUA_UP_baseline,1),size(ipsi_V1_MUA_UP_baseline,1));
    boot_ipsi_V1_MUA_UP_baseline(iBoot,:) = mean(ipsi_V1_MUA_UP_baseline(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_V1_MUA_UP_baseline,1),size(contra_V1_MUA_UP_baseline,1));
    boot_contra_V1_MUA_UP_baseline(iBoot,:) = mean(contra_V1_MUA_UP_baseline(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_HPC_MUA_UP_baseline,1),size(ipsi_HPC_MUA_UP_baseline,1));
    boot_ipsi_HPC_MUA_UP_baseline(iBoot,:) = mean(ipsi_HPC_MUA_UP_baseline(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_HPC_MUA_UP_baseline,1),size(contra_HPC_MUA_UP_baseline,1));
    boot_contra_HPC_MUA_UP_baseline(iBoot,:) = mean(contra_HPC_MUA_UP_baseline(event_id,:),1,'omitnan');
   
end



parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot);
    
    % DOWN states
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_V1_MUA_DOWN,1),size(ipsi_V1_MUA_DOWN,1));
    boot_ipsi_V1_MUA_DOWN(iBoot,:) = mean(ipsi_V1_MUA_DOWN(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_V1_MUA_DOWN,1),size(contra_V1_MUA_DOWN,1));
    boot_contra_V1_MUA_DOWN(iBoot,:) = mean(contra_V1_MUA_DOWN(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_HPC_MUA_DOWN,1),size(ipsi_HPC_MUA_DOWN,1));
    boot_ipsi_HPC_MUA_DOWN(iBoot,:) = mean(ipsi_HPC_MUA_DOWN(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_HPC_MUA_DOWN,1),size(contra_HPC_MUA_DOWN,1));
    boot_contra_HPC_MUA_DOWN(iBoot,:) = mean(contra_HPC_MUA_DOWN(event_id,:),1,'omitnan');
    
    % Baselines DOWN
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_V1_MUA_DOWN_baseline,1),size(ipsi_V1_MUA_DOWN_baseline,1));
    boot_ipsi_V1_MUA_DOWN_baseline(iBoot,:) = mean(ipsi_V1_MUA_DOWN_baseline(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_V1_MUA_DOWN_baseline,1),size(contra_V1_MUA_DOWN_baseline,1));
    boot_contra_V1_MUA_DOWN_baseline(iBoot,:) = mean(contra_V1_MUA_DOWN_baseline(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(ipsi_HPC_MUA_DOWN_baseline,1),size(ipsi_HPC_MUA_DOWN_baseline,1));
    boot_ipsi_HPC_MUA_DOWN_baseline(iBoot,:) = mean(ipsi_HPC_MUA_DOWN_baseline(event_id,:),1,'omitnan');
    
    s = RandStream('mrg32k3a','Seed',iBoot);
    event_id = datasample(s,1:size(contra_HPC_MUA_DOWN_baseline,1),size(contra_HPC_MUA_DOWN_baseline,1));
    boot_contra_HPC_MUA_DOWN_baseline(iBoot,:) = mean(contra_HPC_MUA_DOWN_baseline(event_id,:),1,'omitnan');
end

% Save MUA bootstrap results
ipsi_contra_MUA_psth = struct();
ipsi_contra_MUA_psth.timebins = x_mua;

ipsi_contra_MUA_psth.ipsi_V1_MUA_UP_boot = boot_ipsi_V1_MUA_UP;
ipsi_contra_MUA_psth.contra_V1_MUA_UP_boot = boot_contra_V1_MUA_UP;
ipsi_contra_MUA_psth.ipsi_HPC_MUA_UP_boot = boot_ipsi_HPC_MUA_UP;
ipsi_contra_MUA_psth.contra_HPC_MUA_UP_boot = boot_contra_HPC_MUA_UP;

ipsi_contra_MUA_psth.ipsi_V1_MUA_UP_baseline_boot = boot_ipsi_V1_MUA_UP_baseline;
ipsi_contra_MUA_psth.contra_V1_MUA_UP_baseline_boot = boot_contra_V1_MUA_UP_baseline;
ipsi_contra_MUA_psth.ipsi_HPC_MUA_UP_baseline_boot = boot_ipsi_HPC_MUA_UP_baseline;
ipsi_contra_MUA_psth.contra_HPC_MUA_UP_baseline_boot = boot_contra_HPC_MUA_UP_baseline;

ipsi_contra_MUA_psth.ipsi_V1_MUA_DOWN_boot = boot_ipsi_V1_MUA_DOWN;
ipsi_contra_MUA_psth.contra_V1_MUA_DOWN_boot = boot_contra_V1_MUA_DOWN;
ipsi_contra_MUA_psth.ipsi_HPC_MUA_DOWN_boot = boot_ipsi_HPC_MUA_DOWN;
ipsi_contra_MUA_psth.contra_HPC_MUA_DOWN_boot = boot_contra_HPC_MUA_DOWN;

ipsi_contra_MUA_psth.ipsi_V1_MUA_DOWN_baseline_boot = boot_ipsi_V1_MUA_DOWN_baseline;
ipsi_contra_MUA_psth.contra_V1_MUA_DOWN_baseline_boot = boot_contra_V1_MUA_DOWN_baseline;
ipsi_contra_MUA_psth.ipsi_HPC_MUA_DOWN_baseline_boot = boot_ipsi_HPC_MUA_DOWN_baseline;
ipsi_contra_MUA_psth.contra_HPC_MUA_DOWN_baseline_boot = boot_contra_HPC_MUA_DOWN_baseline;

save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ipsi_contra_MUA_psth.mat'),'ipsi_contra_MUA_psth');



% Sort events by duration for heatmaps
L_UP_duration = merged_event_info.UP_ints(L_UP_idx,2) - merged_event_info.UP_ints(L_UP_idx,1);
R_UP_duration = merged_event_info.UP_ints(R_UP_idx,2) - merged_event_info.UP_ints(R_UP_idx,1);
UP_duration_combined = [L_UP_duration; R_UP_duration];
[~, sorted_UP_index] = sort(UP_duration_combined);

L_DOWN_duration = merged_event_info.DOWN_ints(L_DOWN_idx,2) - merged_event_info.DOWN_ints(L_DOWN_idx,1);
R_DOWN_duration = merged_event_info.DOWN_ints(R_DOWN_idx,2) - merged_event_info.DOWN_ints(R_DOWN_idx,1);
DOWN_duration_combined = [L_DOWN_duration; R_DOWN_duration];
[~, sorted_DOWN_index] = sort(DOWN_duration_combined);

%%%%%%%% Plotting Figure 1: Ipsi vs Contra Ripple distribution
%%%%%%%% DOWN UP
event_averaging_scale = 10;
colour_lines = [0,90,50; 74,20,134]/256; % Green (Ipsi), Purple (Contra)


time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

num_bins = 50+20;
Range=[-0.2 1.2];
bin_edges = linspace(Range(1), Range(2), num_bins+1);
x_norm = bin_edges(1:end-1) + diff(bin_edges)/2;


fig1 = figure();
fig1.Position = [100 100 1525 635];
fig1.Name = 'Ipsi vs Contra Ripple distribution around UP/DOWN transitions';

% Subplot 1: Ipsi UP Ripple absolute time heatmap
subplot(2,4,1)
temp_matrix = event_averaging_scale * movmean(ipsi_ripples_UP(sorted_UP_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(h, 'AlphaData', ~isnan(temp_matrix));
set(gca, 'Color', 'w');
hold on;
duration_in_bins = UP_duration_combined(sorted_UP_index) / 0.02;
duration_bin_position = 51 + duration_in_bins;
plot(flip(duration_bin_position), flip(1:numel(UP_duration_combined)), 'r--', 'LineWidth', 1)
xticks([0.5 13 25.5 38 50.5 62.5 75 87.5 100.5])
xticklabels([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xline(50.5,'r',LineWidth=1)
xlim([35 100])
clim([0 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by UP duration')
title('Ipsi Ripples UP')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 2: Contra UP Ripple absolute time heatmap
subplot(2,4,2)
temp_matrix = event_averaging_scale * movmean(contra_ripples_UP(sorted_UP_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(h, 'AlphaData', ~isnan(temp_matrix));
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(UP_duration_combined)), 'r--', 'LineWidth', 1)
xticks([0.5 13 25.5 38 50.5 62.5 75 87.5 100.5])
xticklabels([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xline(50.5,'r',LineWidth=1)
xlim([35 100])
clim([0 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by UP duration')
title('Contra Ripples UP')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 3: Absolute time bootstrap line plot (Ipsi vs Contra vs Baseline)
subplot(2,4,3)
clear ERROR_SHADE

% Ipsi
binnedArray = temp_ipsi_UP;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra
binnedArray = temp_contra_UP;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline
binnedArray = (temp_ipsi_UP_baseline + temp_contra_UP_baseline) ./ 2;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
ylim([0 0.06])
xlim([-0.3 1])
xticks([-0.25:0.25:1])
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Bootstrapped Ripples UP')

% Subplot 4: Normalised time bootstrap line plot
subplot(2,4,4)
clear ERROR_SHADE

% Ipsi norm
binnedArray = temp_ipsi_UP_norm;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_norm,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x_norm fliplr(x_norm)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra norm
binnedArray = temp_contra_UP_norm;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_norm,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x_norm fliplr(x_norm)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline norm
binnedArray = ipsi_contra_ripple_probability.norm_ipsi_baseline_boot;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_norm,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x_norm fliplr(x_norm)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
ylim([0 0.065])
xlim([0 1])
xticks([0 0.25 0.5 0.75 1])
xlabel('Normalised UP Duration')
ylabel('Probability')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Normalised Ripples UP')

% Subplot 5: Ipsi DOWN Ripple absolute time heatmap
subplot(2,4,5)
temp_matrix = event_averaging_scale * movmean(ipsi_ripples_DOWN(sorted_DOWN_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(h, 'AlphaData', ~isnan(temp_matrix));
set(gca, 'Color', 'w');
hold on;
duration_in_bins = DOWN_duration_combined(sorted_DOWN_index) / 0.02;
duration_bin_position = 51 + duration_in_bins;
plot(flip(duration_bin_position), flip(1:numel(DOWN_duration_combined)), 'r--', 'LineWidth', 1)
xticks([0.5 13 25.5 38 50.5 62.5 75 87.5 100.5])
xticklabels([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xline(50.5,'r',LineWidth=1)
xlim([35 100])
clim([0 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
title('Ipsi Ripples DOWN')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 6: Contra DOWN Ripple absolute time heatmap
subplot(2,4,6)
temp_matrix = event_averaging_scale * movmean(contra_ripples_DOWN(sorted_DOWN_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(h, 'AlphaData', ~isnan(temp_matrix));
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(DOWN_duration_combined)), 'r--', 'LineWidth', 1)
xticks([0.5 13 25.5 38 50.5 62.5 75 87.5 100.5])
xticklabels([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xline(50.5,'r',LineWidth=1)
xlim([35 100])
clim([0 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
title('Contra Ripples DOWN')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 7: Absolute time bootstrap line plot (Ipsi vs Contra vs Baseline)
subplot(2,4,7)
clear ERROR_SHADE

% Ipsi
binnedArray = temp_ipsi_DOWN;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra
binnedArray = temp_contra_DOWN;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline
binnedArray = (temp_ipsi_DOWN_baseline + temp_contra_DOWN_baseline) ./ 2;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
ylim([0 0.081])
xlim([-0.3 1])
xticks([-0.25:0.25:1])
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Probability')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Bootstrapped Ripples DOWN')

% Subplot 8: Normalised time bootstrap line plot
subplot(2,4,8)
clear ERROR_SHADE

% Ipsi norm
binnedArray = temp_ipsi_DOWN_norm;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_norm,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x_norm fliplr(x_norm)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra norm
binnedArray = temp_contra_DOWN_norm;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_norm,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x_norm fliplr(x_norm)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline norm
binnedArray = ipsi_contra_ripple_probability.norm_ipsi_baseline_boot;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_norm,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x_norm fliplr(x_norm)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
ylim([0 0.065])
xlim([0 1])
xticks([0 0.25 0.5 0.75 1])
xlabel('Normalised DOWN Duration')
ylabel('Probability')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Normalised Ripples DOWN')

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH'),[],'ContentType','vector')


%%%
% Plotting Figure 2: MUA around UP transition
fig2 = figure();
fig2.Position = [150 150 1400 800];
fig2.Name = 'Ipsi vs Contra MUA around UP transition';

event_averaging_scale = 30;
% Row 1: V1 MUA
% Subplot 1: Ipsi V1 MUA UP Heatmap
subplot(2,3,1)
temp_matrix = movmean(ipsi_V1_MUA_UP(sorted_UP_index,:), event_averaging_scale, 1, 'omitnan');
% temp_matrix = ipsi_V1_MUA_UP;
h = imagesc(temp_matrix);
colormap(flipud(gray))
set(gca, 'Color', 'w');
hold on;
duration_in_bins = UP_duration_combined(sorted_UP_index) / (x_mua(2) - x_mua(1));
duration_bin_position = 101 + duration_in_bins;
plot(flip(duration_bin_position), flip(1:numel(UP_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-2 2])
colorbar;
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Events sorted by UP duration')
title('Ipsi V1 MUA UP')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 2: Contra V1 MUA UP Heatmap
subplot(2,3,2)
temp_matrix = movmean(contra_V1_MUA_UP(sorted_UP_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(UP_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-2 2])
colorbar; colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Events sorted by UP duration')
title('Contra V1 MUA UP')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 3: V1 MUA UP Bootstrap line plot
subplot(2,3,3)
clear ERROR_SHADE

% Ipsi
binnedArray = boot_ipsi_V1_MUA_UP;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra
binnedArray = boot_contra_V1_MUA_UP;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline
binnedArray = (boot_ipsi_V1_MUA_UP_baseline + boot_contra_V1_MUA_UP_baseline) ./ 2;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
xlim([-0.5 0.5])
xlabel('Time relative to transition (s)')
ylabel('MUA Activity')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Bootstrapped V1 MUA UP')

% Row 2: HPC MUA
% Subplot 4: Ipsi HPC MUA UP Heatmap
subplot(2,3,4)
temp_matrix = movmean(ipsi_HPC_MUA_UP(sorted_UP_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(UP_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-1 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Events sorted by UP duration')
title('Ipsi HPC MUA UP')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 5: Contra HPC MUA UP Heatmap
subplot(2,3,5)
temp_matrix = movmean(contra_HPC_MUA_UP(sorted_UP_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(UP_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-1 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Events sorted by UP duration')
title('Contra HPC MUA UP')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 6: HPC MUA UP Bootstrap line plot
subplot(2,3,6)
clear ERROR_SHADE

% Ipsi
binnedArray = boot_ipsi_HPC_MUA_UP;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra
binnedArray = boot_contra_HPC_MUA_UP;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline
binnedArray = (boot_ipsi_HPC_MUA_UP_baseline + boot_contra_HPC_MUA_UP_baseline) ./ 2;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
xlim([-0.3 1])
xticks(-0.25:0.25:1)
ylim([-0.2 0.05])
xlabel('Time relative to transition (s)')
ylabel('MUA Activity')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Bootstrapped HPC MUA UP')


% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_MUA_PSTH_ipsi_contra'),[],'ContentType','vector')

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH'),[],'ContentType','vector')

% Plotting Figure 3: MUA around DOWN transition
fig3 = figure('Color','w');
fig3.Position = [200 200 1400 800];
fig3.Name = 'Ipsi vs Contra MUA around DOWN transition';

% Row 1: V1 MUA
% Subplot 1: Ipsi V1 MUA DOWN Heatmap
subplot(2,3,1)
temp_matrix = movmean(ipsi_V1_MUA_DOWN(sorted_DOWN_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
duration_in_bins = DOWN_duration_combined(sorted_DOWN_index) / (x_mua(2) - x_mua(1));
duration_bin_position = 101 + duration_in_bins;
plot(flip(duration_bin_position), flip(1:numel(DOWN_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-2 2])
colorbar; colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Events sorted by DOWN duration')
title('Ipsi V1 MUA DOWN')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 2: Contra V1 MUA DOWN Heatmap
subplot(2,3,2)
temp_matrix = movmean(contra_V1_MUA_DOWN(sorted_DOWN_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(DOWN_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-2 2])
colorbar; colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Events sorted by DOWN duration')
title('Contra V1 MUA DOWN')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 3: V1 MUA DOWN Bootstrap line plot
subplot(2,3,3)
clear ERROR_SHADE

% Ipsi
binnedArray = boot_ipsi_V1_MUA_DOWN;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra
binnedArray = boot_contra_V1_MUA_DOWN;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline
binnedArray = (boot_ipsi_V1_MUA_DOWN_baseline + boot_contra_V1_MUA_DOWN_baseline) ./ 2;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
xlim([-0.3 1])
xticks(-0.25:0.25:1)
xlabel('Time relative to transition (s)')
ylabel('MUA Activity')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Bootstrapped V1 MUA DOWN')

% Row 2: HPC MUA
% Subplot 4: Ipsi HPC MUA DOWN Heatmap
subplot(2,3,4)
temp_matrix = movmean(ipsi_HPC_MUA_DOWN(sorted_DOWN_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(DOWN_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-1 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Events sorted by DOWN duration')
title('Ipsi HPC MUA DOWN')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 5: Contra HPC MUA DOWN Heatmap
subplot(2,3,5)
temp_matrix = movmean(contra_HPC_MUA_DOWN(sorted_DOWN_index,:), event_averaging_scale, 1, 'omitnan');
h = imagesc(temp_matrix);
set(gca, 'Color', 'w');
hold on;
plot(flip(duration_bin_position), flip(1:numel(DOWN_duration_combined)), 'r--', 'LineWidth', 1)
xticks([1 50 101 150 200])
xticklabels([-1 -0.5 0 0.5 1])
xline(101,'r',LineWidth=1)
clim([-1 1])
colorbar; colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Events sorted by DOWN duration')
title('Contra HPC MUA DOWN')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)

% Subplot 6: HPC MUA DOWN Bootstrap line plot
subplot(2,3,6)
clear ERROR_SHADE

% Ipsi
binnedArray = boot_ipsi_HPC_MUA_DOWN;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(1,:), 'LineWidth', 1.5); hold on;
ERROR_SHADE(1) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.2','LineStyle','none');

% Contra
binnedArray = boot_contra_HPC_MUA_DOWN;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'Color',colour_lines(2,:), 'LineWidth', 1.5);
ERROR_SHADE(2) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.2','LineStyle','none');

% Baseline
binnedArray = (boot_ipsi_HPC_MUA_DOWN_baseline + boot_contra_HPC_MUA_DOWN_baseline) ./ 2;
y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);
PLOT = plot(x_mua,y,'k', 'LineWidth', 1.5);
ERROR_SHADE(3) = patch([x_mua fliplr(x_mua)],[UCI fliplr(LCI)],'k','FaceAlpha','0.2','LineStyle','none');

xline(0,'r',LineWidth=1)
xlim([-0.3 1])
xticks(-0.25:0.25:1)
xlabel('Time relative to transition (s)')
ylabel('MUA Activity')
legend(ERROR_SHADE, {'Ipsi','Contra','Baseline'}, 'box', 'off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',10)
title('Bootstrapped HPC MUA DOWN')

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP_DOWN_ripples_PSTH'),[],'ContentType','vector')
