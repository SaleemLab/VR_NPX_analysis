
function plot_ipsi_contra_UP_DOWN_spindle_probability(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_baseline.mat'));
probability_normalised_whole_baseline = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_baseline.mat'));
probability_normalised_baseline = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_baseline.mat'));
probability_psth_baseline = probability;


load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised.mat'));
probability_normalised = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));
probability_psth = probability;


%% Ripple probability during UP-DOWN onset
probability= probability_psth_whole;
probability_baseline = probability_psth_whole_baseline;
event_averaging_scale = 30;

% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'ipsi-contra ripple distribution around DOWN-UP transition';

nprobe = 1;
L_ripples = probability(nprobe).L_ripples_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
imagesc(30*movmean(L_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral ripples')


nprobe = 2;
R_ripples = probability(nprobe).R_ripples_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(50*movmean(R_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral ripples')


nprobe = 1;
R_ripples = probability(nprobe).R_ripples_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(30*movmean(R_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral ripples')



nprobe = 2;
L_ripples = probability(nprobe).L_ripples_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(30*movmean(L_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral ripples')


nexttile
nprobe = 1;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability(nprobe).L_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).R_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during left DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
nprobe = 2;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(probability(nprobe).R_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).L_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);


PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability_baseline(nprobe).R_ripples_UP_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during Right DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%% Ripple probability during DOWN-UP transition
probability= probability_psth_whole;
probability_baseline = probability_psth_whole_baseline;
event_averaging_scale = 30;
% colour_lines = [44,123,182;215,25,28]/256;
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'ipsi-contra ripple distribution around UP-DOWN transition';

nprobe = 1;
L_ripples = probability(nprobe).L_ripples_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,1));
imagesc(event_averaging_scale*movmean(L_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by UP duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral ripples')

nprobe = 2;
R_ripples = probability(nprobe).R_ripples_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(R_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral ripples')


nprobe = 1;
R_ripples = probability(nprobe).R_ripples_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(R_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral ripples')


nprobe = 2;
L_ripples = probability(nprobe).L_ripples_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability(nprobe).UP_all_index,1));
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(L_ripples(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral ripples')





nexttile
nprobe = 1;
all_UP_no = length(probability(nprobe).UP_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability(nprobe).L_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).R_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(probability_baseline(nprobe).L_ripples_UP_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during left UP')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
nprobe = 2;
all_DOWN_no = length(probability(nprobe).UP_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability(nprobe).R_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).L_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability_baseline(nprobe).R_ripples_UP_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during Right UP')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%% Ripple probability during Normalised UP and DOWN
probability= probability_normalised_whole;
probability_baseline = probability_normalised_whole_baseline;


num_bins = size(probability(1).L_ripples_DOWN,2);

colour_lines = [44,123,182;215,25,28]/256;
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'ipsi-contra ripple distribution during normalised UP';

%%%%% DOWN
nexttile
nprobe = 1;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = 1:num_bins
y = mean(probability(nprobe).L_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).R_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during Left DOWN')
xlabel('Normalised duration of DOWN')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
nprobe = 2;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
y = mean(probability(nprobe).R_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);


PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).L_ripples_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);
PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during Right DOWN')
xlabel('Normalised duration of DOWN')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


%%%%% UP
nexttile
nprobe = 1;
all_DOWN_no = length(probability(nprobe).UP_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = 1:num_bins
y = mean(probability(nprobe).L_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).R_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability_baseline(nprobe).L_ripples_UP_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during Left UP')
xlabel('Normalised duration of UP')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
nprobe = 2;
all_DOWN_no = length(probability(nprobe).UP_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
y = mean(probability(nprobe).R_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);


PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability(nprobe).L_ripples_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);
PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(probability_baseline(nprobe).R_ripples_UP_bootstrap,'omitnan');
LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi ripples','contra ripples','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of ripples during Right UP')
xlabel('Normalised duration of UP')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])



%%%%%%%
% Phase-amplitude coupling
nBins = 18;
edges = linspace(-pi, pi, nBins+1);

for nchannel = 1:size(SO_phase_ripples,1)
    for mchannel = 1:size(ripple_peak_amplitude,1)
        % Phase-amplitude coupling
        SO_phase_ripples(:,nevent) = ripple_peak_amplitude(tidx,:);
        spindle_phase_ripples(:,nevent) = ripple_peak_amplitude(tidx,:);


        [~,~,binIdx] = histcounts(SO_phase_ripples(nchannel,:), edges);
        % Mean amplitude in each phase bin
        ampByPhase = accumarray(binIdx(binIdx>0), rippleAmp(binIdx>0), [nBins 1], @mean);
        p = ampByPhase / sum(ampByPhase); % normalize to get probability distribution
    end
end
if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])


% SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));