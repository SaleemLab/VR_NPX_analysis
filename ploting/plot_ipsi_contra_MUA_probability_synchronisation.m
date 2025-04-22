function plot_ipsi_contra_MUA_probability_synchronisation(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,event_info,sessions_to_process)

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


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'));
% probability_SO_SO = probability;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'));
% probability_SO_SO_contralateral = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'));
probability_SO_SO_contralateral = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'));
probability_SO_SO = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability.mat'));
probability_SO_SO_contralateral_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_whole.mat'));
probability_SO_SO_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability_whole.mat'));
probability_ripples_ripples_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;

colour_lines = [];
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red for R
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue for L
colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral

colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

colour_lines = [44,123,182;215,25,28]/256;


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


%% Left and Right ipsi-contra DOWN probability

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_whole_baseline.mat'));
probability_SO_SO_whole_baseline = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_whole.mat'));
probability_SO_SO_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability_baseline.mat'));
probability_SO_SO_contralatral_whole_baseline = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability.mat'));
probability_SO_SO_contralatral_whole = probability;

probability= probability_SO_SO_whole;
% probability_baseline = probability_SO_SO_whole_baseline;


% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'ipsi-contra DOWN distribution around UP-DOWN transition';

nprobe = 1;
event_averaging_scale = 5;
L_UP = probability_SO_SO_whole(nprobe).DOWN_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,1));
imagesc(30*movmean(L_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Ipsi DOWN probability')


nprobe = 2;
R_UP = probability_SO_SO_whole(nprobe).DOWN_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(50*movmean(R_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Ipsi DOWN probability')


nprobe = 1;
R_UP = probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(R_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Contra DOWN probabilty')



nprobe = 2;
L_UP = probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability_SO_SO_whole(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(L_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Contra DOWN probabilty')


nexttile
nprobe = 1;
all_DOWN_no = length(probability_SO_SO_whole(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability_SO_SO_whole(nprobe).DOWN_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_whole(nprobe).DOWN_DOWN_bootstrap,2.5);
UCI = prctile(probability_SO_SO_whole(nprobe).DOWN_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(probability_SO_SO_contralatral_whole_baseline(nprobe).DOWN_DOWN_bootstrap,'omitnan');
LCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).DOWN_DOWN_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).DOWN_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of DOWN during left DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
nprobe = 2;
all_DOWN_no = length(probability_SO_SO_whole(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability_SO_SO_whole(nprobe).DOWN_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_whole(nprobe).DOWN_DOWN_bootstrap,2.5);
UCI = prctile(probability_SO_SO_whole(nprobe).DOWN_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralateral_whole(nprobe).DOWN_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(probability_SO_SO_contralatral_whole_baseline(nprobe).DOWN_DOWN_bootstrap,'omitnan');
LCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).DOWN_DOWN_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).DOWN_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
title('Probability of DOWN during Right DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])

%% Left and Right ipsi-contra UP probability
% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'ipsi-contra UP distribution around DOWN-UP transition';

nprobe = 1;
event_averaging_scale = 5;
L_UP = probability_SO_SO_whole(nprobe).UP_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,1));
imagesc(30*movmean(L_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Ipsi UP probability')


nprobe = 2;
R_UP = probability_SO_SO_whole(nprobe).UP_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,1));% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(50*movmean(R_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Ipsi UP probability')


nprobe = 1;
R_UP = probability_SO_SO_contralateral_whole(nprobe).UP_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,1));% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(R_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Contra UP probabilty')



nprobe = 2;
L_UP = probability_SO_SO_contralateral_whole(nprobe).UP_UP;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,2)-slow_waves_all(nprobe).UP_ints(probability_SO_SO_whole(nprobe).UP_all_index,1));% imagesc(movmean(50*movmean(L_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(L_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
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
title('Contra UP probabilty')


nexttile
nprobe = 1;
all_DOWN_no = length(probability_SO_SO_whole(nprobe).UP_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability_SO_SO_whole(nprobe).UP_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_whole(nprobe).UP_UP_bootstrap,2.5);
UCI = prctile(probability_SO_SO_whole(nprobe).UP_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(probability_SO_SO_contralateral_whole(nprobe).UP_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_contralateral_whole(nprobe).UP_UP_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralateral_whole(nprobe).UP_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(probability_SO_SO_contralatral_whole_baseline(nprobe).UP_UP_bootstrap,'omitnan');
LCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).UP_UP_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).UP_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
% xline(0,'r')
title('Probability of UP during left UP')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
nprobe = 2;
all_DOWN_no = length(probability_SO_SO_whole(nprobe).UP_all_index);
time_wondows = [-1 1];
time_bin = 0.02;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(probability_SO_SO_whole(nprobe).UP_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_whole(nprobe).UP_UP_bootstrap,2.5);
UCI = prctile(probability_SO_SO_whole(nprobe).UP_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(probability_SO_SO_contralateral_whole(nprobe).UP_UP,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(probability_SO_SO_contralateral_whole(nprobe).UP_UP_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralateral_whole(nprobe).UP_UP_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(probability_SO_SO_contralatral_whole_baseline(nprobe).UP_UP_bootstrap,'omitnan');
LCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).UP_UP_bootstrap,2.5);
UCI = prctile(probability_SO_SO_contralatral_whole_baseline(nprobe).UP_UP_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
title('Probability of UP during Right UP')
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Probability')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%%%%%
%%%%%
%%%%%
%%%%%
%% Grabbing ipsi and contra values
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');


probability = probability_psth_whole;
nprobe = 1;
event_averaging_scale = 10;
% extract plv, amp corr and lag (latency)
for nprobe=1:2
    mprobe = abs(nprobe-3);

    % UP and DOWN
    ipsi_lag_DU{nprobe} = [];
    contra_lag_DU{nprobe} = [];
    ipsi_lag_UD{nprobe} = [];
    contra_lag_UD{nprobe} = [];

    ipsi_corr_DU{nprobe} = [];
    contra_corr_DU{nprobe} = [];
    ipsi_corr_UD{nprobe} = [];
    contra_corr_UD{nprobe} = [];

    ipsi_plv_DU{nprobe} = [];
    contra_plv_DU{nprobe} = [];
    ipsi_plv_UD{nprobe} = [];
    contra_plv_UD{nprobe} = [];

    % ripples
    ipsi_lag_ripples{nprobe} = [];
    contra_lag_ripples{nprobe} = [];

    ipsi_corr_ripples{nprobe} = [];
    contra_corr_ripples{nprobe} = [];

    ipsi_plv_ripples{nprobe} = [];
    contra_plv_ripples{nprobe} = [];

    % spindles
    ipsi_lag_spindles{nprobe} = [];
    contra_lag_spindles{nprobe} = [];

    ipsi_corr_spindles{nprobe} = [];
    contra_corr_spindles{nprobe} = [];

    ipsi_plv_spindles{nprobe} = [];
    contra_plv_spindles{nprobe} = [];

    for nsession = 1:max(ripples_all(1).session_count)

        %%%%%% DOWN -> UP
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession)),probability(nprobe).UP_all_index);

        ipsi_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        contra_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);

        ipsi_lag_DU{nprobe} = [ipsi_lag_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_DU{nprobe} = [contra_lag_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_DU{nprobe} = [ipsi_plv_DU{nprobe} squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_DU{nprobe} = [contra_plv_DU{nprobe} squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_DU{nprobe} = [ipsi_corr_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_DU{nprobe} = [contra_corr_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        % ipsi_lag_DU{nprobe} = [ipsi_lag_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_lag_DU{nprobe} = [contra_lag_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_plv_DU{nprobe} = [ipsi_plv_DU{nprobe} min(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_plv_DU{nprobe} = [contra_plv_DU{nprobe} min(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_corr_DU{nprobe} = [ipsi_corr_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_corr_DU{nprobe} = [contra_corr_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        %%%%%% UP -> DOWN
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession)),probability(nprobe).DOWN_all_index);

        ipsi_lag_UD{nprobe} = [ipsi_lag_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_UD{nprobe} = [contra_lag_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_UD{nprobe} = [ipsi_plv_UD{nprobe} squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_UD{nprobe} = [contra_plv_UD{nprobe} squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_UD{nprobe} = [ipsi_corr_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_UD{nprobe} = [contra_corr_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        % ipsi_lag_UD{nprobe} = [ipsi_lag_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % contra_lag_UD{nprobe} = [contra_lag_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_plv_UD{nprobe} = [ipsi_plv_UD{nprobe} min(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_plv_UD{nprobe} = [contra_plv_UD{nprobe} min(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_corr_UD{nprobe} = [ipsi_corr_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_corr_UD{nprobe} = [contra_corr_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        %%%%% Spindles
        [C,ia,ib] = intersect(find(spindles_all(nprobe).session_count == sessions_to_process(nsession)),find(spindles_all(nprobe).session_count == sessions_to_process(nsession) & spindles_all(nprobe).SWS_index == 1));

        if ~isempty(ia)
            ipsi_lag_spindles{nprobe} = [ipsi_lag_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            contra_lag_spindles{nprobe} = [contra_lag_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

            ipsi_plv_spindles{nprobe} = [ipsi_plv_spindles{nprobe} squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            contra_plv_spindles{nprobe} = [contra_plv_spindles{nprobe} squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

            ipsi_corr_spindles{nprobe} = [ipsi_corr_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            contra_corr_spindles{nprobe} = [contra_corr_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];
        end

        %%%%% Ripples
        [C,ia,ib] = intersect(find(ripples_all(nprobe).session_count == sessions_to_process(nsession)),find(ripples_all(nprobe).session_count == sessions_to_process(nsession) & ripples_all(nprobe).SWS_index == 1));
        ipsi_shank = find(ripples_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==HPC_ref_shank(nsession,nprobe))=[];
        contra_shank = find(ripples_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);


        ipsi_lag_ripples{nprobe} = [ipsi_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_ripples{nprobe} = [contra_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_ripples{nprobe} = [ipsi_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_ripples{nprobe} = [contra_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_ripples{nprobe} = [ipsi_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_ripples{nprobe} = [contra_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];
    end
end

%% Plot DU transition averaged MUA
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_event_info.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

ipsi_V1_MUA = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA = [PSTH_MUA(2).L_V1_UP; PSTH_MUA(1).R_V1_UP];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(2).L_V1_UP; PSTH_MUA_baseline(1).R_V1_UP];

ipsi_probability = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_probability = [PSTH_MUA(2).L_V1_UP; PSTH_MUA(1).R_V1_UP];

lag_diff = merged_event_info.UP_lag_diff;
plv_diff = merged_event_info.UP_plv_diff;
corr_diff = merged_event_info.UP_corr_diff;
ipsi_lag = [ipsi_lag_DU{1} ipsi_lag_DU{2}]';
contra_lag = [contra_lag_DU{1} contra_lag_DU{2}]';
group_id = merged_event_info.UP_group_id;

sync_threshold = mean(abs([event_info(1).UP_lag_threshold_low event_info(2).UP_lag_threshold_low event_info(1).UP_lag_threshold_high event_info(2).UP_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Contra
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;


event_idx = [];
event_idx{1} = {intersect(find((group_id == 3 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < 0.1 & lag_diff > -0.1 ),merged_event_info.DOWN_index),intersect(find((group_id == 1 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),merged_event_info.DOWN_index)};

event_idx{2} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.UP_index),intersect(find(group_id == 3 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.UP_index),intersect(find(group_id == 1 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.UP_index)};

event_idx{3} = {intersect(find(group_id == 1 & lag_diff <prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.UP_index),intersect(find(group_id == 1 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 3 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.UP_index),intersect(find(group_id == 3 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.UP_index)};

event_idx{4} = {intersect(find(group_id ~= 2 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.UP_index),intersect(find(group_id ~= 2 &  lag_diff <0 &lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff <sync_threshold & lag_diff >-sync_threshold),merged_event_info.UP_index),...
    intersect(find(group_id ~= 2 & lag_diff>0 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.UP_index),intersect(find(group_id ~= 2 & lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.UP_index)}

event_idx{5} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),merged_event_info.UP_index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1 & corr_diff <prctile(corr_diff(group_id == 1),50) ),merged_event_info.UP_index),intersect(find(group_id == 1 & corr_diff >prctile(corr_diff(group_id == 1),50) ),merged_event_info.UP_index)};

event_idx{6} = {intersect(find(group_id == 3 & lag_diff <0 ),merged_event_info.UP_index),intersect(find(group_id == 3 & lag_diff >0 ),merged_event_info.UP_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1 & lag_diff <0 ),merged_event_info.UP_index),intersect(find(group_id == 1 & lag_diff >0 ),merged_event_info.UP_index)};


group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% non-dominant clusters
group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};

group_name{6} = {'Ipsi dominant ipsi leading','Ipsi dominant contra leading','Bilaterally synchronised','Contra dominant ipsi leading','Contra dominant contra leading','Shuffled'};% dominant clusters

title_names = {'Ipsi-contra DOWN-UP V1 MUA by three clusters', 'Ipsi-contra DOWN-UP V1 MUA by dominant clusters + 50% lag diff','Ipsi-contra DOWN-UP V1 MUA by non-dominant clusters + 50% lag diff',...
    'Ipsi-contra DOWN-UP V1 MUA by cluster 1 and 3 merged + 50% lag diff','Ipsi-contra DOWN-UP V1 MUA by corr diff','Ipsi-contra DOWN-UP V1 MUA by 5 clusters'}
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [70 300 1700 500];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-1 0.5])
    title('Ipsi MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = contra_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(contra_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(contra_baseline_bootstrap,2.5);
    UCI = prctile(contra_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])

    % xline(0,'r')
    title('Contra MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_contra_diff_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_contra_diff_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_contra_diff_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])


    % xline(0,'r')
    title('Ipsi-contra MUA diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    
end


%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])




%% Plot DU transition averaged MUA
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_event_info.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

ipsi_V1_MUA = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_V1_MUA = [PSTH_MUA(2).L_V1_DOWN; PSTH_MUA(1).R_V1_DOWN];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_DOWN; PSTH_MUA_baseline(2).R_V1_DOWN];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(2).L_V1_DOWN; PSTH_MUA_baseline(1).R_V1_DOWN];

ipsi_probability = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_probability = [PSTH_MUA(2).L_V1_DOWN; PSTH_MUA(1).R_V1_DOWN];

lag_diff = merged_event_info.DOWN_lag_diff;
plv_diff = merged_event_info.DOWN_plv_diff;
corr_diff = merged_event_info.DOWN_corr_diff;
ipsi_lag = [ipsi_lag_UD{1} ipsi_lag_UD{2}]';
contra_lag = [contra_lag_UD{1} contra_lag_UD{2}]';
group_id = merged_event_info.UP_group_id;

sync_threshold = mean(abs([event_info(1).DOWN_lag_threshold_low event_info(2).DOWN_lag_threshold_low event_info(1).DOWN_lag_threshold_high event_info(2).DOWN_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Contra
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;


event_idx = [];
event_idx{1} = {intersect(find((group_id == 3 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < 0.1 & lag_diff > -0.1 ),merged_event_info.DOWN_index),intersect(find((group_id == 1 & lag_diff < -0.1)|(group_id == 1 & lag_diff > 0.1) ),merged_event_info.DOWN_index)};

event_idx{2} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 3 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 1 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 1 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.DOWN_index)};

event_idx{3} = {intersect(find(group_id == 1 & lag_diff <prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 1 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 3 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 3 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.DOWN_index)};

event_idx{4} = {intersect(find(group_id ~= 2 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.DOWN_index),intersect(find(group_id ~= 2 &  lag_diff <0 &lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff <sync_threshold & lag_diff >-sync_threshold),merged_event_info.DOWN_index),...
    intersect(find(group_id ~= 2 & lag_diff>0 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.DOWN_index),intersect(find(group_id ~= 2 & lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.DOWN_index)}

event_idx{5} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 1 & corr_diff <prctile(corr_diff(group_id == 1),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 1 & corr_diff >prctile(corr_diff(group_id == 1),50) ),merged_event_info.DOWN_index)};

event_idx{6} = {intersect(find(group_id == 3 & lag_diff <0 ),merged_event_info.DOWN_index),intersect(find(group_id == 3 & lag_diff >0 ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 1 & lag_diff <0 ),merged_event_info.DOWN_index),intersect(find(group_id == 1 & lag_diff >0 ),merged_event_info.DOWN_index)};



group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% non-dominant clusters
group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};

group_name{6} = {'Ipsi dominant ipsi leading','Ipsi dominant contra leading','Bilaterally synchronised','Contra dominant ipsi leading','Contra dominant contra leading','Shuffled'};% dominant clusters

title_names = {'Ipsi-contra UP_DOWN V1 MUA by three clusters', 'Ipsi-contra UP_DOWN V1 MUA by dominant clusters + 50% lag diff','Ipsi-contra UP_DOWN V1 MUA by non-dominant clusters + 50% lag diff',...
    'Ipsi-contra UP_DOWN V1 MUA by cluster 1 and 3 merged + 50% lag diff','Ipsi-contra UP_DOWN V1 MUA by corr diff','Ipsi-contra UP_DOWN V1 MUA by 5 clusters'}
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [70 300 1700 500];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-1 0.5])
    title('Ipsi MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = contra_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(contra_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(contra_baseline_bootstrap,2.5);
    UCI = prctile(contra_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])

    % xline(0,'r')
    title('Contra MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_contra_diff_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_contra_diff_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_contra_diff_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])


    % xline(0,'r')
    title('Ipsi-contra MUA diff')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    
end



%%
%%%%%
%%%%%
%%%%%
%%%%%

%% Plot ripples averaged MUA
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_event_info.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_ripples; PSTH_MUA(2).R_HPC_ripples];
contra_HPC_MUA = [PSTH_MUA(2).L_HPC_ripples; PSTH_MUA(1).R_HPC_ripples];

ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_ripples; PSTH_MUA_baseline(2).R_HPC_ripples];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(2).L_HPC_ripples; PSTH_MUA_baseline(1).R_HPC_ripples];
% 
% ipsi_probability = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
% contra_probability = [PSTH_MUA(2).L_V1_UP; PSTH_MUA(1).R_V1_UP];

lag_diff = merged_event_info.ripples_lag_diff;
plv_diff = merged_event_info.ripples_plv_diff;
corr_diff = merged_event_info.ripples_corr_diff;
ipsi_lag = [ipsi_lag_ripples{1} ipsi_lag_ripples{2}]';
contra_lag = [contra_lag_ripples{1} contra_lag_ripples{2}]';
group_id = merged_event_info.ripples_group_id;

sync_threshold = mean(abs([event_info(1).ripples_lag_threshold_low event_info(2).ripples_lag_threshold_low event_info(1).ripples_lag_threshold_high event_info(2).ripples_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Contra
binnedArray = ipsi_HPC_MUA_baseline-contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;


event_idx = [];
event_idx{1} = {find(group_id==1),find(group_id==2),find(group_id==3),find(group_id==4)};



group_name=[];
group_name{1} = {'contra dominant','low plv diff','high plv diff','ipsi dominant','Shuffled'};

title_names = {'Ipsi-contra ripples HPC MUA by four clusters'}
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [400 150 1300 500];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.5 0.5])
    title('Ipsi MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = contra_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(contra_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(contra_baseline_bootstrap,2.5);
    UCI = prctile(contra_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.5 0.5])


    % xline(0,'r')
    title('Contra MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        index =event_idx{ngroup}{i};
        binnedArray = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
        temp=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        y = mean(binnedArray,'omitnan');
        LCI = prctile(temp,2.5);
        UCI = prctile(temp,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline

    y = mean(ipsi_contra_diff_baseline_bootstrap,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(ipsi_contra_diff_baseline_bootstrap,2.5);
    UCI = prctile(ipsi_contra_diff_baseline_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.5 0.5])


    % xline(0,'r')
    title('Ipsi-contra MUA diff')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    
end