function plot_UP_DOWN_ripple_MUA_PSTH

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
% load(fullfile(analysis_folder,'ripples_all_POST.mat'))
% load(fullfile(analysis_folder,'spindles_all_POST.mat'))
% load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));

%% MUA during DOWN
colour_lines = [44,123,182;215,25,28]/256;
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'Left UP-DOWN transition MUA activity';

nprobe = 1;
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

nexttile
[~,sorted_index] = sort(slow_waves_all(1).DOWN_ints(probability(1).DOWN_all_index,2)-slow_waves_all(1).DOWN_ints(probability(1).DOWN_all_index,1));
imagesc(movmean(PSTH_MUA(1).L_V1_DOWN(sorted_index,:),50,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Left V1 MUA spiking')

nexttile
% [~,sorted_index] = sort(slow_waves_all(1).DOWN_ints(probability(1).DOWN_all_index,2)-slow_waves_all(1).DOWN_ints(probability(1).DOWN_all_index,1));
imagesc(movmean(PSTH_MUA(1).L_HPC_DOWN(sorted_index,:),50,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Left HPC MUA spiking')


nexttile
% [~,sorted_index] = sort(slow_waves_all(1).UP_ints(probability(1).UP_all_index,2)-slow_waves_all(1).UP_ints(probability(1).UP_all_index,1));
imagesc(movmean(PSTH_MUA(1).R_V1_DOWN(sorted_index,:),50,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Right V1 MUA spiking')


nexttile
% [~,sorted_index] = sort(slow_waves_all(1).UP_ints(probability(1).UP_all_index,2)-slow_waves_all(1).UP_ints(probability(1).UP_all_index,1));
imagesc(movmean(PSTH_MUA(1).R_HPC_DOWN(sorted_index,:),50,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Right HPC MUA spiking')



nexttile
x = PSTH_MUA(nprobe).timebins;
y = mean(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap);
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
x = PSTH_MUA(nprobe).timebins;
y = mean(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap);
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
x = PSTH_MUA(nprobe).timebins;
y = mean(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap);
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
x = PSTH_MUA(nprobe).timebins;
y = mean(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap);
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
sgtitle('Left V1 UP-DOWN transition')


if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%% MUA during UP
% probability = probability_psth;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;


% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left DOWN-UP transition with shuffle','MUA activity during right DOWN-UP transition with shuffle'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [1355 332 700 600];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_UP_no = size(PSTH_MUA(nprobe).L_V1_UP,1);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left V1 
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_V1_UP_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).L_V1_UP_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_V1_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).L_V1_UP_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).L_V1_UP_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_V1_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Left V1')
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% Right V1
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_V1_UP_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).R_V1_UP_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_V1_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(4,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(4,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).R_V1_UP_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).R_V1_UP_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_V1_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Right V1')
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Left HPC
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_HPC_UP_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).L_HPC_UP_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Left HPC')
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%% Right HPC
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_HPC_UP_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(3,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(3,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).R_HPC_UP_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Right HPC')
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end



%% MUA during DOWN
% probability = probability_psth;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;


% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left UP-DOWN transition with shuffle','MUA activity during right UP-DOWN transition with shuffle'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [1355 332 700 600];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left V1 
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).L_V1_DOWN_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Left V1')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% Right V1
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(4,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(4,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).R_V1_DOWN_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Right V1')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Left HPC
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).L_HPC_DOWN_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Left HPC')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%% Right HPC
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(3,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(3,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).R_HPC_DOWN_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Right HPC')
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


%% MUA during ripples
% probability = probability_psth;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;


% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left ripples with shuffles','MUA activity during right ripples with shuffles'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [1355 332 700 600];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left V1 
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_V1_ripples_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).L_V1_ripples_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_V1_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).L_V1_ripples_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).L_V1_ripples_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_V1_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Left V1')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% Right V1
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_V1_ripples_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).R_V1_ripples_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_V1_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(4,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(4,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).R_V1_ripples_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).R_V1_ripples_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_V1_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Right V1')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Left HPC
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_HPC_ripples_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).L_HPC_ripples_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_HPC_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).L_HPC_ripples_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).L_HPC_ripples_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).L_HPC_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Left HPC')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%% Right HPC
    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_HPC_ripples_bootstrap);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(PSTH_MUA(nprobe).R_HPC_ripples_bootstrap,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_HPC_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(3,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(3,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(PSTH_MUA(nprobe).R_HPC_ripples_shuffled);
    LCI = prctile(PSTH_MUA(nprobe).R_HPC_ripples_shuffled,2.5);
    UCI = prctile(PSTH_MUA(nprobe).R_HPC_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Right HPC')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


%% MUA during UP
% probability = probability_psth;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;


% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left DOWN-UP transition','MUA activity during right DOWN-UP transition'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
fig(nprobe)=figure;
fig(nprobe).Position = [500 150 980 760];
fig(nprobe).Name = probe_hemisphere_texts{nprobe};

for nprobe = 1:2

    all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left V1
    nexttile
    no_events = size(PSTH_MUA(nprobe).L_V1_UP,1);
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_V1_UP);
    LCI = y - std(PSTH_MUA(nprobe).L_V1_UP)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).L_V1_UP)./sqrt(no_events);

    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).L_V1_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).L_V1_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_V1_UP);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).R_V1_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).R_V1_UP_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).R_V1_UP)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).R_V1_UP)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(4,:));hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(4,:),'FaceAlpha','0.3','LineStyle','none');

    legend(ERROR_SHADE(1:2),{'Left V1','Right V1'},'Box','off','Location','southeast')
    % xline(0,'r')
    title(probe_hemisphere_texts{nprobe})
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_HPC_UP);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).L_HPC_UP)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).L_HPC_UP)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
    ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_HPC_UP);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).R_HPC_UP)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).R_HPC_UP)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(3,:));hold on;
    ERROR_SHADE(4) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(3,:),'FaceAlpha','0.3','LineStyle','none');


    legend(ERROR_SHADE(3:4),{'Left HPC','Right HPC'},'Box','off','Location','southeast')
    % xline(0,'r')
    title(probe_hemisphere_texts{nprobe})
    xlabel('Time relative to UP onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    % nexttile
end


%% MUA during DOWN
% probability = probability_psth;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;


% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left UP-DOWN transition','MUA activity during right UP-DOWN transition'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
fig(nprobe)=figure;
fig(nprobe).Position = [500 150 980 760];
fig(nprobe).Name = probe_hemisphere_texts{nprobe};

for nprobe = 1:2
% 
    % all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left V1
    nexttile
    no_events = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_V1_DOWN);
    LCI = y - std(PSTH_MUA(nprobe).L_V1_DOWN)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).L_V1_DOWN)./sqrt(no_events);

    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).L_V1_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).L_V1_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_V1_DOWN);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).R_V1_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).R_V1_UP_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).R_V1_DOWN)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).R_V1_DOWN)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(4,:));hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(4,:),'FaceAlpha','0.3','LineStyle','none');

    legend(ERROR_SHADE(1:2),{'Left V1','Right V1'},'Box','off','Location','southeast')
    % xline(0,'r')
    title(probe_hemisphere_texts{nprobe})
    xlabel('Time relative to DOWN onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_HPC_DOWN);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).L_HPC_UP_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).L_HPC_DOWN)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).L_HPC_DOWN)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
    ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_HPC_DOWN);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).R_HPC_UP_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).R_HPC_DOWN)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).R_HPC_DOWN)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(3,:));hold on;
    ERROR_SHADE(4) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(3,:),'FaceAlpha','0.3','LineStyle','none');


    legend(ERROR_SHADE(3:4),{'Left HPC','Right HPC'},'Box','off','Location','southeast')
    % xline(0,'r')
    title(probe_hemisphere_texts{nprobe})
    xlabel('Time relative to DOWN onset (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    % nexttile
end

%% MUA during ripples
% probability = probability_psth;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;


% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left ripples','MUA activity during right ripples'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
fig(nprobe)=figure;
fig(nprobe).Position = [500 150 980 760];
fig(nprobe).Name = probe_hemisphere_texts{nprobe};

for nprobe = 1:2

    all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left V1 
    nexttile
    no_events = size(PSTH_MUA(nprobe).L_V1_ripples,1);
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_V1_ripples);
    LCI = y - std(PSTH_MUA(nprobe).L_V1_ripples)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).L_V1_ripples)./sqrt(no_events);

    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).L_V1_ripples_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).L_V1_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_V1_ripples);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).R_V1_ripples_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).R_V1_ripples_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).R_V1_ripples)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).R_V1_ripples)./sqrt(no_events);

    PLOT = plot(x,y,'Color',colour_lines(4,:));hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(4,:),'FaceAlpha','0.3','LineStyle','none');

    legend(ERROR_SHADE(1:2),{'Left V1','Right V1'},'Box','off')
    % xline(0,'r')
    title(probe_hemisphere_texts{nprobe})
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).L_HPC_ripples);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).L_HPC_ripples_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).L_HPC_ripples_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).L_HPC_ripples)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).L_HPC_ripples)./sqrt(no_events);

    PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
    ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

    x = PSTH_MUA(nprobe).timebins;
    y = mean(PSTH_MUA(nprobe).R_HPC_ripples);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    % LCI = prctile(PSTH_MUA(nprobe).R_HPC_ripples_bootstrap,2.5);
    % UCI = prctile(PSTH_MUA(nprobe).R_HPC_ripples_bootstrap,97.5);
    LCI = y - std(PSTH_MUA(nprobe).R_HPC_ripples)./sqrt(no_events);
    UCI = y + std(PSTH_MUA(nprobe).R_HPC_ripples)./sqrt(no_events);
    PLOT = plot(x,y,'Color',colour_lines(3,:));hold on;
    ERROR_SHADE(4) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(3,:),'FaceAlpha','0.3','LineStyle','none');


    legend(ERROR_SHADE(3:4),{'Left HPC','Right HPC'},'Box','off')
    % xline(0,'r')
    title(probe_hemisphere_texts{nprobe})
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activty (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    % nexttile
end


%%




if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

