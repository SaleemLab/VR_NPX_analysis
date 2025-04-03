function plot_bilateral_synchronisation()

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

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));


load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'));
probability_SO_SO = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'));
probability_SO_SO_contralateral = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability.mat'));
probability_ripples_SO = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability.mat'));
probability_ripples_ripples = probability;

%% Plotting distribution of ripple during normalised duration of UP  (peaktimes)
probability = probability_normalised;
probability_baseline = probability_normalised_baseline;

time_wondows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;
probe_hemisphere_texts = {'Probability of ripples during left V1 normalised UP-DOWN duration','Probability of ripples during right V1 normalised UP-DOWN duration'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};

colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])




%% plotting Probability of SWR relative to UP and DOWN onset (peaktimes)
probability = probability_psth;
probability_baseline = probability_psth_baseline;

time_wondows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of ripples during left V1 UP-DOWN','Probability of ripples during right V1 UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
colour_lines = [215,25,28;44,123,182]/256;

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')
    xlim([-0.5 0.5])

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])


    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])
end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% Plotting distribution of ripple during normalised duration of UP  (WHOLE)
probability = probability_normalised_whole;
probability_baseline = probability_normalised_whole_baseline;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;
% probe_hemisphere_texts = {'Probability of ripples during left V1 normalised UP-DOWN duration','Probability of ripples during right V1 normalised UP-DOWN duration'};
probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};
colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])




%% plotting Probability of SWR relative to UP and DOWN onset (WHOLE)
probability = probability_psth_whole;
probability_baseline = probability_psth_whole_baseline;

time_wondows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
% probe_hemisphere_texts = {'Probability of ripples during left V1 UP-DOWN','Probability of ripples during right V1 UP-DOWN'};
probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};


for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')
    xlim([-0.5 0.5])

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])


    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])
%% plotting Probability of UP and/or DOWN relative to Ripple peaktime
probability = probability_ripples_SO;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left UP-DOWN during ripples','Probability of right UP-DOWN during ripples'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_UP_no = length(probability(nprobe).UP_all_index);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    all_ripple_no = probability(nprobe).L_ripple_no;

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).L_ripples_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_DOWN)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).L_ripples_DOWN_shuffled);
    LCI = prctile(probability(nprobe).L_ripples_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during left ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).L_ripples_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_UP)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).L_ripples_UP_shuffled);
    LCI = prctile(probability(nprobe).L_ripples_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during left ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    %%%% DOWN
    all_ripple_no = probability(nprobe).R_ripple_no;

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).R_ripples_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_DOWN)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).R_ripples_DOWN_shuffled);
    LCI = prctile(probability(nprobe).R_ripples_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during right ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).R_ripples_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_UP)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).R_ripples_UP_shuffled);
    LCI = prctile(probability(nprobe).R_ripples_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during right ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% plotting Probability of UP and/or DOWN relative to Ipsilateral UP DOWN
probability = probability_SO_SO;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left UP-DOWN relative to ipsilateral UP-DOWN','Probability of right UP-DOWN relative to ipsilateral UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% DOWN during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_UP_shuffled);
    LCI = prctile(probability(nprobe).DOWN_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% DOWN during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_DOWN_shuffled);
    LCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% UP during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_DOWN_shuffled);
    LCI = prctile(probability(nprobe).UP_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_UP_shuffled);
    LCI = prctile(probability(nprobe).UP_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during UP')
    xlabel('Time relative to UP (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% plotting Probability of UP and/or DOWN relative to contralateral UP DOWN
probability = probability_SO_SO_contralateral;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left UP-DOWN relative to contralateral UP-DOWN','Probability of right UP-DOWN relative to contralateral UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    if nprobe ==1
        mprobe = 2;
    else
        mprobe = 1;
    end

    all_UP_no = length(probability(mprobe).UP_all_index);
    all_DOWN_no = length(probability(mprobe).DOWN_all_index);

    %%%% DOWN during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_UP_shuffled);
    LCI = prctile(probability(nprobe).DOWN_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% DOWN during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_DOWN_shuffled);
    LCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% UP during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_DOWN_shuffled);
    LCI = prctile(probability(nprobe).UP_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_UP_shuffled);
    LCI = prctile(probability(nprobe).UP_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during UP')
    xlabel('Time relative to UP (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])



%% plotting Probability of ripple relative to Ripple peaktime
probability = probability_ripples_ripples;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left ripples during ripples','Probability of right ripples during ripples'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_ripple_no = length(probability(nprobe).L_ripple_no);

    %%%% Left ripples plotting
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).L_ripples_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).L_ripples_shuffled);
    LCI = prctile(probability(nprobe).L_ripples_shuffled,2.5);
    UCI = prctile(probability(nprobe).L_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of ripples during left ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    %%%% DOWN
    all_ripple_no = probability(nprobe).R_ripple_no;

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).R_ripples_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).R_ripples_shuffled);
    LCI = prctile(probability(nprobe).R_ripples_shuffled,2.5);
    UCI = prctile(probability(nprobe).R_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of ripples during right ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])
