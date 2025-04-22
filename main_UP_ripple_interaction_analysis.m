%% all events for sleep analysis (spiking but also event info)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar 
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]

%% Loading extracted SO ripples spindles info for plotting and analysis 
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
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
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));

%% Plotting basic temporal probability of UP DOWN and ripple
% plot_UP_DOWN_ripple_probability
plot_bilateral_synchronisation(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)

plot_ipsi_contra_UP_DOWN_ripple_spindle_coupling(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)

plot_ipsi_contra_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)
plot_ipsi_contra_UP_DOWN_spindle_probability(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)


% plot_UP_DOWN_ripple_probability_baseline


plot_ipsilateral_contralateral_UP_DOWN_spindles_coupling(slow_waves_all,ripples_all,spindles_all)
plot_ipsilateral_contralateral_UP_DOWN_spindles_coupling(slow_waves_all,ripples_all,spindles_all)

plot_ipsi_contra_MUA_probability_synchronisation(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,event_info,sessions_to_process)

% plot_UP_DOWN_spindle_probability
plot_UP_DOWN_spindle_probability_baseline
%% Plotting and calculating MUA relative to UP DOWN and ripple
all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

UP_DOWN_ripple_PSTH_MUA = calculate_UP_DOWN_ripple_PSTH...
    (slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,'option','MUA','time_option','absolute');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');

UP_DOWN_relative_PSTH_MUA = calculate_UP_DOWN_relative_PSTH...
    (slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,'option','MUA');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_relative_PSTH_MUA.mat'),'UP_DOWN_relative_PSTH_MUA');

plot_UP_DOWN_ripple_MUA_PSTH

%% The effect of ripple timing and power on DOWN UP transition and MUA response
%%%% P(ripples) during UP DOWN
clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
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
all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');


load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised.mat'));
probability_normalised = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));
probability_psth = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));

% plot_ripple_timing_UP(slow_waves_all,ripples_all,event_info,probability_psth,probability_normalised,probability_psth_whole,probability_normalised_whole,UP_DOWN_ripple_PSTH_MUA);
% plot_ripple_power_UP;


plot_ipsi_contra_ripple_power_UP


%% Analyse UP - DOWN transition
analyse_UP_DOWN_transition(event_info,ripples_all,spindles_all,slow_waves_all)

analyse_DOWN_UP_transition(event_info,ripples_all,spindles_all)


%% distribution of DOWN duration with high vs low ripple

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

hemisphere_text = {'L','R'}
time_index = 6:45;
time_index = reshape(time_index, 5, [])';
time_edges = -0.5:0.02:0.5;
time_windows = time_edges(1)+0.02/2:0.02:time_edges(end)-0.02/2;
% time_windows = reshape(time_windows, 5, [])';
nfigure = 1;
for nprobe = 1:2
    fig(nfigure)=figure;
    fig(nfigure).Name=sprintf('%s DOWN duration with and without L ripples',hemisphere_text{nprobe});
    fig(nfigure).Position = [680 150 1200 800]
    for n = 1:size(time_index,1)
        subplot(2,4,n)

        histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).L_ripples_DOWN(:,time_index(n,:)),2)>0)),0:0.01:0.4,'Normalization','probability');
        hold on;histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).L_ripples_DOWN(:,1:26),2)==0)),0:0.01:0.4,'Normalization','probability')
        ylim([0 0.1])
        title(sprintf('Ripple from %.2f to %.2f',time_windows(time_index(n,1)),time_windows(time_index(n,end))))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlabel('DOWN duration (s)');
        ylabel('probability');

        % han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
        % han.XLabel.Visible = 'on'; % Turn on visibility for x-label
        % han.YLabel.Visible = 'on'; % Turn on visibility for y-label
        % xlabel(han, 'Time (s)');
        % ylabel(han, 'T1/T2 bias');
    end
    legend('with L ripples','without L ripples','box','off')

    nfigure = nfigure+1;
    time_index = 6:45;
    time_index = reshape(time_index, 5, [])';
    time_edges = -0.5:0.02:0.5;
    time_windows = time_edges(1)+0.02/2:0.02:time_edges(end)-0.02/2;
    % time_windows = reshape(time_windows, 5, [])';
    fig(nfigure)=figure;
    fig(nfigure).Name=sprintf('%s DOWN duration with and without R ripples',hemisphere_text{nprobe});
    fig(nfigure).Position = [680 150 1200 800]
    for n = 1:size(time_index,1)
        subplot(2,4,n)
        histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).R_ripples_DOWN(:,time_index(n,:)),2)>0)),0:0.01:0.4,'Normalization','probability');
        hold on;histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).R_ripples_DOWN(:,1:26),2)==0)),0:0.01:0.4,'Normalization','probability')
        ylim([0 0.1])
        title(sprintf('Ripple from %.2f to %.2f',time_windows(time_index(n,1)),time_windows(time_index(n,end))))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlabel('DOWN duration (s)');
        ylabel('probability');

        % han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
        % han.XLabel.Visible = 'on'; % Turn on visibility for x-label
        % han.YLabel.Visible = 'on'; % Turn on visibility for y-label
        % xlabel(han, 'Time (s)');
        % ylabel(han, 'T1/T2 bias');
    end
    legend('with R ripples','without R ripples','box','off')
    nfigure = nfigure+1;
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% distribution of DOWN duration with and without ripples
hemisphere_text = {'L','R'}
time_index = 6:45;
time_index = reshape(time_index, 5, [])';
time_edges = -0.5:0.02:0.5;
time_windows = time_edges(1)+0.02/2:0.02:time_edges(end)-0.02/2;
% time_windows = reshape(time_windows, 5, [])';
nfigure = 1;
for nprobe = 1:2
    fig(nfigure)=figure;
    fig(nfigure).Name=sprintf('%s DOWN duration with and without L ripples',hemisphere_text{nprobe});
    fig(nfigure).Position = [680 150 1200 800]
    for n = 1:size(time_index,1)
        subplot(2,4,n)
        histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).L_ripples_DOWN(:,time_index(n,:)),2)>0)),0:0.01:0.4,'Normalization','probability');
        hold on;histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).L_ripples_DOWN(:,1:26),2)==0)),0:0.01:0.4,'Normalization','probability')
        ylim([0 0.1])
        title(sprintf('Ripple from %.2f to %.2f',time_windows(time_index(n,1)),time_windows(time_index(n,end))))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlabel('DOWN duration (s)');
        ylabel('probability');

        % han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
        % han.XLabel.Visible = 'on'; % Turn on visibility for x-label
        % han.YLabel.Visible = 'on'; % Turn on visibility for y-label
        % xlabel(han, 'Time (s)');
        % ylabel(han, 'T1/T2 bias');
    end
    legend('with L ripples','without L ripples','box','off')

    nfigure = nfigure+1;
    time_index = 6:45;
    time_index = reshape(time_index, 5, [])';
    time_edges = -0.5:0.02:0.5;
    time_windows = time_edges(1)+0.02/2:0.02:time_edges(end)-0.02/2;
    % time_windows = reshape(time_windows, 5, [])';
    fig(nfigure)=figure;
    fig(nfigure).Name=sprintf('%s DOWN duration with and without R ripples',hemisphere_text{nprobe});
    fig(nfigure).Position = [680 150 1200 800]
    for n = 1:size(time_index,1)
        subplot(2,4,n)
        histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).R_ripples_DOWN(:,time_index(n,:)),2)>0)),0:0.01:0.4,'Normalization','probability');
        hold on;histogram(probability(nprobe).DOWN_duration(find(sum(probability(nprobe).R_ripples_DOWN(:,1:26),2)==0)),0:0.01:0.4,'Normalization','probability')
        ylim([0 0.1])
        title(sprintf('Ripple from %.2f to %.2f',time_windows(time_index(n,1)),time_windows(time_index(n,end))))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlabel('DOWN duration (s)');
        ylabel('probability');

        % han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
        % han.XLabel.Visible = 'on'; % Turn on visibility for x-label
        % han.YLabel.Visible = 'on'; % Turn on visibility for y-label
        % xlabel(han, 'Time (s)');
        % ylabel(han, 'T1/T2 bias');
    end
    legend('with R ripples','without R ripples','box','off')
    nfigure = nfigure+1;
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])
%% distribution of UP duration with and without ripples
hemisphere_text = {'L','R'}
time_index = 6:45;
time_index = reshape(time_index, 5, [])';
time_edges = -0.5:0.02:0.5;
time_windows = time_edges(1)+0.02/2:0.02:time_edges(end)-0.02/2;
% time_windows = reshape(time_windows, 5, [])';
nfigure = 1;
for nprobe = 1:2
    fig(nfigure)=figure;
    fig(nfigure).Name=sprintf('%s UP duration with and without L ripples',hemisphere_text{nprobe});
    fig(nfigure).Position = [680 150 1200 800]
    for n = 1:size(time_index,1)
        subplot(2,4,n)
        histogram(probability(nprobe).UP_duration(find(sum(probability(nprobe).L_ripples_UP(:,time_index(n,:)),2)>0)),0:0.01:0.4,'Normalization','probability');
        hold on;histogram(probability(nprobe).UP_duration(find(sum(probability(nprobe).L_ripples_UP(:,1:26),2)==0)),0:0.01:0.4,'Normalization','probability')
        ylim([0 0.1])
        title(sprintf('Ripple from %.2f to %.2f',time_windows(time_index(n,1)),time_windows(time_index(n,end))))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlabel('DOWN duration (s)');
        ylabel('probability');

        % han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
        % han.XLabel.Visible = 'on'; % Turn on visibility for x-label
        % han.YLabel.Visible = 'on'; % Turn on visibility for y-label
        % xlabel(han, 'Time (s)');
        % ylabel(han, 'T1/T2 bias');
    end
    legend('with L ripples','without L ripples','box','off')

    nfigure = nfigure+1;
    time_index = 6:45;
    time_index = reshape(time_index, 5, [])';
    time_edges = -0.5:0.02:0.5;
    time_windows = time_edges(1)+0.02/2:0.02:time_edges(end)-0.02/2;
    % time_windows = reshape(time_windows, 5, [])';
    fig(nfigure)=figure;
    fig(nfigure).Name=sprintf('%s UP duration with and without R ripples',hemisphere_text{nprobe});
    fig(nfigure).Position = [680 150 1200 800]
    for n = 1:size(time_index,1)
        subplot(2,4,n)
        histogram(probability(nprobe).UP_duration(find(sum(probability(nprobe).R_ripples_UP(:,time_index(n,:)),2)>0)),0:0.01:0.4,'Normalization','probability');
        hold on;histogram(probability(nprobe).UP_duration(find(sum(probability(nprobe).R_ripples_UP(:,1:26),2)==0)),0:0.01:0.4,'Normalization','probability')
        ylim([0 0.1])
        title(sprintf('Ripple from %.2f to %.2f',time_windows(time_index(n,1)),time_windows(time_index(n,end))))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlabel('DOWN duration (s)');
        ylabel('probability');

        % han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
        % han.XLabel.Visible = 'on'; % Turn on visibility for x-label
        % han.YLabel.Visible = 'on'; % Turn on visibility for y-label
        % xlabel(han, 'Time (s)');
        % ylabel(han, 'T1/T2 bias');
    end
    legend('with R ripples','without R ripples','box','off')
    nfigure = nfigure+1;
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])


%%
clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

%%%%%%%%%% SWR and UP DOWN relationship
time_wondows = [-0.5 0.5];
time_bin = 0.01;
num_bins=50; % divide one UP event into 20 bins

probability=[];
all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;


for nprobe = 1:2
    %     fig(nprobe)=figure;
    %     fig(nprobe).Position = [982 50 700 950];
    probability(nprobe).L_ripples_DOWN_distribution = [];
    probability(nprobe).L_ripples_DOWN_distribution_shuffled = [];
    probability(nprobe).L_ripples_UP_distribution = [];
    probability(nprobe).L_ripples_UP_distribution_shuffled = [];
    % DOWN state
    for nsession = 1:10
        UP_index = find(slow_waves_all(nprobe).UP_session_count == nsession);
        DOWN_index = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        ripples_index = find(ripples_all(1).session_count == nsession& ripples_all(1).SWS_index == 1);
        UP_duration = slow_waves_all(nprobe).UP_ints(UP_index,2)-slow_waves_all(nprobe).UP_ints(UP_index,1);
        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index(UP_duration<10),:);
        % UP_duration<2

        %         [~,event_index,relative_duration] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripples_all(1).peaktimes(ripples_index),num_bins,0);
        [~,event_index,relative_duration] = calculate_relative_event_probability(UP_ints,ripples_all(1).peaktimes(ripples_index),num_bins,0);
        temp = event_index(~isnan(event_index));
        event_index(~isnan(event_index))= UP_index(temp); % Convert to event id across all sessions.
        probability(nprobe).L_ripples_UP_distribution = [probability(nprobe).L_ripples_UP_distribution [event_index; relative_duration]];

        [~,event_index,relative_duration] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripples_all(1).peaktimes(ripples_index),num_bins,0);
        temp = event_index(~isnan(event_index));
        event_index(~isnan(event_index))= DOWN_index(temp); % Convert to event id across all sessions.
        probability(nprobe).L_ripples_DOWN_distribution = [probability(nprobe).L_ripples_DOWN_distribution [event_index; relative_duration]];

        [~,event_index,relative_duration] = calculate_relative_event_probability(UP_ints,ripples_all(1).peaktimes(ripples_index),num_bins,1);
        temp = event_index(~isnan(event_index));
        event_index(~isnan(event_index))= UP_index(temp); % Convert to event id across all sessions.
        probability(nprobe).L_ripples_UP_distribution_shuffled = [probability(nprobe).L_ripples_UP_distribution_shuffled [event_index; relative_duration]];

        [~,event_index,relative_duration] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripples_all(1).peaktimes(ripples_index),num_bins,1);
        temp = event_index(~isnan(event_index));
        event_index(~isnan(event_index))= DOWN_index(temp); % Convert to event id across all sessions.
        probability(nprobe).L_ripples_DOWN_distribution_shuffled = [probability(nprobe).L_ripples_DOWN_distribution_shuffled [event_index; relative_duration]];


        % probability(nprobe).L_ripples_UP(nsession,:) = calculate_relative_event_probability(UP_ints,ripples_all(1).peaktimes(ripples_index),num_bins,0);
        % probability(nprobe).L_ripples_DOWN_shuffled(nsession,:) = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripples_all(1).peaktimes(ripples_index),num_bins,1);
        % probability(nprobe).L_ripples_UP_shuffled(nsession,:) = calculate_relative_event_probability(UP_ints,ripples_all(1).peaktimes(ripples_index),num_bins,1);
    end
end


for nprobe = 1:2
    figure
    for nsession = 1:10
        nexttile
        histogram(log10(slow_waves_all(nprobe).UP_ints(slow_waves_all(nprobe).UP_session_count ==nsession,2)-slow_waves_all(nprobe).UP_ints(slow_waves_all(nprobe).UP_session_count ==nsession,1)),100,'EdgeColor','none');
        hold on;
        histogram(log10(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).UP_session_count ==nsession,2)-slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).UP_session_count ==nsession,1)),100,'EdgeColor','none');

        xlabel('Event duration (log10 sec)')
        ylabel('event counts')

        title(sprintf('Session %i',nsession))
    end
    legend('UP','DOWN','Box','off')
end

probe_hemisphere_texts = {'left','right'};

for nprobe = 1:2
%     ripples_peaktimes = ripples_all(1).peaktimes(ripples_all(1).SWS_index);
    ripples_zscore = ripples_all(1).peak_zscore(ripples_all(1).SWS_index)';
    ripples_duration = ripples_all(1).offset(ripples_all(1).SWS_index)'-ripples_all(1).onset(ripples_all(1).SWS_index)';
    [N,Xedges,Yedges] = histcounts2(ripples_zscore,probability(nprobe).L_ripples_UP_distribution(2,:),linspace(5,prctile(ripples_zscore,99),30),0:0.05:1);
   
    DOWN_duration = slow_waves_all(nprobe).DOWN_ints(:,2)-slow_waves_all(nprobe).DOWN_ints(:,1);
    UP_duration = slow_waves_all(nprobe).UP_ints(:,2)-slow_waves_all(nprobe).UP_ints(:,1);
    
    event_id = 1:length(slow_waves_all(nprobe).UP_session_count);
    % event_id = event_id(UP_duration<10);

    ripples_UP_index = event_id(ismember(event_id,[probability(nprobe).L_ripples_UP_distribution(1,:)]));
    ripples_UP_index(end) = [];
    UP_only_index = event_id(~ismember(event_id,[probability(nprobe).L_ripples_UP_distribution(1,:)]));
    UP_only_index(end) = [];


    ripples_DOWN_index = event_id(ismember(event_id,[probability(nprobe).L_ripples_DOWN_distribution(1,:)]));
    ripples_DOWN_index(end) = [];
    DOWN_only_index = event_id(~ismember(event_id,[probability(nprobe).L_ripples_DOWN_distribution(1,:)]));
    DOWN_only_index(end) = [];

    event_session_count = slow_waves_all(nprobe).DOWN_session_count(ripples_UP_index)';
    next_DOWN_duration = slow_waves_all(nprobe).DOWN_ints(ripples_UP_index,2)-slow_waves_all(nprobe).DOWN_ints(ripples_UP_index,1);
    next_DOWN_peaks = slow_waves_all(nprobe).SWpeakmag(ripples_UP_index);

    last_ripple_zscore_UP = [];
    peak_ripple_zscore_UP=[];
    ripples_number_UP = [];
    cum_ripples_duration_UP = [];
    UP_end_ripples =[];
    last_ripple_UP=[];
    first_ripple_UP=[];
    mid_ripple_UP=[];
    PSD_slope_UP=[];
    first_ripple_zscore_UP=[];
    median_ripple_zscore_UP=[];

    for nevent = 1:length(ripples_UP_index)
        relative_times = probability(nprobe).L_ripples_UP_distribution(2,(probability(nprobe).L_ripples_UP_distribution(1,:) == ripples_UP_index(nevent)));
        ripples_zscore_this_event = ripples_zscore(probability(nprobe).L_ripples_UP_distribution(1,:) == ripples_UP_index(nevent));
        last_ripple_UP(nevent) = relative_times(end);
        first_ripple_UP(nevent) = relative_times(1);
        mid_ripple_UP(nevent) = median(relative_times);
        PSD_slope_UP(nevent)  = slow_waves_all(nprobe).DOWN_PSD_slope(ripples_UP_index(nevent));

        last_ripple_zscore_UP(nevent) = ripples_zscore_this_event(end);
        peak_ripple_zscore_UP(nevent) = max(ripples_zscore_this_event);
        first_ripple_zscore_UP(nevent) = ripples_zscore_this_event(1);
        median_ripple_zscore_UP(nevent) = median(ripples_zscore_this_event);
        ripples_number_UP(nevent) = sum(probability(nprobe).L_ripples_UP_distribution(1,:) == ripples_UP_index(nevent));
        cum_ripples_duration_UP(nevent) = sum(ripples_duration(probability(nprobe).L_ripples_UP_distribution(1,:) == ripples_UP_index(nevent)));
        if relative_times(end)>0.80
            UP_end_ripples(nevent) = 1;
        else 
             UP_end_ripples(nevent) = 0;
        end
    end

    last_ripple_zscore_DOWN = [];
    peak_ripple_zscore_DOWN=[];
    ripples_number_DOWN = [];
    cum_ripples_duration_DOWN = [];
    DOWN_end_ripples =[];
    last_ripple_DOWN=[];
    first_ripple_DOWN=[];
    mid_ripple_DOWN=[];
    PSD_slope_DOWN=[];
    first_ripple_zscore_DOWN=[];
    median_ripple_zscore_DOWN=[];

    for nevent = 1:length(ripples_DOWN_index)
        relative_times = probability(nprobe).L_ripples_DOWN_distribution(2,(probability(nprobe).L_ripples_DOWN_distribution(1,:) == ripples_DOWN_index(nevent)));
        ripples_zscore_this_event = ripples_zscore(probability(nprobe).L_ripples_DOWN_distribution(1,:) == ripples_DOWN_index(nevent));
        last_ripple_DOWN(nevent) = relative_times(end);
        first_ripple_DOWN(nevent) = relative_times(1);
        mid_ripple_DOWN(nevent) = median(relative_times);
        PSD_slope_DOWN(nevent)  = slow_waves_all(nprobe).DOWN_PSD_slope(ripples_DOWN_index(nevent));

        last_ripple_zscore_DOWN(nevent) = ripples_zscore_this_event(end);
        peak_ripple_zscore_DOWN(nevent) = max(ripples_zscore_this_event);
        first_ripple_zscore_DOWN(nevent) = ripples_zscore_this_event(1);
        median_ripple_zscore_DOWN(nevent) = median(ripples_zscore_this_event);
        ripples_number_DOWN(nevent) = sum(probability(nprobe).L_ripples_DOWN_distribution(1,:) == ripples_DOWN_index(nevent));
        cum_ripples_duration_DOWN(nevent) = sum(ripples_duration(probability(nprobe).L_ripples_DOWN_distribution(1,:) == ripples_DOWN_index(nevent)));
        if relative_times(end)>0.80
            DOWN_end_ripples(nevent) = 1;
        else 
             DOWN_end_ripples(nevent) = 0;
        end
    end

    scatter(ripples_number_DOWN,peak_ripple_zscore_DOWN,'k','filled','MarkerFaceAlpha',.01,'MarkerEdgeAlpha',0);
    xlabel('Ripple number per DOWN')
    ylabel('Peak ripple power dring DOWN')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    
    for nsession = 1:10
        nexttile
        hold on

        X = log10(UP_duration(UP_only_index(ismember(UP_only_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        Y = log10(DOWN_duration(UP_only_index(ismember(UP_only_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        %         X = log10(UP_duration(ripples_UP_index(ismember(ripples_UP_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        %         Y = log10(DOWN_duration(ripples_UP_index(ismember(ripples_UP_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
        mdl = fitlm(X',Y');
        [pval,F_stat,~] = coefTest(mdl);
        R2 = mdl.Rsquared.Adjusted;
        x =[min(X) max(X)];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);
        PLOT(1)= plot(x,y_est,':','Color','k','LineWidth',2)
        %     end
        text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');


        hold on

        %         X = log10(UP_duration(UP_only_index(ismember(UP_only_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        %         Y = log10(DOWN_duration(UP_only_index(ismember(UP_only_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        X = log10(UP_duration(ripples_UP_index(ismember(ripples_UP_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        Y = log10(DOWN_duration(ripples_UP_index(ismember(ripples_UP_index,find(slow_waves_all(nprobe).UP_session_count == nsession)))));
        scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
        mdl = fitlm(X',Y');
        [pval,F_stat,~] = coefTest(mdl);
        R2 = mdl.Rsquared.Adjusted;
        x =[min(X) max(X)];
        b = mdl.Coefficients.Estimate';
        y_est = polyval(fliplr(b),x);
        PLOT(2)= plot(x,y_est,':','Color','r','LineWidth',2)
        text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');
        %     end
        xlabel('Current UP duration (log10 sec)')
        ylabel('Next Down duration (log10 sec)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        legend([PLOT(1) PLOT(2)],{'UP only','UP with ripples'},'Box','off','Location','best')
        title(sprintf('session %i',nsession))
    end

    %%%%%%%%%%%%%%%%% UP and DOWN relationship
    %
    counter = counter + 1;
    fig(counter)=figure;
    fig(counter).Position = [982 50 700 950];
    fig(counter).Name = sprintf('UP and DOWN relationship %s V1',probe_hemisphere_texts{nprobe})
    next_DOWN_duration_no_ripples = slow_waves_all(nprobe).DOWN_ints(UP_only_index,2)-slow_waves_all(nprobe).DOWN_ints(UP_only_index,1);

    nexttile
    hold on
    X = log10(UP_duration(UP_only_index));
    Y = log10(DOWN_duration(UP_only_index));
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','k','LineWidth',2)
    %     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');

%     [N,xbin,bin]=histcounts(X,100);
% 
%     plot(max(y_est)+histcounts(X,100)/max(histcounts(X,100)))

    
    X = log10(UP_duration(ripples_UP_index));
    Y = log10(DOWN_duration(ripples_UP_index));
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2)= plot(x,y_est,':','Color','r','LineWidth',2)
    %     end
    text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');
    xlabel('Current UP duration (log10 sec)')
    ylabel('Next Down duration (log10 sec)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    legend([PLOT(1) PLOT(2)],{'UP only','UP with ripples'},'Box','off','Location','best')



    nexttile
    hold on
    X = log10(UP_duration(UP_only_index+1));
    Y = log10(DOWN_duration(UP_only_index));
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','k','LineWidth',2)
    %     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');

%     [N,xbin,bin]=histcounts(X,100);
% 
%     plot(max(y_est)+histcounts(X,100)/max(histcounts(X,100)))


    X = log10(UP_duration(ripples_UP_index+1));
    Y = log10(DOWN_duration(ripples_UP_index));
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2)= plot(x,y_est,':','Color','r','LineWidth',2)
    %     end
    text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');
    xlabel('Next UP duration (log10 sec)')
    ylabel('Current Down duration (log10 sec)')
    legend([PLOT(1) PLOT(2)],{'UP only','UP with ripples'},'Box','off','Location','best')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%%%%%%% UP without ripple is shorter than UP with ripples
    UP_only_duration_dist=[];
    UP_ripples_duration_dist=[];
    for nsession = 1:all_sessions
        session_index =intersect(UP_only_index,find(slow_waves_all(nprobe).UP_session_count==nsession));
        [N,xedge,~]=histcounts(log10(UP_duration(session_index)),-2:0.01:1);
        UP_only_duration_dist(nsession,:) = N/sum(N);

        session_index =intersect(ripples_UP_index,find(slow_waves_all(nprobe).UP_session_count==nsession));
        [N,xedge,~]=histcounts(log10(UP_duration(session_index)),-2:0.01:1);
        UP_ripples_duration_dist(nsession,:) = N/sum(N);
    end

    nexttile
    x = xedge(1)+mean(diff(xedge))/2:mean(diff(xedge)):xedge(end)-mean(diff(xedge))/2;
    y = mean(cumsum(UP_only_duration_dist,2));
    SE = std(cumsum(UP_only_duration_dist,2))./sqrt(all_sessions);
    PLOT(1) = plot(x,y,'Color','k');hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],'k','FaceAlpha','0.3','LineStyle','none');

    x = xedge(1)+mean(diff(xedge))/2:mean(diff(xedge)):xedge(end)-mean(diff(xedge))/2;
    y = mean(cumsum(UP_ripples_duration_dist,2));
    SE = std(cumsum(UP_ripples_duration_dist,2))./sqrt(all_sessions);
    PLOT(2) = plot(x,y,'Color','r');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],'r','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'})
    % xline(0,'r')
    % title('Cumulative prob of reft ripples during UP')
    xlabel('UP duration (log10 second)')
    ylabel('Cumulative probability')
    legend([PLOT(1) PLOT(2)],{'UP only','UP with ripples'},'Box','off','Location','best')
    ylim([0 1])
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    [p,h,stats] = ranksum(UP_duration(UP_only_index),UP_duration(ripples_UP_index),'tail','left');
    if p < 0.05
        text(gca,.8,0.6,'*','Units','Normalized','FontName','Arial','Color','r');
    else
        text(gca,.8,0.6,'*','Units','Normalized','FontName','Arial','Color','k');
    end
    title('UP duration with and without ripples')
  
    
    DOWN_only_duration_dist=[];
    DOWN_ripples_duration_dist=[];
    for nsession = 1:all_sessions
        session_index =intersect(UP_only_index,find(slow_waves_all(nprobe).UP_session_count==nsession));
        [N,xedge,~]=histcounts(log10(DOWN_duration(session_index)),-2:0.01:1);
        DOWN_only_duration_dist(nsession,:) = N/sum(N);

        session_index =intersect(ripples_UP_index,find(slow_waves_all(nprobe).UP_session_count==nsession));
        [N,xedge,~]=histcounts(log10(DOWN_duration(session_index)),-2:0.01:1);
        DOWN_ripples_duration_dist(nsession,:) = N/sum(N);
    end

    nexttile
    x = xedge(1)+mean(diff(xedge))/2:mean(diff(xedge)):xedge(end)-mean(diff(xedge))/2;
    y = mean(cumsum(DOWN_only_duration_dist,2));
    SE = std(cumsum(DOWN_only_duration_dist,2))./sqrt(all_sessions);
    PLOT(1) = plot(x,y,'Color','k');hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],'k','FaceAlpha','0.3','LineStyle','none');

    x = xedge(1)+mean(diff(xedge))/2:mean(diff(xedge)):xedge(end)-mean(diff(xedge))/2;
    y = mean(cumsum(DOWN_ripples_duration_dist,2));
    SE = std(cumsum(DOWN_ripples_duration_dist,2))./sqrt(all_sessions);
    PLOT(2) = plot(x,y,'Color','r');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],'r','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'})
    % xline(0,'r')
    % title('Cumulative prob of reft ripples during UP')
    xlabel('DOWN duration (log10 second)')
    ylabel('Cumulative probability')
    legend([PLOT(1) PLOT(2)],{'UP only','UP with ripples'},'Box','off','Location','best')
    ylim([0 1])
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    [p,h,stats] = ranksum(DOWN_duration(UP_only_index),DOWN_duration(ripples_UP_index),'tail','left');
    if p < 0.05
        text(gca,.8,0.6,'*','Units','Normalized','FontName','Arial','Color','r');
    else
        text(gca,.8,0.6,'*','Units','Normalized','FontName','Arial','Color','k');
    end
    title('Next DOWN duration with and without ripples')

    DOWN_only_duration_dist=[];
    DOWN_ripples_duration_dist=[];
    for nsession = 1:all_sessions
        session_index =intersect(DOWN_only_index,find(slow_waves_all(nprobe).UP_session_count==nsession));
        [N,xedge,~]=histcounts(log10(DOWN_duration(session_index)),-2:0.01:1);
        DOWN_only_duration_dist(nsession,:) = N/sum(N);

        session_index =intersect(ripples_DOWN_index,find(slow_waves_all(nprobe).UP_session_count==nsession));
        [N,xedge,~]=histcounts(log10(DOWN_duration(session_index)),-2:0.01:1);
        DOWN_ripples_duration_dist(nsession,:) = N/sum(N);
    end

    nexttile
    x = xedge(1)+mean(diff(xedge))/2:mean(diff(xedge)):xedge(end)-mean(diff(xedge))/2;
    y = mean(cumsum(DOWN_only_duration_dist,2));
    SE = std(cumsum(DOWN_only_duration_dist,2))./sqrt(all_sessions);
    PLOT(1) = plot(x,y,'Color','k');hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],'k','FaceAlpha','0.3','LineStyle','none');

    x = xedge(1)+mean(diff(xedge))/2:mean(diff(xedge)):xedge(end)-mean(diff(xedge))/2;
    y = mean(cumsum(DOWN_ripples_duration_dist,2));
    SE = std(cumsum(DOWN_ripples_duration_dist,2))./sqrt(all_sessions);
    PLOT(2) = plot(x,y,'Color','r');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],'r','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'})
    % xline(0,'r')
    % title('Cumulative prob of reft ripples during UP')
    xlabel('DOWN duration (log10 second)')
    ylabel('Cumulative probability')
    legend([PLOT(1) PLOT(2)],{'UP only','UP with ripples'},'Box','off','Location','best')
    ylim([0 1])
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    [p,h,stats] = ranksum(DOWN_duration(DOWN_only_index),DOWN_duration(ripples_DOWN_index),'tail','left');
    if p < 0.05
        text(gca,.8,0.6,'*','Units','Normalized','FontName','Arial','Color','r');
    else
        text(gca,.8,0.6,'*','Units','Normalized','FontName','Arial','Color','k');
    end
    title('Current DOWN duration with and without ripples')

    %%%%%%%%%%%%%%%% DOWN duration 

%     X = log10(cum_ripples_duration_UP);
%     Y = log10(UP_duration(ripples_UP_index));
%     scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
% 
%     X = last_ripple_zscore_UP;
%     Y = log10(UP_duration(ripples_UP_index));
%     scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
% 
% 

    X = peak_ripple_zscore_UP;
    Y = log10(UP_duration(ripples_UP_index));
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    hold on
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','r','LineWidth',2)
%     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');


    hold on;
    X = peak_ripple_zscore_UP;
    Y = log10(DOWN_duration(ripples_UP_index));
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);

    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2) =plot(x,y_est,':','Color','k','LineWidth',2)

    text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');
    xlabel('Ripple power per event')
    ylabel('Event duration (log10 second)')
    legend([PLOT(1) PLOT(2)],{'Current UP','Next DOWN'},'Box','off')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title('Peak Ripple power during UP vs Event duration')

    %%%%%%%%%%%%%%%%% 

    %%%%%%%%%%%%%%%%% cumulative ripple duration during UP vs DOWN and ripple zscore
    counter = counter + 1;
    fig(counter)=figure;
    fig(counter).Position = [522 46 1160 950];
    fig(counter).Name = sprintf('cumulative  %s ripple duration vs ripple zscore during UP and DOWN',probe_hemisphere_texts{nprobe})

    nexttile
    X = cum_ripples_duration_UP;
    Y = peak_ripple_zscore_UP;
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    hold on
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','r','LineWidth',2)
%     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');

    hold on;
    X = cum_ripples_duration_DOWN;
    Y = peak_ripple_zscore_DOWN;
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);

    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2) =plot(x,y_est,':','Color','k','LineWidth',2)

    text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');
    xlabel('Cumulative ripple duration per event')
    ylabel('Peak ripple power per event')
    legend([PLOT(1) PLOT(2)],{'UP','DOWN'},'Box','off')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title('Peak ripple power')

    nexttile
    X = cum_ripples_duration_UP;
    Y = median_ripple_zscore_UP;
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    hold on
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','r','LineWidth',2)
%     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');

    hold on;
    X = cum_ripples_duration_DOWN;
    Y = median_ripple_zscore_DOWN;
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);

    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2) =plot(x,y_est,':','Color','k','LineWidth',2)

    text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');
    xlabel('Cumulative ripple duration per event')
    ylabel('Peak ripple power per event')
    legend([PLOT(1) PLOT(2)],{'UP','DOWN'},'Box','off')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title('Median ripple power')

    nexttile
    X = cum_ripples_duration_UP;
    Y = last_ripple_zscore_UP;
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    hold on
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','r','LineWidth',2)
%     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');

    hold on;
    X = cum_ripples_duration_DOWN;
    Y = last_ripple_zscore_DOWN;
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);

    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2) =plot(x,y_est,':','Color','k','LineWidth',2)

    text(gca,.7,0.2,sprintf('p = %.2d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');
    xlabel('Cumulative ripple duration per event')
    ylabel('Peak ripple power per event')
    legend([PLOT(1) PLOT(2)],{'UP','DOWN'},'Box','off')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title('Last ripple power')

    nexttile
    X = cum_ripples_duration_UP;
    Y = first_ripple_zscore_UP;
    scatter(X,Y,'r','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
    hold on
    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(1)= plot(x,y_est,':','Color','r','LineWidth',2)
%     end
    text(gca,.7,0.1,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','r');

    hold on;
    X = cum_ripples_duration_DOWN;
    Y = first_ripple_zscore_DOWN;
    scatter(X,Y,'k','filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);

    mdl = fitlm(X',Y');
    [pval,F_stat,~] = coefTest(mdl);
    R2 = mdl.Rsquared.Adjusted;
    x =[min(X) max(X)];
    b = mdl.Coefficients.Estimate';
    y_est = polyval(fliplr(b),x);
    PLOT(2) =plot(x,y_est,':','Color','k','LineWidth',2)

    text(gca,.7,0.2,sprintf('p = %.3d \nR2 = %.3d',pval,R2),'Units','Normalized','FontName','Arial','Color','k');
    xlabel('Cumulative ripple duration per event')
    ylabel('Peak ripple power per event')
    legend([PLOT(1) PLOT(2)],{'UP','DOWN'},'Box','off')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title('First ripple power')




    % ripple duration and ripple zscore
    event_index = find(isnan(probability(nprobe).L_ripples_UP_distribution(2,:)));
    scatter(ripples_duration,ripples_zscore,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    scatter(ripples_duration(event_index),ripples_zscore(event_index),'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);



    event_index = find(~isnan(probability(nprobe).L_ripples_UP_distribution(2,:)));

    scatter(ripples_duration(intersect(event_index,find(ripples_duration>0.05))),ripples_zscore(intersect(event_index,find(ripples_duration>0.05))),'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
%     scatter(cum_ripples_duration_UP,last_ripple_zscore_UP,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);



     [N,Xedges,Yedges] = histcounts2(peak_ripple_zscore_UP',next_DOWN_duration,0:1:20,0:0.05:1);
%      imagesc(N')
     imagesc(N'./sum(N'))
     xticks(1:20+0.5)
     xticklabels(Xedges)
     yticks(1:length(Yedges)+0.5)
     yticklabels(Yedges)

     %     [~,index]=max(N);
     colorbar
     colormap(flipud(gray))


    scatter(last_ripple_zscore_UP(UP_end_ripples==1),next_DOWN_duration(UP_end_ripples==1),'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Last ripple ripple power dring UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    scatter(next_DOWN_duration,next_DOWN_peaks,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Last ripple ripple power dring UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    scatter(last_ripple_zscore_UP,next_DOWN_duration,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Last ripple ripple power dring UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(probability(nprobe).L_ripples_UP_distribution(2,:),ripples_zscore,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    xlabel('Ripple times (normalised UP duration)')
    ylabel('Ripples zscore')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(last_ripple_zscore_UP,next_DOWN_duration,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Last ripple ripple power dring UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(peak_ripple_zscore_UP,next_DOWN_duration,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Peak ripple ripple power dring UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(last_ripple_zscore_UP,next_DOWN_peaks,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Last ripple ripple power dring UP')
    ylabel('Next DOWN peak amplitude')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(peak_ripple_zscore_UP,next_DOWN_peaks,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Peak ripple ripple power dring UP')
    ylabel('Next DOWN peak amplitude')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(ripples_number_UP,next_DOWN_duration,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Ripple number per UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(ripples_number_UP,next_DOWN_peaks,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Ripple number per UP')
    ylabel('Next DOWN peak amplitude')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    scatter(cum_ripples_duration_UP,next_DOWN_duration,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Cumulative ripple durations per UP')
    ylabel('Next DOWN duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    
    nexttile
    scatter(cum_ripples_duration_UP,next_DOWN_peaks,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Cumulative ripple durations per UP')
    ylabel('Next DOWN peak amplitude')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    scatter(cum_ripples_duration_UP,UP_duration(ripples_UP_index),'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlabel('Cumulative ripple durations per UP')
    ylabel('Current UP duration')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    scatter(UP_duration,slow_waves_all(nprobe).DOWN_PSD_slope,'k','filled','MarkerFaceAlpha',.02,'MarkerEdgeAlpha',0);
    xlim([0 2])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    next_DOWN_duration = slow_waves_all(nprobe).DOWN_ints(UP_only_index,2)-slow_waves_all(nprobe).DOWN_ints(UP_only_index,1);
    scatter(UP_duration(UP_only_index),next_DOWN_duration,'k','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);
    xlim([0 0.2])
    ylim([0 1])

%     beeswarm(4*UP_end_ripples',last_ripple_zscore_UP')

    for nBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',nBoot); % Set random seed for resampling
        resampled_data = datasample(s,last_ripple_zscore_UP(UP_end_ripples==0),sum(UP_end_ripples==1));
        A(nBoot) = mean(resampled_data);
        resampled_data = datasample(s,last_ripple_zscore_UP(UP_end_ripples==1),sum(UP_end_ripples==1));
        B(nBoot) = mean(resampled_data);

        p = ranksum(resampled_data,last_ripple_zscore_UP(UP_end_ripples==1));
        p_value(nBoot) = p;
    end
    prctile(A,2.5)
    prctile(A,97.5)
    prctile(B,2.5)
    prctile(B,97.5)

    average_duration=[];
    for n = unique(ripples_number_UP)
        for nsession = unique(event_session_count)
            average_duration(n,nsession) = mean(next_DOWN_duration(event_session_count == nsession & ripples_number_UP==n));
        end
    end
    plot(average_duration)

    % 
    next_DOWN_duration = slow_waves_all(nprobe).DOWN_ints(UP_only_index+1,2)-slow_waves_all(nprobe).DOWN_ints(UP_only_index+1,1);
    scatter(UP_duration(UP_only_index),next_DOWN_duration,'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
    xlim([0 5])
    ylim([0 1])
    ripples_zscore


     [N,Xedges,Yedges] = histcounts2(ripples_zscore,probability(nprobe).L_ripples_UP_distribution(2,:),linspace(5,prctile(ripples_zscore,99),20),0:0.05:1);
    imagesc(N./sum(N))

    [~,index]=max(N);
    colorbar
    colormap(flipud(gray))
end






%%%%%%%%%% V1 and HPC MUA z-score relative to all ripples, UP - DOWN ripples,
%%%%%%%%%% DOWN - UP ripples, UP ripples and DOWN ripples 

for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
    %%% Get spike counts
    last_spike = slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}(end);
    tvec = 0:0.01:last_spike;
    % mobility_interp = interp1(selected_clusters.tvec{1},double(selected_clusters.mobility_thresholded{1}),tvec,'previous');
    tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];
    if ~isempty(behavioural_state_merged_all.SWS{nsession})
        sleep_tvec = Restrict(tvec,behavioural_state_merged_all.SWS{nsession});
        sleep_index = ismember(tvec,sleep_tvec);
    end
    HPC_spike_counts=[];
    V1_spike_counts=[];
    w = gausswin(0.03*1/mean(diff(tvec)));
    w = w / sum(w);

    HPC_MUA_sleep=[];
    V1_MUA_sleep=[];
    % speed= filtfilt(w,1,zscore(histcounts(spike_times_sleep,tvec_edges))')';
    for nprobe = 1:length(slow_waves_all) % Get distribution of sleep spike counts
        temp_spike_count =[];
        spike_times = slow_waves_all(nprobe).V1_MUA_spiketimes{nsession};
        % temp_spike_count = filtfilt(w,1,histcounts(spike_times,tvec_edges)')';
        temp_spike_count = histcounts(spike_times,tvec_edges);
        V1_MUA_sleep{nprobe} = temp_spike_count(sleep_index==1);

        spike_times = slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession};
        % temp_spike_count = filtfilt(w,1,histcounts(spike_times,tvec_edges)')';
        temp_spike_count = histcounts(spike_times,tvec_edges);
        HPC_MUA_sleep{nprobe} = temp_spike_count(sleep_index==1);
    end



    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}, event_times, [-0.5 0.5], 0.01);
    binnedArray = filtfilt(w,1,binnedArray')';
    zscored_psth= reshape(binnedArray,1,[]);
    zscored_psth = (zscored_psth-mean(V1_MUA_sleep{nprobe}))/(std(V1_MUA_sleep{nprobe}));
    zscored_psth = reshape(zscored_psth,length(event_times),[]);

end


nprobe = 2;
mprobe = 1;


[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}, event_times, [-0.5 0.5], 0.01);
binnedArray = filtfilt(w,1,binnedArray')';

zscored_psth= reshape(binnedArray,1,[]);
zscored_psth = (zscored_psth-mean(V1_MUA_sleep{nprobe}))/(std(V1_MUA_sleep{nprobe}));
zscored_psth = reshape(zscored_psth,length(event_times),[]);
imagesc(zscored_psth)
colorbar
colormap(flipud(gray))

plot_perievent_spiketime_histogram



%%%%%%%%%%%%%%%% event veiwer
event_times1 = slow_waves(1).ints.UP;
event_times2 = slow_waves(2).ints.UP;

[tvec_epoches,idx] = RestrictInts(tvec',behavioural_state_merged.SWS); %Replace with InInterval
plot(tvec,HPC_spike_counts{3}');hold on
spacer = 1+max(HPC_spike_counts{3}(idx));
plot(tvec,V1_spike_counts{1}'+ spacer,'r');
spacer= 1+max(HPC_spike_counts{3}(idx))+1+max(V1_spike_counts{1}(idx));
plot(tvec,V1_spike_counts{2}' + spacer,'b');
for n = 1:30
    %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
    %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
    %             xline(min(tvec(tvec>event_times(n,2))),'b')
    %             xline(max(tvec(tvec<event_times(n,1))),'r')
    rectangle('Position',[min(tvec(tvec>event_times1(n,1))) 0.5...
        min(tvec(tvec>event_times1(n,2)))-max(tvec(tvec<event_times1(n,1))),...
        20],...
        'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.2)
end

for n = 1:30
    %             xline(min(tvec(tvec>ripples(1).UP_DOWN_transition_onset(n))),'b')
    %             xline(max(tvec(tvec<ripples(1).UP_DOWN_transition_onset(n))),'r')
    %             xline(min(tvec(tvec>event_times(n,2))),'b')
    %             xline(max(tvec(tvec<event_times(n,1))),'r')
    rectangle('Position',[min(tvec(tvec>event_times2(n,1))) 0.5...
        min(tvec(tvec>event_times2(n,2)))-max(tvec(tvec<event_times2(n,1))),...
        20],...
        'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',0.2)
end

% HPC frames
for nprobe = 1:2
    probe_no = session_info(n).probe(nprobe).probe_id+1;
    if ~isempty(behavioural_state_merged.SWS)
        % SWS_no_ripples = SubtractIntervals(behavioural_state_merged.SWS,[ripples(probe_no).onset ripples(probe_no).offset]);
        [HPC_frames] = detect_candidate_frames_masa(tvec,HPC_spikes{probe_no}(:,2),behavioural_state_merged.SWS,0.01,options);
    end
end

%%%%%%%%%%% test detection
nprobe = 1;
if ~isempty(behavioural_state(probe_no).SWS)
    best_channel = find(LFP(nprobe).best_V1_channel==slow_waves(nprobe).best_channel);
    temp_SW = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters,'NREMInts',SWS);
    [temp_spindles] = FindSpindles_masa(LFP(probe_no).best_V1(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
elseif isfield(LFP(nprobe),'best_V1_high_freq')
    [~,best_channel] = max(LFP(nprobe).best_V1_high_freq);
    [temp_spindles] = FindSpindles_masa(LFP(probe_no).best_V1_high_freq(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
else
    [~,best_channel] = max(LFP(nprobe).L5_power(:,7));
    [temp_spindles] = FindSpindles_masa(LFP(probe_no).L5(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
end

