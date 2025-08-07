%% Behavioural analysis

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


clear all


SUBJECTS = {'M24016','M24017','M24018','M24062','M24064','M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS,'bilateral');
Stimulus_type = 'RUN1'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

session_behaviour_summary_training = struct();  % New struct to store the extracted values
session_behaviour_summary_ephys = struct();  % New struct to store the extracted values


for iSub = 1:length(SUBJECTS)
    mouseName = SUBJECTS{iSub};
    all_session_paths = get_mouse_session_paths(mouseName);

    % loading training data
    load(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','training_behaviour.mat'),'behaviour_all')

    mkdir(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','training'));
    % Training
    for iDay = 1:length(behaviour_all.lap_ID_all)

        fields = fieldnames(behaviour_all);
        behaviour = struct();  % New struct to store the extracted values

        for f = 1:numel(fields)
            fieldname = fields{f};
            behaviour.(fieldname) = behaviour_all.(fieldname){iDay};
        end
        if isempty(behaviour.track_ID_all)
            continue
        end
        session_behaviour_summary = plot_session_behaviour_summary(behaviour,'session_name',sprintf('%s Training Day %i',mouseName,iDay));

        mkdir(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','training',sprintf('Day %i',iDay)));
        save_all_figures(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','training',sprintf('Day %i',iDay)),[]);

        fields = fieldnames(session_behaviour_summary);

        for f = 1:numel(fields)
            fieldname = fields{f};
            session_behaviour_summary_training.(fieldname){iSub}{iDay} = session_behaviour_summary.(fieldname);
        end
    end
    

    experiment_info = subject_session_stimuli_mapping({SUBJECTS{iSub}},'bilateral');

    VR_day_id = [];
    familiar_day_id = [];
    for nsession = 1:length(experiment_info)
        VR_day_id = [VR_day_id sum(contains(experiment_info(nsession).StimulusName,'RUN1')>0)];

        if sum(contains(experiment_info(nsession).StimulusName,'RUN1'))==0
            familiar_day_id = [familiar_day_id 0];

        elseif sum(contains(experiment_info(nsession).StimulusName,'RUN1'))>1
            stimuli_index = find(contains(experiment_info(nsession).StimulusName,'RUN1'));
            stimuli_index = stimuli_index(1);
            familiar_day_id = [familiar_day_id contains(experiment_info(nsession).session(stimuli_index).probe(1).trial_type,'Familiar')];
        else

            familiar_day_id = [familiar_day_id contains(experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,'RUN1')).probe(1).trial_type,'Familiar')];
        end

    end

    VR_day_id = find(VR_day_id==1);
    familiar_day_id = find(familiar_day_id==1);
    % Ephys
    for iDay = 1:length(familiar_day_id)
      
        load(fullfile(experiment_info(familiar_day_id(iDay)).ANALYSIS_DATAPATH,'Masa2tracks','extracted_task_info_RUN1.mat'))
        load(fullfile(experiment_info(familiar_day_id(iDay)).ANALYSIS_DATAPATH,'Masa2tracks','extracted_behaviour_RUN1.mat'))

        session_behaviour_summary = plot_session_behaviour_summary(Behaviour,'task_info',Task_info,'session_name',sprintf('%s ephys Day %i',mouseName,iDay));

        mkdir(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','ephys_behaviour',sprintf('Day %i',iDay)));
        save_all_figures(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouseName,'analysis','ephys_behaviour',sprintf('Day %i',iDay)),[]);

        fields = fieldnames(session_behaviour_summary);

        for f = 1:numel(fields)
            fieldname = fields{f};
            session_behaviour_summary_ephys.(fieldname){iSub}{iDay} = session_behaviour_summary.(fieldname);
        end

    end

end
session_behaviour_summary_training.x_bins = session_behaviour_summary.x_bins;
session_behaviour_summary_training.x_bin_edges = session_behaviour_summary.x_bin_edges;

session_behaviour_summary_ephys.x_bins = session_behaviour_summary.x_bins;
session_behaviour_summary_ephys.x_bin_edges = session_behaviour_summary.x_bin_edges;



if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
save(fullfile(analysis_folder,'V1-HPC behaviour','session_behaviour_summary_training.mat'),'session_behaviour_summary_training')
save(fullfile(analysis_folder,'V1-HPC behaviour','session_behaviour_summary_ephys.mat'),'session_behaviour_summary_ephys')


%% Behaviour summary
clear all

SUBJECTS = {'M24016','M24017','M24018','M24062','M24064','M24065'};
% experiment_info = subject_session_stimuli_mapping(SUBJECTS,'bilateral');
Stimulus_type = 'RUN1'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'V1-HPC behaviour','session_behaviour_summary_training.mat'))
load(fullfile(analysis_folder,'V1-HPC behaviour','session_behaviour_summary_ephys.mat'))


for iSub = 1:length(SUBJECTS)
    early_days = 1:floor(length(session_behaviour_summary_training.L_first_lick{iSub})/3);
    mid_days = early_days(end)+1:2*floor(length(session_behaviour_summary_training.L_first_lick{iSub})/3);
    late_days = mid_days(end):round(length(session_behaviour_summary_training.L_first_lick{iSub}));
    ephys_days = 1:round(length(session_behaviour_summary_ephys.L_first_lick{iSub}));

    for track_id = 1:2
        early_days_T1_lick{iSub}{track_id} = [];
        early_days_T2_lick{iSub}{track_id}= [];
        mid_days_T1_lick{iSub}{track_id}= [];
        mid_days_T2_lick{iSub}{track_id}= [];
        late_days_T1_lick{iSub}{track_id}= [];
        late_days_T2_lick{iSub}{track_id}= [];
        ephys_days_T1_lick{iSub}{track_id}= [];
        ephys_days_T2_lick{iSub}{track_id}= [];
    end


    for track_id = 1:2
        early_days_speed{iSub}{track_id} = [];
        mid_days_speed{iSub}{track_id}= [];
        late_days_speed{iSub}{track_id}= [];
        ephys_days_speed{iSub}{track_id}= [];

    end

    for track_id = 1:2
        early_days_passive{iSub}{track_id} = [];
        mid_days_passive{iSub}{track_id}= [];
        late_days_passive{iSub}{track_id}= [];
        ephys_days_passive{iSub}{track_id} = [];
    end

    for track_id = 1:2
        early_days_manual{iSub}{track_id} = [];
        mid_days_manual{iSub}{track_id}= [];
        late_days_manual{iSub}{track_id}= [];
        ephys_days_manual{iSub}{track_id} = [];
    end


    for iDay = early_days
        if ~isempty(session_behaviour_summary_training.L_first_lick{iSub}{iDay})
            manual_t1 = [session_behaviour_summary_training.manual_trial{iSub}{iDay}{1}];
            manual_t1(manual_t1==1) = nan;

            manual_t2 = [session_behaviour_summary_training.manual_trial{iSub}{iDay}{2}];
            manual_t2(manual_t2==1) = nan;

            early_days_manual{iSub}{1} = [early_days_manual{iSub}{1} manual_t1];
            early_days_manual{iSub}{2} = [early_days_manual{iSub}{2} manual_t2];
            
            early_days_T1_lick{iSub}{1} = [early_days_T1_lick{iSub}{1} session_behaviour_summary_training.L_first_lick{iSub}{iDay}{1}];
            early_days_T2_lick{iSub}{1} = [early_days_T2_lick{iSub}{1} session_behaviour_summary_training.L_first_lick{iSub}{iDay}{2}];
        end
        if ~isempty(session_behaviour_summary_training.R_first_lick{iSub}{iDay})

            early_days_T1_lick{iSub}{2} = [early_days_T1_lick{iSub}{2}  session_behaviour_summary_training.R_first_lick{iSub}{iDay}{1}];
            early_days_T2_lick{iSub}{2} = [early_days_T2_lick{iSub}{2}  session_behaviour_summary_training.R_first_lick{iSub}{iDay}{2}];
        end

        if ~isempty(session_behaviour_summary_training.lap_speed{iSub}{iDay})
            early_days_speed{iSub}{1} = [early_days_speed{iSub}{1}; session_behaviour_summary_training.lap_speed{iSub}{iDay}{1}];
            early_days_speed{iSub}{2} = [early_days_speed{iSub}{2}; session_behaviour_summary_training.lap_speed{iSub}{iDay}{2}];
        end
    end


    for iDay = mid_days
        if ~isempty(session_behaviour_summary_training.L_first_lick{iSub}{iDay})
            manual_t1 = [session_behaviour_summary_training.manual_trial{iSub}{iDay}{1}];
            manual_t1(manual_t1==1) = nan;

            manual_t2 = [session_behaviour_summary_training.manual_trial{iSub}{iDay}{2}];
            manual_t2(manual_t2==1) = nan;

            mid_days_manual{iSub}{1} = [mid_days_manual{iSub}{1} manual_t1];
            mid_days_manual{iSub}{2} = [mid_days_manual{iSub}{2} manual_t2];

            mid_days_T1_lick{iSub}{1} = [mid_days_T1_lick{iSub}{1}  session_behaviour_summary_training.L_first_lick{iSub}{iDay}{1}];
            mid_days_T2_lick{iSub}{1} = [mid_days_T2_lick{iSub}{1}  session_behaviour_summary_training.L_first_lick{iSub}{iDay}{2}];
        end
        if ~isempty(session_behaviour_summary_training.R_first_lick{iSub}{iDay})

            mid_days_T1_lick{iSub}{2} = [mid_days_T1_lick{iSub}{2}  session_behaviour_summary_training.R_first_lick{iSub}{iDay}{1}];
            mid_days_T2_lick{iSub}{2} = [mid_days_T2_lick{iSub}{2}  session_behaviour_summary_training.R_first_lick{iSub}{iDay}{2}];
        end
        if ~isempty(session_behaviour_summary_training.lap_speed{iSub}{iDay})
            mid_days_speed{iSub}{1} = [mid_days_speed{iSub}{1}; session_behaviour_summary_training.lap_speed{iSub}{iDay}{1}];
            mid_days_speed{iSub}{2} = [mid_days_speed{iSub}{2}; session_behaviour_summary_training.lap_speed{iSub}{iDay}{2}];
        end
    end



    for iDay = late_days
        if ~isempty(session_behaviour_summary_training.L_first_lick{iSub}{iDay})
            manual_t1 = [session_behaviour_summary_training.manual_trial{iSub}{iDay}{1}];
            manual_t1(manual_t1==1) = nan;

            manual_t2 = [session_behaviour_summary_training.manual_trial{iSub}{iDay}{2}];
            manual_t2(manual_t2==1) = nan;

            late_days_manual{iSub}{1} = [late_days_manual{iSub}{1} manual_t1];
            late_days_manual{iSub}{2} = [late_days_manual{iSub}{2} manual_t2];

            late_days_T1_lick{iSub}{1} = [late_days_T1_lick{iSub}{1}  session_behaviour_summary_training.L_first_lick{iSub}{iDay}{1}];
            late_days_T2_lick{iSub}{1} = [late_days_T2_lick{iSub}{1}  session_behaviour_summary_training.L_first_lick{iSub}{iDay}{2}];
        end
        if ~isempty(session_behaviour_summary_training.R_first_lick{iSub}{iDay})

            late_days_T1_lick{iSub}{2} = [late_days_T1_lick{iSub}{2}  session_behaviour_summary_training.R_first_lick{iSub}{iDay}{1}];
            late_days_T2_lick{iSub}{2} = [late_days_T2_lick{iSub}{2}  session_behaviour_summary_training.R_first_lick{iSub}{iDay}{2}];
        end
        if ~isempty(session_behaviour_summary_training.lap_speed{iSub}{iDay})
            late_days_speed{iSub}{1} = [late_days_speed{iSub}{1}; session_behaviour_summary_training.lap_speed{iSub}{iDay}{1}];
            late_days_speed{iSub}{2} = [late_days_speed{iSub}{2}; session_behaviour_summary_training.lap_speed{iSub}{iDay}{2}];
        end
    end



    for iDay = ephys_days
        if ~isempty(session_behaviour_summary_ephys.L_first_lick{iSub}{iDay})

            ephys_days_T1_lick{iSub}{1} = [ephys_days_T1_lick{iSub}{1}  session_behaviour_summary_ephys.L_first_lick{iSub}{iDay}{1}];
            ephys_days_T2_lick{iSub}{1} = [ephys_days_T2_lick{iSub}{1}  session_behaviour_summary_ephys.L_first_lick{iSub}{iDay}{2}];
        end
        
        
        if iDay<=length(session_behaviour_summary_ephys.manual_trial{iSub})
            manual_t1 = [session_behaviour_summary_ephys.manual_trial{iSub}{iDay}{1}];
            manual_t1(manual_t1==1) = nan;

            manual_t2 = [session_behaviour_summary_ephys.manual_trial{iSub}{iDay}{2}];
            manual_t2(manual_t2==1) = nan;

            ephys_days_manual{iSub}{1} = [ephys_days_manual{iSub}{1} manual_t1];
            ephys_days_manual{iSub}{2} = [ephys_days_manual{iSub}{2} manual_t2];
        else
            manual_t1 = nan(1,length(session_behaviour_summary_ephys.L_first_lick{iSub}{iDay}{1}));
            manual_t2 = nan(1,length(session_behaviour_summary_ephys.L_first_lick{iSub}{iDay}{2}));

            ephys_days_manual{iSub}{1} = [ephys_days_manual{iSub}{1} manual_t1];
            ephys_days_manual{iSub}{2} = [ephys_days_manual{iSub}{2} manual_t2];
        end


        if ~isempty(session_behaviour_summary_ephys.R_first_lick{iSub}{iDay})

            ephys_days_T1_lick{iSub}{2} = [ephys_days_T1_lick{iSub}{2}  session_behaviour_summary_ephys.R_first_lick{iSub}{iDay}{1}];
            ephys_days_T2_lick{iSub}{2} = [ephys_days_T2_lick{iSub}{2}  session_behaviour_summary_ephys.R_first_lick{iSub}{iDay}{2}];
        end
        if ~isempty(session_behaviour_summary_ephys.lap_speed{iSub}{iDay})
            ephys_days_speed{iSub}{1} = [ephys_days_speed{iSub}{1}; session_behaviour_summary_ephys.lap_speed{iSub}{iDay}{1}];
            ephys_days_speed{iSub}{2} = [ephys_days_speed{iSub}{2}; session_behaviour_summary_ephys.lap_speed{iSub}{iDay}{2}];
        end
    end

end



%%%%%%%% Acorss animal running behaviour

% Ratio of first licks in the correct track
early_days_T1_lick{iSub}{1}

early_days_T1_lick{iSub}{2}



late_days_T1_lick{iSub}{1}

late_days_T1_lick{iSub}{2}



% Speed not consitently different across tracks

% Track L speed
mean_speed = [];

colorlines = [ ...
    0.2, 0.4, 0.8;    % blue
    0.85, 0.2, 0.2;   % red
];

for iSub = 1:length(early_days_speed)
    fig = figure;
    fig.Name = sprintf('%s running speed',SUBJECTS{iSub})
    
    nexttile
    histogram(mean(early_days_speed{iSub}{1},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(mean(early_days_speed{iSub}{2},2,'omitnan'),20,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));

    legend('Track Left','Track R','box','off')
    xlabel('Mean lap running speed (cm/s)')
    xlabel('Proportion of laps')
    set(gca,'TickDir','out','Box','off','FontSize',12)

    nexttile
    histogram(mean(mid_days_speed{iSub}{1},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(mean(mid_days_speed{iSub}{2},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));
    legend('Track Left','Track R','box','off')
    xlabel('Mean lap running speed (cm/s)')
    xlabel('Proportion of laps')
    set(gca,'TickDir','out','Box','off','FontSize',12)

    nexttile
    histogram(mean(late_days_speed{iSub}{1},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(mean(late_days_speed{iSub}{2},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));
    legend('Track Left','Track R','box','off')
    xlabel('Mean lap running speed (cm/s)')
    xlabel('Proportion of laps')
    set(gca,'TickDir','out','Box','off','FontSize',12)

    nexttile
    histogram(mean(ephys_days_speed{iSub}{1},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(mean(ephys_days_speed{iSub}{2},2,'omitnan'),0:2:45,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));
    legend('Track Left','Track R','box','off')
    xlabel('Mean lap running speed (cm/s)')
    xlabel('Proportion of laps')
    set(gca,'TickDir','out','Box','off','FontSize',12)

end




colorlines = [ ...
    0.2, 0.4, 0.8;    % blue
    0.85, 0.2, 0.2;   % red
];
%%%% Lick 
for iSub = 1:length(early_days_speed)
    fig = figure;
    fig.Name = sprintf('%s first lick distribution',SUBJECTS{iSub})
    fig.Position= [43 55 494 891];

    nexttile
    histogram(early_days_T1_lick{iSub}{1}(~isnan(early_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(early_days_T1_lick{iSub}{2}(~isnan(early_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track L','Track R','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Early days Track L')

    nexttile
    histogram(early_days_T2_lick{iSub}{2}(~isnan(early_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    histogram(early_days_T2_lick{iSub}{1}(~isnan(early_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track R','Track L','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Early days Track R')


    
    nexttile
    histogram(mid_days_T1_lick{iSub}{1}(~isnan(mid_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(mid_days_T1_lick{iSub}{2}(~isnan(mid_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track L','Track R','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Mid days Track L')

    nexttile
    histogram(mid_days_T2_lick{iSub}{2}(~isnan(mid_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    histogram(mid_days_T2_lick{iSub}{1}(~isnan(mid_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track R','Track L','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Mid days Track R')




    nexttile
    histogram(late_days_T1_lick{iSub}{1}(~isnan(late_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(late_days_T1_lick{iSub}{2}(~isnan(late_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track L','Track R','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Late days Track L')

    nexttile
    histogram(late_days_T2_lick{iSub}{2}(~isnan(late_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    histogram(late_days_T2_lick{iSub}{1}(~isnan(late_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track R','Track L','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Late days Track R')




    nexttile
    histogram(ephys_days_T1_lick{iSub}{1}(~isnan(ephys_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    histogram(ephys_days_T1_lick{iSub}{2}(~isnan(ephys_days_manual{iSub}{1})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track L','Track R','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Ephys days Track L')

    nexttile
    histogram(ephys_days_T2_lick{iSub}{2}(~isnan(ephys_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(2,:));hold on;
    histogram(ephys_days_T2_lick{iSub}{1}(~isnan(ephys_days_manual{iSub}{2})),0:2:140,'Normalization','probability','FaceAlpha',0.5,'FaceColor',colorlines(1,:));hold on;
    xline(100,'k','LineWidth',2)
    legend('Track R','Track L','box','off')
    xlabel('Position (cm)')
    xlabel('Proportion of first licks')
    set(gca,'TickDir','out','Box','off','FontSize',12)
    title('Ephys days Track R')
end

save_all_figures(fullfile(analysis_folder,'V1-HPC behaviour'),[])



%%%% Lick
prop_anticipatory_correct=[];
p_anticipatory_lick_dist=[];
prop_correct=[];
p_lick_dist=[];
first_side=[];
first_pos = [];

for iSub = 1:length(early_days_speed)
    first_side=[];
    first_pos = [];

    correct_licking = [early_days_T1_lick{iSub}{1} early_days_T2_lick{iSub}{2}];
    incorrect_licking = [early_days_T1_lick{iSub}{2} early_days_T2_lick{iSub}{1}];

    % Identify which lick came first per lap
    first_pos{1}  = nan(1, length(correct_licking));
    first_side{1} = nan(1, length(correct_licking));  % 1 = L, 2 = R

    for i = 1:length(correct_licking)
        if ~isnan(correct_licking(i)) && (isnan(incorrect_licking(i))) || correct_licking(i) < incorrect_licking(i)
            first_pos{1}(i) = correct_licking(i);
            first_side{1}(i) = 1;
        elseif isnan(correct_licking(i)) && (~isnan(incorrect_licking(i)))  || correct_licking(i) >incorrect_licking(i)
            first_pos{1}(i) = incorrect_licking(i);
            first_side{1}(i) = 0;
        end
    end

    
    correct_licking = [mid_days_T1_lick{iSub}{1} mid_days_T2_lick{iSub}{2}];
    incorrect_licking = [mid_days_T1_lick{iSub}{2} mid_days_T2_lick{iSub}{1}];

    first_pos{2}  = nan(1, length(correct_licking));
    first_side{2} = nan(1, length(correct_licking));  % 1 = L, 2 = R

    for i = 1:length(correct_licking)
        if ~isnan(correct_licking(i)) && (isnan(incorrect_licking(i))) || correct_licking(i) < incorrect_licking(i)
            first_pos{2}(i) = correct_licking(i);
            first_side{2}(i) = 1;
        elseif isnan(correct_licking(i)) && (~isnan(incorrect_licking(i)))  || correct_licking(i) >incorrect_licking(i)
            first_pos{2}(i) = incorrect_licking(i);
            first_side{2}(i) = 0;
        end
    end



    correct_licking = [late_days_T1_lick{iSub}{1} late_days_T2_lick{iSub}{2}];
    incorrect_licking = [late_days_T1_lick{iSub}{2} late_days_T2_lick{iSub}{1}];

    first_pos{3}  = nan(1, length(correct_licking));
    first_side{3} = nan(1, length(correct_licking));  % 1 = L, 2 = R

    for i = 1:length(correct_licking)
        if ~isnan(correct_licking(i)) && (isnan(incorrect_licking(i))) || correct_licking(i) < incorrect_licking(i)
            first_pos{3}(i) = correct_licking(i);
            first_side{3}(i) = 1;
        elseif isnan(correct_licking(i)) && (~isnan(incorrect_licking(i)))  || correct_licking(i) >incorrect_licking(i)
            first_pos{3}(i) = incorrect_licking(i);
            first_side{3}(i) = 0;
        end
    end


    correct_licking = [ephys_days_T1_lick{iSub}{1} ephys_days_T2_lick{iSub}{2}];
    incorrect_licking = [ephys_days_T1_lick{iSub}{2} ephys_days_T2_lick{iSub}{1}];

    first_pos{4}  = nan(1, length(correct_licking));
    first_side{4} = nan(1, length(correct_licking));  % 1 = L, 2 = R

    for i = 1:length(correct_licking)
        if ~isnan(correct_licking(i)) && (isnan(incorrect_licking(i))) || correct_licking(i) < incorrect_licking(i)
            first_pos{4}(i) = correct_licking(i);
            first_side{4}(i) = 1;
        elseif isnan(correct_licking(i)) && (~isnan(incorrect_licking(i)))  || correct_licking(i) >incorrect_licking(i)
            first_pos{4}(i) = incorrect_licking(i);
            first_side{4}(i) = 0;
        end
    end
    % Calculate hit and FA rates

    % N_valid = sum(is_valid);
    
    for ntype = 1:4
        for iBoot = 1:1000

            s = RandStream('philox4x32_10', 'Seed', iBoot);

            resampled_idx = datasample(s,1:length(first_side{ntype}), length(first_side{ntype}), 'Replace', true);
            is_valid = ~isnan(first_side{ntype}(resampled_idx));  % only laps with at least one lick

            anticip = first_pos{ntype}(resampled_idx) < 100;

            N_valid = sum(is_valid & anticip);

            hits = sum(first_side{ntype}(resampled_idx) == 1 & anticip& is_valid);           % context-matched licks
            fas  = sum(first_side{ntype}(resampled_idx) ~= 1 & anticip& is_valid);  % contex-mismatched licks

            hit_rate = (hits + 0.5) / (N_valid + 1);
            fa_rate  = (fas + 0.5) / (N_valid + 1);

     
            % D prime
            d_prime(iSub,iBoot,ntype) = norminv(hit_rate) - norminv(fa_rate);
        end
    end



    prop_anticipatory_correct(iSub,1) =(sum(early_days_T2_lick{iSub}{2}(~isnan(early_days_manual{iSub}{2}) & ~(early_days_T2_lick{iSub}{1} <= 100)) <= 100) ...
        + sum(early_days_T1_lick{iSub}{1}(~isnan(early_days_manual{iSub}{1}) & ~(early_days_T1_lick{iSub}{2} <= 100)) <= 100))/ ...
        (length(early_days_T2_lick{iSub}{1})+length(early_days_T1_lick{iSub}{1}));

    prop_anticipatory_correct(iSub,2) =(sum(mid_days_T2_lick{iSub}{2}(~isnan(mid_days_manual{iSub}{2}) & ~(mid_days_T2_lick{iSub}{1} <= 100)) <= 100) ...
        + sum(mid_days_T1_lick{iSub}{1}(~isnan(mid_days_manual{iSub}{1}) & ~(mid_days_T1_lick{iSub}{2} <= 100)) <= 100))/ ...
        (length(mid_days_T2_lick{iSub}{1})+length(mid_days_T1_lick{iSub}{1}));

    prop_anticipatory_correct(iSub,3) =(sum(late_days_T2_lick{iSub}{2}(~isnan(late_days_manual{iSub}{2}) & ~(late_days_T2_lick{iSub}{1} <= 100)) <= 100) ...
        + sum(late_days_T1_lick{iSub}{1}(~isnan(late_days_manual{iSub}{1}) &  ~(late_days_T1_lick{iSub}{2} <= 100)) <= 100))/ ...
        (length(late_days_T2_lick{iSub}{1})+length(late_days_T1_lick{iSub}{1}));

    prop_anticipatory_correct(iSub,4) =(sum(ephys_days_T2_lick{iSub}{2}(~isnan(ephys_days_manual{iSub}{2}) & ~(ephys_days_T2_lick{iSub}{1} <= 100)) <= 100) ...
        + sum(ephys_days_T1_lick{iSub}{1}(~isnan(ephys_days_manual{iSub}{1}) & ~(ephys_days_T1_lick{iSub}{2} <= 100)) <= 100))/ ...
        (length(ephys_days_T2_lick{iSub}{1})+length(ephys_days_T1_lick{iSub}{1}));



    [h,p,ks2stat] = kstest2([mid_days_T1_lick{iSub}{1}(~isnan(mid_days_manual{iSub}{1})) mid_days_T2_lick{iSub}{2}(~isnan(mid_days_manual{iSub}{2}))],...
        [early_days_T1_lick{iSub}{1}(~isnan(early_days_manual{iSub}{1})) early_days_T2_lick{iSub}{2}(~isnan(early_days_manual{iSub}{2}))],'Tail','larger');
    p_anticipatory_lick_dist(iSub,1) = p;

    [h,p,ks2stat] = kstest2([late_days_T1_lick{iSub}{1}(~isnan(late_days_manual{iSub}{1})) late_days_T2_lick{iSub}{2}(~isnan(late_days_manual{iSub}{2}))],...
        [early_days_T1_lick{iSub}{1}(~isnan(early_days_manual{iSub}{1})) early_days_T2_lick{iSub}{2}(~isnan(early_days_manual{iSub}{2}))],'Tail','larger');
    p_anticipatory_lick_dist(iSub,2) = p;

    [h,p,ks2stat] = kstest2([ephys_days_T1_lick{iSub}{1}(~isnan(ephys_days_manual{iSub}{1})) ephys_days_T2_lick{iSub}{2}(~isnan(ephys_days_manual{iSub}{2}))],...
        [early_days_T1_lick{iSub}{1}(~isnan(early_days_manual{iSub}{1})) early_days_T2_lick{iSub}{2}(~isnan(early_days_manual{iSub}{2}))],'Tail','larger');
    p_anticipatory_lick_dist(iSub,3) = p;

end


colorlines = [ ...
    213,  62,  79;   % red
    244, 109,  67;   % orange-red
    253, 174,  97;   % light orange
    171, 221, 164;   % light green
    102, 194, 165;   % teal green
     50, 136, 189    % blue
] / 255;




fig = figure;
fig.Name = 'Proportion of context-matched first lick';
fig.Position= [43 55 494 891];
nexttile
for iSub = 1:6
    plot(prop_anticipatory_correct(iSub,:),'Marker','.','MarkerSize',40,'Color',colorlines(iSub,:));hold on;
end
legend(SUBJECTS,'box','off')
xticks([1 2 3 4])
xticklabels({'Early','Mid','Late','Recording'})
ylabel('proportion of context-matched anticipatory licking')
set(gca,'TickDir','out','Box','off','FontSize',12)
xlim([0 5])

nexttile

for iSub = 1:6
    for ntype = 1:4
        mu(ntype) = mean(d_prime(iSub,:,ntype));
        CI(ntype,:) =  prctile(d_prime(iSub,:,ntype), [2.5, 97.5]);

        % plot(,'Marker','.','MarkerSize',40,);hold on;
    end
    hold on
    E(iSub)=errorbar(1+iSub*0.2:2:8+iSub*0.2, mu, mu - CI(:,1)', CI(:,2)' - mu, 'o', 'CapSize', 5,'LineWidth', 1.5,'Color',colorlines(iSub,:));
end
legend(SUBJECTS,'box','off')
xticks(1.8:2:8.8)
xticklabels({'Early','Mid','Late','Recording'})
ylabel('d prime')
set(gca,'TickDir','out','Box','off','FontSize',12)
xlim([0 10])

save_all_figures(fullfile(analysis_folder,'V1-HPC behaviour'),[])



prop_correct


