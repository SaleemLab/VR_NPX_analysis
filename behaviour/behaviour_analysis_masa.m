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
save(fullfile(analysis_folder,'V1-HPC behaviour','session_behaviour_summary_training.mat'))
save(fullfile(analysis_folder,'V1-HPC behaviour','session_behaviour_summary_ephys.mat'))
