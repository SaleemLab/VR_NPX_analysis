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
session_count = 0;

slow_waves_TF_all = struct();
ripples_TF_all = struct();
spindles_TF_all = struct();

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));


    if length(stimulus_name)>1 % Based on if POST or PRE
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
        else

            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));

            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
            else
                n =1;
            end
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
    options = session_info(n).probe(1);

    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
    if isempty(DIR)
        continue
    end

    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
    DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

    if ~isempty(DIR)
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        session_clusters_RUN=session_clusters;
        clear session_clusters
    end

    if ~isempty(DIR1)
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        session_clusters_RUN=session_clusters;
        clear session_clusters
    end

    if contains(stimulus_name{n},'Masa2tracks')
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        % %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         clusters=clusters_ks4;
    elseif contains(stimulus_name{n},'Sleep')
        %         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    else
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
        %         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    end

end

