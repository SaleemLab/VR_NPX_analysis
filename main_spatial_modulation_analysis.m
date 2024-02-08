%% Main place cell and V1 spatial modulation analysis code

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

%% Spatial modulation GLM analysis
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(nprobe);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        
        clusters = clusters_ks3;

        for nprobe = 1:length(session_info(n).probe)
            metric_param = [];
            metric_param.isi_violations_ratio = @(x) x<=0.1;
            metric_param.amplitude_cutoff = @(x) x<=0.1;
            metric_param.amplitude_median = @(x) x>50;
            metric_param.peak_depth = @(x) x>max(clusters(nprobe).peak_depth)-2000 & x<max(clusters(nprobe).peak_depth-1200); % metric_param.depth_range = [] -- full range?
            [selected_clusters good_cell_index] = select_clusters(clusters(nprobe),metric_param);

            spatial_modulation_GLM_analysis(selected_clusters(nprobe),Behaviour,Task_info);
        end
    end
end



% SUBJECTS = {'M23017'}
SUBJECT = 'M23028';
SESSION = '20230706';
Stimulus_type = 'Masa2tracks';
if contains(Stimulus_type,'Masa2tracks')
    session_files = dir(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info*.mat'));
    for n = 1:length(session_files) % May have PRE RUN and POST sessions rather than just one
        load(fullfile(session_files(n).folder, session_files(n).name))
        extract_and_preprocess_NPX(session_info,Stimulus_type)
    end
else
    load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info.mat'))
    extract_and_preprocess_NPX(session_info,Stimulus_type)

end


