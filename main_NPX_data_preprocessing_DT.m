%% Putative pipeline for processing NPX1 data recorded from hippocampus and V1
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)

% main_NPX_data_preprocessing

% Go to SUA_analysis_masa for kilosort + cell explorer cell classification
% Go to CellExplorerTest_masa_adapted for cell exploerer pipeline


%% Set the data folders and processing parameters
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))

addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
% addpath('Z:\ibn-vision\USERS\Masa\code\Masa_utility')
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\NPXAnalysis\NPXAnalysis2022'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\visual_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'));

%     ROOTPATH = 'X:\ibn-vision';
ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive

all_SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};

for n = 1:length(all_SUBJECTS)
    % extract information about this animal
    SUBJECTS = {all_SUBJECTS{n}};
    experiment_info = subject_session_stimuli_mapping(SUBJECTS,'V1-MEC');

    if exist(fullfile(ROOTPATH,'DATA','SUBJECTS',all_SUBJECTS{n},'analysis')) == 0
        mkdir(fullfile(ROOTPATH,'DATA','SUBJECTS',all_SUBJECTS{n},'analysis'))
    end
    % save experiment info into analysis folder
    save(fullfile(ROOTPATH,'DATA','SUBJECTS',all_SUBJECTS{n},'analysis','experiment_info'),'experiment_info')
    
    % For each session, loop through all stimuli
    for nsession = 1:length(experiment_info)
        for nstimuli = 1:length(experiment_info(nsession).session)
            clear session_info
%             for nprobe = 1:length(experiment_info(nsession).stimuli_type(nstimuli).probe)
%                 session_info(nprobe) = experiment_info(nsession).stimuli_type(nstimuli).probe(nprobe);
%             end

            session_info = experiment_info(nsession).session(nstimuli);

            if exist(session_info.probe(1).ANALYSIS_DATAPATH) == 0
                mkdir(session_info.probe(1).ANALYSIS_DATAPATH)
            end
            
            stimulus_name = session_info.probe(1).StimulusName;
            if contains(stimulus_name,'Masa2tracks')
                stimulus_name = 'Masa2tracks';
            end

            for nprobe = 1:length(session_info.probe)
                DIR = dir(fullfile(session_info.probe(nprobe).EPHYS_DATAPATH,'*.meta'));
                for nfile = 1:length(DIR)
                    copyfile(fullfile(DIR(nfile).folder,DIR(nfile).name),fullfile(session_info.probe(nprobe).ANALYSIS_DATAPATH,DIR(nfile).name))
                end
            end

            if contains(session_info.probe(1).StimulusName,'Masa2tracks')
                save(fullfile(session_info.probe(1).ANALYSIS_DATAPATH,sprintf('session_info%s',erase(session_info.probe(1).StimulusName,'Masa2tracks'))),'session_info')
            else
                save(fullfile(session_info.probe(1).ANALYSIS_DATAPATH,'session_info'),'session_info')
            end

        end
    end
end

%% import and align and store Bonsai and cluster spike data

addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
%%%%%% Option 1 use subject_session_stimuli_mapping for all animals you
%%%%%% want to process. 
SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};
options = 'V1-MEC';
ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive

% Stimulus_type = 'Masa2tracks';
Stimulus_type = 'Checkerboard';

experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);
% All_stimuli = {'Masa2tracks','SparseNoise_fullscreen','Checkerboard','StaticGratings'}
All_stimuli = {'SparseNoise_fullscreen','Checkerboard','StaticGratings'}
% experiment_info = experiment_info(2)
for n = 1:length(All_stimuli)
    extract_and_preprocess_NPX_batch(experiment_info,All_stimuli{n})
end



%%%%%% Option 2 go to specific animal folder to do specific session(s) you
%%%%%% want to process.
% SUBJECTS = {'M23017'}
SUBJECT = ['M23031'];

SESSION = {['20230714']};

for iMouse = 1
    load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT(iMouse,:),'analysis','experiment_info.mat'))
    for iSession = 1:size(SESSION{iMouse},1)
        Stimulus_type = experiment_info(iSession).StimulusName;
        for iStim = 1:length(Stimulus_type)
            if ~contains(Stimulus_type{iStim},'StaticGrating')
                load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT(iMouse,:),'analysis',SESSION{iMouse}(iSession,:),Stimulus_type{iStim},'session_info.mat'))
                extract_and_preprocess_NPX(session_info,Stimulus_type{iStim})
            end
        end
    end
end
%% PSD analysis and LFP profile
ROOTPATH = 'Z:\ibn-vision';

addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
% Single session
SUBJECT = ['M23037'];
SESSION = {['20230810';'20230811';'20230812';'20230813']};
% SUBJECT = ['M23032';'M23034';'M23038'];
% 
% SESSION = {['20230718';'20230719';'20230720';'20230721';'20230722'];
%     ['20230804';'20230805';'20230806'];
%     ['20230816';'20230817']};
options = 'bilateral';
Stimulus_type = 'Checkerboard';
% Stimulus_type = 'OpenField';
for iMouse = 1:3
    load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT(iMouse,:),'analysis','experiment_info.mat'))
    for iSession = 1:size(SESSION{iMouse},1)
if contains(Stimulus_type,'Masa2tracks')
    session_files = dir(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info*.mat'));
    for n = 1:length(session_files) % May have PRE RUN and POST sessions rather than just one
        load(fullfile(session_files(n).folder, session_files(n).name))
        extract_PSD_profile(session_info,Stimulus_type)
    end
else
    load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT(iMouse,:),'analysis',SESSION{iMouse}(iSession,:),Stimulus_type,'session_info.mat'))
    extract_PSD_profile(session_info,Stimulus_type)
end
    end
end

% Batch PSD analysis
Stimulus_type = 'Checkerboard'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
% SUBJECTS = {'M23028'};
% SUBJECTS = {'M23087'};
all_SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};
% SUBJECTS = {'M23029','M23087','M23153'};
experiment_info = subject_session_stimuli_mapping(all_SUBJECTS,'V1-MEC');
% experiment_info = experiment_info(end);
extract_PSD_profile_batch(experiment_info,Stimulus_type);

%% Determine L4 of V1 based on checkerboard (require manual updating)
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

% Single session checkerboard
ROOTPATH = 'Z:\ibn-vision';
SUBJECT = ['M23037';'M23038'];

SESSION = {
    
    ['20230812';'20230813'];
    ['20230816';'20230817']};

Stimulus_type = 'Checkerboard';
for iMouse = 1:4
    
    for iSession = 1:size(SESSION{iMouse},1)
load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT(iMouse,:),'analysis',SESSION{iMouse}(iSession,:),Stimulus_type,'session_info.mat'))

        for nprobe = 2:length(session_info.probe) % For each session, how many probes
            options= session_info.probe(nprobe);
            %             options.ROOTPATH = ROOTPATH;
            options.importMode = 'KS';
            options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            
            [lfpAvg(options.probe_no).column,csd(options.probe_no).column,PSD,best_channels] = checkerboard_CSD_profile(options);
        
            save_all_figures(options.ANALYSIS_DATAPATH,[]);
            close all
        
            save(fullfile(options.ANALYSIS_DATAPATH,"checkerboard_CSD.mat"),'lfpAvg','csd');
        
            [LF_FILE imecMeta chan_config ~] = extract_NPX_channel_config(options,[]);% Since it is LF
            [best_channels{options.probe_no}] = update_best_channels(options,chan_config);
        
            save(fullfile(options.ANALYSIS_DATAPATH,'..',"best_channels.mat"),'best_channels')
            
            power = [];
            xcoord = [];
            ycoord = [];
        
            for nchannel = 1:size(PSD{options.probe_no},2)
                power(nchannel,:) = PSD{options.probe_no}(nchannel).mean_power;
                xcoord(nchannel) = PSD{options.probe_no}(nchannel).xcoord;
                ycoord(nchannel) = PSD{options.probe_no}(nchannel).ycoord;
            end
        
            % Replot based on updated channels
            for col = 1:length(lfpAvg(options.probe_no).column)
                xcoord_avaliable = lfpAvg(options.probe_no).column(col).xcoord;
                plot_perievent_CSD_LFP(lfpAvg(options.probe_no).column(col),csd(options.probe_no).column(col),power(xcoord == xcoord_avaliable,:),chan_config,chan_config(xcoord == xcoord_avaliable,:),best_channels{options.probe_no},options)
            end
            save_all_figures(options.ANALYSIS_DATAPATH,[]);
        end
    end
end
% Checkerboard CSD batch
% SUBJECTS = {'M23017','M23028','M23029'};
SUBJECTS = {'M23087'};
options = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);
Stimulus_type = 'Checkerboard';
% determine_best_channels
calculate_checkerboard_CSD_profile_batch(experiment_info,Stimulus_type)

{'surface','L4','L5','CA1'}