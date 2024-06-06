%% Putative pipeline for processing chronic NPX1 and NPX2s data recorded from hippocampus and V1
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)

% main_NPX_data_preprocessing

% Go to SUA_analysis_masa for kilosort + cell explorer cell classification
% Go to CellExplorerTest_masa_adapted for cell exploerer pipeline


%% Set the data folders and processing parameters
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
% addpath('Z:\ibn-vision\USERS\Masa\code\Masa_utility')
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\NPXAnalysis\NPXAnalysis2022'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\visual_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'));

%     ROOTPATH = 'X:\ibn-vision';
ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive

all_SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};

all_SUBJECTS = {'M24017'};
for n = 1:length(all_SUBJECTS)
    % extract information about this animal
    SUBJECTS = {all_SUBJECTS{n}};
    experiment_info = subject_session_stimuli_mapping(SUBJECTS,'bilateral');

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

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

%%%%%% Option 1 use subject_session_stimuli_mapping for all animals you
%%%%%% want to process. 
SUBJECTS = {'M23017','M23029','M23087','M23153'};
SUBJECTS = {'M23028','M23087','M23153'};
SUBJECTS = {'M24017'};
options = 'bilateral';
ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive

% Stimulus_type = 'Masa2tracks';
% Stimulus_type = 'Checkerboard';

experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);
% All_stimuli = {'FullScreenFlash'}
% All_stimuli = {'SparseNoise_fullscreen','Checkerboard','StaticGratings'}
% experiment_info = experiment_info(4)

All_stimuli = {'Masa2tracks','SparseNoise','Checkerboard','StaticGratings'};
for n = 1:length(All_stimuli)
    extract_and_preprocess_NPX_batch(experiment_info,All_stimuli{n})
end

SUBJECTS = {'M24017'};
options = 'bilateral';
Stimulus_type = 'Checkerboard';
% Stimulus_type = 'SparseNoise';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);
% experiment_info = experiment_info(1);
extract_and_preprocess_NPX_batch(experiment_info,Stimulus_type)

Stimulus_type = 'Masa2tracks';
%%%%%% Option 2 go to specific animal folder to do specific session(s) you
%%%%%% want to process.
% SUBJECTS = {'M23017'}
SUBJECT = 'M23028';
SESSION = '20230706';

% SUBJECT = 'M23087';
SUBJECT = 'M23153';
% SESSION = '20231212';
SESSION = '20231212';
options = 'bilateral';
% Stimulus_type = 'Masa2tracks';
Stimulus_type = 'Checkerboard';
Stimulus_type = 'SparseNoise_fullscreen';
% Stimulus_type = 'OpenFieldChronic';
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


%% PSD analysis and LFP profile
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
clear all
ROOTPATH = 'Z:\ibn-vision';
% Single session
SUBJECT = 'M24019';
SESSION = '202405';
options = 'bilateral';
Stimulus_type = 'Checkerboard';
% Stimulus_type = 'OpenField';
if contains(Stimulus_type,'Masa2tracks')
    session_files = dir(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info*.mat'));
    for n = 1:length(session_files) % May have PRE RUN and POST sessions rather than just one
        load(fullfile(session_files(n).folder, session_files(n).name))
        extract_PSD_profile(session_info,Stimulus_type)
    end
else
    load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info.mat'))
    extract_PSD_profile(session_info,Stimulus_type)
end


% Batch PSD analysis
Stimulus_type = 'Checkerboard'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
% SUBJECTS = {'M23028'};
% SUBJECTS = {'M23087'};
SUBJECTS = {'M23017','M23028','M23029','M23087'};
SUBJECTS = {'M23087'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS,'bilateral');
% experiment_info = experiment_info(1);
Stimulus_type= 'Checkerboard_sh1'; 
extract_PSD_profile_batch(experiment_info,Stimulus_type);
Stimulus_type= 'Checkerboard_sh2'; 
extract_PSD_profile_batch(experiment_info,Stimulus_type);
Stimulus_type= 'Checkerboard_sh3'; 
extract_PSD_profile_batch(experiment_info,Stimulus_type);
Stimulus_type= 'Checkerboard_sh4'; 
extract_PSD_profile_batch(experiment_info,Stimulus_type);


%% Determine L4 of V1 based on checkerboard (require manual updating)
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))

clear all

% Single session checkerboard
ROOTPATH = 'Z:\ibn-vision';
SUBJECT = 'M23028';
SESSION = '20230703';
options = 'bilateral';

% Stimulus_type = 'FullScreenFlash_2';
Stimulus_type = 'Checkerboard';
% Stimulus_type = 'RUN';
load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info.mat'))

for nprobe = 1:length(session_info.probe) % For each session, how many probes
    options= session_info.probe(nprobe);
    %             options.ROOTPATH = ROOTPATH;
    options.importMode = 'KS';
    options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
    
%     DIR = dir(fullfile(options.ANALYSIS_DATAPATH,"checkerboard_CSD.mat"))
% 
%     if ~isempty(DIR)
%         load(fullfile(options.ANALYSIS_DATAPATH,"checkerboard_CSD.mat"),'lfpAvg','csd');
%     end

    [lfpAvg(options.probe_no).column,csd(options.probe_no).column,PSD,best_channels] = checkerboard_CSD_profile(options);
    
    save_all_figures(options.ANALYSIS_DATAPATH,[]);
    close all

    save(fullfile(options.ANALYSIS_DATAPATH,"checkerboard_CSD.mat"),'lfpAvg','csd');

    [LF_FILE imecMeta chan_config ~] = extract_NPX_channel_config(options,[]);% Since it is LF
    [best_channels{options.probe_no}] = update_best_channels(options,chan_config);

    save(fullfile(options.ANALYSIS_DATAPATH,'..',"best_channels.mat"),'best_channels')

    %     power = [];
    %     xcoord = [];
    %     ycoord = [];
    %
    %     for nchannel = 1:size(PSD{options.probe_no},2)
    %         power(nchannel,:) = PSD{options.probe_no}(nchannel).mean_power;
    %         xcoord(nchannel) = PSD{options.probe_no}(nchannel).xcoord;
    %         ycoord(nchannel) = PSD{options.probe_no}(nchannel).ycoord;
    %     end
    %
    %     % sort channel according to y coordinate
    %     [ycoord idx] = sort(ycoord,'ascend');
    %     xcoord = xcoord(idx);
    %     power = power(idx,:);
    %     chan_config = chan_config(idx,:);
    %
    %     % Replot based on updated channels
    %     for col = 1:length(lfpAvg(options.probe_no).column)
    %         xcoord_avaliable = lfpAvg(options.probe_no).column(col).xcoord;
    %         plot_perievent_CSD_LFP(lfpAvg(options.probe_no).column(col),csd(options.probe_no).column(col),power(xcoord == xcoord_avaliable,:),chan_config,chan_config(xcoord == xcoord_avaliable,:),best_channels{options.probe_no},options)
    %     end
    checkerboard_CSD_profile(options);
    save_all_figures(options.ANALYSIS_DATAPATH,[]);
end

% Checkerboard CSD batch
% SUBJECTS = {'M23017','M23028','M23029'};
SUBJECTS = {'M23087'};
options = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);
Stimulus_type = 'Checkerboard';
% determine_best_channels
calculate_checkerboard_CSD_profile_batch(experiment_info,Stimulus_type)


%% Visual tuning based on Static Gratings
Stimulus_type = 'StaticGratings';
for nsession =10:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for nprobe = 1:length(session_info.probe) % For each session, how many probes
        options = session_info.probe(nprobe);
        options.BinWidth = 1/60;
        options.importMode = 'KS'; % LF or MUA or KS
        % options.importMode = 'LF'; % LF or MUA or KS
        options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
        options.stim_dur = 0.1;
        options.AnalysisTimeWindow = [0 2];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
        options.ks_unitType = 'good'; % 'mua', 'good' or ''
        options.PD_FLAG = 1;
        options.paradigm = 'SG';
        options.gFileNum = gFileNum;
        options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

        [resps,otherData,stimData,~,wheelData,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
        StaticGratings.probe(options.probe_no).resps=resps;
        StaticGratings.probe(options.probe_no).otherData=otherData;
        StaticGratings.probe(options.probe_no).stimData=stimData;
        StaticGratings.probe(options.probe_no).wheelData=wheelData;
        StaticGratings.probe(options.probe_no).photodiodeData=photodiodeData;
        StaticGratings.probe(options.probe_no).timeVector=timeVector;
        StaticGratings.probe(options.probe_no).options=options;
        switch(stimulus_name{1,1})
            case {'StaticGratings_short','StaticGratings'}
                stim_orientation = readmatrix('X:\ibn-vision\CODE\DEV\BONSAI\Diao\dome_dual_DT\Grating_trials_short.CSV');
            case {'StaticGratings_long'}
                stim_orientation = readmatrix('X:\ibn-vision\CODE\DEV\BONSAI\Diao\dome_dual_DT\Grating_trials.CSV');
        end
        orientation_angles = unique(stim_orientation);
        stim_orientation_tmp = stim_orientation;
        stim_orientation_tmp(wheelData.staticgrating_idx_error) = -1;
        stim_orientation_wo_error = find(stim_orientation_tmp>-1);
        stim_index_wo_error = zeros(200,1);
        stim_index_wo_error(stim_orientation_wo_error) = 1;
        grating_response = cell(length(orientation_angles),1);
        avg_grating_response = zeros([size(resps,[1 2]),length(orientation_angles)]);
        for iAngle = 1:length(orientation_angles)
            orientation_idx = find(stim_orientation == orientation_angles(iAngle));
            grating_response{iAngle,1} = resps(:,:,orientation_idx);
            avg_grating_response(:,:,iAngle) = mean(grating_response{iAngle,1},3);
        end
        grating_response_wo_error = cell(length(orientation_angles),1);
        avg_grating_response_wo_error = zeros([size(resps,[1 2]),length(orientation_angles)]);
        for iAngle = 1:length(orientation_angles)
            orientation_idx = find(stim_orientation == orientation_angles(iAngle));
            iangle_idx = zeros(size(stim_orientation,1),1);
            iangle_idx(orientation_idx) = 1;
            iangle_idx_wo_error = iangle_idx & stim_index_wo_error;
            grating_response_wo_error{iAngle,1} = resps(:,:,iangle_idx_wo_error);
            avg_grating_response_wo_error(:,:,iAngle) = mean(grating_response_wo_error{iAngle,1},3);
        end
        StaticGratings.probe(options.probe_no).stim_orientation=stim_orientation;
        StaticGratings.probe(options.probe_no).stim_index_wo_error=stim_index_wo_error;
        StaticGratings.probe(options.probe_no).grating_response=grating_response;
        StaticGratings.probe(options.probe_no).grating_response_wo_error=grating_response_wo_error;
        tuning_curve = plot_tuning_curve(grating_response,0,0);
        StaticGratings.probe(options.probe_no).tuning_curve=tuning_curve;



    end
    cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
    save('StaticGrating.mat','StaticGratings')
end

%% Static Gratings Tuning Curve Plotting
cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
load('StaticGratings.mat')
tuning_curve = plot_tuning_curve(StaticGratings.probe(1).grating_response,1,0);

