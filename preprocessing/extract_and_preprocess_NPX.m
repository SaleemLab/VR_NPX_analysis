function extract_and_preprocess_NPX(session_info,Stimulus_type)
% Function to import and align bonsai data and spike data for one session
% one stimulus

% Inputs:
% 1. session_info
%
% 2. Stimulus_type
% Specifiy which stimulus type to process -> go through different
% extraction function.

stimulus_name = session_info.probe(1).StimulusName;

% Part 1 Align two NPX probes if there are two probes
if length(session_info.probe) >1
    align_probes_NX1(session_info.probe);
end

% Part 2 Import and align bonsai data to Spikeglx time (always to probe 1)
options = session_info.probe(1);
% Search behaviour files
if contains(Stimulus_type,'Masa2tracks')
    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,sprintf("extracted_behaviour%s.mat",erase(stimulus_name,Stimulus_type))));
else
    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,"extracted_behaviour*.mat"));
end

if isempty(DIR) % skip behaviour extraction if already saved
    if contains(Stimulus_type,'OpenField')
        [Behaviour,~] = import_and_align_Bonsai_OpenField(stimulus_name,options);
    elseif contains(Stimulus_type,'Masa2tracks')
        [Behaviour,Task_info,Peripherals] = import_and_align_Masa_VR_Bonsai(stimulus_name,options);

    elseif contains(Stimulus_type,'Track')
        [Behaviour,Task_info,Peripherals] = import_and_align_Masa_VR_Bonsai(stimulus_name,options);
    elseif contains(Stimulus_type,'Edd')

    else % Else just standard visual stimuli such as Sparse Noise, checkerboard and static grating etc
        [Behaviour,Task_info,Peripherals]  = import_and_align_visual_stimuli_Bonsai(stimulus_name,options);
    end

    if exist(options.ANALYSIS_DATAPATH) == 0
        mkdir(options.ANALYSIS_DATAPATH)
    end

    if contains(Stimulus_type,'Masa2tracks')
        % If Masa2tracks, PRE, RUN and/or POST saved in one folder
        save(fullfile(options.ANALYSIS_DATAPATH,...
            sprintf('extracted_behaviour%s.mat',erase(stimulus_name,Stimulus_type))),'Behaviour')
        save(fullfile(options.ANALYSIS_DATAPATH,...
            sprintf('extracted_task_info%s.mat',erase(stimulus_name,Stimulus_type))),'Task_info')
        save(fullfile(options.ANALYSIS_DATAPATH,...
            sprintf('extracted_peripherals%s.mat',erase(stimulus_name,Stimulus_type))),'Peripherals')
    else
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'),'Task_info')
%           save(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info_without_pd.mat'),'Task_info')
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_peripherals.mat'),'Peripherals')
    end
end

% Part 3 Extract spike data
options = session_info.probe(1); % load behaviour variable first
if contains(Stimulus_type,'Masa2tracks')
    load(fullfile(options.ANALYSIS_DATAPATH,...
        sprintf('extracted_behaviour%s.mat',erase(stimulus_name,Stimulus_type))))
else
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))
end


for nprobe = 1:length(session_info.probe)
    options = session_info.probe(nprobe);
    
    DIR_SORTER = dir(options.SORTER_DATAPATH);
    DIR_KS = dir(options.KS_DATAPATH);
    if ~isempty(DIR_SORTER) % if spike interface sorter folder is present
        [clusters_ks2(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS2','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
        [clusters_ks3(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS3','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
    elseif ~isempty(DIR_KS)% elseif original KS3 folder is present
        [clusters(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','off','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
    end%     [all_clusters chan_config sorted_config] = extract_clusters_NPX(options,'group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
end

% spikes = clusters;
% fields_to_remove = {'spike_count_raw','spike_count_smoothed','zscore_smoothed'};
% clusters = rmfield(clusters,fields_to_remove);


if contains(Stimulus_type,'Masa2tracks')
    % If Masa2tracks, PRE, RUN and/or POST saved in one folder
    if ~isempty(DIR_SORTER) % if spike interface sorter folder is present
        save(fullfile(options.ANALYSIS_DATAPATH,...
            sprintf('extracted_clusters_ks2%s.mat',erase(stimulus_name,Stimulus_type))),'clusters_ks2')
        save(fullfile(options.ANALYSIS_DATAPATH,...
            sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name,Stimulus_type))),'clusters_ks3')
    elseif ~isempty(DIR_KS)% elseif original KS3 folder is present
        save(fullfile(options.ANALYSIS_DATAPATH,...
            sprintf('extracted_clusters%s.mat',erase(stimulus_name,Stimulus_type))),'clusters')
    end
    %             save(fullfile(options.ANALYSIS_DATAPATH,...
    %                 sprintf('extracted_spikes%s.mat',erase(stimulus_name{n},Stimulus_type))),'spikes')
else
    if ~isempty(DIR_SORTER) % if spike interface sorter folder is present
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks2.mat'),'clusters_ks2')
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'),'clusters_ks3')
    elseif ~isempty(DIR_KS)% elseif original KS3 folder is present
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters.mat'),'clusters')
    end
    %             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spikes.mat'),'spikes')
end

end




