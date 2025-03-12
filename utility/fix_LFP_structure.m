function fix_LFP_structure(session_info,stimulus_name,best_channels)

options = session_info.probe(1);
% load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
if isempty(DIR)
    disp('No extracted clusters')
    return
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
if contains(stimulus_name,'Masa2tracks')
    % load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP')
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name,'Masa2tracks'))))
    % load(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state_merged%s.mat',erase(stimulus_name,'Masa2tracks'))))
    % load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name,'Masa2tracks'))))

else
    % load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')
    % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
    % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))
    % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
    % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
    % load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'),'behavioural_state_merged')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
    % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));% clusters for MUA
    % clusters = clusters_ks4;
    % load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name,'Chronic'))),'session_clusters'); % Session clusters for SUA
end

% for nprobe = 1:length(session_info.probe)
for nprobe = 1:2
    options = session_info.probe(nprobe);

    options.importMode = 'KS'
    [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);
    % scatter(chan_config.Ks_xcoord,chan_config.Ks_ycoord)
    probe_no = session_info.probe(nprobe).probe_id + 1;

    for nchannel = 1:length(PSD{nprobe})
        power{nprobe}(nchannel,:) = PSD{nprobe}(nchannel).mean_power;
    end

    for nshank = 1:4
        % channel_id = 48;%40
        % nshank = 1;%2
        channels_selected = [PSD{nprobe}.channel];
        index = find(channels_selected == LFP(probe_no).best_HPC_channel(nshank));
        if nshank ~= PSD{nprobe}(index).shank
            disp('shank id of LFP.best_HPC_shank_id not consistent')
            continue
        end

        % LFP(probe_no).best_HPC(nshank,:) =  raw_LFP(channel_id,:);
        % LFP(probe_no).best_HPC_channel(nshank) = PSD{nprobe}(IB(channel_id)).channel;
        LFP(probe_no).best_HPC_depth(nshank) = PSD{nprobe}(index).ycoord;
        LFP(probe_no).best_HPC_xcoord(nshank) = PSD{nprobe}(index).xcoord;
        % LFP(probe_no).best_HPC_power(nshank,:) =  power{nprobe}(index,:);
        % LFP(probe_no).best_HPC_shank_id(nshank) =  nshank;
    end
    if nprobe ==2
        LFP(probe_no).best_HPC(4,:) = [];
        LFP(probe_no).best_HPC_channel(4) = [];
        LFP(probe_no).best_HPC_depth(4) = [];
        LFP(probe_no).best_HPC_xcoord(4) = [];
        LFP(probe_no).best_HPC_power(4,:) =  [];
        LFP(probe_no).best_HPC_shank_id(4) =  [];

        % LFP(probe_no).best_CA1(nshank,:) =  raw_LFP(channel_id,:);
        % LFP(probe_no).best_CA1_channel(nshank) = PSD{nprobe}(IB(channel_id)).channel;
        % LFP(probe_no).best_CA1_depth(nshank) = PSD{nprobe}(IB(channel_id)).ycoord;
        % LFP(probe_no).best_CA1_xcoord(nshank) = PSD{nprobe}(IB(channel_id)).xcoord;
        % LFP(probe_no).best_CA1_power(nshank,:) =  power{nprobe}(IB(channel_id),:);
        % LFP(probe_no).best_CA1_shank_id(nshank) =  nshank;
        %
    end
end

if contains(stimulus_name,'Masa2tracks')
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
    save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP','-v7.3')
else
    % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')
end
