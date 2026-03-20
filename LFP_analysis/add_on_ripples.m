function add_on_ripples(session_info,stimulus_name,best_channels)

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
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP')
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name,'Masa2tracks'))))
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('behavioural_state_merged%s.mat',erase(stimulus_name,'Masa2tracks'))))
    load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name,'Masa2tracks'))))

else
    % load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
    load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'),'behavioural_state_merged')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'),'PSD','power')% save PSD for the sleep session
    % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));% clusters for MUA
    clusters = clusters_ks4;
    load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name,'Chronic'))),'session_clusters'); % Session clusters for SUA
end


clear ripples_temp
% for nprobe = 1:length(session_info.probe)
for nprobe =1:2
    options = session_info.probe(nprobe);
    probe_no = session_info.probe(nprobe).probe_id + 1;
    options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

    tvec = LFP(probe_no).tvec;
    if contains(stimulus_name,'Sleep')
        speed = double(session_clusters.mobility_thresholded{1});
        Behaviour.mobility_zscore=session_clusters.mobility_zscore{1};
    else
        speed = Behaviour.speed;
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
    end

    for best_channel = 1:size(LFP(nprobe).best_HPC,1)
        [ripples_temp] = FindRipples_masa(LFP(nprobe).best_HPC(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'minDuration',30,'durations',[30 200],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
            'noise',[],'passband',[125 300],'thresholds',[2 5],'show','on','best_channel',best_channel);
        if length(ripples(nprobe).onset) == length(ripples_temp.onset)
            if sum(ripples(nprobe).onset == ripples_temp.onset) == length(ripples(nprobe).onset)
                ripples(nprobe).best_channel= ripples_temp.best_channel;
                continue
            end
        end
    end

end

if ~isempty(ripples(1).best_channel) & ~isempty(ripples(2).best_channel)
    if contains(stimulus_name,'Masa2tracks')
        save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'ripples')
        %     save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name,'Masa2tracks'))),'LFP','-v7.3')
    else
        save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'),'ripples')
        %     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')
    end
else
    disp([options.ANALYSIS_DATAPATH,' ripples best channel not uploaded!!!!'])
end
% save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name,'Masa2tracks'))),'spindles')

% save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'),'spindles')
% save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'),'slow_waves')
% save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP','-v7.3')
% save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP','-v7.3')

close all