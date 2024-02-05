function extract_PSD_profile_batch(experiment_info,Stimulus_type)

% program to extract PSD profile from all sessions (specified by experiment_info)
% suitable for batch preprocessing

% Related:
% extract_PSD_profile(session_info,Stimulus_type)


for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    
    if isempty(session_info)
        continue
    end

    PSD = [];
    power = [];
    best_channels = [];
    for nprobe = 1:length(session_info.probe)
%         session_info.probe(nprobe).importMode = 'LF';
        options = session_info.probe(nprobe);

        column = 1;
        if nprobe ~= 1
            % Information about AP band probe 1 sample size to align probe
            % 2 LFP traces.
            session_info.probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
            [~, imecMeta, chan_config, ~] = extract_NPX_channel_config(session_info.probe(1),[]);
            [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX(options,[],'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp,'selected_channels',chan_config.Channel);
        else
            [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX(options,[]);
        end
        
        [PSD{nprobe} power{nprobe} best_channels{nprobe}] = calculate_channel_PSD(raw_LFP,SR,chan_config,options,'plot_option',1)
%         [PSD{nprobe} power{nprobe} best_channels{nprobe}] = calculate_channel_PSD(raw_LFP,SR,sorted_config,options,'plot_option',1)

        %         title(sprintf('%s %s PSD profile probe %i',options.SUBJECT,options.SESSION,nprobe))
        %         filename = sprintf('%s %s PSD profile probe %i.pdf',options.SUBJECT,options.SESSION,nprobe)
        %         saveas(gcf,filename)
        %         filename = sprintf('%s %s PSD profile probe %i.fig',options.SUBJECT,options.SESSION,nprobe)
        %         saveas(gcf,filename)
    end
    

    save_dir = fullfile(options.ANALYSIS_DATAPATH, '..'); % find parent dir

    save(fullfile(save_dir,"best_channels.mat"),'best_channels')
    save(fullfile(save_dir,"extracted_PSD.mat"),'PSD','power')
end