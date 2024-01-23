function extract_PSD_profile(session_info,Stimulus_type)

% program to extract PSD profile from one session


% Related:
% extract_PSD_profile_batch(experiment_info,Stimulus_type)


PSD = [];
power = [];
best_channels = [];
for nprobe = 1:length(session_info.probe)
%     session_info.probe(nprobe).importMode = 'LF';
    options = session_info.probe(nprobe);

    column = 1;
    if nprobe ~= 1
        % Information about AP band probe 1 sample size to align probe
        % 2 LFP traces.
        session_info.probe(1).importMode = 'KS'; % need total sample number at AP band for later probe aligned LFP
        [~, imecMeta, ~, ~] = extract_NPX_channel_config(session_info.probe(1),column);
        [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column,'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp);
    else
        [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column);
    end

    [PSD{nprobe} power{nprobe} best_channels(nprobe)] = calculate_channel_PSD(raw_LFP,SR,sorted_config,options,'plot_option',1)

    %         title(sprintf('%s %s PSD profile probe %i',options.SUBJECT,options.SESSION,nprobe))
    %         filename = sprintf('%s %s PSD profile probe %i.pdf',options.SUBJECT,options.SESSION,nprobe)
    %         saveas(gcf,filename)
    %         filename = sprintf('%s %s PSD profile probe %i.fig',options.SUBJECT,options.SESSION,nprobe)
    %         saveas(gcf,filename)
end
% For now saving PSD and best channels in analysis folder rather than in a
% specific stimuli folder because of the same session (day), the profile
% should be usually more or less consistent

save(fullfile(erase(options.ANALYSIS_DATAPATH,['\',Stimulus_type]),"best_channels.mat"),'best_channels')
save(fullfile(erase(options.ANALYSIS_DATAPATH,['\',Stimulus_type]),"extracted_PSD.mat"),'PSD','power')