function [lfpAvg,csd,power,best_channels] = determine_best_channels(options,nprobe)
best_channels = [];
power = [];
lfpAvg = [];
csd = [];

options.BinWidth = 1/1250;
options.importMode = 'LF'; % LF or MUA or KS
% options.importMode = 'LF'; % LF or MUA or KS
options.AnalysisTimeWindow = [-0.1 0.5];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
%             options.ks_unitType = 'good'; % 'mua', 'good' or ''
options.paradigm = 'photodiode_timestamp';

[resps,otherData,stimData,~,wheelData,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
timestamps = linspace(options.AnalysisTimeWindow(1),options.AnalysisTimeWindow(2),size(resps,2));
column = 1;
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);

[csd.all, lfpAvg.all ] = perievent_CSD_LFP_amplitude_phase(permute(resps(sorted_config.Channel,:,:), [2, 1, 3]), timestamps,[],'twin',[0.1 0.5]);


load(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis',"best_channels.mat"))
load(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis',"extracted_PSD.mat"))

lfpAvg.filter_type = {'all'};
lfpAvg.event_group = {'Checkerboard'};
lfpAvg.probe_no = nprobe;
cd(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
%             plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power{nprobe},chan_config,sorted_config,best_channels{nprobe})

plot_perievent_CSD_LFP(lfpAvg,csd,power{nprobe},chan_config,sorted_config,best_channels{nprobe},options)
