function [clusters] = extract_clusters(options)
%% load Ephys data
% ephys_path = 'D:\Neuropixel_recording\M22069\20221130\kilosort';
% KS_DATAPATH =  'D:\Neuropixel_recording\M22069\20221130\kilosort';
KS_DATAPATH = '/research/DATA/SUBJECTS/M22069/ephys/20221201/kilosort';
% KS_DATAPATH = 'X:/ibn-vision/DATA/SUBJECTS/M22069/ephys/20221201/kilosort';
% KS_DATAPATH = '/research/DATA/SUBJECTS/M22069/ephys/20221130/kilosort';

cd(KS_DATAPATH)
gfileNum = 0; % load g1 file (Track 1 Track 2 running session)
KS_CATGT_FNAME = 'CatGT_M22069_20221201.log';
% KS_CATGT_FNAME = 'CatGT_M22069_20221130.log';
SampleRate = 30000;
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(KS_DATAPATH,gfileNum,KS_CATGT_FNAME,SampleRate)
mean_waveforms = readNPY([KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform
cluster_id = cluster_id + 1;
for id = 1:length(cluster_id)
    peak_channel_waveforms(id,:) = squeeze(mean_waveforms(cluster_id(id),peakChannel(id),:));
end

% imec = Neuropixel.ImecDataset('D:\Neuropixel_recording\M22069\20221130\M22069_20221130_g1\M22069_20221130_g1_imec0\M22069_20221130_g1_t0.imec0.lf');
% imec.inspectAP_timeWindow(timeWindow, ...)

count = 1;
% Find MUA and SUA activity
unit_id = [];
unit_type = [];
SUA_spike_time = [];
SUA_spike_id = [];
SUA_peak_channel_waveforms = [];
SUA_peak_channel = [];
MUA_peak_channel = [];
MUA_peak_channel_waveforms = [];
MUA_spike_time = [];
MUA_spike_id = [];

for unit = 1:length(nominal_KSLabel)
    if peakChannel(unit) > 100 & peakChannel(unit) < 140

        if nominal_KSLabel(unit) == 'good'
            SUA_peak_channel(count) = peakChannel(unit);
            SUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(unit,:);
            SUA_spike_time = [SUA_spike_time; these_spike_times{unit}];
            SUA_spike_id = [SUA_spike_id; count*ones(length(these_spike_times{unit}),1)];
            id_conversion(count,:) = [count cluster_id(unit)];
            count = count + 1;
        elseif nominal_KSLabel(unit) == 'mua'

            MUA_peak_channel(count) = peakChannel(unit);
            MUA_peak_channel_waveforms(count,:) = peak_channel_waveforms(unit,:);
            MUA_spike_time = [MUA_spike_time; these_spike_times{unit}];
            MUA_spike_id = [MUA_spike_id; cluster_id(unit)*ones(length(these_spike_times{unit}),1)];

            count = count + 1;
        end
    end
end


% Sort spike time from start to end
[~,index] = sort(SUA_spike_time);
data.spike_id = SUA_spike_id(index);
data.spike_times = SUA_spike_time(index);
data.peak_channel = SUA_peak_channel;
data.peak_channel_waveforms = SUA_peak_channel_waveforms;
data.id_conversion = id_conversion;

[~,index] = sort(MUA_spike_time);
data.MUA.spike_id = MUA_spike_id(index);
data.MUA.spike_times = MUA_spike_time(index);

% data.MUA.peak_channel = MUA_peak_channel;
% data.MUA.peak_channel_waveforms = MUA_peak_channel_waveforms;

clusters = data;
analysis_folder = '/research/DATA/SUBJECTS/M22069/ephys/20221201/analysis';
% analysis_folder = 'X:\ibn-vision\DATA\SUBJECTS\M22069\ephys\20221201\analysis';
cd(analysis_folder)
% save extracted_clusters_CA3 clusters -v7.3
% 70 220
% save extracted_thalamus_clusters clusters -v7.3
% 0 50
% save extracted_V1_clusters clusters -v7.3
%320 250
