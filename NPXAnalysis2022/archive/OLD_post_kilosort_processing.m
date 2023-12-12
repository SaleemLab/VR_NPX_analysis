%% working dirs and filepaths
% 

SUBJECT = 'M22002';
SESSION = '20220215';

%for M22002, 20220214


%for M22002, 20220215
gs_all = {[0 ],[2]};
StimulusName_all = {['SparseNoise'],['StaticGratings']};
AnalysisTimeWindow_all = {[0 0.5], [-0.25 1.25]};
%for M22002, 20220216
% gs_all = {[0 ],[1]};
% StimulusName_all = {['SparseNoise'],['StaticGratings']};
% AnalysisTimeWindow_all = {[0 0.5], [-0.25 1.25]};
%for M22003, 20220221
% gs_all = {[2],[4]}; 
% StimulusName_all ={['Textures2D'],['StaticGratings']};%Sparse Noise was not kilosorted(can't be concatenated by catGT)Drifiting gratings two behaviour files, dirdots also something wrong with the bonsai eye data.
% AnalysisTimeWindow_all = {[-0.25 1.75], [-0.25 1.25]};
%for M22003, 20220222
% gs_all = {[0 1],[2 3 4], [7]};
% StimulusName_all = {['SparseNoise'],['StaticGratings'],['Textures2D']};
% AnalysisTimeWindow_all = {[0 0.5],[-0.25 1.25], [-0.25 1.75]};

% for M22003, 20220223
% gs_all = {[0],[4]};  %{[0],[1],[2],[3],[4]}; 
% StimulusName_all ={['SparseNoise'],['StaticGratings']};%{['SparseNoise'],['DriftingGratings'],['Textures2D'],['DirDots'],['StaticGratings']};
% AnalysisTimeWindow_all = {[0 0.5],[-0.25 1.25]};%{[0 0.5],[-0.25 1.75],[-0.25 1.75],[-0.25 1.75], [-0.25 1.25]};


session_info = {['X:\ibn-vision\DATA\SUBJECTS\',SUBJECT,'\bonsai\',SESSION] ,...
    ['X:\ibn-vision\DATA\SUBJECTS\',SUBJECT,'\ephys\',SESSION,'\']};

% g_number = 0:2; 
% out_start_smp= [14976842, 113176838];
% inp_gamp_smp= [680532,10113816];
% out_zeros_smp = [15000, 15000];
% save([GLXDir,'\catGT_number'],'g_number','out_start_smp','inp_gamp_smp','out_zeros_smp');    
%BonsaiDir = 'X:\ibn-vision\DATA\SUBJECTS\M21192\bonsai\20210817';
%GLXDir = 'X:\ibn-vision\DATA\SUBJECTS\M21192\ephys\20210817\M21192_20210817_H0__g0\M21192_20210817_H0__g0_imec0';

sampleRate = 30000;


for isession = 1:size(session_info,1)
    isession
    BonsaiDir = session_info{isession,1};
%     GLXDir = session_info{isession,2};
%% import GLX async pulse times
% syncPulseUniqueText = 'M22003_20220222_g0_tcat.imec0.ap.SY_384_6_500.txt';
% syncTimes_sglx = getGLXSyncTimes(GLXDir, syncPulseUniqueText);
load([GLXDir,'kilosort3\catGT_number.mat']);
%% import table that maps g_idx to stimulus name and get start and end times for each
% gkey = import_process_gkey(BonsaiDir, GLXDir, syncTimes_sglx, sampleRate);
% 
% %% import and process Bonsai Data (trials, peripherals, eye)
% %StimulusNames = {'SpeedTextures', 'Textures2D', 'Asym', 'SpeedDiscrim'};
% StimulusNames = {'SpeedTuning'};
% 
% [trials, wheel] = import_process_BonsaiData(BonsaiDir, StimulusNames, gkey, syncTimes_sglx);
% 

%% load neural data from Phy/Kilosort
% units = import_ks_data(GLXDir,sampleRate);
% session(isession).units = units;
%good_idx = find(strcmp([units.KSLabel], 'good'));


%% Speed textures analysis
% session(isession).units = AnalyseAsymFlow(trials, units);
% 
% session(isession).units = AnalyseSpeedTextures(trials, units);


end


%%load kilosrot data directly
channel_map=readNPY([GLXDir,'\kilosort3\channel_map.npy']);
channel_positions=readNPY([GLXDir,'\kilosort3\channel_positions.npy']);
similar_templates = readNPY([GLXDir,'\kilosort3\similar_templates.npy']);
templates = readNPY([GLXDir,'\kilosort3\templates.npy']);
templates_ind = readNPY([GLXDir,'\kilosort3\templates_ind.npy']);
whitening_mat = readNPY([GLXDir,'\kilosort3\whitening_mat.npy']);
whitening_mat_inv = readNPY([GLXDir,'\kilosort3\whitening_mat_inv.npy']);
amplitudes = readNPY([GLXDir,'\kilosort3\amplitudes.npy']);
spike_clusters = readNPY([GLXDir,'\kilosort3\spike_clusters.npy']);
spike_templates = readNPY([GLXDir,'\kilosort3\spike_templates.npy']);
spike_times = readNPY([GLXDir,'\kilosort3\spike_times.npy']);
cluster_Amplitude = tdfread([GLXDir,'\kilosort3\cluster_Amplitude.tsv']);
cluster_ContamPct = tdfread([GLXDir,'\kilosort3\cluster_ContamPct.tsv']);
cluster_group = tdfread([GLXDir,'\kilosort3\cluster_group.tsv']);
cluster_KSLabel = tdfread([GLXDir,'\kilosort3\cluster_KSLabel.tsv']);
%fetch the timestamps of each recording in the concatnated file
g_no = length(g_number);
sampleStart =zeros(g_no -1,1); sampleEnd = zeros(g_no - 1,1);
sampleStart(1)= 1; sampleEnd(1) = out_start_smp (1) - 1;
for ig= 1: g_no - 2
    sampleStart(ig + 1) = out_start_smp (ig) + out_zeros_smp (ig);
    sampleEnd(ig+1) = out_start_smp(ig+1)-1;
end
sampleStart(end+1) = out_start_smp (end) + out_zeros_smp(end);
sampleEnd(end+1) = max(spike_times);
spike_times_exp= cell (g_no,1);
spike_times_id = cell(g_no,1);
spike_templates_exp = cell(g_no,1);
spike_clusters_exp = cell(g_no,1);
spike_times_exp_sec = cell (g_no,1);
no_unit = length(cluster_group.cluster_id);
unit_spike = cell (no_unit, g_no);
unit_spike_sec = cell(no_unit, g_no); 
for g_x = 1:g_no
 spike_times_id{g_x} = spike_times >= sampleStart(g_x) & spike_times <= sampleEnd(g_x);
 spike_times_tmp = spike_times(spike_times_id{g_x});
 spike_times_exp{g_x} = spike_times_tmp - sampleStart(g_x) +1;
 spike_times_exp_sec{g_x} = double(spike_times_exp{g_x})/sampleRate;
 spike_templates_exp{g_x} = spike_templates(spike_times_id{g_x});
 spike_clusters_exp{g_x} = spike_clusters(spike_times_id{g_x}); 

 for unit_id = 1:no_unit
     unit_spike{unit_id,g_x} = spike_times_exp{g_x}(spike_templates_exp{g_x} == (unit_id - 1));
     unit_spike_sec{unit_id,g_x} = spike_times_exp_sec{g_x}(spike_templates_exp{g_x} == (unit_id - 1));
 end

 
end


if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
end

for iStim = 1:length(gs_all)

gs = gs_all {iStim};
StimulusName = StimulusName_all{iStim};
options.importMode = 'LF'; % LF or MUA
options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
options.AnalysisTimeWindow = AnalysisTimeWindow_all{iStim};% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])

DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);
% Step 2: Find csv files associated with desired stimulus
[BehaviourDataFiles,EyeDataFiles,TrialParamsFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);



sum_stim =0;
for thisFile = 1:length(BehaviourDataFiles)
    % Set bonsai data paths
    PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{thisFile});
    EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{thisFile});
    TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{thisFile});
    options.PERIPHERALS_DATAPATH = PERIPHERALS_DATAPATH; % full path to BehaviourData bonsai logfile
    options.EYEDATA_DATAPATH = EYEDATA_DATAPATH; % full path to EyeTrackingData bonsai logfile
    options.TRIALDATA_DATAPATH = TRIALDATA_DATAPATH; % full path to TrialData bonsai logfile
    
    % Set ephys data path (temporary, we want to automate this sometime)
    T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,sprintf('%s_%s_g%01d',SUBJECT,SESSION,gs(thisFile)),sprintf('%s_%s_g%01d_imec0',SUBJECT,SESSION,gs(thisFile)));
    T_EPHYS_DATAPATH = fullfile(T_EPHYS_DATAPATH,sprintf('%s_%s_g%01d_t0.imec0',SUBJECT,SESSION,gs(thisFile)));
    options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording

    [tresps,totherData,tstimData{thisFile},~,~,photodiodeData,ttimeVector] = extractAndCollateNP1Data(options);
    stim_window = options.AnalysisTimeWindow;
    photodiode_sglx{thisFile} = photodiodeData.stim_on.sglxTime;
    no_stim = length(photodiode_sglx{thisFile});
    stimData_tmp = photodiodeData;

    for iunit = 1:no_unit
        for iStim = 1: no_stim
            spike_times_tmp = unit_spike_sec{iunit,gs(thisFile)+1};
            stim_spike_times = spike_times_tmp(spike_times_tmp > (stimData_tmp.stim_on.sglxTime(iStim) + stim_window(1)) & ...
            spike_times_tmp < (stimData_tmp.stim_on.sglxTime(iStim) + stim_window(2)));
            
            unit_stim_spike_sec{iunit,iStim+sum_stim} = stim_spike_times - stimData_tmp.stim_on.sglxTime(iStim);
            unit_stim_spike_sglx{iunit,iStim+sum_stim} = stim_spike_times;
        end
    end
    sum_stim = no_stim + sum_stim;
    if thisFile == 1
        resps = tresps;
        otherData = totherData;
    else
        resps = cat(3,resps,tresps);
        otherData = cat(3,otherData,totherData);
    end
end
post_kilosort_data.unit_stim_spike_sec = unit_stim_spike_sec;
post_kilosort_data.unit_stim_spike_sglx = unit_stim_spike_sglx;
post_kilosort_data.resps = resps;
post_kilosort_data.otherData = otherData;
post_kilosort_data.photodiode_sglx = photodiode_sglx;
post_kilosort_data.tstimData = tstimData;
post_kilosort_data.stim_window = stim_window;
post_kilosort_data.unit_spike = unit_spike;
post_kilosort_data.unit_spike_sec = unit_spike_sec;
post_kilosort_data.timeVector = ttimeVector;
post_kilosort_data.BinWidth = options.BinWidth;
post_kilosort_data.spike_templates_exp =spike_templates_exp;
post_kilosort_data.spike_clusters_exp = spike_clusters_exp;
post_kilosort_data.spike_times_exp_sec = spike_times_exp_sec;
post_kilosort_data.spike_times_exp = spike_times_exp;
post_kilosort_data.out_start_smp = out_start_smp;
post_kilosort_data.inp_gamp_smp= inp_gamp_smp;
post_kilosort_data.out_zeros_smp = out_zeros_smp;
post_kilosort_data.gs= gs;
post_kilosort_data.StimulusName = StimulusName;
post_kilosort_data.importMode = options.importMode;
post_kilosort_data.channel_map=channel_map;
post_kilosort_data.channel_positions=channel_positions;
post_kilosort_data.similar_templates =similar_templates;
post_kilosort_data.templates = templates;
post_kilosort_data.templates_ind = templates_ind;
post_kilosort_data.whitening_mat = whitening_mat;
post_kilosort_data.whitening_mat_inv = whitening_mat_inv;
post_kilosort_data.amplitudes = amplitudes;
post_kilosort_data.spike_clusters = spike_clusters;
post_kilosort_data.spike_templates = spike_templates;
post_kilosort_data.spike_times = spike_times;
post_kilosort_data.cluster_Amplitude = cluster_Amplitude;
post_kilosort_data.cluster_ContamPct = cluster_ContamPct;
post_kilosort_data.cluster_group = cluster_group;
post_kilosort_data.cluster_KSLabel = cluster_KSLabel;
post_kilosort_data.SUBJECT = SUBJECT;
post_kilosort_data.SESSION = SESSION;
post_kilosort_data.options = options;
save(['C:\Users\adam.tong\Documents\CorticalCommunication','\',SUBJECT,'_',SESSION,'_',StimulusName,'_post_kilosort_data'],'post_kilosort_data','-v7.3'); 
end
return
% Temporary example of how to calculate sparsenoise RFs
% options.figName = sprintf('%s :: %s (file %01d)',    thisAnimal, thisSession,  thisFileNum);
% initMap = sparseNoiseAnalysis(stim_matrix,lfp_data,wheel_data,eye_data,options);



% catGT_table= readCatGTlog('C:\Users\adam.tong\Documents\CorticalCommunication\CatGT.log');

