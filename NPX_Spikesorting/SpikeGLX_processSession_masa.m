%% Script to pre-process ephys data recorded using neuropixels/spike glx

% - Written by EH, first push (xx/05/2022)
% - Edited by FR, 09/05/2022
% - EH 26/05/22. Save kilosort rez variable as a .mat file
% - EH 16/06/22. Fixed g0,gN input arg generation for catGT to work with
% gNums > 9. Prevous version didn't pick up these values. Also prevented
% deletion of temp_wh.dat as needed for phy to take account of drift during
% manual curation.
% - EH 08/07/22 - new function to generate CatGT command to make more
% flexible for future use where we may want to change parameters. EH fixed a
% related bug (22/07/22)
% 


%% To do:
% - this needs converting to a function. 
% Subject/session args 
% 0) check cmd out messages. If fails, report a message
% 1) Need to work on integration with manual curation spike sorting 
% -> i.e. after manual curation, the mean_waveforms/quality_metrics would 
% need to be re-run, and the cluster_table updated. A separate script
% should be written for this
% 2) There's something weird with mean_waveforms, cluster_labels,
% noise_labels -> they seem to cover the original un-merged clusters
% (kilosort automerges some clusters as the final step). I guess these
% files are accessing the original un-merged space? Might need to update
% kilosort outputs? 

%% 0. Choose session and setup directories

function SpikeGLX_processSession_masa(SUBJECT,SESSION)

% SUBJECT = 'M22022';
% SESSION = '20220504';

totalRunTimeStart = tic;

isbeast = true;

addpath(genpath('/home/lab/NPXAnalysis')); % path to where this script is -> change to NPX analysis github repo.
% tempSpikeSortingPath = fullfile('/home/masa/NPX_analysis_code/temp_spikesorting', SUBJECT, SESSION); % where to run pipeline locally.
tempSpikeSortingPath = fullfile('/home/lab/temp_spikesorting', SUBJECT, SESSION); % where to run pipeline locally.
addpath(genpath('/home/lab/Kilosort-main')) % path to main kilosort matlab folder
pathToYourConfigFile = '/home/lab/Kilosort-main/configFiles'; % basic kilosort config file (from original github)
addpath(genpath('/home/lab/SpikeSorting/npy-matlab-master')) % for converting to Phy
path_to_runitsh = '/home/lab/SpikeSorting/CatGT-linux/runit.sh'; % path to CatGT runit.sh file

% rootpath for server
if isbeast
    ROOTPATH = '/research/';
    %ROOTPATH = '/mnt/pfs09/'; % temp until lab gets 'research' access
else
    ROOTPATH = 'X:\ibn-vision\';
end

% Get path to ephys session main parent directory on the server
% Report files found in the folder to the command window
ephysPath = fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'ephys',SESSION);

FilesFound = dir(ephysPath);
nFiles = numel(FilesFound)-2;
disp([num2str(nFiles), ' files found in ' ephysPath, ':'])
dir(ephysPath)

if any(strcmp({FilesFound.name}, 'kilosort'))
    error('Kilsoort directory already present on server for this session, please move or delete\n%s', fullfile(FilesFound(1).folder, 'kilosort'))
end

%% 1. copy files from server to /home/ (fast SSD).
startCopyTime = tic;
disp(['Copying files to local drive: ', tempSpikeSortingPath])
status = copyfile(ephysPath, tempSpikeSortingPath)
fprintf('Time to copy files to local drive: %s\n', duration([0, 0, toc(startCopyTime)]))


% create 'kilosort' folder in local ephys directory
localEphysPath = fullfile(tempSpikeSortingPath);
cd(localEphysPath); 
mkdir(fullfile(localEphysPath, 'kilosort'));
kilosortDir = fullfile(localEphysPath, 'kilosort');


%% 2. Catenate files using CatGT
% This generates the necessary parameter string to run CatGT for this
% subject and session.
% It then runs CatGT using a simple system() command
% Finally it moves the output *.tcat* and CatGT.log files to the kilosort
% directory.

% generate cmd line string to run catGT
[sysString, runName, firstG] = generate_CatGT_command(localEphysPath, path_to_runitsh);

% move to parent ehys folder for session to run CatGT 
cd(localEphysPath)
% run CatGT
disp('Running CatGT...')
catGTstart = tic;
[status,cmdout] = system(sysString);
disp('CatGT complete')
fprintf('Time for CatGT to run: %s\n', duration([0, 0, toc(catGTstart)]))

% find catGT.log file, rename and move to kilosort dir
cd(localEphysPath);
newCatGTLogFileName = ['CatGT', '_', SUBJECT, '_', SESSION, '.log'];
movefile('CatGT.log', [kilosortDir, '/', newCatGTLogFileName]); % change to full file?

% get catenated ephys files and move to kilosort dir
% CatGT stores catenated 'tcat' files in g0 by default (or just the
% first?), so we look for them here using the spikeGlx/imec naming 
% convention of recordings
%cd([localEphysPath, '/', runName,'_g0', '/', runName,'_g0_imec0']) 
% this might change to:
%cd([localEphysPath, '/', runName, '_', firstG{end}, '/', runName, '_', firstG{end}, '_imec0']) % change to fullfile?
tp = [fullfile(localEphysPath, runName), '_g',num2str(firstG)];
tp = [fullfile(tp, runName), '_g', num2str(firstG), '_imec0'];
cd(tp)
tcatFiles = dir('*_tcat.*');

% move files to new kilosort directory
for ifile = 1:numel(tcatFiles)
    movefile(tcatFiles(ifile).name, fullfile(kilosortDir,tcatFiles(ifile).name))
end

%% 3. Run Kilosort

% set current dir as kilosort dir
cd(kilosortDir)

% set necessary options for kilosort
% generate channel map for kilosort
metaName = dir('*tcat.imec0.ap.meta');
binName = dir('*tcat.imec0.ap.bin');

saveDir = kilosortDir;
% modified function points directly to .meta file and saves in specified
% directory (kilosort dir).
SGLXMetaToCoords_ChannelMap_EHversion(metaName,saveDir)

rootZ = kilosortDir; % the raw data binary file is in this folder
rootH = kilosortDir; % path to temporary binary file (same size as data, should be on fast SSD)

% get the basic config settings from this default file
run(fullfile(pathToYourConfigFile, 'configFile384.m'))

% location of temp file
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD

fprintf('Looking for data inside %s \n', rootZ)

% main parameter changes from Kilosort2 to v2.5
ops.trange    = [0 Inf]; % time range to sort
ops.NchanTOT  = 385; % total number of channels in your recording
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
% main parameter changes from Kilosort2.5 to v3.0
ops.Th       = [9 9];

% get channel map
fs = dir(fullfile(rootZ, '*kilosortChanMap.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = dir(fullfile(rootZ, '*tcat.imec0.ap.bin')); % get tcat binary file
ops.fbinary = fullfile(rootZ, fs(1).name);

% run kilosort
disp('running kilosort3...')
ks_start_time = tic;
rez                = preprocessDataSub(ops);
rez                = datashift2(rez, 1);
[rez, st3, tF]     = extract_spikes(rez);
rez                = template_learning(rez, tF, st3);
[rez, st3, tF]     = trackAndSort(rez);
rez                = final_clustering(rez, tF, st3);
rez                = find_merges(rez, 1);
rezToPhy2(rez, rootZ);
disp('kilosort complete!')
fprintf('Time for kilsoort to run: %s\n', duration([0, 0, toc(ks_start_time)]))
% save kilosort figures
cd(kilosortDir)

% get all figures open
h = findall(0,'Type', 'figure');
% save the last 3 figures (generated by kilosort)
savefig(h(end),[SUBJECT, '_', SESSION, '_estDriftTraces']);
savefig(h(end-1),[SUBJECT, '_', SESSION, '_driftMap']);
savefig(h(end-2),[SUBJECT, '_', SESSION, '_unitInfo']);

close all
save([SUBJECT, '_' SESSION, '_ks_rez.mat'], 'rez', '-v7.3', '-nocompression'); % save kilosort rez file.
% delete temp_wh.dat file - no longer deleting as needed for manual
% curation!
%delete(fullfile(rootH, 'temp_wh.dat'))

%% 4) Run allen modules 
% - ks post-processing (remove double counted spikes)
% - label noise clusters
% - mean_waveform 
% - quality_metrics
cd(kilosortDir)
% postgrad139
%cmd = 'cd C:\Users\edward.hrrocks\Documents\Github\JC_ecephys_spike_sorting && pipenv run C:\Users\edward.horrocks\Documents\GitHub\JC_ecephys_spike_sorting\.venv\Scripts\python C:\Users\edward.horrocks\Documents\GitHub\JC_ecephys_spike_sorting\ecephys_spike_sorting\scripts\sglx_runAllenPipeline_eh.py -sub M22008 -sesh 20220408 -ksdir C:\Users\edward.horrocks\Documents\Github\JC_ecephys_spike_sorting -df X:\ibn-vision\DATA\SUBJECTS\M22008\ephys\20220408\kilosort\M22008_20220408_g0_tcat.imec0.ap.bin'
ephysPath = fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'ephys', SESSION);

% beast
% cmd = 'cd /home/edd/Kilosort-main/JC_ecephys_spike_sorting && pipenv run /home/edd/Kilosort-main/JC_ecephys_spike_sorting/.venv/Scripts/python /home/edd/Kilosort-main/JC_ecephys_spike_sorting/ecephys_spike_sorting/scripts/sglx_runAllenPipeline_eh.py';
% cmd = [cmd, ' -sub ', SUBJECT, ' -sesh ', SESSION, ' -ksdir ', kilosortDir, ' -df ', fullfile(binName.folder, binName.name)]
% cmd = ['source /home/edd/anaconda3/etc/profile.d/conda.sh && conda activate AllenPipelineEnv && '];

% edd user
% cmd = ['bash /home/edd/anaconda3/bin/activate AllenPipelineEnv && '];
% cmd = [cmd, 'cd /home/edd/SpikeSorting/ecephys_spike_sorting-master/ && pipenv run python '];
% cmd = [cmd, '"/home/edd/SpikeSorting/ecephys_spike_sorting-master/ecephys_spike_sorting/scripts/sglx_runAllenPipeline_eh.py" '];
%cmd = [cmd, ' -sub M22008 -sesh 20220408 -ksdir ', kilosortDir, ' -df ', fullfile(binName.folder, binName.name)]

% lab user
%cmd = ['bash /home/lab/anaconda3/bin/activate allenPipelineEnv && '];
cmd = ['bash /opt/anaconda3/bin/activate allenPipelineEnv && '];
cmd = [cmd, 'cd /home/lab/SpikeSorting/ecephys_spike_sorting-master/ && pipenv run python '];
cmd = [cmd, '"/home/lab/SpikeSorting/ecephys_spike_sorting-master/ecephys_spike_sorting/scripts/sglx_runAllenPipeline_eh.py" '];
cmd = [cmd, ' -sub ', SUBJECT, ' -sesh ', SESSION, ' -ksdir ', kilosortDir, ' -df ', fullfile(binName.folder, binName.name)]


% need to make sure this saves into the main kilosort directory
%cmd = ['/home/edd/Kilosort-main/JC_ecephys_spike_sorting.venv\Scripts\python.exe'];

allen_pipeline_start = tic;
disp('Running allen pipeline...')
[a, b] = system(cmd)
fprintf('Time for allen pipeline modules to run: %s\n', duration([0, 0, toc(allen_pipeline_start)]))


%% 5) Collate information about clusters into a table and save as.mat file
cd(kilosortDir)
cluster_table =  process_kilosort_and_allen_data(kilosortDir);
save('cluster_table.mat', 'cluster_table') % ~100 MB.


% 6) move relevent files back to server
disp('Copying files back to server...')
copy2serverTime = tic;
copyfile(kilosortDir,fullfile(ephysPath,'kilosort'))
fprintf('Time to copy files to server: %s\n', duration([0, 0, toc(copy2serverTime)]))

fprintf('Total pipeline run time: %s\n', duration([0, 0, toc(totalRunTimeStart)]))
                                               
