function [systemString, runName, firstG] = preprocess_and_save_LFP(options)

% This function generates and run the necessary command string to run CatGT
% and Tprime for data syncrhonisation

%%% example call for reference %%%
%>runit.sh '-dir=ephysPath -run=run_name -g=ga,gb -t=ta,tb <which streams> [ options ]'

%% Find if LF files already there
DIR = dir(fullfile(options.EPHYS_DATAPATH,'*tcat*'));
% cd(options.ANALYSIS_DATAPATH)
if ~isempty(DIR)
    return
end

EphysPath = cd(fullfile(options.EPHYS_DATAPATH,'..','..'));

%% get g numbers of individual recordings for this session (e.g g0-g5)
DIR = dir('C:/Users/masah/Documents/GitHub/VR_NPX_analysis/preprocessing/CatGT-win');
if ~isempty(DIR)
    path_to_runitsh = 'C:/Users/masah/Documents/GitHub/VR_NPX_analysis/preprocessing/CatGT-win/runit.bat'; % path to CatGT runit.bat file
else
    DIR = dir('C:/Users/masahiro.takigawa/Documents/GitHub/VR_NPX_analysis/preprocessing/CatGT-win');
    if ~isempty(DIR)
        path_to_runitsh = 'C:/Users/masahiro.takigawa/Documents/GitHub/VR_NPX_analysis/preprocessing/CatGT-win/runit.bat'; % path to CatGT runit.bat file
    else
        disp('CatGT folder not found!')
    end
end

if isfield(options,'gFileNum')
%     gNum = options.gFileNum;
    gNum= num2str(options.gFileNum);
else
    gNum_index = strfind(options.EPHYS_DATAPATH, '_g'); % find indices with _g
    gNum = options.EPHYS_DATAPATH(gNum_index(end)+2);
end

gNamesTemp = dir(fullfile(EphysPath,['*','_g',gNum,'*']));  % this gets all files with a '_g' in their name.

% get run name (basename of session, without gindices, e.g. M22030_20220615, older sessions may also have _h0 for example.)
runName_temp = strsplit(gNamesTemp(1).name,'_');
runName = runName_temp{1};
for istr = 2:numel(runName_temp)-1
    runName = [runName, '_', runName_temp{istr}];
end

% % generate the system command as a string to run runit.sh (catGT)
% sysString = [' ''-dir=',EphysPath,' -run=',runName,' -g=',gNum];
% sysString = [sysString, ' -t=0 -prb=0 -zerofillmax=500 -prb_fld'];
% sysString = [sysString, ' -ap -apfilter=butter,12,300,0 -gblcar -gfix=0,0.1,0.02 -xd=2,0,-1,6,500'''];

% EphysPath = strrep(EphysPath, '\', '/');

% EphysPath = fullfile('\\live.rd.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS',runName_temp{1},runName_temp{2})

% generate the system command as a string to run runit.bat (catGT)
% sysString = [' ''-dir=',fullfile('\\live.rd.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision',extractAfter(EphysPath,'ibn-vision')),' -run=',runName,' -g=',gNum];

% For windows it is important to ensure 'space' before each '-'
sysString = [' -dir=',EphysPath,' -run=',runName,' -g=',gNum];
sysString = [sysString, ' -t=0 -prb=0,1 -zerofillmax=500 -prb_fld'];
% sysString = [sysString, ' -xa=0,0 -lffilter=butter,12,0,500'''];
sysString = [sysString, ' -gblcar -lf -lffilter=butter,12,0,500'];


% move to parent ehys folder for session to run CatGT 
% cd(options.EPHYS_DATAPATH)
cd(DIR(1).folder)
% prefix runit.sh location to system command
sysString = ['"',path_to_runitsh,'"', sysString];
systemString = sysString;

% run CatGT
disp('Running CatGT to preprocess LFP...')
catGTstart = tic;
[status,cmdout] = system(sysString);
disp('CatGT complete')
fprintf('Time for CatGT to run: %s\n', duration([0, 0, toc(catGTstart)]))