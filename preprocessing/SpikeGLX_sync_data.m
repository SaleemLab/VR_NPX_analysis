function [systemString, runName, firstG] = SpikeGLX_sync_data(EphysPath, options)

% This function generates and run the necessary command string to run CatGT
% and Tprime for data syncrhonisation

%%% example call for reference %%%
%>runit.sh '-dir=ephysPath -run=run_name -g=ga,gb -t=ta,tb <which streams> [ options ]'


if nargin < 3 % no options struct included, use default settings...
    options = [];
end

%% get g numbers of individual recordings for this session (e.g g0-g5)
path_to_runitsh = '/home/lab/SpikeSorting/CatGT-linux/runit.bat'; % path to CatGT runit.bat file

if isfield(options,'gFileNum')
    gNum = options.gFileNum;
else
gNum_index = strfind(EphysPath, '_g'); % find indices with _g
gNum = str2double(EphysPath(gNum_index(end)+2));
end

gNamesTemp = dir(fullfile(EphysPath,['*','_g','*']));  % this gets all files with a '_g' in their name.

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

% generate the system command as a string to run runit.bat (catGT)
sysString = [' ''-dir=',EphysPath,' -run=',runName,' -g=',gNum];
sysString = [sysString, ' -t=0 -prb=0,1 -zerofillmax=500 -prb_fld'];
sysString = [sysString, ' -xa=0,0 -lffilter=butter,12,0,500'''];
sysString = [sysString, ' -lf -lffilter=butter,12,0,500'''];


% move to parent ehys folder for session to run CatGT 
cd(EphysPath)
% prefix runit.sh location to system command
sysString = ['"',path_to_runitsh,'"', sysString];
systemString = sysString;