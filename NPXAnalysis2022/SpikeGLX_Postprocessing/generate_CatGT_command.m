function [systemString, runName, firstG] = generate_CatGT_command(localEphysPath, path_to_runitsh, options)

% This function generates the necessary command string to run CatGT for this
% subject and session.

%%% example call for reference %%%
%>runit.sh '-dir=ephysPath -run=run_name -g=ga,gb -t=ta,tb <which streams> [ options ]'


if nargin < 3 % no options struct included, use default settings...
    options = [];
end

%% get g numbers of individual recordings for this session (e.g g0-g5)
gNamesTemp = dir(fullfile(localEphysPath,['*','_g','*']));  % this gets all files with a '_g' in their name.

% check these are actually gN filenames (i.e. not kilosort etc).
for ig = 1:numel(gNamesTemp)    
    gNames{ig}=gNamesTemp(ig).name;
    tempG = strsplit(gNames{ig},'_'); % split the file name by '_'
    if ~isempty(regexp(tempG{end}, 'g\d')) % check if this file has a g\d part in it's name (i.e. is a recording)
        gVal_temp = [regexp(gNames{ig},'g\d*','Match')]; % this finds the g** part of the string.
        gNum(ig) = str2num(cell2mat(erase(gVal_temp,'g'))); % this removes the g and converts to a double
    else
        gNum(ig) = [];
    end
end    

gNum = sort(gNum); % sort gNs by ascending order
    
firstG = min(gNum);
lastG = max(gNum);
gRange = [num2str(firstG), ',' , num2str(lastG)]; % generate g string (e.g. 'G0,G12')

% get run name (basename of session, without gindices, e.g. M22030_20220615, older sessions may also have _h0 for example.)
runName_temp = strsplit(gNamesTemp(1).name,'_');
runName = runName_temp{1};
for istr = 2:numel(runName_temp)-1
    runName = [runName, '_', runName_temp{istr}];
end

% generate the system command as a string to run runit.sh (catGT)
sysString = [' ''-dir=',localEphysPath,' -run=',runName,' -g=',gRange];
sysString = [sysString, ' -t=0 -prb=0 -zerofillmax=500 -prb_fld'];
sysString = [sysString, ' -ap -aphipass=300 -gbldmx -gfix=0,0.1,0.02 -SY=0,-1,6,500'''];

% move to parent ehys folder for session to run CatGT 
cd(localEphysPath)
% prefix runit.sh location to system command
sysString = ['"',path_to_runitsh,'"', sysString];
systemString = sysString;