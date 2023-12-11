%% read training data


function behaviour = readTrainingData(path,date)

rewardPath = dir(fullfile(path,['*Reward',date,'*.csv']));
lickPerformancePath = dir(fullfile(path,['*Lick_Performance',date,'*.csv']));
wheelLogPath = dir(fullfile(path,['*WheelLog',date,'*.csv']));
reward = readtable([path,rewardPath(1).name]);
lickPerformance = readtable([path,lickPerformancePath(1).name]); 
if ~isempty(wheelLogPath)
    wheelLog = readtable([path,wheelLogPath(1).name]);
end

%% reward type, 1 = passive, 10 = active, non if fail the trial (not licking for active, aborted both the same)
behaviour.rewardType = table2cell(reward(:,"Var1")); 
behaviour.rewardType = cellfun(@(x) x(2:end), behaviour.rewardType, 'UniformOutput',false);
behaviour.rewardType = cellfun(@str2double,behaviour.rewardType);

%%  Reward Delivery Time Based on Eye Frame Count, not Time
behaviour.rewardDeliveryEyeFrameCount = table2array(reward(:,"Var2")); 

%% Track 1 count
behaviour.track1Count = table2array(reward(:,"Var3"));

%% Track 2 count
behaviour.track2Count = table2array(reward(:,"Var4"));

%% Track 1 performance
behaviour.track1Performance = table2array(reward(:,"Var5"));

%% Track 2 performance
behaviour.track2Performance = table2array(reward(:,"Var6"));

%% reward delivery position
behaviour.rewardPosition = table2cell(reward(:,"Var7")); 
behaviour.rewardPosition = cellfun(@(x) x(1:end-1), behaviour.rewardPosition, 'UniformOutput',false);
behaviour.rewardPosition = cellfun(@str2double,behaviour.rewardPosition);

%% LickState
behaviour.lickState = table2cell(lickPerformance(:,"Var1")); 
behaviour.lickState = cellfun(@(x) x(2:end), behaviour.lickState, 'UniformOutput',false);
behaviour.lickState = cellfun(@str2double,behaviour.lickState);

%% Onset Eye Frame Count of Lick
behaviour.lickEyeFrameCount =  table2array(lickPerformance(:,"Var2")); 

%% Track ID respective to lick
behaviour.lickTrackID = table2array(lickPerformance(:,"Var3"));

%% lick track 1 lap count
behaviour.lickTrack1Count = table2array(lickPerformance(:,"Var4"));

%% lick track 2 lap count
behaviour.lickTrack2Count = table2array(lickPerformance(:,"Var5"));

%% reward delivery position
behaviour.lickPosition = table2cell(lickPerformance(:,"Var6")); 
behaviour.lickPosition = cellfun(@(x) x(1:end-1), behaviour.lickPosition, 'UniformOutput',false);
behaviour.lickPosition = cellfun(@str2double,behaviour.lickPosition);

%% wheel log
if ~isempty(wheelLogPath)
%% wheel frame
behaviour.wheelFrame = table2array(wheelLog(:,"Var1"));

%% wheel Time
behaviour.wheelTime = table2array(wheelLog(:,"Var2"));

%% wheel computer time
behaviour.wheelComputerTime = table2cell(wheelLog(:,"Var3"));

%% wheel raw input
behaviour.wheelRawInput = table2cell(wheelLog(:,"Var4"));
behaviour.wheelRawInput = cellfun(@(x) x(8:end), behaviour.wheelRawInput, 'UniformOutput',false);
behaviour.wheelRawInput = cellfun(@str2double,behaviour.wheelRawInput);

%% wheel left lick
behaviour.wheelLickL = table2cell(wheelLog(:,"Var5"));
behaviour.wheelLickL = cellfun(@(x) x(7:end), behaviour.wheelLickL, 'UniformOutput',false);
behaviour.wheelLickL = cellfun(@str2double,behaviour.wheelLickL);

%% wheel right lick
behaviour.wheelLickR = table2cell(wheelLog(:,"Var6"));
behaviour.wheelLickR = cellfun(@(x) x(7:end), behaviour.wheelLickR, 'UniformOutput',false);
behaviour.wheelLickR = cellfun(@str2double,behaviour.wheelLickR);

%% wheel quad state
behaviour.wheelQuadState = table2cell(wheelLog(:,"Var7"));
behaviour.wheelQuadState = cellfun(@(x) x(11:end), behaviour.wheelQuadState, 'UniformOutput',false);
behaviour.wheelQuadState = cellfun(@str2double,behaviour.wheelQuadState);
%% wheel camera frame count
behaviour.wheelEyeFrameCount = table2cell(wheelLog(:,"Var8"));
behaviour.wheelEyeFrameCount = cellfun(@(x) x(15:end), behaviour.wheelEyeFrameCount, 'UniformOutput',false);
behaviour.wheelEyeFrameCount = cellfun(@str2double,behaviour.wheelEyeFrameCount);

%% wheel timestamp of serial string read (time of day, total millisecond)
behaviour.wheelArduinoReadTime = table2cell(wheelLog(:,"Var9"));
behaviour.wheelArduinoReadTime = cellfun(@(x) x(6:end), behaviour.wheelArduinoReadTime, 'UniformOutput',false);
behaviour.wheelArduinoReadTime = cellfun(@str2double,behaviour.wheelArduinoReadTime);

%% wheel position
behaviour.wheelPosition = table2cell(wheelLog(:,"Var10"));
behaviour.wheelPosition = cellfun(@(x) x(10:end), behaviour.wheelPosition, 'UniformOutput',false);
behaviour.wheelPosition = cellfun(@(x) x(1:end-1), behaviour.wheelPosition, 'UniformOutput',false);
behaviour.wheelPosition = cellfun(@str2double,behaviour.wheelPosition);

%% Lick Time
for n = 1:length(behaviour.lickEyeFrameCount)
    wheel_index = find(behaviour.lickEyeFrameCount(n) == behaviour.wheelEyeFrameCount);
    behaviour.lickTime(n) = behaviour.wheelTime(wheel_index(1));
end

end
end