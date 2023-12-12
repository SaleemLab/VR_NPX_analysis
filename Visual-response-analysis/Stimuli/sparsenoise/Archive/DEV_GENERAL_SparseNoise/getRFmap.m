%% Tomaso Muzzu & Mai Morimoto - UCL - 04/09/2019


function ES = getRFmap(mouseName,iexp,plotSingleUnits)

[~,hostName] = system('hostname'); hostName = hostName(1:end-1);
if ~strcmp(hostName, 'saleem12')
    addpath('X:\CODE\STABLE\OpenEphys_analysis_tools');
    SubjectDir = ['X:' filesep 'DATA' filesep 'SUBJECTS' filesep];
else
    addpath('X:\ibn-vision\CODE\STABLE\OpenEphys_analysis_tools');
    SubjectDir = ['X:\ibn-vision' filesep 'DATA' filesep 'SUBJECTS' filesep];
end
addpath([cd filesep 'SparseNoiseSyncR']);

% SubjectDir = ['X:' filesep 'DATA' filesep 'SUBJECTS' filesep];

% check if sparse noise stimulus has already synced with ephys
ProcessedFiles = dir([SubjectDir mouseName filesep 'Processed']);
f_n = 1; f_n_i = 0;
while f_n<=length(ProcessedFiles)
    if strfind(ProcessedFiles(f_n).name,'SN')
        d_i = strfind(ProcessedFiles(f_n).name,['20' iexp(1:2)]);
        if isempty(d_i)
            d_i = strfind(ProcessedFiles(f_n).name,['ks_'])+1;
        end
        if strfind(ProcessedFiles(f_n).name(d_i+5:d_i+6),iexp(3:4))
            if strfind(ProcessedFiles(f_n).name(d_i+8:d_i+9),iexp(5:6))
                if f_n_i == 0 
                    f_n_i = f_n;
                else
                    f_n_i = [f_n_i f_n];
                end
            end
        end
    end
    f_n = f_n + 1;
end

if f_n_i(1) == 0
    % do you want to sync the file now?
    SyncNow = input('Recording does not exist or is not synced yet. Would you like to sync it now? 1-yes 0-no\n');
    if SyncNow
       ES = SparseNoiseSyncR(mouseName);
    end
else
    if length(f_n_i)> 1
        % multiple sparse noise signals were used on the day
        fprintf('\n\nMultiple SN recordings found: \n')
        for r = 1:length(f_n_i)
            fprintf([num2str(r) ' = ' ProcessedFiles(f_n_i(r)).name '\n'])
        end
        Rec_OI = input('Which one do you want to look at? \n');
        load([ProcessedFiles(f_n_i(Rec_OI)).folder filesep ProcessedFiles(f_n_i(Rec_OI)).name]);
    else
        load([ProcessedFiles(f_n_i).folder filesep ProcessedFiles(f_n_i).name]);
    end
end




m = ReverseCorrelationSTA(ES,plotSingleUnits);

% m = ReverseCorrelationSTC(ES,m,plotSingleUnits);

m = ForwardCorrelationOFFMap(ES,m,plotSingleUnits);

m = ForwardCorrelationONMap(ES,m,plotSingleUnits);

% m = MultiMeasureRF_MUA(ES,m);

% m = MultiMeasureRF_MUA_smooth(ES,m);

ES.m = m;
% if plotSingleUnits
%     temp_sel = input('Would you like to look at specific units?\n')
%     if temp_sel == 0
%         selected = 1:length(ES.StimSpiketimes);
%     else
%         selected = temp_sel;
%     end
%     % selected=[1 7 15 16 17 22 27 28];
%     m = MultiMeasureRF(ES,m,selected)
% end
% EOF
end




