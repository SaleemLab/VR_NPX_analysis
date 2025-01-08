clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\Testing\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
%% LFP preprocess -> Catgt

SUBJECT = 'M24072';
all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys','20*'));

for ndate = 1:length(all_DIR)
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))

    session_DIR= dir('20*');

    Error_sessions = [];
    Error_idx = [];
    Error_timestamp=[];
    Error_acquisition = [];

    for nsession = 1:length(session_DIR)
        cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
        DIR = dir('*_g*'); % For each acquisition, loop through all g files -> Extract sync pulse and error signal and preprocess LFP

        for nfolder = 1:length(DIR)
            cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))

            % Probe 1
            options= [];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']);

                EPHYS_DIR = dir(fullfile(options.EPHYS_DATAPATH,'..','*syncpulseTimes.mat'));
                EPHYS_DIR = [];
                if ~isempty(EPHYS_DIR)
                    load(fullfile(options.EPHYS_DATAPATH,'..',EPHYS_DIR.name))
                else
                    [ap_file,lf_file] = findImecBinFile(options.EPHYS_DATAPATH);
                    % if ~isempty(lf_file)
                    %     binpath = fullfile(options.EPHYS_DATAPATH,lf_file);
                    % else
                    fprintf('\n'); warning('Extracting sync-pulses from AP file...this takes about 5 minutes or more'); fprintf('\n');
                    binpath = fullfile(options.EPHYS_DATAPATH,ap_file);
                    % end
                    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
                    parseDate = date;
                    fname = extractBefore(ap_file,'.imec')
                    save(fullfile(options.EPHYS_DATAPATH,'..',[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate','-v7.3'); % used to be saved inside the ephys probe folder, now just same place where niqd is saved
                end

                NIDQ_DIR = dir(fullfile(options.EPHYS_DATAPATH,'..','*nidq*'));
                if ~isempty(NIDQ_DIR)
                    extract_and_align_nidq_signals(options)
                end
                preprocess_and_save_LFP(options)
            end

            % Probe 2
            options=[];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']);
                % extract_and_align_nidq_signals(options)
                preprocess_and_save_LFP(options)
            end
            
            if ~isempty(syncTimes_ephys.Error_on_idx)
                acquisition_id = str2double(DIR(nfolder).folder(end));

                Error_acquisition = [Error_acquisition; acquisition_id];
                Error_sessions = [Error_sessions; {DIR(nfolder).name}];
                Error_idx = [Error_idx; syncTimes_ephys.Error_on_idx];
                Error_timestamp = [Error_timestamp; syncTimes_ephys.Error_on_timestamp];
                
            end
    

        end
    end
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))
    if ~isempty(Error_acquisition)
        T = table(Error_acquisition,Error_sessions, Error_idx, Error_timestamp, 'VariableNames', {'acquisition','session','error_sample', 'error_timestamp'});
        writetable(T,'error_info.csv','Delimiter',',');
    else
        T = table('Size', [0, 4],'VariableNames', {'acquisition','session','error_sample', 'error_timestamp'},'VariableTypes', {'double','double', 'double', 'double'});
        writetable(T,'error_info.csv','Delimiter',',');
    end
end


%%% Rename to have animal id in ephys folder
SUBJECT = 'M24064';
all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys','20*'));

for ndate = 1:length(all_DIR)
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))

    session_DIR= dir('20*');

    Error_sessions = [];
    Error_idx = [];
    Error_timestamp=[];
    Error_acquisition = [];

    for nsession = 1:length(session_DIR)
        cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
        DIR = dir('*_g*'); % For each acquisition, loop through all g files -> Extract sync pulse and error signal and preprocess LFP

        for nfolder = 1:length(DIR)
            cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))

            % Probe 1
            options= [];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']);
                cd(options.EPHYS_DATAPATH)

                gfile_DIR = dir('*imec0*');

                for iFile = 1:length(gfile_DIR)
                    if ~contains(gfile_DIR(iFile).name,SUBJECT)

                        movefile(gfile_DIR(iFile).name,[SUBJECT,'_',gfile_DIR(iFile).name])
                    end

                end
            end

            % Probe 2
            options= [];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']);
                cd(options.EPHYS_DATAPATH)

                gfile_DIR = dir('*imec1*');

                for iFile = 1:length(gfile_DIR)
                    if ~contains(gfile_DIR(iFile).name,SUBJECT)

                        movefile(gfile_DIR(iFile).name,[SUBJECT,'_',gfile_DIR(iFile).name])
                    end

                end
            end

            % Rename all the folders and files other than ephys (Nidq and etc)
            cd(fullfile(DIR(nfolder).folder,DIR(nfolder).name))
            all_DIR_this_folder = dir('*g*');
            for iFile = 1:length(all_DIR_this_folder)
                if ~contains(all_DIR_this_folder(iFile).name,SUBJECT)
                    movefile(all_DIR_this_folder(iFile).name,[SUBJECT,'_',all_DIR_this_folder(iFile).name])
                end
            end

            % Rename aquisition folder
            cd(fullfile(DIR(nfolder).folder,DIR(nfolder).name,'..'))
            if ~contains(DIR(nfolder).name,SUBJECT)
                movefile(DIR(nfolder).name,[SUBJECT,'_',DIR(nfolder).name])
            end
        end
    end
end



SUBJECT = 'M24064';
all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys','20*'));

for ndate = 12:length(all_DIR)
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))

    session_DIR= dir('20*');

    Error_sessions = [];
    Error_idx = [];
    Error_timestamp=[];
    Error_acquisition = [];

    for nsession = 1:length(session_DIR)
        cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
        DIR = dir('*_g*'); % For each acquisition, loop through all g files -> Extract sync pulse and error signal and preprocess LFP

        for nfolder = 1:length(DIR)
            syncTimes_ephys=[];
            cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))

            % Probe 1
            options= [];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']);

                EPHYS_DIR = dir(fullfile(options.EPHYS_DATAPATH,'..','*syncpulseTimes.mat'));
                % EPHYS_DIR = [];
                if ~isempty(EPHYS_DIR)
                    load(fullfile(options.EPHYS_DATAPATH,'..',EPHYS_DIR.name))
                else
                    [ap_file,lf_file] = findImecBinFile(options.EPHYS_DATAPATH);
                    % if ~isempty(lf_file)
                    %     binpath = fullfile(options.EPHYS_DATAPATH,lf_file);
                    % else
                    fprintf('\n'); warning('Extracting sync-pulses from AP file...this takes about 5 minutes or more'); fprintf('\n');
                    binpath = fullfile(options.EPHYS_DATAPATH,ap_file);
                    % end
                    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
                    parseDate = date;
                    fname = extractBefore(ap_file,'.imec')
                    save(fullfile(options.EPHYS_DATAPATH,'..',[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate','-v7.3'); % used to be saved inside the ephys probe folder, now just same place where niqd is saved
                end

                NIDQ_DIR = dir(fullfile(options.EPHYS_DATAPATH,'..','*nidq*'));
                if ~isempty(NIDQ_DIR)
                    extract_and_align_nidq_signals(options)
                end
                preprocess_and_save_LFP(options)
            end

            % Probe 2
            options=[];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']);
                % extract_and_align_nidq_signals(options)
                preprocess_and_save_LFP(options)
            end
            
            if ~isempty(syncTimes_ephys.Error_on_idx)
                acquisition_id = str2double(DIR(nfolder).folder(end));

                Error_acquisition = [Error_acquisition; acquisition_id];
                Error_sessions = [Error_sessions; {DIR(nfolder).name}];
                Error_idx = [Error_idx; syncTimes_ephys.Error_on_idx];
                Error_timestamp = [Error_timestamp; syncTimes_ephys.Error_on_timestamp];
                
            end
    

        end
    end
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))
    if ~isempty(Error_acquisition)
        T = table(Error_acquisition,Error_sessions, Error_idx, Error_timestamp, 'VariableNames', {'acquisition','session','error_sample', 'error_timestamp'});
        writetable(T,'error_info.csv','Delimiter',',');
    else
        T = table('Size', [0, 4],'VariableNames', {'acquisition','session','error_sample', 'error_timestamp'},'VariableTypes', {'double','double', 'double', 'double'});
        writetable(T,'error_info.csv','Delimiter',',');
    end
end


% SUBJECT = 'M24072';
% all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'analysis','20*'));
% 
% for ndate = 1:length(all_DIR)
%     cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))
% 
%     allFilesAndFolders = dir;
% 
%     % Extract only the folders
%     session_DIR = allFilesAndFolders([allFilesAndFolders.isdir]);
% 
%     for nsession = 1:length(session_DIR)
%         if ~strcmp(session_DIR(nsession).name, '.') && ~strcmp(session_DIR(nsession).name, '..')
% 
%             cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
%             old_metafiles = dir('202*.meta');
% 
%             for n = 1:length(old_metafiles)
%                 movefile(old_metafiles(n).name, [SUBJECT,'_',old_metafiles(n).name]);
%                 % delete(old_metafiles(n).name)
%             end
%         end
% 
%     end
% end
% 

