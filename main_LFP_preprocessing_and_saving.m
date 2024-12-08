%% LFP preprocess -> Catgt
% SUBJECT = 'M23087';
SUBJECT = 'M24016';
% SESSION = '20231212';
% date = '20240706';
SESSION = '20240706/20240706_0';


EPHYS_DATAPATH= fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys',SESSION);
cd(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys',SESSION))
DIR = dir('*_g*');

for nfolder = 1:length(DIR)
    options= [];
    temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']));
    if ~isempty(temp_DIR)
        options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']);
        preprocess_and_save_LFP(options)
    end

    temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']));
    if ~isempty(temp_DIR)
        options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']);
        preprocess_and_save_LFP(options)
    end
end

SUBJECT = 'M24016';
all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys','20*'));

for ndate = 1:length(all_DIR)
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))

    session_DIR= dir('20*');
    for nsession = 1:length(session_DIR)
        cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
        DIR = dir('*_g*');
        for nfolder = 1:length(DIR)
            cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
            options= [];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']);
                preprocess_and_save_LFP(options)
            end

%             temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']));
%             if ~isempty(temp_DIR)
%                 options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']);
%                 preprocess_and_save_LFP(options)
%             end


        end
    end
end


SUBJECT = 'M24062';
all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys','20*'));

for ndate = 8:length(all_DIR)
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))

    session_DIR= dir('20*');
    for nsession = 1:length(session_DIR)
        cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
        DIR = dir('*_g*');
        for nfolder = 1:length(DIR)
            cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
            options= [];
            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec0']);
                preprocess_and_save_LFP(options)
                DIR_ephys_sync = dir(fullfile(options.EPHYS_DATAPATH,'..','*syncpulseTimes.mat'))
                %     DIR = [];
                if ~isempty(DIR_ephys_sync)
                    load(fullfile(options.EPHYS_DATAPATH,'..',DIR_ephys_sync.name))
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
                    save(fullfile(options.EPHYS_DATAPATH,'..',[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate','-v7.3'); % used to be saved inside the ephys probe folder, now just same place where niqd is saved
                end
                extract_and_align_nidq_signals(options)
            end

            temp_DIR = dir(fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']));
            if ~isempty(temp_DIR)
                options.EPHYS_DATAPATH = fullfile(DIR(nfolder).folder,DIR(nfolder).name,[DIR(nfolder).name,'_imec1']);
                preprocess_and_save_LFP(options)
                DIR_ephys_sync = dir(fullfile(options.EPHYS_DATAPATH,'..','*syncpulseTimes.mat'))
                %     DIR = [];

                if ~isempty(DIR_ephys_sync)
                    load(fullfile(options.EPHYS_DATAPATH,'..',DIR_ephys_sync.name))
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
                    save(fullfile(options.EPHYS_DATAPATH,'..',[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate','-v7.3'); % used to be saved inside the ephys probe folder, now just same place where niqd is saved
                end

                extract_and_align_nidq_signals(options)
            end


        end
    end
end