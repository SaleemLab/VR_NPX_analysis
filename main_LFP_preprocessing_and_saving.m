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


all_DIR= dir(fullfile('Z:\ibn-vision\DATA\SUBJECTS',SUBJECT,'ephys','20*'));

for ndate = 1:length(all_DIR)
    cd(fullfile(all_DIR(ndate).folder,all_DIR(ndate).name))

    session_DIR= dir('20*');
    for nsession = 1:length(session_DIR)
        cd(fullfile(session_DIR(nsession).folder,session_DIR(nsession).name))
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
    end
end