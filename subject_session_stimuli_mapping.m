function experiment_info = subject_session_stimuli_mapping(SUBJECTS,options)
% Function used to create basic information about each recording sessions
% recorded from the animals specified

% 05/01/2024 M.T updating experiment info to include open field data 
% (Currently open field ephys data is saved as a separate folder 
% (e.g. 20231215_OpenField) to avoid conflict with 20231215)
% We need to find a good way to save multiple ephys
% folders for the same day. Since spikeglx can not just concat all data
% from multiple sessions, currently, we need to process each folder separately... 
% Maybe with spikeinterface it can be sovled?


ROOTPATH = 'Z:\ibn-vision';

% if contains(options,'bilateral')
%     EXPERIMENT_INFO_PATH = 'Z:\ibn-vision\USERS\Masa\recording_info';
% else
%     
% end



experiment_info = [];
% SUBJECTS = {'M23017','M23028','M23029'};
nexperiment = 1;

for nsubject = 1:length(SUBJECTS)
    if contains(options,'bilateral')
        stimuli_info = readtable('Z:\ibn-vision\USERS\Masa\recording_info\session_stimuli_table','Sheet',SUBJECTS{nsubject});
    elseif contains(options,'V1-MEC')
           stimuli_info = readtable(['Z:\ibn-vision\USERS\Diao\ephys_recording_info\session_stimuli_table_DT'],'Sheet',SUBJECTS{nsubject});
    end

    all_sessions_this_animal = unique(stimuli_info.date);
    if ~iscell(stimuli_info.probe_number)
        stimuli_info.probe_number = cellstr(num2str(stimuli_info.probe_number));
        if contains(options,'bilateral')
            stimuli_info.probe_hemisphere = cellstr(num2str(stimuli_info.probe_hemisphere));
        elseif contains(options,'V1-MEC')
            stimuli_info.V1 = cellstr(num2str(stimuli_info.V1));
            stimuli_info.MEC = cellstr(num2str(stimuli_info.MEC));
        end
    end

    for nsession = 1:length(all_sessions_this_animal)

        experiment_info(nexperiment).experiment_ID = [SUBJECTS{nsubject},'_',num2str(all_sessions_this_animal(nsession))];
        experiment_info(nexperiment).subject = SUBJECTS{nsubject};
        experiment_info(nexperiment).session = all_sessions_this_animal(nsession);


        DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
        EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECTS{nsubject},'ephys',num2str(all_sessions_this_animal(nsession)));
        EPHYS_DATAPATH_OpenField = fullfile(DATAPATH,SUBJECTS{nsubject},'ephys',[num2str(all_sessions_this_animal(nsession)),'_OpenField']);

        experiment_info(nexperiment).ANALYSIS_DATAPATH = fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_sessions_this_animal(nsession)));

        experiment_info(nexperiment).BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECTS{nsubject},'stimuli',num2str(all_sessions_this_animal(nsession)));

        experiment_info(nexperiment).gFileNum = stimuli_info.gs_number(stimuli_info.date == all_sessions_this_animal(nsession));
        StimulusName = stimuli_info.stimulus_type(stimuli_info.date == all_sessions_this_animal(nsession));
        experiment_info(nexperiment).StimulusName = StimulusName;
        experiment_info(nexperiment).probe_number = stimuli_info.probe_number(stimuli_info.date == all_sessions_this_animal(nsession));
        if contains(options,'bilateral')
            experiment_info(nexperiment).probe_hemisphere = stimuli_info.probe_hemisphere(stimuli_info.date == all_sessions_this_animal(nsession));
        end
        
        file_counter = 1;

        for nstimuli = 1:length(StimulusName)
            % Bonsai data path
            if contains(StimulusName{nstimuli},'PRE') | contains(StimulusName{nstimuli},'RUN') | contains(StimulusName{nstimuli},'POST')
                bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['Masa2tracks_WheelLog', '*']));
                [~,idx] = sort([bonsai_files.datenum]);
                
                file_time = erase(bonsai_files(idx(file_counter)).name,'Masa2tracks_WheelLog');
                file_time = file_time(1:end-6); % sometimes the exact seconds can be off for different bonsai files, search based on hour and minute (may need optimisation);
                bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',file_time,'*.csv']));

                file_counter = file_counter + 1;
            elseif contains(StimulusName{nstimuli},'StaticGratings')
                this_stimulus = 'StaticGratings';
                bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,[this_stimulus, '*.csv']));

            else 
                this_stimulus = StimulusName{nstimuli};
                if isequal(this_stimulus , 'Track')
                    bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus, '_*.csv']));
                else
                    bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus, '*.csv']));
                end
                
                
                if ~isempty(dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,[this_stimulus, '*.bin'])))
                    bin_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,[this_stimulus, '*.bin']));
                    bonsai_files(size(bonsai_files,1)+1) = bin_files;
                end
            end

            probe_id = str2num(experiment_info(nexperiment).probe_number{nstimuli});
            if contains(options,'bilateral')
                probe_hemisphere = str2num(experiment_info(nexperiment).probe_hemisphere{nstimuli});
            elseif contains(options,'V1-MEC')
                probe_V1 = str2num(experiment_info(nexperiment).probe_V1{nstimuli});
                probe_MEC = str2num(experiment_info(nexperiment).probe_MEC{nstimuli});
            end

            % If it is normal 
            if contains(StimulusName{nstimuli},'OpenField')
                EPHYS_DATAPATH_temp =EPHYS_DATAPATH_OpenField;
            else
                EPHYS_DATAPATH_temp =EPHYS_DATAPATH;
            end

            if length(probe_id) == 1

                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).SUBJECT = SUBJECTS{nsubject};
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).SESSION = num2str(all_sessions_this_animal(nsession));
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).StimulusName = StimulusName{nstimuli};
                folderName = findGFolder(EPHYS_DATAPATH_temp,experiment_info(nexperiment).gFileNum(nstimuli));
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,folderName,[folderName,'_imec',num2str(probe_id)]);
                
                if contains(StimulusName{nstimuli},'Masa2tracks')
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).ANALYSIS_DATAPATH =...
                        fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_sessions_this_animal(nsession)),'Masa2tracks');
                else
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).ANALYSIS_DATAPATH =...
                        fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_sessions_this_animal(nsession)),StimulusName{nstimuli});
                end

                if ~exist(fullfile(EPHYS_DATAPATH_temp,'kilosort'), 'dir')
                    KS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,['kilosort_probe_',num2str(probe_id+1)]); % e.g. imec0 is in kilosort_probe_1 folder

                else

                    KS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,'kilosort');
                end

                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).KS_DATAPATH = KS_DATAPATH;
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).gFileNum = experiment_info(nexperiment).gFileNum(nstimuli);
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).probe_id = probe_id;

                if contains(options,'bilateral')
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).probe_hemisphere = probe_hemisphere;
                elseif contains(options,'V1-MEC')
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).probe_V1 = probe_V1;
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).probe_MEC = probe_MEC;
                end

                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).MAP_FILE = ...
                    fullfile(KS_DATAPATH,[SUBJECTS{nsubject},'_',num2str(all_sessions_this_animal(nsession)),'_g0','_tcat.imec',num2str(probe_id),'.ap_kilosortChanMap.mat']);
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).KS_CATGT_FNAME =...
                    fullfile(['CatGT_',SUBJECTS{nsubject},'_',num2str(all_sessions_this_animal(nsession)),'.log']);

                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).BONSAI_DATAPATH = experiment_info(nexperiment).BONSAI_DATAPATH;
                experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).bonsai_files_names = {bonsai_files.name};
            else
                for nprobe = 1:length(probe_id)
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).SUBJECT = SUBJECTS{nsubject};
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).SESSION = num2str(all_sessions_this_animal(nsession));
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).StimulusName = StimulusName{nstimuli};

                    if contains(StimulusName{nstimuli},'Masa2tracks')
                        experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).ANALYSIS_DATAPATH =...
                            fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_sessions_this_animal(nsession)),'Masa2tracks');
                    else
                        experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).ANALYSIS_DATAPATH =...
                            fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_sessions_this_animal(nsession)),StimulusName{nstimuli});
                    end

                    folderName = findGFolder(EPHYS_DATAPATH_temp,experiment_info(nexperiment).gFileNum(nstimuli));
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,folderName,[folderName,'_imec',num2str(probe_id(nprobe))]);
                    KS_DATAPATH =  fullfile(EPHYS_DATAPATH_temp,['kilosort_probe_',num2str(nprobe)]);
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).KS_DATAPATH = KS_DATAPATH;
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).gFileNum = experiment_info(nexperiment).gFileNum(nstimuli);
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).probe_id = probe_id(nprobe);

                    if contains(options,'bilateral')
                        experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).probe_hemisphere = probe_hemisphere(nprobe);
                    elseif contains(options,'V1-MEC')
                        experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).probe_V1 = probe_V1;
                        experiment_info(nexperiment).stimuli_type(nstimuli).probe(1).probe_MEC = probe_MEC;
                    end

                    fullfile(experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).KS_DATAPATH...
                        ,[SUBJECTS{nsubject},'_',num2str(all_sessions_this_animal(nsession)),'_g',num2str(probe_id(nprobe)),'_tcat.imec0.ap_kilosortChanMap.mat']);

                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).MAP_FILE = ...
                        fullfile(KS_DATAPATH,[SUBJECTS{nsubject},'_',num2str(all_sessions_this_animal(nsession)),'_g0','_tcat.imec',num2str(probe_id(nprobe)),'.ap_kilosortChanMap.mat']);
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).KS_CATGT_FNAME =...
                        fullfile(['CatGT_',SUBJECTS{nsubject},'_',num2str(all_sessions_this_animal(nsession)),'.log']);

                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).BONSAI_DATAPATH = experiment_info(nexperiment).BONSAI_DATAPATH;
                    experiment_info(nexperiment).stimuli_type(nstimuli).probe(nprobe).bonsai_files_names = {bonsai_files.name};
                end
            end
        end

        nexperiment = nexperiment + 1;
    end
end


end