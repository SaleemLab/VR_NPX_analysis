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
        dataCell = readcell('Z:\ibn-vision\USERS\Masa\recording_info\session_stimuli_table','Sheet',SUBJECTS{nsubject}); % this can create missing object when session id is missing.
        % Convert the cell array to a table
        stimuli_info_table = cell2table(dataCell(2:end,:),"VariableNames",stimuli_info.Properties.VariableNames); % not the best way but works
        
    elseif contains(options,'V1-MEC')
        stimuli_info = readtable(['Z:\ibn-vision\USERS\Diao\ephys_recording_info\session_stimuli_table_DT'],'Sheet',SUBJECTS{nsubject});
        dataCell = readcell('Z:\ibn-vision\USERS\Diao\ephys_recording_info\session_stimuli_table_DT','Sheet',SUBJECTS{nsubject});
        % Convert the cell array to a table
        stimuli_info_table = cell2table(dataCell(2:end,:),"VariableNames",stimuli_info.Properties.VariableNames);
    end

    if contains('session_id',stimuli_info.Properties.VariableNames)
        stimuli_info.session_id = stimuli_info_table.session_id;
    end
    clear stimuli_info_table

    all_dates_this_animal = unique(stimuli_info.date);
   
    if contains('session_id',stimuli_info.Properties.VariableNames) % if session number is specifeid (since chronic recording)
        all_sessions_this_animal = stimuli_info.session_id; % this only works if there is at least something in session id column

        if ~iscell(all_sessions_this_animal) % if not in cell structure, convert to cell
            all_sessions_this_animal = num2cell(all_sessions_this_animal);
        end

        for nsession = 1:length(all_sessions_this_animal)
            if ismissing(all_sessions_this_animal{nsession}) % if missing, becomes nan
                all_sessions_this_animal{nsession} = nan;
            elseif ~isstring(all_sessions_this_animal{nsession}) % if not string, becomes string
                all_sessions_this_animal{nsession} = num2str(all_sessions_this_animal{nsession});
            end
        end

    else
        all_sessions_this_animal = cell(length(stimuli_info.date),1);
        all_sessions_this_animal(:) = {nan};
    end


    if ~iscell(stimuli_info.probe_number)
        stimuli_info.probe_number = cellstr(num2str(stimuli_info.probe_number));
        if contains(options,'bilateral')
            stimuli_info.probe_hemisphere = cellstr(num2str(stimuli_info.probe_hemisphere));
        elseif contains(options,'V1-MEC')

        end
    end

    if contains(options,'V1-MEC')
        if ~iscell(stimuli_info.V1)
            stimuli_info.V1 = cellstr(num2str(stimuli_info.V1));
        end
        if ~iscell(stimuli_info.MEC)
            stimuli_info.MEC = cellstr(num2str(stimuli_info.MEC));
        end
    end
        
    for nsession = 1:length(all_dates_this_animal)
        
        experiment_info(nexperiment).experiment_ID = [SUBJECTS{nsubject},'_',num2str(all_dates_this_animal(nsession))];
        experiment_info(nexperiment).date = all_dates_this_animal(nsession);
        experiment_info(nexperiment).subject = SUBJECTS{nsubject};

        DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');

%         EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECTS{nsubject},'ephys',num2str(all_sessions_this_animal(nsession)));
%         EPHYS_DATAPATH_OpenField = fullfile(DATAPATH,SUBJECTS{nsubject},'ephys',[num2str(all_sessions_this_animal(nsession)),'_OpenField']);

        experiment_info(nexperiment).ANALYSIS_DATAPATH = fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_dates_this_animal(nsession)));
        experiment_info(nexperiment).BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECTS{nsubject},'stimuli',num2str(all_dates_this_animal(nsession)));

        experiment_info(nexperiment).gFileNum = stimuli_info.gs_number(stimuli_info.date == all_dates_this_animal(nsession));
        StimulusName = stimuli_info.stimulus_type(stimuli_info.date == all_dates_this_animal(nsession));

        experiment_info(nexperiment).StimulusName = StimulusName;
        experiment_info(nexperiment).probe_number = stimuli_info.probe_number(stimuli_info.date == all_dates_this_animal(nsession));

        [stimulus_counter,nameCountMap] = assignNumbersBasedOnCount(StimulusName);
        all_session_id_this_date = all_sessions_this_animal(stimuli_info.date == all_dates_this_animal(nsession)); %The session id of each stimulus (in caase many recordings of the same stimuli)

        if contains(options,'bilateral')
            experiment_info(nexperiment).probe_hemisphere = stimuli_info.probe_hemisphere(stimuli_info.date == all_dates_this_animal(nsession));
        elseif contains(options,'V1-MEC')
            experiment_info(nexperiment).probe_V1 = stimuli_info.V1(stimuli_info.date == all_dates_this_animal(nsession));
            experiment_info(nexperiment).probe_MEC = stimuli_info.MEC(stimuli_info.date == all_dates_this_animal(nsession));
        end
        
        for nstimuli = 1:length(StimulusName)

            % Bonsai data path
            if contains(StimulusName{nstimuli},'StaticGrating')
                this_stimulus = 'StaticGrating';
                bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus, '*.csv']));

            else
                this_stimulus = StimulusName{nstimuli};
                if isequal(this_stimulus , 'Track')
                    bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus, '_*.csv']));
                elseif contains(this_stimulus,'Masa2tracks')
                    bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*','Masa2tracks', '*.csv']));
                else
                    bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus, '*.csv']));
                end
            end

            [~,idx] = sort([bonsai_files.datenum]);

            file_time = extractPattern(bonsai_files(idx(stimulus_counter(nstimuli))).name); % sometimes the exact seconds can be off for different bonsai files, so search is currently based on hour and minute (may need optimisation);
            bonsai_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',file_time,'*.csv']));


            if ~isempty(dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus,'*.bin'])))

                % Unfortunately... sometimes bin file time stamp can be 1 or 2
                % mins later than others....
                bin_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',this_stimulus,'*.bin']));
                [~,idx] = sort([bin_files.datenum]);

                file_time1 = extractPattern(bin_files(idx(stimulus_counter(nstimuli))).name); % sometimes the exact seconds can be off for different bonsai files, so search is currently based on hour and minute (may need optimisation);
                bin_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',file_time1,'*.bin']));

                bonsai_files(size(bonsai_files,1)+1) = bin_files;
            end


            % Add face data (saved as mat file) if exist
            if ~isempty(dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',file_time,'*proc.mat'])))
                mat_files = dir(fullfile(experiment_info(nexperiment).BONSAI_DATAPATH,['*',file_time,'*proc.mat']));
                bonsai_files(size(bonsai_files,1)+1) = mat_files;
            end

            probe_id = str2num(experiment_info(nexperiment).probe_number{nstimuli});

            if contains(options,'bilateral')
                probe_hemisphere = str2num(experiment_info(nexperiment).probe_hemisphere{nstimuli});
            elseif contains(options,'V1-MEC')
                probe_V1 = str2num(experiment_info(nexperiment).probe_V1{nstimuli});
                probe_MEC = str2num(experiment_info(nexperiment).probe_MEC{nstimuli});
            end

            % If it is normal
            if isnan(all_session_id_this_date{nstimuli})
                EPHYS_DATAPATH_temp =fullfile(DATAPATH,SUBJECTS{nsubject},'ephys',num2str(all_dates_this_animal(nsession)));
            else
                EPHYS_DATAPATH_temp =fullfile(DATAPATH,SUBJECTS{nsubject},'ephys',[num2str(all_dates_this_animal(nsession)),'_',all_session_id_this_date{nstimuli}]);
            end

            for nprobe = 1:length(probe_id)
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).SUBJECT = SUBJECTS{nsubject};
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).SESSION = num2str(all_dates_this_animal(nsession));
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).StimulusName = StimulusName{nstimuli};

                if contains(StimulusName{nstimuli},'Masa2tracks')
                    experiment_info(nexperiment).session(nstimuli).probe(nprobe).ANALYSIS_DATAPATH =...
                        fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_dates_this_animal(nsession)),'Masa2tracks');
                else
                    if  nameCountMap(StimulusName{nstimuli}) == 1 % if this stimulus only happend once
                        experiment_info(nexperiment).session(nstimuli).probe(nprobe).ANALYSIS_DATAPATH =...
                            fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_dates_this_animal(nsession)),StimulusName{nstimuli});
                    else % if multiple recording on the same day
                        experiment_info(nexperiment).session(nstimuli).probe(nprobe).ANALYSIS_DATAPATH =...
                            fullfile(DATAPATH,SUBJECTS{nsubject},'analysis',num2str(all_dates_this_animal(nsession)),[StimulusName{nstimuli},'_',num2str(stimulus_counter(nstimuli))]);
                    end
                end

                folderName = findGFolder(EPHYS_DATAPATH_temp,experiment_info(nexperiment).gFileNum(nstimuli));
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,folderName,[folderName,'_imec',num2str(probe_id(nprobe))]);

                SORTER_DATAPATH =  fullfile(EPHYS_DATAPATH_temp,['probe',num2str(probe_id(nprobe))]); % NEW spikeinterface sorter datapath
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).SORTER_DATAPATH = SORTER_DATAPATH;

                if ~exist(fullfile(EPHYS_DATAPATH_temp,'kilosort'), 'dir') % Old KS datapath
                    KS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,['kilosort_probe_',num2str(probe_id(nprobe)+1)]); % e.g. imec0 is in kilosort_probe_1 folder

                else
                    KS_DATAPATH = fullfile(EPHYS_DATAPATH_temp,'kilosort');
                end

                experiment_info(nexperiment).session(nstimuli).probe(nprobe).KS_DATAPATH = KS_DATAPATH;
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).gFileNum = experiment_info(nexperiment).gFileNum(nstimuli);
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).probe_id = probe_id(nprobe);

                if contains(options,'bilateral')
                    experiment_info(nexperiment).session(nstimuli).probe(nprobe).probe_hemisphere = probe_hemisphere(nprobe);
                elseif contains(options,'V1-MEC')
                    experiment_info(nexperiment).session(nstimuli).probe(1).probe_V1 = probe_V1;
                    experiment_info(nexperiment).session(nstimuli).probe(1).probe_MEC = probe_MEC;
                end

                experiment_info(nexperiment).session(nstimuli).probe(nprobe).MAP_FILE = ...
                    fullfile(KS_DATAPATH,[SUBJECTS{nsubject},'_',num2str(all_dates_this_animal(nsession)),'_g0','_tcat.imec',num2str(probe_id(nprobe)),'.ap_kilosortChanMap.mat']);

                if contains(StimulusName{nstimuli},'OpenField')
                    experiment_info(nexperiment).session(nstimuli).probe(nprobe).KS_CATGT_FNAME =...
                        fullfile(['CatGT_',SUBJECTS{nsubject},'_',num2str(all_dates_this_animal(nsession)),'_OpenField.log']);
                else
                    experiment_info(nexperiment).session(nstimuli).probe(nprobe).KS_CATGT_FNAME =...
                        fullfile(['CatGT_',SUBJECTS{nsubject},'_',num2str(all_dates_this_animal(nsession)),'.log']);
                end

                experiment_info(nexperiment).session(nstimuli).probe(nprobe).segment_frames =...
                    fullfile(SORTER_DATAPATH,'sorters','segment_frames.csv');

                experiment_info(nexperiment).session(nstimuli).probe(nprobe).BONSAI_DATAPATH = experiment_info(nexperiment).BONSAI_DATAPATH;
                experiment_info(nexperiment).session(nstimuli).probe(nprobe).bonsai_files_names = {bonsai_files.name};

            end
        end

        nexperiment = nexperiment + 1;
    end
end


end


function [output,nameCountMap] = assignNumbersBasedOnCount(names)
    % Create a map to store counts
    nameCountMap = containers.Map;

    % Initialize output array
    output = zeros(size(names));

    % Iterate through each stimulus
    for i = 1:length(names)
        % Get the current name

        currentName = names{i};
        if contains(currentName,'Masa2tracks')% if it is Masa2tracks
            currentName = 'Masa2tracks';
        end

        % Check if the stimulus is already in the map
        if isKey(nameCountMap, currentName)
            % Increment the count if the stimulus is present
            count = nameCountMap(currentName);
            nameCountMap(currentName) = count + 1;
        else
            % Add the stimulus to the map if it's not present
            nameCountMap(currentName) = 1;
        end

        % Assign the count to the output array
        output(i) = nameCountMap(currentName);
    end
end


function extractedPattern = extractPattern(filename)
    % Define the regular expression pattern
    pattern = '(\d{4}-\d{2}-\d{2}T\d{2}_\d{2})';

    % Find the match in the filename using regular expression
    match = regexp(filename, pattern, 'match');

    % Check if a match is found
    if ~isempty(match)
        extractedPattern = match{1};
    else
        extractedPattern = 'No match found';
    end
end