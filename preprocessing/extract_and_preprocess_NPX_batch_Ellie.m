function extract_and_preprocess_NPX_batch_Ellie(experiment_info,Stimulus_type)
% Function to import and align bonsai data and spike data
% This is a batch processing version of extract_and_preprocess_NPX function

% Inputs:
% 1. experiment_info
% Based on loaded experiment_info, it may process sessions from just one animal or
% sessions from multiple animals.
% 2. Stimulus_type
% Specifiy which stimulus type to process -> go through different
% extraction function.

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(stimulus_name) % If stimulus not existing for this session
        continue
    end
    
    for n = 1:length(session_info) % just in case there might be multiple recording for the same stimulus type (e.g. PRE RUN POST Masa2tracks)

        
        if str2num(session_info(n).probe(1).SESSION) < 20240401
            % Old Part 1 Align two NPX probes if there are two probes
            if length(session_info(n).probe) >1
                align_probes_NX1(session_info(n).probe);
            end
        else % New Part 1 Extract Nidq sync pulse and event times and ephys sync pulse
            % Currently not using CatGT and Tprime
            options = session_info(n).probe(1);
            extract_and_align_nidq_signals(options);
            EPHYS_parent_folder = cd(fullfile(options.EPHYS_DATAPATH,'..','..'));
            preprocess_and_save_LFP(options)

            for nprobe = 1:length(session_info(n).probe)
                DIR = dir(fullfile(session_info(n).probe(nprobe).EPHYS_DATAPATH,'*tcat*.meta'));
                for nfile = 1:length(DIR)
                    copyfile(fullfile(DIR(nfile).folder,DIR(nfile).name),fullfile(session_info(n).probe(nprobe).ANALYSIS_DATAPATH,DIR(nfile).name))
                end
            end
            disp('Meta files of preprocessed LFP transfered')
        end

        % Part 2 Import and align bonsai data to Spikeglx time (always to probe 1)
        options = session_info(n).probe(1);
        if contains(Stimulus_type,'Masa2tracks')
            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,sprintf("extracted_behaviour%s.mat",erase(stimulus_name{n},Stimulus_type))));
        else
            DIR = dir(fullfile(options.ANALYSIS_DATAPATH,"extracted_behaviour*.mat"));
        end
        % DIR = [];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if isempty(DIR)
            disp('process and extract behavioural data')
            if contains(Stimulus_type,'Sleep')
               
                [Behaviour] = import_and_align_Bonsai_Sleep_Ellie(stimulus_name{n},session_info(n).probe);
          
            else % Else just standard visual stimuli such as Sparse Noise, checkerboard and static grating etc

                [Behaviour,Task_info,Peripherals]  = import_and_align_visual_stimuli_Bonsai_Ellie(stimulus_name{n},options);
            end

            if exist(options.ANALYSIS_DATAPATH) == 0
                mkdir(options.ANALYSIS_DATAPATH)
            end



            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')

            if exist('Task_info', 'var')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'),'Task_info')
            end

            if exist('Peripherals', 'var')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_peripherals.mat'),'Peripherals')
            end
            
        end
% 
        % Part 3 extract spike data
        options = session_info(n).probe(1);
        if contains(Stimulus_type,'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,...
                sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},Stimulus_type))))
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))
        end


        DIR_SORTER = dir(options.SORTER_DATAPATH);
        %DIR_KS = dir(options.KS_DATAPATH);

        if isempty(DIR_SORTER) % if spike interface sorter folder is not present, skip
              continue
        end
        
        disp(['Checking path: ', options.segment_frames])
        if isempty(dir(options.segment_frames))
            options.segment_frames = fullfile(options.EPHYS_DATAPATH,'..','..','..',['probe',num2str(options.probe_id),'segment_frames.csv']);
        end
        
        segment_frames_table = readtable(options.segment_frames);
        segment_frames = table2array(segment_frames_table(:,1));

        session_id=extractAfter(options.EPHYS_DATAPATH,['\',options.SESSION,'\',options.SESSION,'_']);
        if isempty(session_id)
            session_id=extractAfter(options.EPHYS_DATAPATH,['\',options.SESSION,'\',options.SUBJECT,'_',options.SESSION,'_']);
        end
        session_id=str2num(session_id(1));
        
        if sum(contains(segment_frames,[num2str(session_id),'_g',num2str(options.gFileNum)]))==0
            disp('Session without spike sorting is skipped')
            continue
        end

        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            DIR_SORTER = dir(options.SORTER_DATAPATH);
            %DIR_KS = dir(options.KS_DATAPATH);

            if ~isempty(DIR_SORTER) % if spike interface sorter folder is present

                temp = dir(fullfile(options.SORTER_DATAPATH,'waveform','kilosort2'));
                if ~isempty(temp)
                    [clusters_ks2(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS2','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
                end
                temp = dir(fullfile(options.SORTER_DATAPATH,'waveform','kilosort3'));
                if ~isempty(temp)
                    [clusters_ks3(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS3_original','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
                end

                temp = dir(fullfile(options.SORTER_DATAPATH,'waveform','kilosort4'));
                if ~isempty(temp)
                    [clusters_ks4(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS4_original','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
                end

            elseif ~isempty(DIR_KS)% elseif original KS3 folder is present
                [clusters(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','off','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
            end
            %     [all_clusters chan_config sorted_config] = extract_clusters_NPX(options,'group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
        end

        %         spikes = clusters;
        %         fields_to_remove = {'spike_count_raw','spike_count_smoothed','zscore_smoothed'};
        %         clusters = rmfield(clusters,fields_to_remove);
        %
        if contains(Stimulus_type,'Masa2tracks') 
            % If Masa2tracks, PRE, RUN and/or POST saved in one folder
            if ~isempty(DIR_SORTER) % if spike interface sorter folder is present


                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort2'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,...
                        sprintf('extracted_clusters_ks2%s.mat',erase(stimulus_name{n},Stimulus_type))),'clusters_ks2')
                end

                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort3'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,...
                        sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},Stimulus_type))),'clusters_ks3')
                end

                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort4'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,...
                        sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},Stimulus_type))),'clusters_ks4')
                end
            elseif ~isempty(DIR_KS)% elseif original KS3 folder is present
                save(fullfile(options.ANALYSIS_DATAPATH,...
                    sprintf('extracted_clusters%s.mat',erase(stimulus_name{n},Stimulus_type))),'clusters')
            end

        else
            if ~isempty(DIR_SORTER) % if spike interface sorter folder is present

                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort2'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks2.mat'),'clusters_ks2')
                end

                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort3'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'),'clusters_ks3')
                end

                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort4'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'),'clusters_ks4')
                end

             elseif ~isempty(DIR_KS)% elseif original KS3 folder is present
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters.mat'),'clusters')
            end

        end
        
        clear clusters_ks2 clusters_ks3

    end
    close all
end
