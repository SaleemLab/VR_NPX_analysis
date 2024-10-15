function extract_and_preprocess_NPX_batch(experiment_info,Stimulus_type)
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
%         DIR = [];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if isempty(DIR)
            disp('process and extract behavioural data')
            if contains(Stimulus_type,'OpenField')|contains(Stimulus_type,'Sleep')
                if contains(Stimulus_type,'Sleep')
                    [Behaviour] = import_and_align_Bonsai_Sleep(stimulus_name{n},session_info(n).probe);
                else
                    [Behaviour] = import_and_align_Bonsai_OpenField(stimulus_name{n},session_info(n).probe);
                end
            elseif contains(Stimulus_type,'Masa2tracks') |contains(Stimulus_type,'Track') 
%                 [Behaviour,Task_info,Peripherals] = import_and_align_Masa_VR_Bonsai(stimulus_name{n},options);

                if str2num(options(1).SESSION) < 20240401 % if session before 2024/04/01
                    [Behaviour,Task_info,Peripherals]  = import_and_align_Masa_VR_Bonsai(stimulus_name{n},options);
                else
                    [Behaviour,Task_info,Peripherals]  = import_and_align_Masa_VR_Bonsai_MatrixRig(stimulus_name{n},options);
                end
                
            elseif contains(Stimulus_type,'Diao')

            elseif contains(Stimulus_type,'Edd')

            else % Else just standard visual stimuli such as Sparse Noise, checkerboard and static grating etc

                if str2num(options(1).SESSION) < 20240401 % if session before 2024/04/01
                    [Behaviour,Task_info,Peripherals]  = import_and_align_visual_stimuli_Bonsai(stimulus_name{n},options);
                else
                    [Behaviour,Task_info,Peripherals]  = import_and_align_visual_stimuli_Bonsai_MatrixRig(stimulus_name{n},options);
                end
            end

            if exist(options.ANALYSIS_DATAPATH) == 0
                mkdir(options.ANALYSIS_DATAPATH)
            end

            if contains(Stimulus_type,'Masa2tracks')
                % If Masa2tracks, PRE, RUN and/or POST saved in one folder
                save(fullfile(options.ANALYSIS_DATAPATH,...
                    sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},Stimulus_type))),'Behaviour')
                save(fullfile(options.ANALYSIS_DATAPATH,...
                    sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},Stimulus_type))),'Task_info')
                save(fullfile(options.ANALYSIS_DATAPATH,...
                    sprintf('extracted_peripherals%s.mat',erase(stimulus_name{n},Stimulus_type))),'Peripherals')
            elseif contains(Stimulus_type,'OpenField')|contains(Stimulus_type,'Sleep')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')
%                 save(fullfile(options.ANALYSIS_DATAPATH,'extracted_peripherals.mat'),'Peripherals')
            else
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')
                save(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'),'Task_info')
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
        DIR_KS = dir(options.KS_DATAPATH);

        if isempty(DIR_SORTER) % if spike interface sorter folder is not present, skip
              continue
        end

        segment_frames = readtable(options.segment_frames);
        session_id=extractAfter(options.EPHYS_DATAPATH,['\',options.SESSION,'\',options.SESSION,'_']);
        session_id=str2num(session_id(1));

        if sum(table2array(segment_frames(:,3))==session_id)==0
            disp('Session without sipike sorting is skipped due to imFlagError')
            continue
        end

        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            DIR_SORTER = dir(options.SORTER_DATAPATH);
            DIR_KS = dir(options.KS_DATAPATH);
            
            
            if ~isempty(DIR_SORTER) % if spike interface sorter folder is present
                temp = dir(fullfile(options.SORTER_DATAPATH,'waveform','kilosort3','extensions','quality_metrics'));
                if ~isempty(temp)
                    [clusters_ks3_original(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS3_original','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
                end
%                 
                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort2'));
                if ~isempty(temp)
                    [clusters_ks2(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS2','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
                end
                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort3'));
                if ~isempty(temp)
                    [clusters_ks3(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS3','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
                end

                temp = dir(fullfile(options.SORTER_DATAPATH,'sorters','kilosort4'));
                if ~isempty(temp)
                    [clusters_ks4(nprobe) chan_config sorted_config] = extract_clusters_NPX(options,'sorter','KS4','group','all clusters','tvec',Behaviour.tvec,'SR',mean(1./diff(Behaviour.tvec)));
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
%             save(fullfile(options.ANALYSIS_DATAPATH,...
%                 sprintf('extracted_spikes%s.mat',erase(stimulus_name{n},Stimulus_type))),'spikes')
        else
            if ~isempty(DIR_SORTER) % if spike interface sorter folder is present
                temp = dir(fullfile(options.SORTER_DATAPATH,'waveform','kilosort3','extensions','quality_metrics'));
                if ~isempty(temp)
                    save(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3_original.mat'),'clusters_ks3_original')
                end

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
%             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spikes.mat'),'spikes')
        end
        
        clear clusters_ks2 clusters_ks3

    end
    close all
end
