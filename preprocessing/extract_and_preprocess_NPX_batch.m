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
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(stimulus_name) % If stimulus not existing for this session
        continue
    end
    
    for n = 1:length(session_info) % just in case there might be multiple recording for the same stimulus type (e.g. PRE RUN POST Masa2tracks)

        % Part 1 Align two NPX probes if there are two probes
        if length(session_info(n).probe) >1
            align_probes_NX1(session_info(n).probe);
        end

        % Part 2 Import and align bonsai data to Spikeglx time (always to probe 1)
        options = session_info(n).probe(1);
        if contains(Stimulus_type,'OpenField')
            [Behaviour,~] = import_and_align_Bonsai_OpenField(stimulus_name{n},session_info(n).probe);
        elseif contains(Stimulus_type,'Masa2tracks')
           [Behaviour,Task_info,Peripherals] = import_and_align_Masa_VR_Bonsai(stimulus_name{n},options);

        elseif contains(Stimulus_type,'Diao')
            
        elseif contains(Stimulus_type,'Edd')
            
        else % Else just standard visual stimuli such as Sparse Noise, checkerboard and static grating etc
          [Behaviour,Task_info,Peripherals]  = import_and_align_visual_stimuli_Bonsai(stimulus_name{n},options);
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
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'),'Task_info')
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_peripherals.mat'),'Peripherals')
        end

        % Part 3 extract spike data
        for nprobe = 1:length(session_info(n).probe) >1
            options = session_info(n).probe(nprobe);
            if contains(Stimulus_type,'Masa2tracks')
                load(fullfile(options.ANALYSIS_DATAPATH,...
                    sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},Stimulus_type))))
            else
                load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'))
            end
            [clusters chan_config sorted_config] = extract_clusters_NPX(options,varargin)
        end

%         save(sprintf('extracted_position%s.mat',erase(stimulus_name{n},Stimulus_type)),'position')
        
%         if contains(stimulus_name{n},'RUN')
%             lap_times = extract_laps_masa(1,Behaviour,position)
%             save extracted_laps lap_times
%         end


        % figure
        % hold on
        % plot(MousePos.sglxTime,MousePos.pos)
        % scatter(MousePos.stimuli_onset(MousePos.stimuli_track == 1),1000*MousePos.stimuli_track(MousePos.stimuli_track == 1),'r')
        % scatter(MousePos.stimuli_onset(MousePos.stimuli_track == 2),100*MousePos.stimuli_track(MousePos.stimuli_track == 2),'b')
    end

end
