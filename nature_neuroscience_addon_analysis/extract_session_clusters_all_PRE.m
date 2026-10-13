function extract_session_clusters_all_PRE()
    addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
    addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

    clear all
    SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
    option = 'bilateral';
    experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
    experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
    
    Stimulus_type = 'SleepChronic_PRE';

    session_clusters_all = struct();
    session_count = 0;

    for nsession =1:length(experiment_info)
        session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        
        if isempty(stimulus_name)
            continue
        end

        if length(stimulus_name) > 1
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2')); 
            if isempty(n)
                n = 1;
            end
        else
            n = 1;
        end

        options = session_info(n).probe(1);

        SUBJECT_experiment_info = subject_session_stimuli_mapping({options.SUBJECT},option);
        iDate = find([SUBJECT_experiment_info(:).date] == str2double(options.SESSION));
        
        % We also need RUN1 spatial cell IDs to filter
        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        if isempty(DIR1)
            continue;
        end
        
        session_count = session_count + 1;
        
        % Load RUN1 to get spatial_cell_id
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'), 'session_clusters');
        session_clusters_RUN = session_clusters;
        clear session_clusters;

        spatial_cell_id = session_clusters_RUN.cluster_id((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
            | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));

        % Load PRE
        try
            % Often it's 'session_clusters_SleepPRE.mat' or similar
            % If erase(stimulus_name{n},'Chronic') is used in previous script:
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        catch
            fprintf('Could not load session_clusters for PRE session %d\n', nsession);
            continue;
        end
        
        % From cell structure back to spike times and spike id
        session_clusters.spike_id_all=vertcat(session_clusters.spike_id{:});
        session_clusters.spike_times_all=vertcat(session_clusters.spike_times{:});
        [session_clusters.spike_times_all,index] = sort(session_clusters.spike_times_all);
        session_clusters.spike_id_all=session_clusters.spike_id_all(index);
        
        % Filter only spatial cells
        spike_index = ismember(session_clusters.spike_id_all, spatial_cell_id);
        spike_times = session_clusters.spike_times_all(spike_index);
        spike_id = session_clusters.spike_id_all(spike_index);
        
        [spike_times,sort_idx] = sort(spike_times);
        spike_id = spike_id(sort_idx);

        % Save lightweight data
        session_clusters_all.subject(session_count,:) = options.SUBJECT;
        session_clusters_all.session(session_count,:) = options.SESSION;
        session_clusters_all.session_day(session_count) = iDate;

        session_clusters_all.spike_id{nsession} = spike_id;
        session_clusters_all.spike_times{nsession} = spike_times;
        session_clusters_all.spatial_cell_id{nsession} = spatial_cell_id;
    end

    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end

    save(fullfile(analysis_folder,'session_clusters_all_PRE.mat'),'session_clusters_all','-v7.3')
    fprintf('Saved session_clusters_all_PRE.mat successfully!\n');
end
