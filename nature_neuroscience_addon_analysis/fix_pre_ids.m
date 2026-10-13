function fix_pre_ids()
    addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
    SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
    option = 'bilateral';
    experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
    experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
    
    % Find which of the 22 sessions have PRE
    pre_map = zeros(1, 12);
    subj_map = zeros(1, 12);
    pre_count = 0;
    
    for nsession = 1:22
        stims = experiment_info(nsession).StimulusName;
        if sum(contains(stims, 'SleepChronic_PRE')) > 0
            pre_count = pre_count + 1;
            pre_map(pre_count) = nsession;
            
            % Get subject string (e.g. 'M24062') and convert to unique 1-6 id
            subj_str = experiment_info(nsession).session(1).probe(1).SUBJECT;
            subj_num = str2double(subj_str(end-1:end));
            subj_map(pre_count) = subj_num;
        end
    end
    
    % The unique subject nums in experiment_info are (16, 17, 18, 62, 64, 65).
    % In POST ripple_info, they are mapped to 1:6 via `unique`.
    % Let's do the exact same `unique` mapping so they match perfectly!
    all_subj_nums = zeros(1, 22);
    for n = 1:22
        subj_str = experiment_info(n).session(1).probe(1).SUBJECT;
        all_subj_nums(n) = str2double(subj_str(end-1:end));
    end
    [~, ~, unique_subj_ids] = unique(all_subj_nums);
    
    % Now pre_subj_ids is the 1..6 mapping for the PRE sessions
    pre_subj_ids = unique_subj_ids(pre_map);
    
    analysis_folder = 'D:\corticohippocampal_replay';
    
    % Update ripple_info_PRE.mat
    fprintf('Updating ripple_info_PRE.mat...\n');
    f_info = fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'ripple_info_PRE.mat');
    load(f_info, 'ripple_info_PRE');
    
    % Map session_id
    old_session_ids = ripple_info_PRE.session_id;
    new_session_ids = zeros(size(old_session_ids));
    new_subject_ids = zeros(size(old_session_ids));
    
    for i = 1:12
        idx = (old_session_ids == i);
        new_session_ids(idx) = pre_map(i);
        new_subject_ids(idx) = pre_subj_ids(i);
    end
    
    ripple_info_PRE.session_id = new_session_ids;
    ripple_info_PRE.subject_id = new_subject_ids;
    
    save(f_info, 'ripple_info_PRE', '-v7.3');
    
    % Copy to revision folder
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    copyfile(f_info, fullfile(rev_dir, 'ripple_info_PRE.mat'));
    
    fprintf('Fixed ripple_info_PRE.mat session and subject IDs!\n');
    
    % Check if ripples_all_PRE needs the same? 
    % The user said "and ripples_all_PRE should also be copied to the same folder" (I already did that).
    % But they also said: "did i extract LFP already for pre? if so, maybe worth having ripple_info_PRE too and ripples_all_PRE should also be copied to the same folder" - this implies they just want them available. 
    % But "can you make sure that session id and animal id is refered to 22 session with POST"
    % Let's also add a mapped field to ripples_all_PRE just in case, WITHOUT breaking original session_count.
    
    fprintf('Updating ripples_all_PRE.mat...\n');
    f_rip = fullfile(analysis_folder, 'ripples_all_PRE.mat');
    load(f_rip, 'ripples_all');
    
    for nprobe = 1:2
        old_counts = ripples_all(nprobe).session_count;
        new_counts = zeros(size(old_counts));
        new_subjs = zeros(size(old_counts));
        for i = 1:12
            idx = (old_counts == i);
            new_counts(idx) = pre_map(i);
            new_subjs(idx) = pre_subj_ids(i);
        end
        % Add NEW fields so we don't break downstream scripts that expect 1:12 max index
        ripples_all(nprobe).mapped_session_id = new_counts;
        ripples_all(nprobe).mapped_subject_id = new_subjs;
    end
    
    save(f_rip, 'ripples_all', '-v7.3');
    copyfile(f_rip, fullfile(rev_dir, 'ripples_all_PRE.mat'));
    
    fprintf('All done!\n');
end
