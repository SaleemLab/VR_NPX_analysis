function extract_awake_data()
    % Define analysis folder
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    
    out_dir = fullfile(analysis_folder, 'V1-HPC sleep reactivation');
    
    % --- 1. Create KDE_reactivation_ripples_PSTH_awake.mat ---
    fprintf('Loading KDE_reactivation_ripples_PSTH.mat...\n');
    load(fullfile(out_dir, 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');
    
    KDE_reactivation_ripples_PSTH_awake = struct();
    regions = {'HPC', 'V1'};
    for r = 1:2
        region = regions{r};
        KDE_reactivation_ripples_PSTH_awake.([region '_ripples'])             = KDE_reactivation_ripples_PSTH.([region '_awake_ripples']);
        KDE_reactivation_ripples_PSTH_awake.([region '_z_ripples'])           = KDE_reactivation_ripples_PSTH.([region '_z_awake_ripples']);
        KDE_reactivation_ripples_PSTH_awake.([region '_logodds_ripples'])     = KDE_reactivation_ripples_PSTH.([region '_logodds_awake_ripples']);
        KDE_reactivation_ripples_PSTH_awake.([region '_z_logodds_ripples'])   = KDE_reactivation_ripples_PSTH.([region '_z_logodds_awake_ripples']);
    end
    KDE_reactivation_ripples_PSTH_awake.nan_mask = KDE_reactivation_ripples_PSTH.nan_mask_awake;
    
    save(fullfile(out_dir, 'KDE_reactivation_ripples_PSTH_awake.mat'), 'KDE_reactivation_ripples_PSTH_awake', '-v7.3');
    fprintf('Saved KDE_reactivation_ripples_PSTH_awake.mat\n');

    % --- 2. Create ripple_info_awake.mat ---
    fprintf('Loading ripples_all_POST.mat to extract ripple_info_awake...\n');
    % We need cortex_ref_shank for SO_phase indexing. Let's load it from workspace or just assume best channel.
    % Actually we can just load ripples_all_best_V1_SO_POST.mat or ripples_all_POST.mat
    try
        load(fullfile(analysis_folder, 'ripples_all_POST.mat'), 'ripples_all');
    catch
        load(fullfile(analysis_folder, 'ripples_all.mat'), 'ripples_all');
    end
    
    % The user wants SO phase L and R, spindle power L and R (also percentile) and ripple power L and R
    % For awake ripples: ripples_all.awake_index == 1
    % Wait, `awake_index` might be `awake_index==1` or `SWS_index==0`. Let's use awake_index==1
    
    % Find cortex ref shank for each session (often it's just 1, but let's check)
    % Usually we take the first column or the best_channel.
    % Let's see if SO_phase_ripple_peaktime is already 1D per session or 2D (shank x ripple).
    % In the original code: `ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)`
    % If we just take the first row for each session, since we don't have cortex_ref_shank here,
    % wait, the user's `ripple_info` was built in `corticohippocampal_reactivation_distribution.m` which defined `cortex_ref_shank`.
    % Let's load the original `ripple_info.mat` to see if we can get the dimensions right.
    % Actually, if I look at `corticohippocampal_reactivation_distribution.m`, `cortex_ref_shank` was computed.
    % To be safe, I'll extract it using the first row (shank 1) or best_channel.
    % But `ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}` is a cell array.
    % Let's inspect it dynamically.

    % We also need to define `event_ids_first`. In `corticohippocampal_reactivation_distribution.m`:
    % `[~, event_ids_first, ~] = unique(merged_event_info.ripples_peaktimes);`
    % But wait! If we don't need to merge UP/DOWN, we just want a flat list of awake ripples for L and R?
    % The user said: "i only need SO phase L and R and spindle power L and R (also percentile) and ripple power L and R. ripple_info_awake.mat"
    % This means `ripple_info_awake.SO_phase` should be `[N x 2]` where N is the number of awake ripples.
    % But wait! Do we want to align to HPC ripples or merged ripples?
    % The PSTH `KDE_reactivation_ripples_PSTH_awake.HPC_ripples` has dimensions `[N1+N2 x 200]`.
    % It's literally just concatenated probe 1 and probe 2 ripples!
    % `KDE_reactivation_ripples_PSTH.([region '_awake_ripples']) = [ripple_bias_masked_nonSWS.(region).bias{1}; ripple_bias_masked_nonSWS.(region).bias{2}];`
    % So `ripple_info_awake` should have the EXACT SAME rows as `KDE_reactivation_ripples_PSTH_awake.HPC_ripples`.
    % So it's just Probe 1 awake ripples followed by Probe 2 awake ripples.
    
    ripple_info_awake = struct();
    
    % Determine which ripples were used for PSTH
    % In PSTH: `state_mask = ripples_all(nprobe).SWS_index == 0;`
    % Let's use SWS_index == 0 to match exactly.
    
    % Ripple power
    r_pow_1 = ripples_all(1).peak_zscore(ripples_all(1).SWS_index == 0);
    r_pow_2 = ripples_all(2).peak_zscore(ripples_all(2).SWS_index == 0);
    ripple_info_awake.ripple_power = [r_pow_1; r_pow_2];
    
    % Event ID
    event_id_1 = find(ripples_all(1).SWS_index == 0);
    event_id_2 = find(ripples_all(2).SWS_index == 0);
    ripple_info_awake.event_id = [event_id_1; event_id_2]; % Absolute index in ripples_all
    
    % Subject & Session
    sess_1 = ripples_all(1).session_count(ripples_all(1).SWS_index == 0);
    sess_2 = ripples_all(2).session_count(ripples_all(2).SWS_index == 0);
    ripple_info_awake.session_id = [sess_1; sess_2];
    
    % For SO phase and spindle, it's stored in `{nsession}`
    SO_phase_1 = []; SO_phase_2 = [];
    spindle_amp_1 = []; spindle_amp_2 = [];
    spindle_pct_1 = []; spindle_pct_2 = [];
    
    % Load slow_waves_all to calculate cortex_ref_shank
    load(fullfile(analysis_folder, 'slow_waves_all_POST.mat'), 'slow_waves_all');
    cortex_ref_shank = zeros(max(ripples_all(1).session_count), 2);
    for nsession = 1:max(ripples_all(1).session_count)
        for probe_no = 1:2
            ch_idx = slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession);
            best_shank = slow_waves_all(probe_no).shank{nsession}(ch_idx);
            cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == best_shank(1), 1);
        end
    end

    for nprobe = 1:2
        all_SO_phase = [];
        all_spindle_amp = [];
        all_spindle_pct = [];
        
        for nsession = 1:max(ripples_all(nprobe).session_count)
            all_SO_phase = [all_SO_phase, ripples_all(nprobe).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
            all_spindle_amp = [all_spindle_amp, ripples_all(nprobe).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
            if isfield(ripples_all(nprobe), 'spindle_percentile') && length(ripples_all(nprobe).spindle_percentile) >= nsession
                all_spindle_pct = [all_spindle_pct, ripples_all(nprobe).spindle_percentile{nsession}(cortex_ref_shank(nsession,:),:)];
            else
                all_spindle_pct = [all_spindle_pct, nan(2, size(ripples_all(nprobe).SO_phase_ripple_peaktime{nsession}, 2))];
            end
        end
        
        mask = ripples_all(nprobe).SWS_index == 0;
        if nprobe == 1
            SO_phase_1 = all_SO_phase(:, mask)';
            spindle_amp_1 = all_spindle_amp(:, mask)';
            spindle_pct_1 = all_spindle_pct(:, mask)';
        else
            SO_phase_2 = all_SO_phase(:, mask)';
            spindle_amp_2 = all_spindle_amp(:, mask)';
            spindle_pct_2 = all_spindle_pct(:, mask)';
        end
    end
    
    ripple_info_awake.SO_phase = [SO_phase_1; SO_phase_2];
    ripple_info_awake.spindle_amplitude = [spindle_amp_1; spindle_amp_2];
    ripple_info_awake.spindle_percentile = [spindle_pct_1; spindle_pct_2];
    
    save(fullfile(out_dir, 'ripple_info_awake.mat'), 'ripple_info_awake', '-v7.3');
    fprintf('Saved ripple_info_awake.mat\n');

    % --- 3. Copy main variables to D:\corticohippocampal_replay\nature_neuroscience_revision_analysis ---
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    if ~exist(rev_dir, 'dir')
        mkdir(rev_dir);
    end
    fprintf('Copying files to revision directory...\n');
    
    files_to_copy = {
        'session_clusters_all_POST.mat',
        'session_clusters_all_RUN1.mat',
        'session_clusters_all_PRE.mat',
        'ripples_all_POST.mat',
        fullfile('V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'),
        fullfile('V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH_PRE.mat'),
        fullfile('V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH_awake.mat'),
        fullfile('V1-HPC sleep reactivation', 'ripple_info.mat'),
        fullfile('V1-HPC sleep reactivation', 'ripple_info_awake.mat')
    };

    for i = 1:length(files_to_copy)
        src_file = fullfile(analysis_folder, files_to_copy{i});
        [~, name, ext] = fileparts(src_file);
        dst_file = fullfile(rev_dir, [name ext]);
        if exist(src_file, 'file')
            copyfile(src_file, dst_file);
            fprintf('Copied %s\n', name);
        else
            fprintf('Warning: File not found %s\n', src_file);
        end
    end
    
    fprintf('All done!\n');
end
