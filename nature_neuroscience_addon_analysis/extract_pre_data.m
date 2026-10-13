function extract_pre_data()
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    
    out_dir = fullfile(analysis_folder, 'V1-HPC sleep reactivation');
    
    fprintf('Loading ripples_all_PRE.mat and slow_waves_all_PRE.mat...\n');
    load(fullfile(analysis_folder, 'ripples_all_PRE.mat'), 'ripples_all');
    load(fullfile(analysis_folder, 'slow_waves_all_PRE.mat'), 'slow_waves_all');
    
    % Calculate cortex_ref_shank for PRE
    cortex_ref_shank = zeros(max(ripples_all(1).session_count), 2);
    for nsession = 1:max(ripples_all(1).session_count)
        for probe_no = 1:2
            % Some PRE sessions might not have slow_waves if LFP was bad or missing, let's be careful
            if length(slow_waves_all(probe_no).channel) >= nsession && ~isempty(slow_waves_all(probe_no).channel{nsession})
                ch_idx = slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession);
                best_shank = slow_waves_all(probe_no).shank{nsession}(ch_idx);
                cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == best_shank(1), 1);
            else
                % Default to 1 if no slow_waves exist for this session/probe
                cortex_ref_shank(nsession,probe_no) = 1;
            end
        end
    end

    ripple_info_PRE = struct();
    
    % Get power and indices
    r_pow_1 = ripples_all(1).peak_zscore(ripples_all(1).SWS_index == 1);
    r_pow_2 = ripples_all(2).peak_zscore(ripples_all(2).SWS_index == 1);
    ripple_info_PRE.ripple_power = [r_pow_1; r_pow_2];
    
    ripple_info_PRE.event_id = [find(ripples_all(1).SWS_index == 1); find(ripples_all(2).SWS_index == 1)];
    ripple_info_PRE.session_id = [ripples_all(1).session_count(ripples_all(1).SWS_index == 1); ripples_all(2).session_count(ripples_all(2).SWS_index == 1)];

    for nprobe = 1:2
        all_SO_phase = [];
        all_spindle_amp = [];
        all_spindle_pct = [];
        
        for nsession = 1:max(ripples_all(nprobe).session_count)
            % Wait, PRE only has 12 sessions in ripples_all!
            % We mapped PRE sessions into 1..12 in ripples_all.session_count.
            if length(ripples_all(nprobe).SO_phase_ripple_peaktime) >= nsession && ~isempty(ripples_all(nprobe).SO_phase_ripple_peaktime{nsession})
                shank_idx = cortex_ref_shank(nsession,:);
                % If the data doesn't have 2 shanks (e.g. unilateral), guard against index out of bounds
                max_shank = size(ripples_all(nprobe).SO_phase_ripple_peaktime{nsession}, 1);
                shank_idx(shank_idx > max_shank) = max_shank;
                shank_idx(shank_idx < 1) = 1;
                
                all_SO_phase = [all_SO_phase, ripples_all(nprobe).SO_phase_ripple_peaktime{nsession}(shank_idx,:)];
                all_spindle_amp = [all_spindle_amp, ripples_all(nprobe).spindle_amplitude_ripple_peaktime{nsession}(shank_idx,:)];
                
                if isfield(ripples_all(nprobe), 'spindle_percentile') && length(ripples_all(nprobe).spindle_percentile) >= nsession && ~isempty(ripples_all(nprobe).spindle_percentile{nsession})
                    all_spindle_pct = [all_spindle_pct, ripples_all(nprobe).spindle_percentile{nsession}(shank_idx,:)];
                else
                    all_spindle_pct = [all_spindle_pct, nan(2, size(ripples_all(nprobe).SO_phase_ripple_peaktime{nsession}, 2))];
                end
            end
        end
        
        mask = ripples_all(nprobe).SWS_index == 1;
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
    
    ripple_info_PRE.SO_phase = [SO_phase_1; SO_phase_2];
    ripple_info_PRE.spindle_amplitude = [spindle_amp_1; spindle_amp_2];
    ripple_info_PRE.spindle_percentile = [spindle_pct_1; spindle_pct_2];
    
    save(fullfile(out_dir, 'ripple_info_PRE.mat'), 'ripple_info_PRE', '-v7.3');
    fprintf('Saved ripple_info_PRE.mat\n');

    % Copy files to revision folder
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
    if ~exist(rev_dir, 'dir')
        mkdir(rev_dir);
    end
    
    files_to_copy = {
        'ripples_all_PRE.mat',
        fullfile('V1-HPC sleep reactivation', 'ripple_info_PRE.mat')
    };

    fprintf('Copying PRE files to revision directory...\n');
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
