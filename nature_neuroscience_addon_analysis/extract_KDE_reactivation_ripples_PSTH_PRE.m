function extract_KDE_reactivation_ripples_PSTH_PRE()

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic_PRE';

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end

load(fullfile(analysis_folder,'ripples_all_PRE.mat'),'ripples_all');

ripple_window = [-1, 1];
psth_step = 0.01;
nbins = round(diff(ripple_window) / psth_step);
metric_names = {'bias', 'zscored_bias', 'log_odds', 'zscored_log_odds'};
regions = ["HPC", "V1"];

ripple_bias_masked_SWS = struct();

% Initialize
for nprobe = 1:2
    for region = regions
        for metric = metric_names
            metric_str = metric{1};
            ripple_bias_masked_SWS.(region).(metric_str){nprobe} = [];
        end
    end
end
nan_mask = cell(1,2); 

session_count = 0;

for nsession = 1:length(experiment_info)
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    
    if length(stimulus_name) > 1
        n = find(contains(stimulus_name,'_2'));
        if isempty(n), n = 1; end
    else
        n = 1;
    end
    
    options = session_info(n).probe(1);
    
    pre_folder = fullfile(options.ANALYSIS_DATAPATH, '..', sprintf('session_clusters_%s.mat', erase(stimulus_name{n},'Chronic')));
    if ~exist(pre_folder, 'file') && ~exist(fullfile(options.ANALYSIS_DATAPATH, '..', sprintf('session_clusters_%s.mat', stimulus_name{n})), 'file')
        continue;
    end

    session_count = session_count + 1;
    fprintf('Processing session %d (count %d)\n', nsession, session_count);
    
    kde_hpc_file = fullfile(options.ANALYSIS_DATAPATH, 'KDE_reactivation_ripples.mat');
    kde_v1_file = fullfile(options.ANALYSIS_DATAPATH, 'KDE_reactivation_V1_ripples.mat');
    
    fprintf('Checking HPC file: %s\n', kde_hpc_file);
    has_hpc = exist(kde_hpc_file, 'file');
    has_v1 = exist(kde_v1_file, 'file');
    
    if has_hpc
        fprintf('Loading HPC file!\n');
        load(kde_hpc_file, 'KDE_reactivation_ripples');
        KDE_reactivation = KDE_reactivation_ripples;
    else
        fprintf('HPC file MISSING!\n');
        KDE_reactivation = [];
    end
    
    if has_v1
        load(kde_v1_file, 'KDE_reactivation_V1_ripples');
        KDE_reactivation_V1 = KDE_reactivation_V1_ripples;
    else
        KDE_reactivation_V1 = [];
    end

    for nprobe = 1:2
        % Find ripples for this session
        session_mask = ripples_all(nprobe).session_count == session_count;
        ripple_idx_all = find(session_mask); % absolute indices of this session
        state_mask = ripples_all(nprobe).SWS_index == 1;
        % IA (ripple_idx_session) gives the relative index 1, 2, 3... within this session's ripples
        [ripple_idx_abs, ripple_idx_session] = intersect(ripple_idx_all, find(state_mask));
        
        if isempty(ripple_idx_session)
            continue; 
        end
        
        % Ripple onset times for the SWS ripples in this session
        ripple_times = [ripples_all(nprobe).onset(ripple_idx_abs), ripples_all(nprobe).offset(ripple_idx_abs)];
        
        nan_mask_probe = [];
        
        for region = regions
            if region == "V1"
                src = KDE_reactivation_V1;
            else
                src = KDE_reactivation;
            end
            
            if isempty(src) || length(src) < nprobe || isempty(src(nprobe).event_bias)
                continue;
            end
            
            % Compute whole-session stats
            real_bias = src(nprobe).event_bias(:)';
            log_odds = log(src(nprobe).event_T1_probability ./ src(nprobe).event_T2_probability);
            log_odds = log_odds(:)';
            
            mean_bias = mean(real_bias(~isinf(real_bias)), 'omitnan');
            std_bias = std(real_bias(~isinf(real_bias)), 'omitnan');
            z_bias = (real_bias - mean_bias) / std_bias;

            % Keep +/-Inf (same as POST); capped at 0.5/99.5 prctile at analysis time
            mean_lo = mean(log_odds(~isinf(log_odds)), 'omitnan');
            std_lo = std(log_odds(~isinf(log_odds)), 'omitnan');
            z_lo = (log_odds - mean_lo) / std_lo;
            
            metric_vals = struct();
            metric_vals.bias = real_bias;
            metric_vals.zscored_bias = z_bias;
            metric_vals.log_odds = log_odds;
            metric_vals.zscored_log_odds = z_lo;
            
            time_bins = src(nprobe).event_bins(:,1)';
            event_ids = src(nprobe).event_id;
            
            for metric = metric_names
                metric_str = metric{1};
                metric_all_session = [];
                vals = metric_vals.(metric_str);
                
                for i = 1:length(ripple_idx_session)
                    ripple_id = ripple_idx_session(i);
                    t_start = ripple_times(i, 1);
                    rel_time = time_bins - t_start;
                    
                    % Find bins for this specific ripple within the -1 to 1 window
                    % We match only on event_ids since PRE timestamp offsets can be messy
                    idx = find(event_ids == ripple_id);
                    
                    ripple_vec = nan(1, nbins);
                    if ~isempty(idx)
                        if length(idx) > nbins
                            % If there are 201 bins (e.g. -1 to 1 inclusive), take first 200
                            idx = idx(1:nbins);
                        end
                        ripple_vec(1:length(idx)) = vals(idx);
                    end
                    metric_all_session = [metric_all_session; ripple_vec];
                    
                    % Compute masking only once per ripple (using HPC bias as trigger)
                    if strcmp(region, 'HPC') && strcmp(metric_str, 'bias')
                        nan_row = zeros(1, nbins);
                        
                        all_onsets = ripple_times(:,1);
                        all_offsets = ripple_times(:,2);
                        
                        % Mask next ripple
                        next_ripples = all_onsets > t_start & all_onsets < (t_start + 1);
                        if any(next_ripples)
                            t_next = min(all_onsets(next_ripples));
                            mask_start = round((t_next - t_start) / psth_step) + 101;
                            if mask_start <= nbins
                                nan_row(mask_start:end) = NaN;
                            end
                        end
                        
                        % Mask previous ripple
                        prev_ripples = all_offsets < t_start & all_offsets > (t_start - 1);
                        if any(prev_ripples)
                            t_prev = max(all_offsets(prev_ripples));
                            mask_end = round((t_prev - t_start) / psth_step) + 101;
                            if mask_end >= 1
                                nan_row(1:mask_end) = NaN;
                            end
                        end
                        
                        nan_mask_probe = [nan_mask_probe; nan_row];
                    end
                end
                
                % Append to master array for this probe
                ripple_bias_masked_SWS.(region).(metric_str){nprobe} = [ripple_bias_masked_SWS.(region).(metric_str){nprobe}; metric_all_session];
            end
        end
        
        nan_mask{nprobe} = [nan_mask{nprobe}; nan_mask_probe];
    end
end

% Assemble final struct
KDE_reactivation_ripples_PSTH_PRE = struct();
for region = regions
    r = char(region);
    KDE_reactivation_ripples_PSTH_PRE.([r '_ripples'])             = [ripple_bias_masked_SWS.(r).bias{1};             ripple_bias_masked_SWS.(r).bias{2}];
    KDE_reactivation_ripples_PSTH_PRE.([r '_z_ripples'])           = [ripple_bias_masked_SWS.(r).zscored_bias{1};     ripple_bias_masked_SWS.(r).zscored_bias{2}];
    KDE_reactivation_ripples_PSTH_PRE.([r '_logodds_ripples'])     = [ripple_bias_masked_SWS.(r).log_odds{1};         ripple_bias_masked_SWS.(r).log_odds{2}];
    KDE_reactivation_ripples_PSTH_PRE.([r '_z_logodds_ripples'])   = [ripple_bias_masked_SWS.(r).zscored_log_odds{1}; ripple_bias_masked_SWS.(r).zscored_log_odds{2}];
end

KDE_reactivation_ripples_PSTH_PRE.nan_mask = [nan_mask{1}; nan_mask{2}];

% Cap +/-Inf log odds at 99.5 / 0.5 percentile (same as POST)
KDE_reactivation_ripples_PSTH_PRE = cap_log_odds_PSTH(KDE_reactivation_ripples_PSTH_PRE);

% Save
out_dirs = {fullfile(analysis_folder, 'V1-HPC sleep reactivation'), ...
    fullfile(analysis_folder, 'nature_neuroscience_revision_analysis')};
for i = 1:length(out_dirs)
    if ~exist(out_dirs{i}, 'dir')
        mkdir(out_dirs{i});
    end
    save(fullfile(out_dirs{i}, 'KDE_reactivation_ripples_PSTH_PRE.mat'), 'KDE_reactivation_ripples_PSTH_PRE', '-v7.3');
    fprintf('Saved KDE_reactivation_ripples_PSTH_PRE.mat to %s\n', out_dirs{i});
end

end
