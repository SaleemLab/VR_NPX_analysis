addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);
%% Extract key information for DOWN UP, ripples and spindles
 
%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end


%%%% Grouping of DOWN-UP and UP-DOWN and ripples based on detetcion lags

%%%%%%%%%%%%% Overlaps
% Initialize output
L_overlap_idx = [];
R_overlap_idx = [];
overlap_idx = [];
non_overlap_idx = [];
merged_idx = [];
all_lags = [];
all_overlap_idx=[];
% for nsession = 1:max(slow_waves_all.UP_session_count)

%     L_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 1 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
%     R_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 2 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
event_types = {'UP','DOWN','ripples','spindles'};
for n = 1:4
 
    L_idx = find(merged_event_info.(sprintf('%s_hemisphere_id',event_types{n})) == 1);
    R_idx = find(merged_event_info.(sprintf('%s_hemisphere_id',event_types{n})) == 2);
    L_ints = merged_event_info.(sprintf('%s_ints',event_types{n}))(L_idx,:);
    R_ints = merged_event_info.(sprintf('%s_ints',event_types{n}))(R_idx,:);

    windows_threshold = [0.05,0.1,0.2];

    for ngroup = 1:4

        % Track indices into merged_event_info.UP_ints
        ipsi_overlap_idx{ngroup} = [];     % leading events involved in overlaps
        contra_overlap_idx{ngroup} = [];     % lagging events involved in overlaps
        unique_lags{ngroup} = [];
        all_lags{ngroup} = [];

        L_overlap_idx{ngroup} = [];
        R_overlap_idx{ngroup} =[];
        L_lags{ngroup} = [];
        R_lags{ngroup} =[];
        lags = [];

        all_overlap_idx{ngroup} = [];
        overlap_idx{ngroup} = [];
        non_overlap_idx{ngroup} = [];

        for iL = 1:size(L_ints,1)
            % Get current L interval

            L_start = L_ints(iL,1);
            if ngroup == 4
                L_end   = L_ints(iL,2);
            else
                L_end   = L_ints(iL,1)+windows_threshold(ngroup);
            end
            % Check overlap with all R intervals
            if ngroup == 4
                is_overlap = (R_ints(:,1) <= L_end) & (R_ints(:,2) >= L_start);
            else
                is_overlap = (R_ints(:,1) <= L_end) & (R_ints(:,1)+windows_threshold(ngroup) >= L_start);
            end

            if any(is_overlap)
                overlapping_R = find(is_overlap);

                for j = 1:length(overlapping_R)
                    iR = overlapping_R(j);

                    % Store overlapping pair indices
                    L_overlap_idx{ngroup}(end+1) = L_idx(iL);
                    R_overlap_idx{ngroup}(end+1) = R_idx(iR);
                    L_lags{ngroup}(end+1) = L_ints(iL,1) - R_ints(iR,1);
                    R_lags{ngroup}(end+1) = R_ints(iR,1) - L_ints(iL,1) ;

                end
            end
        end

        % Unique overlap winners
        lags = [L_lags{ngroup} R_lags{ngroup}];
        all_overlap_idx{ngroup} = [L_overlap_idx{ngroup} R_overlap_idx{ngroup}];
        [overlap_idx{ngroup},~] = unique(all_overlap_idx{ngroup});
        % unique_overlap_idx{ngroup} = lags(index);
        all_lags{ngroup} = lags;

        % Non-overlapping L and R
        nonoverlap_L_idx = setdiff(L_idx, unique(L_overlap_idx{ngroup}));
        nonoverlap_R_idx = setdiff(R_idx, unique(R_overlap_idx{ngroup}));

        non_overlap_idx{ngroup} = unique([nonoverlap_L_idx; nonoverlap_R_idx]);
    end


    merged_event_info.(sprintf('%s_overlap_idx_all',event_types{n})) = all_overlap_idx;
    % (sprintf('%s_lags_all',event_types{n}))
    merged_event_info.(sprintf('%s_non_overlap_idx',event_types{n})) = non_overlap_idx;
    merged_event_info. (sprintf('%s_overlap_idx',event_types{n})) = overlap_idx;
    merged_event_info.(sprintf('%s_lags_all',event_types{n})) = all_lags;
end

%%%% lag based on LFP
probability = probability_psth_whole;
nprobe = 1;
event_averaging_scale = 10;
% extract plv, amp corr and lag (latency)
for nprobe=1:2
    mprobe = abs(nprobe-3);

  
    % ripples
    ipsi_lag_ripples{nprobe} = [];
    contra_lag_ripples{nprobe} = [];

    ipsi_corr_ripples{nprobe} = [];
    contra_corr_ripples{nprobe} = [];

    ipsi_plv_ripples{nprobe} = [];
    contra_plv_ripples{nprobe} = [];

    for nsession = 1:max(ripples_all(1).session_count)

  
        %%%%% Ripples
        [C,ia,ib] = intersect(find(ripples_all(nprobe).session_count == sessions_to_process(nsession)),find(ripples_all(nprobe).session_count == sessions_to_process(nsession) & ripples_all(nprobe).SWS_index == 1));
        ipsi_shank = find(ripples_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==HPC_ref_shank(nsession,nprobe))=[];
        contra_shank = find(ripples_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);


        ipsi_lag_ripples{nprobe} = [ipsi_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_ripples{nprobe} = [contra_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_ripples{nprobe} = [ipsi_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_ripples{nprobe} = [contra_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_ripples{nprobe} = [ipsi_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_ripples{nprobe} = [contra_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];
    end
end


%%% Ripples and spindles features and UP DOWN features during UP/DOWN
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole.mat'));
spindle_probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% probability_normalised_whole = probability;


%%%%%%% Get DU and UP delta power
%%%%%%% DOWN and UP info
ipsi_Delta_peaks_zscore_UD = [];
ipsi_Delta_peaks_zscore_DU = [];
SWpeakmag_UD = [];
SWpeakmag_DU = [];

contra_Delta_peaks_zscore_UD = [];
contra_Delta_peaks_zscore_DU = [];

index = [];

for nprobe = 1:2
    for nsession = 1:max(slow_waves_all(1).UP_session_count)

        % slow_waves_all(1).SWpeakmag(slow_waves_all(1).UP_session_count == nsession)
        % Find DOWN
        [C,ia,ib]  =intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index);

        ipsi_Delta_peaks_zscore_UD = [ipsi_Delta_peaks_zscore_UD;...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,nprobe),ia)'];
        SWpeakmag_UD = [SWpeakmag_UD;...
            slow_waves_all(nprobe).SWpeakmag(C)];

        % Find DOWN before UP
        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        ipsi_Delta_peaks_zscore_DU = [ipsi_Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,nprobe),ia)'];
        index = [index; intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index)];

        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);

        SWpeakmag_DU = [SWpeakmag_DU;...
            slow_waves_all(nprobe).SWpeakmag(temp(ia))];



        % slow_waves_all(1).SWpeakmag(slow_waves_all(1).UP_session_count == nsession)
        % Find DOWN
        [C,ia,ib]  =intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index);

        contra_Delta_peaks_zscore_UD = [contra_Delta_peaks_zscore_UD;...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,abs(nprobe-3)),ia)'];

        % Find DOWN before UP
        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        contra_Delta_peaks_zscore_DU = [contra_Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,abs(nprobe-3)),ia)'];
        index = [index; intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession),probability(nprobe).DOWN_all_index)];

        [C,ia,ib]  = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession,2), slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession),probability(nprobe).UP_all_index),1));

        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
    end
end

% Calculate ripples probability according to unilateral vs bilateral events
%%%%% Ripple info
ripples_times =  merged_event_info.ripples_ints;
ripples_group_info = merged_event_info.ripples_group_id;
ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id;
ripples_original_index = [find(ripples_all(1).SWS_index); find(ripples_all(2).SWS_index)];


ripples_lag_diff = [contra_lag_ripples{1} contra_lag_ripples{2}]';
ripples_all_overlap_idx = merged_event_info.ripples_overlap_idx_all{end};
ripples_non_overlap_idx = merged_event_info.ripples_non_overlap_idx{end};
ripples_all_lags =merged_event_info.ripples_lags_all{end};

% merged_event_info.ripples_index_sorted

%%%%%%%%%%%%%% Grab ripple info during UP
hemi_labels = {'L', 'R'};
region_labels = {'L_V1','R_V1','L_HPC','R_HPC'};
varnames = {
    'time_from_last_ripples'
    'first_ripples_power'
    'last_ripples_power'
    'first_ripples_lag'
    'last_ripples_lag'
    'first_ripples_lag_diff'
    'last_ripples_lag_diff'
    'ripple_counts'
    'ripples_duration'
    'first_ripples_duration'
    'last_ripples_duration'
};

% Init data_struct
for h = 1:2
    hemi = hemi_labels{h};
    for v = 1:length(varnames)
        data_struct.(hemi).(varnames{v}) = cell(1,2);
    end
    for r = 1:length(region_labels)
        reg = region_labels{r};
        data_struct.(hemi).(['ripple_' reg '_MUA_cumulative']) = cell(1,2);
        data_struct.(hemi).(['normalised_ripple_' reg '_MUA_cumulative']) = cell(1,2);
        data_struct.(hemi).(['first_ripple_' reg '_MUA_peak']) = cell(1,2);
        data_struct.(hemi).(['last_ripple_' reg '_MUA_peak']) = cell(1,2);
    end
end

for r = 1:length(region_labels)
    reg = region_labels{r};
    data_struct.([reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['normalised_' reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['first_half_' reg '_MUA_cumulative']) = cell(1,2);
    data_struct.(['second_half_' reg '_MUA_cumulative']) = cell(1,2);
end


varnames = fieldnames(data_struct.L);

% Main loop
for nprobe = 1:2
    UP_indices = probability_psth_whole(nprobe).UP_all_index;
    for h = 1:2
        hemi = hemi_labels{h};
        ripple_idx = h;
        ripple_data = ripples_all(ripple_idx);
        ripple_norm_dur = event_info(nprobe).([hemi '_ripple_normalised_UP_duration']);

        for nevent = 1:length(UP_indices)
            up_idx = UP_indices(nevent);

            % --- Always store UP-wide MUA metrics ---
            up_duration = event_info(nprobe).UP_duration(nevent);
            for r = 1:length(region_labels)
                reg = region_labels{r};
                trace = event_info(nprobe).([reg '_MUA_UP']){nevent};
                npoints = length(trace);
                half_idx = floor(npoints/2);
                cum_sum = sum(trace);
                data_struct.([reg '_MUA_cumulative']){nprobe}(nevent) = cum_sum;
                data_struct.(['normalised_' reg '_MUA_cumulative']){nprobe}(nevent) = cum_sum / up_duration;
                data_struct.(['first_half_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(trace(1:half_idx)) / (up_duration/2);
                data_struct.(['second_half_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(trace(half_idx+1:end)) / (up_duration/2);
            end

            % --- Ripple-based metrics (only if ripples exist in this UP) ---
            event_index = find(up_idx == ripple_norm_dur(:,2));
            if ~isempty(event_index)
                ripples_index = ripple_norm_dur(event_index, 1);
                onset = ripple_data.onset(ripples_index);
                offset = ripple_data.offset(ripples_index);
                duration = offset - onset;

                data_struct.(hemi).first_ripples_duration{nprobe}(nevent) = offset(1) - onset(1);
                data_struct.(hemi).last_ripples_duration{nprobe}(nevent) = offset(end) - onset(end);
                data_struct.(hemi).ripples_duration{nprobe}(nevent) = sum(duration);
                data_struct.(hemi).ripple_counts{nprobe}(nevent) = length(ripples_index);
                data_struct.(hemi).time_from_last_ripples{nprobe}(nevent) = ...
                    slow_waves_all(nprobe).UP_ints(up_idx,2) - ripple_data.peaktimes(ripples_index(end));

                data_struct.(hemi).first_ripples_power{nprobe}(nevent) = ripple_data.peak_zscore(ripples_index(1));
                data_struct.(hemi).last_ripples_power{nprobe}(nevent) = ripple_data.peak_zscore(ripples_index(end));

                idx1 = find(ripples_original_index == ripples_index(1) & ripples_hemisphere_id == ripple_idx);
                idx2 = find(ripples_original_index == ripples_index(end) & ripples_hemisphere_id == ripple_idx);
                data_struct.(hemi).first_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(idx1);
                data_struct.(hemi).last_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(idx2);

                [C, ia] = intersect(ripples_all_overlap_idx, idx1);
                if ~isempty(C)
                    data_struct.(hemi).first_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
                end
                [C, ia] = intersect(ripples_all_overlap_idx, idx2);
                if ~isempty(C)
                    data_struct.(hemi).last_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
                end

                % MUA from event_info
                for r = 1:length(region_labels)
                    reg = region_labels{r};  % 'L_V1', 'R_HPC', etc.
                    region_hemi = strcmp(reg(1), 'R') + 1;  % 1 = L, 2 = R
                    region_type = reg(3:end);              % 'V1' or 'HPC'

                    data_struct.(hemi).(['ripple_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(event_info(nprobe).([hemi '_ripple_' region_type '_MUA_cumulative_UP'])(event_index, region_hemi));
                    data_struct.(hemi).(['normalised_ripple_' reg '_MUA_cumulative']){nprobe}(nevent) = sum(event_info(nprobe).([hemi '_ripple_' region_type '_MUA_cumulative_UP'])(event_index, region_hemi)) ./ sum(duration);

                    data_struct.(hemi).(['first_ripple_' reg '_MUA_peak']){nprobe}(nevent) = event_info(nprobe).([hemi '_ripple_' region_type '_MUA_peak_UP'])(event_index(1), region_hemi);
                    data_struct.(hemi).(['last_ripple_' reg '_MUA_peak']){nprobe}(nevent) = event_info(nprobe).([hemi '_ripple_' region_type '_MUA_peak_UP'])(event_index(end), region_hemi);
                end
            else
                data_struct.(hemi).time_from_last_ripples{nprobe}(nevent) = 0;
                for v = 1:length(varnames)
                    data_struct.(hemi).(varnames{v}){nprobe}(nevent) = nan;
                end
                for r = 1:length(region_labels)
                    reg = region_labels{r};
                    for suffix = {'ripple_', 'normalised_ripple_', 'first_ripple_', 'last_ripple_'}
                        fname = [suffix{1} reg '_MUA_cumulative'];
                        if contains(fname, 'peak')
                            fname = strrep(fname, '_cumulative', '_peak');
                        end
                        data_struct.(hemi).(fname){nprobe}(nevent) = nan;
                    end
                end
            end
        end
    end
end

% Build UP_DOWN_info
UP_DOWN_info = struct();
varnames = fieldnames(data_struct.L);
for i = 1:length(varnames)
    varname = varnames{i};
    if ~contains(varname,'MUA')
        L_var = data_struct.L.(varname);
        R_var = data_struct.R.(varname);
        UP_DOWN_info.(['ipsi_' varname '_UP']) = [L_var{1}, R_var{2}];
        UP_DOWN_info.(['contra_' varname '_UP']) = [R_var{1}, L_var{2}];
    end
end

% ripple_<hemi>_<region>_MUA_<type>
ripple_regions = {'V1', 'HPC'};
ripple_mua_suffixes = {
    'ripple_%s_MUA_cumulative'
    'normalised_ripple_%s_MUA_cumulative'
    'first_ripple_%s_MUA_peak'
    'last_ripple_%s_MUA_peak'
    'first_ripple_%s_MUA_cumulative'
    'last_ripple_%s_MUA_cumulative'
};


for i = 1:length(ripple_regions)
    region = ripple_regions{i};
    for j = 1:length(ripple_mua_suffixes)
        suffix_template = ripple_mua_suffixes{j};

        % Build field names for L and R
        L_field = sprintf(suffix_template, ['L_' region]);
        R_field = sprintf(suffix_template, ['R_' region]);

        % Build output variable name
        var_base = strrep(suffix_template, '%s', region);
        ipsi_name = ['ipsi_' var_base '_UP'];
        contra_name = ['contra_' var_base '_UP'];

        % Assign
        UP_DOWN_info.(ipsi_name) = [data_struct.L.(L_field){1}, data_struct.R.(R_field){2}];
        UP_DOWN_info.(contra_name) = [data_struct.R.(R_field){1}, data_struct.L.(L_field){2}];
    end
end


% <prefix>_<hemi>_<region>_MUA_<type>
global_regions = {'V1', 'HPC'};
global_prefixes = {'', 'normalised_', 'first_half_', 'second_half_'};

for i = 1:length(global_regions)
    region = global_regions{i};
    for j = 1:length(global_prefixes)
        prefix = global_prefixes{j};

        L_field = [prefix 'L_' region '_MUA_cumulative'];
        R_field = [prefix 'R_' region '_MUA_cumulative'];

        base_name = [prefix region '_MUA_cumulative_UP'];
        UP_DOWN_info.(['ipsi_' base_name]) = [data_struct.(L_field){1}, data_struct.(R_field){2}];
        UP_DOWN_info.(['contra_' base_name]) = [data_struct.(R_field){1}, data_struct.(L_field){2}];
    end
end


% Cleanup
clear hemi_labels region_labels varnames data_struct
clear h r v f i C ia ripple_idx
clear hemi reg suffix fname field
clear nprobe nevent up_idx up_duration npoints half_idx cum_sum trace
clear UP_indices ripple_norm_dur ripple_data
clear event_index ripples_index onset offset duration vals MUA_cum MUA_peak
clear idx1 idx2
clear ripple_regions ripple_types global_regions global_prefixes
clear region prefix base_name L_field R_field ipsi_name contra_name



%%%%% Grab DOWN/UP info
varnames = {
    'ipsi_Delta_peaks_zscore_UD'
    'ipsi_Delta_peaks_zscore_DU'
    'SWpeakmag_UD'
    'SWpeakmag_DU'
    'contra_Delta_peaks_zscore_UD'
    'contra_Delta_peaks_zscore_DU'
    };

% Loop through each variable name
for i = 1:length(varnames)
    varname = varnames{i};
    temp =eval(varname);
    UP_DOWN_info.(varname) = temp;
end
UP_DOWN_info.UP_duration = [event_info(1).UP_duration;event_info(2).UP_duration];
UP_DOWN_info.previous_DOWN_duration = [event_info(1).previous_DOWN_duration;event_info(2).previous_DOWN_duration];
UP_DOWN_info.next_DOWN_duration = [event_info(1).next_DOWN_duration;event_info(2).next_DOWN_duration];


%%%%%%%%%% Spindles
%%%%%%%%%%%%%%%%%%%%
UP_DOWN_info.ipsi_spindles_UP = [spindle_probability_psth_whole(1).L_spindles_UP; spindle_probability_psth_whole(2).R_spindles_UP];
UP_DOWN_info.contra_spindles_UP = [spindle_probability_psth_whole(1).R_spindles_UP; spindle_probability_psth_whole(2).L_spindles_UP];
UP_DOWN_info.ipsi_spindles_DOWN = [spindle_probability_psth_whole(1).L_spindles_DOWN; spindle_probability_psth_whole(2).R_spindles_DOWN];
UP_DOWN_info.contra_spindles_DOWN = [spindle_probability_psth_whole(1).R_spindles_DOWN; spindle_probability_psth_whole(2).L_spindles_DOWN];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% ALL
bins_to_use = bin_centers>0 & bin_centers<0.1;
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);


[bV1_mean_T1, LCI1_V1, UCI1_V1] = calculate_bootstrap_CI(z_bias_V1', T1_events);
[bV1_mean_T2, LCI2_V1, UCI2_V1] = calculate_bootstrap_CI(z_bias_V1', T2_events);
[b_mean_T1, LCI1, UCI1] = calculate_bootstrap_CI(z_bias', T1_events);
[b_mean_T2, LCI2, UCI2] = calculate_bootstrap_CI(z_bias', T2_events);


nfig = figure;
nfig.Name = 'Bayesian Bias PSTH Track 1 vs Track 2 All'; 
nfig.Position = [1070 110 800 860];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

% Colors
track1_color = [0 0 1]; % blue
track2_color = [1 0 0]; % red

% Tile 5: z_bias mean with CI

nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI1 fliplr(LCI1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2 fliplr(LCI2)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
legend('Track 1', 'Track 2','box','off');
xlabel('Time (s)')
ylabel('Bayesian Bias (z)')
title('z\_bias CI');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color, 'LineWidth', 2);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color, 'LineWidth', 2);
legend('Track 1', 'Track 2','box','off');
title('z\_bias\_V1 CI');
xlabel('Time (s)')
ylabel('Bayesian Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
