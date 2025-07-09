%%%%%% Main script for mixed effect and cox regression analysis to examine
%%%%%% predictive relationship between ripples and DOWN UP states and
%%%%%% spindles

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
% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
% load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));
probability_psth = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;
% %%%%%%%
% % Phase-amplitude coupling
% nBins = 18;
% edges = linspace(-pi, pi, nBins+1);
% 
% for nchannel = 1:size(SO_phase_ripples,1)
%     for mchannel = 1:size(ripple_peak_amplitude,1)
%         % Phase-amplitude coupling
%         SO_phase_ripples(:,nevent) = ripple_peak_amplitude(tidx,:);
%         spindle_phase_ripples(:,nevent) = ripple_peak_amplitude(tidx,:);
% 
% 
%         [~,~,binIdx] = histcounts(SO_phase_ripples(nchannel,:), edges);
%         % Mean amplitude in each phase bin
%         ampByPhase = accumarray(binIdx(binIdx>0), rippleAmp(binIdx>0), [nBins 1], @mean);
%         p = ampByPhase / sum(ampByPhase); % normalize to get probability distribution
%     end
% end
% if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
%     mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
% end
% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])


% SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
% SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));


% Find reference channel/shank
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


%% Grouping of DOWN-UP and UP-DOWN and ripples based on detetcion lags
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded.mat'),'MUA_PSTH_merged_thresholded');


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


%% Ripples and spindles features and UP DOWN features during UP/DOWN
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_probability_whole.mat'));
spindle_probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% probability_normalised_whole = probability;


% Get DU and UP delta power
ipsi_Delta_peaks_zscore_UD = [];
ipsi_Delta_peaks_zscore_DU = [];
SWpeakmag_UD = [];
SWpeakmag_DU = [];

contra_Delta_peaks_zscore_UD = [];
contra_Delta_peaks_zscore_DU = [];

index = [];
%%%%% UP info
% UP_index_all = [probability(1).UP_index; probability(2).UP_index];
% DOWN_index_all = [probability(1).DOWN_index; probability(2).DOWN_index];

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
    % 'ripple_%s_MUA_cumulative'
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

%% Predict Ripples from DOWN-UP info
%%%%% UP info
event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.UP_overlap_idx_all{end};
non_overlap_idx = merged_event_info.UP_non_overlap_idx{end};
lags =merged_event_info.UP_lags_all{end};



UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%%%%%%%%%%
%%%%% V1 MUA
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA = [PSTH_MUA(1).R_V1_UP; PSTH_MUA(2).L_V1_UP];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];

%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];

ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];

%%%%%%%%%%% ripples
ipsi_probability = [probability_psth_whole(1).L_ripples_UP; probability_psth_whole(2).R_ripples_UP];
contra_probability = [probability_psth_whole(1).R_ripples_UP; probability_psth_whole(2).L_ripples_UP];

%%%%%%%%%% Predict HPC MUA and ripples during DOWN-UP transition based on
%%%%%%%%%% V1 DOWN UP bilateral synchrony and magnitude
%%%%%%%%% 0 to 100ms ripples
% 
% index = all_overlap_idx(abs(lags)<=0.15);
% lag_index = (abs(lags)<=0.15);
% output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_output.mat'),'output');
% 
% 
% 
% output = predict_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA, contra_HPC_MUA,ipsi_probability, contra_probability,UP_DOWN_info,'subject_id',subject_id);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_output.mat'),'output');
% 
% 
% %%%%%%%%%%%%%%%%%%% Plotting
% 
index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);

load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_output.mat'),'output');
output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),...
    ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'output',output,'plot_option',1);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'),[],'ContentType','image')

% 
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_output.mat'),'output');
output = predict_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA, contra_HPC_MUA,ipsi_probability, contra_probability,UP_DOWN_info,'subject_id',subject_id,'output',output,'plot_option',1);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'),[],'ContentType','image')






%%%%%%%%%% 100 - 200ms ripples
%%%%%%%%%% Predict HPC MUA and ripples during DOWN-UP transition based on
%%%%%%%%%% V1 DOWN UP bilateral synchrony and magnitude

index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),...
    ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),...
    'time_window', [0.15 0.25]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_250ms_output.mat'),'output');



output = predict_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA, contra_HPC_MUA,ipsi_probability, contra_probability,UP_DOWN_info,'subject_id',subject_id,...
    'time_window', [0.15 0.25]);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_250ms_output.mat'),'output');


%%%%%%%%%%%%%%%%%%% Plotting

index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_output.mat'),'output');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_250ms_output.mat'),'output');
output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),...
    ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'output',output,'plot_option',1,'time_window',[0.15 0.25]);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression 250ms ripples (150ms windows)'),[],'ContentType','image')

% 
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_output.mat'),'output');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_250ms_output.mat'),'output');
output = predict_ripples_by_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA, contra_HPC_MUA,ipsi_probability, contra_probability,UP_DOWN_info,'subject_id',subject_id,'output',output,'plot_option',1,'time_window',[0.15 0.25]);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression 250ms ripples (full windows)'),[],'ContentType','image')


%%


%%%%% DOWN info
event_times = merged_event_info.DOWN_ints;
hemisphere_id = merged_event_info.DOWN_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.DOWN_overlap_idx_all{end};
non_overlap_idx = merged_event_info.DOWN_non_overlap_idx{end};
lags =merged_event_info.DOWN_lags_all{end};

DOWN_session_count = [slow_waves_all(1).DOWN_session_count(probability(1).DOWN_all_index); slow_waves_all(2).DOWN_session_count(probability(2).DOWN_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(DOWN_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%%%%%%%%%%%%% Grab ripple info during UP DOWN
%%%%%%%%%%%
%%%%% V1 MUA
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_V1_MUA = [PSTH_MUA(1).R_V1_DOWN; PSTH_MUA(2).L_V1_DOWN];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_DOWN; PSTH_MUA_baseline(2).R_V1_DOWN];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_DOWN; PSTH_MUA_baseline(2).L_V1_DOWN];

%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_DOWN; PSTH_MUA(2).R_HPC_DOWN];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_DOWN; PSTH_MUA(2).L_HPC_DOWN];

ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_DOWN; PSTH_MUA_baseline(2).R_HPC_DOWN];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_DOWN; PSTH_MUA_baseline(2).L_HPC_DOWN];

%%%%%%%%%%% ripples
ipsi_probability = [probability_psth_whole(1).L_ripples_DOWN; probability_psth_whole(2).R_ripples_DOWN];
contra_probability = [probability_psth_whole(1).R_ripples_DOWN; probability_psth_whole(2).L_ripples_DOWN];

ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_DOWN; probability_psth_whole_baseline(2).R_ripples_DOWN];
contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_DOWN; probability_psth_whole_baseline(2).L_ripples_DOWN];

% 
% %%%%%%%%%%%%%%%% -0.1 to 0 ripples
% index = 1:size(ipsi_V1_MUA,1);
% % lag_index = (abs(lags)<=2);
% output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'subject_id',subject_id(index));
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
% 
% 
% index = all_overlap_idx(abs(lags)<=0.15);
% lag_index = (abs(lags)<=0.15);
% output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
% 
% 
% 
% %%%%%%%%%%%% Plot
index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'subject_id',subject_id(index));
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'),[],'ContentType','image')


index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'),[],'ContentType','image')


% 



%%%%%%%%%%%%%%%% -0.2 to 0 ripples
index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'subject_id',subject_id(index),'time_window', -0.2);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_200ms_output.mat'),'output');


index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),...
    UP_DOWN_info,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'time_window', -0.2);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_200ms_output.mat'),'output');



%%%%%%%%%%%% Plot
index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'subject_id',subject_id(index),'time_window', -0.2);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression 200ms ripples (full windows)'),[],'ContentType','image')


index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
    'output',output,'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'time_window', -0.2);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression 200ms ripples (150ms windows)'),[],'ContentType','image')


%%%%%%%%%%%% Control for lag and delta power by using perceding UP
%%%%%%%%%%%% transition lag and delta power
UP_DOWN_info.SWpeakmag_UD1 = UP_DOWN_info.SWpeakmag_UD ;
UP_DOWN_info.SWpeakmag_UD = UP_DOWN_info.SWpeakmag_DU;

all_overlap_idx_UP = merged_event_info.UP_overlap_idx_all{end};
lags =merged_event_info.UP_lags_all{end};
index = all_overlap_idx_UP(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);

output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
    'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output_based_on_UP.mat'),'output');


index = 1:size(ipsi_V1_MUA,1);
output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'subject_id',subject_id(index),'time_window', -0.1);
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_200ms_output.mat'),'output');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output_based_on_UP.mat'),'output');


index = 1:size(ipsi_V1_MUA,1);
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output_based_on_UP.mat'),'output');
output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,'subject_id',subject_id(index),'time_window', -0.1,'plot_option',1,'output',output);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows UP lag control)'),[],'ContentType','image')



lags =merged_event_info.UP_lags_all{end};
index = all_overlap_idx_UP(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output_based_on_UP.mat'),'output');
output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
    'plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index),'output',output);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows UP lag control)'),[],'ContentType','image')




%% Generate CSV files

model_indices = [1 2 7 8 17 18 5 6 11 12 21 22 3 4 9 10 19 20];
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression 250ms ripples (full windows)'));
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_250ms_output.mat'),'output');
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'DU_ripples_Summary.csv');


model_indices = [1 2 7 8 17 18 5 6 11 12 21 22 3 4 9 10 19 20];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_Ripples_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'DU_ripples_Summary.csv');


model_indices = [1 2 9 10 25 26 3 4 11 12 27 28 7 8 15 16 31 32 5 6 13 14 29 30];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','DU_synchrony_Ripples_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'DU_synchrony Ripples_Summary.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output_based_on_UP.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows UP lag control)'));
model_indices = [13 14 29 30 5 6 2 3 25 26 18 19 21 22 32 33 7 23 15];
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1synchrony_Summary.csv');


model_indices = [13 14 29 30 5 6 2 3 25 26 18 19 21 22 32 33 7 23 15];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1synchrony_Summary.csv');

% [10 11 13 14 21 22 5 6 2 3 17 18 7 15]
model_indices = [10 11 13 14 21 22 5 6 2 3 17 18 7 15];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output_based_on_UP.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows UP lag control)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1depression_Summary.csv');


model_indices = [10 11 13 14 21 22 5 6 2 3 17 18 7 15];
load(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');
cd(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'));
% load('Ripples_V1synchrony_output.mat');  % Load your 'output' variable
summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, 'Ripples_V1depression_Summary.csv');



% output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),UP_DOWN_info,...
%     plot_option',1,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output_based_on_UP.mat'),'output');

%% Effect of ripples and cumulative HPC activities  on UP survival probability
UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);
UP_DOWN_info.subject_id = subject_id;

%%%%%%%%%%%%% Ripple power predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[], ...
    'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition','subject_id',subject_id...
);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis'),[])

% scatter([UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP], ...
%     [UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP],...
%     'filled','MarkerFaceAlpha','0.02')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% ylim([0 25])
% xlim([0 0.3])


% histogram(UP_DOWN_info.ipsi_last_ripples_power_UP)
% plot_UP_survival_probability( ...
%     {UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP}, ...
%     {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',1,'subject_id',subject_id, ...
%     'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition (1 ripple)'...
% );

% 
% % Ripple MUA peak predicts UP probability
% plot_UP_survival_probability( ...
%     {UP_DOWN_info.ipsi_last_ripple_V1_MUA_peak_UP, UP_DOWN_info.contra_last_ripple_V1_MUA_peak_UP}, ...
%     {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
%     'event_option',[],'title_name', 'Last ripple V1 peak MUA and last ripple to UP-DOWN transition'...
% );


%%%%%%%%%%%%% cumulative ripple HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_ripple_HPC_MUA_cumulative_UP, UP_DOWN_info.contra_ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Ripple cumulative HPC MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);


% cumulative ripple HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_HPC_MUA_cumulative_UP, UP_DOWN_info.contra_normalised_ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Ripple normalised cumulative HPC MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);

% cumulative ripple V1 MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_ripple_V1_MUA_cumulative_UP, UP_DOWN_info.contra_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Ripple cumulative V1 MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);


% cumulative ripple HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP, UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Ripple normalised cumulative V1 MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);

%%%%%%%%%%%%%%%%% ripple bilateral synchrony does not predict UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripples_lag_diff_UP, UP_DOWN_info.contra_last_ripples_lag_diff_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[], ...
    'event_option',[],'title_name', 'Last ripple synchrony and last ripple to UP-DOWN transition','subject_id',subject_id...
);



%%%%%%%%%%%%%%%% Normalised cumulative HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_HPC_MUA_cumulative_UP, UP_DOWN_info.contra_normalised_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);


% Normalsied cumulative V1 MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_V1_MUA_cumulative_UP, UP_DOWN_info.contra_normalised_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);

%%%%%%%%%%%%%%%% 1st half Normalised cumulative HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_first_half_HPC_MUA_cumulative_UP, UP_DOWN_info.contra_first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);


% 1st half Normalsied cumulative V1 MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 1st and last ripple to UP-DOWN transition','timebin', 0.015...
);

% 2nd half Normalised cumulative HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_HPC_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd half and last ripple to UP-DOWN transition','timebin', 0.015...
);


% 2nd half Normalsied cumulative V1 MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd and last ripple to UP-DOWN transition','timebin', 0.015...
);

%%%%%%%%%%%%%%% 2nd half - 1st half

% 2nd half Normalised cumulative HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.ipsi_first_half_HPC_MUA_cumulative_UP,...
    UP_DOWN_info.contra_second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.contra_first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd-1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);


% 2nd half Normalsied cumulative V1 MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP-UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP,...
    UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP-UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd-1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis'),[])


%%%%%%%%%%%%%%%% cumulative HPC MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_HPC_MUA_cumulative_UP, UP_DOWN_info.contra_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative HPC MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);


% cumulative V1 MUA predicts UP probability
plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_V1_MUA_cumulative_UP, UP_DOWN_info.contra_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative V1 MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis'),[])

%% Ripple predicted by spindle and SO phase and amplitude


%%%%%%%%%%%%%%%%%% Ripple info (based in ipsi or contra)
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'));
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%% ripple power
ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];

%%% spindle co-occurance
[~,spindle_index,~,index] =RestrictInts(merged_event_info.ripples_ints,merged_event_info.spindles_ints);
ripple_info.spindle_presence = spindle_index;
ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.spindle_presence_hemi(find(spindle_index)) = merged_event_info.spindles_hemisphere_id(index);

%%% spindle phase
spindle_phase=[];
spindle_phase1=[];
for probe_no = 1:2
    spindle_phase1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase1{probe_no} = [spindle_phase1{probe_no} ripples_all(probe_no).spindle_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end


spindle_phase(1,:) = [spindle_phase1{1}(1,ripples_all(1).SWS_index==1) spindle_phase1{2}(2,ripples_all(2).SWS_index==1)];
spindle_phase(2,:) = [spindle_phase1{1}(2,ripples_all(1).SWS_index==1) spindle_phase1{2}(1,ripples_all(2).SWS_index==1)];


ripple_info.spindle_phase = spindle_phase';

%%% spindle power
spindle_amplitude=[];
spindle_amplitude1=[];
for probe_no = 1:2
    spindle_amplitude1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_amplitude1{probe_no} = [spindle_amplitude1{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_amplitude(1,:) = [spindle_amplitude1{1}(1,ripples_all(1).SWS_index==1) spindle_amplitude1{2}(2,ripples_all(2).SWS_index==1)];
spindle_amplitude(2,:) = [spindle_amplitude1{1}(2,ripples_all(1).SWS_index==1) spindle_amplitude1{2}(1,ripples_all(2).SWS_index==1)];

ripple_info.spindle_amplitude = spindle_amplitude';

%%% SO phase
SO_phase=[];
SO_phase1=[];
for probe_no = 1:2
    SO_phase1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase1{probe_no} = [SO_phase1{probe_no} ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_phase(1,:) = [SO_phase1{1}(1,ripples_all(1).SWS_index==1) SO_phase1{2}(2,ripples_all(2).SWS_index==1)];
SO_phase(2,:) = [SO_phase1{1}(2,ripples_all(1).SWS_index==1) SO_phase1{2}(1,ripples_all(2).SWS_index==1)];
ripple_info.SO_phase = SO_phase';

%%% SO power
SO_amplitude=[];
SO_amplitude1=[];
for probe_no = 1:2
    SO_amplitude1{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_amplitude1{probe_no} = [SO_amplitude1{probe_no} ripples_all(probe_no).SO_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_amplitude(1,:) = [SO_amplitude1{1}(1,ripples_all(1).SWS_index==1) SO_amplitude1{2}(2,ripples_all(2).SWS_index==1)];
SO_amplitude(2,:) = [SO_amplitude1{1}(2,ripples_all(1).SWS_index==1) SO_amplitude1{2}(1,ripples_all(2).SWS_index==1)];

ripple_info.SO_amplitude = SO_amplitude';


%%% early UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,1)+0.2]);
ripple_info.early_UP_index = UP_index;
ripple_info.early_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.early_UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);


%%% late UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1)-0.2 merged_event_info.DOWN_ints(:,1)]);
ripple_info.late_UP_index = UP_index;
ripple_info.late_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.late_UP_index_hemi(find(UP_index)) = merged_event_info.DOWN_hemisphere_id(index);

%%% UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,2)]);
ripple_info.UP_index = UP_index;
ripple_info.UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);

%%% DOWN transition co-occurance
[~,DOWN_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1) merged_event_info.DOWN_ints(:,2)]);
ripple_info.DOWN_index = DOWN_index;
ripple_info.DOWN_index_hemi = zeros(size(DOWN_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.DOWN_index_hemi(find(DOWN_index)) = merged_event_info.DOWN_hemisphere_id(index);


%%%%%% Ripple session count and subject_id
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));
singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));


%%%%%%%%%%%%% Mixed effect model (Does spindle and SO predict ripple power?)

output = predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_output.mat'),'output');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_output.mat'),'output')
predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id,'output',output,'plot_option',1);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')



%%%%%%%%%%%%% UP transition
output = predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id);
save(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_UP_transition_output.mat'),'output');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_spindles_amp_phase_Ripples_DOWN_transition_output.mat'),'output')
predict_ripples_by_SO_spindles(ripple_info,'subject_id',subject_id,'output',output,'plot_option',1);

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


predict_ripples_by_SO_spindles_TF(ripple)




% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis'),[])
