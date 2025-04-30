function plot_ipsi_contra_MUA_synchronisation_detection_overlap_control(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,event_info,sessions_to_process)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
% load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
% load(fullfile(analysis_folder,'ripples_all_POST.mat'))
% load(fullfile(analysis_folder,'spindles_all_POST.mat'))
% load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'));
% probability_SO_SO = probability;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'));
% probability_SO_SO_contralateral = probability;
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'));
probability_SO_SO_contralateral = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'));
probability_SO_SO = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability.mat'));
probability_SO_SO_contralateral_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability_whole.mat'));
probability_SO_SO_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability_whole.mat'));
probability_ripples_ripples_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;

colour_lines = [];
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red for R
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue for L
colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral

colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

colour_lines = [44,123,182;215,25,28]/256;


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


%%%%%
%%%%%
%%%%%
%%%%%
%% Grabbing ipsi and contra values
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged.mat'),'MUA_PSTH_merged');

MUA_PSTH_merged_overlap_control = MUA_PSTH_merged;


probability = probability_psth_whole;
nprobe = 1;
event_averaging_scale = 10;
% extract plv, amp corr and lag (latency)
for nprobe=1:2
    mprobe = abs(nprobe-3);

    % UP and DOWN

    ipsi_amp_UD{nprobe} = [];
    contra_amp_UD{nprobe} = [];

    ipsi_lag_DU{nprobe} = [];
    contra_lag_DU{nprobe} = [];
    ipsi_lag_UD{nprobe} = [];
    contra_lag_UD{nprobe} = [];

    ipsi_corr_DU{nprobe} = [];
    contra_corr_DU{nprobe} = [];
    ipsi_corr_UD{nprobe} = [];
    contra_corr_UD{nprobe} = [];

    ipsi_plv_DU{nprobe} = [];
    contra_plv_DU{nprobe} = [];
    ipsi_plv_UD{nprobe} = [];
    contra_plv_UD{nprobe} = [];

    % ripples
    ipsi_lag_ripples{nprobe} = [];
    contra_lag_ripples{nprobe} = [];

    ipsi_corr_ripples{nprobe} = [];
    contra_corr_ripples{nprobe} = [];

    ipsi_plv_ripples{nprobe} = [];
    contra_plv_ripples{nprobe} = [];

    % spindles
    ipsi_lag_spindles{nprobe} = [];
    contra_lag_spindles{nprobe} = [];

    ipsi_corr_spindles{nprobe} = [];
    contra_corr_spindles{nprobe} = [];

    ipsi_plv_spindles{nprobe} = [];
    contra_plv_spindles{nprobe} = [];

    for nsession = 1:max(ripples_all(1).session_count)

        %%%%%% DOWN -> UP
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession)),probability(nprobe).UP_all_index);

        ipsi_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        contra_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);

        ipsi_lag_DU{nprobe} = [ipsi_lag_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_DU{nprobe} = [contra_lag_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_DU{nprobe} = [ipsi_plv_DU{nprobe} squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_DU{nprobe} = [contra_plv_DU{nprobe} squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_DU{nprobe} = [ipsi_corr_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_DU{nprobe} = [contra_corr_DU{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        % ipsi_lag_DU{nprobe} = [ipsi_lag_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_lag_DU{nprobe} = [contra_lag_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_plv_DU{nprobe} = [ipsi_plv_DU{nprobe} min(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_plv_DU{nprobe} = [contra_plv_DU{nprobe} min(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_corr_DU{nprobe} = [ipsi_corr_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_corr_DU{nprobe} = [contra_corr_DU{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        %%%%%% UP -> DOWN
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession)),probability(nprobe).DOWN_all_index);

        ipsi_amp_UD{nprobe} = [ipsi_amp_UD{nprobe} squeeze(slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(ipsi_shank,ia))];
        contra_amp_UD{nprobe} = [contra_amp_UD{nprobe} squeeze(slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(contra_shank,ia))];

        ipsi_lag_UD{nprobe} = [ipsi_lag_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_lag_UD{nprobe} = [contra_lag_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_plv_UD{nprobe} = [ipsi_plv_UD{nprobe} squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_UD{nprobe} = [contra_plv_UD{nprobe} squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_UD{nprobe} = [ipsi_corr_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_UD{nprobe} = [contra_corr_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

        % ipsi_lag_UD{nprobe} = [ipsi_lag_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % contra_lag_UD{nprobe} = [contra_lag_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_plv_UD{nprobe} = [ipsi_plv_UD{nprobe} min(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_plv_UD{nprobe} = [contra_plv_UD{nprobe} min(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % 
        % ipsi_corr_UD{nprobe} = [ipsi_corr_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        % contra_corr_UD{nprobe} = [contra_corr_UD{nprobe} min(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        %%%%% Spindles
        [C,ia,ib] = intersect(find(spindles_all(nprobe).session_count == sessions_to_process(nsession)),find(spindles_all(nprobe).session_count == sessions_to_process(nsession) & spindles_all(nprobe).SWS_index == 1));

        if ~isempty(ia)
            ipsi_lag_spindles{nprobe} = [ipsi_lag_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            contra_lag_spindles{nprobe} = [contra_lag_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

            ipsi_plv_spindles{nprobe} = [ipsi_plv_spindles{nprobe} squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            contra_plv_spindles{nprobe} = [contra_plv_spindles{nprobe} squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];

            ipsi_corr_spindles{nprobe} = [ipsi_corr_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            contra_corr_spindles{nprobe} = [contra_corr_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];
        end

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

%% Clustering of DOWN-UP transition lag vs corr vs plv distribution based on overlapping events
% customColors = [
%     0.698, 0.094, 0.169;  % Deep Red
%     0.957, 0.647, 0.510;  % Light Red-Orange
%     0.529, 0.529, 0.529;  % Gray
%     0.102, 0.102, 0.102;  % Near Black
%     ];


%%%%%%%%%%%%% Overlaps
% Initialize output
L_overlap_idx = [];
R_overlap_idx = [];
overlap_idx = [];
non_overlap_idx = [];
merged_idx = [];
% for nsession = 1:max(slow_waves_all.UP_session_count)

%     L_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 1 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
%     R_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 2 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
L_idx = find(merged_event_info.UP_hemisphere_id == 1);
R_idx = find(merged_event_info.UP_hemisphere_id == 2);
L_ints = merged_event_info.UP_ints(L_idx,:);
R_ints = merged_event_info.UP_ints(R_idx,:);

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

%% Plot DU transition averaged MUA
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;


%%%%% V1 MUA
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
contra_V1_MUA = [PSTH_MUA(1).R_V1_UP; PSTH_MUA(2).L_V1_UP];
ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_UP; PSTH_MUA_baseline(2).R_V1_UP];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_UP; PSTH_MUA_baseline(2).L_V1_UP];

lag_diff = merged_event_info.UP_lag_diff;
plv_diff = merged_event_info.UP_plv_diff;
corr_diff = merged_event_info.UP_corr_diff;
ipsi_lag = [ipsi_lag_DU{1} ipsi_lag_DU{2}]';
contra_lag = [contra_lag_DU{1} contra_lag_DU{2}]';
group_id = merged_event_info.UP_group_id;

sync_threshold = mean(abs([event_info(1).UP_lag_threshold_low event_info(2).UP_lag_threshold_low event_info(1).UP_lag_threshold_high event_info(2).UP_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];



%%%%%%%%%%%%% Grouping events
event_idx = [];
% Ipsi-contra UP_DOWN V1 MUA by three clusters
% [index,ia] = intersect(all_overlap_idx{end},merged_event_info.UP_index);
index = intersect(merged_event_info.UP_index_sorted,all_overlap_idx{end});

event_idx{1} = {intersect(find((group_id == 3 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),index),...
    intersect(find(group_id == 2 & lag_diff < 0.1 & lag_diff > -0.1 ),index),intersect(find((group_id == 1 & lag_diff < -0.1)|(group_id == 1 & lag_diff > 0.1) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by ipsi dominant cluster (different lag diff)
event_idx{2} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),index),intersect(find(group_id == 3 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff <0),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),index),...
    intersect(find(group_id == 3 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 3& lag_diff >0),50) ),index),intersect(find(group_id == 3 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff >0),50) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by contra dominant clusters (different lag diff)
event_idx{3} = {intersect(find(group_id == 1 & lag_diff <prctile(lag_diff(group_id == 1& lag_diff <0),50) ),index),intersect(find(group_id == 1 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff <0),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),index),...
    intersect(find(group_id == 1 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 1& lag_diff >0),50) ),index),intersect(find(group_id == 1 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by cluster 1 and 3 merged (different lag diff)
event_idx{4} = {intersect(find(group_id ~= 2 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),index),intersect(find(group_id ~= 2 &  lag_diff <0 &lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),index),...
    intersect(find(group_id == 2 & lag_diff <sync_threshold & lag_diff >-sync_threshold),index),...
    intersect(find(group_id ~= 2 & lag_diff>0 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),index),intersect(find(group_id ~= 2 & lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),index)}

% Ipsi-contra UP_DOWN V1 MUA by ipsi leading (different corr diff)
index = intersect(find(lag_diff <prctile(lag_diff(group_id ~= 2& lag_diff <0),50)),intersect(merged_event_info.UP_index_sorted,all_overlap_idx{end}));

event_idx{5} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index_sorted),...
    intersect(find(group_id == 1& corr_diff <prctile(corr_diff(group_id == 1),50) ),index),intersect(find(group_id == 1& corr_diff >prctile(corr_diff(group_id == 1),50) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by contra leading (different corr diff)
index = intersect(find(lag_diff <prctile(lag_diff(group_id ~= 2& lag_diff <0),50)),intersect(merged_event_info.UP_index_sorted,all_overlap_idx{end}));

event_idx{6} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index_sorted),...
    intersect(find(group_id == 1& corr_diff <prctile(corr_diff(group_id == 1),50) ),index),intersect(find(group_id == 1 & corr_diff >prctile(corr_diff(group_id == 1),50) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by 3 clusters leading only+
index = intersect(merged_event_info.UP_index_sorted,all_overlap_idx{end});

event_idx{7} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),index),...
    intersect(find(group_id == 1 & lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50)),index)};


group_name=[];
group_name{1} = {'Ipsi dominant','Bilaterally synchronised','Contra dominant','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% ipsi dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% contra dominant clusters
group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra DOWN_UP MUA by three clusters (overlapping)', 'Ipsi-contra DOWN_UP MUA by ipsi dominant cluster (different lag diff overlapping)','Ipsi-contra DOWN_UP MUA by contra dominant clusters (different lag diff overlapping)',...
    'Ipsi-contra DOWN_UP MUA by cluster 1 and 3 merged (different lag diff overlapping)','Ipsi-contra DOWN_UP MUA by ipsi leading (diff corr diff overlapping)','Ipsi-contra DOWN_UP MUA by contra leading (diff corr diff overlapping)','Ipsi-contra DOWN_UP MUA by three clusters (leading only overlapping)'}

% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged_overlap_control.x = x;
MUA_PSTH_merged_overlap_control.ipsi_DU_V1 = [];
MUA_PSTH_merged_overlap_control.contra_DU_V1 = [];
MUA_PSTH_merged_overlap_control.ipsi_contra_diff_DU_V1 = [];

MUA_PSTH_merged_overlap_control.ipsi_DU_HPC = [];
MUA_PSTH_merged_overlap_control.contra_DU_HPC = [];
MUA_PSTH_merged_overlap_control.ipsi_contra_diff_DU_HPC = [];

for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged_overlap_control.ipsi_DU_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_overlap_control.contra_DU_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_overlap_control.ipsi_contra_diff_DU_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_overlap_control.ipsi_DU_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_overlap_control.contra_DU_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_overlap_control.ipsi_contra_diff_DU_HPC{ngroup}{i} = temp6;
    end
end


MUA_PSTH_merged_overlap_control.DU_groups = [title_names {'all DOWN_UP'}];
MUA_PSTH_merged_overlap_control.DU_index = [event_idx (1:size(ipsi_V1_MUA,1))'];

%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1|ngroup ==7
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_DU_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_DU_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-1 0.5])
    title('V1 Ipsi MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_DU_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_DU_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_DU_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_DU_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-1 0.5])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_DU_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_DU_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.1])
    title('HPC Ipsi MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_DU_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_DU_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.1])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_DU_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_DU_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.1])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end



%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])






%% Plot UD transition averaged MUA
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_event_info.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

ipsi_V1_MUA = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
contra_V1_MUA = [PSTH_MUA(1).R_V1_DOWN; PSTH_MUA(2).L_V1_DOWN];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_DOWN; PSTH_MUA_baseline(2).R_V1_DOWN];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_DOWN; PSTH_MUA_baseline(2).L_V1_DOWN];

% ipsi_probability = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
% contra_probability = [PSTH_MUA(2).L_V1_DOWN; PSTH_MUA(1).R_V1_DOWN];

lag_diff = merged_event_info.DOWN_lag_diff;
plv_diff = merged_event_info.DOWN_plv_diff;
corr_diff = merged_event_info.DOWN_corr_diff;
ipsi_lag = [ipsi_lag_UD{1} ipsi_lag_UD{2}]';
contra_lag = [contra_lag_UD{1} contra_lag_UD{2}]';
group_id = merged_event_info.DOWN_group_id;

sync_threshold = mean(abs([event_info(1).DOWN_lag_threshold_low event_info(2).DOWN_lag_threshold_low event_info(1).DOWN_lag_threshold_high event_info(2).DOWN_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Contra
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;

%%%%%%%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_DOWN; PSTH_MUA(2).R_HPC_DOWN];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_DOWN; PSTH_MUA(2).L_HPC_DOWN];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_DOWN; PSTH_MUA_baseline(2).R_HPC_DOWN];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_DOWN; PSTH_MUA_baseline(2).L_HPC_DOWN];

%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap1 = temp;

%%% Contra
binnedArray = contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap1 = temp;

%%% Ipsi-contra difference
binnedArray = ipsi_HPC_MUA_baseline-contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap1 = temp;

%%%%%%% Grouping events
event_idx = [];
% Ipsi-contra UP_DOWN V1 MUA by three clusters
event_idx{1} = {intersect(find((group_id == 3 & lag_diff < -0.1)|(group_id == 3 & lag_diff > 0.1) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < 0.1 & lag_diff > -0.1 ),merged_event_info.DOWN_index),intersect(find((group_id == 1 & lag_diff < -0.1)|(group_id == 1 & lag_diff > 0.1) ),merged_event_info.DOWN_index)};

% Ipsi-contra UP_DOWN V1 MUA by ipsi dominant cluster (different lag diff)
event_idx{2} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 3 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 3 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 3 &lag_diff >prctile(lag_diff(group_id == 3& lag_diff >0),50) ),merged_event_info.DOWN_index)};

% Ipsi-contra UP_DOWN V1 MUA by contra dominant clusters (different lag diff)
event_idx{3} = {intersect(find(group_id == 1 & lag_diff <prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 1 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 1 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.DOWN_index),intersect(find(group_id == 1 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50) ),merged_event_info.DOWN_index)};

% Ipsi-contra UP_DOWN V1 MUA by cluster 1 and 3 merged (different lag diff)
event_idx{4} = {intersect(find(group_id ~= 2 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.DOWN_index),intersect(find(group_id ~= 2 &  lag_diff <0 &lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff <sync_threshold & lag_diff >-sync_threshold),merged_event_info.DOWN_index),...
    intersect(find(group_id ~= 2 & lag_diff>0 & lag_diff <prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.DOWN_index),intersect(find(group_id ~= 2 & lag_diff >prctile(lag_diff(group_id ~= 2 & lag_diff >0),50) ),merged_event_info.DOWN_index)};

% Ipsi-contra UP_DOWN V1 MUA by ipsi leading (different corr diff)
index = intersect(find(lag_diff <prctile(lag_diff(group_id ~= 2& lag_diff <0),50)),merged_event_info.DOWN_index);

event_idx{5} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1& corr_diff <prctile(corr_diff(group_id == 1),50) ),index),intersect(find(group_id == 1& corr_diff >prctile(corr_diff(group_id == 1),50) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by contra leading (different corr diff)
index = intersect(find(lag_diff >prctile(lag_diff(group_id ~= 2& lag_diff >0),50)),merged_event_info.DOWN_index);

event_idx{6} = {intersect(find(group_id == 3 & corr_diff >prctile(corr_diff(group_id == 3),50) ),index),intersect(find(group_id == 3 & corr_diff <prctile(corr_diff(group_id == 3),50) ),index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.UP_index),...
    intersect(find(group_id == 1& corr_diff <prctile(corr_diff(group_id == 1),50) ),index),intersect(find(group_id == 1 & corr_diff >prctile(corr_diff(group_id == 1),50) ),index)};

% Ipsi-contra UP_DOWN V1 MUA by 5 clusters
event_idx{7} = {intersect(find(group_id == 3 & lag_diff <prctile(lag_diff(group_id == 3& lag_diff <0),50) ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 2 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),merged_event_info.DOWN_index),...
    intersect(find(group_id == 1 & lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50)),merged_event_info.DOWN_index)};



group_name=[];
group_name{1} = {'Ipsi dominant','Bilaterally synchronised','Contra dominant','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% ipsi dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% contra dominant clusters
group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra UP_DOWN MUA by three clusters', 'Ipsi-contra UP_DOWN MUA by ipsi dominant cluster (different lag diff)','Ipsi-contra UP_DOWN MUA by contra dominant clusters (different lag diff)',...
    'Ipsi-contra UP_DOWN MUA by cluster 1 and 3 merged (different lag diff)','Ipsi-contra UP_DOWN MUA by ipsi leading (diff corr diff)','Ipsi-contra UP_DOWN MUA by contra leading (diff corr diff)','Ipsi-contra UP_DOWN MUA by three clusters (leading only)'}
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged.x = x;
MUA_PSTH_merged.ipsi_UD_V1 = [];
MUA_PSTH_merged.contra_UD_V1 = [];
MUA_PSTH_merged.ipsi_contra_diff_UD_V1 = [];

MUA_PSTH_merged.ipsi_UD_HPC = [];
MUA_PSTH_merged.contra_UD_HPC = [];
MUA_PSTH_merged.ipsi_contra_diff_UD_HPC = [];

for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged.ipsi_UD_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged.contra_UD_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged.ipsi_contra_diff_UD_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged.ipsi_UD_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged.contra_UD_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged.ipsi_contra_diff_UD_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged.ipsi_UD_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged.contra_UD_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged.ipsi_contra_diff_UD_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged.ipsi_UD_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged.contra_UD_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged.ipsi_contra_diff_UD_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged.UD_groups = [title_names {'all UP_DOWN'}];
MUA_PSTH_merged.UD_index = [event_idx (1:size(ipsi_V1_MUA,1))'];


%%%%%%%%% Plot UD transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1|ngroup ==7
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_UD_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_UD_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-1 0.5])
    title('V1 Ipsi MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_UD_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_UD_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.5])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-1 0.5])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_UD_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_UD_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.2 0.3])
    title('HPC Ipsi MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_UD_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_UD_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
ylim([-0.2 0.3])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.2 0.3])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end


%%%%%%%%% All DOWN-UP plots

clear ERROR_SHADE
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name ='All ipsi-contra UP-DOWN MUA';

colour_lines = [0,90,50;74,20,134]/256; % dark purple, meganta, light green, dark green

nexttile
binnedArray = MUA_PSTH_merged.ipsi_UD_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)


binnedArray = MUA_PSTH_merged.contra_UD_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged.ipsi_UD_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-1 0.5])


% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('V1 MUA')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline

binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-1 0.5])

% xline(0,'r')
legend([ERROR_SHADE(1:2)],{'ipsi-contra diff','shuffled'},'box','off')
title('V1 MUA difference')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('ipsi-contra MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
binnedArray = MUA_PSTH_merged.ipsi_UD_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

binnedArray = MUA_PSTH_merged.contra_UD_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged.ipsi_UD_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.2 0.3])

% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('HPC MUA')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray =  MUA_PSTH_merged.ipsi_contra_diff_UD_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline

binnedArray = MUA_PSTH_merged.ipsi_contra_diff_UD_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.2 0.3])

% xline(0,'r')
legend([ERROR_SHADE(1:2)],{'ipsi-contra diff','shuffled'},'box','off')
title('HPC MUA difference')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('ipsi-contra MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])

%%
%%%%%
%%%%%
%%%%%
%%%%%

%% Plot ripples averaged MUA
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_event_info.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;
%%%%%%%%%%%%%%%% HPC spikes
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_ripples; PSTH_MUA(2).R_HPC_ripples];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_ripples; PSTH_MUA(2).L_HPC_ripples];

ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_ripples; PSTH_MUA_baseline(2).R_HPC_ripples];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(2).L_HPC_ripples; PSTH_MUA_baseline(1).R_HPC_ripples];
% 
% ipsi_probability = [PSTH_MUA(1).L_V1_UP; PSTH_MUA(2).R_V1_UP];
% contra_probability = [PSTH_MUA(2).L_V1_UP; PSTH_MUA(1).R_V1_UP];

lag_diff = merged_event_info.ripples_lag_diff;
plv_diff = merged_event_info.ripples_plv_diff;
corr_diff = merged_event_info.ripples_corr_diff;
ipsi_lag = [ipsi_lag_ripples{1} ipsi_lag_ripples{2}]';
contra_lag = [contra_lag_ripples{1} contra_lag_ripples{2}]';
group_id = merged_event_info.ripples_group_id;

sync_threshold = mean(abs([event_info(1).ripples_lag_threshold_low event_info(2).ripples_lag_threshold_low event_info(1).ripples_lag_threshold_high event_info(2).ripples_lag_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Magenta, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% ipsi contra diff
binnedArray = ipsi_HPC_MUA_baseline-contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;

%%%%%%%%%%%%%%%% V1 spikes
ipsi_V1_MUA = [PSTH_MUA(1).L_V1_ripples; PSTH_MUA(2).R_V1_ripples];
contra_V1_MUA = [PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).L_V1_ripples];
ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_ripples; PSTH_MUA_baseline(2).R_V1_ripples];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_ripples; PSTH_MUA_baseline(2).L_V1_ripples];

%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap1 = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap1 = temp;

%%% ipsi contra diff
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap1 = temp;



%%%%%%% Event grouping based on ipsi-contra difference features
event_idx = [];
% Ipsi-contra ripples HPC MUA by four different clusters
event_idx{1} = {find(group_id==4),find(group_id==3),find(group_id==2),find(group_id==1)};

% Ipsi-contra ripples HPC MUA by ipsi dominant cluster (different lag diff)
index = merged_event_info.ripples_index;
event_idx{2} = {intersect(find(group_id == 4 & lag_diff <prctile(lag_diff(group_id == 4& lag_diff <0),50) ),index),intersect(find(group_id == 4 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 4& lag_diff <0),50) ),index),...
    intersect(find(group_id == 2 | group_id == 3 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),index),...
    intersect(find(group_id == 4 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id ==4& lag_diff >0),50) ),index),intersect(find(group_id == 4 &lag_diff >prctile(lag_diff(group_id == 4& lag_diff >0),50) ),index)};

% Ipsi-contra ripples HPC MUA by contra dominant clusters (different lag diff)
event_idx{3} = {intersect(find(group_id == 1 & lag_diff <prctile(lag_diff(group_id == 1& lag_diff <0),50) ),index),intersect(find(group_id == 1 & lag_diff <0 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff <0),50) ),index),...
    intersect(find(group_id == 2 | group_id == 3 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),index),...
    intersect(find(group_id == 1 & lag_diff>0 &lag_diff <prctile(lag_diff(group_id ==1& lag_diff >0),50) ),index),intersect(find(group_id == 1 &lag_diff >prctile(lag_diff(group_id == 1& lag_diff >0),50) ),index)};

% Ipsi-contra ripples HPC MUA by bilateral clusters (different lag diff)
index = intersect(merged_event_info.ripples_index,find(group_id == 2 | group_id == 3));
lag_thresholds = prctile(lag_diff(index),[0 25 50 75 100]);

for n = 1:length(lag_thresholds)-1
    event_idx{4}{n} =intersect(find(lag_diff>lag_thresholds(n)&lag_diff <lag_thresholds(n+1)),index);
end

% Ipsi-contra ripples HPC MUA by lag diff
index = merged_event_info.ripples_index;
lag_thresholds = prctile(lag_diff(index),[0 20 40 60 80 100]);

for n = 1:length(lag_thresholds)-1
    event_idx{5}{n} =intersect(find(lag_diff>lag_thresholds(n)&lag_diff <lag_thresholds(n+1)),index);
end

group_name=[];
group_name{1} = {'Ipsi dominant','Bilateral high plv diff','Bilateral low plv diff','Contra dominant','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% ipsi dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% contra dominant clusters
group_name{4} = {'Top 0-25% ipsi leading','Top 25-50% ipsi leading','Top 50-75% ipsi leading','Top 75-100% ipsi leading','Shuffled'};% purely based on lags
group_name{5} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags

title_names = {'Ipsi-contra ripples MUA by four clusters', 'Ipsi-contra ripples MUA by ipsi dominant cluster (different lag diff)','Ipsi-contra ripples MUA by contra dominant clusters (different lag diff)',...
    'Ipsi-contra ripples MUA by bilateral clusters (different lag diff)','Ipsi-contra ripples MUA by lag diff'};

% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged.x = x;
MUA_PSTH_merged.ipsi_ripples_V1 = [];
MUA_PSTH_merged.contra_ripples_V1 = [];
MUA_PSTH_merged.ipsi_contra_diff_ripples_V1 = [];

MUA_PSTH_merged.ipsi_ripples_HPC = [];
MUA_PSTH_merged.contra_ripples_HPC = [];
MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC = [];

for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged.ipsi_ripples_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged.contra_ripples_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged.ipsi_contra_diff_ripples_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged.ipsi_ripples_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged.contra_ripples_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged.ipsi_ripples_V1_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged.contra_ripples_V1_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged.ipsi_contra_diff_ripples_V1_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged.ipsi_ripples_HPC_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged.contra_ripples_HPC_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged.ripples_groups = [title_names {'all ripples'}];
MUA_PSTH_merged.ripples_index = [event_idx (1:size(ipsi_V1_MUA,1))'];


for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1 | ngroup ==4
        colour_lines = [0,90,50;254,145,198;228,42,168;74,20,134]/256; % Dark Green , light Magenta, dark mageta, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_ripples_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_ripples_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.2 0.4])
    xlim([-0.2 0.2])
    
    % xline(0,'r')
    title('V1 Ipsi MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_ripples_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_ripples_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.2 0.4])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.2 0.4])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_ripples_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_ripples_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.5 5.6])
    xlim([-0.2 0.2])
    
    % xline(0,'r')
    title('HPC Ipsi MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_ripples_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_ripples_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.5 5.6])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.5 0.85])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])




%%%%%%%%% All ripples plots

clear ERROR_SHADE
fig = figure('Color','w');
% fig.Position = [350 59 1100 465];
fig.Position = [350 59 1100 930];
fig.Name ='ipsi-contra ripples MUA';

colour_lines = [0,90,50;74,20,134]/256; % dark purple, meganta, light green, dark green

nexttile
binnedArray = MUA_PSTH_merged.ipsi_ripples_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)


binnedArray = MUA_PSTH_merged.contra_ripples_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged.ipsi_ripples_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

xlim([-0.5 0.5])
ylim([-0.2 0.4])


% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('V1 MUA')
xlabel('Time relative to ripple peaktimes (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline

binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

xlim([-0.5 0.5])
ylim([-0.2 0.4])

% xline(0,'r')
legend([ERROR_SHADE(1:2)],{'ipsi-contra diff','shuffled'},'box','off')
title('V1 MUA difference')
xlabel('Time relative to ripple peaktimes (s)')
ylabel('ipsi-contra MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
binnedArray = MUA_PSTH_merged.ipsi_ripples_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

binnedArray = MUA_PSTH_merged.contra_ripples_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged.ipsi_ripples_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

xlim([-0.5 0.5])
ylim([-0.5 5.6])

% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('HPC MUA')
xlabel('Time relative to ripple peaktimes (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray =  MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline

binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

xlim([-0.5 0.5])
ylim([-0.2 0.4])

% xline(0,'r')
legend([ERROR_SHADE(1:2)],{'ipsi-contra diff','shuffled'},'box','off')
title('HPC MUA difference')
xlabel('Time relative to ripple peaktimes (s)')
ylabel('ipsi-contra MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%% Calculate Unilateral-bilateral diff in response
index = merged_event_info.ripples_index;
bilateral_index = intersect(find(group_id == 2 | group_id == 3 & lag_diff < sync_threshold & lag_diff > -sync_threshold ),index);
group_name=[];
group_name{1} = {'Ipsi dominant','Bilateral high plv diff','Bilateral low plv diff','Contra dominant','Shuffled'};
group_name{2} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% ipsi dominant clusters
group_name{3} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% contra dominant clusters
group_name{4} = {'Top 0-25% ipsi leading','Top 25-50% ipsi leading','Top 50-75% ipsi leading','Top 75-100% ipsi leading','Shuffled'};% purely based on lags
group_name{5} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags

title_names = {'Ipsi-contra ripples unilateral-bilateral MUA by four clusters', 'Ipsi-contra ripples unilateral-bilateral MUA by ipsi dominant cluster (different lag diff)','Ipsi-contra ripples unilateral-bilateral MUA by contra dominant clusters (different lag diff)',...
    'Ipsi-contra ripples unilateral-bilateral MUA by bilateral clusters (different lag diff)','Ipsi-contra ripples unilateral-bilateral MUA by lag diff'};

% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

for ngroup = 1:length(group_name)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1 | ngroup ==4
        colour_lines = [0,90,50;254,145,198;228,42,168;74,20,134]/256; % Dark Green , light Magenta, dark mageta, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end



    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(ipsi_V1_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_ripples_V1{ngroup}{i}-binnedArray1;

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    PLOT = plot([-1 1],[0 0],'k');hold on;

    ylim([-0.15 0.15])
    xlim([-0.2 0.2])
    
    % xline(0,'r')
    title('V1 Ipsi MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Unilateral-Bilateral MUA difference (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(contra_V1_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_ripples_V1{ngroup}{i}-binnedArray1;

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    PLOT = plot([-1 1],[0 0],'k');hold on;
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.15 0.15])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Unilateral-Bilateral MUA difference (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(ipsi_V1_MUA(bilateral_index,:)-contra_V1_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_V1{ngroup}{i}-binnedArray1;

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    PLOT = plot([-1 1],[0 0],'k');hold on;
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.15 0.15])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Unilateral-Bilateral MUA difference (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(ipsi_HPC_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_ripples_HPC{ngroup}{i}-binnedArray1;

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    PLOT = plot([-1 1],[0 0],'k');hold on;

    %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1.5 0.7])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('HPC Ipsi MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Unilateral-Bilateral MUA difference (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(contra_HPC_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_ripples_HPC{ngroup}{i}-binnedArray1;

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    PLOT = plot([-1 1],[0 0],'k');hold on;
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1.5 0.7])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Unilateral-Bilateral MUA difference (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(ipsi_HPC_MUA(bilateral_index,:)-contra_HPC_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_ripples_HPC{ngroup}{i}-binnedArray1;

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    PLOT = plot([-1 1],[0 0],'k');hold on;
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.85 0.85])
    xlim([-0.2 0.2])

    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Unilateral-Bilateral MUA difference (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])

%%
%%%%%
%%%%%
%%%%%
%%%%%

%% Plot spindles averaged MUA
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_event_info.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

ipsi_V1_MUA = [PSTH_MUA(1).L_V1_spindles; PSTH_MUA(2).R_V1_spindles];
contra_V1_MUA = [PSTH_MUA(1).R_V1_spindles; PSTH_MUA(2).L_V1_spindles];

ipsi_V1_MUA_baseline = [PSTH_MUA_baseline(1).L_V1_spindles; PSTH_MUA_baseline(2).R_V1_spindles];
contra_V1_MUA_baseline = [PSTH_MUA_baseline(1).R_V1_spindles; PSTH_MUA_baseline(2).L_V1_spindles];

% ipsi_probability = [PSTH_MUA(1).L_V1_DOWN; PSTH_MUA(2).R_V1_DOWN];
% contra_probability = [PSTH_MUA(2).L_V1_DOWN; PSTH_MUA(1).R_V1_DOWN];

lag_diff = merged_event_info.spindles_lag_diff;
plv_diff = merged_event_info.spindles_plv_diff;
corr_diff = merged_event_info.spindles_corr_diff;
ipsi_lag = [ipsi_lag_spindles{1} ipsi_lag_spindles{2}]';
contra_lag = [contra_lag_spindles{1} contra_lag_spindles{2}]';
group_id = merged_event_info.spindles_group_id;

sync_threshold = mean(abs([event_info(1).spindles_threshold_low event_info(2).spindles_threshold_low event_info(1).spindles_threshold_high event_info(2).spindles_threshold_high]));

colour_lines=[];
colour_lines = [0,90,50;74,20,134;228,42,168]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 


%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap = temp;

%%% Contra
binnedArray = contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap = temp;

%%% Contra
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;


%%%%%%%%%%%%%% HPC spikes
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_spindles; PSTH_MUA(2).R_HPC_spindles];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_spindles; PSTH_MUA(2).L_HPC_spindles];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_spindles; PSTH_MUA_baseline(2).R_HPC_spindles];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_spindles; PSTH_MUA_baseline(2).L_HPC_spindles];

%%%%% calculate shuffled MUA baseline
%%% Ipsi
binnedArray = ipsi_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_baseline_bootstrap1 = temp;

%%% Contra
binnedArray = contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

contra_baseline_bootstrap1 = temp;



%%%%%%%%%%%%%%%% grouping events based on ipsi-contra feature difference
event_idx = [];
% Ipsi-contra spindles V1 MUA by four clusters
index = merged_event_info.spindles_index;
event_idx{1} = {intersect(find(group_id == 4),index),intersect(find(group_id == 3),index),intersect(find(group_id == 2),index),intersect(find(group_id == 1),index)};

% Ipsi-contra spindles V1 MUA by lag diff
index = merged_event_info.spindles_index;

lag_thresholds = prctile(lag_diff(index),[0 20 40 60 80 100]);

for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =intersect(find(lag_diff>lag_thresholds(n)&lag_diff <lag_thresholds(n+1)),index);
end

group_name=[];
group_name{1} = {'Ipsi dominant high lag diff','Ipsi dominant low lag diff','Bilaterally synchronised','Contra dominant','Shuffled'};
group_name{2} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags

title_names = {'Ipsi-contra spindles V1 MUA by four clusters', 'Ipsi-contra spindles V1 MUA by lag diff'};
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged.x = x;
MUA_PSTH_merged.ipsi_spindles_V1 = [];
MUA_PSTH_merged.contra_spindles_V1 = [];
MUA_PSTH_merged.ipsi_contra_diff_spindles_V1 = [];

MUA_PSTH_merged.ipsi_spindles_HPC = [];
MUA_PSTH_merged.contra_spindles_HPC = [];
MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC = [];

for ngroup = 1:length(event_idx)+1

    group_index = [];
    if ngroup <= length(event_idx)
        group_index = event_idx{ngroup};

    else
        group_index{1} = (1:size(ipsi_V1_MUA,1))';
    end

    for i = 1:length(group_index)
        index =group_index{i};

        binnedArray1 = ipsi_V1_MUA(index,:);
        binnedArray2 = contra_V1_MUA(index,:);
        binnedArray3 = ipsi_V1_MUA(index,:)-contra_V1_MUA(index,:);
        binnedArray4 = ipsi_HPC_MUA(index,:);
        binnedArray5 = contra_HPC_MUA(index,:);
        binnedArray6 = ipsi_HPC_MUA(index,:)-contra_HPC_MUA(index,:);

        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp5=[];
        temp6=[];
        parfor iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(binnedArray1,1),size(binnedArray1,1));
            temp1(iBoot,:) =  mean(binnedArray1(event_id,:),'omitnan');
            temp2(iBoot,:) =  mean(binnedArray2(event_id,:),'omitnan');
            temp3(iBoot,:) =  mean(binnedArray3(event_id,:),'omitnan');
            temp4(iBoot,:) =  mean(binnedArray4(event_id,:),'omitnan');
            temp5(iBoot,:) =  mean(binnedArray5(event_id,:),'omitnan');
            temp6(iBoot,:) =  mean(binnedArray6(event_id,:),'omitnan');
            % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
        end

        MUA_PSTH_merged.ipsi_spindles_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged.contra_spindles_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged.ipsi_contra_diff_spindles_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged.ipsi_spindles_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged.contra_spindles_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged.ipsi_spindles_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged.contra_spindles_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged.ipsi_contra_diff_spindles_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged.ipsi_spindles_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged.contra_spindles_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged.spindles_groups = [title_names {'all spindles'}];
MUA_PSTH_merged.spindles_index = [event_idx (1:size(ipsi_V1_MUA,1))'];


for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;254,145,198;228,42,168;74,20,134]/256; % Dark Green , light Magenta, dark mageta, dark purple
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_spindles_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_spindles_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.6])
    title('V1 Ipsi MUA')
    xlabel('Time relative to spindle onset (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_spindles_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_spindles_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.6])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to spindle onset (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_V1{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.6])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to spindle onset (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_spindles_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_spindles_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    ylim([-0.3 0.6])
    title('HPC Ipsi MUA')
    xlabel('Time relative to spindle onset (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.contra_spindles_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.contra_spindles_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    ylim([-0.3 0.6])


    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to spindle onset (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC{ngroup}{i};

        % nprobe = 1;
        time_wondows = [-1 1];
        time_bin = 0.01;
        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


        y = mean(binnedArray,'omitnan');
        LCI = prctile(binnedArray,2.5);
        UCI = prctile(binnedArray,97.5);

        PLOT = plot(x,y,'Color',colour_lines(i,:));hold on;
        ERROR_SHADE(i) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(i,:),'FaceAlpha','0.3','LineStyle','none');
        xline(0,'r',LineWidth=1)
    end

    % baseline
    binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.6])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to spindle onset (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end



%%%%%%%%% All spindles plots

clear ERROR_SHADE
fig = figure('Color','w');
% fig.Position = [350 59 1100 465];
fig.Position = [350 59 1100 930];
fig.Name ='ipsi-contra spindles V1 MUA';


colour_lines = [0,90,50;74,20,134]/256; % dark purple, meganta, light green, dark green

nexttile
binnedArray = MUA_PSTH_merged.ipsi_spindles_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)


binnedArray = MUA_PSTH_merged.contra_spindles_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;


y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged.ipsi_spindles_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.3 0.6])


% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('V1 MUA')
xlabel('Time relative to spindle onset (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_V1{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline

binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_V1_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.3 0.6])

% xline(0,'r')
legend([ERROR_SHADE(1:2)],{'ipsi-contra diff','shuffled'},'box','off')
title('V1 MUA difference')
xlabel('Time relative to spindle onset (s)')
ylabel('ipsi-contra MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
binnedArray = MUA_PSTH_merged.ipsi_spindles_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

binnedArray = MUA_PSTH_merged.contra_spindles_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline
binnedArray = MUA_PSTH_merged.ipsi_spindles_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.3 0.6])

% xline(0,'r')
legend([ERROR_SHADE(1:end)],{'ipsi','contra','shuffled'},'box','off')
title('HPC MUA')
xlabel('Time relative to spindle onset (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binnedArray =  MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC{end}{1};
% nprobe = 1;
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

y = mean(binnedArray,'omitnan');
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)

% baseline

binnedArray = MUA_PSTH_merged.ipsi_contra_diff_spindles_HPC_baseline;
y = mean(binnedArray,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(binnedArray,2.5);
UCI = prctile(binnedArray,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');

ylim([-0.3 0.6])

% xline(0,'r')
legend([ERROR_SHADE(1:2)],{'ipsi-contra diff','shuffled'},'box','off')
title('HPC MUA difference')
xlabel('Time relative to spindle onset (s)')
ylabel('ipsi-contra MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)




%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])




save(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_overlap_control.mat'),'MUA_PSTH_merged_overlap_control');

