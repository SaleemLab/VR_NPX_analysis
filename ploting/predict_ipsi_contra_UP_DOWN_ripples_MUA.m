function predict_ipsi_contra_UP_DOWN_ripples_MUA(slow_waves_all,ripples_all,spindles_all)

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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_baseline.mat'));
probability_normalised_whole_baseline = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

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
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded.mat'),'MUA_PSTH_merged_thresholded');


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


%% Grouping DOWN/UP events and ripples based on lags
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_baseline.mat'));
probability_psth_whole_baseline = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability;

% Calculate ripples probability according to unilateral vs bilateral events
%%%%% Ripple info
ripples_times =  merged_event_info.ripples_ints;
ripples_group_info = merged_event_info.ripples_group_id;
ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id;
ripples_original_index = [find(ripples_all(1).SWS_index); find(ripples_all(2).SWS_index)];


ripples_lag_diff = merged_event_info.ripples_lag_diff;
ripples_all_overlap_idx = merged_event_info.ripples_overlap_idx_all{end};
ripples_non_overlap_idx = merged_event_info.ripples_non_overlap_idx{end};
ripples_all_lags =merged_event_info.ripples_lags_all{end};

% merged_event_info.ripples_index_sorted

%%%%%%%%%%%%%% Grab ripple info during UP
L_first_ripples_power = [];
L_first_ripples_lag = [];
L_first_ripples_lag_diff = [];

L_last_ripples_power = [];
L_last_ripples_lag = [];
L_last_ripples_lag_diff = [];

L_last_ripples_lag_diff = [];
R_last_ripples_lag_diff = [];

L_time_from_last_ripples=[];
R_time_from_last_ripples=[];

L_ripple_counts = [];
R_ripple_counts = [];

R_first_ripples_power = [];
R_first_ripples_lag = [];
R_first_ripples_lags_diff = [];

R_last_ripples_power = [];
R_last_ripples_lag = [];
R_last_ripples_lags_diff = [];

for nprobe = 1:2


    unique(event_info(nprobe).L_ripple_normalised_UP_duration(:,2));
    for nevent = 1:length(probability_psth_whole(nprobe).UP_all_index)
        ripples_index = find(probability_psth_whole(nprobe).UP_all_index(nevent)==event_info(nprobe).L_ripple_normalised_UP_duration(:,2));
        if isempty(ripples_index) ==0
            L_time_from_last_ripples{nprobe}(nevent) =slow_waves_all(nprobe).UP_ints(probability_psth_whole(nprobe).UP_all_index(nevent),2) - ripples_all(1).peaktimes(event_info(nprobe).L_ripple_normalised_UP_duration(ripples_index(end),1));

            L_ripple_counts{nprobe}(nevent) = length(ripples_index);

            L_first_ripples_power{nprobe}(nevent) = ripples_all(1).peak_zscore(event_info(nprobe).L_ripple_normalised_UP_duration(ripples_index(end),1));
            L_last_ripples_power{nprobe}(nevent) = ripples_all(1).peak_zscore(event_info(nprobe).L_ripple_normalised_UP_duration(ripples_index(end),1));


            ripples_index = event_info(nprobe).L_ripple_normalised_UP_duration(ripples_index,1);
            L_first_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(find(ripples_original_index == ripples_index(1) & ripples_hemisphere_id == 1));
            L_last_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(find(ripples_original_index == ripples_index(end) & ripples_hemisphere_id == 1));

            [C,ia,ib] = intersect(ripples_all_overlap_idx,find(ripples_original_index == ripples_index(1) & ripples_hemisphere_id == 1));
            if ~isempty(C)
                L_first_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
            end

            [C,ia,ib] = intersect(ripples_all_overlap_idx,find(ripples_original_index == ripples_index(end) & ripples_hemisphere_id == 1));
            if ~isempty(C)
                L_last_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
            end
        else
            L_time_from_last_ripples{nprobe}(nevent) = 0;
            L_ripple_counts{nprobe}(nevent) = 0;
            L_first_ripples_lag_diff{nprobe}(nevent) = nan;
            L_first_ripples_lag{nprobe}(nevent) = nan;
            L_first_ripples_power{nprobe}(nevent) = nan;

            L_last_ripples_lag_diff{nprobe}(nevent) = nan;
            L_last_ripples_lag{nprobe}(nevent) = nan;
            L_last_ripples_power{nprobe}(nevent) = nan;
        end
    end

    for nevent = 1:length(probability_psth_whole(nprobe).UP_all_index)
        ripples_index = find(probability_psth_whole(nprobe).UP_all_index(nevent)==event_info(nprobe).R_ripple_normalised_UP_duration(:,2));
        if isempty(ripples_index) ==0
            R_time_from_last_ripples{nprobe}(nevent) = slow_waves_all(nprobe).UP_ints(probability_psth_whole(nprobe).UP_all_index(nevent),2) - ripples_all(2).peaktimes(event_info(nprobe).R_ripple_normalised_UP_duration(ripples_index(1),1)) ;

            R_ripple_counts{nprobe}(nevent) = length(ripples_index);

            R_first_ripples_power{nprobe}(nevent) = ripples_all(2).peak_zscore(event_info(nprobe).R_ripple_normalised_UP_duration(ripples_index(1),1));
            R_last_ripples_power{nprobe}(nevent) = ripples_all(2).peak_zscore(event_info(nprobe).R_ripple_normalised_UP_duration(ripples_index(end),1));

            ripples_index = event_info(nprobe).R_ripple_normalised_UP_duration(ripples_index,1);
            R_first_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(find(ripples_original_index == ripples_index(1) & ripples_hemisphere_id == 2));
            R_last_ripples_lag_diff{nprobe}(nevent) = ripples_lag_diff(find(ripples_original_index == ripples_index(end) & ripples_hemisphere_id == 2));

            [C,ia,ib] = intersect(ripples_all_overlap_idx,find(ripples_original_index == ripples_index(1) & ripples_hemisphere_id == 2));
            if ~isempty(C)
                R_first_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
            end

            [C,ia,ib] = intersect(ripples_all_overlap_idx,find(ripples_original_index == ripples_index(end) & ripples_hemisphere_id == 2));
            if ~isempty(C)
                R_last_ripples_lag{nprobe}(nevent) = ripples_all_lags(ia);
            end
        else
            R_time_from_last_ripples{nprobe}(nevent)= 0;
            R_ripple_counts{nprobe}(nevent) = 0;
            R_first_ripples_lag_diff{nprobe}(nevent) = nan;
            R_first_ripples_lag{nprobe}(nevent) = nan;
            R_first_ripples_power{nprobe}(nevent) = nan;

            R_last_ripples_lag_diff{nprobe}(nevent) = nan;
            R_last_ripples_lag{nprobe}(nevent) = nan;
            R_last_ripples_power{nprobe}(nevent) = nan;
        end
    end
end

%%%%% Grab ripples info
ripples_info=[];
varnames = {
    'time_from_last_ripples'
    'first_ripples_power'
    'last_ripples_power'
    'first_ripples_lag'
    'last_ripples_lag'
    'first_ripples_lag_diff'
    'last_ripples_lag_diff'
    'ripple_counts'
    };

% Loop through each variable name
for i = 1:length(varnames)
    varname = varnames{i};

    % Build L and R variable names dynamically
    L_var = eval(['L_' varname]);
    R_var = eval(['R_' varname]);
    
    % Create ipsi and contra versions
    ripples_info.(['ipsi_' varname '_UP']) = [L_var{1}, R_var{2}];
    ripples_info.(['contra_' varname '_UP']) = [R_var{1}, L_var{2}];
end


%%%%% UP info
event_times = merged_event_info.UP_ints;
hemisphere_id = merged_event_info.UP_hemisphere_id;
% lag_diff = merged_event_info.UP_lag_diff;
% group_id = merged_event_info.UP_group_id;
all_overlap_idx = merged_event_info.UP_overlap_idx_all{end};
non_overlap_idx = merged_event_info.UP_non_overlap_idx{end};
all_lags =merged_event_info.UP_lags_all{end};

%%%%%%%%%% Based on detection lag
probability_merged = [];

index=[];

for nprobe = 1:2
    index{1}{nprobe} = intersect(ripples_all_overlap_idx(ripples_all_lags<-0.005& ripples_all_lags>-0.02),find(ripples_hemisphere_id == nprobe))';% leading
    index{2}{nprobe} = intersect(ripples_all_overlap_idx(ripples_all_lags>0.005 & ripples_all_lags<0.02),find(ripples_hemisphere_id == nprobe))';% lagging
    index{3}{nprobe} = ripples_all_overlap_idx(ripples_all_lags<=0&ripples_all_lags>-0.005); % bilaterally synchronised
    index{4}{nprobe} = [intersect(ripples_all_overlap_idx(ripples_all_lags<-0.02 | ripples_all_lags>0.02),find(ripples_hemisphere_id == nprobe))' intersect(ripples_non_overlap_idx,find(ripples_hemisphere_id == nprobe))']; % unilateral
end
group_name = {'ipsi_leading_ripple','contra_leading_ripple','bilateral_ripple','unilateral_ripple'};



%%%%%%%%% LFP based lag

index=[];
probability_merged = [];
for nprobe = 1:2
    lag_thresholds = (prctile(ripples_lag_diff,[12.5 37.5]));
    index{1}{nprobe} = intersect(find(ripples_lag_diff>lag_thresholds(1) & ripples_lag_diff<lag_thresholds(2)),find(ripples_hemisphere_id == nprobe))';% leading
    lag_thresholds = (prctile(ripples_lag_diff,[62.5 87.5]));
    index{2}{nprobe} = intersect(find(ripples_lag_diff>lag_thresholds(1) & ripples_lag_diff<lag_thresholds(2)),find(ripples_hemisphere_id == nprobe))';% leading
    lag_thresholds = (prctile(ripples_lag_diff,[37.5 62.5]));
    index{3}{nprobe} = find(ripples_lag_diff>lag_thresholds(1) & ripples_lag_diff<lag_thresholds(2)); % bilaterally synchronised
    lag_thresholds = (prctile(ripples_lag_diff,[12.5 87.5]));
    index{4}{nprobe} = intersect(find(ripples_lag_diff<lag_thresholds(1) | ripples_lag_diff>lag_thresholds(2)),find(ripples_hemisphere_id == nprobe))';% leading
end
group_name = {'ipsi_leading_ripple','contra_leading_ripple','bilateral_ripple','unilateral_ripple'};



% UP DOWN based on detection lags 
% Putatively use 0.05s as threshold for bilateral shared
% 0.2s as threshold for bilaterally coordinated but with unilaterally biased
% > 0.2s and/or non overlapping events 

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged_LFP_lag.mat'),'probability_merged');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_merged_LFP_lag.mat'),'probability_normalised_merged');
probability_merged_LFP_lag = probability_merged;
probability_normalised_merged_LFP_lag = probability_normalised_merged;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_merged.mat'),'probability_merged');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole_merged.mat'),'probability_normalised_merged');


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

ipsi_probability_baseline = [probability_psth_whole_baseline(1).L_ripples_UP; probability_psth_whole_baseline(2).R_ripples_UP];
contra_probability_baseline = [probability_psth_whole_baseline(1).R_ripples_UP; probability_psth_whole_baseline(2).L_ripples_UP];


%%%%%%%%%% Predict HPC MUA and ripples during DOWN-UP transition based on
%%%%%%%%%% V1 DOWN UP bilateral synchrony and magnitude
% index = all_overlap_idx(abs(lags)<=0.2);
% lag_index = (abs(lags)<=0.2);
% outputStruct = predict_HPC_MUA_DOWN_UP(ipsi_V1_MUA(index), contra_V1_MUA(index), ipsi_HPC_MUA(index), contra_HPC_MUA(index));

% 
% index = all_overlap_idx(abs(lags)<=0.15);
% lag_index = (abs(lags)<=0.15);
% output = predict_ripples_by_DOWN_UP_V1_MUA(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),ripples_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','V1excitation_HPCexcitation_output.mat'),'output');
% 


index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),ripples_info,'DOWN_UP_index',index,'DOWN_UP_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
save(fullfile(analysis_folder,'V1-HPC sleep interaction','V1synchrony_HPCexcitation_output.mat'),'output');


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



index = 1:size(ipsi_V1_MUA,1);
% lag_index = (abs(lags)<=2);
output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),ripples_info,'subject_id',subject_id(index));
save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1depression_output.mat'),'output');


index = all_overlap_idx(abs(lags)<=0.15);
lag_index = (abs(lags)<=0.15);
output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA(index,:), contra_V1_MUA(index,:), ipsi_HPC_MUA(index,:), contra_HPC_MUA(index,:),ipsi_probability(index,:), contra_probability(index,:),ripples_info,'UP_DOWN_index',index,'UP_DOWN_lag',abs(lags(lag_index)'),'subject_id',subject_id(index));
save(fullfile(analysis_folder,'V1-HPC sleep interaction','Ripples_V1synchrony_output.mat'),'output');


%%
%%%%%%%%%%%%%% Survival
[b, logl, H, stats] = coxphfit(abs(lags(lag_index)'), ripples_info.ipsi_time_from_last_ripples_UP(index)', 'Strata', subject_id);




colour_lines = [74,20,134;158,154,200]/256;
colour_lines = [0,90,50;116,196,118]/256;
all_overlap_idx
[C,index,ib] = intersect(all_overlap_idx(abs(lags)<=1),find(~isnan(ripples_info.ipsi_last_ripples_power_UP))');
[C,index,ib] = intersect(non_overlap_idx,find(~isnan(ripples_info.ipsi_last_ripples_power_UP))');

% % index = find(~isnan(ripples_info.ipsi_last_ripples_power_UP))';
data = ripples_info.ipsi_last_ripples_power_UP(index);
% data = UP_DOWN_lag(index);
timeToTransition_ripples = ripples_info.ipsi_time_from_last_ripples_UP(index);
ripple_threshold1 = prctile(data,25);
ripple_threshold2 = prctile(data,75);

[y,x,LCI,UCI] = ecdf(timeToTransition_ripples(data<ripple_threshold1),'Function', 'survivor','Alpha',0.05,'Bounds','on');
y(isnan(LCI)) = [];x(isnan(LCI)) = [];UCI(isnan(LCI)) = [];LCI(isnan(LCI)) = [];
y = y(x >= 0);LCI = LCI(x >= 0);UCI = UCI(x >= 0);x = x(x >= 0);
PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x' fliplr(x')],[UCI' fliplr(LCI')],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

  hold on;
  [y,x,LCI,UCI] = ecdf(timeToTransition_ripples(data>ripple_threshold2),'Function', 'survivor','Alpha',0.05,'Bounds','on');
  y(isnan(LCI)) = [];x(isnan(LCI)) = [];UCI(isnan(LCI)) = [];LCI(isnan(LCI)) = [];
  y = y(x >= 0);LCI = LCI(x >= 0);UCI = UCI(x >= 0);x = x(x >= 0);
  PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
  ERROR_SHADE(2) = patch([x' fliplr(x')],[UCI' fliplr(LCI')],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
  legend(ERROR_SHADE(1:2),{'top 25% ripple power','bottom 25% ripple power'},'box','off')

%   title(sprintf('%s V1 UP and %s ripples',hemisphere{nprobe},hemisphere{mprobe}))
xlabel('Time (s)');
ylabel('Survival Probability of UP');
xlim([0 0.2])
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


%%%%%%%%% Survival
binEdges = 0:0.01:0.3;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;


% [C,index,ib] = intersect(non_overlap_idx,find(~isnan(ripples_info.ipsi_last_ripples_power_UP))');
index = intersect(find(~isnan(ripples_info.ipsi_last_ripples_power_UP))',find(ripples_info.ipsi_ripple_counts_UP>0)');
data = ripples_info.ipsi_last_ripples_power_UP(index);

% index = intersect(find(~isnan(ripples_info.ipsi_last_ripples_lag_UP))',find(ripples_info.ipsi_ripple_counts_UP>0)');
% data = ripples_info.ipsi_last_ripples_lag_UP(index);
ripple_threshold1 = prctile(data,25);
ripple_threshold2 = prctile(data,75);
event_id = find(data<ripple_threshold1);
timeToTransition_ripples = ripples_info.ipsi_time_from_last_ripples_UP(index);


% colour_lines = [158,154,200;74,20,134]/256;
colour_lines = [116,196,118;0,90,50]/256;
subplot(2,2,1)
y_boot=[];
for iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling

    [temp_y,temp_x] = ecdf(timeToTransition_ripples(datasample(s,event_id,length(event_id))),'Function', 'survivor');
    temp_x(isnan(temp_y)) = [];temp_y(isnan(temp_y)) = [];
    [temp_x,idx]= unique(temp_x);
    temp_y = temp_y(idx);
    y_boot(iBoot,:) = interp1(temp_x, temp_y, binCenters, 'previous', 'extrap')';  % 'previous' for step behavior like ECDF
end

% yyaxis left
x = binCenters';  % 'previous' for step behavior like ECDF
y = mean(y_boot)';
LCI = prctile(y_boot,2.5)';
UCI = prctile(y_boot,97.5)';
%         y = y(x >= 0);LCI = LCI(x >= 0);UCI = UCI(x >= 0);x = x(x >= 0);
PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x' fliplr(x')],[UCI' fliplr(LCI')],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


subplot(2,2,2)
hold on;
x = binCenters';  % 'previous' for step behavior like ECDF
y = mean(diff(y_boot')')';
LCI = prctile(diff(y_boot')',2.5)';
UCI = prctile(diff(y_boot')',97.5)';
y = [0; y];LCI = [0; LCI];UCI = [0; UCI];
PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x' fliplr(x')],[UCI' fliplr(LCI')],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');



event_id = find(data>ripple_threshold2);
subplot(2,2,1)
y_boot=[];
for iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling

    [temp_y,temp_x] = ecdf(timeToTransition_ripples(datasample(s,event_id,length(event_id))),'Function', 'survivor');
    temp_x(isnan(temp_y)) = [];temp_y(isnan(temp_y)) = [];
    [temp_x,idx]= unique(temp_x);
    temp_y = temp_y(idx);
    y_boot(iBoot,:) = interp1(temp_x, temp_y, binCenters, 'previous', 'extrap')';  % 'previous' for step behavior like ECDF
end
hold on;
% yyaxis left
x = binCenters';  % 'previous' for step behavior like ECDF
y = mean(y_boot)';
LCI = prctile(y_boot,2.5)';
UCI = prctile(y_boot,97.5)';
%         y = y(x >= 0);LCI = LCI(x >= 0);UCI = UCI(x >= 0);x = x(x >= 0);
PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(1) = patch([x' fliplr(x')],[UCI' fliplr(LCI')],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


subplot(2,2,2)
hold on;
x = binCenters';  % 'previous' for step behavior like ECDF
y = mean(diff(y_boot')')';
LCI = prctile(diff(y_boot')',2.5)';
UCI = prctile(diff(y_boot')',97.5)';
y = [0; y];LCI = [0; LCI];UCI = [0; UCI];
PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(1) = patch([x' fliplr(x')],[UCI' fliplr(LCI')],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


end