function plot_ipsi_contra_MUA_synchronisation_detection_thresholded(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,event_info,sessions_to_process)

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
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'));
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;


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


%% Grabbing ipsi and contra values
load(fullfile(analysis_folder,'V1-HPC sleep interaction','k_cluster_ipsi_contra_events.mat'),'k_cluster');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');


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

%%


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
all_lags = [];
all_overlap_idx=[];
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



ints = merged_event_info.UP_ints;
%%%%%%%%%%%%%%%%
fig = figure('Color','w');
fig.Position = [30 60 1200 460];
fig.Name = 'DOWN-UP transition lag distribution based on detection';
hemisphere_texts = {'Left DOWN-UP transition','Right DOWN-UP transition'}

% for ngroup = 1:4
nexttile
lags = all_lags{end};
% lags = L_lags{ngroup}
histogram(lags,-2:0.01:2,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('probability');
title('DOWN-UP transition lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
% lags = R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1);  % R onset - L onset
histogram(lags,-2:0.01:2,'EdgeColor','none','Normalization','cdf'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('cumulative probability');
title('DOWN-UP transition lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



customColors = [74,20,134;228,42,168;0,90,50]/256; % dark purple, magenta and dark green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'DOWN-UP transition lag vs corr vs plv distribution with overlapping or non-overlapping events';
hemisphere_texts = {'Left DOWN-UP transition','Right DOWN-UP transition'}



colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:4
    nexttile
    X = [merged_event_info.UP_lag_diff(overlap_idx{ngroup}) merged_event_info.UP_corr_diff(overlap_idx{ngroup}) merged_event_info.UP_plv_diff(overlap_idx{ngroup})];
    idx = merged_event_info.UP_group_id(overlap_idx{ngroup});

    for k = 1:3
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            5, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.24 0.24])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('overlapping DOWN-UP transition (%.2f sec win)',windows_threshold(ngroup)))
    else
        title('overlapping DOWN-UP transition (full win)')
    end


    nexttile
    X = [merged_event_info.UP_lag_diff(non_overlap_idx{ngroup}) merged_event_info.UP_corr_diff(non_overlap_idx{ngroup}) merged_event_info.UP_plv_diff(non_overlap_idx{ngroup})];
    idx = merged_event_info.UP_group_id(non_overlap_idx{ngroup});
    for k = 1:3
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            5, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.24 0.24])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('non-overlapping DOWN-UP transition (%.2f sec win)',windows_threshold(ngroup)))
    else
        title('non-overlapping DOWN-UP transition (full win)')
    end
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


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

%%% Ipsi-Contra baseline
binnedArray = ipsi_V1_MUA_baseline-contra_V1_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap = temp;

%%%%% HPC MUA
ipsi_HPC_MUA = [PSTH_MUA(1).L_HPC_UP; PSTH_MUA(2).R_HPC_UP];
contra_HPC_MUA = [PSTH_MUA(1).R_HPC_UP; PSTH_MUA(2).L_HPC_UP];
ipsi_HPC_MUA_baseline = [PSTH_MUA_baseline(1).L_HPC_UP; PSTH_MUA_baseline(2).R_HPC_UP];
contra_HPC_MUA_baseline = [PSTH_MUA_baseline(1).R_HPC_UP; PSTH_MUA_baseline(2).L_HPC_UP];

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



%%% Ipsi-Contra baseline
binnedArray = ipsi_HPC_MUA_baseline-contra_HPC_MUA_baseline;
temp=[];
parfor iBoot = 1:1000
    s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
    event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
    temp(iBoot,:) =  mean(binnedArray(event_id,:),'omitnan');
    % temp(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
end

ipsi_contra_diff_baseline_bootstrap1 = temp;




%%%%%%%%%%%%% Grouping events
event_idx = [];
% Ipsi-contra UP_DOWN V1 MUA by three clusters
% lags = [R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1)];  % R onset - L onset
lags =all_lags{end};

event_idx{1} = {all_overlap_idx{end}(lags<-0.05),all_overlap_idx{end}(lags<0.05&lags>-0.05),all_overlap_idx{end}(lags>0.05)};

% event_idx{2} = {overlap_idx{end}(lags<prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<-0.05&lags>prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<0.05&lags>-0.05),overlap_idx{end}(lags>0.05&lags<prctile(lags(lags>0.05),50)),overlap_idx{end}(lags>prctile(lags(lags>0.05),50))};

lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{2};
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{3}{n} =(all_overlap_idx{2}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lag_thresholds = [-0.2 -0.1 -0.02 0.02 0.1 0.2];
for n = 1:length(lag_thresholds)-1
    event_idx{4}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
event_idx{5} = {all_overlap_idx{end}(lags>-0.2&lags<-0.05),all_overlap_idx{end}(lags<0.02&lags>-0.05),all_overlap_idx{end}(lags>0.05&lags<0.2),non_overlap_idx{end}};


lags =all_lags{end}(all_lags{end}>=-0.15 & all_lags{end}<=0.15);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lags =all_lags{end};
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{6}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
% group_name{4} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{4} = {'-0.2 to -0.1s','-0.1 to -0.02s','-0.02 to 0.02s','0.02 to 0.1s','0.1 to 0.2s','Shuffled'};% purely based on lags
group_name{5} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Non-overlapping','Shuffled'};
group_name{6} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
% group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
% group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra DOWN_UP MUA by three lags (full windows)', 'Ipsi-contra DOWN_UP MUA by bilateral thresholded lags (full windows)',...
    'Ipsi-contra DOWN_UP MUA by 5 lags (100ms windows)','Ipsi-contra DOWN_UP MUA by 5 lags (200ms windows)','Ipsi-contra DOWN_UP MUA with non-overlapping','Ipsi-contra DOWN_UP MUA by 5 lags (150ms windows)'}
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged_thresholded.x = x;
MUA_PSTH_merged_thresholded.ipsi_DU_V1 = [];
MUA_PSTH_merged_thresholded.contra_DU_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_DU_HPC = [];
MUA_PSTH_merged_thresholded.contra_DU_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_HPC = [];

for ngroup = 1:length(event_idx)

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

        MUA_PSTH_merged_thresholded.ipsi_DU_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_DU_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_DU_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_DU_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.ipsi_DU_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged_thresholded.contra_DU_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_DU_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.contra_DU_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.DU_groups = [title_names];
MUA_PSTH_merged_thresholded.DU_index = [event_idx];



colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple

% nexttile
% scatter(ipsi_lag(all_overlap_idx{end}),all_lags{end},'MarkerFaceColor',colour_lines(1, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% nexttile
% scatter(contra_lag(all_overlap_idx{end}),all_lags{end},'MarkerFaceColor',colour_lines(3, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% nexttile
% scatter(lag_diff(all_overlap_idx{end}),all_lags{end},'MarkerFaceColor',colour_lines(2, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% 
% 
% nexttile
% scatter(corr_diff(all_overlap_idx{end}),all_lags{end},'MarkerFaceColor',colour_lines(2, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% nexttile
% scatter(plv_diff(all_overlap_idx{end}),all_lags{end},'MarkerFaceColor',colour_lines(2, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% nexttile
% scatter(lag_diff(all_overlap_idx{end}),all_lags{end},'MarkerFaceColor',colour_lines(2, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% 

%%%%%%%%%%%%% Distribution of detection thresholded events groups

customColors = [74,20,134;228,42,168;0,90,50]/256; % dark purple, magenta and dark green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'DOWN-UP transition lag vs corr vs plv distribution with overlapping or non-overlapping event clusters';
hemisphere_texts = {'Left DOWN-UP transition','Right DOWN-UP transition'}



% colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:length(event_idx)
    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear PLOT
    for k = 1:length(event_idx{ngroup})
        X = [lag_diff(event_idx{ngroup}{k}) corr_diff(event_idx{ngroup}{k}) plv_diff(event_idx{ngroup}{k})];
        hold on
        PLOT(k) = scatter3( ...
            X(:, 1), X(:, 2), X(:, 3), ...
            5, 'MarkerFaceColor',colour_lines(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.24 0.24])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title(title_names{ngroup})
    legend(PLOT(1:end),{group_name{ngroup}{1:end-1}},'box','off')
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


%%%%%%%%% Plot DU transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_DU_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_DU_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-1 1])
    title('V1 Ipsi MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_DU_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_DU_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 1])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-1 1])


    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_DU_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_DU_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.2])
    title('HPC Ipsi MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_DU_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_DU_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
%     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.2])

    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_DU_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.3 0.2])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to DOWN-UP transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')

%% Clustering of UP-DOWN transition lag vs corr vs plv distribution based on overlapping events
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
all_lags= [];
all_overlap_idx = [];
% for nsession = 1:max(slow_waves_all.UP_session_count)

%     L_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 1 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
%     R_ints = merged_event_info.UP_ints((merged_event_info.UP_ints(:,1) - nsession * 1000000) > 0 &merged_event_info.UP_hemisphere_id == 2 & (merged_event_info.UP_ints(:,1) - nsession * 1000000) < 1000000,:);
L_idx = find(merged_event_info.DOWN_hemisphere_id == 1);
R_idx = find(merged_event_info.DOWN_hemisphere_id == 2);
L_ints = merged_event_info.DOWN_ints(L_idx,:);
R_ints = merged_event_info.DOWN_ints(R_idx,:);

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



ints = merged_event_info.DOWN_ints;


%%%%%%%%%%%%%%%%
fig = figure('Color','w');
fig.Position = [30 60 1200 460];
fig.Name = 'UP-DOWN transition lag distribution based on detection';
hemisphere_texts = {'Left UP-DOWN transition','Right UP-DOWN transition'}

% for ngroup = 1:4
nexttile
lags = all_lags{end};
histogram(lags,-0.5:0.01:0.5,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('probability');
title('UP-DOWN transition lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
% lags = R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1);  % R onset - L onset
histogram(lags,-0.5:0.01:0.5,'EdgeColor','none','Normalization','cdf'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('cumulative probability');
title('UP-DOWN transition lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



customColors = [74,20,134;228,42,168;0,90,50]/256; % dark purple, magenta and dark green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'UP-DOWN transition lag vs corr vs plv distribution with overlapping or non-overlapping events';
% hemisphere_texts = {'Left DOWN-UP transition','Right DOWN-UP transition'}

colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:4
    nexttile
    X = [merged_event_info.DOWN_lag_diff(overlap_idx{ngroup}) merged_event_info.DOWN_corr_diff(overlap_idx{ngroup}) merged_event_info.DOWN_plv_diff(overlap_idx{ngroup})];
    idx = merged_event_info.DOWN_group_id(overlap_idx{ngroup});

    for k = 1:3
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            5, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.24 0.24])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('overlapping UP-DOWN transition (%.2f sec windows)',windows_threshold(ngroup)))
    else
        title('overlapping UP-DOWN transition (full windows)')
    end


    nexttile
    X = [merged_event_info.DOWN_lag_diff(non_overlap_idx{ngroup}) merged_event_info.DOWN_corr_diff(non_overlap_idx{ngroup}) merged_event_info.DOWN_plv_diff(non_overlap_idx{ngroup})];
    idx = merged_event_info.DOWN_group_id(non_overlap_idx{ngroup});

    for k = 1:3
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            5, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.24 0.24])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('non-overlapping  UP-DOWN transition (%.2f sec windows)',windows_threshold(ngroup)))
    else
        title('non-overlapping UP-DOWN transition (full windows)')
    end
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')




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
% group_id = merged_event_info.UP_group_id;

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

%%%%%%%%%%%%% Grouping events
event_idx = [];
% Ipsi-contra UP_DOWN V1 MUA by three clusters
% lags = [R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1)];  % R onset - L onset
lags =all_lags{end};

event_idx{1} = {all_overlap_idx{end}(lags<-0.05),all_overlap_idx{end}(lags<0.05&lags>-0.05),all_overlap_idx{end}(lags>0.05)};

% event_idx{2} = {overlap_idx{end}(lags<prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<-0.05&lags>prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<0.05&lags>-0.05),overlap_idx{end}(lags>0.05&lags<prctile(lags(lags>0.05),50)),overlap_idx{end}(lags>prctile(lags(lags>0.05),50))};

lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{2};
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{3}{n} =(all_overlap_idx{2}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lag_thresholds = [-0.2 -0.1 -0.02 0.02 0.1 0.2];
for n = 1:length(lag_thresholds)-1
    event_idx{4}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
event_idx{5} = {all_overlap_idx{end}(lags>-0.2&lags<-0.05),all_overlap_idx{end}(lags<0.02&lags>-0.05),all_overlap_idx{end}(lags>0.05&lags<0.2),non_overlap_idx{end}};


lags =all_lags{end}(all_lags{end}<=0.2 & all_lags{end}>=-0.2);
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);

lags =all_lags{end};
for n = 1:length(lag_thresholds)-1
    event_idx{6}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
% group_name{4} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{4} = {'-0.2 to -0.1s','-0.1 to -0.02s','-0.02 to 0.02s','0.02 to 0.1s','0.1 to 0.2s','Shuffled'};% purely based on lags
group_name{5} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Non-overlapping','Shuffled'};
group_name{6} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
% group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
% group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra UP_DOWN MUA by three lags (full windows)', 'Ipsi-contra UP_DOWN MUA by bilateral thresholded lags (full windows)',...
    'Ipsi-contra UP_DOWN MUA by 5 lags (100ms windows)','Ipsi-contra UP_DOWN MUA by 5 lags (200ms windows)','Ipsi-contra UP_DOWN MUA with non-overlapping','Ipsi-contra UP_DOWN MUA by 5 lags (150ms windows)'}

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged_thresholded.x = x;
MUA_PSTH_merged_thresholded.ipsi_UD_V1 = [];
MUA_PSTH_merged_thresholded.contra_UD_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_UD_HPC = [];
MUA_PSTH_merged_thresholded.contra_UD_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_HPC = [];

for ngroup = 1:length(event_idx)

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

        MUA_PSTH_merged_thresholded.ipsi_UD_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_UD_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_UD_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_UD_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.ipsi_UD_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged_thresholded.contra_UD_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_UD_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.contra_UD_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.UD_groups = [title_names];
MUA_PSTH_merged_thresholded.UD_index = [event_idx];

%%%%%%%%%%%%% Distribution of detection thresholded events groups

customColors = [74,20,134;228,42,168;0,90,50]/256; % dark purple, magenta and dark green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'UP-DOWN transition lag vs corr vs plv distribution with overlapping or non-overlapping event clusters';
hemisphere_texts = {'Left UP-DOWN transition','Right UP-DOWN transition'}



% colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:length(event_idx)
    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear PLOT
    for k = 1:length(event_idx{ngroup})
        X = [lag_diff(event_idx{ngroup}{k}) corr_diff(event_idx{ngroup}{k}) plv_diff(event_idx{ngroup}{k})];
        hold on
        PLOT(k) = scatter3( ...
            X(:, 1), X(:, 2), X(:, 3), ...
            5, 'MarkerFaceColor',colour_lines(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.24 0.24])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title(title_names{ngroup})
    legend(PLOT(1:end),{group_name{ngroup}{1:end-1}},'box','off')
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


%%%%%%%%% Plot UD transition
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930]; 
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_UD_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_UD_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-1 0.7])
    title('V1 Ipsi MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_UD_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_UD_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);
    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-1 0.7])

    % xline(0,'r')
    title('V1 Contra MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_V1_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-1 0.7])

    % xline(0,'r')
    title('V1 Ipsi-contra MUA diff')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%%%%% HPC MUA
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_UD_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_UD_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})

    % xline(0,'r')
    ylim([-0.3 0.7])
    title('HPC Ipsi MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.contra_UD_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_UD_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    %     legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}})
    ylim([-0.3 0.7])
    % xline(0,'r')
    title('HPC Contra MUA')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_UD_HPC_baseline;
    y = mean(binnedArray,'omitnan');
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(binnedArray,2.5);
    UCI = prctile(binnedArray,97.5);


    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(length(ERROR_SHADE)+1) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend([ERROR_SHADE(1:end)],{group_name{ngroup}{1:end}},'box','off')
    ylim([-0.2 0.7])


    % xline(0,'r')
    title('HPC Ipsi-contra MUA diff')
    xlabel('Time relative to UP-DOWN transition (s)')
    ylabel('MUA activity (z)')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end


save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])



%% Clustering of ripples lag vs corr vs plv distribution based on overlapping events
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
L_idx = find(merged_event_info.ripples_hemisphere_id == 1);
R_idx = find(merged_event_info.ripples_hemisphere_id == 2);
L_ints = merged_event_info.ripples_ints(L_idx,:);
R_ints = merged_event_info.ripples_ints(R_idx,:);

windows_threshold = [0.01,0.05,0.1];

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


% ints = merged_event_info.UP_ints;
%%%%%%%%%%%%%%%%
fig = figure('Color','w');
fig.Position = [30 60 1200 460];
fig.Name = 'ripples lag distribution based on detection';
hemisphere_texts = {'Left ripples','Right ripples'}

% for ngroup = 1:4
nexttile
lags = all_lags{end}; 
histogram(lags,-0.1:0.002:0.1,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('probability');
title('ripples lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
% lags = R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1);  % R onset - L onset
histogram(lags,-0.1:0.002:0.1,'EdgeColor','none','Normalization','cdf'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('cumulative probability');
title('ripples lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



customColors = [74,20,134;228,42,168;254,145,198;0,90,50]/256; % purple, dark magenta, light magenta, green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'Ripples transition lag vs corr vs plv distribution with overlapping or non-overlapping events';
hemisphere_texts = {'Left Ripples','Right Ripples'}



colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:4
    nexttile
    X = [merged_event_info.ripples_lag_diff(overlap_idx{ngroup}) merged_event_info.ripples_corr_diff(overlap_idx{ngroup}) merged_event_info.ripples_plv_diff(overlap_idx{ngroup})];
    idx = merged_event_info.ripples_group_id(overlap_idx{ngroup});

    for k = 1:4
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            5, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.05 0.05])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('overlapping ripples (%.2f sec win)',windows_threshold(ngroup)))
    else
        title('overlapping ripples (full win)')
    end


    nexttile
    X = [merged_event_info.ripples_lag_diff(non_overlap_idx{ngroup}) merged_event_info.ripples_corr_diff(non_overlap_idx{ngroup}) merged_event_info.ripples_plv_diff(non_overlap_idx{ngroup})];
    idx = merged_event_info.ripples_group_id(non_overlap_idx{ngroup});
    for k = 1:4
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            5, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.05 0.05])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('non-overlapping Ripples transition (%.2f sec win)',windows_threshold(ngroup)))
    else
        title('non-overlapping ripples (full win)')
    end
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')

% 
% ipsi_lag_ripples
% nexttile
% scatter(all_lags{end},merged_event_info.ripples_lag_diff(all_overlap_idx{end}),5, 'MarkerFaceColor',customColors(1, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% nexttile
% scatter(all_lags{end},contra_lag(all_overlap_idx{end}),5, 'MarkerFaceColor',customColors(2, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% nexttile
% scatter(all_lags{end},contra_corr_ripples(all_overlap_idx{end}),5, 'MarkerFaceColor',customColors(3, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% 
% scatter(all_lags{end},ipsi_lag(all_overlap_idx{end}),5, 'MarkerFaceColor',customColors(3, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)

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

% merged_event_info.ripples_lag_diff(ll_overlap_idx{end})


%%%%%%% Event grouping based on ipsi-contra difference features
event_idx = [];
% Ipsi-contra ripples HPC MUA by four different clusters

lags =all_lags{end};

event_idx{1} = {all_overlap_idx{end}(lags<-0.005),all_overlap_idx{end}(lags<0.005&lags>-0.005),all_overlap_idx{end}(lags>0.005)};

% event_idx{2} = {overlap_idx{end}(lags<prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<-0.05&lags>prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<0.05&lags>-0.05),overlap_idx{end}(lags>0.05&lags<prctile(lags(lags>0.05),50)),overlap_idx{end}(lags>prctile(lags(lags>0.05),50))};

lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{2};
lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{3}{n} =(all_overlap_idx{2}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lag_thresholds = [-0.05 -0.02 -0.005 0.005 0.02 0.05];
for n = 1:length(lag_thresholds)-1
    event_idx{4}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
event_idx{5} = {all_overlap_idx{end}(lags>-0.2&lags<-0.005),all_overlap_idx{end}(lags<0.005&lags>-0.005),all_overlap_idx{end}(lags>0.005&lags<0.2),non_overlap_idx{end}};


group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
% group_name{4} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{4} = {'-0.05 to -0.02s','-0.02 to -0.005s','-0.005 to 0.005s','0.005 to 0.02s','0.02 to 0.05s','Shuffled'};% purely based on lags
group_name{5} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Non-overlapping','Shuffled'};

% group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
% group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra ripples MUA by three lags (full windows)', 'Ipsi-contra ripples MUA by bilateral thresholded lags (full windows)',...
    'Ipsi-contra ripples MUA by 5 lags percentile (100ms windows)','Ipsi-contra ripples MUA by 5 lags (50ms windows)','Ipsi-contra ripples MUA with non-overlapping'}

% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged_thresholded.x = x;
MUA_PSTH_merged_thresholded.ipsi_ripples_V1 = [];
MUA_PSTH_merged_thresholded.contra_ripples_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_ripples_HPC = [];
MUA_PSTH_merged_thresholded.contra_ripples_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_HPC = [];

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

        MUA_PSTH_merged_thresholded.ipsi_ripples_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_ripples_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_ripples_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_ripples_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.ipsi_ripples_V1_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.contra_ripples_V1_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_V1_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_ripples_HPC_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged_thresholded.contra_ripples_HPC_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_HPC_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ripples_groups = [title_names];
MUA_PSTH_merged_thresholded.ripples_index = [event_idx];


%%%%%%%%%%%%% Distribution of detection thresholded events groups

customColors = [74,20,134;228,42,168;0,90,50]/256; % dark purple, magenta and dark green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'ripples lag vs corr vs plv distribution with overlapping or non-overlapping event clusters';

% colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:length(event_idx)
    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear PLOT
    for k = 1:length(event_idx{ngroup})
        X = [lag_diff(event_idx{ngroup}{k}) corr_diff(event_idx{ngroup}{k}) plv_diff(event_idx{ngroup}{k})];
        hold on
        PLOT(k) = scatter3( ...
            X(:, 1), X(:, 2), X(:, 3), ...
            5, 'MarkerFaceColor',colour_lines(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.05 0.05])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title(title_names{ngroup})
    legend(PLOT(1:end),{group_name{ngroup}{1:end-1}},'box','off')
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


%%%%%%%%% Plot ripples MUA
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end


    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_ripples_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_ripples_V1_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.contra_ripples_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_ripples_V1_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_V1_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_ripples_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_ripples_HPC_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.contra_ripples_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_ripples_HPC_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_HPC_baseline;
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


%% Calculate Unilateral-bilateral diff in response
index = merged_event_info.ripples_index;
bilateral_index = all_overlap_idx{end}(lags<0.005&lags>-0.005);

group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{3} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
% group_name{4} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{4} = {'-0.05 to -0.02s','-0.02 to -0.005s','-0.005 to 0.005s','0.005 to 0.02s','0.02 to 0.05s','Shuffled'};% purely based on lags
group_name{5} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Non-overlapping','Shuffled'};

% group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
% group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra ripples unilateral-bilateral MUA by three lags (full windows)', 'Ipsi-contra ripples unilateral-bilateral MUA by bilateral thresholded lags (full windows)',...
    'Ipsi-contra ripples unilateral-bilateral MUA by 5 lags percentile (100ms windows)','Ipsi-contra ripples unilateral-bilateral MUA by 5 lags (50ms windows)','Ipsi-contra ripples unilateral-bilateral MUA with non-overlapping'}

% colour_lines = [0,90,50;74,20,134]/256; % Green Purple

% colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral
% 
% colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

for ngroup = 1:length(group_name)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end


    nexttile
    clear ERROR_SHADE
    binnedArray1 = mean(ipsi_V1_MUA(bilateral_index,:),'omitnan');
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_ripples_V1{ngroup}{i}-binnedArray1;

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

        binnedArray = MUA_PSTH_merged_thresholded.contra_ripples_V1{ngroup}{i}-binnedArray1;

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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_V1{ngroup}{i}-binnedArray1;

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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_ripples_HPC{ngroup}{i}-binnedArray1;

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
    ylim([-2 0.7])
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

        binnedArray = MUA_PSTH_merged_thresholded.contra_ripples_HPC{ngroup}{i}-binnedArray1;

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
    ylim([-2 0.7])
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_ripples_HPC{ngroup}{i}-binnedArray1;

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

%% Clustering of spindles lag vs corr vs plv distribution based on overlapping events
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
L_idx = find(merged_event_info.spindles_hemisphere_id == 1);
R_idx = find(merged_event_info.spindles_hemisphere_id == 2);
L_ints = merged_event_info.spindles_ints(L_idx,:);
R_ints = merged_event_info.spindles_ints(R_idx,:);

windows_threshold = [0.01,0.05,0.1];

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


% ints = merged_event_info.UP_ints;
%%%%%%%%%%%%%%%%
fig = figure('Color','w');
fig.Position = [30 60 1200 460];
fig.Name = 'spindles lag distribution based on detection';
hemisphere_texts = {'Left spindles','Right spindles'}

% for ngroup = 1:4
nexttile
lags = all_lags{end}; 
histogram(lags,-1:0.01:1,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('probability');
title('spindles lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
% lags = R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1);  % R onset - L onset
histogram(lags,-1:0.01:1,'EdgeColor','none','Normalization','cdf'); % You can change 50 to control bin numbers
xlabel('Lag (seconds)');  % or whatever unit your times are
ylabel('cumulative probability');
title('spindles lag distribution based on detection');
grid on;
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')


customColors = [74,20,134;228,42,168;254,145,198;0,90,50]/256; % purple, dark magenta, light magenta, green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'Spindles transition lag vs corr vs plv distribution with overlapping or non-overlapping events';
hemisphere_texts = {'Left spindles','Right spindles'}



colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:4
    nexttile
    X = [merged_event_info.spindles_lag_diff(overlap_idx{ngroup}) merged_event_info.spindles_corr_diff(overlap_idx{ngroup}) merged_event_info.spindles_plv_diff(overlap_idx{ngroup})];
    idx = merged_event_info.spindles_group_id(overlap_idx{ngroup});

    for k = 1:4
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            20, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.8 0.8])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('overlapping spindles (%.2f sec win)',windows_threshold(ngroup)))
    else
        title('overlapping spindles (full win)')
    end


    nexttile
    X = [merged_event_info.spindles_lag_diff(non_overlap_idx{ngroup}) merged_event_info.spindles_corr_diff(non_overlap_idx{ngroup}) merged_event_info.spindles_plv_diff(non_overlap_idx{ngroup})];
    idx = merged_event_info.spindles_group_id(non_overlap_idx{ngroup});
    for k = 1:4
        hold on
        scatter3( ...
            X(idx == k, 1), X(idx == k, 2), X(idx == k, 3), ...
            20, 'MarkerFaceColor',customColors(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.8 0.8])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    if ngroup <4
        title(sprintf('non-overlapping spindles transition (%.2f sec win)',windows_threshold(ngroup)))
    else
        title('non-overlapping spindles (full win)')
    end
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')

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



%%%%%%%%%%%%% Grouping events
event_idx = [];
% Ipsi-contra UP_DOWN V1 MUA by three clusters
% lags = [R_ints(R_overlap_idx{end},1) - L_ints(L_overlap_idx{end},1)];  % R onset - L onset
lags =all_lags{end};

event_idx{1} = {all_overlap_idx{end}(lags<-0.1),all_overlap_idx{end}(lags<0.1&lags>-0.1),all_overlap_idx{end}(lags>0.1)};

% event_idx{2} = {overlap_idx{end}(lags<prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<-0.05&lags>prctile(lags(lags<-0.05),50)),overlap_idx{end}(lags<0.05&lags>-0.05),overlap_idx{end}(lags>0.05&lags<prctile(lags(lags>0.05),50)),overlap_idx{end}(lags>prctile(lags(lags>0.05),50))};

lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
% lag_thresholds = [-0.2]
for n = 1:length(lag_thresholds)-1
    event_idx{2}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end


event_idx{3} = {all_overlap_idx{end},non_overlap_idx{end}};


lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
lag_thresholds = [-0.2 -0.1 -0.02 0.02 0.1 0.2];
for n = 1:length(lag_thresholds)-1
    event_idx{4}{n} =(all_overlap_idx{end}(lags>lag_thresholds(n)&lags <lag_thresholds(n+1)));
end

lags =all_lags{end};
% lag_thresholds = prctile(lags,[0 20 40 60 80 100]);
event_idx{5} = {all_overlap_idx{end}(lags>-0.2&lags<-0.1),all_overlap_idx{end}(lags<0.1&lags>-0.1),all_overlap_idx{end}(lags>0.1&lags<0.2),non_overlap_idx{end}};


group_name=[];
group_name{1} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Shuffled'};
group_name{2} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{3} = {'overlapping','non-overlapping','Shuffled'};% purely based on lags
% group_name{4} = {'Top 0-20% ipsi leading','Top 20-40% ipsi leading','Top 40-60% ipsi leading','Top 60-80% ipsi leading','Top 80-100% ipsi leading','Shuffled'};% purely based on lags
group_name{4} = {'-0.2 to -0.1s','-0.1 to -0.02s','-0.02 to 0.02s','0.02 to 0.1s','0.1 to 0.2s','Shuffled'};% purely based on lags
group_name{5} = {'Ipsi leading','Bilaterally synchronised','Contra leading','Non-overlapping','Shuffled'};

% group_name{4} = {'Top 50% ipsi leading','Bottom 50% ipsi leading','Bilaterally synchronised','Bottom 50% contra leading','Top 50% contra leading','Shuffled'};% purely based on lags exlcuding cluster 2
% group_name{5} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{6} = {'Top 50% ipsi dominant','Bottom 50% ipsi dominant','Bilaterally synchronised','Bottom 50% contra dominant','Top 50% contra dominant','Shuffled'};% clusters corr
% group_name{7} = {'Ipsi dominant ipsi leading','Bilaterally synchronised','Contra dominant contra leading','Shuffled'};

title_names = {'Ipsi-contra spindles MUA by three lags (full windows)', 'Ipsi-contra spindles MUA by bilateral thresholded lags (full windows)',...
    'Ipsi-contra spindles MUA by overlapping and non-overlapping','Ipsi-contra spindles MUA by 5 lags (200ms windows)','Ipsi-contra spindles MUA with non-overlapping'}

%%%%%%%%%% Calculate bootstrapped MUA

time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;

MUA_PSTH_merged_thresholded.x = x;
MUA_PSTMUA_PSTH_merged_thresholdedH_merged.ipsi_spindles_V1 = [];
MUA_PSTH_merged_thresholded.contra_spindles_V1 = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_V1 = [];

MUA_PSTH_merged_thresholded.ipsi_spindles_HPC = [];
MUA_PSTH_merged_thresholded.contra_spindles_HPC = [];
MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_HPC = [];

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

        MUA_PSTH_merged_thresholded.ipsi_spindles_V1{ngroup}{i} = temp1;
        MUA_PSTH_merged_thresholded.contra_spindles_V1{ngroup}{i} = temp2;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_V1{ngroup}{i} = temp3;
        MUA_PSTH_merged_thresholded.ipsi_spindles_HPC{ngroup}{i} = temp4;
        MUA_PSTH_merged_thresholded.contra_spindles_HPC{ngroup}{i} = temp5;
        MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_HPC{ngroup}{i} = temp6;
    end
end

MUA_PSTH_merged_thresholded.ipsi_spindles_V1_baseline = ipsi_baseline_bootstrap;
MUA_PSTH_merged_thresholded.contra_spindles_V1_baseline = contra_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_V1_baseline = ipsi_contra_diff_baseline_bootstrap;
MUA_PSTH_merged_thresholded.ipsi_spindles_HPC_baseline = ipsi_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.contra_spindles_HPC_baseline = contra_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_HPC_baseline = ipsi_contra_diff_baseline_bootstrap1;
MUA_PSTH_merged_thresholded.spindles_groups = [title_names];
MUA_PSTH_merged_thresholded.spindles_index = [event_idx];

%%%%%%%%%%%%% Distribution of detection thresholded events groups

customColors = [74,20,134;228,42,168;0,90,50]/256; % dark purple, magenta and dark green


fig = figure('Color','w');
fig.Position = [30 60 1880 900];
fig.Name = 'spindles transition lag vs corr vs plv distribution with overlapping or non-overlapping event clusters';
hemisphere_texts = {'Left spindles transition','Right spindles transition'}



% colour_lines = [0,90,50;74,20,134]/256; % Green Purple
for ngroup = 1:length(event_idx)
    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    elseif ngroup ==3
        colour_lines = [228,42,168;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end

    nexttile
    clear PLOT
    for k = 1:length(event_idx{ngroup})
        X = [lag_diff(event_idx{ngroup}{k}) corr_diff(event_idx{ngroup}{k}) plv_diff(event_idx{ngroup}{k})];
        hold on
        PLOT(k) = scatter3( ...
            X(:, 1), X(:, 2), X(:, 3), ...
            20, 'MarkerFaceColor',colour_lines(k, :),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1 ...
            );
    end
    xlim([-0.8 0.8])
    view([-20 60]); % <-- Add this line for better 3D perspective
    grid on; % <-- Add this line
    xlabel('ipsi-contra lag')
    ylabel('ipsi-contra corr diff')
    zlabel('ipsi-contra plv diff')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    title(title_names{ngroup})
    legend(PLOT(1:end),{group_name{ngroup}{1:end-1}},'box','off')
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')

%%%%%%%%%%%%%% plot spindles
for ngroup = 1:length(event_idx)
    fig = figure('Color','w');
    fig.Position = [350 59 1650 930];
    fig.Name =title_names{ngroup};

    if ngroup ==1
        colour_lines = [0,90,50;228,42,168;74,20,134]/256; % Dark Green, Magenta, dark purple
    elseif ngroup ==5
        colour_lines = [0,90,50;228,42,168;74,20,134;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    elseif ngroup ==3
        colour_lines = [228,42,168;82,82,82]/256; % Dark Green, Magenta, dark purple and gray
    else
        colour_lines = [0,90,50;65,171,93;228,42,168;128,125,186;74,20,134]/256; % Dark Green, Light Geen, Magenta, light purple, dark purple
    end
    nexttile
    clear ERROR_SHADE
    for i = 1:length(event_idx{ngroup})

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_spindles_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_spindles_V1_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.contra_spindles_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_spindles_V1_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_V1{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_V1_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_spindles_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_spindles_HPC_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.contra_spindles_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.contra_spindles_HPC_baseline;
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

        binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_HPC{ngroup}{i};

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
    binnedArray = MUA_PSTH_merged_thresholded.ipsi_contra_diff_spindles_HPC_baseline;
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



%%%%%%%%%%%%%%%%%%%%%%
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])




save(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_PSTH_merged_thresholded.mat'),'MUA_PSTH_merged_thresholded');

