%% Corticohippocampal reactivation during ripples

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
% load(fullfile(analysis_folder,'ripples_all_best_V1_SO_POST.mat'))

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


sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);


load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));


%% Extract key information for DOWN UP, ripples and spindles
 
% load((fullfile(analysis_folder,'cortex_SO_ref_shank_all.mat')),'cortex_SO_ref_shank_all');
%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);

        % cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession}==cortex_SO_ref_shank_all(nsession,probe_no) ...
        %     & slow_waves_all(probe_no).probe_hemisphere{nsession}==probe_no);

        % cortex_ref_shank(nsession,probe_no) = ripples_all;
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end


%%%%%%%%%%%%%%%%%% Ripple info

% if merge events to bilateral events as one event
[event_ids_first,event_ids_second] = merge_bilateral_ripple_events(merged_event_info.ripples_hemisphere_id,merged_event_info.ripples_peaktimes,0.05);


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'));
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%% ripple power
ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];
ripple_info.ripple_power = mean([ripple_info.ripple_power(event_ids_first) ripple_info.ripple_power(event_ids_second)],2);


%%% spindle co-occurance
[~,spindle_index,~,index] =RestrictInts([merged_event_info.ripples_ints(:,1) merged_event_info.ripples_ints(:,2)],merged_event_info.spindles_ints(merged_event_info.spindles_hemisphere_id==1,:));
ripple_info.spindle_presence = zeros(size(spindle_index,1),2);
ripple_info.spindle_presence(spindle_index,1) = 1;

[~,spindle_index,~,index] =RestrictInts([merged_event_info.ripples_ints(:,1) merged_event_info.ripples_ints(:,2)],merged_event_info.spindles_ints(merged_event_info.spindles_hemisphere_id==2,:));
% ripple_info.spindle_presence = zeros(size(spindle_index,1),2);
ripple_info.spindle_presence(spindle_index,2) = 1;
ripple_info.spindle_presence = ripple_info.spindle_presence(event_ids_first,:);


ripple_info.spindle_presence_PRE = zeros(size(merged_event_info.ripples_ints,1),2);
for hemi = 1:2
    ripple_interval = [merged_event_info.ripples_ints(:,1)-0.2 merged_event_info.ripples_ints(:,1)];
    spindle_interval = [merged_event_info.spindles_ints(merged_event_info.spindles_hemisphere_id==hemi,:)];

    for ii = 1:size(ripple_interval, 1)
        % Find the first index in restrictints that overlaps with the current int
        idx = find(ripple_interval(ii,1) <= spindle_interval(:,2) & ripple_interval(ii,2) >= spindle_interval(:,1), 1, 'first');

        if ~isempty(idx)
            ripple_info.spindle_presence_PRE(ii,hemi) = idx;
        end
    end
end
ripple_info.spindle_presence_PRE = ripple_info.spindle_presence_PRE(event_ids_first,:);


ripple_info.spindle_presence_POST = zeros(size(merged_event_info.ripples_ints,1),2);
for hemi = 1:2
    ripple_interval = [merged_event_info.ripples_ints(:,1) merged_event_info.ripples_ints(:,1)+0.2];
    spindle_interval = [merged_event_info.spindles_ints(merged_event_info.spindles_hemisphere_id==hemi,:)];

    for ii = 1:size(ripple_interval, 1)
        % Find the first index in restrictints that overlaps with the current int
        idx = find(ripple_interval(ii,1) <= spindle_interval(:,2) & ripple_interval(ii,2) >= spindle_interval(:,1), 1, 'first');

        if ~isempty(idx)
            ripple_info.spindle_presence_POST(ii,hemi) = idx;
        end
    end
end
ripple_info.spindle_presence_POST = ripple_info.spindle_presence_POST(event_ids_first,:);

% ripple_info.spindle_presence_hemi = ripple_info.spindle_presence_hemi(event_ids_first);

%%% spindle phase
spindle_phase=[];
for probe_no = 1:2
    spindle_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_phase = [spindle_phase{1}(:,ripples_all(1).SWS_index==1) spindle_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_phase = spindle_phase';
ripple_info.spindle_phase = ripple_info.spindle_phase(event_ids_first,:);

%%% spindle power
spindle_amplitude=[];
for probe_no = 1:2
    spindle_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_amplitude{probe_no} = [spindle_amplitude{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_amplitude = [spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_amplitude = spindle_amplitude';
ripple_info.spindle_amplitude = ripple_info.spindle_amplitude(event_ids_first,:);

%%% SO phase
SO_phase=[];
for probe_no = 1:2
    SO_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
        % SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_phase = [SO_phase{1}(:,ripples_all(1).SWS_index==1) SO_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.SO_phase = SO_phase';
ripple_info.SO_phase = ripple_info.SO_phase(event_ids_first,:);

%%% SO power
SO_amplitude=[];
for probe_no = 1:2
    SO_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_amplitude{probe_no} = [SO_amplitude{probe_no} ripples_all(probe_no).SO_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_amplitude = [SO_amplitude{1}(:,ripples_all(1).SWS_index==1) SO_amplitude{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.SO_amplitude = SO_amplitude';
ripple_info.SO_amplitude = ripple_info.SO_amplitude(event_ids_first,:);



%%%%%%
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));

% merged_event_info.ripples_hemisphere_id
% singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%% RRR
%% RRR 
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','RRR_reactivation_ripples_PSTH.mat'))
timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;

z_bias = RRR_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = RRR_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';
z_bias = z_bias(:,event_ids_first);
z_bias_V1 = z_bias_V1(:,event_ids_first);
z_bias_RRR = z_bias;
z_bias_V1_RRR = z_bias_V1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.2;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


subject_id = subject_id(event_ids_first);
session_count = session_count(event_ids_first);

singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));


%% Save into tables RRR

%%%%%% mixed effect model (SO phase)
bins_to_use = bin_centers>0 & bin_centers<0.1;
% bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_select = bin_centers>0 & bin_centers<0.2;

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>-0.1 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.1;
[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

mean_bias_V1 = mean(z_bias_V1(bin_centers>0 & bin_centers<0.2, :), 'omitnan');
mean_bias_V1_PRE = mean(z_bias_V1(bin_centers>-0.2 & bin_centers<0,:), 'omitnan');

% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds_PRE = prctile(abs(mean_bias_V1_PRE), 0:10:100);
thresholds_HPC = prctile(abs(mean_bias), 0:10:100);

thresholds = thresholds(1:end-1);
thresholds_PRE = thresholds_PRE(1:end-1);
thresholds_HPC = thresholds_HPC(1:end-1);

% spindle_power_log_odds

mean_beta = [];
SO_phase_modulation_lme = struct();


ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];
% ripple_info.ripple_power = mean([ripple_info.ripple_power(event_ids_first) ripple_info.ripple_power(event_ids_second)],2);
ripple_info.ripple_power = max([ripple_info.ripple_power(event_ids_first) ripple_info.ripple_power(event_ids_second)]')';

for i = 1:1


    th = thresholds_HPC(i);
    t1 = find(mean_bias >= th);
    t2 = find(mean_bias <= -th);
        selected_events = [t1 t2];
%     selected_events1 = [t1 t2];
%     t11 = t1;
%     t22 = t2;

    th = thresholds_PRE(i);
    t1 = find(mean_bias_V1_PRE >= th);
    t2 = find(mean_bias_V1_PRE <= -th);
    selected_events = [t1 t2];

    th = thresholds(i);
    t1 = find(mean_bias_V1 >= th);
    t2 = find(mean_bias_V1 <= -th);
    selected_events = [t1 t2];


%     [~, indices] = ismember(selected_events, selected_events1);

    % %

    tbl = table(normalize([ripple_info.ripple_power(t1); ripple_info.ripple_power(t2)],'zscore'),...
        [ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
        [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
        normalize([ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],'zscore'),...
        normalize([ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],'zscore'),...
        [ripple_info.spindle_presence(t1,2); ripple_info.spindle_presence(t2,1)],...
        [ripple_info.spindle_presence(t1,1); ripple_info.spindle_presence(t2,2)],...
        normalize(mean_bias_V1_PRE(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),normalize(mean_bias(selected_events)','zscore'),categorical([ones(length(t1),1); 2*ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),categorical(session_count(selected_events)),'VariableNames',{'RipplePower','SOPhase_Match','SOPhase_NonMatch','SpindlePower_Match','SpindlePower_NonMatch',...
        'Spindle_Match','Spindle_NonMatch','V1_logodds_PRE','V1_logodds','HPC_logodds','V1_trackID','AnimalID','SessionID'});

    tbl.Match_trough = tbl.SOPhase_Match < -pi/2 | tbl.SOPhase_Match > pi/2 ;
    tbl.NonMatch_trough = tbl.SOPhase_NonMatch < -pi/2 | tbl.SOPhase_NonMatch > pi/2 ;
    tbl.JointState = (2 * tbl.Match_trough + tbl.NonMatch_trough);
    tbl.JointState = categorical(tbl.JointState);
    %         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_PRE1.csv');
    %         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_POST1.csv');

   
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_RRR.mat'),'tbl');
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_PRE1_RRR.mat'),'tbl');
    save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR.mat'),'tbl');
%     load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.mat'),'tbl');

%     writetable(tbl, 'V1_HPC_reactivation_coherence_lme1.csv');

    tbl = table([ripple_info.ripple_power(t1); ripple_info.ripple_power(t2)], ...
        [ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
        [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
        [ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],...
        [ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],...
        [ripple_info.spindle_presence(t1,2); ripple_info.spindle_presence(t2,1)],...
        [ripple_info.spindle_presence(t1,1); ripple_info.spindle_presence(t2,2)],...
        mean_bias_V1_PRE(selected_events)',mean_bias_V1(selected_events)',mean_bias(selected_events)',categorical([ones(length(t1),1); 2*ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),categorical(session_count(selected_events)),'VariableNames',{'RipplePower','SOPhase_Match','SOPhase_NonMatch','SpindlePower_Match','SpindlePower_NonMatch',...
        'Spindle_Match','Spindle_NonMatch','V1_logodds_PRE','V1_logodds','HPC_logodds','V1_trackID','AnimalID','SessionID'});

    tbl.Match_trough = tbl.SOPhase_Match < -pi/2 | tbl.SOPhase_Match > pi/2 ;
    tbl.NonMatch_trough = tbl.SOPhase_NonMatch < -pi/2 | tbl.SOPhase_NonMatch > pi/2 ;
    tbl.JointState = (2 * tbl.Match_trough + tbl.NonMatch_trough);
    tbl.JointState = categorical(tbl.JointState);
%         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_PRE1_raw.csv');
%         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_POST1_raw.csv');
%         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_HPC1_raw.csv');
%     writetable(tbl, 'V1_HPC_reactivation_coherence_lme1_raw.csv');
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_RRR_raw.mat'),'tbl');
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_PRE1_RRR_raw.mat'),'tbl');
    save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR_raw.mat'),'tbl');
end

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_RRR_raw.mat'));tbl1 = tbl;
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR_raw.mat'));
% tbl1 = load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_raw.mat'));
[~,indices]=ismember(tbl1.V1_logodds,tbl.V1_logodds);
% [~,indices1]=ismember(tbl1.HPC_logodds,tbl.HPC_logodds);
% tbl1.V1_logodds
tbl.SpindlePower_Match(indices)= tbl1.SpindlePower_Match;
tbl.SpindlePower_NonMatch(indices)= tbl1.SpindlePower_NonMatch;
tbl.Spindle_Match(indices)= tbl1.Spindle_Match;
tbl.Spindle_NonMatch(indices)= tbl1.Spindle_NonMatch;
% writetable(tbl, 'V1_HPC_reactivation_coherence_lme1_raw.csv');
writetable(tbl, fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR_raw.csv'));
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR_raw.mat'),'tbl');
tbl1 = tbl;


% tbl1 = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1.csv'));tbl1 = tbl;
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR.mat'));
% tbl1 = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.csv'));

tbl.SpindlePower_Match= normalize(tbl1.SpindlePower_Match);
tbl.SpindlePower_NonMatch= normalize(tbl1.SpindlePower_NonMatch);
tbl.Spindle_Match= tbl1.Spindle_Match;
tbl.Spindle_NonMatch= tbl1.Spindle_NonMatch;
writetable(tbl, fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR.csv'));
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR.mat'),'tbl');


%% Bayesian 
%%%%%%%%%%%%%%%%%%
load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))

timebin = 0.02;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;

nan_mask = KDE_reactivation_ripples_PSTH.nan_mask(:,1:2:end)'; %KDE nan mask 10ms bin, make it 20ms bin

z_bias = [bayesian_reactivation_all(1).z_log_odds(:,ripples_all(1).SWS_index==1) bayesian_reactivation_all(2).z_log_odds(:,ripples_all(2).SWS_index==1)];
z_bias_V1 = [bayesian_reactivation_V1_all(1).z_log_odds(:,ripples_all(1).SWS_index==1) bayesian_reactivation_V1_all(2).z_log_odds(:,ripples_all(2).SWS_index==1)];

z_bias = z_bias + nan_mask;
z_bias_V1 = z_bias_V1 + nan_mask;

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';
z_bias = z_bias(:,event_ids_first);
z_bias_V1 = z_bias_V1(:,event_ids_first);
z_bias_bayesian = z_bias;
z_bias_V1_bayesian = z_bias_V1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.2;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


subject_id = subject_id(event_ids_first);
session_count = session_count(event_ids_first);

singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));


%% Save into tables Bayesian

%%%%%% mixed effect model (SO phase)
bins_to_use = bin_centers>0 & bin_centers<0.1;
% bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_select = bin_centers>0 & bin_centers<0.2;

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>-0.1 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.1;
[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

mean_bias_V1 = mean(z_bias_V1(bin_centers>0 & bin_centers<0.2, :), 'omitnan');
mean_bias_V1_PRE = mean(z_bias_V1(bin_centers>-0.2 & bin_centers<0,:), 'omitnan');

% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds_PRE = prctile(abs(mean_bias_V1_PRE), 0:10:100);
thresholds_HPC = prctile(abs(mean_bias), 0:10:100);

thresholds = thresholds(1:end-1);
thresholds_PRE = thresholds_PRE(1:end-1);
thresholds_HPC = thresholds_HPC(1:end-1);

% spindle_power_log_odds

mean_beta = [];
SO_phase_modulation_lme = struct();


ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];
% ripple_info.ripple_power = mean([ripple_info.ripple_power(event_ids_first) ripple_info.ripple_power(event_ids_second)],2);
ripple_info.ripple_power = max([ripple_info.ripple_power(event_ids_first) ripple_info.ripple_power(event_ids_second)]')';

for i = 1:1


    th = thresholds_HPC(i);
    t1 = find(mean_bias >= th);
    t2 = find(mean_bias <= -th);
    selected_events = [t1 t2];
    %     selected_events1 = [t1 t2];
    %     t11 = t1;
    %     t22 = t2;

    th = thresholds_PRE(i);
    t1 = find(mean_bias_V1_PRE >= th);
    t2 = find(mean_bias_V1_PRE <= -th);
    selected_events = [t1 t2];

    th = thresholds(i);
    t1 = find(mean_bias_V1 >= th);
    t2 = find(mean_bias_V1 <= -th);
    selected_events = [t1 t2];


%     [~, indices] = ismember(selected_events, selected_events1);

    % %

    tbl = table(normalize([ripple_info.ripple_power(t1); ripple_info.ripple_power(t2)],'zscore'),...
        [ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
        [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
        normalize([ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],'zscore'),...
        normalize([ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],'zscore'),...
        [ripple_info.spindle_presence(t1,2); ripple_info.spindle_presence(t2,1)],...
        [ripple_info.spindle_presence(t1,1); ripple_info.spindle_presence(t2,2)],...
        normalize(mean_bias_V1_PRE(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),normalize(mean_bias(selected_events)','zscore'),categorical([ones(length(t1),1); 2*ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),categorical(session_count(selected_events)),'VariableNames',{'RipplePower','SOPhase_Match','SOPhase_NonMatch','SpindlePower_Match','SpindlePower_NonMatch',...
        'Spindle_Match','Spindle_NonMatch','V1_logodds_PRE','V1_logodds','HPC_logodds','V1_trackID','AnimalID','SessionID'});

    tbl.Match_trough = tbl.SOPhase_Match < -pi/2 | tbl.SOPhase_Match > pi/2 ;
    tbl.NonMatch_trough = tbl.SOPhase_NonMatch < -pi/2 | tbl.SOPhase_NonMatch > pi/2 ;
    tbl.JointState = (2 * tbl.Match_trough + tbl.NonMatch_trough);
    tbl.JointState = categorical(tbl.JointState);
    %         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_PRE1.csv');
    %         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_POST1.csv');

   
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_bayesian.mat'),'tbl');
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_PRE1_bayesian.mat'),'tbl');
    save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_bayesian.mat'),'tbl');
%     load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.mat'),'tbl');

%     writetable(tbl, 'V1_HPC_reactivation_coherence_lme1.csv');

    tbl = table([ripple_info.ripple_power(t1); ripple_info.ripple_power(t2)], ...
        [ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
        [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
        [ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],...
        [ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],...
        [ripple_info.spindle_presence(t1,2); ripple_info.spindle_presence(t2,1)],...
        [ripple_info.spindle_presence(t1,1); ripple_info.spindle_presence(t2,2)],...
        mean_bias_V1_PRE(selected_events)',mean_bias_V1(selected_events)',mean_bias(selected_events)',categorical([ones(length(t1),1); 2*ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),categorical(session_count(selected_events)),'VariableNames',{'RipplePower','SOPhase_Match','SOPhase_NonMatch','SpindlePower_Match','SpindlePower_NonMatch',...
        'Spindle_Match','Spindle_NonMatch','V1_logodds_PRE','V1_logodds','HPC_logodds','V1_trackID','AnimalID','SessionID'});

    tbl.Match_trough = tbl.SOPhase_Match < -pi/2 | tbl.SOPhase_Match > pi/2 ;
    tbl.NonMatch_trough = tbl.SOPhase_NonMatch < -pi/2 | tbl.SOPhase_NonMatch > pi/2 ;
    tbl.JointState = (2 * tbl.Match_trough + tbl.NonMatch_trough);
    tbl.JointState = categorical(tbl.JointState);
%         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_PRE1_raw.csv');
%         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_POST1_raw.csv');
%         writetable(tbl, 'V1_HPC_reactivation_coherence_lme_HPC1_raw.csv');
%     writetable(tbl, 'V1_HPC_reactivation_coherence_lme1_raw.csv');
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_bayesian_raw.mat'),'tbl');
%     save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_PRE1_bayesian_raw.mat'),'tbl');
    save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_bayesian_raw.mat'),'tbl');
end

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_bayesian_raw.mat'));tbl1 = tbl;
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_bayesian_raw.mat'));
% tbl1 = load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_raw.mat'));
[~,indices]=ismember(tbl1.V1_logodds,tbl.V1_logodds);
% [~,indices1]=ismember(tbl1.HPC_logodds,tbl.HPC_logodds);
% tbl1.V1_logodds
tbl.SpindlePower_Match(indices)= tbl1.SpindlePower_Match;
tbl.SpindlePower_NonMatch(indices)= tbl1.SpindlePower_NonMatch;
tbl.Spindle_Match(indices)= tbl1.Spindle_Match;
tbl.Spindle_NonMatch(indices)= tbl1.Spindle_NonMatch;
% writetable(tbl, 'V1_HPC_reactivation_coherence_lme1_raw.csv');
writetable(tbl, fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian_raw.csv'));
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian_raw.mat'),'tbl');
tbl1 = tbl;


% tbl1 = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1.csv'));tbl1 = tbl;
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_bayesian.mat'));
% tbl1 = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.csv'));

tbl.SpindlePower_Match= normalize(tbl1.SpindlePower_Match);
tbl.SpindlePower_NonMatch= normalize(tbl1.SpindlePower_NonMatch);
tbl.Spindle_Match= tbl1.Spindle_Match;
tbl.Spindle_NonMatch= tbl1.Spindle_NonMatch;
writetable(tbl, fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian.csv'));
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian.mat'),'tbl');


%% RRR Ripple power
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_RRR.mat'),'tbl');

tbl.Match_trough = categorical(tbl.Match_trough);
% tbl.Match_peak = categorical(tbl.Match_peak);
tbl.JointState = categorical(tbl.JointState);
tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
% tbl.HPC_trackID = categorical(tbl.HPC_trackID);
% 
% formula = 'HPC_logodds ~ V1_logodds * RipplePower    + (1| SessionID) +(1 | AnimalID) ';
formula = 'V1_logodds ~ HPC_logodds * RipplePower    + (RipplePower| SessionID) +(RipplePower | AnimalID) ';
% 
% % formula = 'V1_logodds ~ HPC_trackID * RipplePower    + (HPC_trackID| SessionID) +(HPC_trackID | AnimalID) ';
% % lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
lme = fitlme(tbl, formula,'DummyVarCoding','reference','CovariancePattern', {'Diagonal', 'Diagonal'});
lme.Coefficients
% 


models = {'V1_logodds_PRE ~ HPC_logodds  * RipplePower + (RipplePower| SessionID) +(RipplePower | AnimalID)  ',
    'V1_logodds ~ HPC_logodds * RipplePower + (RipplePower| SessionID)+(RipplePower | AnimalID) '};

for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

%     lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme.Coefficients

    ripple_power_modulation_lme(i).model = lme;
    ripple_power_modulation_lme(i).variable = [lme.Coefficients.Name];
    ripple_power_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    ripple_power_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    ripple_power_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    ripple_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    ripple_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    ripple_power_modulation_lme(i).marginal_R2 = marginal_R2;
end


% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
% VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
%     226, 132, 187;   % interpolated 2/3
%     212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

fig = figure('Name','Ripple power RRR reactivation coherence','Position',[482 111 665 851]);
sgtitle('Ripple power RRR reactivation coherence')
% fig = figure('Name','Ripple power KDE reactivation coherence (HPC Ripple Track)','Position',[482 111 665 851]);
% sgtitle('Ripple power KDE reactivation coherence (HPC Ripple Track)')

VariableName = {'RipplePower:HPC_'}
% VariableName = {'PRE','HPC'};
titles = {'PRE V1 (RipplePowerXHPC)','POST V1 (RipplePowerXHPC)'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    clear VariableIdx
    for n = 1:length(VariableName)
        VariableIdx(n) = find(startsWith(ripple_power_modulation_lme(i).variable, VariableName{n}));
    end
    %     end
    clear b_CI b_mean
    b_CI = [ripple_power_modulation_lme(i).b_CI];
    b_mean = [ripple_power_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0+n*0.01,sprintf('%.3e',ripple_power_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0 0.03])
    title(titles{i})
end



models = {'V1_logodds_PRE ~ HPC_logodds + (1| SessionID) +(1 | AnimalID)  ',
    'V1_logodds ~ HPC_logodds + (1| SessionID)+(1 | AnimalID) '};
titles = {'PRE V1 (HPC log odds)','POST V1 (HPC log odds)'};
for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

%     lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme.Coefficients

    log_odds_modulation_lme(i).model = lme;
    log_odds_modulation_lme(i).variable = [lme.Coefficients.Name];
    log_odds_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    log_odds_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    log_odds_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    log_odds_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    log_odds_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    log_odds_modulation_lme(i).marginal_R2 = marginal_R2;
end


VariableName = {'HPC_'}
% VariableName = {'PRE','HPC'};
titles = {'PRE V1','POST V1'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    clear VariableIdx
    for n = 1:length(VariableName)
        VariableIdx(n) = find(startsWith(log_odds_modulation_lme(i).variable, VariableName{n}));
    end
    %     end
    clear b_CI b_mean
    b_CI = [log_odds_modulation_lme(i).b_CI];
    b_mean = [log_odds_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i+4)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0+n*0.01,sprintf('%.3e',log_odds_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([0 0.08])
    title(titles{i})
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','RRR_log_odds_modulation_lme.mat'),'log_odds_modulation_lme');

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_RRR.mat'),'ripple_power_modulation_lme');

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_RRR.mat'),'ripple_power_modulation_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])



%% Bayesian Ripple power
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_bayesian.mat'),'tbl');
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HCP1_bayesian.csv'));

tbl.Match_trough = categorical(tbl.Match_trough);
% tbl.Match_peak = categorical(tbl.Match_peak);
tbl.JointState = categorical(tbl.JointState);
tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
% tbl.HPC_trackID = categorical(tbl.HPC_trackID);
% 
% formula = 'HPC_logodds ~ V1_logodds * RipplePower    + (1| SessionID) +(1 | AnimalID) ';
formula = 'V1_logodds ~ HPC_logodds * RipplePower    + (RipplePower| SessionID) +(RipplePower | AnimalID) ';
% 
% % formula = 'V1_logodds ~ HPC_trackID * RipplePower    + (HPC_trackID| SessionID) +(HPC_trackID | AnimalID) ';
% % lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
lme = fitlme(tbl, formula,'DummyVarCoding','reference','CovariancePattern', {'Diagonal', 'Diagonal'});
lme.Coefficients
% 


models = {'V1_logodds_PRE ~ HPC_trackID  * RipplePower + (1| SessionID) +(1 | AnimalID)  ',
    'V1_logodds ~ HPC_trackID * RipplePower + (1| SessionID)+(1 | AnimalID) '};

for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

%     lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme.Coefficients

    ripple_power_modulation_lme(i).model = lme;
    ripple_power_modulation_lme(i).variable = [lme.Coefficients.Name];
    ripple_power_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    ripple_power_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    ripple_power_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    ripple_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    ripple_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    ripple_power_modulation_lme(i).marginal_R2 = marginal_R2;
end


% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
% VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
%     226, 132, 187;   % interpolated 2/3
%     212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

fig = figure('Name','Ripple power Bayesian reactivation coherence','Position',[482 111 665 851]);
sgtitle('Ripple power Bayesian reactivation coherence')
% fig = figure('Name','Ripple power KDE reactivation coherence (HPC Ripple Track)','Position',[482 111 665 851]);
% sgtitle('Ripple power KDE reactivation coherence (HPC Ripple Track)')

VariableName = {'RipplePower:HPC_'}
% VariableName = {'PRE','HPC'};
titles = {'PRE V1 (RipplePowerXHPC)','POST V1 (RipplePowerXHPC)'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    clear VariableIdx
    for n = 1:length(VariableName)
        VariableIdx(n) = find(startsWith(ripple_power_modulation_lme(i).variable, VariableName{n}));
    end
    %     end
    clear b_CI b_mean
    b_CI = [ripple_power_modulation_lme(i).b_CI];
    b_mean = [ripple_power_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0+n*0.01,sprintf('%.3e',ripple_power_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.01 0.03])
    title(titles{i})
end



models = {'V1_logodds_PRE ~ HPC_logodds + (1| SessionID) +(1 | AnimalID)  ',
    'V1_logodds ~ HPC_logodds + (1| SessionID)+(1 | AnimalID) '};
titles = {'PRE V1 (HPC log odds)','POST V1 (HPC log odds)'};
for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

%     lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
    lme.Coefficients

    log_odds_modulation_lme(i).model = lme;
    log_odds_modulation_lme(i).variable = [lme.Coefficients.Name];
    log_odds_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    log_odds_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    log_odds_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    log_odds_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    log_odds_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    log_odds_modulation_lme(i).marginal_R2 = marginal_R2;
end


VariableName = {'HPC_'}
% VariableName = {'PRE','HPC'};
titles = {'PRE V1','POST V1'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    clear VariableIdx
    for n = 1:length(VariableName)
        VariableIdx(n) = find(startsWith(log_odds_modulation_lme(i).variable, VariableName{n}));
    end
    %     end
    clear b_CI b_mean
    b_CI = [log_odds_modulation_lme(i).b_CI];
    b_mean = [log_odds_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i+4)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0+n*0.01,sprintf('%.3e',log_odds_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.01 0.07])
    title(titles{i})
end
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_bayesian.mat'),'ripple_power_modulation_lme');

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_bayesian.mat'),'ripple_power_modulation_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','bayesian_log_odds_modulation_lme.mat'),'log_odds_modulation_lme');

%% RRR SO phase state

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR.mat'));

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR.mat'));
tbl.Match_trough = categorical(tbl.Match_trough);
% tbl.Match_peak = categorical(tbl.Match_peak);
tbl.JointState = categorical(tbl.JointState);
tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);

models = {'HPC_logodds ~ V1_logodds_PRE * JointState + (1| SessionID) + (1 | AnimalID)',
    'HPC_logodds ~ V1_logodds * JointState + (1| SessionID)+ (1 | AnimalID) ',
    'HPC_logodds ~ V1_logodds_PRE * Match_trough + ( 1| SessionID) + (1 | AnimalID)',
    'HPC_logodds ~ V1_logodds * Match_trough + ( 1| SessionID) + (1 | AnimalID)'}

% models = {'HPC_logodds ~ V1_logodds_PRE * JointState + (V1_logodds_PRE | SessionID) + (V1_logodds_PRE +JointState | AnimalID)',
%     'HPC_logodds ~ V1_logodds * JointState + (V1_logodds | SessionID) + (V1_logodds +JointState | AnimalID)',
%     'HPC_logodds ~ V1_logodds_PRE * Match_trough + ( V1_logodds_PRE + Match_trough| SessionID) + (V1_logodds_PRE+ Match_trough | AnimalID)',
%     'HPC_logodds ~ V1_logodds * Match_trough + ( V1_logodds + Match_trough| SessionID) + (V1_logodds+ Match_trough | AnimalID)'}

% CovariancePattern = {{'Diagonal', 'Diagonal'},{},{},{}};
SO_phase_modulation_lme = struct();
for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
%     lme = fitlme(tbl, formula,'DummyVarCoding','effects');
%     lme = fitlme(tbl, formula);

%     lme.Coefficients

    % 1. Get the Fixed Effects and their Covariance Matrix
    beta_vec = lme.fixedEffects;
    vcov_mat = lme.CoefficientCovariance; % The Var-Cov matrix of fixed effects

    % 2. Identify indices for States 0, 1, and 2 interactions
    names = lme.CoefficientNames;
    if sum(contains(names,'JointState'))>0
        idx0 = find(contains(names, ':JointState_0'));
        idx1 = find(contains(names, ':JointState_1'));
        idx2 = find(contains(names, ':JointState_2'));
        %         idx =find(contains(names, ':JointState'));

        % 3. Calculate Beta 3 (The negative sum)
        beta_3 = -(beta_vec(idx0) + beta_vec(idx1) + beta_vec(idx2));

        % 4. Calculate the Variance for Beta 3 using the Covariance Matrix
        % Formula: Var(H*beta) = H * Vcov * H'
        H = zeros(1, length(names));
        H([idx0, idx1, idx2]) = -1;
        var_3 = H * vcov_mat * H';

    else
        idx0 = find(contains(names, ':Match_trough'));
        % 3. Calculate Beta 3 (The negative sum)
        beta_3 = -(beta_vec(idx0));

        % 4. Calculate the Variance for Beta 3 using the Covariance Matrix
        % Formula: Var(H*beta) = H * Vcov * H'
        H = zeros(1, length(names));
        H([idx0]) = -1;
        var_3 = H * vcov_mat * H';
    end

    % 5. Get SE and CI
    ci_3 = [beta_3 - 1.9599*sqrt(var_3), beta_3 + 1.9599*sqrt(var_3)];
    t_stat = beta_3 / sqrt(var_3);
    p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));


    SO_phase_modulation_lme(i).model = lme;
    SO_phase_modulation_lme(i).variable = [lme.Coefficients.Name; {[lme.Coefficients.Name{end}(1:end-1) '3']}];
    SO_phase_modulation_lme(i).p = [lme.Coefficients.pValue(:); p_value];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    SO_phase_modulation_lme(i).b = [lme.Coefficients.Estimate(:); beta_3];
    SO_phase_modulation_lme(i).t = [lme.Coefficients.tStat(:); t_stat];
    SO_phase_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:); ci_3];
    SO_phase_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    SO_phase_modulation_lme(i).marginal_R2 = marginal_R2;
end

% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
VariableIdx = find(contains(SO_phase_modulation_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

fig = figure('Name','SO phase state RRR reactivation coherence','Position',[482 111 665 851]);
sgtitle('SO phase state RRR reactivation coherence')
% fig = figure('Name','SO phase state KDE reactivation coherence (complex random models)','Position',[482 111 665 851]);
% sgtitle('SO phase state KDE reactivation coherence (complex random models)')

titles = {'PRE V1 log odds','POST V1 log odds','PRE V1 log odds','POST V1 log odds'}
for i = 1:4
    if sum(contains(SO_phase_modulation_lme(i).variable, ':JointState'))>0
        VariableIdx = find(contains(SO_phase_modulation_lme(i).variable, ':JointState'));
    else
        VariableIdx = find(contains(SO_phase_modulation_lme(i).variable, ':Match_trough'));
    end
    clear b_CI b_mean
    b_CI = [SO_phase_modulation_lme(i).b_CI];
    b_mean = [SO_phase_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0.02+n*0.02,sprintf('%.3e',SO_phase_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
%     ylim([-0.07 0.1])
    ylim([-0.07 0.06])
    title(titles{i})
end

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_complex_random.mat'),'SO_phase_modulation_lme');
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_RRR.mat'),'SO_phase_modulation_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])


%% Bayesian SO phase state

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_bayesian.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_bayesian.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR.mat'));
tbl.Match_trough = categorical(tbl.Match_trough);
% tbl.Match_peak = categorical(tbl.Match_peak);
tbl.JointState = categorical(tbl.JointState);
tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);

% models = {'HPC_logodds ~ V1_logodds_PRE * JointState + (1| SessionID) + (1 | AnimalID)',
%     'HPC_logodds ~ V1_logodds * JointState + (1| SessionID)+ (1 | AnimalID) ',
%     'HPC_logodds ~ V1_logodds_PRE * Match_trough + ( 1| SessionID) + (1 | AnimalID)',
%     'HPC_logodds ~ V1_logodds * Match_trough + ( 1| SessionID) + (1 | AnimalID)'}

models = {'HPC_logodds ~ V1_logodds_PRE * JointState + (V1_logodds_PRE | SessionID) + (V1_logodds_PRE +JointState | AnimalID)',
    'HPC_logodds ~ V1_logodds * JointState + (V1_logodds | SessionID) + (V1_logodds +JointState | AnimalID)',
    'HPC_logodds ~ V1_logodds_PRE * Match_trough + ( V1_logodds_PRE + Match_trough| SessionID) + (V1_logodds_PRE+ Match_trough | AnimalID)',
    'HPC_logodds ~ V1_logodds * Match_trough + ( V1_logodds + Match_trough| SessionID) + (V1_logodds+ Match_trough | AnimalID)'}

% CovariancePattern = {{'Diagonal', 'Diagonal'},{},{},{}};
SO_phase_modulation_lme = struct();
for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
%     lme = fitlme(tbl, formula,'DummyVarCoding','effects');
%     lme = fitlme(tbl, formula);

%     lme.Coefficients

    % 1. Get the Fixed Effects and their Covariance Matrix
    beta_vec = lme.fixedEffects;
    vcov_mat = lme.CoefficientCovariance; % The Var-Cov matrix of fixed effects

    % 2. Identify indices for States 0, 1, and 2 interactions
    names = lme.CoefficientNames;
    if sum(contains(names,'JointState'))>0
        idx0 = find(contains(names, ':JointState_0'));
        idx1 = find(contains(names, ':JointState_1'));
        idx2 = find(contains(names, ':JointState_2'));
        %         idx =find(contains(names, ':JointState'));

        % 3. Calculate Beta 3 (The negative sum)
        beta_3 = -(beta_vec(idx0) + beta_vec(idx1) + beta_vec(idx2));

        % 4. Calculate the Variance for Beta 3 using the Covariance Matrix
        % Formula: Var(H*beta) = H * Vcov * H'
        H = zeros(1, length(names));
        H([idx0, idx1, idx2]) = -1;
        var_3 = H * vcov_mat * H';

    else
        idx0 = find(contains(names, ':Match_trough'));
        % 3. Calculate Beta 3 (The negative sum)
        beta_3 = -(beta_vec(idx0));

        % 4. Calculate the Variance for Beta 3 using the Covariance Matrix
        % Formula: Var(H*beta) = H * Vcov * H'
        H = zeros(1, length(names));
        H([idx0]) = -1;
        var_3 = H * vcov_mat * H';
    end

    % 5. Get SE and CI
    ci_3 = [beta_3 - 1.9599*sqrt(var_3), beta_3 + 1.9599*sqrt(var_3)];
    t_stat = beta_3 / sqrt(var_3);
    p_value = 2 * (1 - tcdf(abs(t_stat), lme.DFE));


    SO_phase_modulation_lme(i).model = lme;
    SO_phase_modulation_lme(i).variable = [lme.Coefficients.Name; {[lme.Coefficients.Name{end}(1:end-1) '3']}];
    SO_phase_modulation_lme(i).p = [lme.Coefficients.pValue(:); p_value];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    SO_phase_modulation_lme(i).b = [lme.Coefficients.Estimate(:); beta_3];
    SO_phase_modulation_lme(i).t = [lme.Coefficients.tStat(:); t_stat];
    SO_phase_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:); ci_3];
    SO_phase_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    SO_phase_modulation_lme(i).marginal_R2 = marginal_R2;
end

% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
VariableIdx = find(contains(SO_phase_modulation_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;
% fig = figure('Name','SO phase state Bayesian reactivation coherence','Position',[482 111 665 851]);
% sgtitle('SO phase state Bayesian reactivation coherence')
fig = figure('Name','SO phase state Bayesian reactivation coherence HPC','Position',[482 111 665 851]);
sgtitle('SO phase state Bayesian reactivation coherence HPC')
% fig = figure('Name','SO phase state KDE reactivation coherence (complex random models)','Position',[482 111 665 851]);
% sgtitle('SO phase state KDE reactivation coherence (complex random models)')

titles = {'PRE V1 log odds','POST V1 log odds','PRE V1 log odds','POST V1 log odds'}
for i = 1:4
    if sum(contains(SO_phase_modulation_lme(i).variable, ':JointState'))>0
        VariableIdx = find(contains(SO_phase_modulation_lme(i).variable, ':JointState'));
    else
        VariableIdx = find(contains(SO_phase_modulation_lme(i).variable, ':Match_trough'));
    end
    clear b_CI b_mean
    b_CI = [SO_phase_modulation_lme(i).b_CI];
    b_mean = [SO_phase_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0.02+n*0.02,sprintf('%.3e',SO_phase_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
%     ylim([-0.07 0.1])
    ylim([-0.07 0.06])
    title(titles{i})
end

% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_complex_random.mat'),'SO_phase_modulation_lme');
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_bayesian.mat'),'SO_phase_modulation_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_bayesian_HPC.mat'),'SO_phase_modulation_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])


%% RRR Spindle power
%%%%%%
%%%%%%
%%%%%% PRE and HPC -> POST
spindle_power_modulation_lme = struct();
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1.csv'));
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.csv'));

% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1.csv'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_RRR.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.mat'));
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_PRE1.csv'));


tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
tbl.Spindle_Match = categorical(tbl.Spindle_Match);
tbl.Spindle_NonMatch = categorical(tbl.Spindle_NonMatch);
% tbl.HPC_trackID = categorical(tbl.HPC_trackID);
% tbl.V1_trackID_PRE = categorical(tbl.V1_trackID_PRE);

models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + (1| SessionID) + (1| AnimalID)',
    'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + (1 | SessionID) + (1 | AnimalID)'};

% models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + (V1_logodds_PRE * (SpindlePower_Match + SpindlePower_NonMatch)| SessionID) + (V1_logodds_PRE * (SpindlePower_Match + SpindlePower_NonMatch)| AnimalID)',
%     'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + (V1_logodds * (SpindlePower_Match + SpindlePower_NonMatch) | SessionID) + (V1_logodds * (SpindlePower_Match + SpindlePower_NonMatch) | AnimalID)'};

% models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + (V1_logodds_PRE + SpindlePower_Match + SpindlePower_NonMatch| SessionID) + (V1_logodds_PRE + SpindlePower_Match + SpindlePower_NonMatch| AnimalID)',
%     'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + (V1_logodds +SpindlePower_Match + SpindlePower_NonMatch | SessionID) + (V1_logodds +SpindlePower_Match + SpindlePower_NonMatch | AnimalID)'};

for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
%     lme = fitlme(tbl, formula)
    lme.Coefficients
% [Psi, ~, ~] = covarianceParameters(lme);
% Psi_matrix = Psi{1}
    spindle_power_modulation_lme(i).model = lme;
    spindle_power_modulation_lme(i).variable = [lme.Coefficients.Name];
    spindle_power_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    spindle_power_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    spindle_power_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    spindle_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    spindle_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    spindle_power_modulation_lme(i).marginal_R2 = marginal_R2;
end



% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
% VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
%     226, 132, 187;   % interpolated 2/3
%     212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

fig = figure('Name','Spindle power RRR reactivation coherence','Position',[482 111 665 851]);
sgtitle('Spindle power RRR reactivation coherence')

VariableName = {'SpindlePower_NonMatch:V1','SpindlePower_Match:V1'}
% VariableName = {'PRE','HPC'};
titles = {'PRE V1','POST V1'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    clear VariableIdx
    for n = 1:length(VariableName)
    VariableIdx(n) = find(contains(spindle_power_modulation_lme(i).variable, VariableName{n}));
    end
    %     end
    clear b_CI b_mean
    b_CI = [spindle_power_modulation_lme(i).b_CI];
    b_mean = [spindle_power_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0+n*0.01,sprintf('%.3e',spindle_power_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.025 0.025])
    title(titles{i})
end


save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_RRR.mat'),'spindle_power_modulation_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])



%% Bayesian Spindle power
%%%%%%
%%%%%%
%%%%%% PRE and HPC -> POST
spindle_power_modulation_lme = struct();
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1.csv'));
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.csv'));

% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1.csv'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_POST1_bayesian.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.mat'));
% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_PRE1.csv'));
% 

tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
tbl.Spindle_Match = categorical(tbl.Spindle_Match);
tbl.Spindle_NonMatch = categorical(tbl.Spindle_NonMatch);


models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + (1| SessionID) + (1| AnimalID)',
    'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + (1 | SessionID) + (1 | AnimalID)'};

% models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + (V1_logodds_PRE * (SpindlePower_Match + SpindlePower_NonMatch)| SessionID) + (V1_logodds_PRE * (SpindlePower_Match + SpindlePower_NonMatch)| AnimalID)',
%     'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + (V1_logodds * (SpindlePower_Match + SpindlePower_NonMatch) | SessionID) + (V1_logodds * (SpindlePower_Match + SpindlePower_NonMatch) | AnimalID)'};

% models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + (V1_logodds_PRE + SpindlePower_Match + SpindlePower_NonMatch| SessionID) + (V1_logodds_PRE + SpindlePower_Match + SpindlePower_NonMatch| AnimalID)',
%     'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + (V1_logodds +SpindlePower_Match + SpindlePower_NonMatch | SessionID) + (V1_logodds +SpindlePower_Match + SpindlePower_NonMatch | AnimalID)'};

for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
%     lme = fitlme(tbl, formula)
    lme.Coefficients
% [Psi, ~, ~] = covarianceParameters(lme);
% Psi_matrix = Psi{1}
    spindle_power_modulation_lme(i).model = lme;
    spindle_power_modulation_lme(i).variable = [lme.Coefficients.Name];
    spindle_power_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    spindle_power_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    spindle_power_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    spindle_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    spindle_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    spindle_power_modulation_lme(i).marginal_R2 = marginal_R2;
end



% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
% VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
%     226, 132, 187;   % interpolated 2/3
%     212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

fig = figure('Name','Spindle power Bayesian reactivation coherence','Position',[482 111 665 851]);
sgtitle('Spindle power Bayesian reactivation coherence')

VariableName = {'SpindlePower_NonMatch:V1','SpindlePower_Match:V1'}
% VariableName = {'PRE','HPC'};
titles = {'PRE V1','POST V1'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    clear VariableIdx
    for n = 1:length(VariableName)
    VariableIdx(n) = find(contains(spindle_power_modulation_lme(i).variable, VariableName{n}));
    end
    %     end
    clear b_CI b_mean
    b_CI = [spindle_power_modulation_lme(i).b_CI];
    b_mean = [spindle_power_modulation_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0+n*0.01,sprintf('%.3e',spindle_power_modulation_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.025 0.025])
    title(titles{i})
end


save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_bayesian.mat'),'spindle_power_modulation_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])


%% RRR SO and spindle interaction
SO_spindle_modulation_lme = struct();

% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.csv'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_RRR.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_RRR.mat'));
tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
tbl.Spindle_Match = categorical(tbl.Spindle_Match);
tbl.Spindle_NonMatch = categorical(tbl.Spindle_NonMatch);
tbl.Match_trough = categorical(tbl.Match_trough);
% tbl.Match_peak = categorical(tbl.Match_peak);
tbl.JointState = categorical(tbl.JointState);

% models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + V1_logodds_PRE * JointState + (1| SessionID) + (1| AnimalID)',
%     'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + V1_logodds * JointState + (1| SessionID) + (1| AnimalID)'};
models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match * Match_trough + (1| SessionID) + (1| AnimalID)',
    'HPC_logodds ~ V1_logodds * SpindlePower_Match * Match_trough + (1| SessionID) + (1| AnimalID)'};

for i = 1:length(models)
    formula = models{i};

%     lme = fitlme(tbl, formula,'CovariancePattern', {'Diagonal', 'Diagonal'});
    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
%     lme = fitlme(tbl, formula)
    lme.Coefficients

    SO_spindle_modulation_lme(i).model = lme;
    SO_spindle_modulation_lme(i).variable = [lme.Coefficients.Name];
    SO_spindle_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    SO_spindle_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    SO_spindle_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    SO_spindle_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    SO_spindle_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    SO_spindle_modulation_lme(i).marginal_R2 = marginal_R2;
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_spindle_modulation_lme_RRR.mat'),'SO_spindle_modulation_lme');

%% Bayesian SO and spindle interaction
SO_spindle_modulation_lme = struct();

% tbl = readtable(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1.csv'));
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme_HPC1_bayesian.mat'));
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','V1_HPC_reactivation_coherence_lme1_bayesian.mat'));
tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
tbl.Spindle_Match = categorical(tbl.Spindle_Match);
tbl.Spindle_NonMatch = categorical(tbl.Spindle_NonMatch);
tbl.Match_trough = categorical(tbl.Match_trough);
% tbl.Match_peak = categorical(tbl.Match_peak);
tbl.JointState = categorical(tbl.JointState);

% models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match + V1_logodds_PRE * SpindlePower_NonMatch + V1_logodds_PRE * JointState + (1| SessionID) + (1| AnimalID)',
%     'HPC_logodds ~ V1_logodds * SpindlePower_Match + V1_logodds * SpindlePower_NonMatch + V1_logodds * JointState + (1| SessionID) + (1| AnimalID)'};
models = {'HPC_logodds ~ V1_logodds_PRE * SpindlePower_Match * Match_trough + (1| SessionID) + (1| AnimalID)',
    'HPC_logodds ~ V1_logodds * SpindlePower_Match * Match_trough + (1| SessionID) + (1| AnimalID)'};

for i = 1:length(models)
    formula = models{i};

%     lme = fitlme(tbl, formula,'CovariancePattern', {'Diagonal', 'Diagonal'});
    lme = fitlme(tbl, formula,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
%     lme = fitlme(tbl, formula)
    lme.Coefficients

    SO_spindle_modulation_lme(i).model = lme;
    SO_spindle_modulation_lme(i).variable = [lme.Coefficients.Name];
    SO_spindle_modulation_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    SO_spindle_modulation_lme(i).b = [lme.Coefficients.Estimate(:)];
    SO_spindle_modulation_lme(i).t = [lme.Coefficients.tStat(:)];
    SO_spindle_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    SO_spindle_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    SO_spindle_modulation_lme(i).marginal_R2 = marginal_R2;
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_spindle_modulation_lme.mat'),'SO_spindle_modulation_lme');

%% PRE and HPC log odds -> POST
%%%%%%
%%%%%%
%%%%%% PRE and HPC -> POST
PRE_POST_V1_log_odds_lme = struct();
% models = {['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
%         '+ (1 | AnimalID) + (V1_logodds_PRE - 1 | AnimalID) + (JointState - 1 | AnimalID) ' ...
%         '+ (1 | SessionID) + (V1_logodds_PRE - 1 | SessionID) + (JointState - 1 | SessionID)'],...
%         ['HPC_logodds ~ V1_logodds * JointState ' ...
%         '+ (1 | AnimalID) + (V1_logodds - 1 | AnimalID) + (JointState - 1 | AnimalID) ' ...
%         '+ (1 | SessionID) + (V1_logodds - 1 | SessionID) + (JointState - 1 | SessionID)']};

% models = {['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
%         '+ (V1_logodds_PRE | AnimalID)' ...
%         '+ (V1_logodds_PRE| SessionID)'],...
%         ['HPC_logodds ~ V1_logodds * JointState ' ...
%         '+ (V1_logodds  | AnimalID) ' ...
%         '+ (V1_logodds | SessionID)']};
tbl = readtable('V1_HPC_reactivation_coherence_lme_HPC1.csv');

tbl.HPC_trackID = nan(size(tbl.AnimalID));
tbl.V1_trackID_PRE = nan(size(tbl.AnimalID));
for nsession = 1:max(tbl.SessionID)
    idx = find(tbl.SessionID == nsession);
    HPC_thresholds = prctile(tbl.HPC_logodds(idx),[25 75]);
    tbl.HPC_trackID(idx(tbl.HPC_logodds(idx) > HPC_thresholds(2))) = 1;
    tbl.HPC_trackID(idx(tbl.HPC_logodds(idx) < HPC_thresholds(1))) = 2;

    V1_thresholds = prctile(tbl.V1_logodds_PRE(idx),[25 75]);
    tbl.V1_trackID_PRE(idx(tbl.V1_logodds_PRE(idx) > V1_thresholds(2))) = 1;
    tbl.V1_trackID_PRE(idx(tbl.V1_logodds_PRE(idx) < V1_thresholds(1))) = 2;
end

tbl.AnimalID = categorical(tbl.AnimalID);
tbl.SessionID = categorical(tbl.SessionID);
tbl.HPC_trackID = categorical(tbl.HPC_trackID);
tbl.V1_trackID_PRE = categorical(tbl.V1_trackID_PRE);



% 
% formula = 'V1_logodds ~ V1_logodds_PRE + HPC_logodds + (1| SessionID) + (V1_logodds_PRE - 1| SessionID) + (HPC_logodds - 1| SessionID) + (1| AnimalID) + (V1_logodds_PRE - 1| AnimalID) + (HPC_logodds - 1| AnimalID)';
% formula = 'V1_logodds ~ V1_logodds_PRE + HPC_logodds +  (V1_logodds_PRE + HPC_logodds| SessionID) +  (V1_logodds_PRE + HPC_logodds| AnimalID)';
% 
% formula = 'V1_logodds ~ V1_logodds_PRE + HPC_trackID + (1| SessionID) + (V1_logodds_PRE - 1| SessionID) + (HPC_trackID - 1| SessionID) + (1| AnimalID) + (V1_logodds_PRE - 1| AnimalID) + (HPC_trackID - 1| AnimalID)';
% 
% formula = 'V1_logodds ~ V1_logodds_PRE + HPC_trackID +  (V1_logodds_PRE + HPC_trackID| SessionID) +  (V1_logodds_PRE + HPC_trackID| AnimalID)';
% 
% formula = 'V1_logodds ~ V1_trackID_PRE + HPC_trackID +  (V1_trackID_PRE + HPC_trackID| SessionID) +  (V1_trackID_PRE + HPC_trackID| AnimalID)';
formula = 'V1_logodds ~ V1_trackID_PRE + HPC_trackID + (1| SessionID) + (V1_trackID_PRE - 1| SessionID) + (HPC_trackID - 1| SessionID) + (1| AnimalID) + (V1_trackID_PRE - 1| AnimalID) + (HPC_trackID - 1| AnimalID)';

lme = fitlme(tbl, formula,'DummyVarCoding','effects');
%     lme = fitlme(tbl, formula)
lme.Coefficients
[Psi, ~, ~] = covarianceParameters(lme);
Psi_matrix = Psi{1}


models = {'V1_logodds ~ V1_logodds_PRE + HPC_logodds +  (V1_logodds_PRE + HPC_logodds| SessionID) +  (V1_logodds_PRE + HPC_logodds| AnimalID) ',
    'V1_logodds ~ V1_trackID_PRE + HPC_trackID +  (V1_trackID_PRE + HPC_trackID| SessionID) +  (V1_trackID_PRE + HPC_trackID| AnimalID)'};

for i = 1:length(models)
    formula = models{i};

    % formula = ['HPC_logodds ~ V1_logodds_PRE * JointState ' ...
    %            '+ (V1_logodds_PRE + JointState| SessionID) + (V1_logodds_PRE + JointState | AnimalID)'];

    lme = fitlme(tbl, formula,'DummyVarCoding','effects');
%     lme = fitlme(tbl, formula)
    lme.Coefficients

    PRE_POST_V1_log_odds_lme(i).model = lme;
    PRE_POST_V1_log_odds_lme(i).variable = [lme.Coefficients.Name];
    PRE_POST_V1_log_odds_lme(i).p = [lme.Coefficients.pValue(:)];
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    PRE_POST_V1_log_odds_lme(i).b = [lme.Coefficients.Estimate(:)];
    PRE_POST_V1_log_odds_lme(i).t = [lme.Coefficients.tStat(:)];
    PRE_POST_V1_log_odds_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    PRE_POST_V1_log_odds_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    PRE_POST_V1_log_odds_lme(i).marginal_R2 = marginal_R2;
end

% for i = 1:2
clear VariableIdx
% Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}
% VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
% nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
%     226, 132, 187;   % interpolated 2/3
%     212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

fig = figure('Name','PRE and HPC log odds predicting V1 log odds KDE reactivation coherence','Position',[482 111 665 851]);
sgtitle('SO phase state KDE reactivation coherence')

VariableName = {'PRE','HPC'};
titles = {'HPC log odds','HPC reactivation track'};
for i = 1:2
    %     if sum(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'))>0
    %         VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, ':JointState'));
    %     else
    VariableIdx = find(contains(PRE_POST_V1_log_odds_lme(i).variable, VariableName));
    %     end
    clear b_CI b_mean
    b_CI = [PRE_POST_V1_log_odds_lme(i).b_CI];
    b_mean = [PRE_POST_V1_log_odds_lme(i).b];

    % figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

    % fig = figure('Name','Spindle power KDE reactivation coherence PRE','Position',[482 111 665 851]);
    % sgtitle('Spindle power KDE reactivation coherence PRE')
    %

    lci_boot = [b_mean - b_CI(:,1)];
    uci_boot = [b_CI(:,2)-b_mean ];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    for n = 1:length(VariableIdx)
        bar(1+group_offset*(n-1),b_mean(VariableIdx(n)),bar_width, 'FaceColor', colour_lines(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot 95% CI Error Bar
        errorbar(1+group_offset*(n-1), b_mean(VariableIdx(n)), uci_boot(VariableIdx(n)), uci_boot(VariableIdx(n)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
        text(1+group_offset*(n-1),0.02+n*0.02,sprintf('%.3e',PRE_POST_V1_log_odds_lme(i).p(VariableIdx(n))))
    end

    ylabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.01 0.35])
    title(titles{i})
end
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','PRE_POST_V1_log_odds_lme.mat'),'PRE_POST_V1_log_odds_lme');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])








%% RRR 
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','RRR_reactivation_ripples_PSTH.mat'))
timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;

z_bias = RRR_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = RRR_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';

z_bias_RRR = z_bias(:,event_ids_first);
z_bias_V1_RRR = z_bias_V1(:,event_ids_first);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.2;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


subject_id = subject_id(event_ids_first);
session_count = session_count(event_ids_first);

singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));


%%% Ripple power
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



%%%%%% mixed effect model (ripple power)
bins_to_use = bin_centers>0 & bin_centers<0.1;

bins_to_select = bin_centers>-0.2& bin_centers<0;

bins_to_select = bin_centers>0 & bin_centers<0.2;

[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
mean_bias_V1 = mean(z_bias_V1_RRR(bins_to_select, :), 'omitnan');
% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias_RRR(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds = thresholds(1:end-1);

% spindle_power_log_odds

mean_beta = [];
ripple_power_modulation_lme = struct();

for i = 1:10

    th = thresholds(i);
    t1 = find(mean_bias >= th);
    t2 = find(mean_bias <= -th);

    selected_events = [t1 t2];
    %
    tbl = table([ripple_info.ripple_power(t1); ripple_info.ripple_power(t2)],...
        mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),'VariableNames',{'RipplePower','V1_logodds','HPC_logodds','HPC_trackID','AnimalID'});


    % lme = fitlme(tbl, 'V1_logodds ~ HPC_trackID + (1 | AnimalID)');
    lme = fitlme(tbl, 'V1_logodds ~ HPC_trackID * RipplePower + (1 | AnimalID)');
    % lme = fitlme(tbl, 'V1_logodds ~ HPC_logodds * RipplePower + (1 | AnimalID)');
    lme.Coefficients
    %
    ripple_power_modulation_lme(i).variable = lme.Coefficients.Name;
    ripple_power_modulation_lme(i).p = lme.Coefficients.pValue(:);
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    ripple_power_modulation_lme(i).b = lme.Coefficients.Estimate(:);
    ripple_power_modulation_lme(i).t = lme.Coefficients.tStat(:);
    ripple_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    ripple_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    ripple_power_modulation_lme(i).marginal_R2 = marginal_R2;
end

% log_odds_lme = ripple_power_modulation_lme;
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme_PRE.mat'),'log_odds_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme.mat'),'log_odds_lme');


save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_RRR.mat'),'ripple_power_modulation_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_PRE_RRR.mat'),'ripple_power_modulation_lme');



%%%%% Plot glme beta CI and p value
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_PRE_RRR.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10

    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end

fig = figure('Name','Ripple power RRR reactivation coherence','Position',[482 111 665 851]);
sgtitle('Ripple power RRR reactivation coherence')

subplot(2,2,1)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,4),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % bar(1+group_offset*2,b_mean(i,6), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ i*group_offset, b_mean(i,4), uci_boot(i,4), uci_boot(i,4), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+2*group_offset, b_mean(i,6), uci_boot(i,6), uci_boot(i,6), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(4)))
    % text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(6)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
title('-0.2 -0s')

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_RRR.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10

    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end


subplot(2,2,2)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,4),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % bar(1+group_offset*2,b_mean(i,6), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ i*group_offset, b_mean(i,4), uci_boot(i,4), uci_boot(i,4), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+2*group_offset, b_mean(i,6), uci_boot(i,6), uci_boot(i,6), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(4)))
    % text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(6)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
title('0 -0.2s')

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])





%%%%% Plot glme beta CI and p value
clear match_beta_CI non_match_beta_CI mean_beta non_match_mean_beta

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme_PRE.mat'),'log_odds_lme');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_PRE_RRR.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10

    % b_CI(i,:,:) = [log_odds_lme(i).b_CI];
    % b_mean(i,:) = [log_odds_lme(i).b];
    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
    % spindle_power_modulation_lme(i).non_match_b_shuffle
end

% figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

fig = figure('Name','log odds RRR reactivation coherence (ripple power model)','Position',[482 111 665 660]);
sgtitle('log odds RRR reactivation coherence')
subplot(2,2,1)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1+ i*group_offset,b_mean(i,2),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,3),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);


    % Plot 95% CI Error Bar
    errorbar(1+ i*group_offset, b_mean(i,3), uci_boot(i,3), uci_boot(i,3), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+ i*group_offset, b_mean(i,2), uci_boot(i,2), uci_boot(i,2), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    % text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(2)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(3)))


    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
ylim([0 2.7])
title('-0.2 -0s')

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme.mat'),'log_odds_lme');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_RRR.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10
    % b_CI(i,:,:) = [log_odds_lme(i).b_CI];
    % b_mean(i,:) = [log_odds_lme(i).b];
    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end


subplot(2,2,2)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % bar(1+ i*group_offset,b_mean(i,2),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,3),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+ i*group_offset, b_mean(i,2), uci_boot(i,2), uci_boot(i,2), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ i*group_offset, b_mean(i,3), uci_boot(i,3), uci_boot(i,3), 'k', 'LineWidth', 1.5, 'CapSize', 10);

    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(3)))
    % text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(6)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
title('0 -0.2s')
ylim([0 2.7])

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])



%% RRR spindle power effect

% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2



load(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'));
% load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));
periripple_LFP_info_V1 = periripple_LFP_info_V1_best_SO;

LFP_tvec = periripple_LFP_info_V1(1).tvec;
% nan_mask = interp1(bin_centers, KDE_reactivation_ripples_PSTH.nan_mask',LFP_tvec, 'previous');
% spindle_amplitude1 = [periripple_LFP_info_V1(1).spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{1}(:,ripples_all(2).SWS_index==1)] + nan_mask;
% spindle_amplitude2 = [periripple_LFP_info_V1(1).spindle_amplitude{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)]+ nan_mask;
% % 
spindle_amplitude1 = [periripple_LFP_info_V1(1).spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{1}(:,ripples_all(2).SWS_index==1)] ;
spindle_amplitude2 = [periripple_LFP_info_V1(1).spindle_amplitude{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)];
spindle_amplitude = nan([size(spindle_amplitude1),2]);
spindle_amplitude(:,:,1) = spindle_amplitude1;
spindle_amplitude(:,:,2) = spindle_amplitude2;

ripple_info.spindle_amplitude_temporal = spindle_amplitude;
ripple_info.spindle_amplitude_temporal = ripple_info.spindle_amplitude_temporal(:,event_ids_first,:);



nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



%%%%%% mixed effect model (spindle power)
bins_to_use = bin_centers>0 & bin_centers<0.1;

bins_to_select = bin_centers>-0.2& bin_centers<0;

bins_to_select = bin_centers>0 & bin_centers<0.2;

[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% mean_bias_V1 = mean(z_bias_V1_original(bins_to_select, :), 'omitnan');
mean_bias_V1 = mean(z_bias_V1_RRR(bins_to_select, :), 'omitnan');

% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias_RRR(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds = thresholds(1:end-1);

% spindle_power_log_odds

mean_beta = [];
spindle_power_modulation_lme = struct();

% s

for i = 1:10

    th = thresholds(i);
    t1 = find(mean_bias_V1 >= th);
    t2 = find(mean_bias_V1 <= -th);
    % %
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0);
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2) >0);

    % t1 = t1(ripple_info.spindle_amplitude(t1,2) >0 & ripple_info.spindle_amplitude(t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude(t2,1) > 0 & ripple_info.spindle_amplitude(t2,2) >0);
    selected_events = [t1 t2];
    %
    % tbl = table([ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],...
    %     [ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],...
    %     mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % tbl = table([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1)]',...
    %     [ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2)]',...
    %     mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});
    % 
    % tbl = table(normalize([ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],'zscore'),...
    %     normalize([ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],'zscore'),...
    %     normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});
    % % 
    tbl = table(normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1)],'zscore')',...
        normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2)],'zscore')',...
        normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});


    lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match + V1_trackID * SpindlePower_NonMatch + (1 | AnimalID)');
    % lme = fitlme(tbl, 'HPC_logodds ~V1_trackID*SpindlePower_Match + (1 | AnimalID)');
    % lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match + V1_trackID * SpindlePower_NonMatch + SpindlePower_NonMatch*SpindlePower_Match + (1 | AnimalID)');
    % lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match * SpindlePower_NonMatch + (1 | AnimalID)');

    spindle_power_modulation_lme(i).variable = lme.Coefficients.Name;
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    spindle_power_modulation_lme(i).b = lme.Coefficients.Estimate(:);
    spindle_power_modulation_lme(i).t = lme.Coefficients.tStat(:);
    spindle_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    spindle_power_modulation_lme(i).p = lme.Coefficients.pValue(:);
    spindle_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    spindle_power_modulation_lme(i).marginal_R2 = marginal_R2;

    % 
    % 
    % tbl = table(normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1)],'zscore')',...
    %     normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2)],'zscore')',...
    %     normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});
    % 
    % R2_shuffled = [];
    % R2_shuffled = [];
    % for iBoot = 1:1000
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     % idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);
    %     idx = datasample(1:length(selected_events), length(selected_events), 'Replace', false);
    %     % tbl.SpindlePower_NonMatch = tbl.SpindlePower_NonMatch(idx);
    % 
    %     targetPredictor = 'SpindlePower_Match';
    %     shuffledTbl = tbl;
    %     shuffledTbl.SpindlePower_Match = tbl.SpindlePower_Match(idx);
    % 
    %     animals = unique(tbl.AnimalID);
    % 
    %     for i = 1:numel(animals)
    %         % Find indices for the current animal
    %         animalIdx = find(tbl.AnimalID == animals(i));
    % 
    %         % Use your specific sampling logic:
    %         % datasample without replacement on indices is equivalent to a shuffle
    %         shuffled_indices = datasample(s, animalIdx, length(animalIdx), 'Replace', false);
    % 
    %         % Assign the original values to the shuffled positions for this animal
    %         shuffledTbl.(targetPredictor)(animalIdx) = tbl.(targetPredictor)(shuffled_indices);
    %     end
    % 
    % 
    %     lme = fitlme(shuffledTbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match   + (1 | AnimalID)');
    %        % lme = fitlme(shuffledTbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match* SpindlePower_NonMatch   + (1 | AnimalID)');
    %     [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    %     % spindle_power_modulation_lme(i).marginal_R2 = marginal_R2;
    %     R2_shuffled(1,iBoot) = marginal_R2;
    % end
    % spindle_power_modulation_lme(i).R2_shuffled(1,:) = R2_shuffled;
    % 
    % spindle_power_modulation_lme(i).marginal_R2 -  R2_shuffled(1,iBoot)

end


% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_RRR.mat'),'spindle_power_modulation_lme');
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_PRE_RRR.mat'),'spindle_power_modulation_lme');

%%%%% Plot glme beta CI and p value
clear match_beta_CI non_match_beta_CI mean_beta non_match_mean_beta

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_PRE_RRR.mat'),'spindle_power_modulation_lme');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_RRR.mat'),'spindle_power_modulation_lme');

clear b_CI b_mean
for i = 1:10
    b_CI(i,:,:) = [spindle_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [spindle_power_modulation_lme(i).b];
end

clear VariableIdx
Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}

for i = 1:length(Variablenames)
    VariableIdx(i) = find(contains(spindle_power_modulation_lme(1).variable,Variablenames{i}));
end

% figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

% fig = figure('Name','Spindle power RRR reactivation coherence PRE','Position',[482 111 665 851]);
% sgtitle('Spindle power RRR reactivation coherence PRE')
% 
fig = figure('Name','Spindle power RRR reactivation coherence','Position',[482 111 665 851]);
sgtitle('Spindle power RRR reactivation coherence')
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ group_offset,b_mean(i,VariableIdx(1)),bar_width, 'FaceColor', colour_lines(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+group_offset*2,b_mean(i,VariableIdx(2)), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ group_offset, b_mean(i,VariableIdx(1)), uci_boot(i,VariableIdx(1)), uci_boot(i,VariableIdx(1)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+2*group_offset, b_mean(i,VariableIdx(2)), uci_boot(i,VariableIdx(2)), uci_boot(i,VariableIdx(2)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ group_offset,0.03,sprintf('%.3e',spindle_power_modulation_lme(i).p(VariableIdx(1))))
    text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(VariableIdx(2))))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    % ylim([-0.03 0.05])

    ylim([-0.03 0.07])
    title(i)
end


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])

    % spindle_power_modulation_lme(i).b_CI 






%%
%%% RRR SO effect

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Moving spindle power binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2


load(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'));
% load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));
periripple_LFP_info_V1 = periripple_LFP_info_V1_best_SO;

LFP_tvec = periripple_LFP_info_V1(1).tvec;
% nan_mask = interp1(bin_centers, KDE_reactivation_ripples_PSTH.nan_mask',LFP_tvec, 'previous');
% SO_phase1 = [periripple_LFP_info_V1(1).SO_phase{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{1}(:,ripples_all(2).SWS_index==1)]+nan_mask;
% SO_phase2 = [periripple_LFP_info_V1(1).SO_phase{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{2}(:,ripples_all(2).SWS_index==1)] + nan_mask;

SO_phase1 = [periripple_LFP_info_V1(1).SO_phase{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{1}(:,ripples_all(2).SWS_index==1)];
SO_phase2 = [periripple_LFP_info_V1(1).SO_phase{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{2}(:,ripples_all(2).SWS_index==1)];
SO_phase = nan([size(SO_phase1),2]);
SO_phase(:,:,1) = SO_phase1;
SO_phase(:,:,2) = SO_phase2;

ripple_info.SO_phase_temporal = SO_phase;
ripple_info.SO_phase_temporal = ripple_info.SO_phase_temporal(:,event_ids_first,:);


nBoot = 1000;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



%%%%%% mixed effect model (SO phase)
bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>-0.1 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.1;
[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

mean_bias_V1 = mean(z_bias_V1_RRR(bins_to_select, :), 'omitnan');

% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias_RRR(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds = thresholds(1:end-1);

% spindle_power_log_odds

mean_beta = [];
SO_phase_modulation_lme = struct();

for i = 1:10

    th = thresholds(i);
    t1 = find(mean_bias_V1 >= th);
    t2 = find(mean_bias_V1 <= -th);
    % %
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0);
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2) >0);

    % t1 = t1(ripple_info.spindle_amplitude(t1,2) >0 & ripple_info.spindle_amplitude(t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude(t2,1) > 0 & ripple_info.spindle_amplitude(t2,2) >0);
    selected_events = [t1 t2];

    tbl = table([ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
        [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
        normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % 
    % tbl = table([ripple_info.SO_phase_temporal(LFP_bin,t1,2) ripple_info.SO_phase_temporal(LFP_bin,t2,1)]',...
    %     [ripple_info.SO_phase_temporal(LFP_bin,t1,1) ripple_info.SO_phase_temporal(LFP_bin,t2,2)]',...
    %     normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % tbl = table([ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
    % [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
    % mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    % categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % 
    % tbl = table([ripple_info.SO_phase_temporal(LFP_bin,t1,2) ripple_info.SO_phase_temporal(LFP_bin,t2,1)]',...
    %     [ripple_info.SO_phase_temporal(LFP_bin,t1,1) ripple_info.SO_phase_temporal(LFP_bin,t2,2)]',...
    %     mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});



    tbl.Match_trough = tbl.SOPhase_Match < -pi/2 | tbl.SOPhase_Match > pi/2 ;
    tbl.NonMatch_trough = tbl.SOPhase_NonMatch < -pi/2 | tbl.SOPhase_NonMatch > pi/2 ;

    % 2. Combine into a categorical 'State' variable
    % tbl.JointState = categorical(2 * tbl.Match_High + tbl.NonMatch_High);

    tbl.JointState = (2 * tbl.Match_trough + tbl.NonMatch_trough);
    % tbl.JointState(tbl.JointState >2) = 2;
    tbl.JointState = categorical(tbl.JointState);
    % 0 is when both were peak, 1 is when only non matching was trough, 2 is
    % when only matching was trough and 3 is when both were trough.
    lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * JointState + (1 | AnimalID)');
    H = [0,0, 0, 0, 0, 0, 1 , -1];
    lme.Coefficients
    [pVal, FStat] = coefTest(lme, H);


    SO_phase_modulation_lme(i).variable = lme.Coefficients.Name;
    SO_phase_modulation_lme(i).p = lme.Coefficients.pValue(:);
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    SO_phase_modulation_lme(i).b = lme.Coefficients.Estimate(:);
    SO_phase_modulation_lme(i).t = lme.Coefficients.tStat(:);
    SO_phase_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];


    SO_phase_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    SO_phase_modulation_lme(i).marginal_R2 = marginal_R2;

    % 
    % VariableNames = tbl.Properties.VariableNames;
    % beta_boot=[];
    % beta_shuffled=[];
    % t_boot=[];
    % t_shuffled=[];
    % parfor iBoot = 1:1000
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     idx = datasample(s, 1:length(selected_events), length(selected_events), 'Replace',true);
    % 
    %     tbl_temp = tbl;
    %     for var = 1:length(VariableNames)
    %         tbl_temp.(VariableNames{var}) = tbl.(VariableNames{var})(idx);
    %     end
    % 
    % 
    %     lme = fitlme(tbl_temp, 'HPC_logodds ~ V1_trackID * JointState + (1 | AnimalID)');
    %     beta_boot(iBoot,:) = [lme.Coefficients.Estimate(:)];
    %     t_boot(iBoot,:) = [lme.Coefficients.tStat(:)];
    %     % t_boot(iBoot) = [lme.Coefficients.tStat(5)];
    %     % non_match_t_boot(iBoot) = [lme.Coefficients.tStat(6)];
    % 
    % 
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     idx = datasample(s, 1:length(selected_events), length(selected_events), 'Replace',false);
    % 
    %     tbl_temp = tbl;
    %     for var = 1:length(VariableNames)
    %         tbl_temp.(VariableNames{var}) = tbl.(VariableNames{var})(idx);
    %     end
    %     tbl_temp.HPC_logodds = tbl.HPC_logodds;
    % 
    % 
    %     lme = fitlme(tbl_temp, 'HPC_logodds ~ V1_trackID * JointState + (1 | AnimalID)');
    % 
    %     beta_shuffled(iBoot,:) = [lme.Coefficients.Estimate(:)];
    %     t_shuffled(iBoot,:) = [lme.Coefficients.tStat(:)];
    %     % t_shuffled(iBoot) = [lme.Coefficients.tStat(5)];
    %     % non_match_t_shuffled(iBoot) = [lme.Coefficients.tStat(6)];
    % end


% b_CI(i,:,:) = prctile(beta_boot,[2.5 97.5]);
end


prctile(beta_boot(:,7)-beta_boot(:,8),[1.5 99.5])


% prctile(beta_boot,[2.5 97.5])

% 
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_RRR.mat'),'SO_phase_modulation_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_PRE_RRR.mat'),'SO_phase_modulation_lme');

%%%%% Plot glme beta CI and p value
clear match_beta_CI non_match_beta_CI mean_beta non_match_mean_beta

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_PRE_RRR.mat'),'SO_phase_modulation_lme');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_RRR.mat'),'SO_phase_modulation_lme');

clear b_CI b_mean
for i = 1:10
    % match_beta_CI(1,i,:) = prctile(spindle_power_modulation_lme(i).b_boot,[2.5 97.5]);
    % non_match_beta_CI(1,i,:) = prctile(spindle_power_modulation_lme(i).non_match_b_boot,[2.5 97.5]);
    % 
    % match_beta_CI(2,i,:) = prctile(spindle_power_modulation_lme(i).b_shuffle,[2.5 97.5]);
    % non_match_beta_CI(2,i,:) = prctile(spindle_power_modulation_lme(i).non_match_b_shuffle,[2.5 97.5]);
    b_CI(i,:,:) = [SO_phase_modulation_lme(i).b_CI];
    b_mean(i,:) = [SO_phase_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end

% figure;plot(squeeze(non_match_beta_CI(1,:,:)),'r');hold on; plot(squeeze(non_match_beta_CI(2,:,:)),'k');;yline(0,'k')
% 
% figure;plot(squeeze(match_beta_CI(1,:,:)),'r');hold on; plot(squeeze(match_beta_CI(2,:,:)),'k');yline(0,'k')


% figure;plot(squeeze(b_CI(:,5,:)),'r');hold on; plot(squeeze(b_CI(:,6,:)),'k')

fig = figure('Name','SO phase RRR reactivation coherence PRE','Position',[482 111 665 851]);
sgtitle('SO phase RRR reactivation coherence PRE')

% fig = figure('Name','SO phase RRR reactivation coherence','Position',[482 111 665 851]);
% sgtitle('SO phase RRR reactivation coherence')
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    bar(1,b_mean(i,6),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ group_offset,b_mean(i,8),bar_width, 'FaceColor', colour_lines(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+group_offset*2,b_mean(i,7), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    errorbar(1, b_mean(i,6), uci_boot(i,6), uci_boot(i,6), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ group_offset, b_mean(i,8), uci_boot(i,8), uci_boot(i,8), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+2*group_offset, b_mean(i,7), uci_boot(i,7), uci_boot(i,7), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    text(1,0.02,sprintf('%.3e',SO_phase_modulation_lme(i).p(6)))
    text(1+ group_offset,0.12,sprintf('%.3e',SO_phase_modulation_lme(i).p(8)))
    text(1+ 2*group_offset,0.22,sprintf('%.3e',SO_phase_modulation_lme(i).p(7)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.15 0.28])
    title(i)
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])




%% Bayesian
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
timebin = 0.02;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;

nan_mask = KDE_reactivation_ripples_PSTH.nan_mask(:,1:2:end)'; %KDE nan mask 10ms bin, make it 20ms bin

z_bias = [bayesian_reactivation_all(1).z_log_odds(:,ripples_all(1).SWS_index==1) bayesian_reactivation_all(2).z_log_odds(:,ripples_all(2).SWS_index==1)];
z_bias_V1 = [bayesian_reactivation_V1_all(1).z_log_odds(:,ripples_all(1).SWS_index==1) bayesian_reactivation_V1_all(2).z_log_odds(:,ripples_all(2).SWS_index==1)];

z_bias = z_bias + nan_mask;
z_bias_V1 = z_bias_V1 + nan_mask;

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';

z_bias_bayesian = z_bias(:,event_ids_first);
z_bias_V1_bayesian = z_bias_V1(:,event_ids_first);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.2;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


subject_id = subject_id(event_ids_first);
session_count = session_count(event_ids_first);

singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));

%%
%%% Ripple power
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



%%%%%% mixed effect model (ripple power)
bins_to_use = bin_centers>0 & bin_centers<0.1;

bins_to_select = bin_centers>-0.2& bin_centers<0;

bins_to_select = bin_centers>0 & bin_centers<0.2;

[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
mean_bias_V1 = mean(z_bias_V1_bayesian(bins_to_select, :), 'omitnan');
% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias_bayesian(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds = thresholds(1:end-1);

% spindle_power_log_odds

mean_beta = [];
ripple_power_modulation_lme = struct();

for i = 1:10

    th = thresholds(i);
    t1 = find(mean_bias >= th);
    t2 = find(mean_bias <= -th);

    selected_events = [t1 t2];
    %
    tbl = table([ripple_info.ripple_power(t1); ripple_info.ripple_power(t2)],...
        mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),'VariableNames',{'RipplePower','V1_logodds','HPC_logodds','HPC_trackID','AnimalID'});


    % lme = fitlme(tbl, 'V1_logodds ~ HPC_trackID + (1 | AnimalID)');
    lme = fitlme(tbl, 'V1_logodds ~ HPC_trackID * RipplePower + (1 | AnimalID)');
    % lme = fitlme(tbl, 'V1_logodds ~ HPC_logodds * RipplePower + (1 | AnimalID)');
    lme.Coefficients
    %
    ripple_power_modulation_lme(i).variable = lme.Coefficients.Name;
    ripple_power_modulation_lme(i).p = lme.Coefficients.pValue(:);
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    ripple_power_modulation_lme(i).b = lme.Coefficients.Estimate(:);
    ripple_power_modulation_lme(i).t = lme.Coefficients.tStat(:);
    ripple_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    ripple_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    ripple_power_modulation_lme(i).marginal_R2 = marginal_R2;
end

% log_odds_lme = ripple_power_modulation_lme;
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme_PRE.mat'),'log_odds_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme.mat'),'log_odds_lme');


save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_bayesian.mat'),'ripple_power_modulation_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_PRE_bayesian.mat'),'ripple_power_modulation_lme');



%%%%% Plot glme beta CI and p value
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_PRE_bayesian.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10

    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end

fig = figure('Name','Ripple power Bayesian reactivation coherence','Position',[482 111 665 851]);
sgtitle('Ripple power Bayesian reactivation coherence')

subplot(2,2,1)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,4),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % bar(1+group_offset*2,b_mean(i,6), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ i*group_offset, b_mean(i,4), uci_boot(i,4), uci_boot(i,4), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+2*group_offset, b_mean(i,6), uci_boot(i,6), uci_boot(i,6), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(4)))
    % text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(6)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
title('-0.2 -0s')

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_bayesian.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10

    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end


subplot(2,2,2)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,4),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % bar(1+group_offset*2,b_mean(i,6), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ i*group_offset, b_mean(i,4), uci_boot(i,4), uci_boot(i,4), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+2*group_offset, b_mean(i,6), uci_boot(i,6), uci_boot(i,6), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(4)))
    % text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(6)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
title('0 -0.2s')

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])





%%%%% Plot glme beta CI and p value
clear match_beta_CI non_match_beta_CI mean_beta non_match_mean_beta

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme_PRE.mat'),'log_odds_lme');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_PRE_bayesian.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10

    % b_CI(i,:,:) = [log_odds_lme(i).b_CI];
    % b_mean(i,:) = [log_odds_lme(i).b];
    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
    % spindle_power_modulation_lme(i).non_match_b_shuffle
end

% figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

fig = figure('Name','log odds Bayesian reactivation coherence (ripple power model)','Position',[482 111 665 660]);
sgtitle('log odds Bayesian reactivation coherence')
subplot(2,2,1)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1+ i*group_offset,b_mean(i,2),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,3),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);


    % Plot 95% CI Error Bar
    errorbar(1+ i*group_offset, b_mean(i,3), uci_boot(i,3), uci_boot(i,3), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+ i*group_offset, b_mean(i,2), uci_boot(i,2), uci_boot(i,2), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    % text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(2)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(3)))


    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
ylim([0 2.8])
title('-0.2 -0s')

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','log_odds_lme.mat'),'log_odds_lme');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','ripple_power_modulation_lme_bayesian.mat'),'ripple_power_modulation_lme');

clear b_CI b_mean
for i = 1:10
    % b_CI(i,:,:) = [log_odds_lme(i).b_CI];
    % b_mean(i,:) = [log_odds_lme(i).b];
    b_CI(i,:,:) = [ripple_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [ripple_power_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end


subplot(2,2,2)
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    % subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % bar(1+ i*group_offset,b_mean(i,2),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ i*group_offset,b_mean(i,3),bar_width, 'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % errorbar(1+ i*group_offset, b_mean(i,2), uci_boot(i,2), uci_boot(i,2), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ i*group_offset, b_mean(i,3), uci_boot(i,3), uci_boot(i,3), 'k', 'LineWidth', 1.5, 'CapSize', 10);

    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ i*group_offset,0.03,sprintf('%.3e',ripple_power_modulation_lme(i).p(3)))
    % text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(6)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    % xlim([0.5 2.3])
    % ylim([-0.03 0.05])
    % title(i)
end
title('0 -0.2s')
ylim([0 2.8])

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])


%%
%%% Bayesian spindle power effect

% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2



load(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'));
% load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));
periripple_LFP_info_V1 = periripple_LFP_info_V1_best_SO;

LFP_tvec = periripple_LFP_info_V1(1).tvec;
% nan_mask = interp1(bin_centers, KDE_reactivation_ripples_PSTH.nan_mask',LFP_tvec, 'previous');
% spindle_amplitude1 = [periripple_LFP_info_V1(1).spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{1}(:,ripples_all(2).SWS_index==1)] + nan_mask;
% spindle_amplitude2 = [periripple_LFP_info_V1(1).spindle_amplitude{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)]+ nan_mask;
% % 
spindle_amplitude1 = [periripple_LFP_info_V1(1).spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{1}(:,ripples_all(2).SWS_index==1)] ;
spindle_amplitude2 = [periripple_LFP_info_V1(1).spindle_amplitude{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)];
spindle_amplitude = nan([size(spindle_amplitude1),2]);
spindle_amplitude(:,:,1) = spindle_amplitude1;
spindle_amplitude(:,:,2) = spindle_amplitude2;

ripple_info.spindle_amplitude_temporal = spindle_amplitude;
ripple_info.spindle_amplitude_temporal = ripple_info.spindle_amplitude_temporal(:,event_ids_first,:);



nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



%%%%%% mixed effect model (spindle power)
bins_to_use = bin_centers>0 & bin_centers<0.1;

bins_to_select = bin_centers>-0.2& bin_centers<0;

bins_to_select = bin_centers>0 & bin_centers<0.2;

[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% mean_bias_V1 = mean(z_bias_V1_original(bins_to_select, :), 'omitnan');
mean_bias_V1 = mean(z_bias_V1_bayesian(bins_to_select, :), 'omitnan');

% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias_bayesian(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds = thresholds(1:end-1);

% spindle_power_log_odds

mean_beta = [];
spindle_power_modulation_lme = struct();

% s

for i = 1:10

    th = thresholds(i);
    t1 = find(mean_bias_V1 >= th);
    t2 = find(mean_bias_V1 <= -th);
    % %
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0);
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2) >0);

    % t1 = t1(ripple_info.spindle_amplitude(t1,2) >0 & ripple_info.spindle_amplitude(t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude(t2,1) > 0 & ripple_info.spindle_amplitude(t2,2) >0);
    selected_events = [t1 t2];
    %
    % tbl = table([ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],...
    %     [ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],...
    %     mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % tbl = table([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1)]',...
    %     [ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2)]',...
    %     mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});
    % 
    % tbl = table(normalize([ripple_info.spindle_amplitude(t1,2); ripple_info.spindle_amplitude(t2,1)],'zscore'),...
    %     normalize([ripple_info.spindle_amplitude(t1,1); ripple_info.spindle_amplitude(t2,2)],'zscore'),...
    %     normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});
    % % 
    tbl = table(normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1)],'zscore')',...
        normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2)],'zscore')',...
        normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});


    lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match + V1_trackID * SpindlePower_NonMatch + (1 | AnimalID)');
    % lme = fitlme(tbl, 'HPC_logodds ~V1_trackID*SpindlePower_Match + (1 | AnimalID)');
    % lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match + V1_trackID * SpindlePower_NonMatch + SpindlePower_NonMatch*SpindlePower_Match + (1 | AnimalID)');
    % lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match * SpindlePower_NonMatch + (1 | AnimalID)');

    spindle_power_modulation_lme(i).variable = lme.Coefficients.Name;
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    spindle_power_modulation_lme(i).b = lme.Coefficients.Estimate(:);
    spindle_power_modulation_lme(i).t = lme.Coefficients.tStat(:);
    spindle_power_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    spindle_power_modulation_lme(i).p = lme.Coefficients.pValue(:);
    spindle_power_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    spindle_power_modulation_lme(i).marginal_R2 = marginal_R2;

    % 
    % 
    % tbl = table(normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1)],'zscore')',...
    %     normalize([ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2)],'zscore')',...
    %     normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SpindlePower_Match','SpindlePower_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});
    % 
    % R2_shuffled = [];
    % R2_shuffled = [];
    % for iBoot = 1:1000
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     % idx = datasample(1:length(subject_id), length(subject_id), 'Replace', true);
    %     idx = datasample(1:length(selected_events), length(selected_events), 'Replace', false);
    %     % tbl.SpindlePower_NonMatch = tbl.SpindlePower_NonMatch(idx);
    % 
    %     targetPredictor = 'SpindlePower_Match';
    %     shuffledTbl = tbl;
    %     shuffledTbl.SpindlePower_Match = tbl.SpindlePower_Match(idx);
    % 
    %     animals = unique(tbl.AnimalID);
    % 
    %     for i = 1:numel(animals)
    %         % Find indices for the current animal
    %         animalIdx = find(tbl.AnimalID == animals(i));
    % 
    %         % Use your specific sampling logic:
    %         % datasample without replacement on indices is equivalent to a shuffle
    %         shuffled_indices = datasample(s, animalIdx, length(animalIdx), 'Replace', false);
    % 
    %         % Assign the original values to the shuffled positions for this animal
    %         shuffledTbl.(targetPredictor)(animalIdx) = tbl.(targetPredictor)(shuffled_indices);
    %     end
    % 
    % 
    %     lme = fitlme(shuffledTbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match   + (1 | AnimalID)');
    %        % lme = fitlme(shuffledTbl, 'HPC_logodds ~ V1_trackID * SpindlePower_Match* SpindlePower_NonMatch   + (1 | AnimalID)');
    %     [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    %     % spindle_power_modulation_lme(i).marginal_R2 = marginal_R2;
    %     R2_shuffled(1,iBoot) = marginal_R2;
    % end
    % spindle_power_modulation_lme(i).R2_shuffled(1,:) = R2_shuffled;
    % 
    % spindle_power_modulation_lme(i).marginal_R2 -  R2_shuffled(1,iBoot)

end


% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_bayesian.mat'),'spindle_power_modulation_lme');
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_PRE_bayesian.mat'),'spindle_power_modulation_lme');

%%%%% Plot glme beta CI and p value
clear match_beta_CI non_match_beta_CI mean_beta non_match_mean_beta

% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_PRE_bayesian.mat'),'spindle_power_modulation_lme');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','spindle_power_modulation_lme_bayesian.mat'),'spindle_power_modulation_lme');

clear b_CI b_mean
for i = 1:10
    b_CI(i,:,:) = [spindle_power_modulation_lme(i).b_CI];
    b_mean(i,:) = [spindle_power_modulation_lme(i).b];
end

clear VariableIdx
Variablenames = {'SpindlePower_NonMatch:V1_trackID_2','SpindlePower_Match:V1_trackID_2'}

for i = 1:length(Variablenames)
    VariableIdx(i) = find(contains(spindle_power_modulation_lme(1).variable,Variablenames{i}));
end

% figure;plot(squeeze(b_CI(:,6,:)),'r');hold on; plot(squeeze(b_CI(:,7,:)),'k');yline(0,'k')

fig = figure('Name','Spindle power Bayesian reactivation coherence PRE','Position',[482 111 665 851]);
sgtitle('Spindle power Bayesian reactivation coherence PRE')
% 
% fig = figure('Name','Spindle power Bayesian reactivation coherence','Position',[482 111 665 851]);
% sgtitle('Spindle power Bayesian reactivation coherence')
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    % bar(1,b_mean(i,5),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ group_offset,b_mean(i,VariableIdx(1)),bar_width, 'FaceColor', colour_lines(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+group_offset*2,b_mean(i,VariableIdx(2)), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    % errorbar(1, b_mean(i,5), uci_boot(i,5), uci_boot(i,5), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ group_offset, b_mean(i,VariableIdx(1)), uci_boot(i,VariableIdx(1)), uci_boot(i,VariableIdx(1)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+2*group_offset, b_mean(i,VariableIdx(2)), uci_boot(i,VariableIdx(2)), uci_boot(i,VariableIdx(2)), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    % text(1,0.02,sprintf('%.3e',spindle_power_modulation_lme(i).p(5)))
    text(1+ group_offset,0.03,sprintf('%.3e',spindle_power_modulation_lme(i).p(VariableIdx(1))))
    text(1+ 2*group_offset,0.04,sprintf('%.3e',spindle_power_modulation_lme(i).p(VariableIdx(2))))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    % ylim([-0.03 0.05])

    ylim([-0.03 0.07])
    title(i)
end


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])

    % spindle_power_modulation_lme(i).b_CI 






%%
%%% Bayesian SO effect

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Moving spindle power binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2


load(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'));
% load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));
periripple_LFP_info_V1 = periripple_LFP_info_V1_best_SO;

LFP_tvec = periripple_LFP_info_V1(1).tvec;
% nan_mask = interp1(bin_centers, KDE_reactivation_ripples_PSTH.nan_mask',LFP_tvec, 'previous');
% SO_phase1 = [periripple_LFP_info_V1(1).SO_phase{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{1}(:,ripples_all(2).SWS_index==1)]+nan_mask;
% SO_phase2 = [periripple_LFP_info_V1(1).SO_phase{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{2}(:,ripples_all(2).SWS_index==1)] + nan_mask;

SO_phase1 = [periripple_LFP_info_V1(1).SO_phase{1}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{1}(:,ripples_all(2).SWS_index==1)];
SO_phase2 = [periripple_LFP_info_V1(1).SO_phase{2}(:,ripples_all(1).SWS_index==1) periripple_LFP_info_V1(2).SO_phase{2}(:,ripples_all(2).SWS_index==1)];
SO_phase = nan([size(SO_phase1),2]);
SO_phase(:,:,1) = SO_phase1;
SO_phase(:,:,2) = SO_phase2;

ripple_info.SO_phase_temporal = SO_phase;
ripple_info.SO_phase_temporal = ripple_info.SO_phase_temporal(:,event_ids_first,:);


nBoot = 1000;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



%%%%%% mixed effect model (SO phase)
bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>-0.1 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.1;
[~,LFP_bin] = min(abs(LFP_tvec - mean(bin_centers(bins_to_select))));
% for npower = 1:nBins

% bins_to_select = bin_centers>-0.2 & bin_centers<0;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

mean_bias_V1 = mean(z_bias_V1_bayesian(bins_to_select, :), 'omitnan');

% mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
mean_bias = mean(z_bias_bayesian(bins_to_use, :), 'omitnan');
% selected_events = length(mean_bias);

thresholds = prctile(abs(mean_bias_V1), 0:10:100);
thresholds = thresholds(1:end-1);

% spindle_power_log_odds

mean_beta = [];
SO_phase_modulation_lme = struct();

for i = 1:10

    th = thresholds(i);
    t1 = find(mean_bias_V1 >= th);
    t2 = find(mean_bias_V1 <= -th);
    % %
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0);
    % t1 = t1(ripple_info.spindle_amplitude_temporal(LFP_bin,t1,2) >0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude_temporal(LFP_bin,t2,1) > 0 & ripple_info.spindle_amplitude_temporal(LFP_bin,t2,2) >0);

    % t1 = t1(ripple_info.spindle_amplitude(t1,2) >0 & ripple_info.spindle_amplitude(t1,1) >0);
    % t2 = t2(ripple_info.spindle_amplitude(t2,1) > 0 & ripple_info.spindle_amplitude(t2,2) >0);
    selected_events = [t1 t2];

    tbl = table([ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
        [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
        normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
        categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % 
    % tbl = table([ripple_info.SO_phase_temporal(LFP_bin,t1,2) ripple_info.SO_phase_temporal(LFP_bin,t2,1)]',...
    %     [ripple_info.SO_phase_temporal(LFP_bin,t1,1) ripple_info.SO_phase_temporal(LFP_bin,t2,2)]',...
    %     normalize(mean_bias(selected_events)','zscore'),normalize(mean_bias_V1(selected_events)','zscore'),categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % tbl = table([ripple_info.SO_phase(t1,2); ripple_info.SO_phase(t2,1)],...
    % [ripple_info.SO_phase(t1,1); ripple_info.SO_phase(t2,2)],...
    % mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    % categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});

    % 
    % tbl = table([ripple_info.SO_phase_temporal(LFP_bin,t1,2) ripple_info.SO_phase_temporal(LFP_bin,t2,1)]',...
    %     [ripple_info.SO_phase_temporal(LFP_bin,t1,1) ripple_info.SO_phase_temporal(LFP_bin,t2,2)]',...
    %     mean_bias(selected_events)',mean_bias_V1(selected_events)',categorical([2*ones(length(t1),1); ones(length(t2),1)]),...
    %     categorical(subject_id(selected_events)),'VariableNames',{'SOPhase_Match','SOPhase_NonMatch','V1_logodds','HPC_logodds','V1_trackID','AnimalID'});



    tbl.Match_trough = tbl.SOPhase_Match < -pi/2 | tbl.SOPhase_Match > pi/2 ;
    tbl.NonMatch_trough = tbl.SOPhase_NonMatch < -pi/2 | tbl.SOPhase_NonMatch > pi/2 ;

    % 2. Combine into a categorical 'State' variable
    % tbl.JointState = categorical(2 * tbl.Match_High + tbl.NonMatch_High);

    tbl.JointState = (2 * tbl.Match_trough + tbl.NonMatch_trough);
    % tbl.JointState(tbl.JointState >2) = 2;
    tbl.JointState = categorical(tbl.JointState);
    % 0 is when both were peak, 1 is when only non matching was trough, 2 is
    % when only matching was trough and 3 is when both were trough.
    lme = fitlme(tbl, 'HPC_logodds ~ V1_trackID * JointState + (1 | AnimalID)');
    H = [0,0, 0, 0, 0, 0, 1 , -1];
    lme.Coefficients
    [pVal, FStat] = coefTest(lme, H);


    SO_phase_modulation_lme(i).variable = lme.Coefficients.Name;
    SO_phase_modulation_lme(i).p = lme.Coefficients.pValue(:);
    % spindle_power_modulation_lme(i).non_match_p = lme.Coefficients.pValue(:);
    SO_phase_modulation_lme(i).b = lme.Coefficients.Estimate(:);
    SO_phase_modulation_lme(i).t = lme.Coefficients.tStat(:);
    SO_phase_modulation_lme(i).b_CI = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];


    SO_phase_modulation_lme(i).R2 = lme.Rsquared.Adjusted;
    % spindle_power_modulation_lme(i).marginal_R2 = [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    [marginal_R2,~] = calculate_marginal_R2(tbl,lme);
    SO_phase_modulation_lme(i).marginal_R2 = marginal_R2;

    % 
    % VariableNames = tbl.Properties.VariableNames;
    % beta_boot=[];
    % beta_shuffled=[];
    % t_boot=[];
    % t_shuffled=[];
    % parfor iBoot = 1:1000
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     idx = datasample(s, 1:length(selected_events), length(selected_events), 'Replace',true);
    % 
    %     tbl_temp = tbl;
    %     for var = 1:length(VariableNames)
    %         tbl_temp.(VariableNames{var}) = tbl.(VariableNames{var})(idx);
    %     end
    % 
    % 
    %     lme = fitlme(tbl_temp, 'HPC_logodds ~ V1_trackID * JointState + (1 | AnimalID)');
    %     beta_boot(iBoot,:) = [lme.Coefficients.Estimate(:)];
    %     t_boot(iBoot,:) = [lme.Coefficients.tStat(:)];
    %     % t_boot(iBoot) = [lme.Coefficients.tStat(5)];
    %     % non_match_t_boot(iBoot) = [lme.Coefficients.tStat(6)];
    % 
    % 
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     idx = datasample(s, 1:length(selected_events), length(selected_events), 'Replace',false);
    % 
    %     tbl_temp = tbl;
    %     for var = 1:length(VariableNames)
    %         tbl_temp.(VariableNames{var}) = tbl.(VariableNames{var})(idx);
    %     end
    %     tbl_temp.HPC_logodds = tbl.HPC_logodds;
    % 
    % 
    %     lme = fitlme(tbl_temp, 'HPC_logodds ~ V1_trackID * JointState + (1 | AnimalID)');
    % 
    %     beta_shuffled(iBoot,:) = [lme.Coefficients.Estimate(:)];
    %     t_shuffled(iBoot,:) = [lme.Coefficients.tStat(:)];
    %     % t_shuffled(iBoot) = [lme.Coefficients.tStat(5)];
    %     % non_match_t_shuffled(iBoot) = [lme.Coefficients.tStat(6)];
    % end


% b_CI(i,:,:) = prctile(beta_boot,[2.5 97.5]);
end


prctile(beta_boot(:,7)-beta_boot(:,8),[1.5 99.5])


% prctile(beta_boot,[2.5 97.5])

% 
save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_bayesian.mat'),'SO_phase_modulation_lme');
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_PRE_bayesian.mat'),'SO_phase_modulation_lme');

%%%%% Plot glme beta CI and p value
clear match_beta_CI non_match_beta_CI mean_beta non_match_mean_beta

load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_PRE_bayesian.mat'),'SO_phase_modulation_lme');
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation\reactivation_coherence_KDE_glme','SO_phase_modulation_lme_bayesian.mat'),'SO_phase_modulation_lme');

clear b_CI b_mean
for i = 1:10
    % match_beta_CI(1,i,:) = prctile(spindle_power_modulation_lme(i).b_boot,[2.5 97.5]);
    % non_match_beta_CI(1,i,:) = prctile(spindle_power_modulation_lme(i).non_match_b_boot,[2.5 97.5]);
    % 
    % match_beta_CI(2,i,:) = prctile(spindle_power_modulation_lme(i).b_shuffle,[2.5 97.5]);
    % non_match_beta_CI(2,i,:) = prctile(spindle_power_modulation_lme(i).non_match_b_shuffle,[2.5 97.5]);
    b_CI(i,:,:) = [SO_phase_modulation_lme(i).b_CI];
    b_mean(i,:) = [SO_phase_modulation_lme(i).b];
% spindle_power_modulation_lme(i).non_match_b_shuffle
end

% figure;plot(squeeze(non_match_beta_CI(1,:,:)),'r');hold on; plot(squeeze(non_match_beta_CI(2,:,:)),'k');;yline(0,'k')
% 
% figure;plot(squeeze(match_beta_CI(1,:,:)),'r');hold on; plot(squeeze(match_beta_CI(2,:,:)),'k');yline(0,'k')


% figure;plot(squeeze(b_CI(:,5,:)),'r');hold on; plot(squeeze(b_CI(:,6,:)),'k')

fig = figure('Name','SO phase Bayesian reactivation coherence PRE','Position',[482 111 665 851]);
sgtitle('SO phase Bayesian reactivation coherence PRE')

% fig = figure('Name','SO phase Bayesian reactivation coherence','Position',[482 111 665 851]);
% sgtitle('SO phase Bayesian reactivation coherence')
for i = 1:10
    lci_boot = [b_mean - squeeze(b_CI(:,:,1))];
    uci_boot = [b_mean - squeeze(b_CI(:,:,1))];

    subplot(3,4,i)
    bar_width = 0.3;      % Width of the bars
    group_offset = 0.4;    % Distance from the center integer (half the gap between bars)
    % Plot Bar
    hold on;
    bar(1,b_mean(i,6),bar_width, 'FaceColor', colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+ group_offset,b_mean(i,8),bar_width, 'FaceColor', colour_lines(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    bar(1+group_offset*2,b_mean(i,7), bar_width,'FaceColor', colour_lines(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot 95% CI Error Bar
    errorbar(1, b_mean(i,6), uci_boot(i,6), uci_boot(i,6), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+ group_offset, b_mean(i,8), uci_boot(i,8), uci_boot(i,8), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    errorbar(1+2*group_offset, b_mean(i,7), uci_boot(i,7), uci_boot(i,7), 'k', 'LineWidth', 1.5, 'CapSize', 10);
    text(1,0.02,sprintf('%.3e',SO_phase_modulation_lme(i).p(6)))
    text(1+ group_offset,0.12,sprintf('%.3e',SO_phase_modulation_lme(i).p(8)))
    text(1+ 2*group_offset,0.22,sprintf('%.3e',SO_phase_modulation_lme(i).p(7)))

    xlabel('Standardised b')
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
    xlim([0.5 2.3])
    ylim([-0.15 0.27])
    title(i)
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','reactivation_coherence_KDE_glme'),[])

