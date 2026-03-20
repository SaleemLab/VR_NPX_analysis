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

load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);


load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));


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
[~,spindle_index,~,index] =RestrictInts(merged_event_info.ripples_ints,merged_event_info.spindles_ints);
ripple_info.spindle_presence = spindle_index;
ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.spindle_presence_hemi(find(spindle_index)) = merged_event_info.spindles_hemisphere_id(index);

ripple_info.spindle_presence = ripple_info.spindle_presence(event_ids_first);
ripple_info.spindle_presence_hemi = ripple_info.spindle_presence_hemi(event_ids_first);

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
        SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
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


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%UP_normalised_durationUP_normalised_durationUP_normalised_duration
%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripple_normalised_duration.mat'));

% add time for each sessions for later sorting
sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);
UP_ints=[];
DOWN_ints=[];
ripple_peaktimes=[];
ripple_ints=[];
SO_ints=[];

for nprobe = 1:2
    UP_ints{nprobe}=slow_waves_all(nprobe).UP_ints;
    DOWN_ints{nprobe}=slow_waves_all(nprobe).DOWN_ints;
    SO_ints{nprobe} = slow_waves_all(nprobe).DOWN_intervals;
    ripple_peaktimes{nprobe}=ripples_all(nprobe).peaktimes;
    ripple_ints{nprobe}=[ripples_all(nprobe).onset ripples_all(nprobe).offset];
    % spindle_peaktimes{nprobe}=spindles_all(nprobe).peaktimes(spindles_all(nprobe).SWS_index == 1);
    % spindle_ints{nprobe}=[spindles_all(nprobe).onset(spindles_all(nprobe).SWS_index == 1) spindles_all(nprobe).offset(spindles_all(nprobe).SWS_index == 1)];

    for nsession = 1:max(slow_waves_all(1).UP_session_count)
        index = find(slow_waves_all(nprobe).DOWN_intervals_session == sessions_to_process(nsession));
        SO_ints{nprobe}(index,:) = SO_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession));
        UP_ints{nprobe}(index,:) = UP_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession));
        DOWN_ints{nprobe}(index,:) = DOWN_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(ripples_all(nprobe).SWS_index == 1 &ripples_all(nprobe).session_count == sessions_to_process(nsession));
        ripple_ints{nprobe}(index,:) = ripple_ints{nprobe}(index,:) + nsession * 1000000;
        ripple_peaktimes{nprobe}(index,:) = ripple_peaktimes{nprobe}(index,:) + nsession * 1000000;

        % [C,ia,ib] = intersect(find(spindles_all(nprobe).session_count == sessions_to_process(nsession)),find(spindles_all(nprobe).SWS_index == 1));
        % spindle_ints{nprobe}(ib,:) = spindle_ints{nprobe}(ib,:) + nsession * 1000000;
    end
end


%%% early and late UP transition
ripple_info.early_UP_index=[];
ripple_info.late_UP_index = [];
% ripple_info.early_UP_index=[];
ripple_info.after_SO_index = [];
ripple_info.before_SO_index = [];
twin = 0.25;

for nprobe = 1:2

    % temp =  merged_event_info.DOWN_hemisphere_id==nprobe;
    [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[UP_ints{nprobe}(:,1) UP_ints{nprobe}(:,1)+twin]);
    % [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[mean(merged_event_info.DOWN_ints(temp,:),2) mean(merged_event_info.DOWN_ints(temp,:),2)+0.25]);
    ripple_info.early_UP_index(:,nprobe) = UP_index;
    % [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[mean(merged_event_info.DOWN_ints(temp,:),2)-0.25 mean(merged_event_info.DOWN_ints(temp,:),2)]);
    [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[DOWN_ints{nprobe}(:,1)-twin DOWN_ints{nprobe}(:,1)]);
    ripple_info.late_UP_index(:,nprobe) = UP_index;


    [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[SO_ints{nprobe}(:,1) SO_ints{nprobe}(:,1)+twin]);
    % [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[mean(merged_event_info.DOWN_ints(temp,:),2) mean(merged_event_info.DOWN_ints(temp,:),2)+0.25]);
    ripple_info.after_SO_index(:,nprobe) = UP_index;
    % [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[mean(merged_event_info.DOWN_ints(temp,:),2)-0.25 mean(merged_event_info.DOWN_ints(temp,:),2)]);
    [~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[SO_ints{nprobe}(:,1)-twin SO_ints{nprobe}(:,1)]);
    ripple_info.before_SO_index(:,nprobe) = UP_index;
end

ripple_info.early_UP_index = ripple_info.early_UP_index(event_ids_first,:);
ripple_info.late_UP_index = ripple_info.late_UP_index(event_ids_first,:);


ripple_info.after_SO_index = ripple_info.after_SO_index(event_ids_first,:);
ripple_info.before_SO_index = ripple_info.before_SO_index(event_ids_first,:);


% 
% figure;
% histogram(ripple_info.SO_phase(ripple_info.early_UP_index(:,1)==1,1),'Normalization','probability')
% hold on;histogram(ripple_info.SO_phase(ripple_info.late_UP_index(:,1)==1,1),'Normalization','probability')
% 
% 
% figure;
% histogram(ripple_info.SO_phase(ripple_info.early_UP_index(:,2)==1,2),'Normalization','probability')
% hold on;histogram(ripple_info.SO_phase(ripple_info.late_UP_index(:,2)==1,2),'Normalization','probability')
%%%% Normalised SO (from DOWN peak to DOWN peak)

%%% First half and second half SO
for nprobe = 1:2

    % event_info (based on DOWN)
    UP_info = [];
    if nprobe == 1
        % UP_info{1} = event_info(1).L_ripple_normalised_DOWN_duration;
        % UP_info{2} = event_info(2).L_ripple_normalised_DOWN_duration;
        UP_info{1} = SO_normalised_duration(1).L_ripple;
        UP_info{2} = SO_normalised_duration(2).L_ripple;
    else
        % UP_info{1} = event_info(1).R_ripple_normalised_DOWN_duration;
        % UP_info{2} = event_info(2).R_ripple_normalised_DOWN_duration;
        UP_info{1} = SO_normalised_duration(1).R_ripple;
        UP_info{2} = SO_normalised_duration(2).R_ripple;
    end


    ripple_index = find(ripples_all(nprobe).SWS_index==1);
    normalised_duration{nprobe}=nan(length(ripple_index),2);
    event_duration{nprobe}=nan(length(ripple_index),2);

    for nevent = 1:length(ripple_index)

        for hemi = 1:2
            this_index = find(UP_info{hemi}(:,1) == ripple_index(nevent));
            if ~isempty(this_index)
                normalised_duration{nprobe}(nevent,hemi) =UP_info{hemi}(this_index,3);
                event_duration{nprobe}(nevent,hemi) =UP_info{hemi}(this_index,4);
            end
        end
    end
end

ripple_info.normalised_SO_duration = [normalised_duration{1}; normalised_duration{2}];
ripple_info.normalised_SO_duration = ripple_info.normalised_SO_duration(event_ids_first,:);

ripple_info.SO_event_duration = [event_duration{1}; event_duration{2}];
ripple_info.SO_event_duration = ripple_info.SO_event_duration(event_ids_first,:);


%%%% UP
%%% First half and second half 
for nprobe = 1:2

    % event_info (based on UP DOWN)
    UP_info = [];
    if nprobe == 1
        % UP_info{1} = event_info(1).L_ripple_normalised_UP_duration;
        % UP_info{2} = event_info(2).L_ripple_normalised_UP_duration;
        UP_info{1} = UP_normalised_duration(1).L_ripple;
        UP_info{2} = UP_normalised_duration(2).L_ripple;
        
    else
        % UP_info{1} = event_info(1).R_ripple_normalised_UP_duration;
        % UP_info{2} = event_info(2).R_ripple_normalised_UP_duration;
        UP_info{1} = UP_normalised_duration(1).R_ripple;
        UP_info{2} = UP_normalised_duration(2).R_ripple;
    end


    ripple_index = find(ripples_all(nprobe).SWS_index==1);
    normalised_duration{nprobe}=nan(length(ripple_index),2);
    event_duration{nprobe}=nan(length(ripple_index),2);

    for nevent = 1:length(ripple_index)

        for hemi = 1:2
            this_index = find(UP_info{hemi}(:,1) == ripple_index(nevent));
            if ~isempty(this_index)
                normalised_duration{nprobe}(nevent,hemi) =UP_info{hemi}(this_index,3);
                event_duration{nprobe}(nevent,hemi) =UP_info{hemi}(this_index,4);

            end
        end
    end
end

ripple_info.normalised_UP_duration = [normalised_duration{1}; normalised_duration{2}];
ripple_info.normalised_UP_duration = ripple_info.normalised_UP_duration(event_ids_first,:);

ripple_info.UP_duration = [event_duration{1}; event_duration{2}];
ripple_info.UP_duration = ripple_info.UP_duration(event_ids_first,:);
ripple_info.UP_duration(ripple_info.SO_event_duration>10)=nan;

%%% First half and second half DOWN
for nprobe = 1:2

    % event_info (based on DOWN)
    UP_info = [];
    if nprobe == 1
        % UP_info{1} = event_info(1).L_ripple_normalised_DOWN_duration;
        % UP_info{2} = event_info(2).L_ripple_normalised_DOWN_duration;
        UP_info{1} = DOWN_normalised_duration(1).L_ripple;
        UP_info{2} = DOWN_normalised_duration(2).L_ripple;
    else
        % UP_info{1} = event_info(1).R_ripple_normalised_DOWN_duration;
        % UP_info{2} = event_info(2).R_ripple_normalised_DOWN_duration;
        UP_info{1} = DOWN_normalised_duration(1).R_ripple;
        UP_info{2} = DOWN_normalised_duration(2).R_ripple;
    end


    ripple_index = find(ripples_all(nprobe).SWS_index==1);
    normalised_duration{nprobe}=nan(length(ripple_index),2);
    event_duration{nprobe}=nan(length(ripple_index),2);

    for nevent = 1:length(ripple_index)

        for hemi = 1:2
            this_index = find(UP_info{hemi}(:,1) == ripple_index(nevent));
            if ~isempty(this_index)
                normalised_duration{nprobe}(nevent,hemi) =UP_info{hemi}(this_index,3);
                event_duration{nprobe}(nevent,hemi) =UP_info{hemi}(this_index,4);
            end
        end
    end
end

ripple_info.normalised_DOWN_duration = [normalised_duration{1}; normalised_duration{2}];
ripple_info.normalised_DOWN_duration = ripple_info.normalised_DOWN_duration(event_ids_first,:);

ripple_info.DOWN_duration = [event_duration{1}; event_duration{2}];
ripple_info.DOWN_duration = ripple_info.DOWN_duration(event_ids_first,:);
ripple_info.DOWN_duration(ripple_info.SO_event_duration>10)=nan;
ripple_info.SO_event_duration(ripple_info.SO_event_duration>10)=nan;


% figure
% temp1 = ripple_info.normalised_DOWN_duration(:,1)<0.5 | ripple_info.normalised_DOWN_duration(:,2)<0.5;
% temp2 = ripple_info.normalised_DOWN_duration(:,1)>0.5 | ripple_info.normalised_DOWN_duration(:,2)>0.5;
% hold on;
% histogram(ripple_info.SO_phase(temp1,1),'Normalization','probability')
% histogram(ripple_info.SO_phase(temp2,1),'Normalization','probability')
% 
% figure
% temp1 = ripple_info.normalised_UP_duration(:,1)<0.33 | ripple_info.normalised_UP_duration(:,2)<0.33;
% temp2 = ripple_info.normalised_UP_duration(:,1)>0.33 & ripple_info.normalised_UP_duration(:,1)<0.66|  ripple_info.normalised_UP_duration(:,2)>0.33 & ripple_info.normalised_UP_duration(:,2)<0.66;
% temp3 = ripple_info.normalised_UP_duration(:,1)>0.66 | ripple_info.normalised_UP_duration(:,2)>0.66;
% hold on;
% histogram(ripple_info.SO_phase(temp1,1),-pi:0.1:pi,'Normalization','probability')
% histogram(ripple_info.SO_phase(temp2,1),-pi:0.1:pi,'Normalization','probability')
% histogram(ripple_info.SO_phase(temp3,1),-pi:0.1:pi,'Normalization','probability')


%%%%%%
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));

% merged_event_info.ripples_hemisphere_id
% singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));



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
%%%%%%%%%%%%%%% KDE reactivation bias 
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'))

timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;


% z_bias1 = KDE_reactivation_content.HPC_z_ripples';
% T1_events = find(nanmean(z_bias1(1:10,:))>0.5);
% T2_events = find(nanmean(z_bias1(1:10,:))<-0.5);

% z_bias = KDE_reactivation_ripples_PSTH.HPC_z_ripples';
% z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_ripples';

z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

% z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples';
% z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

%%%% Only grab unique ripple events
z_bias = z_bias(:,event_ids_first);
z_bias_V1 = z_bias_V1(:,event_ids_first);
% event_ids_first
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.2;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));


%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples with long UP state and short UP state
% 
Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];
event_duration_threshold = [0 1.5 10];

% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:33:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
nBins = 2;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)


                t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)<event_duration_threshold(npower+1) & ripple_info.UP_duration(true_idx(idx),2)>event_duration_threshold(npower) ...
                    & ripple_info.normalised_UP_duration(true_idx(idx),2) > 0) > 1';
                t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)<event_duration_threshold(npower+1) & ripple_info.UP_duration(true_idx(idx),1)>event_duration_threshold(npower) ...
                    & ripple_info.normalised_UP_duration(true_idx(idx),1) > 0) > 1';

                % HPC bias difference between Track 1 and Track 2
                t1_HPC = boot_HPC(t1);
                t2_HPC = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                end

                % --- Shifted pairing (null) ---
                t1s = boot_V1_shift >= th;
                t2s = boot_V1_shift <= -th;

                t1s = t1s' + (ripple_info.UP_duration(true_idx,2)<event_duration_threshold(npower+1) & ripple_info.UP_duration(true_idx,2)>event_duration_threshold(npower) ...
                    & ripple_info.normalised_UP_duration(true_idx,2) > 0) > 1';
                t2s = t2s' + (ripple_info.UP_duration(true_idx,1)<event_duration_threshold(npower+1) & ripple_info.UP_duration(true_idx,1)>event_duration_threshold(npower) ...
                    & ripple_info.normalised_UP_duration(true_idx,1) > 0) > 1';

                t1_HPCs = boot_HPC(t1s);
                t2_HPCs = boot_HPC(t2s);
                if any(t1s) && any(t2s)
                    diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                end


            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_duration.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_duration.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by UP state duration (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

% spindle_thresholds = {'SO trough and power 0-25','SO trough and power 25-50','SO trough and power 50-75','SO trough and power 75-100'};
UP_thresholds = {'UP less than 1s','UP longer than 1s'};

for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(UP_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by short UP vs long UP (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' Short Vs long UP state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);


% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])





%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during UP state or not UP state (less than 1.5 second)
% 
Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

event_duration_threshold = 1;

prctile(Duration,50)

% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:33:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
nBins = 2;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)

                if npower == 1 % UP (1.5s) and UP (not 1.5s or less)
                    t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),2) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx(idx),1)<event_duration_threshold)) > 1';
                    t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),1) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx(idx),2)<event_duration_threshold)) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.UP_duration(true_idx,2)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,2) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx,1)<event_duration_threshold)) > 1';
                    t2s = t2s' + (ripple_info.UP_duration(true_idx,1)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,1) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx,2)<event_duration_threshold)) > 1';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                else
                    t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),2) > 0 & ...
                        ripple_info.normalised_UP_duration(true_idx(idx),1) > 0) > 1';
                    t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),1) > 0 & ...
                        ripple_info.normalised_UP_duration(true_idx(idx),2) > 0) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.UP_duration(true_idx,2)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,2) > 0 & ...
                        ripple_info.UP_duration(true_idx,1)<event_duration_threshold) > 1';
                    t2s = t2s' + (ripple_info.UP_duration(true_idx,1)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,1) > 0 & ...
                        ripple_info.UP_duration(true_idx,2)<event_duration_threshold) > 1';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                end

            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_short_1s.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_short.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by UP states (less than 1s) (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

% spindle_thresholds = {'SO trough and power 0-25','SO trough and power 25-50','SO trough and power 50-75','SO trough and power 75-100'};
UP_thresholds = {'Unilateral UP','Bilateral UP'};

for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(UP_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by Unilateral UP vs Bilateral UP (less than 1s) (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' 1st vs 2nd half UP state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);


% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])






%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during UP state or not UP state (longer than 1.5 second)
% 
Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

event_duration_threshold = prctile(Duration,50);

% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:33:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
nBins = 2;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)

                if npower == 1 % UP (1.5s) and UP (not 1.5s or less)
                    t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),2) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx(idx),1)>event_duration_threshold)) > 1';
                    t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),1) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx(idx),2)>event_duration_threshold)) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.UP_duration(true_idx,2)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,2) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx,1)>event_duration_threshold)) > 1';
                    t2s = t2s' + (ripple_info.UP_duration(true_idx,1)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,1) > 0 & ...
                        ~(ripple_info.UP_duration(true_idx,2)>event_duration_threshold)) > 1';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                else
                    t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),2) > 0 & ...
                        ripple_info.normalised_UP_duration(true_idx(idx),1) > 0) > 1';
                    t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),1) > 0 & ...
                        ripple_info.normalised_UP_duration(true_idx(idx),2) > 0) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.UP_duration(true_idx,2)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,2) > 0 & ...
                        ripple_info.UP_duration(true_idx,1)>event_duration_threshold) > 1';
                    t2s = t2s' + (ripple_info.UP_duration(true_idx,1)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,1) > 0 & ...
                        ripple_info.UP_duration(true_idx,2)>event_duration_threshold) > 1';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                end

            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_long.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_long.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by UP states (longer than 1.5s) (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

% spindle_thresholds = {'SO trough and power 0-25','SO trough and power 25-50','SO trough and power 50-75','SO trough and power 75-100'};
UP_thresholds = {'Unilateral UP','Bilateral UP'};

for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(UP_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by Unilateral UP vs Bilateral UP (longer than 1.5s) (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' Unilateral vs Bilateral UP');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);


% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during first half and second half of normalised UP (less than 1.5s)
% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:50:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

event_duration_threshold = prctile(Duration,50);

duration_thresholds = [0, 1/2 1;0 1/2, 1]';
nBins = length(duration_thresholds)-1 ;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)

                t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),2) > duration_thresholds(npower,2) & ...
                    ripple_info.normalised_UP_duration(true_idx(idx),2) <= duration_thresholds(npower+1,2)) > 1';
                t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),1) > duration_thresholds(npower,1) & ...
                    ripple_info.normalised_UP_duration(true_idx(idx),1) <= duration_thresholds(npower+1,1)) > 1';

                % HPC bias difference between Track 1 and Track 2
                t1_HPC = boot_HPC(t1);
                t2_HPC = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                end

                % --- Shifted pairing (null) ---
                t1s = boot_V1_shift >= th;
                t2s = boot_V1_shift <= -th;


                t1s = t1s' + (ripple_info.UP_duration(true_idx,2)<event_duration_threshold& ripple_info.normalised_UP_duration(true_idx,2) > duration_thresholds(npower,2) & ...
                    ripple_info.normalised_UP_duration(true_idx,2) <= duration_thresholds(npower+1,2)) > 1';
                t2s = t2s' + (ripple_info.UP_duration(true_idx,1)<event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,1) > duration_thresholds(npower,1) & ...
                    ripple_info.normalised_UP_duration(true_idx,1) <= duration_thresholds(npower+1,1)) > 1';
                t1_HPCs = boot_HPC(t1s);
                t2_HPCs = boot_HPC(t2s);
                if any(t1s) && any(t2s)
                    diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_UP_duration_short.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_UP_duration_short.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by normalised UP duration less than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');


for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(sprintf('Normalised UP duration %d (%.2f–%.2f)', npower, ...
        duration_thresholds(npower), duration_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC first half vs second half UP less than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' 1st vs 2nd half UP state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])




%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during first half and second half of normalised UP (more than 1.5s)
% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:50:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

event_duration_threshold = prctile(Duration,50);

duration_thresholds = [0, 1/2 1;0 1/2, 1]';
nBins = length(duration_thresholds)-1 ;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)

                t1 = t1' + (ripple_info.UP_duration(true_idx(idx),2)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),2) > duration_thresholds(npower,2) & ...
                    ripple_info.normalised_UP_duration(true_idx(idx),2) <= duration_thresholds(npower+1,2)) > 1';
                t2 = t2' + (ripple_info.UP_duration(true_idx(idx),1)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx(idx),1) > duration_thresholds(npower,1) & ...
                    ripple_info.normalised_UP_duration(true_idx(idx),1) <= duration_thresholds(npower+1,1)) > 1';

                % HPC bias difference between Track 1 and Track 2
                t1_HPC = boot_HPC(t1);
                t2_HPC = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                end

                % --- Shifted pairing (null) ---
                t1s = boot_V1_shift >= th;
                t2s = boot_V1_shift <= -th;


                t1s = t1s' + (ripple_info.UP_duration(true_idx,2)>event_duration_threshold& ripple_info.normalised_UP_duration(true_idx,2) > duration_thresholds(npower,2) & ...
                    ripple_info.normalised_UP_duration(true_idx,2) <= duration_thresholds(npower+1,2)) > 1';
                t2s = t2s' + (ripple_info.UP_duration(true_idx,1)>event_duration_threshold & ripple_info.normalised_UP_duration(true_idx,1) > duration_thresholds(npower,1) & ...
                    ripple_info.normalised_UP_duration(true_idx,1) <= duration_thresholds(npower+1,1)) > 1';
                t1_HPCs = boot_HPC(t1s);
                t2_HPCs = boot_HPC(t2s);
                if any(t1s) && any(t2s)
                    diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_UP_duration_long.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_UP_duration_long.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by normalised UP duration more than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');


for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(sprintf('Normalised SO duration %d (%.2f–%.2f)', npower, ...
        duration_thresholds(npower), duration_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC first half vs second half UP more than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' 1st vs 2nd half UP state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])



%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during normalised SO 'peak' vs 'trough' (less than 1.5s)

% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:50:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
% duration_thresholds = [0, 1/2 1;0 1/2, 1]';
nBins = 2 ;

Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

% event_duration_threshold = 2;
event_duration_threshold = prctile(Duration,50);

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)
                if npower == 1 % close to DOWN
                    t1 = t1' + (ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.75) > 2';
                    t2 = t2' + (ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.75) > 2';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx,2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,2) > 0.75) > 2';
                    t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx,1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,1) > 0.75) > 2';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                else % Away from DOWN
                    t1 = t1' +(ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.75) > 2';
                    t2 = t2' +(ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.75) > 2';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,2) < 0.75) > 2';
                    t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,1) < 0.75) > 2';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end


                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_2s.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_2s.mat'))
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_short.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_short.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by normalised UP duration peak vs trough less than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

duration_thresholds = {'Normalised SO peak','Normalised SO trough'}
for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(duration_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC normalised SO peak vs trough less than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),1);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),1,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title('Normalised SO peak vs trough');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])



%% Temporal log odds AUC of ripples that occur during unilateral vs bilateral normalised SO 'trough' (less than 1.5s)

nBins = 2 ;

Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

event_duration_threshold = 3;
% event_duration_threshold = prctile(Duration,50);

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)
                if npower == 1 % Unilateral trough
                    t1 = t1' + (ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.75) + (ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.75) + (ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold)> 4';
                    t2 = t2' + (ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.75) + (ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.75) +(ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold)> 4';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx,1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,1) > 0.75) + (ripple_info.normalised_SO_duration(true_idx,2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,2) < 0.75) + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) > 4';
                    t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx,2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,2) > 0.75) + (ripple_info.normalised_SO_duration(true_idx,1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,1) < 0.75) + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold)> 4';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                else % Bilateral SO trough
                    t1 = t1' +(ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.75) + (ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.75) + (ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold) > 4';
                    t2 = t2' +(ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.75) + (ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.75) + (ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold)>4';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,2) < 0.75) + (ripple_info.normalised_SO_duration(true_idx,1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,1) < 0.75) + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) > 4';
                    t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,1) < 0.75) + (ripple_info.normalised_SO_duration(true_idx,2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,2) < 0.75) + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold)> 4';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end


                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_2s.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_2s.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_synchrony.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_synchrony.mat'))
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_synchrony_3s.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_synchrony_3s.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC normalised SO trough synchony less than 3s (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

duration_thresholds = {'Normalised SO peak','Normalised SO trough'}
for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(duration_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC normalised SO unilateral vs bilateral trough less than 3s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),1);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),1,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title('Normalised SO peak vs trough');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);






% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])

%%



%% Temporal log odds AUC of ripples that occur during normalised SO 'peak' vs 'trough' (longer than 1.5s)

% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:50:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
% duration_thresholds = [0, 1/2 1;0 1/2, 1]';
nBins = 2 ;

Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];

event_duration_threshold = prctile(Duration,50);

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)
                if npower == 1 % close to DOWN
                    t1 = t1' + (ripple_info.SO_event_duration(true_idx(idx),2) > event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.75) > 2';
                    t2 = t2' + (ripple_info.SO_event_duration(true_idx(idx),1) > event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.75) > 2';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) > event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx,2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,2) > 0.75) > 2';
                    t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) > event_duration_threshold) + (ripple_info.normalised_SO_duration(true_idx,1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,1) > 0.75) > 2';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                else % Away from DOWN
                    t1 = t1' +(ripple_info.SO_event_duration(true_idx(idx),2) > event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.75) > 2';
                    t2 = t2' +(ripple_info.SO_event_duration(true_idx(idx),1) > event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.75) > 2';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) > event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,2) < 0.75) > 2';
                    t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) > event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,1) < 0.75) > 2';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end


                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_long.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_long.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by normalised UP duration peak vs trough longer than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

duration_thresholds = {'Normalised SO peak','Normalised SO trough'}
for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(duration_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC normalised SO peak vs trough longer than 1.5s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),1);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),1,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title('Normalised SO peak vs trough');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])















%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during first half and second half of SO state 
% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:50:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
duration_thresholds = [0, 1/2 1;0 1/2, 1]';
nBins = length(duration_thresholds)-1 ;


Duration = [ripple_info.UP_duration(:,1); ripple_info.UP_duration(:,2)];
Duration(Duration>10)=[];
event_duration_threshold = 2;
% event_duration_threshold = prctile(Duration,50);

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                t1 = t1' +(ripple_info.SO_event_duration(true_idx(idx),2) < event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),2) > duration_thresholds(npower) & ...
                    ripple_info.normalised_SO_duration(true_idx(idx),2) < duration_thresholds(npower+1)) > 2';
                t2 = t2' +(ripple_info.SO_event_duration(true_idx(idx),1) < event_duration_threshold)+ (ripple_info.normalised_SO_duration(true_idx(idx),1) > duration_thresholds(npower) & ...
                    ripple_info.normalised_SO_duration(true_idx(idx),1) < duration_thresholds(npower+1)) > 2';

                % HPC bias difference between Track 1 and Track 2
                t1_HPC = boot_HPC(t1);
                t2_HPC = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                end


                % --- Shifted pairing (null) ---
                t1s = boot_V1_shift >= th;
                t2s = boot_V1_shift <= -th;

                t1s = t1s' + (ripple_info.SO_event_duration(true_idx,2) < event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,2) > duration_thresholds(npower) & ...
                    ripple_info.normalised_SO_duration(true_idx,2) < duration_thresholds(npower+1)) > 2';
                t2s = t2s' + (ripple_info.SO_event_duration(true_idx,1) < event_duration_threshold) +(ripple_info.normalised_SO_duration(true_idx,1) > duration_thresholds(npower) & ...
                    ripple_info.normalised_SO_duration(true_idx,1) < duration_thresholds(npower+1)) > 2';

                t1_HPCs = boot_HPC(t1s);
                t2_HPCs = boot_HPC(t2s);
                if any(t1s) && any(t2s)
                    diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                end

            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_short.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_short.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by normalised SO duration less than 2s (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');


for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(sprintf('Normalised SO duration %d (%.2f–%.2f)', npower, ...
        duration_thresholds(npower), duration_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC first half vs second half SO less than 2s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' 1st vs 2nd half SO state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])




%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal log odds AUC of ripples that occur during normalised SO 'peak' vs 'trough'

% duration_thresholds = prctile(ripple_info.normalised_UP_duration, 0:50:100);
% duration_thresholds = [0, 1/3,2/3, 1;0, 1/3,2/3, 1]';
% duration_thresholds = [0, 1/2 1;0 1/2, 1]';
nBins = 2 ;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2

                % --- Normalised duration ---
                % Track 1 → use right probe (2), Track 2 → left probe (1)
                if npower == 1 % close to DOWN
                    t1 = t1' + (ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.75) > 1';
                    t2 = t2' + (ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.75) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.normalised_SO_duration(true_idx,2) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,2) > 0.75) > 1';
                    t2s = t2s' + (ripple_info.normalised_SO_duration(true_idx,1) < 0.25 | ...
                        ripple_info.normalised_SO_duration(true_idx,1) > 0.75) > 1';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end

                else % Away from DOWN
                    t1 = t1' + (ripple_info.normalised_SO_duration(true_idx(idx),2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),2) < 0.75) > 1';
                    t2 = t2' + (ripple_info.normalised_SO_duration(true_idx(idx),1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx(idx),1) < 0.75) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;

                    t1s = t1s' + (ripple_info.normalised_SO_duration(true_idx,2) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,2) < 0.75) > 1';
                    t2s = t2s' + (ripple_info.normalised_SO_duration(true_idx,1) > 0.25 & ...
                        ripple_info.normalised_SO_duration(true_idx,1) < 0.75) > 1';

                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end


                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end
% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_normalised_SO_duration_tertiles.mat'))

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by normalised SO duration peak vs trough (0.1s win 0.02s step)','Position',[640 100 400 900/4*2]);
tiledlayout(nBins,1,'TileSpacing','compact');

duration_thresholds = {'Normalised SO peak','Normalised SO trough'}
for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(duration_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC normalised SO peak vs trough (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 2]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),1);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),1,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title('Normalised SO peak vs trough');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])


%% Early vs late UP (0.25 second)

% ripple_info.early_UP_index 
% ripple_info.late_UP_index
% Time windows
nBins = 2;
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);
event_duration_threshold =10;

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2
                
                if npower == 1
                    % Early
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1 = t1' & ripple_info.early_UP_index(true_idx(idx),2)==1 & ripple_info.SO_event_duration(true_idx(idx),2)<event_duration_threshold;
                    t2 = t2' & ripple_info.early_UP_index(true_idx(idx),1)==1 & ripple_info.SO_event_duration(true_idx(idx),1)<event_duration_threshold;
                else
                    % Late
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1 = t1' & ripple_info.late_UP_index(true_idx(idx),2)==1 & ripple_info.SO_event_duration(true_idx(idx),2)<event_duration_threshold;
                    t2 = t2' & ripple_info.late_UP_index(true_idx(idx),1)==1 & ripple_info.SO_event_duration(true_idx(idx),1)<event_duration_threshold;
                end

                % HPC bias difference between Track 1 and Track 2
                t1_HPC = boot_HPC(t1);
                t2_HPC = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                end

                % --- Shifted pairing (null) ---
                t1s = boot_V1_shift >= th;
                t2s = boot_V1_shift <= -th;


                if npower == 1
                    % Early
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1s = t1s' & ripple_info.early_UP_index(true_idx,2)==1 & ripple_info.SO_event_duration(true_idx,2)<event_duration_threshold;
                    t2s = t2s' & ripple_info.early_UP_index(true_idx,1)==1 & ripple_info.SO_event_duration(true_idx,1)<event_duration_threshold;
                else
                    % Late
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1s = t1s' & ripple_info.late_UP_index(true_idx,2)==1& ripple_info.SO_event_duration(true_idx,2)<event_duration_threshold;
                    t2s = t2s' & ripple_info.late_UP_index(true_idx,1)==1& ripple_info.SO_event_duration(true_idx,1)<event_duration_threshold;
                end

                t1_HPCs = boot_HPC(t1s);
                t2_HPCs = boot_HPC(t2s);
                if any(t1s) && any(t2s)
                    diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end

% -------- save --------
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DU_UD_500ms.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DU_UD_500ms.mat'))
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DU_UD_250ms_10s.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DU_UD_250ms_10s.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_250ms.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_250ms.mat'))
% -------- Plot --------
% fig = figure('Name','Temporal HPC log-odds AUC by DOWN-UP transitions 500ms (0.1s win 0.02s step)','Position',[640 100 400 900/2]);
fig = figure('Name','Temporal HPC log-odds AUC by DOWN-UP transitions 250ms 10s (0.1s win 0.02s step)','Position',[640 100 400 900/2]);
% fig = figure('Name','Temporal HPC log-odds AUC by DOWN-UP transitions (0.1s win 0.02s step)','Position',[640 100 400 900/2]);
tiledlayout(nBins,1,'TileSpacing','compact');

duration_thresholds = {'DOWN-UP transition','UP-DOWN transition','all'}
for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title( duration_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end


% -------- Plot --------
% fig = figure('Name','Temporal HPC log-odds AUC by DU and UD transition 500ms (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
fig = figure('Name','Temporal HPC log-odds AUC by DU and UD transition 250ms 10s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
% fig = figure('Name','Temporal HPC log-odds AUC by DU and UD transition 250ms (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 nBins]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' 1st vs 2nd half UP state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])







%% After and before SO peak (0.25 second)

% ripple_info.early_UP_index 
% ripple_info.late_UP_index
% Time windows
nBins = 2;
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;
event_duration_threshold = 10;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;

    % Sliding V1 window (used for event selection)
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % All events included, spindle bin applied later conditional on track side
        event_index = true(1, length(z_bias));

        % Compute mean log-odds
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');       % measure
        total_events = length(mean_bias_V1);
        if total_events < 10, continue; end

        % Quantile thresholds on |V1 bias|
        thresholds = prctile(abs(mean_bias_V1), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);
            true_idx = find(event_index);

            boot_V1  = mean_bias_V1(idx);
            boot_HPC = mean_bias_HPC(idx);

            % “Shifted”: randomise pairing between V1 & HPC
            boot_V1_shift = mean_bias_V1;
            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                t1 = boot_V1 >= th;     % Track 1
                t2 = boot_V1 <= -th;    % Track 2
                
                if npower == 1
                    % Early
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1 = t1' & ripple_info.after_SO_index(true_idx(idx),2)==1 & ripple_info.SO_event_duration(true_idx(idx),2)<event_duration_threshold;
                    t2 = t2' & ripple_info.after_SO_index(true_idx(idx),1)==1 & ripple_info.SO_event_duration(true_idx(idx),1)<event_duration_threshold;
                else
                    % Late
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1 = t1' & ripple_info.before_SO_index(true_idx(idx),2)==1 & ripple_info.SO_event_duration(true_idx(idx),2)<event_duration_threshold;
                    t2 = t2' & ripple_info.before_SO_index(true_idx(idx),1)==1 & ripple_info.SO_event_duration(true_idx(idx),1)<event_duration_threshold;
                end

                % HPC bias difference between Track 1 and Track 2
                t1_HPC = boot_HPC(t1);
                t2_HPC = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                end

                % --- Shifted pairing (null) ---
                t1s = boot_V1_shift >= th;
                t2s = boot_V1_shift <= -th;


                if npower == 1
                    % Early
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1s = t1s' & ripple_info.after_SO_index(true_idx,2)==1 & ripple_info.SO_event_duration(true_idx,2)<event_duration_threshold;
                    t2s = t2s' & ripple_info.after_SO_index(true_idx,1)==1 & ripple_info.SO_event_duration(true_idx,1)<event_duration_threshold;
                else
                    % Late
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1s = t1s' & ripple_info.before_SO_index(true_idx,2)==1& ripple_info.SO_event_duration(true_idx,2)<event_duration_threshold;
                    t2s = t2s' & ripple_info.before_SO_index(true_idx,1)==1& ripple_info.SO_event_duration(true_idx,1)<event_duration_threshold;
                end

                t1_HPCs = boot_HPC(t1s);
                t2_HPCs = boot_HPC(t2s);
                if any(t1s) && any(t2s)
                    diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end

% -------- save --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_250ms_10s.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_250ms_10s.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_500ms.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_500ms.mat'))
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_250ms.mat'),'AUC')
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_DOWN_peak_250ms.mat'))
% -------- Plot --------
% fig = figure('Name','Temporal HPC log-odds AUC by DOWN peak (0.1s win 0.02s step)','Position',[640 100 400 900/2]);
% fig = figure('Name','Temporal HPC log-odds AUC by DOWN peak 500ms (0.1s win 0.02s step)','Position',[640 100 400 900/2]);
fig = figure('Name','Temporal HPC log-odds AUC by DOWN peak 250ms 10s (0.1s win 0.02s step)','Position',[640 100 400 900/2]);
tiledlayout(nBins,1,'TileSpacing','compact');

duration_thresholds = {'After SO peak','Before SO peak','all'}
for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title( duration_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end


% -------- Plot --------
% fig = figure('Name','Temporal HPC log-odds AUC by before vs after DOWN peak 500ms (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
% fig = figure('Name','Temporal HPC log-odds AUC by before vs after DOWN peak (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
fig = figure('Name','Temporal HPC log-odds AUC by before vs after DOWN peak 250ms 10s (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 nBins]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(' 1st vs 2nd half UP state');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])



%% Temporal log odds AUC with UP state + spindle power (V1→HPC)
%%%%%%%%%%%%%%%%%%%%%%%% spindle phase binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/4:99.9);


% spindle_thresholds1 = prctile(ripple_info.spindle_amplitude(ripple_info.spindle_amplitude(:,1)>0,1), 0:99.9/2:99.9);
% spindle_thresholds2 = prctile(ripple_info.spindle_amplitude(ripple_info.spindle_amplitude(:,2)>0,2),  0:99.9/2:99.9);
% spindle_thresholds = [spindle_thresholds1;spindle_thresholds2]';
nBins =size(spindle_thresholds,1)-1;


% bins_to_select = bin_centers>0 & bin_centers<0.2;
% bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

win_size = 0.1;
step_size = 0.02;

time_bins = -1:step_size:1;
nTime = numel(time_bins);
bins_to_use = bin_centers>0 & bin_centers<0.1;

AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

% colour_lines = [ ...
%     241, 182, 218;   % original end (lightest)
%     231, 41, 138    % original start (darkest)
%     ] / 256;
% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    % 226, 132, 187;
    % 212,  78, 156;
    231,  41, 138] / 256;

for t = 1:nTime
    t0 = time_bins(t);
    t1 = t0 + win_size;
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('SO phase: processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);
    tic
    for npower = 1:nBins
        event_index =1:length(z_bias);
        % event_index = find(power_index);


        mean_bias_V1   = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
        % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
        mean_bias_HPC = mean(z_bias(bins_to_use, event_index), 'omitnan');
        mean_bias = mean_bias_V1;

        selected_events = length(mean_bias);

        thresholds = prctile(abs(mean_bias), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = length(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shifted_boot = NaN(nBoot, nThresh);
        % event_phase = ripple_info.SO_phase';

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, selected_events, selected_events, 1);

            true_idx = find(event_index);
            boot_HPC = mean_bias_HPC(idx);

            bb_shift = mean_bias_V1;
            boot_V1 = mean_bias_V1(idx);

            diff_tmp = NaN(1, nThresh);
            % prop_tmp = NaN(1, nThresh);
            diff_tmp_shifted = NaN(1, nThresh);
            % prop_tmp_shifted = NaN(1, nThresh);
            % event_phase = ripple_info.SO_phase(true_idx(idx),:)';
            % event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

            spindle_phase = ripple_info.spindle_phase(true_idx(idx),:)';
            spindle_phase_shifted = ripple_info.spindle_phase(true_idx,:)';

            for i = 1:nThresh
                th = thresholds(i);

                t1 = boot_V1 >= th;
                t2 = boot_V1 <= -th;
                
                % SO UP state with different spindle powers
                t1 = t1' + (ripple_info.normalised_UP_duration(true_idx(idx),2) > 0 & ...
                    ripple_info.spindle_amplitude(true_idx(idx),2) > spindle_thresholds(npower,2) & ripple_info.spindle_amplitude(true_idx(idx),2) < spindle_thresholds(npower+1,2)) > 1';
                t2 = t2' + (ripple_info.normalised_UP_duration(true_idx(idx),1) > 0 & ...
                    ripple_info.spindle_amplitude(true_idx(idx),1) > spindle_thresholds(npower,1) & ripple_info.spindle_amplitude(true_idx(idx),1) < spindle_thresholds(npower+1,1)) > 1';

                t1_V1 = boot_HPC(t1);
                t2_V1 = boot_HPC(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                % total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                %     sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);
                %
                % prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;


                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;

                t1s = t1s' + (ripple_info.normalised_UP_duration(true_idx,2) > 0 & ...
                    ripple_info.spindle_amplitude(true_idx,2) > spindle_thresholds(npower,2) & ripple_info.spindle_amplitude(true_idx,2) < spindle_thresholds(npower+1,2)) > 1';
                t2s = t2s' + (ripple_info.normalised_UP_duration(true_idx,1) > 0 & ...
                    ripple_info.spindle_amplitude(true_idx,1) > spindle_thresholds(npower,1) & ripple_info.spindle_amplitude(true_idx,1) < spindle_thresholds(npower+1,1)) > 1';


                t1_V1 = boot_HPC(t1s);
                t2_V1 = boot_HPC(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                % prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
            end
            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shifted;
        end

        % Quantile-based AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';


        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end

% load
% -------- Plot --------
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_spindle_power.mat'),'AUC')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_UP_spindle_power.mat'))

% spindle_thresholds = {'SO trough and spindle 0-25','SO trough and spindle 25-50','SO trough and spindle 50-75','SO trough and spindle 75-100'};
spindle_thresholds = {'UP and spindle low','UP and spindle high'};

fig = figure('Name','Temporal HPC log-odds AUC by UP with different spindle power (0.1s win 0.02s step)','Position',[640 100 400 900]);
tiledlayout(nBins,1,'TileSpacing','compact');


for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(sprintf('%s', spindle_thresholds{npower}));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end


% Save results
% save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC SO trough with low vs high spindle power (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = [1 nBins]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title('UP with low spindle vs UP with high spindle');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.35])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);


% Save results
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])


