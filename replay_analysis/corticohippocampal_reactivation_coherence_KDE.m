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

%%%%%%%%%%%%%%%%%% Ripple info
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
for probe_no = 1:2
    spindle_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_phase = [spindle_phase{1}(:,ripples_all(1).SWS_index==1) spindle_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_phase = spindle_phase';

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

%%% SO phase
SO_phase=[];
for probe_no = 1:2
    SO_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_phase = [SO_phase{1}(:,ripples_all(1).SWS_index==1) SO_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.SO_phase = SO_phase';

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



%%%%%%
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));
singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));



%%%%% Plot distribution

%%%%%%%%%%%%%
% %%%%%%%%%%% Ripple power distribution
ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];

fig = figure('Color','w');
fig.Position=[158 358 702/2 546];
fig.Name = 'Ripple peaktime ripple power distribution';


% histogram([ripple_info.ripple_power],4:0.5:25,'Normalization','probability')
temp = ripple_info.ripple_power;
temp_thresholds = prctile(temp,[25 50 75]);

% Ipsi amplitude
nexttile

histogram(temp,4:0.5:23,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Ripple power (z)');  % or whatever unit your times are
ylabel('probability');
hold on
xline(temp_thresholds,'r')
title('Ripple power');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

% Contra amplitude
temp = ripple_info.ripple_power;
temp_thresholds = prctile(temp,[20 40 60 80]);
nexttile


histogram(temp,4:0.5:23,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Ripple power (z)');  % or whatever unit your times are
ylabel('probability');
hold on
xline(temp_thresholds,'r')
title('Ripple power');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)





save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')


%%%%%%%%%%%%%
% %%%%%%%%%%%
% %%%%%%%%%%% SO phase coupling
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
%%% SO phase
SO_phase=[];
for probe_no = 1:2
    SO_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

% SO_phase = [SO_phase{1} SO_phase{2}]';



fig = figure
fig.Name = 'Ripple-SO phase-locking'
fig.Position=[158 358 702 546];

temp = [SO_phase{1}(1,ripples_all(1).SWS_index==1) SO_phase{2}(2,ripples_all(2).SWS_index==1)]';

mean_angles = [];
vector_lengths = [];
pvals = [];
% angles_wrapped = mod(spindle_phase, 2*pi);

nexttile
h = polarhistogram(temp, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
rmax = max(h.Values);  % get max count from histogram

mean_length = circ_r(temp);
mean_angle = circ_mean(temp);
pvals(end+1) = circ_rtest(temp);

% Scale mean length to match histogram scale
scaled_length = mean_length * rmax;

% Plot mean vector
hold on
polarplot([mean_angle mean_angle], [0 scaled_length], 'r-', 'LineWidth', 2);
txt = sprintf('angle = %.2f,r = %.2f, p = %.4f', mean_angle,mean_length, pvals);
title(txt);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
% figure
mean_angles = [];
vector_lengths = [];
pvals = [];
for i = 1:max(session_count)
    % nexttile
    phases = temp(session_count == i);
    angles_wrapped = mod(phases, 2*pi);
    % polarhistogram(angles_wrapped, 24);  % 12 bins (adjust as needed)
    % title('Polar Histogram of Mean Angles');

    mean_angles(i) = circ_mean(phases);
    vector_lengths(i) = circ_r(phases);  % mean resultant length
    pvals(i) = circ_rtest(phases);
end
angles_wrapped = mod(mean_angles, 2*pi);
polarhistogram(angles_wrapped, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsi')



temp = [SO_phase{1}(2,ripples_all(1).SWS_index==1) SO_phase{2}(1,ripples_all(2).SWS_index==1)]';

mean_angles = [];
vector_lengths = [];
pvals = [];
% angles_wrapped = mod(spindle_phase, 2*pi);

nexttile
h = polarhistogram(temp, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
rmax = max(h.Values);  % get max count from histogram

mean_length = circ_r(temp);
mean_angle = circ_mean(temp);
pvals(end+1) = circ_rtest(temp);

% Scale mean length to match histogram scale
scaled_length = mean_length * rmax;

% Plot mean vector
hold on
polarplot([mean_angle mean_angle], [0 scaled_length], 'r-', 'LineWidth', 2);
txt = sprintf('angle = %.2f,r = %.2f, p = %.4f', mean_angle,mean_length, pvals);
title(txt);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
% figure
mean_angles = [];
vector_lengths = [];
pvals = [];
for i = 1:max(session_count)
    % nexttile
    phases = temp(session_count == i);
    angles_wrapped = mod(phases, 2*pi);
    % polarhistogram(angles_wrapped, 24);  % 12 bins (adjust as needed)
    % title('Polar Histogram of Mean Angles');

    mean_angles(i) = circ_mean(phases);
    vector_lengths(i) = circ_r(phases);  % mean resultant length
    pvals(i) = circ_rtest(phases);
end
angles_wrapped = mod(mean_angles, 2*pi);
polarhistogram(angles_wrapped, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsi')


save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')


%%%%%%%%%%%%%
% %%%%%%%%%%%
% %%%%%%%%%%% Spindle phase coupling

%%% spindle phase
spindle_phase=[];
for probe_no = 1:2
    spindle_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];

% spindle_phase = [spindle_phase{1} spindle_phase{2}]';


fig = figure
fig.Name = 'Ripple-spindle phase-locking';
fig.Position=[158 358 702 546];
temp = [spindle_phase{1}(1,ripples_all(1).SWS_index==1) spindle_phase{2}(2,ripples_all(2).SWS_index==1)]';


mean_angles = [];
vector_lengths = [];
pvals = [];
% angles_wrapped = mod(spindle_phase, 2*pi);

nexttile
h = polarhistogram(temp, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
rmax = max(h.Values);  % get max count from histogram

mean_length = circ_r(temp);
mean_angle = circ_mean(temp);
pvals(end+1) = circ_rtest(temp);

% Scale mean length to match histogram scale
scaled_length = mean_length * rmax;

% Plot mean vector
hold on
polarplot([mean_angle mean_angle], [0 scaled_length], 'r-', 'LineWidth', 2);
txt = sprintf('angle = %.2f,r = %.2f, p = %.4f', mean_angle,mean_length, pvals);
title(txt);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
% figure
mean_angles = [];
vector_lengths = [];
pvals = [];
for i = 1:max(session_count)
    % nexttile
    phases = temp(session_count == i);
    angles_wrapped = mod(phases, 2*pi);
    % polarhistogram(angles_wrapped, 24);  % 12 bins (adjust as needed)
    % title('Polar Histogram of Mean Angles');

    mean_angles(i) = circ_mean(phases);
    vector_lengths(i) = circ_r(phases);  % mean resultant length
    pvals(i) = circ_rtest(phases);
end
angles_wrapped = mod(mean_angles, 2*pi);
polarhistogram(angles_wrapped, 'BinWidth', pi/8);  % 12 bins (adjust as needed)

set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsi')

temp = [spindle_phase{1}(2,ripples_all(1).SWS_index==1) spindle_phase{2}(1,ripples_all(2).SWS_index==1)]';


mean_angles = [];
vector_lengths = [];
pvals = [];
% angles_wrapped = mod(spindle_phase, 2*pi);

nexttile
h = polarhistogram(temp, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
rmax = max(h.Values);  % get max count from histogram

mean_length = circ_r(temp);
mean_angle = circ_mean(temp);
pvals(end+1) = circ_rtest(temp);

% Scale mean length to match histogram scale
scaled_length = mean_length * rmax;

% Plot mean vector
hold on
polarplot([mean_angle mean_angle], [0 scaled_length], 'r-', 'LineWidth', 2);
txt = sprintf('angle = %.2f,r = %.2f, p = %.4f', mean_angle,mean_length, pvals);
title(txt);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
% figure
mean_angles = [];
vector_lengths = [];
pvals = [];
for i = 1:max(session_count)
    % nexttile
    phases = temp(session_count == i);
    angles_wrapped = mod(phases, 2*pi);
    % polarhistogram(angles_wrapped, 24);  % 12 bins (adjust as needed)
    % title('Polar Histogram of Mean Angles');

    mean_angles(i) = circ_mean(phases);
    vector_lengths(i) = circ_r(phases);  % mean resultant length
    pvals(i) = circ_rtest(phases);
end
angles_wrapped = mod(mean_angles, 2*pi);
polarhistogram(angles_wrapped, 'BinWidth', pi/8);  % 12 bins (adjust as needed)
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contra')


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
spindle_amplitude=[];
for probe_no = 1:2
    spindle_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_amplitude{probe_no} = [spindle_amplitude{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end


fig = figure('Color','w');
fig.Position=[158 358 702/2 546];
fig.Name = 'Ripple peaktime spindle amplitude distribution';


temp = [spindle_amplitude{1}(1,ripples_all(1).SWS_index==1) spindle_amplitude{2}(2,ripples_all(2).SWS_index==1)];
temp_thresholds = prctile(temp,[25 50 75]);

% Ipsi amplitude
nexttile

histogram(temp,-2:0.05:5,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Amplitude (z)');  % or whatever unit your times are
ylabel('probability');
hold on
xline(temp_thresholds,'r')
title('ipsi amplitude');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

% Contra amplitude
temp = [spindle_amplitude{1}(2,ripples_all(1).SWS_index==1) spindle_amplitude{2}(1,ripples_all(2).SWS_index==1)];
temp_thresholds = prctile(temp,[25 50 75]);
nexttile

histogram(temp,-2:0.05:5,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Amplitude (z)');  % or whatever unit your times are
ylabel('probability');
hold on
xline(temp_thresholds,'r')
title('contra amplitude');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)




%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
SO_amplitude=[];
for probe_no = 1:2
    SO_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_amplitude{probe_no} = [SO_amplitude{probe_no} ripples_all(probe_no).SO_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end


fig = figure('Color','w');
fig.Position=[158 358 702/2 546];
fig.Name = 'Ripple peaktime SO amplitude distribution';


temp = [SO_amplitude{1}(1,ripples_all(1).SWS_index==1) SO_amplitude{2}(2,ripples_all(2).SWS_index==1)];
temp_thresholds = prctile(temp,[25 50 75]);

% Ipsi amplitude
nexttile

histogram(temp,-2:0.05:5,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Amplitude (z)');  % or whatever unit your times are
ylabel('probability');
hold on
xline(temp_thresholds,'r')
title('ipsi amplitude');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

% Contra amplitude
temp = [SO_amplitude{1}(2,ripples_all(1).SWS_index==1) SO_amplitude{2}(1,ripples_all(2).SWS_index==1)];
temp_thresholds = prctile(temp,[25 50 75]);
nexttile

histogram(temp,-2:0.05:5,'EdgeColor','none','Normalization','probability'); % You can change 50 to control bin numbers
xlabel('Amplitude (z)');  % or whatever unit your times are
ylabel('probability');
hold on
xline(temp_thresholds,'r')
title('contra amplitude');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')




%%%%%
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

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';



%%%%%%%%%%% ALL
bins_to_use = bin_centers>0 & bin_centers<0.1;

log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

nSamples = min([length(T1_events) length(T2_events)]);

[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', 1:length(z_bias_V1),'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', 1:length(z_bias_V1),'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);



nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 All'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%%%%%% cross_corr_bias between HPC bias and V1 bias
max_lag_bins = 100;
num_shuffles = 1000;
num_boot = 1000;

log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

[lags, real_mean, real_CI, shuffle_mean, shuffle_CI] = ...
    xcorr_bias(z_bias(:,[T1_events T2_events])', z_bias_V1(:,[T1_events T2_events])', max_lag_bins, num_shuffles, num_boot,'use_gpu',true);

time_lags = lags * bin_size_sec;

figure; hold on;

% Real ± 95% bootstrap CI
fill([time_lags fliplr(time_lags)], ...
     [real_CI(2,:) fliplr(real_CI(1,:))], ...
     [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Shuffle ± 95% CI
fill([time_lags fliplr(time_lags)], ...
     [shuffle_CI(2,:) fliplr(shuffle_CI(1,:))], ...
     'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot lines
plot(time_lags, real_mean, 'Color', [231,41,138]/256, 'LineWidth', 2);
plot(time_lags, shuffle_mean, 'k', 'LineWidth', 1.5);

xline(0, '--k');
xlabel('Lag (s)');
ylabel('Spearman correlation');
title('HPC–V1 Reactivation Bias Cross-Correlation');
legend({'Real 95% CI', 'Shuffle 95% CI', 'Real Mean', 'Shuffle Mean'});
grid on;



%%%%%%%%%%%%%%%%%%%%%%% within session temporal reacitvation
% Unique subject IDs
unique_subjects = unique(subject_id);
% unique_subjects = unique(session_count);

for subjIdx = 1:length(unique_subjects)
    subj_id = unique_subjects(subjIdx);

    % Step 1: Filter Events
    bins_to_use = bin_centers > -0 & bin_centers < 0.1;
    power_thresholds = prctile(ripple_info.ripple_power, 0:99.9/4:99.9);


    power_index = ripple_info.ripple_power > power_thresholds(end-1) & ripple_info.ripple_power < power_thresholds(end);
    % event_index = find((singlet_index + power_index) > 1);
    %
    % % Restrict to current subject
    % subj_events = find(subject_id == subj_id);
    % event_index = intersect(event_index, subj_events);
    %
    % T1_events = intersect(T1_all, event_index);
    % T2_events = intersect(T2_all, event_index);

    % power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
    % Define logical arrays for whether phase is in the peak range
    % is_peak_phase_1 = ripple_info.SO_phase(:,1) >= 0 & ripple_info.SO_phase(:,1) <= pi;
    %
    % is_peak_phase_2 = ripple_info.SO_phase(:,2) >= 0 & ripple_info.SO_phase(:,2) <= pi;


    % event_index = find((singlet_index)==1);
    subj_events = find(subject_id == subj_id);
    event_index = subj_events;
    % event_index = find(power_index);
    % event_index = intersect(event_index, subj_events);

    % hemi_index = ripple_info.spindle_presence_hemi;

    T1_events = find(nanmean(z_bias(bins_to_use, :)) > prctile(nanmean(z_bias(bins_to_use, event_index)),80))';
    T2_events = find(nanmean(z_bias(bins_to_use, :)) < prctile(nanmean(z_bias(bins_to_use, event_index)),20))';

    T1_events = intersect(T1_events,event_index);
    T2_events = intersect(T2_events,event_index);

    % % T1 events: keep those coupled to spindle peak in hemisphere 2
    % T1_events = intersect(T1_events, find(is_peak_phase_2));
    %
    % % T2 events: keep those coupled to spindle peak in hemisphere 1
    % T2_events = intersect(T2_events, find(is_peak_phase_1));

    nSamples = min([length(T1_events), length(T2_events)]);
    if nSamples < 5
        fprintf('Skipping subject %d (not enough events).\n', subj_id);
        continue;
    end

    % Step 2: Bootstrap
    [bV1_mean_T1, LCI1_V1, UCI1_V1, bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events, 'nSamples', nSamples);
    [bV1_mean_T2, LCI2_V1, UCI2_V1, bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events, 'nSamples', nSamples);
    [bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events], 'nSamples', nSamples);
    [~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events], 'nSamples', nSamples, 'nseed', 1000);

    bootV1_diff = bootV1_T1 - bootV1_T2;
    bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;
    bV1_mean_diff = mean(bootV1_diff);
    LCI_V1_diff = prctile(bootV1_diff, 2.5);
    UCI_V1_diff = prctile(bootV1_diff, 97.5);
    bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
    LCI_V1_diff_shuffled = prctile(bootV1_diff_shuffled, 2.5);
    UCI_V1_diff_shuffled = prctile(bootV1_diff_shuffled, 97.5);

    [b_mean_T1, LCI1, UCI1, boot_T1] = calculate_bootstrap_CI(z_bias', T1_events, 'nSamples', nSamples);
    [b_mean_T2, LCI2, UCI2, boot_T2] = calculate_bootstrap_CI(z_bias', T2_events, 'nSamples', nSamples);
    [b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled, boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events], 'nSamples', nSamples);
    [~, ~, ~, boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events], 'nSamples', nSamples, 'nseed', 1000);

    boot_diff = boot_T1 - boot_T2;
    boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;
    b_mean_diff = mean(boot_diff);
    LCI_diff = prctile(boot_diff, 2.5);
    UCI_diff = prctile(boot_diff, 97.5);
    b_mean_diff_shuffled = mean(boot_diff_shuffled);
    LCI_diff_shuffled = prctile(boot_diff_shuffled, 2.5);
    UCI_diff_shuffled = prctile(boot_diff_shuffled, 97.5);

    % Step 3: Plot
    nfig = figure;
    nfig.Name = sprintf('KDE Bias PSTH - Subject %d', subj_id);
    nfig.Position = [500 75 700 850];

    event_index = [T1_events T2_events];
    z_event = nanmean(z_bias(bins_to_use,[T1_events T2_events]));

    [~,sorted_index] = sort(z_event);

    nexttile
    h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
    set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
    xlim([-0.5 0.5])
    colorbar;
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    % imagesc(z_bias_V1);colorbar;
    xlabel('Time (s)')
    ylabel('Ripple event')

    nexttile
    h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
    set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
    xlim([-0.5 0.5])
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
    fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(bin_centers, b_mean_T1, 'Color', track1_color);
    plot(bin_centers, b_mean_T2, 'Color', track2_color);
    plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
    xlim([-0.5 0.5])
    ylim([-3 3])
    legend('Track 1', 'Track 2','Shuffled','box','off');
    xlabel('Time (s)')
    ylabel('Bias (z)')
    title('z bias');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    % Tile 6: z_bias_V1 mean with CI

    nexttile
    hold on;
    fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
    plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
    plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
    ylim([-0.2 0.2])
    xlim([-0.5 0.5])
    legend('Track 1', 'Track 2','Shuffled','box','off');
    title('z bias V1');
    xlabel('Time (s)')
    ylabel('Bias (z)')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



    nexttile
    hold on;

    fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
    plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
    xlim([-0.5 0.5])
    ylim([-1 4.5])
    legend('Track diff','Shuffled','box','off');
    xlabel('Time (s)')
    ylabel('Bias T1-T2 difference (z)')
    title('z bias diff');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    % Tile 6: z_bias_V1 mean with CI

    nexttile
    hold on;
    fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
    plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
    ylim([-0.1 0.3])
    xlim([-0.5 0.5])
    legend('Track diff','Shuffled','box','off');
    title('z bias diff V1');
    xlabel('Time (s)')
    ylabel('Bias T1-T2 difference (z)')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


end


%%%%%%%%%%%%%%%%%%%%%%%% xcorr
% Define folder for saving output
output_dir = 'reactivation_bias_xcorr';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Unique subject IDs
unique_subjects = unique(subject_id);
KDE_bias_xcorr = struct(); % Pre-allocate

% Loop through each subject
for subjIdx = 1:length(unique_subjects)
    subj_id = unique_subjects(subjIdx);
    fprintf('Processing subject %d (%d of %d)\n', subj_id, subjIdx, length(unique_subjects));

    % Filter event indices for this subject
    bins_to_use = bin_centers > 0 & bin_centers < 0.1;
    subj_events = find(subject_id == subj_id);

    % Define Track 1 and Track 2 events based on z_bias
    event_index = subj_events;
    bias_scores = nanmean(z_bias(bins_to_use, :), 1);
    T1_events = find(bias_scores > prctile(bias_scores(event_index), 75))';
    T2_events = find(bias_scores < prctile(bias_scores(event_index), 25))';

    T1_events = intersect(T1_events, event_index);
    T2_events = intersect(T2_events, event_index);
    nSamples = min([length(T1_events), length(T2_events)]);

    if nSamples < 5
        fprintf('Skipping subject %d (not enough events).\n', subj_id);
        continue;
    end

    % --- Bootstrap Bias Traces ---
    [bV1_mean_T1, LCI1_V1, UCI1_V1, bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events, 'nSamples', nSamples);
    [bV1_mean_T2, LCI2_V1, UCI2_V1, bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events, 'nSamples', nSamples);
    [bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = ...
        calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events], 'nSamples', nSamples);
    [~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events], 'nSamples', nSamples, 'nseed', 1000);

    bootV1_diff = bootV1_T1 - bootV1_T2;
    bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;
    bV1_mean_diff = mean(bootV1_diff);
    LCI_V1_diff = prctile(bootV1_diff, 2.5);
    UCI_V1_diff = prctile(bootV1_diff, 97.5);
    bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
    LCI_V1_diff_shuffled = prctile(bootV1_diff_shuffled, 2.5);
    UCI_V1_diff_shuffled = prctile(bootV1_diff_shuffled, 97.5);

    [b_mean_T1, LCI1, UCI1, boot_T1] = calculate_bootstrap_CI(z_bias', T1_events, 'nSamples', nSamples);
    [b_mean_T2, LCI2, UCI2, boot_T2] = calculate_bootstrap_CI(z_bias', T2_events, 'nSamples', nSamples);
    [b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled, boot_T1_shuffled] = ...
        calculate_bootstrap_CI(z_bias', [T1_events; T2_events], 'nSamples', nSamples);
    [~, ~, ~, boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events], 'nSamples', nSamples, 'nseed', 1000);

    boot_diff = boot_T1 - boot_T2;
    boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;
    b_mean_diff = mean(boot_diff);
    LCI_diff = prctile(boot_diff, 2.5);
    UCI_diff = prctile(boot_diff, 97.5);
    b_mean_diff_shuffled = mean(boot_diff_shuffled);
    LCI_diff_shuffled = prctile(boot_diff_shuffled, 2.5);
    UCI_diff_shuffled = prctile(boot_diff_shuffled, 97.5);

    % --- Cross-Correlation (with parallel bootstrap & shuffle) ---
    hpc_bias = z_bias(:, event_index)';     % N_events x T
    v1_bias = z_bias_V1(:, event_index)';   % N_events x T
    max_lag_bins = 50;
    num_shuffles = 1000;
    num_boot = 1000;

    [lags, corr_mean, corr_CI, shuffle_mean, shuffle_CI] = ...
        xcorr_bias(hpc_bias, v1_bias, max_lag_bins, num_shuffles, num_boot);

    % --- Save into structure
    KDE_bias_xcorr(subjIdx).subj_id = subj_id;
    KDE_bias_xcorr(subjIdx).lags = lags;
    KDE_bias_xcorr(subjIdx).corr_mean = corr_mean;
    KDE_bias_xcorr(subjIdx).corr_CI = corr_CI;
    KDE_bias_xcorr(subjIdx).shuffle_mean = shuffle_mean;
    KDE_bias_xcorr(subjIdx).shuffle_CI = shuffle_CI;

    % --- Plot and save cross-correlation
    time_lags = lags * 0.01;
    fig_xcorr = figure('Name', sprintf('HPC-V1 Bias Cross-Corr – Subject %d', subj_id), ...
                       'Position', [200 100 700 450]);
    hold on;
    fill([time_lags fliplr(time_lags)], [corr_CI(2,:) fliplr(corr_CI(1,:))], ...
         [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([time_lags fliplr(time_lags)], [shuffle_CI(2,:) fliplr(shuffle_CI(1,:))], ...
         'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(time_lags, corr_mean, 'Color', [231,41,138]/256, 'LineWidth', 2);
    plot(time_lags, shuffle_mean, 'k', 'LineWidth', 1.5);
    xline(0, '--k');
    xlabel('Lag (s)'); ylabel('Spearman correlation');
    title(sprintf('HPC–V1 Bias Cross-Corr (Subject %d)', subj_id));
    legend({'Real CI', 'Shuffle CI', 'Real', 'Shuffle'});
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    grid on;

    saveas(fig_xcorr, fullfile(output_dir, sprintf('KDE_bias_xcorr_subj%d.png', subj_id)));
    close(fig_xcorr);
end

% --- Save all results
save(fullfile(output_dir, 'KDE_bias_xcorr_all_subjects.mat'), 'KDE_bias_xcorr');



%%
%%%% Low ripple power

log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

bins_to_use = bin_centers>0 & bin_centers<0.1;
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/4:99.9);

power_index = ripple_info.ripple_power > power_thresholds(1) & ripple_info.ripple_power < power_thresholds(2);
event_index = find(power_index ==1);
% event_index = find((singlet_index+power_index)>1);
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 Low ripple power'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%% High ripple power
bins_to_use = bin_centers>0 & bin_centers<0.1;
% power_thresholds = prctile(ripple_info.ripple_power,0:99.9/3:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));
power_index = ripple_info.ripple_power > power_thresholds(end-1) & ripple_info.ripple_power < power_thresholds(end);
% event_index = find((singlet_index+power_index)>1);
event_index = find(power_index ==1);
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 High ripple power'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%% Spindle power
%%%% with low spindle power

% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));
spindle_index1 = ripple_info.spindle_amplitude(:,1) < prctile(ripple_info.spindle_amplitude(:,1),[25]);
spindle_index2 = ripple_info.spindle_amplitude(:,2) < prctile(ripple_info.spindle_amplitude(:,2),[25]);
% spindle_index = (spindle_index1 + spindle_index2) > 0;

% T1_events = intersect(T1_events,find(spindle_index));
% T2_events = intersect(T2_events,find(spindle_index));
T1_events = intersect(T1_events,find(spindle_index2));
T2_events = intersect(T2_events,find(spindle_index1));

% spindle_index = ripple_info.spindle_amplitude < prctile(ripple_info.spindle_amplitude,[25]);
% 
% event_index = find((singlet_index+spindle_index)>1);
% hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with low spindle power'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% with high spindle power
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% spindle_index = ripple_info.spindle_presence == 1;
spindle_index1 = ripple_info.spindle_amplitude(:,1) > prctile(ripple_info.spindle_amplitude(:,1),[75]);
spindle_index2 = ripple_info.spindle_amplitude(:,2) > prctile(ripple_info.spindle_amplitude(:,2),[75]);

% spindle_index = (spindle_index1 + spindle_index2) > 0;

% T1_events = intersect(T1_events,find(spindle_index));
% T2_events = intersect(T2_events,find(spindle_index));
% 
T1_events = intersect(T1_events,find(spindle_index2));
T2_events = intersect(T2_events,find(spindle_index1));

% event_index = find((singlet_index+spindle_index)>1);
% hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 2));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 1));



nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with high spindle power'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%
%%%% Spindle peak

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));


% Define logical arrays for whether phase is in the peak range
% is_peak_phase_1 = ripple_info.spindle_phase(:,1) >= 0 & ripple_info.spindle_phase(:,1) <= pi;
% 
% is_peak_phase_2 = ripple_info.spindle_phase(:,2) >= 0 & ripple_info.spindle_phase(:,2) <= pi;

is_peak_phase_1 = ripple_info.spindle_phase(:,1) >= -pi/2 & ripple_info.spindle_phase(:,1) <= pi/2;

is_peak_phase_2 = ripple_info.spindle_phase(:,2) >= -pi/2 & ripple_info.spindle_phase(:,2) <= pi/2;


% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_peak_phase_2));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_peak_phase_1));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with spindle phase peak'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%% Spindle trough
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
% is_trough_phase_1 = ripple_info.spindle_phase(:,1) >= -pi & ripple_info.spindle_phase(:,1) < 0;
% 
% is_trough_phase_2 = ripple_info.spindle_phase(:,2) >= -pi & ripple_info.spindle_phase(:,2) < 0;

is_trough_phase_1 = ripple_info.spindle_phase(:,1) >= -pi & ripple_info.spindle_phase(:,1) < -pi/2 |  ripple_info.spindle_phase(:,1) <= pi & ripple_info.spindle_phase(:,1) > pi/2;
is_trough_phase_2 = ripple_info.spindle_phase(:,2) >= -pi & ripple_info.spindle_phase(:,2) < -pi/2 |  ripple_info.spindle_phase(:,2) <= pi & ripple_info.spindle_phase(:,2) > pi/2;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_trough_phase_2));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_trough_phase_1));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);

MUA_bin_centers = []
nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with spindle phase trough'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);









%%%% Spindle peak in opposite hemisphere

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
is_peak_phase_1 = ripple_info.spindle_phase(:,1) >= 0 & ripple_info.spindle_phase(:,1) <= pi;

is_peak_phase_2 = ripple_info.spindle_phase(:,2) >= 0 & ripple_info.spindle_phase(:,2) <= pi;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_peak_phase_1));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_peak_phase_2));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with spindle phase peak in opposite hemi'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);







%%

%%%% SO peak

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
% is_peak_phase_1 = ripple_info.SO_phase(:,1) >= -pi & ripple_info.SO_phase(:,1) <= 0;
% 
% is_peak_phase_2 = ripple_info.SO_phase(:,2) >= -pi & ripple_info.SO_phase(:,2) <= 0;
% 

is_peak_phase_1 = ripple_info.SO_phase(:,1) >= -pi/2 & ripple_info.SO_phase(:,1) <= pi/2;

is_peak_phase_2 = ripple_info.SO_phase(:,2) >= -pi/2 & ripple_info.SO_phase(:,2) <= pi/2;


% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_peak_phase_2));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_peak_phase_1));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase peak'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% SO trough
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
is_trough_phase_1 = (ripple_info.SO_phase(:,1) >0 & ripple_info.SO_phase(:,1) <pi);

is_trough_phase_2 = (ripple_info.SO_phase(:,2) >0 & ripple_info.SO_phase(:,2) <pi);

is_trough_phase_1 = ripple_info.SO_phase(:,1) >= -pi & ripple_info.SO_phase(:,1) < -pi/2 |  ripple_info.SO_phase(:,1) <= pi & ripple_info.SO_phase(:,1) > pi/2;
is_trough_phase_2 = ripple_info.SO_phase(:,2) >= -pi & ripple_info.SO_phase(:,2) < -pi/2 |  ripple_info.SO_phase(:,2) <= pi & ripple_info.SO_phase(:,2) > pi/2;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_trough_phase_2));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_trough_phase_1));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase trough'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);






%%%%%%%%%%%% SO synchronised

%%%% SO peak 

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
% is_peak_phase_1 = ripple_info.SO_phase(:,1) >= -pi & ripple_info.SO_phase(:,1) <= 0;
% 
% is_peak_phase_2 = ripple_info.SO_phase(:,2) >= -pi & ripple_info.SO_phase(:,2) <= 0;
% 

is_peak_phase_1 = ripple_info.SO_phase(:,1) >= -pi/2 & ripple_info.SO_phase(:,1) <= pi/2;
is_peak_phase_2 = ripple_info.SO_phase(:,2) >= -pi/2 & ripple_info.SO_phase(:,2) <= pi/2;
is_peak_phase = is_peak_phase_1 + is_peak_phase_2==2;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_peak_phase_2+is_peak_phase==2));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_peak_phase_1+is_peak_phase==2));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase peak (in phase)'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% SO trough
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
% is_trough_phase_1 = (ripple_info.SO_phase(:,1) >0 & ripple_info.SO_phase(:,1) <pi);
% 
% is_trough_phase_2 = (ripple_info.SO_phase(:,2) >0 & ripple_info.SO_phase(:,2) <pi);

is_trough_phase_1 = ripple_info.SO_phase(:,1) >= -pi & ripple_info.SO_phase(:,1) < -pi/2 |  ripple_info.SO_phase(:,1) <= pi & ripple_info.SO_phase(:,1) > pi/2;
is_trough_phase_2 = ripple_info.SO_phase(:,2) >= -pi & ripple_info.SO_phase(:,2) < -pi/2 |  ripple_info.SO_phase(:,2) <= pi & ripple_info.SO_phase(:,2) > pi/2;
is_trough_phase = is_trough_phase_1 + is_trough_phase_2==2;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_trough_phase_2+is_trough_phase==2));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_trough_phase_1+is_trough_phase==2));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase trough (in phase)'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%%%%% SO not synchronised (out of phase)
%%%% SO peak 

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
% is_peak_phase_1 = ripple_info.SO_phase(:,1) >= -pi & ripple_info.SO_phase(:,1) <= 0;
% 
% is_peak_phase_2 = ripple_info.SO_phase(:,2) >= -pi & ripple_info.SO_phase(:,2) <= 0;
% 

is_peak_phase_1 = ripple_info.SO_phase(:,1) >= -pi/2 & ripple_info.SO_phase(:,1) <= pi/2;
is_peak_phase_2 = ripple_info.SO_phase(:,2) >= -pi/2 & ripple_info.SO_phase(:,2) <= pi/2;
is_peak_phase = is_peak_phase_1 + is_peak_phase_2==2;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_peak_phase_2+is_peak_phase==1));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_peak_phase_1+is_peak_phase==1));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase peak (out phase)'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% SO trough
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% Define logical arrays for whether phase is in the peak range
% is_trough_phase_1 = (ripple_info.SO_phase(:,1) >0 & ripple_info.SO_phase(:,1) <pi);
% 
% is_trough_phase_2 = (ripple_info.SO_phase(:,2) >0 & ripple_info.SO_phase(:,2) <pi);

is_trough_phase_1 = ripple_info.SO_phase(:,1) >= -pi & ripple_info.SO_phase(:,1) < -pi/2 |  ripple_info.SO_phase(:,1) <= pi & ripple_info.SO_phase(:,1) > pi/2;
is_trough_phase_2 = ripple_info.SO_phase(:,2) >= -pi & ripple_info.SO_phase(:,2) < -pi/2 |  ripple_info.SO_phase(:,2) <= pi & ripple_info.SO_phase(:,2) > pi/2;
is_trough_phase = is_trough_phase_1 + is_trough_phase_2==2;

% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_trough_phase_2+is_trough_phase==1));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_trough_phase_1+is_trough_phase==1));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase trough (out phase)'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%% Delta power



log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));


% SO_index1 = ripple_info.spindle_amplitude(:,1) > prctile(ripple_info.spindle_amplitude(:,1),[75]);
% SO_index2 = ripple_info.spindle_amplitude(:,2) > prctile(ripple_info.spindle_amplitude(:,2),[75]);
SO_index1 = ripple_info.SO_amplitude(:,1) < prctile(ripple_info.SO_amplitude(:,1),[25]);
SO_index2 = ripple_info.SO_amplitude(:,2) < prctile(ripple_info.SO_amplitude(:,2),[25]);

T1_events = intersect(T1_events,find(SO_index2));
T2_events = intersect(T2_events,find(SO_index1));


% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);



nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with low delta power'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%%% High delta power

log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));


SO_index1 = ripple_info.SO_amplitude(:,1) > prctile(ripple_info.SO_amplitude(:,1),[75]);
SO_index2 = ripple_info.SO_amplitude(:,2) > prctile(ripple_info.SO_amplitude(:,2),[75]);


T1_events = intersect(T1_events,find(SO_index2));
T2_events = intersect(T2_events,find(SO_index1));


% event_index = find((singlet_index)==1);
% % hemi_index = ripple_info.spindle_presence_hemi;
% T1_events = intersect(T1_events,event_index);
% T2_events = intersect(T2_events,event_index);



nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with high delta power'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

%%

%%%% UP transition

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

% event_index = find((singlet_index+ripple_info.UP_index)>1);
hemi_index = ripple_info.UP_index_hemi;

T1_events = intersect(T1_events,find(hemi_index >0));
T2_events = intersect(T2_events,find(hemi_index >0));



nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 in UP'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% DOWN transition

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));


% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

% event_index = find((singlet_index+ripple_info.DOWN_index)>1);
hemi_index = ripple_info.DOWN_index_hemi;

T1_events = intersect(T1_events,find(hemi_index >0));
T2_events = intersect(T2_events,find(hemi_index >0));



nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with DOWN'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% Early UP

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

% event_index = find((singlet_index+ripple_info.early_UP_index)>1);
hemi_index = ripple_info.early_UP_index_hemi;

T1_events = intersect(T1_events,find(hemi_index>0));
T2_events = intersect(T2_events,find(hemi_index >0));


nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with early UP'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% Late UP

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
log_odds_threshold = prctile(nanmean(z_bias(bins_to_use,:)),[20 80]);
T1_events = find(nanmean(z_bias(bins_to_use,:))>log_odds_threshold(2));
T2_events = find(nanmean(z_bias(bins_to_use,:))<log_odds_threshold(1));

% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

% event_index = find((singlet_index+ripple_info.late_UP_index)>1
event_index = find(ripple_info.late_UP_index>0);
hemi_index = ripple_info.late_UP_index_hemi;

T1_events = intersect(intersect(T1_events,event_index),find(hemi_index>0));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index>0));



nSamples = min([length(T1_events) length(T2_events)]);
[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~, bootV1_T2_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
bootV1_diff_shuffled = bootV1_T1_shuffled - bootV1_T2_shuffled;

bV1_mean_diff_shuffled = mean(bootV1_diff_shuffled);
LCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 2.5);
UCI_V1_diff_shuffled =  prctile(bootV1_diff_shuffled, 97.5);

bootV1_diff = bootV1_T1 - bootV1_T2;
bV1_mean_diff = mean(bootV1_diff);
LCI_V1_diff =  prctile(bootV1_diff, 2.5);
UCI_V1_diff =  prctile(bootV1_diff, 97.5);

[b_mean_T1, LCI1, UCI1,boot_T1] = calculate_bootstrap_CI(z_bias', T1_events,'nSamples',nSamples);
[b_mean_T2, LCI2, UCI2,boot_T2] = calculate_bootstrap_CI(z_bias', T2_events,'nSamples',nSamples);
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples);
[~, ~, ~,boot_T2_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events; T2_events],'nSamples',nSamples,'nseed', 1000);
boot_diff_shuffled = boot_T1_shuffled - boot_T2_shuffled;

b_mean_diff_shuffled = mean(boot_diff_shuffled);
LCI_diff_shuffled =  prctile(boot_diff_shuffled, 2.5);
UCI_diff_shuffled =  prctile(boot_diff_shuffled, 97.5);

boot_diff = boot_T1 - boot_T2;
b_mean_diff = mean(boot_diff);
LCI_diff =  prctile(boot_diff, 2.5);
UCI_diff =  prctile(boot_diff, 97.5);


nfig = figure;
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with late UP'; 
nfig.Position = [500 75 700 850];

event_index = [T1_events; T2_events];
z_event = nanmean(z_bias(bins_to_use,[T1_events; T2_events]));

[~,sorted_index] = sort(z_event);

nexttile
h = imagesc(bin_centers,[],z_bias(:,event_index(sorted_index))');clim([-3 3])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
h= imagesc(bin_centers,[],z_bias_V1(:,event_index(sorted_index))');clim([-1 1])
set(h, 'AlphaData', ~isnan(z_bias(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
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
fill([bin_centers fliplr(bin_centers)], [UCI_label_shuffled fliplr(LCI_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(bin_centers, b_mean_T1, 'Color', track1_color);
plot(bin_centers, b_mean_T2, 'Color', track2_color);
plot(bin_centers, b_mean_label_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-3 3])
legend('Track 1', 'Track 2','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias (z)')
title('z bias');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI1_V1 fliplr(LCI1_V1)], track1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI2_V1 fliplr(LCI2_V1)], track2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_label_shuffled fliplr(LCI_V1_label_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_label_shuffled, 'Color', 'k');
plot(bin_centers, bV1_mean_T1, 'Color', track1_color);
plot(bin_centers, bV1_mean_T2, 'Color', track2_color);
ylim([-0.2 0.2])
xlim([-0.5 0.5])
legend('Track 1', 'Track 2','Shuffled','box','off');
title('z bias V1');
xlabel('Time (s)')
ylabel('Bias (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



nexttile
hold on;

fill([bin_centers fliplr(bin_centers)], [UCI_diff fliplr(LCI_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_diff_shuffled fliplr(LCI_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, b_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, b_mean_diff_shuffled, 'Color', 'k');
xlim([-0.5 0.5])
ylim([-1 4.5])
legend('Track diff','Shuffled','box','off');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
title('z bias diff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Tile 6: z_bias_V1 mean with CI

nexttile
hold on;
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff fliplr(LCI_V1_diff)], [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([bin_centers fliplr(bin_centers)], [UCI_V1_diff_shuffled fliplr(LCI_V1_diff_shuffled)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(bin_centers, bV1_mean_diff, 'Color', [231,41,138]/256);
plot(bin_centers, bV1_mean_diff_shuffled, 'Color', 'k');
ylim([-0.1 0.3])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])

