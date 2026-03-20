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
        spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_phase = [spindle_phase{1}(:,ripples_all(1).SWS_index==1) spindle_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_phase = spindle_phase';

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

%%% early UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,1)+0.1]);
ripple_info.early_UP_index = UP_index;
ripple_info.early_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.early_UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);


%%% late UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1)-0.1 merged_event_info.DOWN_ints(:,1)]);
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


singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));


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
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'))
load(fullfile(analysis_folder,'KDE_reactivation_ripples_all_POST.mat'))
load(fullfile(analysis_folder,'KDE_reactivation_V1_ripples_all_POST.mat'))

timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;


% z_bias1 = KDE_reactivation_content.HPC_z_ripples';
% T1_events = find(nanmean(z_bias1(1:10,:))>0.5);
% T2_events = find(nanmean(z_bias1(1:10,:))<-0.5);

z_bias = KDE_reactivation_ripples_PSTH.HPC_z_ripples';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_ripples';
% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';



%%%%%%%%%%% ALL
bins_to_use = bin_centers>0 & bin_centers<0.1;
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);
nSamples = min([length(T1_events) length(T2_events)]);

[bV1_mean_T1, LCI1_V1, UCI1_V1,bootV1_T1] = calculate_bootstrap_CI(z_bias_V1', T1_events,'nSamples',nSamples);
[bV1_mean_T2, LCI2_V1, UCI2_V1,bootV1_T2] = calculate_bootstrap_CI(z_bias_V1', T2_events,'nSamples',nSamples);
[bV1_mean_label_shuffled, LCI_V1_label_shuffled, UCI_V1_label_shuffled, bootV1_T1_shuffled] = calculate_bootstrap_CI(z_bias_V1', [T1_events T2_events],'nSamples',nSamples);
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
[b_mean_label_shuffled, LCI_label_shuffled, UCI_label_shuffled,boot_T1_shuffled] = calculate_bootstrap_CI(z_bias', [T1_events T2_events],'nSamples',nSamples);
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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%% Low ripple power
bins_to_use = bin_centers>0 & bin_centers<0.1;
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);
power_index = ripple_info.ripple_power > power_thresholds(1) & ripple_info.ripple_power < power_thresholds(2);
event_index = find((singlet_index+power_index)>1);
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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% High ripple power
bins_to_use = bin_centers>0 & bin_centers<0.1;
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.3);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.3);
power_index = ripple_info.ripple_power > power_thresholds(2) & ripple_info.ripple_power < power_thresholds(3);
event_index = find((singlet_index+power_index)>1);
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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% Without spindle

% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);
spindle_index = ripple_info.spindle_presence == 0;

event_index = find((singlet_index+spindle_index)>1);
hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));


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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 without spindles'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);
spindle_index = ripple_info.spindle_presence == 1;

event_index = find((singlet_index+spindle_index)>1);
hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 2));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 1));



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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with spindles'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% Spindle peak

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_peak_phase_1 = ripple_info.spindle_phase(:,1) >= 0 & ripple_info.spindle_phase(:,1) <= pi;

is_peak_phase_2 = ripple_info.spindle_phase(:,2) >= 0 & ripple_info.spindle_phase(:,2) <= pi;


event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
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
is_trough_phase_1 = ripple_info.spindle_phase(:,1) >= -pi & ripple_info.spindle_phase(:,1) < 0;

is_trough_phase_2 = ripple_info.spindle_phase(:,2) >= -pi & ripple_info.spindle_phase(:,2) < 0;


event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);









%%%% Spindle peak in opposite hemisphere

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_peak_phase_1 = ripple_info.spindle_phase(:,1) >= 0 & ripple_info.spindle_phase(:,1) <= pi;

is_peak_phase_2 = ripple_info.spindle_phase(:,2) >= 0 & ripple_info.spindle_phase(:,2) <= pi;

event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%% Spindle trough in opposite hemisphere
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_trough_phase_1 = ripple_info.spindle_phase(:,1) >= -pi & ripple_info.spindle_phase(:,1) < 0;

is_trough_phase_2 = ripple_info.spindle_phase(:,2) >= -pi & ripple_info.spindle_phase(:,2) < 0;

event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_trough_phase_1));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_trough_phase_2));


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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with spindle phase trough in opposite hemi'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);








%%%% SO peak

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_peak_phase_1 = ripple_info.SO_phase(:,1) >= 0 & ripple_info.SO_phase(:,1) <= pi;

is_peak_phase_2 = ripple_info.SO_phase(:,2) >= 0 & ripple_info.SO_phase(:,2) <= pi;


event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.15 0.15])
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
ylim([-1 3])
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
ylim([-0.1 0.2])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% SO trough
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_trough_phase_1 = (ripple_info.SO_phase(:,1) >-pi & ripple_info.SO_phase(:,1) <0);

is_trough_phase_2 = (ripple_info.SO_phase(:,2) >-pi & ripple_info.SO_phase(:,2) <0);

event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.15 0.15])
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
ylim([-1 3])
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
ylim([-0.1 0.2])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);






%%%% SO peak with ripple power

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.3);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.3);
power_index = ripple_info.ripple_power > power_thresholds(2) & ripple_info.ripple_power < power_thresholds(3);

% Define logical arrays for whether phase is in the peak range
is_peak_phase_1 = ripple_info.SO_phase(:,1) >= 0 & ripple_info.SO_phase(:,1) <= pi;

is_peak_phase_2 = ripple_info.SO_phase(:,2) >= 0 & ripple_info.SO_phase(:,2) <= pi;

event_index = find((singlet_index + power_index)>1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase peak + high ripple power'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-1 3])
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
ylim([-0.1 0.25])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);









%%%% SO trough + high ripple power
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.3);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.3);
power_index = ripple_info.ripple_power > power_thresholds(2) & ripple_info.ripple_power < power_thresholds(3);
% Define logical arrays for whether phase is in the peak range
is_trough_phase_1 = (ripple_info.SO_phase(:,1) >-pi & ripple_info.SO_phase(:,1) <0);

is_trough_phase_2 = (ripple_info.SO_phase(:,2) >-pi & ripple_info.SO_phase(:,2) <0);

event_index = find((singlet_index+power_index)>1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase trough + high ripple power'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.15 0.15])
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
ylim([-1 3])
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
ylim([-0.1 0.25])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);





%%%% SO peak in opposite hemisphere

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_peak_phase_1 = ripple_info.SO_phase(:,1) >= 0 & ripple_info.SO_phase(:,1) <= pi;

is_peak_phase_2 = ripple_info.SO_phase(:,2) >= 0 & ripple_info.SO_phase(:,2) <= pi;

event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 1
T1_events = intersect(T1_events, find(is_peak_phase_1));

% T2 events: keep those coupled to spindle peak in hemisphere 2
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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase peak in opposite hemi'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.15 0.15])
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
ylim([-1 3])
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
ylim([-0.1 0.2])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% SO trough in opposite hemisphere
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);

% Define logical arrays for whether phase is in the peak range
is_trough_phase_1 = (ripple_info.SO_phase(:,1) >-pi & ripple_info.SO_phase(:,1) <0);

is_trough_phase_2 = (ripple_info.SO_phase(:,2) >-pi & ripple_info.SO_phase(:,2) <0);

event_index = find((singlet_index)==1);
% hemi_index = ripple_info.spindle_presence_hemi;
T1_events = intersect(T1_events,event_index);
T2_events = intersect(T2_events,event_index);

% T1 events: keep those coupled to spindle peak in hemisphere 2
T1_events = intersect(T1_events, find(is_trough_phase_1));

% T2 events: keep those coupled to spindle peak in hemisphere 1
T2_events = intersect(T2_events, find(is_trough_phase_2));


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
nfig.Name = 'KDE Bias PSTH Track 1 vs Track 2 with SO phase trough in opposite hemi'; 
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.15 0.15])
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
ylim([-1 3])
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
ylim([-0.1 0.2])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);






%%%% UP transition

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);


% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

event_index = find((singlet_index+ripple_info.UP_index)>1);
hemi_index = ripple_info.UP_index_hemi;

T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 2));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 1));



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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);




%%%% DOWN transition

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);


% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

event_index = find((singlet_index+ripple_info.DOWN_index)>1);
hemi_index = ripple_info.DOWN_index_hemi;

T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 2));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 1));



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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% Early UP

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);


% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

event_index = find((singlet_index+ripple_info.early_UP_index)>1);
hemi_index = ripple_info.early_UP_index_hemi;

T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 2));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 1));



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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



%%%% Late UP

%%%% with Spindle
% power_thresholds = prctile(ripple_info.spindle_presence,0:99.9/2:99.9);
T1_events = find(nanmean(z_bias(bins_to_use,:))>0.5);
T2_events = find(nanmean(z_bias(bins_to_use,:))<-0.5);


% T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 0));
% T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 0));

event_index = find((singlet_index+ripple_info.late_UP_index)>1);
hemi_index = ripple_info.late_UP_index_hemi;

T1_events = intersect(intersect(T1_events,event_index),find(hemi_index == 2));
T2_events = intersect(intersect(T2_events,event_index),find(hemi_index == 1));



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
nfig.Position = [7 75 1900 910];
nexttile
imagesc(bin_centers,[],z_bias(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(bin_centers,[],z_bias(:,T2_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T1_events)');clim([-3 3])
xlim([-0.5 0.5])
colorbar;
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')

nexttile
imagesc(MUA_bin_centers,[],z_bias_V1(:,T2_events)');clim([-3 3])
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
ylim([-2 2])
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
ylim([-0.1 0.1])
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
ylim([-1 3])
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
ylim([-0.1 0.15])
xlim([-0.5 0.5])
legend('Track diff','Shuffled','box','off');
title('z bias diff V1');
xlabel('Time (s)')
ylabel('Bias T1-T2 difference (z)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])

