function plot_ipsi_contra_MUA_probability_synchronisation(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)

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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_baseline.mat'));
probability_normalised_baseline = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_baseline.mat'));
probability_psth_baseline = probability;


load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
probability_normalised_whole = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised.mat'));
probability_normalised = probability_normalised;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability.mat'));
probability_psth = probability;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_normalised_whole.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));


load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_probability.mat'));
probability_SO_SO = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_probability.mat'));
probability_SO_SO_contralateral = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_SO_contralateral_whole_probability.mat'));
probability_SO_SO_contralateral_whole = probability;


load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability.mat'));
probability_ripples_SO = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability_whole.mat'));
probability_ripples_SO_whole = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability.mat'));
probability_ripples_ripples = probability;


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
nprobe = 1;
event_averaging_scale = 10;
% extract plv, amp corr and lag (latency)
for nprobe=1:2
    mprobe = abs(nprobe-3);

    % UP and DOWN
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

        ipsi_lag_UD{nprobe} = [ipsi_lag_UD{nprobe} squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];
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
            contra_lag_spindles{nprobe} = [contra_lag_spindles{nprobe} squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];

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
        contra_lag_ripples{nprobe} = [contra_lag_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];

        ipsi_plv_ripples{nprobe} = [ipsi_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_plv_ripples{nprobe} = [contra_plv_ripples{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];

        ipsi_corr_ripples{nprobe} = [ipsi_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        contra_corr_ripples{nprobe} = [contra_corr_ripples{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];
    end
end
%% Plot ipsi vs contra latency, amp corr and phase locking


%% MUA during UP-DOWN transition

%%%%% L DOWN
event_averaging_scale = 50;

% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'Left UP-DOWN transition MUA activity';

nprobe = 1;
MUA = PSTH_MUA(nprobe).L_V1_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))

xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral V1')



nexttile
MUA = PSTH_MUA(nprobe).L_HPC_DOWN;
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))

xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral HPC')


MUA = PSTH_MUA(nprobe).R_V1_DOWN;
nexttile
% [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral V1')


MUA = PSTH_MUA(nprobe).R_HPC_DOWN;
nexttile
% [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral HPC')


nexttile
% nprobe = 1;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(PSTH_MUA(nprobe).L_V1_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(PSTH_MUA(nprobe).R_V1_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(PSTH_MUA_baseline(nprobe).L_V1_DOWN,'omitnan');
LCI = prctile(PSTH_MUA_baseline(nprobe).L_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA_baseline(nprobe).L_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
% xline(0,'r')
title('V1 MUA during Left DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
% nprobe = 1;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(PSTH_MUA(nprobe).L_HPC_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(PSTH_MUA(nprobe).R_HPC_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(PSTH_MUA_baseline(nprobe).L_HPC_DOWN,'omitnan');
LCI = prctile(PSTH_MUA_baseline(nprobe).L_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA_baseline(nprobe).L_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
% xline(0,'r')
title('HPC MUA during Left DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%%%% R DOWN
event_averaging_scale = 50;

% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
% colour_lines = [44,123,182;215,25,28]/256;% Blue Red
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
fig = figure('Color','w');
fig.Position = [350 59 1100 930];
fig.Name = 'Right UP-DOWN transition MUA activity';

nprobe = 2;
MUA = PSTH_MUA(nprobe).R_V1_DOWN;

nexttile
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))

xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral V1')


nexttile
MUA = PSTH_MUA(nprobe).R_HPC_DOWN;
[~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))

xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Ipsilateral HPC')

MUA = PSTH_MUA(nprobe).L_V1_DOWN;
nexttile
% [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral V1')

MUA = PSTH_MUA(nprobe).L_HPC_DOWN;
nexttile
% [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(movmean(MUA(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 50.5 100.5 150.5 200.5])
xticklabels([-1 -0.5 0 0.5 1])
xline(100.5,'r',LineWidth=1)
clim([-2 2])
colorbar
colormap(flipud(gray))
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('Event sorted by DOWN duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral HPC')


nexttile
% nprobe = 2;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(PSTH_MUA(nprobe).R_V1_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).R_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(PSTH_MUA(nprobe).L_V1_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).L_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(PSTH_MUA_baseline(nprobe).R_V1_DOWN,'omitnan');
LCI = prctile(PSTH_MUA_baseline(nprobe).R_V1_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA_baseline(nprobe).R_V1_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
% xline(0,'r')
title('V1 MUA during Right DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
% nprobe = 2;
all_DOWN_no = length(probability(nprobe).DOWN_all_index);
time_wondows = [-1 1];
time_bin = 0.01;
x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
y = mean(PSTH_MUA(nprobe).R_HPC_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).R_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
xline(0,'r',LineWidth=1)
y = mean(PSTH_MUA(nprobe).L_HPC_DOWN,'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA(nprobe).L_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(PSTH_MUA_baseline(nprobe).R_HPC_DOWN,'omitnan');
LCI = prctile(PSTH_MUA_baseline(nprobe).R_HPC_DOWN_bootstrap,2.5);
UCI = prctile(PSTH_MUA_baseline(nprobe).R_HPC_DOWN_bootstrap,97.5);

PLOT = plot(x,y,'k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'ipsi','contra','Shuffled'},'Box','off')
% xline(0,'r')
title('HPC MUA during Right DOWN')
xlabel('Time relative to UP-DOWN transition (s)')
ylabel('MUA activity (z)')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)




%% Contra-UP probability


nprobe = 1;
event_averaging_scale = 10;
probability_UP = probability_SO_SO_contralateral_whole(1).UP_UP; 

nexttile
[~,sorted_index] = sort(event_info(nprobe).UP_duration);
% imagesc(movmean(50*movmean(R_ripples(sorted_index,:),50,1,'omitnan'),3,2,'omitnan'))
imagesc(event_averaging_scale*movmean(probability_UP(sorted_index,:),event_averaging_scale,1,'omitnan'))
xticks([1.5 25.5 50.5 75.5 100.5])
% xticklabels([PSTH_MUA(nprobe).timebins([1 50 100 150 200])+mean(diff(PSTH_MUA(nprobe).timebins)/2)])
xticklabels([-1 -0.5 0 0.5 1])
xline(50.5,'r',LineWidth=1)
clim([0 1])
colorbar
colormap(flipud(gray))
xlabel('Time relative to DOWN-UP transition (s)')
ylabel('Event sorted by UP duration')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
title('Contralateral UP during left UP')