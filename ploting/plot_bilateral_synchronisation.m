function plot_bilateral_synchronisation(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)

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


%% Basic UP/DOWN, spindle and ripple properties

% Duration of UP and DOWN
binEdges = -2:0.1:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

angleEdges =linspace(-pi, pi, 18+1);
angelCentre = angleEdges(1:end-1) + diff(angleEdges)/2;


UP_duration_L = [];
DOWN_duration_L = [];
UP_duration_R = [];
DOWN_duration_R = [];

UP_phase_L = [];
DOWN_phase_L = [];
UP_phase_R = [];
DOWN_phase_R = [];

delta_peaks_L = [];
delta_peaks_R = [];
delta_speed_L = [];
delta_speed_R = [];

ripples_duration_L = [];
ripples_duration_R = [];
ripple_amplitude_L = [];
ripple_amplitude_R = [];
ripple_speed_L = [];
ripple_speed_R = [];

spindle_duration_L = [];
spindle_duration_R = [];
spindle_amplitude_L = [];
spindle_amplitude_R = [];

spindle_SO_phase_L = nan(max(ripples_all(1).session_count), length(angelCentre));
spindle_SO_phase_R = nan(max(ripples_all(1).session_count), length(angelCentre));
ripple_spindle_phase_L = [];
ripple_spindle_phase_R = [];

for nsession = 1:max(ripples_all(1).session_count)

    % UP DOWN
    duration_dist = log10(slow_waves_all(1).UP_ints(slow_waves_all(1).UP_session_count == nsession ,2) - slow_waves_all(1).UP_ints(slow_waves_all(1).UP_session_count == nsession ,1));
    UP_duration_L(nsession,:) = histcounts(duration_dist,binEdges)/length(duration_dist);
    
    duration_dist = log10(slow_waves_all(2).UP_ints(slow_waves_all(2).UP_session_count == nsession ,2) - slow_waves_all(2).UP_ints(slow_waves_all(2).UP_session_count == nsession ,1));
    UP_duration_R(nsession,:) = histcounts(duration_dist,binEdges)/length(duration_dist);
    
    duration_dist = log10(slow_waves_all(1).DOWN_ints(slow_waves_all(1).UP_session_count == nsession ,2) - slow_waves_all(1).DOWN_ints(slow_waves_all(1).UP_session_count == nsession ,1));
    DOWN_duration_L(nsession,:) = histcounts(duration_dist,binEdges)/length(duration_dist);
    
    duration_dist = log10(slow_waves_all(2).DOWN_ints(slow_waves_all(2).UP_session_count == nsession ,2) - slow_waves_all(2).DOWN_ints(slow_waves_all(2).UP_session_count == nsession ,1));
    DOWN_duration_R(nsession,:) = histcounts(duration_dist,binEdges)/length(duration_dist);

    UP_phase_L(nsession,:) = histcounts(slow_waves_all(1).mean_phase_UP{nsession}(3,:),angleEdges);
    DOWN_phase_L(nsession,:) = histcounts(slow_waves_all(1).mean_phase_DOWN{nsession}(3,:),angleEdges);
    UP_phase_R(nsession,:) = histcounts(slow_waves_all(2).mean_phase_UP{nsession}(6,:),angleEdges);
    DOWN_phase_R(nsession,:) = histcounts(slow_waves_all(2).mean_phase_DOWN{nsession}(6,:),angleEdges);


    delta_peaks_L(nsession,:) = histcounts(slow_waves_all(1).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,1),:),-1:0.2:5);
    delta_peaks_R(nsession,:) = histcounts(slow_waves_all(2).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession,2),:),-1:0.2:5);

    delta_speed_L(nsession,:)  = histcounts(slow_waves_all(1).cortex_speed_UD(1,slow_waves_all(1).DOWN_session_count == nsession),-100:2:100);
    delta_speed_R(nsession,:)  = histcounts(slow_waves_all(2).cortex_speed_UD(2,slow_waves_all(2).DOWN_session_count == nsession),-100:2:100);
    % 
    % delta_latency_speed_L(nsession,:) = histcounts(mean(diff(0.00025*slow_waves_all(1).shank_id{nsession}(slow_waves_all(1).probe_hemisphere{nsession} == 1))'./diff(slow_waves_all(1).DOWN_peaktimes{nsession}(slow_waves_all(1).probe_hemisphere{nsession}==1,:))),-0.1:0.005:0.1);
    % delta_latency_speed_R(nsession,:)  = histcounts(-1*mean(diff(0.00025*slow_waves_all(2).shank_id{nsession}(slow_waves_all(2).probe_hemisphere{nsession} == 2))'./diff(slow_waves_all(2).DOWN_peaktimes{nsession}(slow_waves_all(2).probe_hemisphere{nsession}==2,:))),-0.1:0.005:0.1);
    % 
    % Ripples
    ripples_duration_L(nsession,:) = histcounts(ripples_all(1).offset(ripples_all(1).session_count == nsession) - ripples_all(1).onset(ripples_all(1).session_count == nsession),0:0.01:0.2);
    ripples_duration_R(nsession,:) = histcounts(ripples_all(2).offset(ripples_all(2).session_count == nsession) - ripples_all(2).onset(ripples_all(2).session_count == nsession),0:0.01:0.2);
    ripple_amplitude_L(nsession,:) = histcounts(ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession),5:0.5:25);
    ripple_amplitude_R(nsession,:) = histcounts(ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession),5:0.5:25);
    ripple_speed_L(nsession,:)  = histcounts(ripples_all(1).HPC_speed(1,ripples_all(1).session_count == nsession),-500:10:500);
    ripple_speed_R(nsession,:)  = histcounts(ripples_all(2).HPC_speed(2,ripples_all(2).session_count == nsession),-500:10:500);
    ripple_SO_phase_L(nsession,:) = histcounts(ripples_all(1).SO_phase_ripple_peaktime{nsession}(3,:),angleEdges);
    ripple_SO_phase_R(nsession,:) = histcounts(ripples_all(2).SO_phase_ripple_peaktime{nsession}(6,:),angleEdges);
    ripple_spindle_phase_L(nsession,:) = histcounts(ripples_all(1).spindle_phase_ripple_peaktime{nsession}(3,:),angleEdges);
    ripple_spindle_phase_R(nsession,:) = histcounts(ripples_all(2).spindle_phase_ripple_peaktime{nsession}(6,:),angleEdges);

    % Spindles
    spindle_duration_L(nsession,:) = histcounts(spindles_all(1).offset(spindles_all(1).session_count == nsession) - spindles_all(1).onset(spindles_all(1).session_count == nsession),0:0.05:2);
    spindle_duration_R(nsession,:) = histcounts(spindles_all(2).offset(spindles_all(2).session_count == nsession) - spindles_all(2).onset(spindles_all(2).session_count == nsession),0:0.05:2);
    spindle_amplitude_L(nsession,:) = histcounts(spindles_all(1).peak_zscore(spindles_all(1).session_count == nsession),1:0.2:15);
    spindle_amplitude_R(nsession,:) = histcounts(spindles_all(2).peak_zscore(spindles_all(2).session_count == nsession),1:0.2:15);
    if ~isempty(spindles_all(1).SO_phase_spindle_onset{nsession})
        spindle_SO_phase_L(nsession,:) = histcounts(spindles_all(1).SO_phase_spindle_onset{nsession}(3,:),angleEdges);
    end

    if ~isempty(spindles_all(1).SO_phase_spindle_onset{nsession})
        spindle_SO_phase_R(nsession,:) = histcounts(spindles_all(2).SO_phase_spindle_onset{nsession}(6,:),angleEdges);
    end
end


% figure
% for nsession = 1:22
%     nexttile
%     plot(angelCentre,(UP_phase_L(nsession,:)./sum(UP_phase_L(nsession,:)))')
%     xticks([-pi -1/2*pi 0 1/2*pi pi])
% xticklabels({'-π','-π/2','0','π/2','π'})
% end
% 
% figure
% for nsession = 1:22
%     nexttile
%     plot(angelCentre,(UP_phase_R(nsession,:)./sum(UP_phase_R(nsession,:)))')
%     xticks([-pi -1/2*pi 0 1/2*pi pi])
% xticklabels({'-π','-π/2','0','π/2','π'})
% end
% 

%%%%%%%%%%%%%%% UP DOWN
fig = figure('Color','w');
fig.Position = [850 58 1020 920];
fig.Name = 'UP DOWN Left and Right Basic properties';
sgtitle('Left and Right V1 UP/DOWN')
colour_lines = [44,123,182;215,25,28]/256;


binEdges = -2:0.1:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;
nexttile
x = binCentre;
y = mean(UP_duration_L./sum(UP_duration_L,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(UP_duration_L./sum(UP_duration_L,2));
UCI = y + std(UP_duration_L./sum(UP_duration_L,2));

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(UP_duration_R./sum(UP_duration_R,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(UP_duration_R./sum(UP_duration_R,2));
UCI = y + std(UP_duration_R./sum(UP_duration_R,2));

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-2 -1 -0.3 0])
xticklabels([0.01 0.1 0.5 1])

xlabel('UP duration (log10 sec)')
ylabel('Proportion of UP events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
x = binCentre;
y = mean(DOWN_duration_L./sum(DOWN_duration_L,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(DOWN_duration_L./sum(DOWN_duration_L,2));
UCI = y + std(DOWN_duration_L./sum(DOWN_duration_L,2));

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(DOWN_duration_R./sum(DOWN_duration_R,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(DOWN_duration_R./sum(DOWN_duration_R,2));
UCI = y + std(DOWN_duration_R./sum(DOWN_duration_R,2));

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-2 -1 -0.3 0])
xticklabels([0.01 0.1 0.5 1])

xlabel('DOWN duration (log10 sec)')
ylabel('Proportion of DOWN events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


% figure
% for nsession = 1:22
%     nexttile
%     plot(angelCentre,(UP_phase_L(nsession,:)./sum(UP_phase_L(nsession,:)))')
%     xticks([-pi -1/2*pi 0 1/2*pi pi])
% xticklabels({'-π','-π/2','0','π/2','π'})
% end
% 
% figure
% for nsession = 1:22
%     nexttile
%     plot(angelCentre,(UP_phase_R(nsession,:)./sum(UP_phase_R(nsession,:)))')
%     xticks([-pi -1/2*pi 0 1/2*pi pi])
% xticklabels({'-π','-π/2','0','π/2','π'})
% end


% 
% figure;
% y = mean(UP_phase_L./sum(UP_phase_L,2));
% LCI = y - std(UP_phase_L./sum(UP_phase_L,2))
% UCI = y + std(UP_phase_L./sum(UP_phase_L,2))
% 
% polarplot([angelCentre angelCentre(1)], [y y(1)],'b', 'LineWidth', 2);
% hold on;
% y = mean(UP_phase_R./sum(UP_phase_R,2));
% %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
% LCI = y - std(UP_phase_R./sum(UP_phase_R,2))
% UCI = y + std(UP_phase_R./sum(UP_phase_R,2))
% polarplot([angelCentre angelCentre(1)], [y y(1)],'r', 'LineWidth', 2);
% 
% 
% 
% 
% figure;
% y = mean(DOWN_phase_L./sum(DOWN_phase_L,2));
% LCI = y - std(DOWN_phase_L./sum(DOWN_phase_L,2))
% UCI = y + std(DOWN_phase_L./sum(DOWN_phase_L,2))
% 
% polarplot([angelCentre angelCentre(1)], [y y(1)],'b', 'LineWidth', 2);
% hold on;
% y = mean(DOWN_phase_R./sum(DOWN_phase_R,2));
% %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
% LCI = y - std(DOWN_phase_R./sum(DOWN_phase_R,2))
% UCI = y + std(DOWN_phase_R./sum(DOWN_phase_R,2))
% polarplot([angelCentre angelCentre(1)], [y y(1)],'r', 'LineWidth', 2);


nexttile

x = angelCentre;
y = mean(UP_phase_L./sum(UP_phase_L,2));
LCI = y - std(UP_phase_L./sum(UP_phase_L,2));
UCI = y + std(UP_phase_L./sum(UP_phase_L,2));

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(UP_phase_R./sum(UP_phase_R,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(UP_phase_R./sum(UP_phase_R,2));
UCI = y + std(UP_phase_R./sum(UP_phase_R,2));

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-pi -1/2*pi 0 1/2*pi pi])
xticklabels({'-π','-π/2','0','π/2','π'})
xlabel('Mean Phase of UP')
ylabel('Proportion of UP events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
x = angelCentre;
y = mean(DOWN_phase_L./sum(DOWN_phase_L,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(DOWN_phase_L./sum(DOWN_phase_L,2));
UCI = y + std(DOWN_phase_L./sum(DOWN_phase_L,2));


PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(DOWN_phase_R./sum(DOWN_phase_R,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(DOWN_phase_R./sum(DOWN_phase_R,2));
UCI = y + std(DOWN_phase_R./sum(DOWN_phase_R,2));

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-pi -1/2*pi 0 1/2*pi pi])
xticklabels({'-π','-π/2','0','π/2','π'})
xlabel('Mean phase of DOWN')
ylabel('Proportion of DOWN events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



binEdges = -1:0.2:5;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

nexttile
x = binCentre;
y = mean(delta_peaks_L./sum(delta_peaks_L,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(delta_peaks_L./sum(delta_peaks_L,2));
UCI = y + std(delta_peaks_L./sum(delta_peaks_L,2));

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(delta_peaks_R./sum(delta_peaks_R,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(delta_peaks_R./sum(delta_peaks_R,2));
UCI = y + std(delta_peaks_R./sum(delta_peaks_R,2));

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xlabel('Peak Delta power (z)')
ylabel('Proportion of DOWN/UP events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
binEdges = -100:2:100;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;
x = binCentre;
y = mean(delta_speed_L./sum(delta_speed_L,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(delta_speed_L./sum(delta_speed_L,2));
UCI = y + std(delta_speed_L./sum(delta_speed_L,2));

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(delta_speed_R./sum(delta_speed_R,2));
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(delta_speed_R./sum(delta_speed_R,2));
UCI = y + std(delta_speed_R./sum(delta_speed_R,2));

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xlim([-30 30])
xlabel('Delta wave propogation speed (mm/s)')
ylabel('Proportion of DOWN/UP events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%%%%%%%%%%%%%%% Spindles

fig = figure('Color','w');
fig.Position = [800 300 630 530];
fig.Name = 'Spindles Left and Right Basic properties';

sgtitle('Left and Right V1 Spindles')
colour_lines = [44,123,182;215,25,28]/256;

nexttile
x = 0+0.05/2:0.05:2-0.05/2;
y = mean(spindle_duration_L./sum(spindle_duration_L,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(spindle_duration_L./sum(spindle_duration_L,2,'omitnan'),'omitnan')
UCI = y + std(spindle_duration_L./sum(spindle_duration_L,2,'omitnan'),'omitnan')

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(spindle_duration_R./sum(spindle_duration_R,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(spindle_duration_R./sum(spindle_duration_R,2,'omitnan'),'omitnan')
UCI = y + std(spindle_duration_R./sum(spindle_duration_R,2,'omitnan'),'omitnan')
PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xticks([0 ])
% xticklabels([0.01 0.1 0.5 1])
xlim([0.2 1.4])
xlabel('spindle duration (s)')
ylabel('Proportion of spindles events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
x = 1+0.2/2:0.2:15-0.2/2;
y = mean(spindle_amplitude_L./sum(spindle_amplitude_L,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(spindle_amplitude_L./sum(spindle_amplitude_L,2,'omitnan'),'omitnan');
UCI = y + std(spindle_amplitude_L./sum(spindle_amplitude_L,2,'omitnan'),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(spindle_amplitude_R./sum(spindle_amplitude_R,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(spindle_amplitude_R./sum(spindle_amplitude_R,2,'omitnan'),'omitnan');
UCI = y + std(spindle_amplitude_R./sum(spindle_amplitude_R,2,'omitnan'),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xlim([1 10])
xlabel('Spindle peak amplitude (z)')
ylabel('Proportion of spindle events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
x = angelCentre;
y = mean(spindle_SO_phase_L./sum(spindle_SO_phase_L,2),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(spindle_SO_phase_L./sum(spindle_SO_phase_L,2),'omitnan');
UCI = y + std(spindle_SO_phase_L./sum(spindle_SO_phase_L,2),'omitnan');
PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y = mean(spindle_SO_phase_R./sum(spindle_SO_phase_R,2),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(spindle_SO_phase_R./sum(spindle_SO_phase_R,2),'omitnan');
UCI = y + std(spindle_SO_phase_R./sum(spindle_SO_phase_R,2),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-pi -1/2*pi 0 1/2*pi pi])
xticklabels({'-π','-π/2','0','π/2','π'})
xlabel('SO phase at spindle onset')
ylabel('Proportion of spindle events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


%%%%%%%%%%%%%% Ripples
fig = figure('Color','w');
fig.Position = [800 300 630 530];
fig.Name = 'Ripples Left and Right Basic properties';
sgtitle('Left and Right Ripples')
colour_lines = [44,123,182;215,25,28]/256;


nexttile
x = 0+0.01/2:0.01:0.2-0.01/2;
y = mean(ripples_duration_L./sum(ripples_duration_L,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripples_duration_L./sum(ripples_duration_L,2,'omitnan'),'omitnan')
UCI = y + std(ripples_duration_L./sum(ripples_duration_L,2,'omitnan'),'omitnan')

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(ripples_duration_R./sum(ripples_duration_R,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripples_duration_R./sum(ripples_duration_R,2,'omitnan'),'omitnan')
UCI = y + std(ripples_duration_R./sum(ripples_duration_R,2,'omitnan'),'omitnan')
PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xticks([0 ])
% xticklabels([0.01 0.1 0.5 1])
% xlim([0.2 1.4])
xlabel('ripple duration (s)')
ylabel('Proportion of ripple events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
x = 5+0.5/2:0.5:25-0.5/2;
y = mean(ripple_amplitude_L./sum(ripple_amplitude_L,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripple_amplitude_L./sum(ripple_amplitude_L,2,'omitnan'),'omitnan')
UCI = y + std(ripple_amplitude_L./sum(ripple_amplitude_L,2,'omitnan'),'omitnan')

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(ripple_amplitude_R./sum(ripple_amplitude_R,2,'omitnan'),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripple_amplitude_R./sum(ripple_amplitude_R,2,'omitnan'),'omitnan')
UCI = y + std(ripple_amplitude_R./sum(ripple_amplitude_R,2,'omitnan'),'omitnan')
PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% xlim([1 10])
xlabel('Ripple peak amplitude')
ylabel('Proportion of ripple events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
x = angelCentre;
y = mean(ripple_SO_phase_L./sum(ripple_SO_phase_L,2),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripple_SO_phase_L./sum(ripple_SO_phase_L,2),'omitnan');
UCI = y + std(ripple_SO_phase_L./sum(ripple_SO_phase_L,2),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(ripple_SO_phase_R./sum(ripple_SO_phase_R,2),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripple_SO_phase_R./sum(ripple_SO_phase_R,2),'omitnan');
UCI = y + std(ripple_SO_phase_R./sum(ripple_SO_phase_R,2),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-pi -1/2*pi 0 1/2*pi pi])
xticklabels({'-π','-π/2','0','π/2','π'})
xlabel('SO phase at ripple peak')
ylabel('Proportion of ripple events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
x = angelCentre;
y = mean(ripple_spindle_phase_L./sum(ripple_spindle_phase_L,2),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripple_spindle_phase_L./sum(ripple_spindle_phase_L,2),'omitnan');
UCI = y + std(ripple_spindle_phase_L./sum(ripple_spindle_phase_L,2),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(1,:),'LineWidth',2);hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


y = mean(ripple_spindle_phase_R./sum(ripple_spindle_phase_R,2),'omitnan');
%     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
LCI = y - std(ripple_spindle_phase_R./sum(ripple_spindle_phase_R,2),'omitnan');
UCI = y + std(ripple_spindle_phase_R./sum(ripple_spindle_phase_R,2),'omitnan');

PLOT = plot(x,y,'Color',colour_lines(2,:),'LineWidth',2);hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
xticks([-pi -1/2*pi 0 1/2*pi pi])
xticklabels({'-π','-π/2','0','π/2','π'})
xlabel('Spindle phase at ripple peak')
ylabel('Proportion of ripple events')
legend(ERROR_SHADE(1:2),{'Left','Right'},'Location','northwest','Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

% if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
%     mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
% end
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])

%% Ipsilateral vs contralateral correlation (UP DOWN)
probability = probability_psth_whole;

for nprobe = 1:2
    ipsi_amp_corr{nprobe} = [];
    contra_amp_corr{nprobe} = [];

    ipsi_plv{nprobe} = [];
    contra_plv{nprobe} = [];

    ipsi_lag{nprobe} = [];
    contra_lag{nprobe} = [];

    ipsi_amp_corr_bootstrap{nprobe} = [];
    contra_amp_corr_bootstrap{nprobe} = [];

    ipsi_plv_bootstrap{nprobe} = [];
    contra_plv_bootstrap{nprobe} = [];

    ipsi_lag_bootstrap{nprobe} = [];
    contra_lag_bootstrap{nprobe} = [];

    mprobe = abs(nprobe-3);
    % ipsi_amp_corr_shuffle{nprobe} = [];
    % contra_amp_corr_shuffle{nprobe} = [];
    % 
    % ipsi_plv_shuffle{nprobe} = [];
    % contra_plv_shuffle{nprobe} = [];

    for nsession = 1:max(ripples_all(nprobe).session_count)
        % ipsi_amp_corr_bootstrap{nprobe}{nsession} = [];
        % contra_amp_corr_bootstrap{nprobe}{nsession} = [];
        % 
        % ipsi_plv_bootstrap{nprobe}{nsession} = [];
        % contra_plv_bootstrap{nprobe}{nsession} = [];

        % Find UP events with less than 2 seconds and followed by a
        % DOWN
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession)),probability(nprobe).DOWN_all_index);

        ipsi_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        contra_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
            % mean_corr_ipsi(nshank) = mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
            % mean_corr_contra(nshank)= mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);

        if length(ipsi_shank)==1
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} (squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
        else
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        end

        if length(contra_shank)==1
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
            contra_plv{nprobe} = [contra_plv{nprobe} (squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
            contra_lag{nprobe} = [contra_lag{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
        else
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_plv{nprobe} = [contra_plv{nprobe} mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_lag{nprobe} = [contra_lag{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        end

        temp1 = [];
        temp2 = [];
        temp3 = [];
        temp4 = [];
        temp5 = [];
        temp6 = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,ia,size(ia,1));
            if length(ipsi_shank)==1
                temp1(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= (squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp5(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            else
                temp1(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp5(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            end

            if length(contra_shank)==1
                temp2(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp4(iBoot,:)=  (squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp6(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            else
                temp2(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp4(iBoot,:)=  mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp6(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_lag_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            end

        end
        % for iBoot = 1:500
        %     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        %     probe_hemisphere = datasample(s,slow_waves_all(nprobe).probe_hemisphere{nsession},size(slow_waves_all(nprobe).probe_hemisphere{nsession},2),'Replace',false);
        % 
        %     ipsi_shank = find(probe_hemisphere == nprobe);
        %     ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        %     contra_shank = find(probe_hemisphere == mprobe);
        % 
        %     ipsi_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        % 
        %     ipsi_plv_shuffle{nprobe}{nsession}(iBoot,:) =mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_plv_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        % 
        % end

        ipsi_amp_corr_bootstrap{nprobe} = [ipsi_amp_corr_bootstrap{nprobe} temp1];
        contra_amp_corr_bootstrap{nprobe} = [contra_amp_corr_bootstrap{nprobe} temp2];
        ipsi_plv_bootstrap{nprobe} = [ipsi_plv_bootstrap{nprobe} temp3];
        contra_plv_bootstrap{nprobe} = [contra_plv_bootstrap{nprobe} temp4];
        ipsi_lag_bootstrap{nprobe} = [ipsi_lag_bootstrap{nprobe} temp5];
        contra_lag_bootstrap{nprobe} = [contra_lag_bootstrap{nprobe} temp6];
    end
end

[p,h,stats] = signrank(ipsi_plv{1},contra_plv{1},'tail','right');
[p,h,stats] = signrank(ipsi_plv{2},contra_plv{2},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{1},contra_amp_corr{1},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{2},contra_amp_corr{2},'tail','right');
[p,h,stats] = signrank(ipsi_lag{1},contra_lag{1});
[p,h,stats] = signrank(ipsi_lag{2},contra_lag{2});

nprobe = 1;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(slow_waves_all(nprobe).UP_session_count(probability(nprobe).UP_all_index),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataL = table(ipsi_plv{1}',contra_plv{1}',ipsi_amp_corr{1}',contra_amp_corr{1}',ipsi_lag{1}',contra_lag{1}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});

nprobe = 2;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(slow_waves_all(nprobe).UP_session_count(probability(nprobe).UP_all_index),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataR = table(ipsi_plv{2}',contra_plv{2}',ipsi_amp_corr{2}',contra_amp_corr{2}',ipsi_lag{2}',contra_lag{2}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});


nfig = figure('Color','w','Name','Left Right V1 UP-DOWN transition ipsilateral-contralateral difference')
nfig.Position = [103 111 1650 840];
orient(nfig,'landscape')

colour_lines = [44,123,182;215,25,28]/256;


%%%%% Left
data = dataL;
nprobe =1;

nexttile
formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Left UP-DOWN transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi delta power correlation')
ylabel('contra delta power correlation')
% set(gca,'FontSize',14)
title('Left UP-DOWN transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.15 0.15])
ylim([-0.15 0.15])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Left UP-DOWN transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('delta power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.15:0.01:0.15;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('lag')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


%%%%% Right
data = dataR;
nprobe =2;

nexttile
colour_lines = [44,123,182;215,25,28]/256;

formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Right UP-DOWN transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi delta power correlation')
ylabel('contra delta power correlation')
% set(gca,'FontSize',14)
title('Right UP-DOWN transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.15 0.15])
ylim([-0.15 0.15])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Right UP-DOWN transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('delta power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.15:0.01:0.15;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('lag')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% colour_lines = [44,123,182;215,25,28]/256;
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
%     x = 1:length(ipsi_plv{nprobe});
%     y = mean(cumsum(ipsi_plv_bootstrap{nprobe}(:,sorted_index) - contra_plv_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('UP-DOWN transition')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in phase-locking value'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
%     x = 1:length(ipsi_amp_corr{nprobe});
%     y = mean(cumsum(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('UP-DOWN transition')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in delta power correlation'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')


%% Ipsilateral vs contralateral correlation (DOWN UP)
probability = probability_psth_whole;
ipsi_lag=[];contra_lag=[];
for nprobe = 1:2
    ipsi_amp_corr{nprobe} = [];
    contra_amp_corr{nprobe} = [];

    ipsi_plv{nprobe} = [];
    contra_plv{nprobe} = [];

    ipsi_lag{nprobe} = [];
    contra_lag{nprobe} = [];

    ipsi_amp_corr_bootstrap{nprobe} = [];
    contra_amp_corr_bootstrap{nprobe} = [];

    ipsi_plv_bootstrap{nprobe} = [];
    contra_plv_bootstrap{nprobe} = [];

    ipsi_lag_bootstrap{nprobe} = [];
    contra_lag_bootstrap{nprobe} = [];

    mprobe = abs(nprobe-3);
    % ipsi_amp_corr_shuffle{nprobe} = [];
    % contra_amp_corr_shuffle{nprobe} = [];
    % 
    % ipsi_plv_shuffle{nprobe} = [];
    % contra_plv_shuffle{nprobe} = [];

    for nsession = 1:max(ripples_all(nprobe).session_count)
        % ipsi_amp_corr_bootstrap{nprobe}{nsession} = [];
        % contra_amp_corr_bootstrap{nprobe}{nsession} = [];
        % 
        % ipsi_plv_bootstrap{nprobe}{nsession} = [];
        % contra_plv_bootstrap{nprobe}{nsession} = [];

        % Find UP events with less than 2 seconds and followed by a
        % DOWN
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession)),probability(nprobe).UP_all_index);

        ipsi_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        contra_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
            % mean_corr_ipsi(nshank) = mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
            % mean_corr_contra(nshank)= mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);
        
        if length(ipsi_shank)==1
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} (squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
        else
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        end

        if length(contra_shank)==1
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
            contra_plv{nprobe} = [contra_plv{nprobe} (squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
            contra_lag{nprobe} = [contra_lag{nprobe} (squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
        else
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_plv{nprobe} = [contra_plv{nprobe} mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_lag{nprobe} = [contra_lag{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        end

        temp1 = [];
        temp2 = [];
        temp3 = [];
        temp4 = [];
        temp5 = [];
        temp6 = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,ia,size(ia,1));
            if length(ipsi_shank)==1
                temp1(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= (squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp5(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            else
                temp1(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp5(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            end

            if length(contra_shank)==1
                temp2(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp4(iBoot,:)=  (squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp6(iBoot,:)= (squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            else
                temp2(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp4(iBoot,:)=  mean(squeeze(slow_waves_all(nprobe).plv_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp6(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_lag_DU{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            end

        end
        % for iBoot = 1:500
        %     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        %     probe_hemisphere = datasample(s,slow_waves_all(nprobe).probe_hemisphere{nsession},size(slow_waves_all(nprobe).probe_hemisphere{nsession},2),'Replace',false);
        % 
        %     ipsi_shank = find(probe_hemisphere == nprobe);
        %     ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        %     contra_shank = find(probe_hemisphere == mprobe);
        % 
        %     ipsi_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        % 
        %     ipsi_plv_shuffle{nprobe}{nsession}(iBoot,:) =mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_plv_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        % 
        % end

        ipsi_amp_corr_bootstrap{nprobe} = [ipsi_amp_corr_bootstrap{nprobe} temp1];
        contra_amp_corr_bootstrap{nprobe} = [contra_amp_corr_bootstrap{nprobe} temp2];
        ipsi_plv_bootstrap{nprobe} = [ipsi_plv_bootstrap{nprobe} temp3];
        contra_plv_bootstrap{nprobe} = [contra_plv_bootstrap{nprobe} temp4];
        ipsi_lag_bootstrap{nprobe} = [ipsi_lag_bootstrap{nprobe} temp5];
        contra_lag_bootstrap{nprobe} = [contra_lag_bootstrap{nprobe} temp6];
    end
end
[p,h,stats] = signrank(ipsi_plv{1},contra_plv{1},'tail','right');
[p,h,stats] = signrank(ipsi_plv{2},contra_plv{2},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{1},contra_amp_corr{1},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{2},contra_amp_corr{2},'tail','right');
[p,h,stats] = signrank(ipsi_lag{1},contra_lag{1});
[p,h,stats] = signrank(ipsi_lag{2},contra_lag{2});

nprobe = 1;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(slow_waves_all(nprobe).UP_session_count(probability(nprobe).UP_all_index),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataL = table(ipsi_plv{1}',contra_plv{1}',ipsi_amp_corr{1}',contra_amp_corr{1}',ipsi_lag{1}',contra_lag{1}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});

nprobe = 2;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(slow_waves_all(nprobe).UP_session_count(probability(nprobe).UP_all_index),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataR = table(ipsi_plv{2}',contra_plv{2}',ipsi_amp_corr{2}',contra_amp_corr{2}',ipsi_lag{2}',contra_lag{2}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});


nfig = figure('Color','w','Name','Left Right V1 DOWN-UP transition ipsilateral-contralateral difference')
nfig.Position = [103 111 1650 840];
orient(nfig,'landscape')

colour_lines = [44,123,182;215,25,28]/256;


%%%%% Left
data = dataL;
nprobe =1;

nexttile
formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Left DOWN-UP transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi delta power correlation')
ylabel('contra delta power correlation')
% set(gca,'FontSize',14)
title('Left DOWN-UP transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.15 0.15])
ylim([-0.15 0.15])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Left DOWN-UP transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('delta power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.15:0.01:0.15;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('lag')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%%%%% Right
data = dataR;
nprobe =2;

nexttile
colour_lines = [44,123,182;215,25,28]/256;

formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Right DOWN-UP transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi delta power correlation')
ylabel('contra delta power correlation')
% set(gca,'FontSize',14)
title('Right DOWN-UP transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.15 0.15])
ylim([-0.15 0.15])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Left DOWN-UP transition');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.15:0.01:0.15;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('lag')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


% nexttile
% colour_lines = [0,90,50;74,20,134]/256; % Green Purple
% binEdges = -1:0.01:1;
% binCentre = binEdges(1:end-1) + diff(binEdges)/2;
% 
% ipsi_histcounts=[];
% contra_histcounts=[];
% 
% for iBoot = 1:1000
%     ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
%     contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
% end
% 
% x = binCentre;
% y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
% LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
% UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
% plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
% ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');
% 
% y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
% LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
% UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
% plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
% ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
% legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
% xlabel('delta power correlation')
% ylabel('cumulative proportion')
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
%     x = 1:length(ipsi_plv{nprobe});
%     y = mean(cumsum(ipsi_plv_bootstrap{nprobe}(:,sorted_index) - contra_plv_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('DOWN-UP transition')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in phase-locking value'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
%     x = 1:length(ipsi_amp_corr{nprobe});
%     y = mean(cumsum(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('DOWN-UP transition')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in delta power correlation'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')

%% Ipsilateral vs contralateral correlation (Ripple)
for nprobe = 1:2
    ipsi_amp_corr{nprobe} = [];
    contra_amp_corr{nprobe} = [];

    ipsi_plv{nprobe} = [];
    contra_plv{nprobe} = [];

    ipsi_lag{nprobe} = [];
    contra_lag{nprobe} = [];

    ipsi_amp_corr_bootstrap{nprobe} = [];
    contra_amp_corr_bootstrap{nprobe} = [];

    ipsi_plv_bootstrap{nprobe} = [];
    contra_plv_bootstrap{nprobe} = [];

    ipsi_lag_bootstrap{nprobe} = [];
    contra_lag_bootstrap{nprobe} = [];

    mprobe = abs(nprobe-3);
    % ipsi_amp_corr_shuffle{nprobe} = [];
    % contra_amp_corr_shuffle{nprobe} = [];
    %
    % ipsi_plv_shuffle{nprobe} = [];
    % contra_plv_shuffle{nprobe} = [];

    for nsession = 1:max(ripples_all(nprobe).session_count)
        % ipsi_amp_corr_bootstrap{nprobe}{nsession} = [];
        % contra_amp_corr_bootstrap{nprobe}{nsession} = [];
        %
        % ipsi_plv_bootstrap{nprobe}{nsession} = [];
        % contra_plv_bootstrap{nprobe}{nsession} = [];

        % Find UP events with less than 2 seconds and followed by a
        % DOWN
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

        if length(ipsi_shank)==1
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} (squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
        else
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} mean(squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        end

        if length(contra_shank)==1
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];
            contra_plv{nprobe} = [contra_plv{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia))'];
            contra_lag{nprobe} = [contra_lag{nprobe} (squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia)))'];
        else
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_plv{nprobe} = [contra_plv{nprobe} mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_lag{nprobe} = [contra_lag{nprobe} mean(squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia)))];
        end

        % if length(contra_amp_corr{nprobe}) ~= length( ipsi_amp_corr{nprobe})
        % 
        %     nsession
        % end

        temp1 = [];
        temp2 = [];
        temp3 = [];
        temp4 = [];
        temp5 = [];
        temp6 = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,ia,size(ia,1));
            if length(ipsi_shank)==1
                temp1(iBoot,:)= squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id));
                temp3(iBoot,:)= squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id));
                temp5(iBoot,:)= (squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            else
                temp1(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp5(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            end

            if length(contra_shank)==1
                temp2(iBoot,:)= squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id));
                temp4(iBoot,:)= squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id));
                temp6(iBoot,:)= (squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id)));
            else
                temp2(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp4(iBoot,:)= mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp6(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_lag{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id)));
            end

        end
        % for iBoot = 1:500
        %     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        %     probe_hemisphere = datasample(s,slow_waves_all(nprobe).probe_hemisphere{nsession},size(slow_waves_all(nprobe).probe_hemisphere{nsession},2),'Replace',false);
        %
        %     ipsi_shank = find(probe_hemisphere == nprobe);
        %     ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        %     contra_shank = find(probe_hemisphere == mprobe);
        %
        %     ipsi_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        %
        %     ipsi_plv_shuffle{nprobe}{nsession}(iBoot,:) =mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_plv_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        %
        % end

        ipsi_amp_corr_bootstrap{nprobe} = [ipsi_amp_corr_bootstrap{nprobe} temp1];
        contra_amp_corr_bootstrap{nprobe} = [contra_amp_corr_bootstrap{nprobe} temp2];
        ipsi_plv_bootstrap{nprobe} = [ipsi_plv_bootstrap{nprobe} temp3];
        contra_plv_bootstrap{nprobe} = [contra_plv_bootstrap{nprobe} temp4];
        ipsi_lag_bootstrap{nprobe} = [ipsi_lag_bootstrap{nprobe} temp5];
        contra_lag_bootstrap{nprobe} = [contra_lag_bootstrap{nprobe} temp6];
    end
end



[p,h,stats] = signrank(ipsi_plv{1},contra_plv{1},'tail','right');
[p,h,stats] = signrank(ipsi_plv{2},contra_plv{2},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{1},contra_amp_corr{1},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{2},contra_amp_corr{2},'tail','right');

% probability = probability_ripples_SO_whole;
nprobe = 1;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(ripples_all(nprobe).session_count(ripples_all(nprobe).SWS_index == 1),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataL = table(ipsi_plv{1}',contra_plv{1}',ipsi_amp_corr{1}',contra_amp_corr{1}',ipsi_lag{1}',contra_lag{1}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});

nprobe = 2;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(ripples_all(nprobe).session_count(ripples_all(nprobe).SWS_index == 1),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataR = table(ipsi_plv{2}',contra_plv{2}',ipsi_amp_corr{2}',contra_amp_corr{2}',ipsi_lag{2}',contra_lag{2}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});



nfig = figure('Color','w','Name','Left Right ripples ipsilateral-contralateral difference')
nfig.Position = [103 111 1650 840];
orient(nfig,'landscape')

colour_lines = [44,123,182;215,25,28]/256;


%%%%% Left
data = dataL;
nprobe =1;

nexttile
formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Left ripple');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi ripple power correlation')
ylabel('contra ripple power correlation')
% set(gca,'FontSize',14)
title('Left ripple');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.1 0.1])
ylim([-0.1 0.1])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Left ripples');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('ripple power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.1:0.01:0.1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('ripple lag (s)')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


%%%%% Right
data = dataR;
nprobe =2;

nexttile
colour_lines = [44,123,182;215,25,28]/256;

formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Right ripple');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi ripple power correlation')
ylabel('contra ripple power correlation')
% set(gca,'FontSize',14)
title('Right ripple');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.1 0.1])
ylim([-0.1 0.1])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Left ripples');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('ripple power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.1:0.01:0.1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('ripple lag (s)')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% colour_lines = [44,123,182;215,25,28]/256;
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(ripples_all(nprobe).peak_zscore(ripples_all(nprobe).SWS_index == 1));
%     x = 1:length(ipsi_plv{nprobe});
%     y = mean(cumsum(ipsi_plv_bootstrap{nprobe}(:,sorted_index) - contra_plv_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left ripples','Right ripples'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('ripple events')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in phase-locking value'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(ripples_all(nprobe).peak_zscore(ripples_all(nprobe).SWS_index == 1));
%     x = 1:length(ipsi_amp_corr{nprobe});
%     y = mean(cumsum(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left ripples','Right ripples'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('ripple events')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in ripple power correlation'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')

%% Ipsilateral vs contralateral correlation (Spindles)
% HPC_ref_shank
for nprobe = 1:2
    ipsi_amp_corr{nprobe} = [];
    contra_amp_corr{nprobe} = [];

    ipsi_plv{nprobe} = [];
    contra_plv{nprobe} = [];

    ipsi_lag{nprobe} = [];
    contra_lag{nprobe} = [];

    ipsi_amp_corr_bootstrap{nprobe} = [];
    contra_amp_corr_bootstrap{nprobe} = [];

    ipsi_plv_bootstrap{nprobe} = [];
    contra_plv_bootstrap{nprobe} = [];

    ipsi_lag_bootstrap{nprobe} = [];
    contra_lag_bootstrap{nprobe} = [];

    mprobe = abs(nprobe-3);
    % ipsi_amp_corr_shuffle{nprobe} = [];
    % contra_amp_corr_shuffle{nprobe} = [];
    %
    % ipsi_plv_shuffle{nprobe} = [];
    % contra_plv_shuffle{nprobe} = [];

    for nsession = 1:max(spindles_all(nprobe).session_count)
        % ipsi_amp_corr_bootstrap{nprobe}{nsession} = [];
        % contra_amp_corr_bootstrap{nprobe}{nsession} = [];
        %
        % ipsi_plv_bootstrap{nprobe}{nsession} = [];
        % contra_plv_bootstrap{nprobe}{nsession} = [];

        % Find UP events with less than 2 seconds and followed by a
        % DOWN
        [C,ia,ib] = intersect(find(spindles_all(nprobe).session_count == sessions_to_process(nsession)),find(spindles_all(nprobe).session_count == sessions_to_process(nsession) & spindles_all(nprobe).SWS_index == 1));

        if isempty(ia)
            continue
        end

        ipsi_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        contra_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == mprobe);

        mean_corr_ipsi = [];mean_corr_contra=[];
        for nshank = 1:length(ipsi_shank)
            mean_corr_ipsi(nshank) = mean(squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank(nshank),ia)));
        end

        for nshank = 1:length(contra_shank)
            mean_corr_contra(nshank)= mean(squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank(nshank),ia)));
        end

        [~,id] = max(mean_corr_ipsi);
        ipsi_shank = ipsi_shank(id);
        [~,id] = max(mean_corr_contra);
        contra_shank = contra_shank(id);

        if length(ipsi_shank)==1
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} (squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))'];
        else
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} mean(squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} mean(squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_lag{nprobe} = [ipsi_lag{nprobe} mean(squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        end

        if length(contra_shank)==1
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];
            contra_plv{nprobe} = [contra_plv{nprobe} squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia))'];
            contra_lag{nprobe} = [contra_lag{nprobe} (squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))'];
        else
            contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} mean(squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_plv{nprobe} = [contra_plv{nprobe} mean(squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
            contra_lag{nprobe} = [contra_lag{nprobe} mean(squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];
        end

        % if length(contra_amp_corr{nprobe}) ~= length( ipsi_amp_corr{nprobe})
        %
        %     nsession
        % end

        temp1 = [];
        temp2 = [];
        temp3 = [];
        temp4 = [];
        temp5 = [];
        temp6 = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,ia,size(ia,1));
            if length(ipsi_shank)==1
                temp1(iBoot,:)= squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id));
                temp3(iBoot,:)= squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id));
                temp5(iBoot,:)= (squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            else
                temp1(iBoot,:)= mean(squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= mean(squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp5(iBoot,:)= mean(squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            end

            if length(contra_shank)==1
                temp2(iBoot,:)= squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id));
                temp4(iBoot,:)= squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id));
                temp6(iBoot,:)= (squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            else
                temp2(iBoot,:)= mean(squeeze(spindles_all(nprobe).xcorr_r{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp4(iBoot,:)= mean(squeeze(spindles_all(nprobe).plv{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
                temp6(iBoot,:)= mean(squeeze(spindles_all(nprobe).xcorr_lag{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            end

        end
        % for iBoot = 1:500
        %     s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        %     probe_hemisphere = datasample(s,slow_waves_all(nprobe).probe_hemisphere{nsession},size(slow_waves_all(nprobe).probe_hemisphere{nsession},2),'Replace',false);
        %
        %     ipsi_shank = find(probe_hemisphere == nprobe);
        %     ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        %     contra_shank = find(probe_hemisphere == mprobe);
        %
        %     ipsi_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_amp_corr_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        %
        %     ipsi_plv_shuffle{nprobe}{nsession}(iBoot,:) =mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,:)));
        %     contra_plv_shuffle{nprobe}{nsession}(iBoot,:) = mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,:)));
        %
        % end

        ipsi_amp_corr_bootstrap{nprobe} = [ipsi_amp_corr_bootstrap{nprobe} temp1];
        contra_amp_corr_bootstrap{nprobe} = [contra_amp_corr_bootstrap{nprobe} temp2];
        ipsi_plv_bootstrap{nprobe} = [ipsi_plv_bootstrap{nprobe} temp3];
        contra_plv_bootstrap{nprobe} = [contra_plv_bootstrap{nprobe} temp4];
        ipsi_lag_bootstrap{nprobe} = [ipsi_lag_bootstrap{nprobe} temp5];
        contra_lag_bootstrap{nprobe} = [contra_lag_bootstrap{nprobe} temp6];
    end
end



[p,h,stats] = signrank(ipsi_plv{1},contra_plv{1},'tail','right');
[p,h,stats] = signrank(ipsi_plv{2},contra_plv{2},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{1},contra_amp_corr{1},'tail','right');
[p,h,stats] = signrank(ipsi_amp_corr{2},contra_amp_corr{2},'tail','right');
[p,h,stats] = signrank(ipsi_lag{1},contra_lag{1});
[p,h,stats] = signrank(ipsi_lag{2},contra_lag{2});


% probability = probability_ripples_SO_whole;
nprobe = 1;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(spindles_all(nprobe).session_count(spindles_all(nprobe).SWS_index == 1),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataL = table(ipsi_plv{1}',contra_plv{1}',ipsi_amp_corr{1}',contra_amp_corr{1}',ipsi_lag{1}',contra_lag{1}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});

nprobe = 2;
subject_id = str2double(cellstr(slow_waves_all(nprobe).subject(spindles_all(nprobe).session_count(spindles_all(nprobe).SWS_index == 1),end-1:end)));
[~, ~, mappedIDs] = unique(subject_id);
dataR = table(ipsi_plv{2}',contra_plv{2}',ipsi_amp_corr{2}',contra_amp_corr{2}',ipsi_lag{2}',contra_lag{2}',mappedIDs,...
    'VariableNames',{'ipsi_plv','contra_plv','ipsi_amp_corr','contra_amp_corr','ipsi_lag','contra_lag','animal_label'});


nfig = figure('Color','w','Name','Left Right spindles ipsilateral-contralateral difference')
nfig.Position = [103 111 1422 840];
orient(nfig,'landscape')

colour_lines = [44,123,182;215,25,28]/256;


%%%%% Left
data = dataL;
nprobe =1;

nexttile
formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Left spindle');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi spindle power correlation')
ylabel('contra spindle power correlation')
% set(gca,'FontSize',14)
title('Left spindle');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end



nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(1,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.1 0.1])
ylim([-0.1 0.1])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Left spindles');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('spindle power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.1:0.01:0.1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('spindle lag (s)')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%%%%% Right
data = dataR;
nprobe =2;

nexttile
colour_lines = [44,123,182;215,25,28]/256;

formula = 'contra_plv~  ipsi_plv + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_plv,data.contra_plv,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([0 1],[0 1],'k')

mdl = fitlm(data.ipsi_plv,data.contra_plv);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_plv) max(data.ipsi_plv)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi phase locking')
ylabel('contra phase locking')
% set(gca,'FontSize',14)
title('Right spindle');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
formula = 'contra_amp_corr~  ipsi_amp_corr + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_amp_corr,data.contra_amp_corr,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-1 1],[-1 1],'k')

mdl = fitlm(data.ipsi_amp_corr,data.contra_amp_corr);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_amp_corr) max(data.ipsi_amp_corr)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi spindle power correlation')
ylabel('contra spindle power correlation')
% set(gca,'FontSize',14)
title('Right spindle');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end


nexttile
formula = 'contra_lag~  ipsi_lag + (1|animal_label)';
mdl2 = fitlme(data,formula);

hold on
% arrayfun(@(x) scatter(awake_rate(x),awake_theta(x),86,new_cls(x,:),'filled','o'),1:length(awake_theta))
scatter(data.ipsi_lag,data.contra_lag,10,colour_lines(2,:),'filled','o','MarkerFaceAlpha',0.05)
plot([-0.15 0.15],[-0.15 0.15],'k')
xlim([-0.1 0.1])
ylim([-0.1 0.1])
mdl = fitlm(data.ipsi_lag,data.contra_lag);
% [pval,F_stat,~] = coefTest(mdl);
% awake_rate_R2 = mdl.Rsquared.Adjusted;
x =[min(data.ipsi_lag) max(data.ipsi_lag)];
b = mdl.Coefficients.Estimate';
y_est = polyval(fliplr(b),x);
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('ipsi lag')
ylabel('contra lag')
% set(gca,'FontSize',14)
title('Right spindles');

f=get(gca,'Children');
% Mind that order is reversed
if mdl2.Coefficients.pValue(2) < 0.05
    plot(x,y_est,':','Color','m','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square

else
    plot(x,y_est,':','Color','k','LineWidth',3)
    text(gca,.5,0.1,sprintf('p = %.2d & R2 = %.3f',mdl2.Coefficients.pValue(2),mdl2.Rsquared.Adjusted),'Units','Normalized','FontName','Arial');
    axis square
end

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = 0:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_plv_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_plv_bootstrap{nprobe}(iBoot,:),binEdges);
end
x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('phase locking value')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -1:0.01:1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_amp_corr_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('spindle power correlation')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
colour_lines = [0,90,50;74,20,134]/256; % Green Purple
binEdges = -0.1:0.01:0.1;
binCentre = binEdges(1:end-1) + diff(binEdges)/2;

ipsi_histcounts=[];
contra_histcounts=[];

for iBoot = 1:1000
    ipsi_histcounts(iBoot,:) = histcounts(ipsi_lag_bootstrap{nprobe}(iBoot,:),binEdges);
    contra_histcounts(iBoot,:) = histcounts(contra_lag_bootstrap{nprobe}(iBoot,:),binEdges);
end

x = binCentre;
y  = cumsum(mean(ipsi_histcounts)/sum(mean(ipsi_histcounts)));
LCI = cumsum(prctile(ipsi_histcounts,2.5),2)/sum(prctile(ipsi_histcounts,2.5),2);
UCI = cumsum(prctile(ipsi_histcounts,97.5),2)/sum(prctile(ipsi_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(1,:)); hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');

y  = cumsum(mean(contra_histcounts)/sum(mean(contra_histcounts)));
LCI = cumsum(prctile(contra_histcounts,2.5),2)/sum(prctile(contra_histcounts,2.5),2);
UCI = cumsum(prctile(contra_histcounts,97.5),2)/sum(prctile(contra_histcounts,97.5),2);
plot(binCentre,y,'Color',colour_lines(2,:)); hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:2),{'ipsi','contra'},'box', 'off')
xlabel('spindle lag (s)')
ylabel('cumulative proportion')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

% 
% colour_lines = [44,123,182;215,25,28]/256;
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(spindles_all(nprobe).peak_zscore(spindles_all(nprobe).SWS_index == 1));
%     x = 1:length(ipsi_plv{nprobe});
%     y = mean(cumsum(ipsi_plv_bootstrap{nprobe}(:,sorted_index) - contra_plv_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left spindles','Right spindles'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('spindle events')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in phase-locking value'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% 
% 
% nexttile
% for nprobe = 1:2
%     [~,sorted_index] = sort(spindles_all(nprobe).peak_zscore(spindles_all(nprobe).SWS_index == 1));
%     x = 1:length(ipsi_amp_corr{nprobe});
%     y = mean(cumsum(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2));
%     LCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),2.5);
%     UCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),97.5);
%     PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
%     ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left spindles','Right spindles'},'Location','northwest','Box','off')
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('spindle events')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in spindle power correlation'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','image')

