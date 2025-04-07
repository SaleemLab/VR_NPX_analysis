function plot_bilateral_synchronisation(slow_waves_all,ripples_all,spindles_all,behavioural_state_merged_all,sessions_to_process)

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
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_SO_probability.mat'));
probability_ripples_SO = probability;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_ripples_probability.mat'));
probability_ripples_ripples = probability;


colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red for R
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue for L
colour_lines{3} = [255,185,205;254,145,198;228,42,168;182,0,140;122,1,119]/256;% 5 megenta for bilateral

colour_lines{4} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
colour_lines{5} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 



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


if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%% Ipsilateral vs contralateral correlation (UP DOWN)
duration_threshold = 2;
for nprobe = 1:2
    ipsi_amp_corr{nprobe} = [];
    contra_amp_corr{nprobe} = [];

    ipsi_plv{nprobe} = [];
    contra_plv{nprobe} = [];

    ipsi_amp_corr_bootstrap{nprobe} = [];
    contra_amp_corr_bootstrap{nprobe} = [];

    ipsi_plv_bootstrap{nprobe} = [];
    contra_plv_bootstrap{nprobe} = [];

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
        [C,ia,ib] = intersect(find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession)),probability_psth_whole(nprobe).DOWN_all_index);

        ipsi_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == nprobe);
        ipsi_shank(ipsi_shank==cortex_ref_shank(nsession,nprobe))=[];
        contra_shank = find(slow_waves_all(nprobe).probe_hemisphere{nsession} == mprobe);


        ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];

        ipsi_plv{nprobe} = [ipsi_plv{nprobe} mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        contra_plv{nprobe} = [contra_plv{nprobe} mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,ia)))];

        temp1 = [];
        temp2 = [];
        temp3 = [];
        temp4 = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,ia,size(ia,1));

            temp1(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            temp2(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).xcorr_r_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
            temp3(iBoot,:)= mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            temp4(iBoot,:)=  mean(squeeze(slow_waves_all(nprobe).plv_UD{nsession}(cortex_ref_shank(nsession,nprobe),contra_shank,event_id)));
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
    end

end


[p,h,stats] = signrank(ipsi_plv{1},contra_plv{1});
[p,h,stats] = signrank(ipsi_plv{2},contra_plv{2});
[p,h,stats] = signrank(ipsi_amp_corr{1},contra_amp_corr{1});
[p,h,stats] = signrank(ipsi_amp_corr{2},contra_amp_corr{2});

nfig = figure('Color','w','Name','Left Right V1 UP/DOWN ipsilateral-contralateral difference')
nfig.Position = [940 100 550 500];
orient(nfig,'landscape')

colour_lines = [44,123,182;215,25,28]/256;

nexttile
histcounts(ipsi_plv{1},0.5:0.001:1); hold on;
histogram(contra_plv{1},0.5:0.001:1)


nexttile
for nprobe = 1:2
    [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
    x = 1:length(ipsi_plv{nprobe});
    y = mean(cumsum(ipsi_plv_bootstrap{nprobe}(:,sorted_index) - contra_plv_bootstrap{nprobe}(:,sorted_index),2));
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),2.5);
    UCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
end
legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')

% xline(0,'r')
% title('Prob of left ripples during DOWN')
xlabel('UP-DOWN transition')
ylabel(sprintf('Cumulative ipsi-contralateral difference\n in phase-locking value'))
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
for nprobe = 1:2
      [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
    x = 1:length(ipsi_amp_corr{nprobe});
    y = mean(cumsum(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2));
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),2.5);
    UCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
end
legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')

% xline(0,'r')
% title('Prob of left ripples during DOWN')
xlabel('UP-DOWN transition')
ylabel(sprintf('Cumulative ipsi-contralateral difference\n in delta power correlation'))
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


% for nprobe = 1:2
%       [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
%     x = 1:length(ipsi_amp_corr{nprobe});
%     y = mean(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index));
%     %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
%     LCI = prctile(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2.5);
%     UCI = prctile(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index),97.5);
% 
%     PLOT = bar(x,y,'EdgeColor','none','FaceColor',colour_lines(nprobe,:),'FaceAlpha',0.5);hold on;
%     % ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
% end
% legend(ERROR_SHADE(1:2),{'Left V1 UP-DOWN','Right V1 UP-DOWN'},'Location','northwest','Box','off')
% 
% % xline(0,'r')
% % title('Prob of left ripples during DOWN')
% xlabel('UP-DOWN transition')
% ylabel(sprintf('Cumulative ipsi-contralateral difference\n in delta power correlation'))
% set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

%% Ipsilateral vs contralateral correlation (Ripple)
duration_threshold = 2;
for nprobe = 1:2
    ipsi_amp_corr{nprobe} = [];
    contra_amp_corr{nprobe} = [];

    ipsi_plv{nprobe} = [];
    contra_plv{nprobe} = [];

    ipsi_amp_corr_bootstrap{nprobe} = [];
    contra_amp_corr_bootstrap{nprobe} = [];

    ipsi_plv_bootstrap{nprobe} = [];
    contra_plv_bootstrap{nprobe} = [];

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

        if length(ipsi_shank)==1
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia))'];
        else
            ipsi_amp_corr{nprobe} = [ipsi_amp_corr{nprobe} mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
            ipsi_plv{nprobe} = [ipsi_plv{nprobe} mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,ia)))];
        end

        contra_amp_corr{nprobe} = [contra_amp_corr{nprobe} mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia)))];
        contra_plv{nprobe} = [contra_plv{nprobe} mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,ia)))];
        % if length(contra_amp_corr{nprobe}) ~= length( ipsi_amp_corr{nprobe})
        % 
        %     nsession
        % end

        temp1 = [];
        temp2 = [];
        temp3 = [];
        temp4 = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,ia,size(ia,1));
            if length(ipsi_shank)==1
                temp1(iBoot,:)= squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id));
                temp3(iBoot,:)= squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id));
            else
                temp1(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
                temp3(iBoot,:)= mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            end

            temp1(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            temp2(iBoot,:)= mean(squeeze(ripples_all(nprobe).xcorr_r{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id)));
            temp3(iBoot,:)= mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),ipsi_shank,event_id)));
            temp4(iBoot,:)=  mean(squeeze(ripples_all(nprobe).plv{nsession}(HPC_ref_shank(nsession,nprobe),contra_shank,event_id)));
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
    end
end


[p,h,stats] = signrank(ipsi_plv{1},contra_plv{1});
[p,h,stats] = signrank(ipsi_plv{2},contra_plv{2});
[p,h,stats] = signrank(ipsi_amp_corr{1},contra_amp_corr{1});
[p,h,stats] = signrank(ipsi_amp_corr{2},contra_amp_corr{2});

nfig = figure('Color','w','Name','Left Right ripples ipsilateral-contralateral difference')
nfig.Position = [940 100 550 500];
orient(nfig,'landscape')

colour_lines = [44,123,182;215,25,28]/256;

nexttile
for nprobe = 1:2
    [~,sorted_index] = sort(ripples_all(nprobe).offset(ripples_all(nprobe).SWS_index==1)-ripples_all(nprobe).onset(ripples_all(nprobe).SWS_index==1));
    % sorted_index = 1:length(ipsi_plv_bootstrap{nprobe});
    x = 1:length(ipsi_plv{nprobe});
    y = mean(cumsum(ipsi_plv_bootstrap{nprobe}(:,sorted_index) - contra_plv_bootstrap{nprobe}(:,sorted_index),2));
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),2.5);
    UCI = prctile(cumsum((ipsi_plv_bootstrap{nprobe}(:,sorted_index)-contra_plv_bootstrap{nprobe}(:,sorted_index)),2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
end
legend(ERROR_SHADE(1:2),{'Left ripples','Right ripples'},'Location','northwest','Box','off')

% xline(0,'r')
% title('Prob of left ripples during DOWN')
xlabel('UP-DOWN transition')
ylabel(sprintf('Cumulative ipsi-contralateral difference\n in phase-locking value'))
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

nexttile
for nprobe = 1:2
      % [~,sorted_index] = sort(slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,2)-slow_waves_all(nprobe).DOWN_ints(probability(nprobe).DOWN_all_index,1));
      % sorted_index = 1:length(ipsi_amp_corr_bootstrap{nprobe});
      [~,sorted_index] = sort(ripples_all(nprobe).offset(ripples_all(nprobe).SWS_index==1)-ripples_all(nprobe).onset(ripples_all(nprobe).SWS_index==1));
    x = 1:length(ipsi_amp_corr{nprobe});
    y = mean(cumsum(ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index) - contra_amp_corr_bootstrap{nprobe}(:,sorted_index),2));
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),2.5);
    UCI = prctile(cumsum((ipsi_amp_corr_bootstrap{nprobe}(:,sorted_index)-contra_amp_corr_bootstrap{nprobe}(:,sorted_index)),2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(nprobe) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');
end
legend(ERROR_SHADE(1:2),{'Left ripples','Right ripples'},'Location','northwest','Box','off')

% xline(0,'r')
% title('Prob of left ripples during DOWN')
xlabel('UP-DOWN transition')
ylabel(sprintf('Cumulative ipsi-contralateral difference\n in ripple power correlation'))
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

if exist(fullfile(analysis_folder,'V1-HPC bilateral interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC bilateral interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[])


%%

nexttile
colour_lines = [74,20,134;0,90,50]/256;
s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
sampled_events = datasample(s,1:length(ipsi_plv{1}),1000);

grp = [ones(length(sampled_events),1);ones(length(sampled_events),1)*2];
tst=[ipsi_plv{1}(sampled_events) contra_plv{1}(sampled_events)]';

% xbe = beeswarm(grp,tst,'sort_style','nosort','colormap',[PP.RUN1T1;PP.RUN1T2;PP.RUN2T1;PP.RUN2T2],'dot_size',2,'corral_style','rand');
xbe = beeswarm(grp,tst,'sort_style','nosort','colormap',colour_lines,'dot_size',2,'overlay_style','sd','corral_style','rand');

% yticks([0:0.02:0.06])
xticks([1:2])
xticklabels({'POST1 T1','POST1 T2','POST2 T1','POST2 T2'})
ylabel('Replay rate (events/sec)')
set(gca,'FontSize',12)
ylim([0 0.07])
hold on

tst=[INTER_T1_rate_events; INTER_T2_rate_events; FINAL_RT1_rate_events; FINAL_RT2_rate_events]';

xbe = reshape(xbe,size(tst,1),size(tst,2));
for i = 1:size(tst,1)
    plot(xbe(i,[1 2]),tst(i,[1 2]),'Color',[0,0,0,0.2])
    plot(xbe(i,[3 4]),tst(i,[3 4]),'Color',[0,0,0,0.2])
end

axis square
title(sprintf('POST replay rate for both exposures (%s)',rest_option));



%% Plot UP DOWN ipsilateral vs contralateral







%% Plotting distribution of ripple during normalised duration of UP  (peaktimes)
probability = probability_normalised;
probability_baseline = probability_normalised_baseline;

time_wondows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;
probe_hemisphere_texts = {'Probability of ripples during left V1 normalised UP-DOWN duration','Probability of ripples during right V1 normalised UP-DOWN duration'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};

colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])




%% plotting Probability of SWR relative to UP and DOWN onset (peaktimes)
probability = probability_psth;
probability_baseline = probability_psth_baseline;

time_wondows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of ripples during left V1 UP-DOWN','Probability of ripples during right V1 UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
colour_lines = [215,25,28;44,123,182]/256;

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')
    xlim([-0.5 0.5])

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])


    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])
end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% Plotting distribution of ripple during normalised duration of UP  (WHOLE)
probability = probability_normalised_whole;
probability_baseline = probability_normalised_whole_baseline;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;
% probe_hemisphere_texts = {'Probability of ripples during left V1 normalised UP-DOWN duration','Probability of ripples during right V1 normalised UP-DOWN duration'};
probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};
colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = linspace(0,1,num_bins);
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = linspace(0,1,num_bins);
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])




%% plotting Probability of SWR relative to UP and DOWN onset (WHOLE)
probability = probability_psth_whole;
probability_baseline = probability_psth_whole_baseline;

time_wondows = [-1 1];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
% probe_hemisphere_texts = {'Probability of ripples during left V1 UP-DOWN','Probability of ripples during right V1 UP-DOWN'};
probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};


for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')
    xlim([-0.5 0.5])

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).L_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])


    %%%% Right ripples
    % all_ripple_no = probability(nprobe).R_ripple_no;

    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Time relative to DOWN onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = cumsum(sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability_baseline(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to UP onset')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = sum(probability_baseline(nprobe).R_ripples_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability_baseline(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Time relative to UP onset')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.5 0.5])
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])
%% plotting Probability of UP and/or DOWN relative to Ripple peaktime
probability = probability_ripples_SO;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left UP-DOWN during ripples','Probability of right UP-DOWN during ripples'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_UP_no = length(probability(nprobe).UP_all_index);
    % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    all_ripple_no = probability(nprobe).L_ripple_no;

    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_DOWN)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).L_ripples_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_DOWN)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).L_ripples_DOWN_shuffled);
    LCI = prctile(probability(nprobe).L_ripples_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).L_ripples_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during left ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples_UP)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).L_ripples_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples_UP)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).L_ripples_UP_shuffled);
    LCI = prctile(probability(nprobe).L_ripples_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).L_ripples_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during left ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    %%%% DOWN
    all_ripple_no = probability(nprobe).R_ripple_no;

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_DOWN)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).R_ripples_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_DOWN)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).R_ripples_DOWN_shuffled);
    LCI = prctile(probability(nprobe).R_ripples_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).R_ripples_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during right ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples_UP)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).R_ripples_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples_UP)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).R_ripples_UP_shuffled);
    LCI = prctile(probability(nprobe).R_ripples_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).R_ripples_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during right ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% plotting Probability of UP and/or DOWN relative to Ipsilateral UP DOWN
probability = probability_SO_SO;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left UP-DOWN relative to ipsilateral UP-DOWN','Probability of right UP-DOWN relative to ipsilateral UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);

    %%%% DOWN during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_UP_shuffled);
    LCI = prctile(probability(nprobe).DOWN_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% DOWN during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_DOWN_shuffled);
    LCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% UP during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_DOWN_shuffled);
    LCI = prctile(probability(nprobe).UP_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_UP_shuffled);
    LCI = prctile(probability(nprobe).UP_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during UP')
    xlabel('Time relative to UP (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

%% plotting Probability of UP and/or DOWN relative to contralateral UP DOWN
probability = probability_SO_SO_contralateral;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left UP-DOWN relative to contralateral UP-DOWN','Probability of right UP-DOWN relative to contralateral UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    if nprobe ==1
        mprobe = 2;
    else
        mprobe = 1;
    end

    all_UP_no = length(probability(mprobe).UP_all_index);
    all_DOWN_no = length(probability(mprobe).DOWN_all_index);

    %%%% DOWN during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_UP_shuffled);
    LCI = prctile(probability(nprobe).DOWN_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% DOWN during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).DOWN_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).DOWN_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).DOWN_DOWN_shuffled);
    LCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).DOWN_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of DOWN during DOWN')
    xlabel('Time relative to DOWN (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% UP during U-D
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_DOWN)./all_DOWN_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_DOWN_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_DOWN_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_DOWN)./all_DOWN_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_DOWN_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_DOWN_shuffled);
    LCI = prctile(probability(nprobe).UP_DOWN_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_DOWN_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during DOWN')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    %%% UP during D-U
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).UP_UP)./all_UP_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).UP_UP_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).UP_UP_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to UP (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).UP_UP)./all_UP_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).UP_UP_bootstrap,2.5);
    UCI = prctile(probability(nprobe).UP_UP_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).UP_UP_shuffled);
    LCI = prctile(probability(nprobe).UP_UP_shuffled,2.5);
    UCI = prctile(probability(nprobe).UP_UP_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of UP during UP')
    xlabel('Time relative to UP (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])



%% plotting Probability of ripple relative to Ripple peaktime
probability = probability_ripples_ripples;

time_wondows = [-0.5 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
% all_sessions = max(slow_waves_all(1).DOWN_session_count);
colour_lines = [215,25,28;253,174,97;171,217,233;44,123,182]/256;


colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'Probability of left ripples during ripples','Probability of right ripples during ripples'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};

for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    all_ripple_no = length(probability(nprobe).L_ripple_no);

    %%%% Left ripples plotting
    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).L_ripples)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).L_ripples_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).L_ripples_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).L_ripples_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of left ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).L_ripples)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).L_ripples_bootstrap,2.5);
    UCI = prctile(probability(nprobe).L_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).L_ripples_shuffled);
    LCI = prctile(probability(nprobe).L_ripples_shuffled,2.5);
    UCI = prctile(probability(nprobe).L_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of ripples during left ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples
    %%%% DOWN
    all_ripple_no = probability(nprobe).R_ripple_no;

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = cumsum(sum(probability(nprobe).R_ripples)./all_ripple_no);
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_bootstrap,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_bootstrap,2),97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(cumsum(probability(nprobe).R_ripples_shuffled,2));
    LCI = prctile(cumsum(probability(nprobe).R_ripples_shuffled,2),2.5);
    UCI = prctile(cumsum(probability(nprobe).R_ripples_shuffled,2),97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Location','northwest','Box','off')

    % xline(0,'r')
    % title('Prob of right ripples during DOWN')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile
    x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
    y = sum(probability(nprobe).R_ripples)./all_ripple_no;
    %     y = mean(cumsum(probability(nprobe).R_ripples_DOWN_bootstrap,2));
    LCI = prctile(probability(nprobe).R_ripples_bootstrap,2.5);
    UCI = prctile(probability(nprobe).R_ripples_bootstrap,97.5);

    PLOT = plot(x,y,'Color',colour_lines(nprobe,:));hold on;
    ERROR_SHADE(1) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines(nprobe,:),'FaceAlpha','0.3','LineStyle','none');

    y = mean(probability(nprobe).R_ripples_shuffled);
    LCI = prctile(probability(nprobe).R_ripples_shuffled,2.5);
    UCI = prctile(probability(nprobe).R_ripples_shuffled,97.5);

    PLOT = plot(x,y,'k');hold on;
    ERROR_SHADE(2) = patch([x fliplr(x)],[UCI fliplr(LCI)],'k','FaceAlpha','0.3','LineStyle','none');
    % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
    % xline(0,'r')
    title('Probability of ripples during right ripples')
    xlabel('Time relative to ripple peaktime (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])
