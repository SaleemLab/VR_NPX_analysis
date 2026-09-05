
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

% % load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info_100ms.mat'),'UP_DOWN_info');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','ripple_info.mat'),'ripple_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');


%%%%% Ripple info (combined)

UP_ints=[];
DOWN_ints=[];
ripple_peaktimes=[];
ripple_ints=[];
SO_ints=[];
prev_down_idx=[];
next_down_idx=[];
V1_MUA_spiketimes = [];%{1} is L and {2} is R
HC_MUA_spiketimes = [];


for nprobe = 1:2
    V1_MUA_spiketimes{nprobe}=[];
    HC_MUA_spiketimes{nprobe}=[];

    UP_ints{nprobe}=slow_waves_all(nprobe).UP_ints;
    DOWN_ints{nprobe}=slow_waves_all(nprobe).DOWN_ints;
    SO_ints{nprobe} = slow_waves_all(nprobe).DOWN_intervals;
    ripple_peaktimes{nprobe}=ripples_all(nprobe).peaktimes(ripples_all(nprobe).SWS_index == 1);
    ripple_ints{nprobe}=[ripples_all(nprobe).onset(ripples_all(nprobe).SWS_index == 1) ripples_all(nprobe).offset(ripples_all(nprobe).SWS_index == 1)];
    % spindle_peaktimes{nprobe}=spindles_all(nprobe).peaktimes(spindles_all(nprobe).SWS_index == 1);
    % spindle_ints{nprobe}=[spindles_all(nprobe).onset(spindles_all(nprobe).SWS_index == 1) spindles_all(nprobe).offset(spindles_all(nprobe).SWS_index == 1)];

    for nsession = 1:max(slow_waves_all(1).UP_session_count)
        index = find(slow_waves_all(nprobe).DOWN_intervals_session == sessions_to_process(nsession));
        SO_ints{nprobe}(index,:) = SO_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession));
        UP_ints{nprobe}(index,:) = UP_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession));
        DOWN_ints{nprobe}(index,:) = DOWN_ints{nprobe}(index,:) + nsession * 1000000;

        index = find(ripples_all(nprobe).session_count(ripples_all(nprobe).SWS_index == 1) == sessions_to_process(nsession));
        ripple_ints{nprobe}(index,:) = ripple_ints{nprobe}(index,:) + nsession * 1000000;
        ripple_peaktimes{nprobe}(index,:) = ripple_peaktimes{nprobe}(index,:) + nsession * 1000000;

        % [C,ia,ib] = intersect(find(spindles_all(nprobe).session_count == sessions_to_process(nsession)),find(spindles_all(nprobe).SWS_index == 1));
        % spindle_ints{nprobe}(ib,:) = spindle_ints{nprobe}(ib,:) + nsession * 1000000;
        % if nprobe ==1 % only need once
        V1_MUA_spiketimes{nprobe} = [V1_MUA_spiketimes{nprobe}; slow_waves_all(nprobe).V1_MUA_spiketimes{nsession} + nsession * 1000000];
        HC_MUA_spiketimes{nprobe} = [HC_MUA_spiketimes{nprobe}; slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession} + nsession * 1000000];

        % slow_waves_all(nprobe).V1_MUA_spiketimes{nsession} = [slow_waves_all(nprobe).V1_MUA_spiketimes{nsession} + nsession * 1000000];
        % slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession} = [slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession} + nsession * 1000000];
        % slow_waves_all(2).V1_MUA_spiketimes{nsession} = [slow_waves_all(2).V1_MUA_spiketimes{nsession} + nsession * 1000000; HC_MUA_spiketimes{2}];
        % end
    end

    nUP = length(UP_ints{nprobe});
    prev_down_idx{nprobe} = nan(nUP, 1);
    next_down_idx{nprobe} = nan(nUP, 1);

    % tol = 1e-6; % Tolerance for floating point comparison

    % 1. Find Previous DOWN: Where DOWN_offset == UP_onset
    % We look for UP(:,1) inside the list of DOWN(:,2)
    [is_prev, prev_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,1), DOWN_ints{nprobe}(:,2), 1e-10);
    prev_down_idx{nprobe}(~is_prev) = nan; % Set zeros to NaN

    % 2. Find Next DOWN: Where DOWN_onset == UP_offset
    % We look for UP(:,2) inside the list of DOWN(:,1)
    [is_next, next_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,2), DOWN_ints{nprobe}(:,1), 1e-10);
    next_down_idx{nprobe}(~is_next) = nan; % Set zeros to NaN
end


merged_event_info.UP_ints = [UP_ints{1}(probability_psth_whole(1).UP_all_index,:); UP_ints{2}(probability_psth_whole(2).UP_all_index,:)];
merged_event_info.DOWN_ints = [DOWN_ints{1}(probability_psth_whole(1).DOWN_all_index,:); DOWN_ints{2}(probability_psth_whole(2).DOWN_all_index,:)];

merged_event_info.UP_hemisphere_id = [ones(length(probability_psth_whole(1).UP_all_index),1); 2*ones(length(probability_psth_whole(2).UP_all_index),1)];
merged_event_info.DOWN_hemisphere_id = [ones(length(probability_psth_whole(1).DOWN_all_index),1); 2*ones(length(probability_psth_whole(2).DOWN_all_index),1)];

merged_event_info.ripples_peaktimes = [ripple_peaktimes{1}; ripple_peaktimes{2}];
merged_event_info.ripples_ints = [ripple_ints{1}; ripple_ints{2}];
merged_event_info.ripples_hemisphere_id = [ones(length(ripple_peaktimes{1}),1); 2*ones(length(ripple_peaktimes{2}),1)];
[event_ids_first,event_ids_second] = merge_bilateral_ripple_events(merged_event_info.ripples_hemisphere_id,merged_event_info.ripples_peaktimes,0.05);

merged_event_info.ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id(event_ids_first);
% merged_event_info.ripples_original_index = merged_event_info.ripples_original_index(event_ids_first);
merged_event_info.ripples_peaktimes = merged_event_info.ripples_peaktimes(event_ids_first,:);
merged_event_info.ripples_ints = merged_event_info.ripples_ints(event_ids_first,:);


ripples_original_index = [find(ripples_all(1).SWS_index); find(ripples_all(2).SWS_index)];
merged_event_info.ripples_original_index = ripples_original_index(event_ids_first);
ripplePower = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index); ripples_all(2).peak_zscore(ripples_all(2).SWS_index)];
ripplePower = mean([ripplePower(event_ids_first) ripplePower(event_ids_second)],2);


V1_spiketimes = [];%{1} is L and {2} is R
HC_spiketimes = [];
slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession}
slow_waves_all(nprobe).V1_MUA_spiketimes{nsession}

L_HPC_MUA = [event_info(1).L_HPC_MUA_UP event_info(2).L_HPC_MUA_UP];
R_HPC_MUA = [event_info(1).R_HPC_MUA_UP event_info(2).R_HPC_MUA_UP];
L_V1_MUA = [event_info(1).L_V1_MUA_UP event_info(2).L_V1_MUA_UP];
R_V1_MUA = [event_info(1).R_V1_MUA_UP event_info(2).R_V1_MUA_UP];



for nevent = 1:length(merged_event_info.UP_ints)
% nevent = nevent + 2;
    tvec_this_event = (merged_event_info.UP_ints(nevent,1)+0.01/2:0.01:merged_event_info.UP_ints(nevent,2)-0.01/2);

    HPC_MUA_this_event = mean([L_HPC_MUA{nevent}; R_HPC_MUA{nevent}],1);
    if merged_event_info.UP_hemisphere_id==1
        ipsi_V1_MUA_this_event = L_V1_MUA{nevent};
        contra_V1_MUA_this_event = R_V1_MUA{nevent};
    else
        ipsi_V1_MUA_this_event = R_V1_MUA{nevent};
        contra_V1_MUA_this_event = L_V1_MUA{nevent};
    end


    Rippletidx=[];
    RippleIDs = [];
    for i = 1:length(merged_event_info.ripples_ints)
        [~,temp,~] = find(tvec_this_event>=merged_event_info.ripples_ints(i,1)-0.005&tvec_this_event<=merged_event_info.ripples_ints(i,2)+0.005);

        Rippletidx = [Rippletidx temp];
        if ~isempty(temp)
            RippleIDs = [RippleIDs i*ones(1,length(temp))];
        end
    end

    ripplePower(unique(RippleIDs))
    % figure;plot(mean([L_HPC_MUA{nevent}; R_HPC_MUA{nevent}],1));xline(tidx);
    % hold;plot(ipsi_V1_MUA_this_event)
    % ylim([0 1])
end





cumRipple
% Cumulative 
for nevent = 1:length(merged_event_info.UP_ints)
    merged_event_info.UP_ints(nevent,:)
    merged_event_info.ripples_ints(:,:)

    

end



    mean([event_info(nprobe).L_HPC_MUA_UP{nevent}; event_info(nprobe).R_HPC_MUA_UP{nevent}],'omitnan')





%% Effect of ripples and cumulative HPC activities  on UP survival probability
% save(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info_100ms.mat'),'UP_DOWN_info');
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info_100ms.mat'),'UP_DOWN_info');

UP_session_count = [slow_waves_all(1).UP_session_count(probability(1).UP_all_index); slow_waves_all(2).UP_session_count(probability(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);
UP_DOWN_info.subject_id = subject_id;
UP_DOWN_info.session_id = UP_session_count;


fig = figure('Name','Survival probability UP feature distribution')
fig.Position = [350 59 1000 750];

nexttile
histogram(UP_DOWN_info.last_ripples_power_UP,5:0.3:20,...
            'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.last_ripples_power_UP,[25 50 75]),'r')
xlabel('Last ripple power (z)')
ylabel('Proportion of UP events')
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

nexttile
histogram(UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP,0:1:110,...
        'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP,[25 50 75]),'r')
xlabel('Normalised cumulative ripple HC MUA')
ylabel('Proportion of UP events')
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

nexttile
histogram(UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP,0:1:110,...
        'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP,[25 50 75]),'r')
xlabel('Normalised cumulative ripple V1 MUA')
ylabel('Proportion of UP events')
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)


nexttile
histogram(UP_DOWN_info.first_half_HPC_MUA_cumulative_UP,0:1:80,...
    'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.first_half_HPC_MUA_cumulative_UP,[25 50 75]),'r')
ylabel('Proportion of UP events')
xlabel('Normalised cumulative HC MUA (1st half of UP)')
ylim([0 0.017])
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

nexttile
histogram(UP_DOWN_info.second_half_HPC_MUA_cumulative_UP,0:1:80,...
    'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.second_half_HPC_MUA_cumulative_UP,[25 50 75]),'r')
ylabel('Proportion of UP events')
xlabel('Normalised cumulative HC MUA (2nd half of UP)')
ylim([0 0.017])
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)


nexttile
histogram(UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP,-40:1:40,...
    'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP,[25 50 75]),'r')
ylabel('Proportion of UP events')
xlabel('Normalised cumulative HC MUA difference between first and second half of UP (2nd - 1st)')
ylim([0 0.015])
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)



nexttile
histogram(UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP,0:1:80,...
    'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP,[25 50 75]),'r')
ylabel('Proportion of UP events')
xlabel('Normalised cumulative V1 MUA (1st half of UP)')
ylim([0 0.01])
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

nexttile
histogram(UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP,0:1:80,...
    'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP,[25 50 75]),'r')
ylabel('Proportion of UP events')
xlabel('Normalised cumulative V1 MUA (2nd half of UP)')
ylim([0 0.01])
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

nexttile
histogram(UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP-UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP,-40:1:40,...
    'Normalization','probability','EdgeAlpha',0)
xline(prctile(UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP-UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP,[25 50 75]),'r')
ylabel('Proportion of UP events')
xlabel('Normalised cumulative V1 MUA difference between first and second half of UP (2nd - 1st)')
ylim([0 0.02])
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%% Ripple power predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_power_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition','subject_id',UP_DOWN_info.session_id...
);
save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_power_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripples_power_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',UP_DOWN_info.session_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_power_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


% 
% %%%%%%%%%%%%% Ripple duration predicts UP probability
% output = plot_UP_survival_probability( ...
%     {UP_DOWN_info.last_ripples_duration_UP}, ...
%     {UP_DOWN_info.time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
%     'event_option',[],'title_name', 'Last ripple duration and last ripple to UP-DOWN transition','subject_id',subject_id...
% );
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_duration_survival.mat'),'output');
% % save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option')
% 


% scatter([UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP], ...
%     [UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP],...
%     'filled','MarkerFaceAlpha','0.02')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% ylim([0 25])
% xlim([0 0.3])


% histogram(UP_DOWN_info.ipsi_last_ripples_power_UP)
% plot_UP_survival_probability( ...
%     {UP_DOWN_info.ipsi_last_ripples_power_UP, UP_DOWN_info.contra_last_ripples_power_UP}, ...
%     {UP_DOWN_info.ipsi_time_from_last_ripples_UP, UP_DOWN_info.contra_time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ipsi_ripple_counts_UP, UP_DOWN_info.contra_ripple_counts_UP},'count_option',1,'subject_id',subject_id, ...
%     'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition (1 ripple)'...
% );

%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ripple MUA peak predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripple_HPC_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id,'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak MUA and last ripple to UP-DOWN transition'...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.last_ripple_HPC_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Last ripple peak MUA and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',UP_DOWN_info.session_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])

%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last ripple V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_V1_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Last ripple ipsi V1 peak MUA and last ripple to UP-DOWN transition'...
);

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.contra_last_ripple_V1_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Last ripple contra V1 peak MUA and last ripple to UP-DOWN transition'...
);

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_V1_MUA_peak_UP - UP_DOWN_info.contra_last_ripple_V1_MUA_peak_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Last ripple ipsi contra V1 peak MUA diff and last ripple to UP-DOWN transition'...
);
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised cumulative ripple HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ripple activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Normalised cumulative ripple activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',UP_DOWN_info.session_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised cumulative ripple V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 ripple activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 ripple activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_V1_MUA_survival.mat'),'output');
% load(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_V1_MUA_survival.mat'),'output');

% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised cumulative ripple V1 MUA ipsi contra diff predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP - UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ipsi contra V1 ripple activity diff and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_ripple_V1_MUA_cumulative_UP - UP_DOWN_info.contra_normalised_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ipsi contra V1 ripple activity diff and last ripple to UP-DOWN transition (Shuffled)','timebin', 0.015...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_ripple_V1_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cumulative ripple HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative ripple activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ripple_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
    'event_option',[],'title_name', 'Cumulative ripple activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'subject_id',subject_id...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','cumulative_ripple_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cumulative ripple V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative ripple V1 activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_ripple_V1_MUA_cumulative_UP,UP_DOWN_info.contra_ripple_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Cumulative ripple V1 activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','cumulative_ripple_V1_MUA_survival.mat'),'output');
load(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','cumulative_ripple_V1_MUA_survival.mat'),'output');

% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])






%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%% Normalised cumulative V1 MUA predicts UP probability
% Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP V1 activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP V1 activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% load(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%%%%% Normalised cumulative HPC MUA predicts UP probability
% Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP HPC activity and last ripple to UP-DOWN transition','timebin', 0.015...
);
% save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','last_ripple_normalised_HPC_MUA_survival.mat'),'output');

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.normalised_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative UP HPC activity and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);

output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','normalised_cumulative_UP_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])



%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% 1st half Normalised cumulative HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','1st_half_normalised_cumulative_UP_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])

%%%%%%%%%%%%%%%% 1st half Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','1st_half_normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% 2nd half Normalised cumulative HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_half_normalised_cumulative_UP_HPC_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% 2nd half Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_half_normalised_cumulative_UP_V1_MUA_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])




%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% 2nd half - 1st half Normalised cumulative HPC MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd - 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative HPC MUA 2nd - 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


output.p_shuffled = output_shuffled.p;
output.b_shuffled = output_shuffled.b;

save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_1st_half_normalised_cumulative_UP_HPC_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%%%%%%%%%%%% 2nd half - 1st half Normalised cumulative V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd - 1st half and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP}, ....
    {UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',UP_DOWN_info.session_id, ...
    'event_option',[],'title_name', 'Normalised cumulative V1 MUA 2nd - 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_1st_half_normalised_cumulative_UP_V1_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])








%%%%%%%%%%%%%%%% Normalised cumulative ipsi-contra V1 MUA predicts UP probability
output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_normalised_V1_MUA_cumulative_UP,UP_DOWN_info.contra_normalised_V1_MUA_cumulative_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP},...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name','Normalised cumulative V1 MUA and last ripple to UP-DOWN transition','timebin', 0.015...
);

output_shuffled = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP, UP_DOWN_info.contra_second_half_V1_MUA_cumulative_UP - UP_DOWN_info.contra_first_half_V1_MUA_cumulative_UP}, ....
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'subject_id',subject_id, ...
    'event_option',[],'title_name', 'Normalised cumulative ipsi-contra V1 MUA 2nd - 1st half and last ripple to UP-DOWN transition (Shuffled)','shuffle_option',1,'timebin', 0.015...
);


for i = 1:length(output)
output(i).p_shuffled = output_shuffled(i).p;
output(i).b_shuffled = output_shuffled(i).b;
end


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','2nd_1st_half_normalised_cumulative_UP_V1_MUA_diff_survival.mat'),'output');
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[],'SVG_option',1)
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])


%%%%%% high Spindle power at the end of UP predicts faster UP termination
%%%%%% (but maybe a bit circular reasoning)

output = plot_UP_survival_probability( ...
    {double(UP_DOWN_info.ipsi_spindles_UP),double(UP_DOWN_info.contra_spindles_UP)}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',[],'timebin', 0.015, ...
    'title_name', 'UP with spindles and last ripple to UP-DOWN transition','event_option','binary','subject_id',subject_id...
    );

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_spindle_power_UP,UP_DOWN_info.contra_last_ripple_spindle_power_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',[],'timebin', 0.015, ...
    'title_name', 'last ripple UP spindle powers and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
    );

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_next_spindle_power_UP,UP_DOWN_info.contra_last_ripple_next_spindle_power_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',[],'timebin', 0.015, ...
    'title_name', 'last ripple next UP spindle power and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
    );

output = plot_UP_survival_probability( ...
    {UP_DOWN_info.ipsi_last_ripple_next_spindle_power_UP-UP_DOWN_info.ipsi_last_ripple_next_spindle_diff_UP,...
    UP_DOWN_info.contra_last_ripple_next_spindle_power_UP-UP_DOWN_info.contra_last_ripple_next_spindle_diff_UP}, ...
    {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
    {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',[],'timebin', 0.015, ...
    'title_name', 'last ripple end of UP spindle power and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
    );


% output = plot_UP_survival_probability( ...
%     {UP_DOWN_info.ipsi_last_ripple_next_spindle_diff_UP,UP_DOWN_info.contra_last_ripple_next_spindle_diff_UP}, ...
%     {UP_DOWN_info.time_from_last_ripples_UP,UP_DOWN_info.time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ripple_counts_UP,UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',[],'timebin', 0.015, ...
%     'title_name', 'last ripple next UP spindle power diff and last ripple to UP-DOWN transition','event_option',[],'subject_id',subject_id...
%     );

%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% 
% UP_DOWN_info.ipsi_last_ripple_spindle_power_UP

%% Multiple ripple features

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','UP_DOWN_info_100ms.mat'),'UP_DOWN_info');

%%%%%%%%%%%%% Ripple power predicts UP probability
% output = plot_UP_survival_probability( ...
%     {UP_DOWN_info.last_ripples_power_UP}, ...
%     {UP_DOWN_info.time_from_last_ripples_UP}, ...
%     {UP_DOWN_info.ripple_counts_UP},'count_option',[],'bilateral_merging',1, ...
%     'event_option',[],'title_name', 'Last ripple peak power and last ripple to UP-DOWN transition','subject_id',subject_id...
% );

%%%%% Ripple power + normalised cumulative ripple activity 
feature = [UP_DOWN_info.last_ripples_power_UP;...
    UP_DOWN_info.normalised_ripple_HPC_MUA_cumulative_UP;...
    UP_DOWN_info.second_half_HPC_MUA_cumulative_UP-UP_DOWN_info.first_half_HPC_MUA_cumulative_UP]';
    % UP_DOWN_info.ipsi_second_half_V1_MUA_cumulative_UP-UP_DOWN_info.ipsi_first_half_V1_MUA_cumulative_UP;]';
subject_id = UP_DOWN_info.session_id;
% subject_id = UP_DOWN_info.subject_id;
ttt = UP_DOWN_info.time_from_last_ripples_UP;
counts = UP_DOWN_info.ripple_counts_UP;
count_option = [];
shuffle_option = [];

% for hemi = 1:length(ripple_feature)  % 1 = ipsi, 2 = contra

% feature = ripple_feature{hemi};
% ttt = timeToTransition{hemi};
% counts = ripple_counts{hemi};
% Get valid entries
if isempty(count_option)
    valid_idx = intersect(find(~isnan(feature)), find(counts > 0));
elseif count_option == 1
    valid_idx = intersect(find(~isnan(feature)), find(counts == 1));
elseif count_option == 2
    valid_idx = intersect(find(~isnan(feature)), find(counts > 1));
end
feature_valid = feature(valid_idx,:);
ttt_valid = ttt(valid_idx);
subject_used =  subject_id(valid_idx);

%         if isempty(event_option)
%             % Compute thresholds
%             low_thresh = prctile(feature_valid, 25);
%             high_thresh = prctile(feature_valid, 75);
%         else
%             low_thresh = 0;
%             high_thresh = 1;
%         end

cox_time = ttt_valid;
cox_event = ones(size(cox_time));  % All UPs transition
feature_used = feature_valid;


% Bootstrap Cox Regression (1000x)
nBoot = 1000;
boot_b = nan(nBoot, size(feature,2));
boot_p = nan(nBoot, size(feature,2));
parfor iBoot = 1:nBoot
    s = RandStream('mrg32k3a','Seed',iBoot);

    boot_idx = datasample(s, 1:length(cox_time), length(cox_time));
    % boot_idx = 1:length(cox_time);
    feat_sample = feature_used(boot_idx,:);
    time_sample = cox_time(boot_idx);
    subj_sample = subject_used(boot_idx);

    try
        [b_tmp, ~, ~, stats_tmp] = coxphfit(feat_sample, time_sample, ...
            'Strata', subj_sample);
        boot_b(iBoot,:) = b_tmp;
        boot_p(iBoot,:) = stats_tmp.p;
    catch
        boot_b(iBoot,:) = NaN;
        boot_p(iBoot,:) = NaN;
    end
end

% Clean and save
boot_b = boot_b;
boot_p = boot_p;
boot_output.b = boot_b;
boot_output.p = boot_p;

% Store summary for bar plot
barData = median(boot_b);
lowerCI = median(boot_b) - prctile(boot_b, 2.5);
upperCI = prctile(boot_b, 97.5) - median(boot_b);
p50 = prctile(boot_p, 50);

% end



fig = figure
fig.Position = [350 59 515 600];


fig.Name = 'HPC combined COX regression';

x_pos = [1 2 3];
for  n = 1:3  % 1 = ipsi, 2 = contra
    hold on
    bar(x_pos(n), barData(n), 0.4, ...
        'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    errorbar(x_pos(n), barData(n), lowerCI(n), upperCI(n), ...
        'k', 'linestyle', 'none', 'linewidth', 1.5);
    text(x_pos(n), barData(n) + upperCI(n) + 0.0001, ...
        sprintf('p_{50%%} = %.3e', p50(n)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlim([0.5 3.5])
xticks(x_pos)
xticklabels({'Ripple power','Cumulative ripple MUA','2nd - 1st UP cumulative MUA'})
ylabel('Cox coefficient (b)')
title('Bootstrap Cox regression by hemisphere')
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)


save(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined','combined_HPC_effect_on_UP_survival.mat'),'boot_output');
save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','survival analysis combined'),[])





