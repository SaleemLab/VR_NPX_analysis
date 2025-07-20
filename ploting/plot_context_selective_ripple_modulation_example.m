
%% Plotting context selecitve ripple modulation
% scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')
load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');


ripple_modulation_percentile = [context_modulation_all.ripple_modulation_percentile{:}];
PRE_ripple_modulation_percentile = [context_modulation_all.modulation_percentile_PRE{:}];
ripple_modulation_percentile = max([ripple_modulation_percentile; PRE_ripple_modulation_percentile])>0.95;
ripple_modulation_id = ripple_modulation_percentile>0.95;

session_id_all = [context_modulation_all.session_id{:}];
regions_all = [context_modulation_all.region{:}];
regions_all(ripple_modulation_id==0)=nan;

FR_track = [context_modulation_all.FR_track{:}];
z_FR_track = [context_modulation_all.z_FR_track{:}];
z_FR_track_diff = z_FR_track(1,:)-z_FR_track(2,:);

FR_track = [context_modulation_all.FR_track{:}];
FR_track_diff = FR_track(1,:)-FR_track(2,:);




colorlines = [ ...
    0.2, 0.4, 0.8    % blue
    0.85, 0.2, 0.2;  % red
    ];


%%%%%%%%%%%%%%%%%% context selective ripple modulation based on different
%%%%%%%%%%%%%%%%%% conditions (ripple power or spindle power or SO phase and etc)

ripple_types = {'ALL','LOW','HIGH'};
PSTH_fields = {'PSTH_diff','PSTH_diff_low_ripple','PSTH_diff_high_ripple'};
% conditions = {'Ripple power',''}
% Preallocate
all_vars_V1 = struct();
all_vars_HPC = struct();

for iType = 1:numel(ripple_types)
    tag = ripple_types{iType};
    PSTH_data = vertcat(context_modulation_all.(PSTH_fields{iType}){:});

    % Compute FR diffs
    POST = mean(PSTH_data(:, context_modulation_all.timebin > 0 & context_modulation_all.timebin < 0.2), 2, 'omitnan');
    PRE = mean(PSTH_data(:, context_modulation_all.timebin > -0.2 & context_modulation_all.timebin < 0), 2, 'omitnan');
    SHIFT = mean(PSTH_data(:, context_modulation_all.timebin > -1 & context_modulation_all.timebin < -0.8), 2, 'omitnan');

    % Index
    isV1 = contains(regions_all,'V1');
    isHPC = contains(regions_all,'HPC');

    % Store for V1
    all_vars_V1.(['POST_' tag]) = normalize(double(POST(isV1)));
    all_vars_V1.(['PRE_' tag]) = normalize(double(PRE(isV1)));
    all_vars_V1.(['SHIFT_' tag]) = normalize(double(SHIFT(isV1)));

    % Store for HPC
    all_vars_HPC.(['POST_' tag]) = normalize(double(POST(isHPC)));
    all_vars_HPC.(['PRE_' tag]) = normalize(double(PRE(isHPC)));
    all_vars_HPC.(['SHIFT_' tag]) = normalize(double(SHIFT(isHPC)));
end

% Add z_FR_track_diff and subjectID
z_FR_V1 = normalize(double(z_FR_track_diff(contains(regions_all,'V1'))'));
z_FR_HPC = normalize(double(z_FR_track_diff(contains(regions_all,'HPC'))'));
subjectID_V1 = categorical(session_id_all(contains(regions_all,'V1'))');
subjectID_HPC = categorical(session_id_all(contains(regions_all,'HPC'))');

% Construct tables
tbl_V1 = table(z_FR_V1, ...
    all_vars_V1.POST_ALL, all_vars_V1.POST_LOW, all_vars_V1.POST_HIGH, ...
    all_vars_V1.PRE_ALL, all_vars_V1.PRE_LOW, all_vars_V1.PRE_HIGH, ...
    all_vars_V1.SHIFT_ALL, all_vars_V1.SHIFT_LOW, all_vars_V1.SHIFT_HIGH, ...
    subjectID_V1, ...
    'VariableNames', {'z_FR_track_diff', ...
    'POST_ripple_FR_diff_ALL','POST_ripple_FR_diff_LOW','POST_ripple_FR_diff_HIGH', ...
    'PRE_ripple_FR_diff_ALL','PRE_ripple_FR_diff_LOW','PRE_ripple_FR_diff_HIGH', ...
    'SHIFT_ripple_FR_diff_ALL','SHIFT_ripple_FR_diff_LOW','SHIFT_ripple_FR_diff_HIGH', ...
    'subjectID'});

tbl_HPC = table(z_FR_HPC, ...
    all_vars_HPC.POST_ALL, all_vars_HPC.POST_LOW, all_vars_HPC.POST_HIGH, ...
    all_vars_HPC.PRE_ALL, all_vars_HPC.PRE_LOW, all_vars_HPC.PRE_HIGH, ...
    all_vars_HPC.SHIFT_ALL, all_vars_HPC.SHIFT_LOW, all_vars_HPC.SHIFT_HIGH, ...
    subjectID_HPC, ...
    'VariableNames', {'z_FR_track_diff', ...
    'POST_ripple_FR_diff_ALL','POST_ripple_FR_diff_LOW','POST_ripple_FR_diff_HIGH', ...
    'PRE_ripple_FR_diff_ALL','PRE_ripple_FR_diff_LOW','PRE_ripple_FR_diff_HIGH', ...
    'SHIFT_ripple_FR_diff_ALL','SHIFT_ripple_FR_diff_LOW','SHIFT_ripple_FR_diff_HIGH', ...
    'subjectID'});

ModelList = {
    'POST_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)';
    'POST_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)';
    'POST_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)';

    'PRE_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)';
    'PRE_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)';
    'PRE_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)';

    'SHIFT_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)';
    'SHIFT_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)';
    'SHIFT_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)';

    'z_FR_track_diff ~ POST_ripple_FR_diff_ALL + PRE_ripple_FR_diff_ALL + SHIFT_ripple_FR_diff_ALL + (1|subjectID)';
    'z_FR_track_diff ~ POST_ripple_FR_diff_LOW + PRE_ripple_FR_diff_LOW + SHIFT_ripple_FR_diff_LOW + (1|subjectID)';
    'z_FR_track_diff ~ POST_ripple_FR_diff_HIGH + PRE_ripple_FR_diff_HIGH + SHIFT_ripple_FR_diff_HIGH + (1|subjectID)';

    'z_FR_track_diff ~ POST_ripple_FR_diff_LOW + POST_ripple_FR_diff_HIGH + (1|subjectID)';
    'z_FR_track_diff ~ PRE_ripple_FR_diff_LOW + PRE_ripple_FR_diff_HIGH + (1|subjectID)';
    'z_FR_track_diff ~ SHIFT_ripple_FR_diff_LOW + SHIFT_ripple_FR_diff_HIGH + (1|subjectID)';
    };


%%%%%%%%%
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
load(fullfile(analysis_folder,'session_clusters_all_POST.mat'))
z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);
clear z_bias1

%%%%%%%%%%%%% ripple PSTH (cell examples)
timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;
bins_to_use = bin_centers>0 & bin_centers<0.1;

psthBinSize = 0.01;
windows = [-1.5 1.5];
sessions_to_process = 1:max(ripples_all(1).session_count);

for nsession = 1:22

    session_id = [ripples_all(1).session_count(ripples_all(1).SWS_index); ripples_all(2).session_count(ripples_all(2).SWS_index)];

    all_clusters = session_clusters_all.spatial_cell_id{nsession};

    event_times = [ripples_all(1).onset(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).onset(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)];
    event_id = [ones(sum(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1),1); ones(sum((ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)),1)];

    ripple_modulation = ripple_modulation_analysis(session_clusters_all.spike_times{nsession},session_clusters_all.spike_id{nsession},windows,psthBinSize,...
        'unit_id',all_clusters,'event_times',event_times,'event_id',event_id,'saving_PSTH',1,'shuffle_option',0);


    amplitudes = [ripples_all(1).peak_zscore(ripples_all(1).session_count == nsession&ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).session_count == nsession&ripples_all(2).SWS_index==1)]';
    mean_bias = mean(z_bias(bins_to_use,session_id == sessions_to_process(nsession)),1,'omitnan');

    power_threshold = prctile(amplitudes,[25 75]);
    log_odds_threshold = prctile(mean_bias,[20 80]);

    T1_index = find(mean_bias > log_odds_threshold(2));
    T2_index = find(mean_bias < log_odds_threshold(1));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_Hz{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH(nCell,T1_index,:)));
        context_modulation_all.PSTH_Hz{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH(nCell,T2_index,:)));
    end


    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes >= power_threshold(end));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes >= power_threshold(end));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_high_ripple_Hz{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH(nCell,T1_index,:)));
        context_modulation_all.PSTH_high_ripple_Hz{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH(nCell,T2_index,:)));
    end


    T1_index = find(mean_bias > log_odds_threshold(2) & amplitudes <= power_threshold(1));
    T2_index = find(mean_bias < log_odds_threshold(1) & amplitudes <= power_threshold(1));

    for nCell = 1:length(all_clusters)
        % Ripple PSTH
        context_modulation_all.PSTH_low_ripple_Hz{nsession}(1,nCell,:) = mean(squeeze(ripple_modulation.PSTH(nCell,T1_index,:)));
        context_modulation_all.PSTH_low_ripple_Hz{nsession}(2,nCell,:) = mean(squeeze(ripple_modulation.PSTH(nCell,T2_index,:)));
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
    all_PSTH = [context_modulation_all.PSTH_high_ripple{:}];
    PSTH_Hz = [context_modulation_all.PSTH_high_ripple_Hz{nsession}];

    % all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
    PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
    shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');

    colorlines = [ ...
        0.2, 0.4, 0.8    % blue
        0.85, 0.2, 0.2;  % red
        ];


    X = double(z_FR_track_diff)';
    Y = double(POST_ripple_FR_diff);
    Z = session_id_all';


    cluster_id = zeros(length(Z),1);
    cluster_id(Z == nsession) = 1:sum(Z == nsession);

    T2_cell = cluster_id(X < -0.3 & Y < -0.1 & Z == nsession & contains(regions_all','V1_L'));
    T2_cell_all = find((X < -0.3 & Y < -0.1 & Z == nsession & contains(regions_all','V1_L')));

    T1_cell = cluster_id(X > 0.3 & Y > 0.1 & Z == nsession & contains(regions_all','V1_R'));
    T1_cell_all = find((X > 0.3 & Y > 0.1 & Z == nsession & contains(regions_all','V1_R')));


    
    PSTH_Hz = [context_modulation_all.PSTH_Hz{nsession}];

    high_PSTH_Hz = [context_modulation_all.PSTH_high_ripple_Hz{nsession}];
    low_PSTH_Hz = [context_modulation_all.PSTH_low_ripple_Hz{nsession}];

    if ~isempty(T2_cell)
        for ncell = 1:length(T2_cell)
            fig = figure;
            fig.Position = [1385 650 470 340];
            fig.Name = sprintf('Session %i V1 L spatial and ripple PSTH %i (%i)', nsession,session_clusters_all.spatial_cell_id{nsession}(T2_cell(ncell)),T2_cell_all(ncell));
            tlo = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % fixed 4-subplot layout

            %%% all ripples
            nexttile
            plot(1:140,mean(session_clusters_all.spatial_response{nsession}{T2_cell(ncell),1}),'Color',colorlines(1,:));hold on;
            plot(1:140,mean(session_clusters_all.spatial_response{nsession}{T2_cell(ncell),2}),'Color',colorlines(2,:))
            % title(sprintf('V1 L %i',session_clusters_all.spatial_cell_id{nsession}(T2_cell(ncell))))
            set(gca,'TickDir','out','Box','off','FontSize',12)
            xlabel('Position (cm)')
            ylabel('Firing rate (Hz)')
            % ylim([0 20])
            xlim([0 140])
            legend('Track L','Track R','box','off')

            nexttile
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(1,T2_cell_all(ncell),:)),'Color',colorlines(1,:));hold on;
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(2,T2_cell_all(ncell),:)),'Color',colorlines(2,:))
            %
            plot(context_modulation_all.timebin,squeeze(PSTH_Hz(1,T2_cell(ncell),:)),'Color',colorlines(1,:));hold on;
            plot(context_modulation_all.timebin,squeeze(PSTH_Hz(2,T2_cell(ncell),:)),'Color',colorlines(2,:))

            xlim([-1 1])
            xlabel('Time (s)')
            ylabel('Firing rate (Hz)')
            xlim([-1 1])
            title('all ripples')
            set(gca,'TickDir','out','Box','off','FontSize',12)


            nexttile
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(1,T2_cell_all(ncell),:)),'Color',colorlines(1,:));hold on;
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(2,T2_cell_all(ncell),:)),'Color',colorlines(2,:))
            %
            plot(context_modulation_all.timebin,squeeze(low_PSTH_Hz(1,T2_cell(ncell),:)),'Color',colorlines(1,:));hold on;
            plot(context_modulation_all.timebin,squeeze(low_PSTH_Hz(2,T2_cell(ncell),:)),'Color',colorlines(2,:))

            xlim([-1 1])
            xlabel('Time (s)')
            ylabel('Firing rate (Hz)')
            xlim([-1 1])
            title('low ripples')
            set(gca,'TickDir','out','Box','off','FontSize',12)

            nexttile
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(1,T2_cell_all(ncell),:)),'Color',colorlines(1,:));hold on;
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(2,T2_cell_all(ncell),:)),'Color',colorlines(2,:))
            %
            plot(context_modulation_all.timebin,squeeze(high_PSTH_Hz(1,T2_cell(ncell),:)),'Color',colorlines(1,:));hold on;
            plot(context_modulation_all.timebin,squeeze(high_PSTH_Hz(2,T2_cell(ncell),:)),'Color',colorlines(2,:))

            xlim([-1 1])
            xlabel('Time (s)')
            ylabel('Firing rate (Hz)')
            xlim([-1 1])
            title('high ripples')
            set(gca,'TickDir','out','Box','off','FontSize',12)

        end
    end

    PSTH_Hz = [context_modulation_all.PSTH_Hz{nsession}];

    high_PSTH_Hz = [context_modulation_all.PSTH_high_ripple_Hz{nsession}];
    low_PSTH_Hz = [context_modulation_all.PSTH_low_ripple_Hz{nsession}];

    if ~isempty(T1_cell)
        for ncell = 1:length(T1_cell)
            fig = figure;
            fig.Position = [1385 650 470 340];
            fig.Name = sprintf('Session %i V1 R spatial and ripple PSTH %i (%i)', nsession,session_clusters_all.spatial_cell_id{nsession}(T1_cell(ncell)),T1_cell_all(ncell));
            tlo = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % fixed 4-subplot layout

            %%% all ripples
            nexttile
            plot(1:140,mean(session_clusters_all.spatial_response{nsession}{T1_cell(ncell),1}),'Color',colorlines(1,:));hold on;
            plot(1:140,mean(session_clusters_all.spatial_response{nsession}{T1_cell(ncell),2}),'Color',colorlines(2,:))
            % title(sprintf('V1 L %i',session_clusters_all.spatial_cell_id{nsession}(T2_cell(ncell))))
            set(gca,'TickDir','out','Box','off','FontSize',12)
            xlabel('Position (cm)')
            ylabel('Firing rate (Hz)')
            % ylim([0 20])
            xlim([0 140])
            legend('Track L','Track R','box','off')

            nexttile
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(1,T2_cell_all(ncell),:)),'Color',colorlines(1,:));hold on;
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(2,T2_cell_all(ncell),:)),'Color',colorlines(2,:))
            %
            plot(context_modulation_all.timebin,squeeze(PSTH_Hz(1,T1_cell(ncell),:)),'Color',colorlines(1,:));hold on;
            plot(context_modulation_all.timebin,squeeze(PSTH_Hz(2,T1_cell(ncell),:)),'Color',colorlines(2,:))

            xlim([-1 1])
            xlabel('Time (s)')
            ylabel('Firing rate (Hz)')
            xlim([-1 1])
            title('all ripples')
            set(gca,'TickDir','out','Box','off','FontSize',12)


            nexttile
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(1,T2_cell_all(ncell),:)),'Color',colorlines(1,:));hold on;
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(2,T2_cell_all(ncell),:)),'Color',colorlines(2,:))
            %
            plot(context_modulation_all.timebin,squeeze(low_PSTH_Hz(1,T1_cell(ncell),:)),'Color',colorlines(1,:));hold on;
            plot(context_modulation_all.timebin,squeeze(low_PSTH_Hz(2,T1_cell(ncell),:)),'Color',colorlines(2,:))

            xlim([-1 1])
            xlabel('Time (s)')
            ylabel('Firing rate (Hz)')
            xlim([-1 1])
            title('low ripples')
            set(gca,'TickDir','out','Box','off','FontSize',12)

            nexttile
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(1,T2_cell_all(ncell),:)),'Color',colorlines(1,:));hold on;
            % plot(context_modulation_all.timebin,squeeze(all_PSTH(2,T2_cell_all(ncell),:)),'Color',colorlines(2,:))
            %
            plot(context_modulation_all.timebin,squeeze(high_PSTH_Hz(1,T1_cell(ncell),:)),'Color',colorlines(1,:));hold on;
            plot(context_modulation_all.timebin,squeeze(high_PSTH_Hz(2,T1_cell(ncell),:)),'Color',colorlines(2,:))

            xlim([-1 1])
            xlabel('Time (s)')
            ylabel('Firing rate (Hz)')
            xlim([-1 1])
            title('high ripples')
            set(gca,'TickDir','out','Box','off','FontSize',12)

        end
    end


    save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','context-selective ripple PSTH'),[],'ContentType','vector')

end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')


%%%%%%%%%%%%% High ripple

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples Track prefering with neuron labelled)';
nfig.Position = [   842   345   954   578];


subplot(2,3,1)
% X = double(z_FR_track_diff(contains(regions_all,'V1')))';

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
is_V1R = contains(regions_all, 'V1')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'V1')&z_FR_track_diff<0;


scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
xline(0,'k--')
yline(0,'k--')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.45 0.45])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('0 to 0.2s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,2)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on

X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'V1')));

scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-0.2 to 0s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,3)

% is_V1R = contains(regions_all, 'V1_R');
% is_V1L = contains(regions_all, 'V1_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

hold on
X = double(z_FR_track_diff(contains(regions_all,'V1')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'V1')));
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
hold on
xline(0,'k--')
yline(0,'k--')

valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
title('-1 to -0.8s relative to ripple')
set(gca,'TickDir','out','Box','off','FontSize',12)


subplot(2,3,4)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

is_V1R = contains(regions_all, 'HPC')&z_FR_track_diff>0;
is_V1L = contains(regions_all, 'HPC')&z_FR_track_diff<0;

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(POST_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-1 1])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



subplot(2,3,5)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');

scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];


X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(PRE_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.7 0.7])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)

% is_V1R = contains(regions_all, 'HPC_R');
% is_V1L = contains(regions_all, 'HPC_L');
scolors = repmat([nan nan nan], length(regions_all), 1);  % red
scolors(is_V1R, :) = repmat(colorlines(1,:), sum(is_V1R), 1);  % red
scolors(is_V1L, :) = repmat(colorlines(2,:), sum(is_V1L), 1);  % blue
scolors(isnan(scolors(:,1)),:) = [];

X = double(z_FR_track_diff(contains(regions_all,'HPC')))';
Y = double(shifted_ripple_FR_diff(contains(regions_all,'HPC')));
% scatter(X,Y,'k','filled','MarkerFaceAlpha', 0.1)
scatter(X,Y,20,scolors,'filled','MarkerFaceAlpha', 0.1)

hold on
xline(0,'k--')
yline(0,'k--')


valid_id = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(valid_id);
Y = Y(valid_id);
coeffs = polyfit(X, Y, 1);
x_fit = linspace(min(X), max(X), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_HIGH ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



%%%% Loop though bins


%%%%%%%%%%% kstest2
ks2stat_all = nan(1000,10);

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 (high ripples bins) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');

% timewindows = -1:0.2:1;
timewindows = -1:0.2:1;
for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
    text(0.4,0.05,sprintf('p = %.3e',pval))
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ylim([0 1])
    xline(0,'k--','LineWidth',2)
    xlabel('Track L - Track R Ripple FR diff (z)')
    ylabel('cum prop of cells')
    legend('V1 Track R prefering','V1 Track L prefering','box','off')
    title(sprintf('%.1f to %.1fs relative to ripples',timewindows(nbin),timewindows(nbin+1)))
    set(gca,'TickDir','out','Box','off','FontSize',12)
end



%%%%%%%%%%% kstest2
ks2stat_all = nan(1000,10);

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 (low ripples bins) (200ms bins)';
nfig.Position = [640         253        1239         725];

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
% shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.4&context_modulation_all.timebin<-0.2),2,'omitnan');

timewindows = -1:0.2:1;

for nbin = 1:length(timewindows)-1
    nexttile

    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(nbin)&context_modulation_all.timebin<timewindows(nbin+1)),2,'omitnan');
    [h,pval,ks2stat] = kstest2(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),'Tail','smaller');


    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(2,:));
    hold on;histogram(POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0),-0.5:0.02:0.5,'Normalization','probability','FaceColor',colorlines(1,:))
    text(0.4,0.05,sprintf('p = %.3e',pval))
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ylim([0 1])
    xline(0,'k--','LineWidth',2)
    xlabel('Track L - Track R Ripple FR diff (z)')
    ylabel('cum prop of cells')
    legend('V1 Track R prefering','V1 Track L prefering','box','off')
    title(sprintf('%.1f to %.1fs relative to ripples',timewindows(nbin),timewindows(nbin+1)))
    set(gca,'TickDir','out','Box','off','FontSize',12)
   

end



%%%%%%%%%%% KS max difference bootstrap
bin_width = 0.1;    % 50 ms
step_size = 0.01;   % 10 ms
t_start = -0.5;
t_end = 0.5;

% Bin centers
bin_centers = t_start:step_size:t_end;

% Bin edges (each bin spans 100 ms centered at bin_centers)
timewindows = [bin_centers - bin_width/2; bin_centers + bin_width/2];
ks2stat_all = nan(1000,length(timewindows)-1);

%%% Low ripples
all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
for nbin = 1:length(timewindows)-1
    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(1,nbin)&context_modulation_all.timebin<timewindows(2,nbin)),2,'omitnan');

    T1_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0);
    T2_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0);
    num_neurons = min([length(T1_FR_dff) length(T2_FR_dff)]);


    parfor iBoot = 1:1000
        temp1 = [];temp2 = [];
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index1 = randsample(s, 1:length(T1_FR_dff), num_neurons, true);
        temp1 = T1_FR_dff(index1);
        index2 = randsample(s, 1:length(T2_FR_dff), num_neurons, true);
        temp2 = T2_FR_dff(index2);
        [h,pval,ks2stat_boot] = kstest2(temp1,temp2);
        ks2stat_all(iBoot,nbin) = ks2stat_boot;
    end
end
ks2stat_low_ripples = ks2stat_all;

%%%% High ripples
all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
ks2stat_all = nan(1000,length(timewindows)-1);

for nbin = 1:length(timewindows)-1
    POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>timewindows(1,nbin)&context_modulation_all.timebin<timewindows(2,nbin)),2,'omitnan');

    T1_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff>0);
    T2_FR_dff = POST_ripple_FR_diff(contains(regions_all,'V1')&z_FR_track_diff<0);
    num_neurons = min([length(T1_FR_dff) length(T2_FR_dff)]);


    parfor iBoot = 1:1000
        temp1 = [];temp2 = [];
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index1 = randsample(s, 1:length(T1_FR_dff), num_neurons, true);
        temp1 = T1_FR_dff(index1);
        index2 = randsample(s, 1:length(T2_FR_dff), num_neurons, true);
        temp2 = T2_FR_dff(index2);
        [h,pval,ks2stat_boot] = kstest2(temp1,temp2);
        ks2stat_all(iBoot,nbin) = ks2stat_boot;
    end
end
ks2stat_high_ripples = ks2stat_all;



nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 KS difference time series (high vs low ripple)';
nfig.Position = [ 1150         322         363         282];
timewindows = [bin_centers - bin_width/2; bin_centers + bin_width/2];
% x = timewindows(1:end-1)+mean(diff(timewindows))/2;
x = mean(timewindows);
x = x(1:end-1)+mean(diff(x))/2;
colour_lines = [158,202,225;33,113,181]/256;% two blue 

nexttile
% x = timewindows(1:end-1)+mean(diff(timewindows))/2;
plot(x,mean(ks2stat_low_ripples),'Color',colour_lines(1,:));hold on;
ci_low  = prctile(ks2stat_low_ripples, 2.5, 1);
ci_high = prctile(ks2stat_low_ripples, 97.5, 1);

% Fill 95% CI
F(1) = fill([x, fliplr(x)], [ci_low, fliplr(ci_high)], ...
    colour_lines(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

plot(x,mean(ks2stat_high_ripples),'Color',colour_lines(end,:));hold on;
ci_low  = prctile(ks2stat_high_ripples, 2.5, 1);
ci_high = prctile(ks2stat_high_ripples, 97.5, 1);

% Fill 95% CI
F(2)=fill([x, fliplr(x)], [ci_low, fliplr(ci_high)], ...
    colour_lines(end,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
set(gca,'TickDir','out','Box','off','FontSize',12)
xticks([-0.5 -0.25 0 0.25 0.5])
xline(0,'r')

legend(F(1:2),{'low ripples','high ripples'},'box','off')
xlabel('Time relative to ripple onset (s)')
ylabel('Maximum empirical cumulative distribution difference');

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])
