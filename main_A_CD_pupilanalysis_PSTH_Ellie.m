%%% For analysis of trial-by-trial pupil size vs "omission responses". ERB 2026
% Depths of V1 L5 and CA1 are both based on best regional power per PSD analysis (for which, use PSD_analysis_UnProcessedLFP_ellie).
% V1 L4 depth is identified via CSD analysis (for which, run CSD_Gratings_afterFILT_ellie or CSD_GAVNIKstims_afterFILT_ellie

%% setting metrics to screen good clusters
clear all
addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

% 1/5 Choose your probe depth of interest 
depth_for_analysis = 'V1'; % choose 'L4' or 'V1' or 'CA1' or 'Sub_CA1' or 
stim_duration = 150;
subjects_dir = 'V:\Ellie\DATA\SUBJECTS';
SUBJECTS = {'M00014'};
params = create_cluster_selection_params('sorting_option','ellie');
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_type = 'GAVNIK_A_CD'; % 'GAVNIK_A_CD' 'GAVNIK_ABCD_1' 
plot_choice = 'aggregate'; % curated 'single_units' or in 'aggregate' or uncurated 'MUA'; MUA includes all clusters from kilosort, unfiltered
plot_type = 'FR'; % 'FR' firing rate or 'raster' or 'struct' (for no plotting but output of metrics).
sliced_plot_option = 'no'; % 'yes' if you want to plot traces by groups of 40 trials to look for changes during the session
z_method = 'per_neuron'; 
z_score_period = 'stim_session'; % 'none' or 'stim_session' excludes variable greyscreen periods before and after stim paradigm
if contains(Stimulus_type, '_A_CD')
    stimulus_color = [1.0,  0.0,  0.0];  %  red - GAVNIK_A_CD  
elseif contains(Stimulus_type, '_E_CD') 
    stimulus_color = [0.0,  1.0,  0.0];  % green - GAVNIK_E_CD  
elseif contains(Stimulus_type, '_ABCD')    
    stimulus_color = [0.0,  0.0,  1.0];  % blue - GAVNIK_ABCD
end  

% 3/5 files will be saved here in the cd
analysis_dir = fullfile(subjects_dir, SUBJECTS{1}, ...
    'analysis', '20250321', Stimulus_type);
cd(analysis_dir);

%% SET THIS 4/5***** For NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
probe_type = 0; % NPX1.0 is type 0, NPX2.0 is type 1.

% 5/5 
for nsession = 9 % row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    
    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        subject_number = session_info(n).probe(1).SUBJECT;
        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters = clusters_ks4;
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'earliest_V1sink_CSD.mat'))
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'depths_from_PSD.mat'))
        files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel
        file_to_load = fullfile(options.EPHYS_DATAPATH, files(1).name); %load() does not accept wildcards like *
        load(file_to_load);

        
        if probe_type == 0
            layerfour_channels = earliest_V1sink_CSD.overall_best_halfmax_channel;
            layerfour_channels = layerfour_channels(:);
            shank_ids = ones(size(layerfour_channels)); % dummy shank ID = 1
        
            assert(isfield(earliest_V1sink_CSD, 'overall_best_halfmax_depth'), ...
                'Missing overall_best_halfmax_depth in earliest_V1sink_CSD');
        
            earliest_V1sink_CSD = struct( ...
                'shank_id', 1, ...
                'best_channel_this_shank', layerfour_channels(1), ...
                'best_depth_this_shank', earliest_V1sink_CSD(1).overall_best_halfmax_depth ...
            );
        elseif probe_type == 1
            % NPX2.0 → multiple shanks
            % earliest_V1sink_CSD is a struct array, extract fields safely
            layerfour_channels = [earliest_V1sink_CSD.best_channel_this_shank]; 
            shank_ids = [earliest_V1sink_CSD.shank_id]; 
            % ensure column vectors
            layerfour_channels = layerfour_channels(:); 
            shank_ids = shank_ids(:);
        end

        
        % --- Identify shanks and their layer 4 channels ---
        [unique_shanks, ia, ~] = unique(shank_ids, 'stable');
        best_channels_per_shank = layerfour_channels(ia); % one per shank
        
        % --- Preallocate shank-specific depth ranges ---
        L4_depth_range = cell(numel(unique_shanks), 1);
        V1_depth_range = cell(numel(unique_shanks), 1);
        CA1_depth_range = cell(numel(unique_shanks), 1);
        Sub_CA1_depth_range = cell(numel(unique_shanks), 1);
               
        % --- Compute per-shank depth ranges using CSD and PSD data ---
        %% --- Normalise PSD struct for NPX1.0 to look like NPX2.0 ---
        if probe_type == 0
            % If PSD fields are not shank-wrapped, wrap them into shank1
            if ~isfield(depths_from_PSD, 'shank1')
                tmp = depths_from_PSD;
                depths_from_PSD = struct();
                depths_from_PSD.shank1 = tmp;
            end
        end
        
        for iShank = 1:numel(unique_shanks)
            this_shank   = unique_shanks(iShank);
            this_channel = best_channels_per_shank(iShank);
            
            this_shank_name = ['shank', num2str(this_shank)];  % creates e.g. 'shank2'
            Brain_surface_depth = depths_from_PSD.(this_shank_name).surface_depth_PSD;
            L4_channel_depth    = earliest_V1sink_CSD(find([earliest_V1sink_CSD.shank_id] == this_shank, 1)).best_depth_this_shank;
            L5_depth            = depths_from_PSD.(this_shank_name).L5_depth_PSD;
            CA1_depth           = depths_from_PSD.(this_shank_name).CA1_depth_PSD;
        
            L4_depth_range{iShank}   = [L4_channel_depth - 60, L4_channel_depth + 60]; % giving electrodes the full extent of 120um inclusive 
            % (hence measuring spiking over depth range greater than 120um. As 60 is divisible by 15 and 10, this will give the same effective range for NPX1.0 (which has staggered electrodes every 10um down the shank) and NPX2.0 (which has electrode rows every 15um down each shank)
            V1_depth_range{iShank}   = [L5_depth - 330, L5_depth + 700];
            CA1_depth_range{iShank}  = [CA1_depth - 150, CA1_depth + 150];
            Sub_CA1_depth_range{iShank} = [min(CA1_depth_range{iShank}) - 1000, min(CA1_depth_range{iShank})];  
        end
             
        all_orientations = unique(Task_info.stim_orientation); % uniqe values sorted in ascending order
               
        params = create_cluster_selection_params('sorting_option','ellie');
        psthBinSize = 0.01; % but use 1ms for raster plots
        
        switch depth_for_analysis
            case 'L4' 
                depth_ranges = L4_depth_range;
            case 'V1'
                depth_ranges = V1_depth_range;
            case 'CA1' 
                depth_ranges = CA1_depth_range;
            case 'Sub_CA1'
                depth_ranges = Sub_CA1_depth_range;
        end

        if contains(z_score_period, 'stim_session') % excludes variable grey screen period before and after the stimulus paradigm ran
            time_edges = (min(Task_info.stim_onset) - 2):psthBinSize:(max(Task_info.stim_onset) + 2);
        end             
        
        
        if (contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD') || contains(Stimulus_type, 'GAVNIK_ABCD')) && contains(plot_choice, 'aggregate') && contains(plot_type, 'FR')
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);

                    cluster_channels = clusters(np).peak_channel(sc.cluster_id);
                    cluster_depths = ycoords(cluster_channels); % get depths of selected clusters from ycoords (peak_depths from SI are 15 microns different..)
    
                    % Map each cluster to its shank
                    cluster_shanks = kcoords(cluster_channels);
    
                    % Initialize container for depth-filtered cluster IDs
                    depth_cluster_ids = [];
            
                    % Loop over each shank
                    for iShank = 1:numel(unique_shanks)
                        this_shank = unique_shanks(iShank);
                        this_range = depth_ranges{iShank};
            
                        % Logical mask for clusters on this shank and within this shank’s depth range
                        shank_mask = cluster_shanks == this_shank;
                        depth_mask = cluster_depths >= min(this_range) & cluster_depths <= max(this_range);
            
                        % Keep only clusters satisfying both conditions
                        keep_mask = shank_mask & depth_mask;
            
                        % Append these cluster IDs
                        depth_cluster_ids = [depth_cluster_ids; sc.cluster_id(keep_mask(:))];
                    end

                    % Keep only spike times and IDs for clusters at selected depth
                    depth_selected_clusters(np).cluster_id = depth_cluster_ids;
                    depth_selected_clusters(np).spike_times = sc.spike_times(ismember(sc.spike_id, depth_cluster_ids));
                    depth_selected_clusters(np).spike_id = sc.spike_id(ismember(sc.spike_id, depth_cluster_ids));
                end
                
                % Combine spike times for all clusters at selected depth (i.e., for curated MUA)
                all_spike_times = [];
                all_spike_ids = [];
                for np = 1:length(depth_selected_clusters)
                    all_spike_times = [all_spike_times; depth_selected_clusters(np).spike_times];
                    all_spike_ids = [all_spike_ids; depth_selected_clusters(np).spike_id];
                end

                if strcmp(z_method, 'per_neuron') && (contains(z_score_period, 'stim_session'))
                    % define each unit's baseline distribution; mean and sd of its firing across the entire session
                    unique_units = unique(all_spike_ids);
                    nUnits = numel(unique_units);              
                    unit_baseline_mean = nan(nUnits,1);
                    unit_baseline_std  = nan(nUnits,1);
                
                    for u = 1:nUnits
                        this_unit = unique_units(u);
                        % Get spikes from this neuron
                        unit_spike_times = all_spike_times(all_spike_ids == this_unit);
                        % Bin this neuron's spikes across the stim session
                        unit_counts = histcounts(unit_spike_times, time_edges); %number of spikes fired in each time bin across the recording session (g number, should be about 540s for Gavornik 150ms protocol)
                
                        unit_baseline_mean(u) = mean(unit_counts);
                        unit_baseline_std(u)  = std(unit_counts);
                    end
                end

                
                ordered_oris = unique(Task_info.stim_orientation, 'stable');               
                % make a figure of the mean trace for sense-checking - ensure it matches up with the corresponding GAVNIK_ABCD vs A_CD vs E_CD
                
                fig = figure;
                fig.Name = sprintf('%s: Aggregate activity: %s', Stimulus_type, depth_for_analysis);
                fig.Position = [114 90 770 650];

                % Define grey intervals
                if (contains(Stimulus_type, 'GAVNIK250_A_CD'))
                    grey_intervals = [-0.5 0; 1.0 2.0];
                    element2_onset = 0.25;
                else
                    grey_intervals = [-0.5 0; 0.6 2.0];
                    element2_onset = 0.15;
                end    

                for i = 1:size(grey_intervals,1)
                    h = xregion(grey_intervals(i,1), grey_intervals(i,2), ...
                        FaceColor=[0.7 0.7 0.7], ...
                        FaceAlpha=0.3);
                    h.HandleVisibility = 'off';
                end 
                
                ori = 1; % plot from onset of first stimulus in sequence
                stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));

                %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(all_spike_times, stim_onsets, [-0.3 1.75], psthBinSize);
                psth_window = [-0.30 1.75];
                
                if strcmp(z_method, 'per_neuron')
                    % express each unit's stimulus-evoked response versus its own baseline mean and sd across the session
                    unique_units = unique(all_spike_ids);
                    nUnits = numel(unique_units);                    
                    
                    % Get the binned data from the first neuron to determine dimensions
                    first_unit = unique_units(1);
                    first_unit_spikes = all_spike_times(all_spike_ids == first_unit);
                
                    [~, bins, ~, ~, ~, first_binnedArray_u] = ...
                        psthAndBA(first_unit_spikes, stim_onsets, psth_window, psthBinSize);
                
                    nTrials = size(first_binnedArray_u, 1);
                    nTimeBins = size(first_binnedArray_u, 2);
                
                    % Preallocate: neurons x trials x time bins
                    z_binnedArray = nan(nUnits, nTrials, nTimeBins);

                    for u = 1:nUnits
                        this_unit = unique_units(u);
                
                        if u == 1
                            binnedArray_u = first_binnedArray_u;
                        else
                            unit_mask = all_spike_ids == this_unit;
                            unit_spike_times = all_spike_times(unit_mask);
                    
                            [~, ~, ~, ~, ~, binnedArray_u] = ...
                                psthAndBA(unit_spike_times, stim_onsets, ...
                                          psth_window, psthBinSize);
                        end
                        mu = unit_baseline_mean(u);
                        sd = unit_baseline_std(u);
                
                        if sd == 0 || isnan(sd)
                            continue
                        end
                        % Z-score every trial using this neuron's own mean and SD
                        z_u = (binnedArray_u - mu) ./ sd; % express each unit's stimulus-evoked response versus its own baseline mean and sd across the session
                        % Store all trials
                        z_binnedArray(u,:,:) = z_u;
                    end               
                end
                
                % Mean z-scored firing rate across neurons for each trial
                mean_z_trace_by_trial = squeeze(mean(z_binnedArray, 1, 'omitnan'));
                
                % Mean z-scored firing rate across trials for sense-checking
                mean_z_trace = mean(mean_z_trace_by_trial, 1, 'omitnan');

                hold on;
                plot(bins, mean_z_trace, 'Color', stimulus_color)
                h = xline(element2_onset, 'k', 'Element 2 onset');
                h.LabelHorizontalAlignment = 'left';
                xline(element2_onset + 0.03, 'k--', '+30 ms');
                xline(element2_onset + 0.08, 'k--', '+80 ms');
                title(sprintf('%s %s: Aggregate activity (mean_z_trace) for sense-checking: %s', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
                hold off;
                % legend(flipud(findobj(gca,'-property','DisplayName')), 'Location', 'northeast', 'Interpreter', 'none');
                % hold on;               
                % xline(0, 'k', (sprintf('A %d%s onset', round(ordered_oris(1)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                % xline(0.15, 'k', (sprintf('grey onset')), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                % xline(0.30, 'k', (sprintf('C %d%s onset', round(ordered_oris(3)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                % xline(0.45, 'k', (sprintf('D %d%s onset', round(ordered_oris(4)), char(176))), 'LabelVerticalAlignment','top', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off', 'FontSize', 24);
                

                % calculate "omission response" for each trial
                % z_binnedArray = neurons x trials x time bins
                
                % Find bins within the 30–80 ms post-element-2 window
                peak_idx = bins >= element2_onset + 0.03 & ...
                           bins <= element2_onset + 0.08;
                
                % Number of neurons and trials
                nUnits  = size(z_binnedArray, 1);
                nTrials = size(z_binnedArray, 2);
                
                % Store element2 "omission" response for each neuron and trial
                element2_response_by_neuron = nan(nUnits, nTrials);
                
                for u = 1:nUnits
                
                    % Extract this neuron's trial x time matrix
                    z_u = squeeze(z_binnedArray(u,:,:));
                
                    % Interpolate firing at the exact element-2 onset for every trial
                    FR_at_element2_onset = nan(nTrials,1);
                
                    for trial = 1:nTrials
                        FR_at_element2_onset(trial) = interp1( ...
                            bins, z_u(trial,:), element2_onset, 'linear');
                    end
                
                    % Find peak response during 30–80 ms after element-2 onset
                    peak_FR = max(z_u(:,peak_idx), [], 2);
                
                    % Element2 "omission" response = peak - firing at element-2 onset
                    element2_response_by_neuron(u,:) = ...
                        peak_FR - FR_at_element2_onset;
                
                end
                
                % Average across neurons so that every neuron contributes equally
                element2_response = mean(element2_response_by_neuron, 1, 'omitnan');
                element2_response = element2_response(:);

                disp(size(z_binnedArray))
                disp(size(mean_z_trace_by_trial))
                disp(size(element2_response))
                disp(length(bins))

                % Pupil data
                pupil_area = Behaviour.pupil_area_sixtyHz;
                pupil_time = Behaviour.sglxTime;
                
                disp(size(pupil_area))
                disp(size(pupil_time))
                disp([pupil_time(1), pupil_time(end)])
                disp([min(stim_onsets), max(stim_onsets)])
                

                % Pupil area at onset of A and element2 for each trial
                pupil_at_A_onset = nan(nTrials,1);
                pupil_at_element2_onset = nan(nTrials,1);
                % interpolate pupil_area from the two values at the two pupil_times that straddle the stim_onset time
                % pupil is only sampled at 60Hz (every 16.7ms)
                for trial = 1:nTrials
                    % Pupil area at A onset
                    pupil_at_A_onset(trial) = interp1( ...
                        pupil_time, pupil_area, ...
                        stim_onsets(trial), ...
                        'linear');
                    % Pupil area at element 2 onset
                    pupil_at_element2_onset(trial) = interp1( ...
                        pupil_time, pupil_area, ...
                        stim_onsets(trial) + element2_onset, ...
                        'linear');
                end
                
                % Change in pupil area from A onset to element 2 onset
                % Negative = pupil became smaller
                % Positive = pupil became larger
                pupil_change = pupil_at_element2_onset - pupil_at_A_onset;

                disp(size(pupil_at_A_onset))
                disp(size(element2_response))
                
                disp([pupil_at_A_onset(1:10), element2_response(1:10)])
                
                %% ---------------- PLOT 1: Pupil at A onset ----------------
                valid_trials = isfinite(element2_response) & isfinite(pupil_at_A_onset);

                x = element2_response(valid_trials);
                y = pupil_at_A_onset(valid_trials);
                
                fig = figure;
                scatter(x, y, 45, stimulus_color, 'filled');
                
                xlabel('Peak element2 spiking vs spiking at element2 onset');
                ylabel('Pupil area at A onset, arb. units');
                sgtitle(sprintf('%s - %s: Pupil area at A onset vs. element2 spiking bounce: %s', subject_number, Stimulus_type, depth_for_analysis), 'Interpreter', 'none');
                xlim([-0.5 5]);
                ylim([0 3000]);

                grid off;
                set(gca, 'FontSize', 24);  % Tick labels font size
                box off;
                
                %% Correlation and line of best fit

                [r, p_corr] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
                
                coeffs = polyfit(x, y, 1);
                x_fit = linspace(min(x), max(x), 100);
                y_fit = polyval(coeffs, x_fit);
                
                hold on;
                plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
                hold off;
                
                fprintf('Pearson correlation: r = %.3f, p = %.4g\n', r, p_corr);
                
                % Print correlation statistics on plot
                n_datapoints = numel(x);
                
                text(0.05, 0.95, ...
                    sprintf('Pearson r = %.2f\np = %.3g\nn = %d', r, p_corr, n_datapoints), ...
                    'Units', 'normalized', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize', 12);
                                
                box off;
                
                % Calculate means and SEMs
                mean_omission = mean(element2_response, 'omitnan');
                %sem_omission  = std(element2_response, 'omitnan') / sqrt(sum(~isnan(element2_response)));
                
                mean_pupil = mean(pupil_at_A_onset, 'omitnan');
                %sem_pupil  = std(pupil_at_A_onset, 'omitnan') / sqrt(sum(~isnan(pupil_at_A_onset)));
                hold on;
                plot(mean_omission, mean_pupil, ...
                    'ko', ...
                    'MarkerFaceColor', 'k', ...
                    'MarkerSize', 8);

                % Also save as .fig
                fig_filename = sprintf('%s %s Pupil at A onset vs element2 response %s.fig', subject_number, depth_for_analysis, Stimulus_type);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
                hold off;   
                %% ---------------- PLOT 2: Pupil at element 2 onset ----------------

                valid_trials = isfinite(element2_response) & ...
                               isfinite(pupil_at_element2_onset);
                
                x = element2_response(valid_trials);
                y = pupil_at_element2_onset(valid_trials);
                
                fig = figure;
                scatter(x, y, 45, stimulus_color, 'filled');
                
                xlabel('Peak element2 spiking vs spiking at element2 onset');
                ylabel('Pupil area at element 2 onset, arb. units');
                xlim([-0.5 5]);
                ylim([0 3000]);

                grid off;
                set(gca, 'FontSize', 24);
                box off;
                
                %% Correlation and line of best fit
                
                [r, p_corr] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
                
                coeffs = polyfit(x, y, 1);
                x_fit = linspace(min(x), max(x), 100);
                y_fit = polyval(coeffs, x_fit);
                
                hold on;
                plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
                
                fprintf('Element 2 pupil: Pearson r = %.3f, p = %.4g\n', r, p_corr);
                
                n_datapoints = numel(x);

                text(0.05, 0.95, ...
                    sprintf('Pearson r = %.2f\np = %.3g\nn = %d', r, p_corr, n_datapoints), ...
                    'Units', 'normalized', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize', 12);
                

                %% Save mouse-level correlation result
        
                element2_pupil_vsSpiking = struct();
                
                % --- Subject/session information ---
                element2_pupil_vsSpiking.subject_number = subject_number;
                element2_pupil_vsSpiking.Stimulus_type = Stimulus_type;
                element2_pupil_vsSpiking.depth_for_analysis = depth_for_analysis;
                element2_pupil_vsSpiking.z_method = z_method;
                element2_pupil_vsSpiking.z_score_period = z_score_period;
                
                % --- Correlation ---
                element2_pupil_vsSpiking.pearson_r = r;
                element2_pupil_vsSpiking.pearson_p = p_corr;
                element2_pupil_vsSpiking.n_trials = n_datapoints;
                
                % --- Pupil measurement ---
                element2_pupil_vsSpiking.pupil_measurement = 'Pupil area at element 2 onset';
                element2_pupil_vsSpiking.pupil_time = element2_onset;
                
                % --- Spiking omission response ---
                element2_pupil_vsSpiking.Spiking_measurement = ...
                    'Peak spiking 30-80 ms after element 2 onset less spiking at element 2 onset';
                
                element2_pupil_vsSpiking.element2_onset = element2_onset;
                
                if stim_duration == 150
                    element2_pupil_vsSpiking.Spiking_peak_window = [0.180 0.230];
                elseif stim_duration == 250
                    element2_pupil_vsSpiking.Spiking_peak_window = [0.280 0.330];
                end
                
                element2_pupil_vsSpiking.pupil_units = 'arbitrary units';
                
                % --- Trial-level data used for the correlation ---
                element2_pupil_vsSpiking.element2_response = x;
                element2_pupil_vsSpiking.pupil_at_element2_onset = y;
                
                % --- Analysis information ---
                element2_pupil_vsSpiking.correlation_type = 'Pearson';
                element2_pupil_vsSpiking.analysis_date = datetime("now");
                
                %% Save result struct
        
                result_filename = sprintf('element2_pupil_vsSpiking_%s_%s.mat', ...
                    subject_number, Stimulus_type);
        
                save(fullfile(subjects_dir, result_filename), 'element2_pupil_vsSpiking');
                       
                fprintf('Saved correlation result to:\n%s\n', ...
                    fullfile(subjects_dir, result_filename));



                % Mean ± SEM
                mean_omission = mean(x, 'omitnan');
                %sem_omission = std(x, 'omitnan') / sqrt(sum(isfinite(x)));
                
                mean_pupil = mean(y, 'omitnan');
                %sem_pupil = std(y, 'omitnan') / sqrt(sum(isfinite(y)));
                
                plot(mean_omission, mean_pupil, ...
                    'ko', ...
                    'MarkerFaceColor', 'k', ...
                    'MarkerSize', 8);
                
                title(sprintf('%s - %s: Pupil area at element 2 onset vs. element2 spiking bounce', ...
                    subject_number, Stimulus_type), ...
                    'Interpreter', 'none');

                % Also save as .fig
                fig_filename = sprintf('%s %s Pupil at element2 onset vs element2 response %s.fig', subject_number, depth_for_analysis, Stimulus_type);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
                hold off; 
                
                
                
                %% ---------------- PLOT 3: Change in pupil area ----------------
                
                valid_trials = isfinite(element2_response) & ...
                               isfinite(pupil_change);
                
                x = element2_response(valid_trials);
                y = pupil_change(valid_trials);
                
                fig = figure;
                scatter(x, y, 45, stimulus_color, 'filled');
                
                xlabel('Peak element2 spiking vs spiking at element2 onset');
                ylabel('\Delta pupil area (element 2 onset - A onset)');
                xlim([-0.5 5]);
                
                grid off;
                set(gca, 'FontSize', 24);
                box off;
                
                %% Correlation and line of best fit
                
                [r, p_corr] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
                
                coeffs = polyfit(x, y, 1);
                x_fit = linspace(min(x), max(x), 100);
                y_fit = polyval(coeffs, x_fit);
                
                hold on;
                plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
                
                fprintf('Pupil change: Pearson r = %.3f, p = %.4g\n', r, p_corr);
                
                n_datapoints = numel(x);

                text(0.05, 0.95, ...
                    sprintf('Pearson r = %.2f\np = %.3g\nn = %d', r, p_corr, n_datapoints), ...
                    'Units', 'normalized', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize', 12);                
                % Mean ± SEM
                mean_omission = mean(x, 'omitnan');
                %sem_omission = std(x, 'omitnan') / sqrt(sum(isfinite(x)));
                
                mean_pupil_change = mean(y, 'omitnan');
                %sem_pupil_change = std(y, 'omitnan') / sqrt(sum(isfinite(y)));
                
                plot(mean_omission, mean_pupil_change, ...
                    'ko', ...
                    'MarkerFaceColor', 'k', ...
                    'MarkerSize', 8);
                
                title(sprintf('%s - %s: Pupil change vs. element2 spiking bounce', ...
                    subject_number, Stimulus_type), ...
                    'Interpreter', 'none');
                % Also save as .fig
                fig_filename = sprintf('%s %s Pupil at element2 onset less at Aonset vs element2 response %s.fig', subject_number, depth_for_analysis, Stimulus_type);
                fig_save_path = fullfile(pwd, fig_filename);
                savefig(fig, fig_save_path);
                hold off;

                %% Sense check: pupil area across entire session

                fig = figure;
                plot(pupil_time, pupil_area, 'k-', 'LineWidth', 1);
                hold on;
                
                % Overlay A stimulus onsets
                for i = 1:length(stim_onsets)
                    xline(stim_onsets(i), 'r-', 'HandleVisibility', 'off');
                end
                
                xlabel('Time (s)');
                ylabel('Pupil area (a.u.)');
                title(sprintf('%s %s: Pupil area across session', subject_number, Stimulus_type), ...
                    'Interpreter', 'none');
                
                set(gca, 'FontSize', 18);
                box off;
            end
        end



        
        if (contains(Stimulus_type, 'GAVNIK_A_CD') || contains(Stimulus_type, 'GAVNIK_E_CD')) &&...
           contains(plot_choice, 'single_units') 
            for nprobe = 1:length(clusters)
                selected_clusters(nprobe) = select_clusters(clusters(nprobe),params); %only look at good clusters, which pass the set parameters
                
                depth_selected_clusters = selected_clusters; % initialize with same structure
                for np = 1:length(clusters)
                    sc = selected_clusters(np);
                    cluster_channels = clusters(np).peak_channel(sc.cluster_id);
                    cluster_depths = ycoords(cluster_channels); % get depths of selected clusters from ycoords (peak_depths from SI are 15 microns different..)
    
                    % Map each cluster to its shank
                    cluster_shanks = kcoords(cluster_channels);
    
                    % Initialize container for depth-filtered cluster IDs
                    depth_cluster_ids = [];
            
                    % Loop over each shank
                    for iShank = 1:numel(unique_shanks)
                        this_shank = unique_shanks(iShank);
                        this_range = depth_ranges{iShank};
            
                        % Logical mask for clusters on this shank and within this shank’s depth range
                        shank_mask = cluster_shanks == this_shank;
                        depth_mask = cluster_depths >= min(this_range) & cluster_depths <= max(this_range);
            
                        % Keep only clusters satisfying both conditions
                        keep_mask = shank_mask & depth_mask;
            
                        % Append these cluster IDs
                        depth_cluster_ids = [depth_cluster_ids; sc.cluster_id(keep_mask(:))];
                    end

                    % Keep only spike times and IDs for clusters at selected depth
                    depth_selected_clusters(np).cluster_id = depth_cluster_ids;
                    depth_selected_clusters(np).spike_times = sc.spike_times(ismember(sc.spike_id, depth_cluster_ids));
                    depth_selected_clusters(np).spike_id = sc.spike_id(ismember(sc.spike_id, depth_cluster_ids));
                end

                cluster_id = depth_selected_clusters(nprobe).cluster_id; % cluster_ids of units which pass the set parameters [NB these are one count higher than per zero-based pythonic SI output cluster IDs...]
                
                baseline_window = [0 30]; % in seconds - first 30s of recording is grey screen - can use for z-scoring

                % Define time_edges depending on z_score_period
                if contains(z_score_period, 'entire_session')
                    time_edges = 0:psthBinSize:max(depth_selected_clusters(nprobe).spike_times);
                elseif contains(z_score_period, 'first30secs')
                    time_edges = baseline_window(1):psthBinSize:baseline_window(2);
                end
                
                ori_response = [];
                z_binnedArray = [];

                for nCluster = 1:length(cluster_id) % loop through each good cluster at selected depth
                
                    spike_times_this_cluster = depth_selected_clusters(nprobe).spike_times(depth_selected_clusters(nprobe).spike_id == cluster_id(nCluster)); % extract the spike times for this cluster
                    %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset, [-0.150 0.30], 0.001); % gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
             
                    % Compute appropriate histogram counts for z-scoring
                    if contains(z_score_period, 'entire_session') || contains(z_score_period, 'stim_session')
                        zscore_counts = histcounts(spike_times_this_cluster, time_edges);
                    elseif contains(z_score_period, 'first30secs')
                        baseline_spikes = spike_times_this_cluster(spike_times_this_cluster >= baseline_window(1) & spike_times_this_cluster <= baseline_window(2));
                        zscore_counts = histcounts(baseline_spikes, time_edges);
                    end

                    fig(nCluster)=figure; % open a figure window
                    fig(nCluster).Name=sprintf('%s Grating responses Cluster %i', Stimulus_type, cluster_id(nCluster)); %overall figure title includes the cluster_id (one count higher than zero-based SI output)
                    fig(nCluster).Position = [114 90 770 650]; % sets the size of the figure window
                
                    tiledlayout(5,1); % vertical layout (5 rows × 1 column)
                    % Extract ordered orientations based on presentation sequence
                    ordered_oris = unique(Task_info.stim_orientation, 'stable'); % radians for TRAIN but degrees for GAVNIK
                    
                    for ori = 1:length(ordered_oris)  

                        if contains(plot_type, 'raster')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.02 0.17], psthBinSize/10); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            nexttile % opens a subplot tile
                            plot(rasterX,rasterY,'k','LineWidth',1) % plot the spiking raster for the current orientation
                            xline(0,'r',LineWidth=1) % put a vertical red line at time zero (stimulus onset)
                            xlim([-0.02 0.17])
                            xticks([0, 0.05, 0.1, 0.15]);
                            ylim([0 sum(Task_info.stim_orientation==ordered_oris(ori))]) % set the max y coord to be the total # trials with this orientation
                            ylabel('Trial');
                            % imagesc(binnedArray)
                            % xticks([1.5 10.5 15.5 20.5 30.5])
                            % xticklabels([-0.150 -0.050 0 0.050 0.150])
                            % xline(15.5,'r',LineWidth=1)
                            % colorbar
                            % colormap(flipud(gray))
                            
                            if ori == 2
                                title('Grey screen'); % in the GAVNIK_A_CD and GAVNIK_E_CD conditions, the second stimulus is absent
                            else
                                title(sprintf('Orientation %d%s', round(ordered_oris(ori)), char(176)));
                            end
                            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',24)
                        end
                        
                        
                        if contains(plot_type, 'FR')
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.02 0.17], psthBinSize); % for this orientation, gets the rasterplot coordinates and binnedArray (matrix of spike counts per timebin), with a time window of -150ms to +150ms around stimulus onset
                            
                            mean_trace = mean(binnedArray, 1);
                            z_trace = (mean_trace - mean(zscore_counts)) / std(zscore_counts);
                            nexttile;
                            plot(bins, z_trace, 'k', 'LineWidth', 1.5);
                            xline(0, 'r', 'LineWidth', 1);
                            xlim([-0.02 0.17]);
                            xticks([0, 0.05, 0.1, 0.15]);
                            ylabel('Z-scored FR');
                            if ori == 2
                                title('Grey screen'); % in the GAVNIK_A_CD and GAVNIK_E_CD conditions, the second stimulus is absent
                            else
                                title(sprintf('Orientation %d%s', round(ordered_oris(ori)), char(176)));
                            end
                            
                            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                        end

                        %[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_this_cluster,  Task_info.stim_onset(Task_info.stim_orientation==ordered_oris(ori)), [-0.150 0.30], psthBinSize);
                        %z_binnedArray{nCluster}{ori} = (binnedArray - mean(spikecounts_cluster))./std(spikecounts_cluster); %z-score normalisation of the binned Array for this orientation: subtracts the mean and divides by the s.d. of the spike count histogram for the whole session
                        % hold on
                        % plot(bins,mean(z_binnedArray));
             
                        % time_selected 
                    
                    end
                    
                    % === Extra tile for post-stimulus activity for 4th orientation ===
                    if length(ordered_oris) >= 4 && contains(plot_type, 'FR')
                        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(4));
                        [~, bins_long, ~, ~, ~, binnedArray_long] = psthAndBA(spike_times_this_cluster, stim_onsets, [-0.02 0.75], psthBinSize);
                        mean_trace_long = mean(binnedArray_long, 1);
                        z_trace_long = (mean_trace_long - mean(zscore_counts)) / std(zscore_counts);

                        nexttile;
                        plot(bins_long, z_trace_long, 'k', 'LineWidth', 1.5);
                        xline(0, 'r', 'LineWidth', 1);
                        xlim([-0.02 0.75]);
                        xticks(0:0.05:0.75);
                        ylabel('Z-scored FR');
                        
                        title(sprintf('Orientation %d%s (Extended to show post-stimulus oscillations)', round(ordered_oris(4)), char(176)));
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                    end
                    
                    if length(ordered_oris) >= 4 && contains(plot_type, 'raster')
                        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(4));
    
                        [~, ~, rasterX_long, rasterY_long, ~, ~] = psthAndBA(spike_times_this_cluster, stim_onsets, [-0.02 0.75], psthBinSize/10);
    
                        nexttile;
                        plot(rasterX_long, rasterY_long, 'k', 'LineWidth', 1);
                        xline(0, 'r', 'LineWidth', 1);
                        xlim([-0.02 0.75]);
                        ylim([0 sum(Task_info.stim_orientation == ordered_oris(4))]);
                        xticks(0:0.05:0.75);
                        
                        title(sprintf('Orientation %d%s (Extended to show post-stimulus oscillations)', round(ordered_oris(4)), char(176)));
                        ylabel('Trial');
                        set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 24);
                    end


                    % Sanitize Stimulus_type for filenames
                    safeStimulusType = regexprep(Stimulus_type, '[:\\/*?"<>| ]', '_');
                    cluster_depth = clusters(nprobe).peak_depth(cluster_id(nCluster));
                    sgtitle(sprintf('%s - %s: Response of %s Cluster %i (%.0f µm)', subject_number, Stimulus_type, depth_for_analysis, cluster_id(nCluster), cluster_depth), 'Interpreter', 'None');

                    if contains(plot_type, 'raster')
                        exportgraphics(fig(nCluster), sprintf('%s_%s_Cluster_%i_raster.pdf', safeStimulusType, depth_for_analysis, cluster_id(nCluster)));
                    end

                    if contains(plot_type, 'FR')
                        exportgraphics(fig(nCluster), sprintf('%s_%s_Cluster_%i_FR.pdf', safeStimulusType, depth_for_analysis, cluster_id(nCluster)));
                    end
                end           
            end
        end
    end
end

