%% 0. Load
load('V1_HPC_reactivation_coherence_lme1_raw.mat');


%% 1. Outlier Removal
% Define the continuous variables that are vulnerable to extreme noise
outlierVars = {'RipplePower', 'SpindlePower_Match', 'SpindlePower_NonMatch', ...
               'HPC_logodds', 'V1_logodds_PRE', 'V1_logodds'};

% Find outliers (trimming the extreme 0.1% tails)
outlierMask = isoutlier(tbl{:, outlierVars}, 'percentiles', [0.01, 99.99]);
rowsWithOutliers = any(outlierMask, 2);

% Create the clean table
cleanTbl = tbl(~rowsWithOutliers, :);
fprintf('Removed %d rows containing extreme outliers.\n', sum(rowsWithOutliers));

%% 2. Z-Score Continuous Predictors (CRITICAL for interactions)
% Now that outliers are gone, we can safely compute means and standard deviations
cleanTbl.RipplePower_Z = zscore(cleanTbl.RipplePower);
cleanTbl.SpindlePower_Match_Z = zscore(cleanTbl.SpindlePower_Match);
cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds);
cleanTbl.V1_logodds_PRE_Z = zscore(cleanTbl.V1_logodds_PRE);
cleanTbl.V1_logodds_Z = zscore(cleanTbl.V1_logodds); % Post-ripple

%% Calculate coherence
% The manuscript specifically uses z-scored log-odds for its LME models
cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds);
cleanTbl.V1_logodds_Z  = zscore(cleanTbl.V1_logodds);     % Post-ripple
cleanTbl.V1_logodds_PRE_Z = zscore(cleanTbl.V1_logodds_PRE); % Pre-ripple

% 2. Calculate Per-Ripple Coherence (Agreement Score)
% This score represents the "coupling strength" of each individual event.
% Positive = Coherent (Matching content)
% Negative = Incoherent (Mismatching content)
% Near Zero = One or both regions are neutral
cleanTbl.Event_Coherence_Post = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z;
cleanTbl.Event_Coherence_Pre  = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z;

% Verify the distribution
fprintf('Mean Event Coherence (Post): %.3f\n', mean(cleanTbl.Event_Coherence_Post));

%% Polar and Linear Plots: Coherence vs. SO Phase
% Choose which coherence metric to plot
targetCoherence = cleanTbl.Event_Coherence_Post; 

% 1. Define Phase Bins (12 bins of 30 degrees each)
numBins = 6;
edges = linspace(-pi, pi, numBins + 1);
centers = edges(1:end-1) + diff(edges)/2;

% To make a continuous line on a polar plot, we duplicate the first point
plotAngles = [centers, centers(1)]; 

% 2. Calculate Mean Coherence per Bin (MATCHED Hemisphere)
[~, ~, binIdx_M] = histcounts(cleanTbl.SOPhase_Match, edges);
meanCoh_M = zeros(1, numBins);
for b = 1:numBins
    meanCoh_M(b) = mean(targetCoherence(binIdx_M == b), 'omitnan');
end
plotCoh_M = [meanCoh_M, meanCoh_M(1)]; % Close the loop for polar

% 3. Calculate Mean Coherence per Bin (NON-MATCHED Hemisphere)
[~, ~, binIdx_NM] = histcounts(cleanTbl.SOPhase_NonMatch, edges);
meanCoh_NM = zeros(1, numBins);
for b = 1:numBins
    meanCoh_NM(b) = mean(targetCoherence(binIdx_NM == b), 'omitnan');
end
plotCoh_NM = [meanCoh_NM, meanCoh_NM(1)]; % Close the loop for polar

%% 4. Plot the 4-Panel Graph
figure('Color', 'w', 'Position', [100, 100, 900, 800]); % Made taller for 2 rows

% ==========================================
% TOP ROW: POLAR PLOTS
% ==========================================
% --- Plot 1: Matched Hemisphere (Polar) ---
ax1 = subplot(2, 2, 1, polaraxes);
polarplot(ax1, plotAngles, plotCoh_M, '-o', 'LineWidth', 2.5, 'Color', '#D95319');
ax1.ThetaZeroLocation = 'top';     
ax1.ThetaDir = 'clockwise';        
ax1.ThetaTick = [0 90 180 270];
ax1.ThetaTickLabel = {'0 (Peak)', '\pi/2', '\pm\pi (Trough)', '-\pi/2'};
title(ax1, {'Matched Hemisphere', '(Polar)'}, 'FontWeight', 'bold');

% --- Plot 2: Non-Matched Hemisphere (Polar) ---
ax2 = subplot(2, 2, 2, polaraxes);
polarplot(ax2, plotAngles, plotCoh_NM, '-o', 'LineWidth', 2.5, 'Color', '#0072BD');
ax2.ThetaZeroLocation = 'top';     
ax2.ThetaDir = 'clockwise';
ax2.ThetaTick = [0 90 180 270];
ax2.ThetaTickLabel = {'0 (Peak)', '\pi/2', '\pm\pi (Trough)', '-\pi/2'};
title(ax2, {'Non-Matched Hemisphere', '(Polar)'}, 'FontWeight', 'bold');

% Equalize radial limits for polar plots
maxRadius = max([plotCoh_M, plotCoh_NM]) * 1.1;
minRadius = min([plotCoh_M, plotCoh_NM, 0]); 
rlim(ax1, [minRadius, maxRadius]);
rlim(ax2, [minRadius, maxRadius]);

% ==========================================
% BOTTOM ROW: LINEAR PLOTS
% ==========================================
% --- Plot 3: Matched Hemisphere (Linear) ---
ax3 = subplot(2, 2, 3);
plot(centers, meanCoh_M, '-o', 'LineWidth', 2.5, 'Color', '#D95319', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 6);
grid on; box off;
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi (Trough)', '-\pi/2', '0 (Peak)', '\pi/2', '\pi (Trough)'});
xlabel('Matched SO Phase', 'FontWeight', 'bold');
ylabel('Mean Coherence', 'FontWeight', 'bold');
title('Matched Hemisphere (Linear)', 'FontSize', 12);

% --- Plot 4: Non-Matched Hemisphere (Linear) ---
ax4 = subplot(2, 2, 4);
plot(centers, meanCoh_NM, '-o', 'LineWidth', 2.5, 'Color', '#0072BD', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 6);
grid on; box off;
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi (Trough)', '-\pi/2', '0 (Peak)', '\pi/2', '\pi (Trough)'});
xlabel('Non-Matched SO Phase', 'FontWeight', 'bold');
ylabel('Mean Coherence', 'FontWeight', 'bold');
title('Non-Matched Hemisphere (Linear)', 'FontSize', 12);

% Equalize Y-axis limits for linear plots
linearMin = min([meanCoh_M, meanCoh_NM]) * 0.9;
linearMax = max([meanCoh_M, meanCoh_NM]) * 1.1;
ylim(ax3, [linearMin, linearMax]);
ylim(ax4, [linearMin, linearMax]);
xlim(ax3, [-pi, pi]);
xlim(ax4, [-pi, pi]);

sgtitle('Mean V1-HC Coherence across Slow Oscillation Phases', 'FontSize', 16, 'FontWeight', 'bold');


%% 4-Panel 2D Heatmap: Coherence, Occupancy, Ripple Power, Spindle Power
targetCoherence = cleanTbl.Event_Coherence_Post; % Focus on Post-Ripple
totalRipples = length(targetCoherence);

% 1. Shift the Phases mathematically from [-pi, pi] to [pi/2, 5pi/2]
shifted_Match = mod(cleanTbl.SOPhase_Match - pi/2, 2*pi) + pi/2;
shifted_NonMatch = mod(cleanTbl.SOPhase_NonMatch - pi/2, 2*pi) + pi/2;

% 2. Define Phase Bins over the new shifted range
numBins = 6; 
edges = linspace(pi/2, 5*pi/2, numBins + 1);
centers = edges(1:end-1) + diff(edges)/2;

% 3. Initialize the 2D Grids for all 4 metrics
coherenceMap   = NaN(numBins, numBins);
occupancyMap   = zeros(numBins, numBins); 
ripplePowerMap = NaN(numBins, numBins);
spindlePowerMap= NaN(numBins, numBins);

% 4. Discretize the phases into bin indices using the SHIFTED phases
[~, ~, matchIdx] = histcounts(shifted_Match, edges);
[~, ~, nonMatchIdx] = histcounts(shifted_NonMatch, edges);

% 5. Calculate Means & Occupancy for every combination
for r = 1:numBins % Row = Non-Matched
    for c = 1:numBins % Col = Matched
        validRipples = (nonMatchIdx == r) & (matchIdx == c);
        
        % Calculate Occupancy (Proportion of total ripples)
        occupancyMap(r, c) = sum(validRipples) / totalRipples;
        
        % Calculate Means for the other metrics (require at least 5 ripples)
        if sum(validRipples) > 5
            coherenceMap(r, c)    = mean(targetCoherence(validRipples), 'omitnan');
            ripplePowerMap(r, c)  = mean(cleanTbl.RipplePower_Z(validRipples), 'omitnan');
            spindlePowerMap(r, c) = mean(cleanTbl.SpindlePower_Match_Z(validRipples), 'omitnan');
        end
    end
end

% 6. Plot the 4 Panels
figure('Color', 'w', 'Position', [100, 100, 1100, 950]); 

% ==========================================
% PANEL 1: COHERENCE (Top Left)
% ==========================================
ax1 = subplot(2, 2, 1);
imagesc(centers, centers, coherenceMap);
axis xy; colormap(ax1, viridis);
c1 = colorbar; c1.Label.String = 'Mean Coherence (Z)'; c1.Label.FontWeight = 'bold';
title('1. Coherence: State Space Competition', 'FontSize', 13);

% ==========================================
% PANEL 2: OCCUPANCY (Top Right)
% ==========================================
ax2 = subplot(2, 2, 2);
imagesc(centers, centers, occupancyMap);
axis xy; colormap(ax2, flipud(gray)); 
c2 = colorbar; c2.Label.String = 'Proportion of Events'; c2.Label.FontWeight = 'bold';
title('2. Occupancy: Distribution of Ripples', 'FontSize', 13);

% ==========================================
% PANEL 3: RIPPLE POWER (Bottom Left)
% ==========================================
ax3 = subplot(2, 2, 3);
imagesc(centers, centers, ripplePowerMap);
axis xy; colormap(ax3, plasma); % Using plasma to distinguish from coherence
c3 = colorbar; c3.Label.String = 'Mean Ripple Power (Z)'; c3.Label.FontWeight = 'bold';
title('3. Mean Ripple Power', 'FontSize', 13);

% ==========================================
% PANEL 4: SPINDLE POWER (Bottom Right)
% ==========================================
ax4 = subplot(2, 2, 4);
imagesc(centers, centers, spindlePowerMap);
axis xy; colormap(ax4, plasma); 
c4 = colorbar; c4.Label.String = 'Mean Matched Spindle Power (Z)'; c4.Label.FontWeight = 'bold';
title('4. Mean Matched Spindle Power', 'FontSize', 13);

sgtitle('Cortico-Hippocampal Circuit: State Space Dynamics', 'FontSize', 16, 'FontWeight', 'bold');

% 7. Apply Formatting to All Axes (The Bug Fix)
axesList = [ax1, ax2, ax3, ax4];
for i = 1:length(axesList)
    ax = axesList(i);
    
    % Set Ticks and Labels
    xticks(ax, [pi/2, pi, 3*pi/2, 2*pi, 5*pi/2]);
    xticklabels(ax, {'\pi/2', '\pi (Trough)', '-\pi/2', '0 (Peak)', '\pi/2'});
    yticks(ax, [pi/2, pi, 3*pi/2, 2*pi, 5*pi/2]);
    yticklabels(ax, {'\pi/2', '\pi (Trough)', '-\pi/2', '0 (Peak)', '\pi/2'});
    
    % Set Axis Labels
    xlabel(ax, 'Matched Hemisphere SO Phase', 'FontWeight', 'bold', 'FontSize', 10);
    ylabel(ax, 'Non-Matched Hemisphere SO Phase', 'FontWeight', 'bold', 'FontSize', 10);
    
    % Add Reference Lines depending on the plot theme
    if ax == ax2 % If it's the grayscale occupancy plot, use dark lines
        lineColor = 'k';
    else % For all the colored heatmaps, use white lines
        lineColor = 'w';
    end
    
    hold(ax, 'on');
    xline(ax, pi, [lineColor ':'], 'LineWidth', 1.5, 'Alpha', 0.6);   % Vertical Trough
    yline(ax, pi, [lineColor ':'], 'LineWidth', 1.5, 'Alpha', 0.6);   % Horizontal Trough
    xline(ax, 2*pi, [lineColor '--'], 'LineWidth', 1.5, 'Alpha', 0.6); % Vertical Peak
    yline(ax, 2*pi, [lineColor '--'], 'LineWidth', 1.5, 'Alpha', 0.6); % Horizontal Peak
end



%% 6-Panel Analysis: Ripple, Match Spindle, and Non-Match Spindle Power
targetCoherence = cleanTbl.Event_Coherence_Post; 
totalRipples = length(targetCoherence);

cleanTbl.SpindlePower_NonMatch_Z = (cleanTbl.SpindlePower_NonMatch - mean(cleanTbl.SpindlePower_NonMatch, 'omitnan')) ./ std(cleanTbl.SpindlePower_NonMatch, 'omitnan');

% 1. Create Bins based on Percentiles (Sextiles for numBins=6)
numBins = 6;
p_edges = linspace(0, 100, numBins + 1);
ripEdges = prctile(cleanTbl.RipplePower_Z, p_edges);
spiMEdges = prctile(cleanTbl.SpindlePower_Match_Z, p_edges);
spiNMEdges = prctile(cleanTbl.SpindlePower_NonMatch_Z, p_edges);

[~, ~, ripIdx]  = histcounts(cleanTbl.RipplePower_Z, ripEdges);
[~, ~, spiMIdx] = histcounts(cleanTbl.SpindlePower_Match_Z, spiMEdges);
[~, ~, spiNMIdx] = histcounts(cleanTbl.SpindlePower_NonMatch_Z, spiNMEdges);

% 2. Calculate 1D Statistics for the 3 Line Plots
% (Pre-allocate arrays)
[rC, rCoh, rErr, sMC, sMCoh, sMErr, sNMC, sNMCoh, sNMErr] = deal(zeros(1, numBins));

for i = 1:numBins
    % Ripple
    mR = (ripIdx == i);
    rC(i) = mean(cleanTbl.RipplePower_Z(mR), 'omitnan');
    rCoh(i) = mean(targetCoherence(mR), 'omitnan');
    rErr(i) = std(targetCoherence(mR), 'omitnan') / sqrt(sum(mR));
    % Match Spindle
    mSM = (spiMIdx == i);
    sMC(i) = mean(cleanTbl.SpindlePower_Match_Z(mSM), 'omitnan');
    sMCoh(i) = mean(targetCoherence(mSM), 'omitnan');
    sMErr(i) = std(targetCoherence(mSM), 'omitnan') / sqrt(sum(mSM));
    % Non-Match Spindle
    mSNM = (spiNMIdx == i);
    sNMC(i) = mean(cleanTbl.SpindlePower_NonMatch_Z(mSNM), 'omitnan');
    sNMCoh(i) = mean(targetCoherence(mSNM), 'omitnan');
    sNMErr(i) = std(targetCoherence(mSNM), 'omitnan') / sqrt(sum(mSNM));
end

% 3. Calculate 2D Statistics for Heatmaps
cohMapM  = NaN(numBins, numBins);
cohMapNM = NaN(numBins, numBins);
occMap   = zeros(numBins, numBins);

for r = 1:numBins
    for c = 1:numBins
        % Match Map & Occupancy (Using Match for occupancy ref)
        maskM = (spiMIdx == r) & (ripIdx == c);
        occMap(r, c) = sum(maskM) / totalRipples;
        if sum(maskM) > 10, cohMapM(r, c) = mean(targetCoherence(maskM), 'omitnan'); end
        
        % Non-Match Map
        maskNM = (spiNMIdx == r) & (ripIdx == c);
        if sum(maskNM) > 10, cohMapNM(r, c) = mean(targetCoherence(maskNM), 'omitnan'); end
    end
end

% 4. Plotting (2 Rows x 3 Columns)
figure('Color', 'w', 'Position', [50, 50, 1400, 850]);
tL = {'1st', '2nd', '3rd', '4th', '5th', '6th'}; % Labels for 6 bins

% --- Top Row: 1D Line Plots ---
axLine = zeros(1,3);
axLine(1) = subplot(2, 3, 1); errorbar(rC, rCoh, rErr, '-o', 'Color', '#D95319', 'LineWidth', 2); title('1. Ripple Power');ylim([0 0.12]);
axLine(2) = subplot(2, 3, 2); errorbar(sMC, sMCoh, sMErr, '-o', 'Color', '#0072BD', 'LineWidth', 2); title('2. Match Spindle Power');ylim([0 0.12]);
axLine(3) = subplot(2, 3, 3); errorbar(sNMC, sNMCoh, sNMErr, '-o', 'Color', '#7E2F8E', 'LineWidth', 2); title('3. Non-Match Spindle Power');ylim([0 0.12]);

% Format Line Plots
allY = [rCoh, sMCoh, sNMCoh];
for i=1:3
    grid(axLine(i), 'on'); box(axLine(i), 'off'); 
    ylabel(axLine(i), 'Coherence (Z)'); xlabel(axLine(i), 'Power (Z)');
end

% --- Bottom Row: 2D Heatmaps ---
% Match Heatmap
ax4 = subplot(2, 3, 4); imagesc(cohMapM); axis xy; colormap(ax4, viridis); colorbar;
title('4. Coherence: Rip x Match Spi'); setupHeatmap(ax4, tL, 'Match Spindle');

% Non-Match Heatmap
ax5 = subplot(2, 3, 5); imagesc(cohMapNM); axis xy; colormap(ax5, viridis); colorbar;
title('5. Coherence: Rip x Non-Match Spi'); setupHeatmap(ax5, tL, 'Non-Match Spindle');

% Occupancy (Match side example)
ax6 = subplot(2, 3, 6); imagesc(occMap); axis xy; colormap(ax6, flipud(gray)); colorbar;
title('6. Occupancy (Rip x Match Spi)'); setupHeatmap(ax6, tL, 'Match Spindle');

% Sync color scales for Coherence heatmaps
cLimits = [min([cohMapM(:); cohMapNM(:)], [], 'omitnan'), max([cohMapM(:); cohMapNM(:)], [], 'omitnan')];
clim(ax4, cLimits); clim(ax5, cLimits);

sgtitle('Post-Ripple Coherence: Match vs Non-Match Dynamics', 'FontSize', 18, 'FontWeight', 'bold');

%% Heatmap: Match vs. Non-Match Spindle Power Interaction
targetCoherence = cleanTbl.Event_Coherence_Post; 
totalRipples = length(targetCoherence);

% 1. Define Bins based on Percentiles (Sextiles)
numBins = 6;
p_edges = linspace(0, 100, numBins + 1);
spiMEdges  = prctile(cleanTbl.SpindlePower_Match_Z, p_edges);
spiNMEdges = prctile(cleanTbl.SpindlePower_NonMatch_Z, p_edges);

% Discretize into bin indices
[~, ~, spiMIdx]  = histcounts(cleanTbl.SpindlePower_Match_Z, spiMEdges);
[~, ~, spiNMIdx] = histcounts(cleanTbl.SpindlePower_NonMatch_Z, spiNMEdges);

% 2. Initialize Grids
cohMap_Inter = NaN(numBins, numBins);
occMap_Inter = zeros(numBins, numBins);

% 3. Calculate Mean Coherence & Occupancy
for r = 1:numBins % Y-axis: Non-Match Spindle Power
    for c = 1:numBins % X-axis: Match Spindle Power
        mask = (spiNMIdx == r) & (spiMIdx == c);
        
        occMap_Inter(r, c) = sum(mask) / totalRipples;
        
        if sum(mask) > 10
            cohMap_Inter(r, c) = mean(targetCoherence(mask), 'omitnan');
        end
    end
end

% 4. Plotting the Comparison
figure('Color', 'w', 'Position', [100, 100, 1100, 500]);
tL = {'1st', '2nd', '3rd', '4th', '5th', '6th'};

% --- PANEL 1: Coherence Interaction ---
ax1 = subplot(1, 2, 1);
imagesc(cohMap_Inter);
axis xy; colormap(ax1, viridis);
c1 = colorbar; c1.Label.String = 'Mean Coherence (Z)';
c1.Label.FontWeight = 'bold';

set(ax1, 'XTick', 1:numBins, 'XTickLabel', tL, 'YTick', 1:numBins, 'YTickLabel', tL);
xlabel('Matched Spindle Power (Sextile)', 'FontWeight', 'bold');
ylabel('Non-Matched Spindle Power (Sextile)', 'FontWeight', 'bold');
title('Coherence: Inter-Hemispheric Interaction', 'FontSize', 14);

% --- PANEL 2: Occupancy (State Density) ---
ax2 = subplot(1, 2, 2);
imagesc(occMap_Inter);
axis xy; colormap(ax2, flipud(gray));
c2 = colorbar; c2.Label.String = 'Proportion of Events';
c2.Label.FontWeight = 'bold';

set(ax2, 'XTick', 1:numBins, 'XTickLabel', tL, 'YTick', 1:numBins, 'YTickLabel', tL);
xlabel('Matched Spindle Power (Sextile)', 'FontWeight', 'bold');
ylabel('Non-Matched Spindle Power (Sextile)', 'FontWeight', 'bold');
title('Occupancy: Bilateral Spindle Coordination', 'FontSize', 14);

sgtitle('Competition Analysis: Do Hemispheres Interfere?', 'FontSize', 16, 'FontWeight', 'bold');


% Helper function for axes
function setupHeatmap(ax, labels, yLab)
    set(ax, 'XTick', 1:length(labels), 'XTickLabel', labels, 'YTick', 1:length(labels), 'YTickLabel', labels);
    xlabel(ax, 'Ripple Power (Sextile)'); ylabel(ax, [yLab, ' (Sextile)']);
end