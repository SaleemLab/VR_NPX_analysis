function output = predict_ripples_by_SO_spindles(ripple_info, varargin)
% Analyse how SO and spindles phase and power at ripple peak modulates
% ripple power.


% --- Input Parser
p = inputParser;
addParameter(p, 'nPerm', 1000);
addParameter(p, 'nBoot', 1000);
addParameter(p, 'output', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
parse(p, varargin{:});

nBoot = p.Results.nBoot;
subject_id = p.Results.subject_id;
ripples_index = find(~isnan(ripple_info.ripple_power));  % All valid events
output = p.Results.output;

% --- Preallocate output
if isempty(output)
 parfor iBoot = 1:nBoot
     tic
     s = RandStream('philox4x32_10', 'Seed', iBoot);
     index = randsample(s, ripples_index, numel(ripples_index), true);

     % Extract base variables
     ripple_power = normalize(ripple_info.ripple_power(index));
     subjectID = categorical(subject_id(index));

     % Spindle
     sp_phase_ipsi   = ripple_info.spindle_phase(index,1);
     sp_phase_contra = ripple_info.spindle_phase(index,2);
     spindle_amp_ipsi   = normalize(ripple_info.spindle_amplitude(index,1));
     spindle_amp_contra = normalize(ripple_info.spindle_amplitude(index,2));
     cos_sp_ipsi = cos(sp_phase_ipsi); sin_sp_ipsi = sin(sp_phase_ipsi);
     cos_sp_contra = cos(sp_phase_contra); sin_sp_contra = sin(sp_phase_contra);

     % SO
     SO_phase_ipsi   = ripple_info.SO_phase(index,1);
     SO_phase_contra = ripple_info.SO_phase(index,2);
     SO_amp_ipsi   = normalize(ripple_info.SO_amplitude(index,1));
     SO_amp_contra = normalize(ripple_info.SO_amplitude(index,2));
     cos_SO_ipsi = cos(SO_phase_ipsi); sin_SO_ipsi = sin(SO_phase_ipsi);
     cos_SO_contra = cos(SO_phase_contra); sin_SO_contra = sin(SO_phase_contra);

     % Full table
     tbl_all = table(ripple_power, subjectID, ...
         cos_sp_ipsi, sin_sp_ipsi, cos_sp_contra, sin_sp_contra, ...
         spindle_amp_ipsi, spindle_amp_contra, ...
         cos_SO_ipsi, sin_SO_ipsi, cos_SO_contra, sin_SO_contra, ...
         SO_amp_ipsi, SO_amp_contra);

     % Model formulas
     modelList = {
         'ripple_power ~ cos_sp_ipsi + sin_sp_ipsi + (1|subjectID)';
         'ripple_power ~ cos_sp_contra + sin_sp_contra + (1|subjectID)';
         'ripple_power ~ cos_sp_ipsi + sin_sp_ipsi + cos_sp_contra + sin_sp_contra + (1|subjectID)';

         'ripple_power ~ spindle_amp_ipsi + (1|subjectID)';
         'ripple_power ~ spindle_amp_contra + (1|subjectID)';
         'ripple_power ~ spindle_amp_ipsi + spindle_amp_contra + (1|subjectID)';

         'ripple_power ~ cos_SO_ipsi + sin_SO_ipsi + (1|subjectID)';
         'ripple_power ~ cos_SO_contra + sin_SO_contra + (1|subjectID)';
         'ripple_power ~ cos_SO_ipsi + sin_SO_ipsi + cos_SO_contra + sin_SO_contra + (1|subjectID)';

         'ripple_power ~ SO_amp_ipsi + (1|subjectID)';
         'ripple_power ~ SO_amp_contra + (1|subjectID)';
         'ripple_power ~ SO_amp_ipsi + SO_amp_contra + (1|subjectID)';

         'ripple_power ~ cos_sp_ipsi + sin_sp_ipsi + spindle_amp_ipsi + (1|subjectID)';
         'ripple_power ~ cos_SO_ipsi + sin_SO_ipsi + SO_amp_ipsi + (1|subjectID)';
         };

     modelNames = {
         'spindle_phase_ipsi';
         'spindle_phase_contra';
         'spindle_phase_both';

         'spindle_amp_ipsi';
         'spindle_amp_contra';
         'spindle_amp_both';

         'SO_phase_ipsi';
         'SO_phase_contra';
         'SO_phase_both';

         'SO_amp_ipsi';
         'SO_amp_contra';
         'SO_amp_both';

         'spindle_phase_amp_ipsi';
         'SO_phase_amp_ipsi';
         };


     % Preallocate local results
     nModels = numel(modelList);
     local_b = cell(nModels,1);
     local_tstat = cell(nModels,1);
     local_pval = cell(nModels,1);
     local_variable = cell(nModels,1);
     local_R2 = zeros(nModels,1);

     for m = 1:nModels
         thisModel = modelList{m};

         % Extract variables from model
         tokens = regexp(thisModel, '[a-zA-Z_][a-zA-Z_0-9]*', 'match');
         modelVars = unique(tokens(2:end));  % skip 'ripple_power'

         % Subset table for this model
         tbl = tbl_all(:, modelVars);
         tbl.ripple_power = ripple_power;

         % Fit mixed-effects model
         glme = fitlme(tbl, thisModel);

         local_b{m} = glme.Coefficients.Estimate(2:end);
         local_tstat{m} = glme.Coefficients.tStat(2:end);
         local_pval{m} = glme.Coefficients.pValue(2:end);
         local_variable{m} = glme.CoefficientNames(2:end);

         if isprop(glme, 'Rsquared') && isfield(glme.Rsquared, 'Adjusted')
             local_R2(m) = glme.Rsquared.Adjusted;
         else
             local_R2(m) = NaN;  % fallback
         end
     end

     % Store in output
     output(iBoot).b = local_b;
     output(iBoot).R2 = local_R2;
     output(iBoot).t_stat = local_tstat;
     output(iBoot).pval = local_pval;
     output(iBoot).variable = local_variable;
     output(iBoot).model = modelList;
     output(iBoot).type = modelNames;
     toc
 end

end

% --- Optional plotting
if p.Results.plot_option

    customColors = [0,90,50;74,20,134; 228,42,168 ]/256; % dark purple, dark green, magenta

    for i = 1:1000
        b1 = output(i).b{1}; % ipsi model: [intercept, cos, sin]
        b2 = output(i).b{2}; % contra model

        beta_cos_ipsi = b1(1);
        beta_sin_ipsi = b1(2);
        beta_cos_contra = b2(1);
        beta_sin_contra = b2(2);

        mod_depth_ipsi(i) = sqrt(beta_cos_ipsi^2 + beta_sin_ipsi^2);
        mod_depth_contra(i) = sqrt(beta_cos_contra^2 + beta_sin_contra^2);

        pref_phase_ipsi(i) = atan2(beta_sin_ipsi, beta_cos_ipsi);
        pref_phase_contra(i) = atan2(beta_sin_contra, beta_cos_contra);

        pval_cos_ipsi(i) = output(i).pval{1}(1);  % cos
        pval_sin_ipsi(i) = output(i).pval{1}(2);  % sin

        pval_cos_contra(i) = output(i).pval{2}(1);
        pval_sin_contra(i) = output(i).pval{2}(2);

        R2_ipsi(i) = output(i).R2(1);
        R2_contra(i) = output(i).R2(2);
    end

    % === Preferred phase and R² ===
    medianR2_ipsi = median(R2_ipsi);
    medianR2_contra = median(R2_contra);

    medianP_cos_ipsi = median(pval_cos_ipsi);
    medianP_sin_ipsi = median(pval_sin_ipsi);
    medianP_cos_contra = median(pval_cos_contra);
    medianP_sin_contra = median(pval_sin_contra);

    pref_ipsi = median(pref_phase_ipsi);
    pref_contra = median(pref_phase_contra);

    % === Mod depth stats ===
    mean_mod = [mean(mod_depth_ipsi), mean(mod_depth_contra)];
    ci_mod = prctile([mod_depth_ipsi; mod_depth_contra]', [2.5 97.5]);
    lowerCI = mean_mod - ci_mod(1,:);
    upperCI = ci_mod(2,:) - mean_mod;

    % === Colors ===
    customColors = [0,90,50;74,20,134]/256;

    % === FIGURE ===
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name = 'Ripple power predicted by spindle phase';
    tiledlayout(2,2, 'TileSpacing', 'tight');

    % === SCATTER 1: IPSI ===
    nexttile
    scatter(ripple_info.spindle_phase(:,1), ripple_info.ripple_power, 10, ...
        'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.03);
    hold on
    yline(mean(ripple_info.ripple_power), 'k--');
    xline(pref_ipsi, 'r-', 'LineWidth', 2);
    text(pi/2, prctile(ripple_info.ripple_power,99.9), ...
        sprintf('\\phi = %.2f rad\nR^2 = %.3f\np_{cos} = %.3e\np_{sin} = %.3e', ...
        pref_ipsi, medianR2_ipsi, medianP_cos_ipsi, medianP_sin_ipsi), ...
        'Color','r','FontSize',10, 'HorizontalAlignment','center');
    xlabel('Spindle Phase (ipsi)'); ylabel('Ripple Power');
    title('Ipsi Spindle Phase Modulation');
    set(gca,'TickDir','out','Box','off','FontSize',12)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    ylim([4 30])
    xlim([-pi pi])


    % === SCATTER 2: CONTRA ===
    nexttile
    scatter(ripple_info.spindle_phase(:,2), ripple_info.ripple_power, 10, ...
        'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.03);
    hold on
    yline(mean(ripple_info.ripple_power), 'k--');
    xline(pref_contra, 'r-', 'LineWidth', 2);
    text(pi/2, prctile(ripple_info.ripple_power,99.9), ...
        sprintf('\\phi = %.2f rad\nR^2 = %.3f\np_{cos} = %.2e\np_{sin} = %.2e', ...
        pref_contra, medianR2_contra, medianP_cos_contra, medianP_sin_contra), ...
        'Color','r','FontSize',10, 'HorizontalAlignment','center');
    xlabel('Spindle Phase (contra)'); ylabel('Ripple Power');
    title('Contra Spindle Phase Modulation');
    set(gca,'TickDir','out','Box','off','FontSize',12)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    ylim([4 30])
    xlim([-pi pi])

    % === BAR PLOT: Mod depth ===
    nexttile
    barColors = customColors;
    x_pos = [1 2];
    bar(x_pos, mean_mod, 0.5, 'FaceColor', 'flat', 'CData', barColors, 'EdgeColor', 'none','FaceAlpha',0.5); hold on;
    errorbar(x_pos, mean_mod, lowerCI, upperCI, ...
        'k', 'LineStyle','none', 'LineWidth',1.5);
    xticks(x_pos); xticklabels({'Ipsi', 'Contra'});
    ylabel('Modulation Depth');
    title('Modulation Depth (95% CI)');
    set(gca,'TickDir','out','Box','off','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spindple amplitude

for nBoot = 1:1000
    b_spindle_ipsi(nBoot) = output(nBoot).b{4};
    b_spindle_contra(nBoot) = output(nBoot).b{5};
    t_spindle_ipsi(nBoot) = output(nBoot).t_stat{4};
    t_spindle_contra(nBoot) = output(nBoot).t_stat{5};
    R2_spindle_ipsi(nBoot) = output(nBoot).R2(4);
    R2_spindle_contra(nBoot) = output(nBoot).R2(5);
    pval_spindle_ipsi(nBoot) = output(nBoot).pval{4};
    pval_spindle_contra(nBoot) = output(nBoot).pval{5};
end



fig = figure('Color','w');
fig.Position = [350 59 800 930/3*2];
fig.Name = 'Ripple power predicted by spindle power';
tiledlayout(2,2, 'TileSpacing', 'tight');
% --- IPSI SCATTER
nexttile
scatter(ripple_info.spindle_amplitude(:,1), ripple_info.ripple_power, 10, ...
    'filled','MarkerFaceColor',customColors(1,:), 'MarkerFaceAlpha',0.03);
hold on
coeffs = polyfit(ripple_info.spindle_amplitude(:,1), ripple_info.ripple_power, 1);
x_fit = linspace(min(ripple_info.spindle_amplitude(:,1)), max(ripple_info.spindle_amplitude(:,1)), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
text(prctile(ripple_info.spindle_amplitude(1,:),99.9), prctile(ripple_info.ripple_power,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', median(R2_spindle_ipsi), median(pval_spindle_ipsi)), ...
    'FontSize',10,'Color','r');
xlabel('Spindle power')
ylabel('Ripple power')
title('Ipsi spindle power')
set(gca,'TickDir','out','Box','off','FontSize',12)
ylim([4 30])
xlim([-2 6])

% --- CONTRA SCATTER
nexttile
scatter(ripple_info.spindle_amplitude(:,2), ripple_info.ripple_power, 10, ...
    'filled','MarkerFaceColor',customColors(2,:), 'MarkerFaceAlpha',0.03);
hold on
coeffs = polyfit(ripple_info.spindle_amplitude(:,2), ripple_info.ripple_power, 1);
x_fit = linspace(min(ripple_info.spindle_amplitude(:,2)), max(ripple_info.spindle_amplitude(:,2)), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
text(prctile(ripple_info.spindle_amplitude(2,:),99.9), prctile(ripple_info.ripple_power,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', median(R2_spindle_contra), median(pval_spindle_contra)), ...
    'FontSize',10,'Color','r');
xlabel('Spindle power')
ylabel('Ripple power')
title('Contra spindle power')
set(gca,'TickDir','out','Box','off','FontSize',12)
ylim([4 30])
xlim([-2 6])

% --- BAR: T-STATISTICS
nexttile
barData = [mean(t_spindle_ipsi), mean(t_spindle_contra)];
ci_ipsi = prctile(t_spindle_ipsi, [2.5 97.5]);
ci_contra = prctile(t_spindle_contra, [2.5 97.5]);
lowerCI = [barData(1)-ci_ipsi(1), barData(2)-ci_contra(1)];
upperCI = [ci_ipsi(2)-barData(1), ci_contra(2)-barData(2)];
barColors = customColors;
x_pos = [1 2];
for i = 1:2
    bar(x_pos(i), barData(i), 0.4, 'FaceColor', barColors(i,:), 'EdgeColor','none','FaceAlpha',0.5); hold on
    errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), 'k', 'LineWidth', 1.5)
end
xticks(x_pos)
xticklabels({'ipsi','contra'})
ylabel('T-statistic')
title('Spindle power effect')
set(gca,'TickDir','out','Box','off','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SO phase

customColors = [0,90,50;74,20,134; 228,42,168 ]/256; % dark purple, dark green, magenta

% === Extract SO model info (model 7 & 8)
for i = 1:1000
    b1 = output(i).b{7}; % model 7: ipsi
    b2 = output(i).b{8}; % model 8: contra

    beta_cos_ipsi = b1(1);
    beta_sin_ipsi = b1(2);
    beta_cos_contra = b2(1);
    beta_sin_contra = b2(2);

    mod_depth_ipsi_SO(i) = sqrt(beta_cos_ipsi^2 + beta_sin_ipsi^2);
    mod_depth_contra_SO(i) = sqrt(beta_cos_contra^2 + beta_sin_contra^2);

    pref_phase_ipsi_SO(i) = atan2(beta_sin_ipsi, beta_cos_ipsi);
    pref_phase_contra_SO(i) = atan2(beta_sin_contra, beta_cos_contra);

    pval_cos_ipsi_SO(i) = output(i).pval{7}(1);  % cos
    pval_sin_ipsi_SO(i) = output(i).pval{7}(2);  % sin
    pval_cos_contra_SO(i) = output(i).pval{8}(1);
    pval_sin_contra_SO(i) = output(i).pval{8}(2);

    R2_ipsi_SO(i) = output(i).R2(7);
    R2_contra_SO(i) = output(i).R2(8);
end

% === Summary stats
medianR2_ipsi = median(R2_ipsi_SO);
medianR2_contra = median(R2_contra_SO);
medianP_cos_ipsi = median(pval_cos_ipsi_SO);
medianP_sin_ipsi = median(pval_sin_ipsi_SO);
medianP_cos_contra = median(pval_cos_contra_SO);
medianP_sin_contra = median(pval_sin_contra_SO);
pref_ipsi = median(pref_phase_ipsi_SO);
pref_contra = median(pref_phase_contra_SO);
mean_mod = [mean(mod_depth_ipsi_SO), mean(mod_depth_contra_SO)];
ci_mod = prctile([mod_depth_ipsi_SO; mod_depth_contra_SO]', [2.5 97.5]);
lowerCI = mean_mod - ci_mod(1,:);
upperCI = ci_mod(2,:) - mean_mod;

% === Colors
customColors = [0,90,50;74,20,134]/256;

% === FIGURE ===
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name = 'Ripple power predicted by SO phase';
    tiledlayout(2,2, 'TileSpacing', 'tight');
tiledlayout(2,2, 'TileSpacing', 'tight');

% === SCATTER: IPSI SO Phase
nexttile
scatter(ripple_info.SO_phase(:,1), ripple_info.ripple_power, 10, ...
    'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.02);
hold on
yline(mean(ripple_info.ripple_power), 'k--');
xline(pref_ipsi, 'r-', 'LineWidth', 2);
text(pi/2, prctile(ripple_info.ripple_power,99.9), ...
    sprintf('\\phi = %.2f rad\nR^2 = %.3f\np_{cos} = %.2e\np_{sin} = %.2e', ...
    pref_ipsi, medianR2_ipsi, medianP_cos_ipsi, medianP_sin_ipsi), ...
    'Color','r','FontSize',10, 'HorizontalAlignment','center');
xlabel('SO Phase (ipsi)'); ylabel('Ripple Power');
title('Ipsi SO Phase Modulation');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlim([-pi pi]); ylim([4 30])
set(gca,'TickDir','out','Box','off','FontSize',12)

% === SCATTER: CONTRA SO Phase
nexttile
scatter(ripple_info.SO_phase(:,2), ripple_info.ripple_power, 10, ...
    'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.02);
hold on
yline(mean(ripple_info.ripple_power), 'k--');
xline(pref_contra, 'r-', 'LineWidth', 2);
text(pi/2, prctile(ripple_info.ripple_power,99.9), ...
    sprintf('\\phi = %.2f rad\nR^2 = %.3f\np_{cos} = %.2e\np_{sin} = %.2e', ...
    pref_contra, medianR2_contra, medianP_cos_contra, medianP_sin_contra), ...
    'Color','r','FontSize',10, 'HorizontalAlignment','center');
xlabel('SO Phase (contra)'); ylabel('Ripple Power');
title('Contra SO Phase Modulation');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlim([-pi pi]); ylim([4 30])
set(gca,'TickDir','out','Box','off','FontSize',12)

% === BAR PLOT: Modulation Depth
nexttile
x_pos = [1 2];
bar(x_pos, mean_mod, 0.5, 'FaceColor', 'flat', 'CData', customColors, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
errorbar(x_pos, mean_mod, lowerCI, upperCI, ...
    'k', 'LineStyle','none', 'LineWidth',1.5);
xticks(x_pos); xticklabels({'Ipsi', 'Contra'});
ylabel('Modulation Depth');
title('SO Phase Modulation Depth (95% CI)');
set(gca,'TickDir','out','Box','off','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So amplitude

% === Extract values for SO amplitude (model 10 = ipsi, model 11 = contra)
for nBoot = 1:1000
    b_SO_ipsi(nBoot) = output(nBoot).b{10};
    b_SO_contra(nBoot) = output(nBoot).b{11};
    t_SO_ipsi(nBoot) = output(nBoot).t_stat{10};
    t_SO_contra(nBoot) = output(nBoot).t_stat{11};
    R2_SO_ipsi(nBoot) = output(nBoot).R2(10);
    R2_SO_contra(nBoot) = output(nBoot).R2(11);
    pval_SO_ipsi(nBoot) = output(nBoot).pval{10};
    pval_SO_contra(nBoot) = output(nBoot).pval{11};
end

fig = figure('Color','w');
fig.Position = [350 59 800 930/3*2];
fig.Name = 'Ripple power predicted by SO power';
tiledlayout(2,2, 'TileSpacing', 'tight');

% --- IPSI SO SCATTER
nexttile
scatter(ripple_info.SO_amplitude(:,1), ripple_info.ripple_power, 10, ...
    'filled','MarkerFaceColor',customColors(1,:), 'MarkerFaceAlpha',0.03);
hold on
coeffs = polyfit(ripple_info.SO_amplitude(:,1), ripple_info.ripple_power, 1);
x_fit = linspace(min(ripple_info.SO_amplitude(:,1)), max(ripple_info.SO_amplitude(:,1)), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
text(prctile(ripple_info.SO_amplitude(:,1),99), prctile(ripple_info.ripple_power,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', median(R2_SO_ipsi), median(pval_SO_ipsi)), ...
    'FontSize',10,'Color','r');
xlabel('SO power')
ylabel('Ripple power')
title('Ipsi SO power')
set(gca,'TickDir','out','Box','off','FontSize',12)
ylim([4 30])
xlim([-2 6])

% --- CONTRA SO SCATTER
nexttile
scatter(ripple_info.SO_amplitude(:,2), ripple_info.ripple_power, 10, ...
    'filled','MarkerFaceColor',customColors(2,:), 'MarkerFaceAlpha',0.03);
hold on
coeffs = polyfit(ripple_info.SO_amplitude(:,2), ripple_info.ripple_power, 1);
x_fit = linspace(min(ripple_info.SO_amplitude(:,2)), max(ripple_info.SO_amplitude(:,2)), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
text(prctile(ripple_info.SO_amplitude(:,2),99), prctile(ripple_info.ripple_power,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', median(R2_SO_contra), median(pval_SO_contra)), ...
    'FontSize',10,'Color','r');
xlabel('SO power')
ylabel('Ripple power')
title('Contra SO power')
set(gca,'TickDir','out','Box','off','FontSize',12)
ylim([4 30])
xlim([-2 6])

% --- BAR: T-STATISTICS
nexttile
barData = [mean(t_SO_ipsi), mean(t_SO_contra)];
ci_ipsi = prctile(t_SO_ipsi, [2.5 97.5]);
ci_contra = prctile(t_SO_contra, [2.5 97.5]);
lowerCI = [barData(1)-ci_ipsi(1), barData(2)-ci_contra(1)];
upperCI = [ci_ipsi(2)-barData(1), ci_contra(2)-barData(2)];
barColors = customColors;
x_pos = [1 2];
for i = 1:2
    bar(x_pos(i), barData(i), 0.4, 'FaceColor', barColors(i,:), 'EdgeColor','none','FaceAlpha',0.5); hold on
    errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), 'k', 'LineWidth', 1.5)
end
xticks(x_pos)
xticklabels({'ipsi','contra'})
ylabel('T-statistic')
title('SO power effect')
set(gca,'TickDir','out','Box','off','FontSize',12)


end
