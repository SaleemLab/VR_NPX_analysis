
%% Plotting context selecitve ripple modulation
% scatter(context_modulation_all.z_FR_track(1,V1_id) - context_modulation_all.z_FR_track(2,V1_id),context_modulation_all.PRE_ripple_FR(V1_id))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','z_V1_population_ripple_PSTH.mat'),'z_V1_population_ripple_PSTH')
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_modulation_all.mat'),'context_modulation_all')
load(fullfile(analysis_folder,'ripple_modulation_PSTH_all_POST.mat'),'ripple_modulation_PSTH_all')


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
    0.85, 0.2, 0.2;  % red
    0.2, 0.4, 0.8    % blue
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

clear output
tbl_all = {tbl_V1,tbl_HPC};
for nregion = 1:2
    parfor iBoot = 1:1000
        local_b = cell(length(ModelList),1);
        local_tstat = cell(length(ModelList),1);
        local_pval = cell(length(ModelList),1);
        local_variable = cell(length(ModelList),1);
        local_R2 = zeros(length(ModelList),1);

        tic


        for m = 1:length(ModelList)
            tbl = tbl_all{nregion};
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            index = randsample(s, 1:height(tbl), height(tbl), true);
            tbl = tbl(index,:);

            glme = fitlme(tbl, ModelList{m});
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

        output(iBoot).b = local_b;
        output(iBoot).R2 = local_R2;
        output(iBoot).t_stat = local_tstat;
        output(iBoot).pval = local_pval;
        output(iBoot).variable = local_variable;
        output(iBoot).model = ModelList;
        %     output(iBoot).type = modelNames;
        toc
    end
    if nregion == 1
        output_V1 = output;
    else
        output_HPC = output;
    end
   
end
 save(fullfile(analysis_folder,'V1-HPC sleep reactivation','context_ripple_modulation_glme.mat'),'output_V1','output_HPC');
% z_FR_track_diff(contains(regions_all,'V1'))





%%%%%%%%%%% ALL ripples

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (Track prefering)';
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
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_ALL ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)



%%%%%%%%%%%%% Low ripple


all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_low_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');

nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (low ripples Track prefering)';
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
glme = fitlme(tbl_V1,'POST_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_V1,'PRE_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_V1,'SHIFT_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_HPC,'POST_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

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
glme = fitlme(tbl_HPC,'PRE_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','r');

xlim([-2 2])
ylim([-0.7 0.7])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)

subplot(2,3,6)
% 
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
glme = fitlme(tbl_HPC,'SHIFT_ripple_FR_diff_LOW ~ z_FR_track_diff + (1|subjectID)');

text(prctile(X,99.9), prctile(Y,99.9), ...
    sprintf('R^2 = %.3f\np = %.2e', glme.Rsquared.Adjusted, glme.Coefficients.pValue(2)), ...
    'FontSize',10,'Color','k');

xlim([-2 2])
ylim([-0.4 0.4])
xlabel('Track FR diff (z)')
ylabel('Ripple FR diff (z)')
set(gca,'TickDir','out','Box','off','FontSize',12)






%%%%%%%%%%%%% High ripple

all_PSTH_diff = vertcat(context_modulation_all.PSTH_diff_high_ripple{:});
POST_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>0&context_modulation_all.timebin<0.2),2,'omitnan');
PRE_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-0.2&context_modulation_all.timebin<0),2,'omitnan');
shifted_ripple_FR_diff = mean(all_PSTH_diff(:,context_modulation_all.timebin>-1&context_modulation_all.timebin<-0.8),2,'omitnan');


nfig = figure;
nfig.Name = 'Context selective ripple modulation in V1 and HPC regression (high ripples Track prefering)';
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
