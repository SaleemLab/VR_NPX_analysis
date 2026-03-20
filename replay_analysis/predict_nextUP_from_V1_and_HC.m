function V1_HC_nextUP_lme = predict_nextUP_from_V1_and_HC(tbl,V1_HC_nextUP_lme)
if isempty(V1_HC_ripple_lme)
    V1_HC_ripple_lme = struct;
    % --- Configuration ---
    nBoot = 1000;
    nShuff = 1000;
    unique_sessions = unique(tbl.session_id);

    % --- Initialize Flat Containers ---
    % 1:Persistence, 2:HPC_Driver, 3:Mediation_Loop, 4:Recruitment
    beta_boot = zeros(nBoot, 4);
    prop_boot = zeros(nBoot, 4);
    prop_mediated_boot = zeros(nBoot,1);

    % Triple Null Containers
    shuff_HPC_b = zeros(nShuff, 4);   shuff_HPC_p = zeros(nShuff, 4);
    shuff_V1PRE_b = zeros(nShuff, 4); shuff_V1PRE_p = zeros(nShuff, 4);
    shuff_LINK_b = zeros(nShuff, 4);  shuff_LINK_p = zeros(nShuff, 4);

    % Pre-group data by session for speed
    session_data = cell(length(unique_sessions), 1);
    for s = 1:length(unique_sessions)
        session_data{s} = tbl(tbl.session_id == unique_sessions(s), :);
    end

    %% --- PART 1: NESTED BOOTSTRAP (Reliability) ---
    fprintf('Running %d Bootstrap iterations...\n', nBoot);
    parfor i = 1:nBoot
        seed = RandStream('philox4x32_10', 'Seed', i);
        resamp_idx = randsample(seed, length(unique_sessions), length(unique_sessions), true);
        b_tbl = vertcat(session_data{resamp_idx});

        [betas, props,~,prop_mediated] = calculate_metrics_local(b_tbl);
        beta_boot(i,:) = betas;
        prop_boot(i,:) = props;
        prop_mediated_boot(i)=prop_mediated;
    end

    %% --- PART 2: TRIPLE SURGICAL SHUFFLE (Significance) ---
    fprintf('Running %d Triple Shuffle iterations...\n', nShuff);
    parfor i = 1:nShuff
        seed = RandStream('philox4x32_10', 'Seed', i);

        % Null 1: Shuffle HPC (Tests Driver)
        s_hpc = tbl;
        for s = 1:length(unique_sessions)
            idx = find(s_hpc.session_id == unique_sessions(s));
            s_hpc.lastRippleHPC(idx) = s_hpc.lastRippleHPC(idx(randperm(seed, length(idx))));
        end
        [b1, p1] = calculate_metrics_local(s_hpc);
        shuff_HPC_b(i,:) = b1; shuff_HPC_p(i,:) = p1;

        % Null 2: Shuffle V1PRE (Tests Persistence)
        s_v1 = tbl;
        for s = 1:length(unique_sessions)
            idx = find(s_v1.session_id == unique_sessions(s));
            s_v1.lastRippleV1PRE(idx) = s_v1.lastRippleV1PRE(idx(randperm(seed, length(idx))));
        end
        [b2, p2] = calculate_metrics_local(s_v1);
        shuff_V1PRE_b(i,:) = b2; shuff_V1PRE_p(i,:) = p2;

        % Null 3: Shuffle Link (Tests Recruitment & Loop)
        s_link = tbl;
        for s = 1:length(unique_sessions)
            idx = find(s_link.session_id == unique_sessions(s));
            s_link.lastRippleV1PRE(idx) = s_link.lastRippleV1PRE(idx(randperm(seed, length(idx))));
        end
        [b3, p3] = calculate_metrics_local(s_link);
        shuff_LINK_b(i,:) = b3; shuff_LINK_p(i,:) = p3;
    end

    %% --- PART 3: PACKAGING (Flat Structure) ---
    [obs_b, obs_p,pval] = calculate_metrics_local(tbl);

    % Metadata & Observed
    V1_HC_ripple_lme.labels = {'Persistence', 'HPC_Driver', 'Mediation_Loop', 'Recruitment'};
    V1_HC_ripple_lme.beta_observed = obs_b;
    V1_HC_ripple_lme.prop_observed = obs_p;
    V1_HC_ripple_lme.pval_observed = pval;

    mFull = fitlme(tbl, 'lastRippleV1 ~ lastRippleV1PRE + lastRippleHPC + (1|subject_id) + (1|session_id)');
    mA    = fitlme(tbl, 'lastRippleHPC ~ lastRippleV1PRE + (1|subject_id)+ (1|session_id)');
    mV1   = fitlme(tbl, 'lastRippleV1 ~ lastRippleV1PRE + (1|subject_id)+ (1|session_id)');
    mHPC  = fitlme(tbl, 'lastRippleV1 ~ lastRippleHPC + (1|subject_id)+ (1|session_id)');
    V1_HC_ripple_lme.models = {mFull,mA,mV1,mHPC};
    % Reliability (Distributions & CIs)
    V1_HC_ripple_lme.beta_boot = beta_boot;
    V1_HC_ripple_lme.prop_boot = prop_boot;
    V1_HC_ripple_lme.beta_CI = quantile(beta_boot, [0.025, 0.975]);
    V1_HC_ripple_lme.prop_CI = quantile(prop_boot, [0.025, 0.975]);

    % Significance (Surgical P-values)
    % 1. Persistence vs V1PRE-shuff
    V1_HC_ripple_lme.p_values_beta(1) = mean(abs(shuff_V1PRE_b(:,1)) >= abs(obs_b(1)));
    V1_HC_ripple_lme.p_values_prop(1) = mean(shuff_V1PRE_p(:,1) >= obs_p(1));

    % 2. Driver vs HPC-shuff
    V1_HC_ripple_lme.p_values_beta(2) = mean(abs(shuff_HPC_b(:,2)) >= abs(obs_b(2)));
    V1_HC_ripple_lme.p_values_prop(2) = mean(shuff_HPC_p(:,2) >= obs_p(2));

    % 3 & 4. Loop & Recruitment vs Link-shuff
    V1_HC_ripple_lme.p_values_beta(3) = mean(abs(shuff_LINK_b(:,3)) >= abs(obs_b(3)));
    V1_HC_ripple_lme.p_values_prop(3) = mean(shuff_LINK_p(:,3) >= obs_p(3));
    V1_HC_ripple_lme.p_values_beta(4) = mean(abs(shuff_LINK_b(:,4)) >= abs(obs_b(4)));
    V1_HC_ripple_lme.p_values_prop(4) = mean(shuff_LINK_p(:,4) >= obs_p(4));

    V1_HC_ripple_lme.shuff_LINK_b = shuff_LINK_b;
    V1_HC_ripple_lme.shuff_HPC_b = shuff_HPC_b;
    V1_HC_ripple_lme.shuff_V1PRE_b = shuff_V1PRE_b;
    V1_HC_ripple_lme.shuff_LINK_prop = shuff_LINK_p;
    V1_HC_ripple_lme.shuff_HPC_prop = shuff_HPC_p;
    V1_HC_ripple_lme.shuff_V1PRE_prop = shuff_V1PRE_p;

    V1_HC_ripple_lme.prop_mediated_boot = prop_mediated_boot;

end

shuffled_b = {V1_HC_ripple_lme.shuff_V1PRE_b,V1_HC_ripple_lme.shuff_HPC_b,V1_HC_ripple_lme.shuff_LINK_b,V1_HC_ripple_lme.shuff_V1PRE_b};
shuffled_prop = {V1_HC_ripple_lme.shuff_V1PRE_prop,V1_HC_ripple_lme.shuff_HPC_prop,V1_HC_ripple_lme.shuff_LINK_prop,V1_HC_ripple_lme.shuff_V1PRE_prop};

res = V1_HC_ripple_lme;
labels = res.labels;
colors = [0.2 0.6 0.8; 0.8 0.4 0.4; 0.6 0.4 0.8; 0.4 0.8 0.4];
gray_col = [0.7 0.7 0.7]; % Neutral gray for shuffles

figure('Color', 'w', 'Position', [100, 100, 1000, 500],'Name','Predict V1 from pre V1 and HC');

%% Panel 1: Beta Coefficients (Observed vs. Shuffled Null)
subplot(1,2,1); hold on;

for i = 1:3
    % --- 1. Main Bar (Observed) ---
    x_obs = i;
    bar(x_obs, res.beta_observed(i), 0.25, 'FaceColor', colors(i,:), 'EdgeColor', 'none');
    % Bootstrap CI
    low_b = res.beta_observed(i) - res.beta_CI(1,i);
    high_b = res.beta_CI(2,i) - res.beta_observed(i);
    errorbar(x_obs, res.beta_observed(i), low_b, high_b, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);

    % --- 2. Shuffle Bar (Null) ---
    x_shuff = i + 0.3;
    % We need the mean and 95% CI of the relevant shuffle distribution
    % (Logic: 1=V1PRE_shuff, 2=HPC_shuff, 3&4=LINK_shuff)
    % Note: You'll need to pass the shuff distributions into this function or extract them
    % For this example, I assume you extracted means/CIs during the main function run
    temp = shuffled_b{i}(:,i);
    shuff_mean = mean(temp); % mean
    shuff_low = prctile(temp,2.5); % Shuffle 2.5th percentile
    shuff_high = prctile(temp',97.5); %  Shuffle 97.5th percentile

    bar(x_shuff, shuff_mean, 0.25, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    errorbar(x_shuff, shuff_mean, shuff_mean-shuff_low, shuff_high-shuff_mean, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 4);

    % --- 3. Significance Stars ---
    p = res.p_values_beta(i);
    if p < 0.05
        txt = '*'; if p < 0.01, txt = '**'; end; if p < 0.001, txt = '***'; end
        txt = sprintf('%.3e',p);
        text(i+0.15, max(res.beta_CI(2,i), shuff_high) + 0.02, txt, ...
            'FontSize', 10);
    end
end

ylabel('Standardized Coefficient (\beta)');
set(gca, 'XTick', 1.15:1:4.15, 'XTickLabel', labels, 'XTickLabelRotation', 30);
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off
title('Effect Size: Observed vs. Shuffled Null');
yline(0, 'k-');

%% Panel 2: Proportions (Information Content)
subplot(1,2,2); hold on;
for i = 1:3
    if i == 1,     s_dist = res.shuff_V1PRE_prop(:,i);
    elseif i == 2, s_dist = res.shuff_HPC_prop(:,i);
    else,          s_dist = res.shuff_LINK_prop(:,i);
    end

    % Observed Bar
    bar(i - 0.15, res.prop_observed(i), 0.25, 'FaceColor', colors(i,:), 'EdgeColor', 'none');
    errorbar(i - 0.15, res.prop_observed(i), ...
        res.prop_observed(i) - res.prop_CI(1,i), ...
        res.prop_CI(2,i) - res.prop_observed(i), 'k', 'linestyle','none','LineWidth',1.2);

    % Shuffle Bar (Gray)
    % s_mean = mean(s_dist);
    % s_ci = quantile(s_dist, [0.025, 0.975]);
    % bar(i + 0.15, s_mean, 0.3, 'FaceColor', gray_col, 'EdgeColor', 'none');
    % errorbar(i + 0.15, s_mean, s_mean - s_ci(1), s_ci(2) - s_mean, ...
    %     'Color', [0.7 0.7 0.7], 'linestyle','none','LineWidth',1.2);
end
xline(3.5,'--')
ylim([-1 100])
ylabel('Prop. Explained Variance (%)');
set(gca, 'XTick', 1:4, 'XTickLabel', labels, 'XTickLabelRotation', 30);
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off
title('Variance Proportions: Observed vs. Shuffle');

sgtitle('V1-HPC Information Pathways: Significance vs. Noise', 'FontSize', 16);


end % End of main function

%% --- Local Helper Function ---
function [betas, props,pval,prop_mediated] = calculate_metrics_local(data)
mFull = fitlme(data, 'lastRippleV1 ~ lastRippleV1PRE + lastRippleHPC + (1|subject_id)+ (1|session_id)');
mA    = fitlme(data, 'lastRippleHPC ~ lastRippleV1PRE + (1|subject_id)+ (1|session_id)');
mV1   = fitlme(data, 'lastRippleV1 ~ lastRippleV1PRE + (1|subject_id)+ (1|session_id)');
mHPC  = fitlme(data, 'lastRippleV1 ~ lastRippleHPC + (1|subject_id)+ (1|session_id)');

a = mA.Coefficients.Estimate(2);  % V1_PRE -> HPC
b = mFull.Coefficients.Estimate(3);  % HPC -> V1
c_p = mFull.Coefficients.Estimate(2); % V1_PRE -> V1 (Direct)

indirect_effect = a * b;
total_effect = c_p + indirect_effect; % Total Effect (c)
if total_effect ~= 0
    prop_mediated = (indirect_effect / total_effect) * 100;
else
    prop_mediated = 0;
end


betas = [c_p, b, a*b, a];
pval = [mFull.Coefficients.pValue(2),mFull.Coefficients.pValue(3),nan,mA.Coefficients.pValue(2)];

rf = mFull.Rsquared.Adjusted; rv = mV1.Rsquared.Adjusted; rh = mHPC.Rsquared.Adjusted;
uv = rf - rh; uh = rf - rv; sh = (rv + rh) - rf;

total = uv + uh + sh;
props = [uv, uh, sh, mA.Rsquared.Adjusted] ./ [total, total, total, 1] * 100;
end