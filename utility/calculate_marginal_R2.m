function [marginal_R2,conditional_R2] = calculate_marginal_R2(tbl,lme)

% 1. Fixed Effects Variance (sigma_f^2)
% Predict using ONLY fixed effects
pred_fixed = predict(lme, tbl, 'Conditional', false); 
var_f = var(pred_fixed,'omitnan');

% 2. Random Effects Variance (sigma_a^2) and Residual Variance (sigma_e^2)
[psi, sigma2_e] = covarianceParameters(lme);
var_a = sum(cellfun(@(x) sum(diag(x)), psi)); % Variance of random intercepts

% 3. Total Variance (The denominator)
var_total = var_f + var_a + sigma2_e;

% 4. Marginal R2 (Fixed only)
marginal_R2 = var_f / var_total;

% 5. Conditional R2 (Fixed + Random)
conditional_R2 = (var_f + var_a) / var_total;