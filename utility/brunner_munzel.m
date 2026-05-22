function [p_value, stat, effect_size] = brunner_munzel(x, y, alpha, tail)
    % BRUNNER_MUNZEL Performs the Brunner-Munzel test for stochastic equality.
    % This test does not assume equal variances or shapes between groups.
    %
    % Inputs:
    %   x          - Vector of data from sample 1 (e.g., Orange)
    %   y          - Vector of data from sample 2 (e.g., Blue)
    %   alpha      - Significance level (default: 0.05)
    %   tail       - Alternative hypothesis:
    %                'both'  : Median(X) ~= Median(Y) (default)
    %                'right' : Median(X) < Median(Y)  (Y is shifted right)
    %                'left'  : Median(X) > Median(Y)  (Y is shifted left)
    %
    % Outputs:
    %   p_value    - Computed p-value based on the specified tail
    %   stat       - The W statistic (follows t-distribution)
    %   effect_size- Relative effect size P(X < Y) + 0.5*P(X = Y)

    % Handle default inputs
    if nargin < 3 || isempty(alpha)
        alpha = 0.05;
    end
    if nargin < 4 || isempty(tail)
        tail = 'both';
    end

    % Ensure columns and remove NaNs
    x = x(~isnan(x)); x = x(:);
    y = y(~isnan(y)); y = y(:);
    
    n = length(x);
    m = length(y);
    
    % 1. Combine data and calculate overall ranks
    combined = [x; y];
    R_combined = tiedrank(combined);
    
    R_x = R_combined(1:n);
    R_y = R_combined(n+1:end);
    
    % Mean ranks
    mean_R_x = mean(R_x);
    mean_R_y = mean(R_y);
    
    % 2. Internal ranks
    R_x_internal = tiedrank(x);
    R_y_internal = tiedrank(y);
    
    % 3. Calculate rank variances (unbiased estimators)
    S_x_sq = (1 / (n - 1)) * sum((R_x - R_x_internal - mean_R_x + (n + 1)/2).^2);
    S_y_sq = (1 / (m - 1)) * sum((R_y - R_y_internal - mean_R_y + (m + 1)/2).^2);
    
    % 4. Relative Effect Size (P(X < Y) + 0.5 * P(X = Y))
    effect_size = (mean_R_y - mean_R_x) / (n + m) + 0.5;
    
    % 5. Test Statistic W
    sigma_x_sq = S_x_sq / (m^2);
    sigma_y_sq = S_y_sq / (n^2);
    pooled_sigma = sqrt((n + m) * (sigma_x_sq / n + sigma_y_sq / m));
    
    if pooled_sigma == 0
        error('Brunner-Munzel test statistic is undefined because variance is zero.');
    end
    
    stat = (mean_R_y - mean_R_x) / (pooled_sigma * (n + m));
    
    % 6. Degrees of Freedom (Welch–Satterthwaite-like approximation for ranks)
    df_numerator = (S_x_sq / m + S_y_sq / n)^2;
    df_denominator = ((S_x_sq / m)^2 / (n - 1)) + ((S_y_sq / n)^2 / (m - 1));
    df = df_numerator / df_denominator;
    
    % 7. Compute p-value based on the selected tail
    switch lower(tail)
        case 'both'
            % Two-tailed test
            p_value = 2 * (1 - tcdf(abs(stat), df));
            
        case 'right'
            % One-tailed test: Is Y (Blue) shifted right compared to X (Orange)?
            % This means we expect stat > 0 and effect_size > 0.5
            p_value = 1 - tcdf(stat, df);
            
        case 'left'
            % One-tailed test: Is Y (Blue) shifted left compared to X (Orange)?
            % This means we expect stat < 0 and effect_size < 0.5
            p_value = tcdf(stat, df);
            
        otherwise
            error('Invalid tail option. Use ''both'', ''right'', or ''left''.');
    end
end