h = findobj(gca, 'Type', 'Line');yData = {h.YData};
n = length(h);
firingRates = zeros(2, n);

% 3. Loop through and extract
for i = 1:n
    % yData(1) is usually the 'Peak' and yData(2) is 'Trough' 
    % (or vice versa depending on plot order)
    firingRates(:, i) = h(i).YData'; 
end

mean(firingRates,2)
std(firingRates')'./sqrt(length(firingRates))


ripple_info.SO_phase(:,1)


[rho, pval] = circ_corrcc(ripple_info.SO_phase(:,1), ripple_info.SO_phase(:,2));


sum(ripple_info.SO_phase(:,1) > -pi/2 & ripple_info.SO_phase(:,1) < pi/2 & (ripple_info.SO_phase(:,2) < -pi/2 | ripple_info.SO_phase(:,2) > pi/2))...
    +sum(ripple_info.SO_phase(:,2) > -pi/2 & ripple_info.SO_phase(:,2) < pi/2 & (ripple_info.SO_phase(:,1) < -pi/2 | ripple_info.SO_phase(:,1) > pi/2))


sum((ripple_info.SO_phase(:,1) < -pi/2 | ripple_info.SO_phase(:,1) > pi/2) & (ripple_info.SO_phase(:,2) < -pi/2 | ripple_info.SO_phase(:,2) > pi/2))


length(ripple_info.SO_phase)

log_odds_modulation_lme
ripple_power_modulation_lme
glme_models = spindle_power_modulation_lme
ripple_modulation_lme
% Initialize the master cell array
finalData = {};
glme_models = PRE_POST_V1_log_odds_lme;


% Loop through each model in the SO_phase_modulation_lme struct
for i = 1:length(glme_models)
    % Extracting data for current model
    vars = glme_models(i).variable; 
    pVals = glme_models(i).p;
    betas = glme_models(i).b;
    tStats = glme_models(i).t;
    CIs = glme_models(i).b_CI; % Nx2 matrix
    r2 = glme_models(i).R2;
    
    % Extract LCI (column 1) and HCI (column 2)
    LCI = CIs(:, 1);
    HCI = CIs(:, 2);
    
    % Repeat the Model index and R2 for each variable row
    modelIdx = repmat(i, length(vars), 1);
    r2Col = repmat(r2, length(vars), 1);
    
    % Concatenate columns into a temporary cell
    % We convert numeric arrays to cells to combine with 'vars' which is a cell
    tempData = [num2cell(modelIdx), vars, num2cell(pVals), ...
                num2cell(betas), num2cell(LCI), num2cell(HCI), ...
                num2cell(tStats), num2cell(r2Col)];
            
    % Append to the main list
    finalData = [finalData; tempData];
end

% Create the formal Table
outputTable = cell2table(finalData, 'VariableNames', ...
    {'ModelNum', 'Variable', 'P_Value', 'Beta', 'LCI', 'HCI', 't_Stat', 'Total_R2'});

% Display top results
head(outputTable)