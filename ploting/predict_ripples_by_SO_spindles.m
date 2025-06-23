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
    figure;
    subplot(3,1,1); histogram(output.spindle.pref_phase); title('Spindle preferred phase');
    subplot(3,1,2); histogram(output.SO.pref_phase); title('SO preferred phase');
    subplot(3,1,3); histogram(output.combined.pref_phase); title('Combined preferred phase');
end
