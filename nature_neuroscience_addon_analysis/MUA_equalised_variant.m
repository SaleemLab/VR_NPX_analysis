function v = MUA_equalised_variant(name)
% Settings of the two versions of the MUA-equalised event control (main_04), each
% saved in its own folder under nature_neuroscience_revision_analysis:
%   'dominant'         : Low and High equalised for HC MUA and dominant-hemisphere
%                        V1 MUA (nQ x nQ joint bins) -> MUA_equalised_event_control
%   'with_nondominant' : also for non-dominant-hemisphere V1 MUA (nQ^3 joint bins)
%                        -> MUA_equalised_event_control_with_nondominant
% v.folder     : output folder
% v.key        : name in plot_AUC_low_vs_high_controls / plot_temporal_AUC_low_vs_high
% v.label      : matched measures, for figure titles
% v.vars       : matched MUA measures (columns of v.matrix), v.var_labels axis labels
% v.nQ         : quantile bins per measure
% v.matrix     : X = v.matrix(HC_MUA, V1_L, V1_R, hemi), events x measures, with
%                dominant hemisphere hemi 1 (V1 L) or 2 (V1 R); NaN where hemi is NaN
    if nargin < 1, name = 'dominant'; end
    switch name
        case 'dominant'
            v.folder = 'MUA_equalised_event_control';
            v.key = 'MUA_equalised';
            v.label = 'within-session HC and dominant V1 MUA matched';
            v.vars = {'HC_MUA', 'V1_MUA_dom'};
            v.var_labels = {'HC MUA 0-0.1 s (z)', 'Dominant V1 MUA (z)'};
        case 'with_nondominant'
            v.folder = 'MUA_equalised_event_control_with_nondominant';
            v.key = 'MUA_equalised_with_nondominant';
            v.label = 'within-session HC, dominant and non-dominant V1 MUA matched';
            v.vars = {'HC_MUA', 'V1_MUA_dom', 'V1_MUA_nondom'};
            v.var_labels = {'HC MUA 0-0.1 s (z)', 'Dominant V1 MUA (z)', 'Non-dominant V1 MUA (z)'};
        otherwise
            error('MUA_equalised_variant:name', 'Unknown variant %s', name);
    end
    v.name = name;
    v.nQ = 5;
    vars = v.vars;
    v.matrix = @(HC_MUA, V1_L, V1_R, hemi) match_matrix(vars, HC_MUA, V1_L, V1_R, hemi);
end

function X = match_matrix(vars, HC_MUA, V1_L, V1_R, hemi)
    dom = NaN(numel(hemi), 1); nondom = dom;
    dom(hemi == 1) = V1_L(hemi == 1);    dom(hemi == 2) = V1_R(hemi == 2);
    nondom(hemi == 1) = V1_R(hemi == 1); nondom(hemi == 2) = V1_L(hemi == 2);
    M = struct('HC_MUA', HC_MUA(:), 'V1_MUA_dom', dom, 'V1_MUA_nondom', nondom);
    X = cell2mat(cellfun(@(f) M.(f), vars, 'UniformOutput', false));
end
