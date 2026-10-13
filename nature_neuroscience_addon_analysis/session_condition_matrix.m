function Y = session_condition_matrix(T, value_field, sleep, window, analysis, labels)
% Sessions (by session_id index) x conditions; NaN where missing
Y = NaN(max(T.session_id), numel(labels));
for iCond = 1:numel(labels)
    Tc = T(T.sleep == sleep & T.V1_window == window & T.analysis == analysis & T.condition == labels(iCond), :);
    Y(Tc.session_id, iCond) = Tc.(value_field);
end
end
