function labels = session_condition_labels(T, analysis)
% Condition labels of an analysis, ordered so the reference condition
% (Low / Peak / Bilateral) comes first: one-sided tests are condition 2 > condition 1
% (High > Low, Trough > Peak, Unilateral Trough > Bilateral Trough).
labels = unique(T.condition(T.analysis == analysis), 'stable');
if numel(labels) == 2 && is_reference_condition(labels(2)) && ~is_reference_condition(labels(1))
    labels = labels([2 1]);
end
end
