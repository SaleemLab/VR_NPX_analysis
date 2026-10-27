function st = ar_stats(y, session, animal, varargin)
% Group statistics of a transition-level measure y (test of mean > 0, or of a
% paired difference when y is already a difference):
%   session means -> signed-rank test across sessions, sessions > 0
%   animal means (mean of its session means) -> animals > 0
%   LME y ~ 1 + (1|animal) + (1|animal:session) on transitions
% Optional 'covariates' (table with the same rows) are added as fixed effects;
% the intercept is still the tested quantity (centre covariates first).
p = inputParser;
addParameter(p, 'covariates', table());
parse(p, varargin{:});
ok = isfinite(y);
y = y(ok); session = session(ok); animal = animal(ok);
[us, ~, is] = unique(session);
sm = accumarray(is, y, [], @mean);
sa = accumarray(is, animal, [], @(x) x(1));
[ua, ~, ia] = unique(sa);
am = accumarray(ia, sm, [], @mean);
st.n_trans = numel(y);
st.sessions = us; st.session_mean = sm; st.session_animal = sa;
st.animals = ua; st.animal_mean = am;
st.mean = mean(sm); st.sem = std(sm) / sqrt(numel(sm));
if numel(sm) > 1, st.p_signrank = signrank(sm); else, st.p_signrank = NaN; end
st.n_sess_pos = sum(sm > 0); st.n_sess = numel(sm);
st.n_anim_pos = sum(am > 0); st.n_anim = numel(am);
tbl = table(y, categorical(session), categorical(animal), 'VariableNames', {'y', 'session', 'animal'});
formula = 'y ~ 1';
cov = p.Results.covariates;
if ~isempty(cov)
    cov = cov(ok, :);
    for v = cov.Properties.VariableNames
        tbl.(v{1}) = cov.(v{1});
        formula = [formula ' + ' v{1}]; %#ok<AGROW>
    end
end
formula = [formula ' + (1|animal) + (1|animal:session)'];
try
    lme = fitlme(tbl, formula);
    c = lme.Coefficients;
    st.lme_est = c.Estimate(1); st.lme_p = c.pValue(1);
    ci = coefCI(lme); st.lme_ci = ci(1, :);
    st.lme = lme;
catch err
    warning('LME failed: %s', err.message);
    st.lme_est = NaN; st.lme_p = NaN; st.lme_ci = [NaN NaN]; st.lme = [];
end
end
