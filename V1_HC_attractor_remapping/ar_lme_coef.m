function st = ar_lme_coef(lme, name)
% Estimate, 95% CI, p and n of one fixed-effect term of a fitted LME
c = lme.Coefficients;
k = strcmp(c.Name, name);
ci = coefCI(lme);
st.lme = lme;
st.est = c.Estimate(k); st.p = c.pValue(k); st.ci = ci(k, :);
st.n = lme.NumObservations;
end
