function PSTH = cap_log_odds_PSTH(PSTH)
% Cap +/-Inf in every log-odds field of a KDE ripple PSTH struct.
% +Inf -> 99.5th prctile, -Inf -> 0.5th prctile of finite values.
% Percentiles are computed after applying the matching nan_mask (same as
% corticohippocampal_reactivation_coherence_KDE.m), but the saved field
% itself is left unmasked so downstream code can still add nan_mask.

fields = fieldnames(PSTH);
for i = 1:length(fields)
    f = fields{i};
    if ~contains(f, 'logodds')
        continue
    end

    if contains(f, 'awake')
        mask = PSTH.nan_mask_awake;
    else
        mask = PSTH.nan_mask;
    end

    vals = PSTH.(f);
    masked = vals + mask;
    finite_vals = masked(isfinite(masked));

    vals(vals >= inf) = prctile(finite_vals, 99.5);
    vals(vals <= -inf) = prctile(finite_vals, 0.5);
    PSTH.(f) = vals;

    fprintf('%s: capped at [%.3f, %.3f]\n', f, prctile(finite_vals, 0.5), prctile(finite_vals, 99.5));
end
end
