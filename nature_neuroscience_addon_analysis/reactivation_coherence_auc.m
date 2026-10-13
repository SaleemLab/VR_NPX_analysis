function r = reactivation_coherence_auc(sel, meas, nBoot, c1, c2)
% Reactivation coherence AUC.
% sel  : selector bias per event (e.g. HPC log odds, or V1 log odds)
% meas : measured bias per event in the other region
% Track 1 events: sel >= th (& c1), Track 2 events: sel <= -th (& c2), for
% thresholds th = 0:10:90 prctile of |sel| over all given events. AUC is the
% normalised trapz of mean(meas | T1) - mean(meas | T2) across thresholds.
% c1 / c2 : optional logical per event (e.g. hemisphere-matched spindle/SO
%           condition: c1 on hemisphere 2, c2 on hemisphere 1, as in POST).
% nBoot = 0 : point estimate on the observed events (no bootstrap).
% nBoot > 0 : bootstrap over events (seeded); shuffled control keeps the
%             selector (and its conditions) un-resampled while resampling
%             meas, which breaks the selector-measure pairing.

sel = sel(:)'; meas = meas(:)';
n_ev = numel(sel);
if nargin < 4 || isempty(c1), c1 = true(1, n_ev); end
if nargin < 5 || isempty(c2), c2 = true(1, n_ev); end
c1 = logical(c1(:)'); c2 = logical(c2(:)');

r = struct('n_events', n_ev, 'AUC_mean', NaN, 'AUC_CI', [NaN NaN], ...
    'AUC_mean_shuffled', NaN, 'AUC_CI_shuffled', [NaN NaN]);
if n_ev < 10
    return
end

th = prctile(abs(sel), 0:10:100); th = th(1:end-1); nT = length(th);

if nBoot == 0
    r.AUC_mean = trapz(th, bias_diff(sel, meas, c1, c2, th)) / (max(th) - min(th));
    return
end

b_d = NaN(nBoot, nT); b_s = NaN(nBoot, nT);
parfor iBoot = 1:nBoot
    s = RandStream('philox4x32_10', 'Seed', iBoot);
    idx = randi(s, n_ev, n_ev, 1);
    b_d(iBoot, :) = bias_diff(sel(idx), meas(idx), c1(idx), c2(idx), th);
    b_s(iBoot, :) = bias_diff(sel, meas(idx), c1, c2, th);
end

au = (trapz(th, b_d') / (max(th) - min(th)))';
aus = (trapz(th, b_s') / (max(th) - min(th)))';
r.AUC_mean = mean(au, 'omitnan');
r.AUC_CI = prctile(au, [2.5 97.5]);
r.AUC_mean_shuffled = mean(aus, 'omitnan');
r.AUC_CI_shuffled = prctile(aus, [2.5 97.5]);
end

function d = bias_diff(sel, meas, c1, c2, th)
d = NaN(1, numel(th));
for j = 1:numel(th)
    t1 = sel >= th(j) & c1;
    t2 = sel <= -th(j) & c2;
    if any(t1) && any(t2)
        d(j) = mean(meas(t1), 'omitnan') - mean(meas(t2), 'omitnan');
    end
end
end
