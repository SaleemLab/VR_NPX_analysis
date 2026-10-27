function pairwise_coupling_oscillation_POST(varargin)
% Pairwise V1-HC coupling vs shared track preference (as pairwise_V1_HC_coupling_POST,
% gain control 'total_V1') computed separately in two event subsets per oscillation
% condition. Spindle and SO conditions use each V1 cell's OWN hemisphere (no dominant
% hemisphere, no decoder): a V1_L cell's events are split by the left-hemisphere feature,
% a V1_R cell's by the right one (ripple_info column 1 = left, 2 = right).
% Conditions (A -> B; B is the condition expected to show stronger coupling):
%   ripple_power     : within-session ripple power Q1 -> Q4 (same events for all cells)
%   spindle          : own-hemisphere spindle_percentile <= 25 -> > 75 (bottom 10% NaN,
%                      so Low = 10-25th, as hemisphere_conditions)
%   SO_phase         : own-hemisphere SO peak [-pi/2, pi/2] -> trough
%   SO_trough_uni_bi : own hemisphere in trough with the other hemisphere also in trough
%                      (bilateral) -> other hemisphere in peak (unilateral)
% For each condition, V1 window (pre -0.2-0 s, post 0-0.2 s; HC 0-0.1 s) and subset:
%   r   : partial correlation per pair (leave-one-out total V1 MUA and HC MUA), computed
%         within the subset's events
%   rho : Spearman correlation across pairs of r with z_V1 * z_HC, minus the mean of its
%         permutation null (preferences permuted within region)
% The same pairs enter A and B: a cell must fire in >= min_active events in both subsets.
% Firing-rate control (binomial thinning): for every cell, spikes of the subset with the
% higher mean count are kept with probability (lower mean / higher mean), so both subsets
% have the same mean count per cell; covariates stay unthinned. Repeated nThin times
% (seeded); rho and its null are averaged over repeats.
% hemisphere_mode 'dominant' (paper-equivalent check): spindle and SO use the event's
% dominant hemisphere instead, from the sign of the original PLS-KDE V1 log odds in the
% same V1 window (V1 > 0 -> right / column 2, V1 < 0 -> left / column 1, as
% hemisphere_conditions); every V1 cell then uses the same event subsets.
% hemisphere_mode 'preference': each V1 cell uses the hemisphere that represents its
% preferred track (z_track_diff > 0, Track 1 / L -> right; < 0 -> left), wherever the
% cell is; identical to 'own' except for cells preferring the ipsilateral track.
% Options: 'sessions' (all), 'nPerm' (1000), 'nThin' (20), 'nPermThin' (50 per repeat),
%   'min_active' (10), 'hemisphere_mode' ('own' | 'dominant' | 'preference'), 'counts_name' (counts
%   subfolder, default 'counts'), 'out_name' (default 'oscillation_effects')
%   'gain_folder' (default 'gain_total_V1'; 'gain_total_V1_all_cells' with counts_all_rates)
%   'measure' : 'correlation' (default; Spearman rho across pairs) or 'covariance' (pair
%               partial covariance in spikes^2 per window; statistic = slope across pairs on
%               z_V1 * z_HC, minus null, which is sensitive to coupling amplitude; thinning
%               matters more here since covariance grows with firing rate)
% Output: pairwise_V1_HC_coupling\<gain_folder>\<out_name>\
%   oscillation_coupling_POST.mat (session table and options), session_stats.csv

p = inputParser;
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'nPerm', 1000, @isnumeric);
addParameter(p, 'nThin', 20, @isnumeric);
addParameter(p, 'nPermThin', 50, @isnumeric);
addParameter(p, 'min_active', 10, @isnumeric);
addParameter(p, 'hemisphere_mode', 'own', @(s) ismember(s, {'own', 'dominant', 'preference'}));
addParameter(p, 'counts_name', 'counts', @ischar);
addParameter(p, 'out_name', 'oscillation_effects', @ischar);
addParameter(p, 'gain_folder', 'gain_total_V1', @ischar);   % e.g. gain_total_V1_all_cells
addParameter(p, 'measure', 'correlation', @(s) ismember(s, {'correlation', 'covariance'}));
parse(p, varargin{:});
P = p.Results;

rev_dir = addon_setup();
base_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling');
out_dir = fullfile(base_dir, P.gain_folder, P.out_name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
wins = {'V1_pre', 'V1_post'};
win_range = {[-0.2 0], [0 0.2]};
V1_bias = [];
if strcmp(P.hemisphere_mode, 'dominant')
    L = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    V1 = (L.KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + L.KDE_reactivation_ripples_PSTH.nan_mask)';
    clear L
    bc = -1 + 0.005:0.01:1;
    V1_bias = cell2mat(cellfun(@(r) mean(V1(bc > r(1) & bc < r(2), :), 'omitnan')', win_range, 'UniformOutput', false));
    clear V1
end

T = table();
for nsession = P.sessions
    load(fullfile(base_dir, P.counts_name, sprintf('session_%02d.mat', nsession)), 'C');
    rows = C.ripple_info_rows;
    hemi = 1 + (C.V1.region == "V1_R");          % 1 = left, 2 = right
    if strcmp(P.hemisphere_mode, 'preference')
        hemi = 1 + (C.V1.z_track_diff > 0);      % hemisphere of the preferred track
    end
    h = C.HPC.counts.HC; H = C.MUA.HPC.HC;
    for w = 1:numel(wins)
        dom = [];
        if strcmp(P.hemisphere_mode, 'dominant')
            dom = NaN(numel(rows), 1);
            dom(V1_bias(rows, w) < 0) = 1; dom(V1_bias(rows, w) > 0) = 2;
        end
        conds = condition_masks(ripple_info, rows, dom);
        v = C.V1.counts.(wins{w});
        V1_tot = C.MUA.V1_L.(wins{w}) + C.MUA.V1_R.(wins{w});
        for c = 1:numel(conds)
            seed = nsession * 1000 + w * 100 + c;
            R = condition_coupling(v, h, V1_tot, H, hemi, C.V1.z_track_diff, C.HPC.z_track_diff, ...
                conds(c).mask, P, seed);
            T = [T; table(nsession, C.subject_id, string(wins{w}), string(conds(c).name), ...
                string(conds(c).labels{1}), string(conds(c).labels{2}), R.n_pairs, R.n_events(1), R.n_events(2), ...
                R.val(1), R.val(2), R.val(2) - R.val(1), R.val_thin(1), R.val_thin(2), R.val_thin(2) - R.val_thin(1), ...
                R.rate_V1(1), R.rate_V1(2), R.rate_HC(1), R.rate_HC(2), ...
                'VariableNames', {'session_id', 'animal_id', 'V1_window', 'condition', 'label_A', 'label_B', ...
                'n_pairs', 'n_events_A', 'n_events_B', 'rho_A', 'rho_B', 'delta', ...
                'rho_A_thinned', 'rho_B_thinned', 'delta_thinned', ...
                'V1_count_A', 'V1_count_B', 'HC_count_A', 'HC_count_B'})]; %#ok<AGROW>
            fprintf('Session %2d %s %-16s: %5d pairs, rho %s %.3f -> %s %.3f | thinned %.3f -> %.3f\n', ...
                nsession, wins{w}, conds(c).name, R.n_pairs, conds(c).labels{1}, R.val(1), ...
                conds(c).labels{2}, R.val(2), R.val_thin(1), R.val_thin(2));
        end
    end
end
save(fullfile(out_dir, 'oscillation_coupling_POST.mat'), 'T', 'P', '-v7.3');
writetable(T, fullfile(out_dir, 'session_stats.csv'));
fprintf('Saved %s\n', out_dir);
end

function conds = condition_masks(ripple_info, rows, dom)
% mask: events x 2 hemispheres x 2 subsets (A, B), logical. With dom (events x 1, 1/2/NaN)
% the features are taken from the dominant hemisphere (column 1) and the other hemisphere
% (column 2) of each event and the event masks are used for both V1 hemispheres.
q = within_session_quartile(ripple_info.ripple_power(rows));
sp = ripple_info.spindle_percentile(rows, :);
ph = ripple_info.SO_phase(rows, :);
if ~isempty(dom)
    ok = ~isnan(dom); idx = (1:numel(rows))';
    d = sub2ind(size(sp), idx(ok), dom(ok)); o = sub2ind(size(sp), idx(ok), 3 - dom(ok));
    sp_d = NaN(numel(rows), 2); ph_d = NaN(numel(rows), 2);
    sp_d(ok, :) = [sp(d) sp(o)]; ph_d(ok, :) = [ph(d) ph(o)];
    sp = repmat(sp_d(:, 1), 1, 2);                   % dominant-hemisphere spindle for all cells
    ph = [ph_d(:, 1) ph_d(:, 1)]; ph_other = [ph_d(:, 2) ph_d(:, 2)];
else
    ph_other = ph(:, [2 1]);
end
peak = ph >= -pi/2 & ph <= pi/2;
trough = (ph >= -pi & ph < -pi/2) | (ph > pi/2 & ph <= pi);
conds = struct('name', {}, 'labels', {}, 'mask', {});
conds(1).name = 'ripple_power'; conds(1).labels = {'Low (Q1)', 'High (Q4)'};
conds(1).mask = cat(3, repmat(q == 1, 1, 2), repmat(q == 4, 1, 2));
conds(2).name = 'spindle'; conds(2).labels = {'Low (Q1)', 'High (Q4)'};
conds(2).mask = cat(3, sp >= 0 & sp <= 25, sp > 75 & sp <= 100);
conds(3).name = 'SO_phase'; conds(3).labels = {'Peak', 'Trough'};
conds(3).mask = cat(3, peak, trough);
conds(4).name = 'SO_trough_uni_bi'; conds(4).labels = {'Bilateral trough', 'Unilateral trough'};
peak_o = ph_other >= -pi/2 & ph_other <= pi/2;
trough_o = (ph_other >= -pi & ph_other < -pi/2) | (ph_other > pi/2 & ph_other <= pi);
conds(4).mask = cat(3, trough & trough_o, trough & peak_o);
end

function R = condition_coupling(v, h, V1_tot, H, hemi, zV1, zHC, mask, P, seed)
% Pair coupling in subsets A and B (V1 cells use the mask column of their hemisphere),
% rho - null per subset, unthinned and thinned
rng(seed, 'twister');
pairs_v = []; pairs_h = []; r = {[], []}; r_thin = cell(P.nThin, 2);
R.n_events = [0 0]; R.rate_V1 = [0 0]; R.rate_HC = [0 0]; n_rate_V1 = 0; n_rate_HC = 0;
for hh = 1:2
    iv = find(hemi == hh);
    E = {mask(:, hh, 1), mask(:, hh, 2)};
    R.n_events = R.n_events + cellfun(@sum, E) / 2;      % mean over hemispheres
    if isempty(iv) || any(cellfun(@sum, E) < P.min_active), continue; end
    ok_v = all(cell2mat(cellfun(@(e) sum(v(e, iv) > 0, 1)' >= P.min_active, E, 'UniformOutput', false)), 2);
    ok_h = all(cell2mat(cellfun(@(e) sum(h(e, :) > 0, 1)' >= P.min_active, E, 'UniformOutput', false)), 2);
    iv = iv(ok_v); ih = find(ok_h);
    if isempty(iv) || isempty(ih), continue; end
    [a, b] = ndgrid(iv, ih);
    pairs_v = [pairs_v; a(:)]; pairs_h = [pairs_h; b(:)]; %#ok<AGROW>
    m_v = cellfun(@(e) mean(v(e, iv), 1), E, 'UniformOutput', false);   % mean count per cell
    m_h = cellfun(@(e) mean(h(e, ih), 1), E, 'UniformOutput', false);
    R.rate_V1 = R.rate_V1 + [sum(m_v{1}) sum(m_v{2})]; n_rate_V1 = n_rate_V1 + numel(iv);
    R.rate_HC = R.rate_HC + [sum(m_h{1}) sum(m_h{2})]; n_rate_HC = n_rate_HC + numel(ih);
    for k = 1:2
        X = v(E{k}, iv); Y = h(E{k}, ih);
        Z1 = V1_tot(E{k}) - X; Z2 = H(E{k}) - Y;              % leave-one-out, unthinned
        rk = coupling_matrix(X, Y, Z1, Z2, P.measure);
        r{k} = [r{k}; rk(:)];
        % thinning: keep probability for subset k = min(1, other mean / this mean)
        pv = min(1, m_v{3 - k} ./ m_v{k}); ph = min(1, m_h{3 - k} ./ m_h{k});
        for t = 1:P.nThin
            Xt = binornd(X, repmat(pv, size(X, 1), 1));
            Yt = binornd(Y, repmat(ph, size(Y, 1), 1));
            rk = coupling_matrix(Xt, Yt, Z1, Z2, P.measure);
            r_thin{t, k} = [r_thin{t, k}; rk(:)];
        end
    end
end
R.n_pairs = numel(pairs_v);
R.rate_V1 = R.rate_V1 / max(n_rate_V1, 1); R.rate_HC = R.rate_HC / max(n_rate_HC, 1);
R.val = [NaN NaN]; R.val_thin = [NaN NaN];
if R.n_pairs < 3, return; end
for k = 1:2
    R.val(k) = rho_minus_null(r{k}, pairs_v, pairs_h, zV1, zHC, P.nPerm, P.measure);
    R.val_thin(k) = mean(cellfun(@(x) rho_minus_null(x, pairs_v, pairs_h, zV1, zHC, P.nPermThin, P.measure), ...
        r_thin(:, k)), 'omitnan');
end
end

function m = coupling_matrix(X, Y, Z1, Z2, measure)
% Pair coupling: partial correlation or partial covariance
[r, c] = pairwise_partial_corr_matrix(X, Y, Z1, Z2);
if strcmp(measure, 'covariance'), m = c; else, m = r; end
end

function val = rho_minus_null(r, pv, ph, zV1, zHC, nPerm, measure)
% Preference statistic over pairs minus its mean with z_track_diff permuted among the V1
% cells and among the HPC cells that enter the pairs:
%   correlation : Spearman rho(r, z_V1 * z_HC) (scale-free)
%   covariance  : least-squares slope of r on z_V1 * z_HC (in the units of r, so it
%                 grows with the amplitude of the coupling)
ok = isfinite(r);
r = r(ok); pv = pv(ok); ph = ph(ok);
if numel(r) < 3, val = NaN; return; end
if strcmp(measure, 'covariance')
    rr = r - mean(r);
    stat = @(pz) sum(rr .* (pz - mean(pz))) / sum((pz - mean(pz)) .^ 2);
else
    rr = tiedrank(r);
    stat = @(pz) corr(rr, tiedrank(pz));
end
rho = stat(zV1(pv) .* zHC(ph));
uv = unique(pv); uh = unique(ph);
null = NaN(nPerm, 1);
for k = 1:nPerm
    zv = zV1; zh = zHC;
    zv(uv) = zV1(uv(randperm(numel(uv))));
    zh(uh) = zHC(uh(randperm(numel(uh))));
    null(k) = stat(zv(pv) .* zh(ph));
end
val = rho - mean(null);
end

function q = within_session_quartile(x)
% Quartile (1-4) within the session (as reactivation_AUC_coherence_PRE)
q = NaN(size(x));
pcts = prctile(x, 0:25:100);
for iBin = 1:4
    if iBin == 1
        in_q = x >= pcts(iBin) & x <= pcts(iBin + 1);
    else
        in_q = x > pcts(iBin) & x <= pcts(iBin + 1);
    end
    q(in_q) = iBin;
end
end
