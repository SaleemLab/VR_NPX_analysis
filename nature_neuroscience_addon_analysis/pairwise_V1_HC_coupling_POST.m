function pairwise_V1_HC_coupling_POST(gain_control, varargin)
% Pairwise V1-HC coupling across POST ripples and its relation to track preference
% (R1 summary: "interactions between V1 and CA1 neurons"; R4.2).
% For every V1 x HPC neuron pair of a session (counts from extract_pairwise_counts_POST):
%   r_raw     : Pearson correlation across ripple events of the V1 count (pre -0.2-0 s or
%               post 0-0.2 s) and the HPC count (0-0.1 s)
%   r_partial : the same after regressing both counts on the gain covariates Z
%               (leave-one-out: the V1 cell's own spikes are removed from the V1 MUA and
%               the HPC cell's own spikes from the HPC MUA)
% gain_control:
%   'total_V1'       : Z = [V1_L + V1_R MUA, HPC MUA]. Removes global cortical state but
%                      keeps the left-right hemisphere difference (the track signal).
%   'own_hemisphere' : Z = [MUA of the V1 cell's own hemisphere, HPC MUA]. Also removes
%                      hemisphere-wide activity, so only within-hemisphere coupling is left.
% Track preference: z_track_diff (> 0 prefers Track 1 / L). Per session and V1 window:
%   rho      : Spearman correlation across pairs of r with z_V1 * z_HC
%   diff     : mean r of same-preference pairs (z_V1 * z_HC > 0) - opposite-preference pairs
%   null     : nPerm permutations of z_track_diff among the V1 cells and among the HPC
%              cells of the session (keeps every pair's coupling, breaks the link to tuning)
% A cell enters a window only if it fires in >= min_active events in that window.
% Options: 'sessions' (default all), 'nPerm' (1000), 'min_active' (10),
%   'cell_set' ('lt10': counts\, cells < 10 Hz (default) | 'all_cells': counts_all_rates\)
% Output: nature_neuroscience_revision_analysis\pairwise_V1_HC_coupling\gain_<gain_control>[_all_cells]\
%   coupling_POST.mat (S: per-session results incl. pairs and nulls), session_stats.csv

p = inputParser;
addRequired(p, 'gain_control', @(s) ismember(s, {'total_V1', 'own_hemisphere'}));
addParameter(p, 'sessions', 1:22, @isnumeric);
addParameter(p, 'nPerm', 1000, @isnumeric);
addParameter(p, 'min_active', 10, @isnumeric);
addParameter(p, 'cell_set', 'lt10', @(s) ismember(s, {'lt10', 'all_cells', 'rate_0p5_5'}));
addParameter(p, 'sleep', 'POST', @(s) ismember(s, {'POST', 'PRE'}));   % PRE: *_PRE folders
parse(p, gain_control, varargin{:});
P = p.Results;

rev_dir = addon_setup();
base_dir = fullfile(rev_dir, 'pairwise_V1_HC_coupling');
[counts_name, out_name, ~, rate_range] = pairwise_cell_set_folders(gain_control, P.cell_set, P.sleep);
P.rate_range = rate_range;
out_dir = fullfile(base_dir, out_name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
v1_windows = {'V1_pre', 'V1_post'};

S = struct([]);
T = table();
for nsession = P.sessions
    f = fullfile(base_dir, counts_name, sprintf('session_%02d.mat', nsession));
    if ~exist(f, 'file')
        fprintf('Session %d: no counts file, skipping\n', nsession);   % e.g. no PRE sleep
        continue
    end
    load(f, 'C');
    C.V1 = select_by_rate(C.V1, rate_range);      % single units only; MUA keeps all units
    C.HPC = select_by_rate(C.HPC, rate_range);
    R = struct('session_count', nsession, 'subject_id', C.subject_id, 'SUBJECT', C.SUBJECT, ...
        'SESSION', C.SESSION, 'n_events', numel(C.onset), 'gain_control', gain_control);
    h = C.HPC.counts.HC;
    h_MUA = C.MUA.HPC.HC;
    for w = 1:numel(v1_windows)
        wn = v1_windows{w};
        v = C.V1.counts.(wn);
        hemi_L = C.V1.region == "V1_L";
        if strcmp(gain_control, 'total_V1')
            v_MUA = repmat(C.MUA.V1_L.(wn) + C.MUA.V1_R.(wn), 1, size(v, 2));
        else
            v_MUA = NaN(size(v));
            v_MUA(:, hemi_L) = repmat(C.MUA.V1_L.(wn), 1, sum(hemi_L));
            v_MUA(:, ~hemi_L) = repmat(C.MUA.V1_R.(wn), 1, sum(~hemi_L));
        end
        v_MUA = v_MUA - v;  % leave-one-out
        ok_v = sum(v > 0, 1)' >= P.min_active;
        ok_h = sum(h > 0, 1)' >= P.min_active;
        iv = find(ok_v); ih = find(ok_h);

        r_raw = NaN(numel(iv), numel(ih));
        r_par = NaN(numel(iv), numel(ih));
        for a = 1:numel(iv)
            x = v(:, iv(a));
            for b = 1:numel(ih)
                y = h(:, ih(b));
                r_raw(a, b) = corr_safe(x, y);
                r_par(a, b) = partial_corr(x, y, [v_MUA(:, iv(a)), h_MUA - y]);
            end
        end
        zv = C.V1.z_track_diff(iv);
        zh = C.HPC.z_track_diff(ih);

        W = struct();
        W.V1_cluster_id = C.V1.cluster_id(iv); W.V1_region = C.V1.region(iv); W.z_V1 = zv;
        W.HPC_cluster_id = C.HPC.cluster_id(ih); W.HPC_region = C.HPC.region(ih); W.z_HC = zh;
        W.r_raw = r_raw; W.r_partial = r_par;
        W.mean_r_raw = mean(r_raw(:), 'omitnan'); W.mean_r_partial = mean(r_par(:), 'omitnan');
        rng(nsession * 10 + w, 'twister');
        W.raw = preference_stats(r_raw, zv, zh, P.nPerm);
        rng(nsession * 10 + w, 'twister');
        W.partial = preference_stats(r_par, zv, zh, P.nPerm);
        R.(wn) = W;

        T = [T; table(nsession, C.subject_id, string(wn), numel(C.onset), numel(iv), sum(hemi_L(iv)), ...
            numel(ih), sum(isfinite(r_par(:))), W.mean_r_raw, W.mean_r_partial, ...
            W.raw.rho, W.raw.rho_null_mean, W.raw.p_rho, W.raw.diff, W.raw.diff_null_mean, W.raw.p_diff, ...
            W.partial.rho, W.partial.rho_null_mean, W.partial.p_rho, W.partial.diff, W.partial.diff_null_mean, W.partial.p_diff, ...
            'VariableNames', {'session_id', 'animal_id', 'V1_window', 'n_events', 'n_V1', 'n_V1_L', ...
            'n_HPC', 'n_pairs', 'mean_r_raw', 'mean_r_partial', ...
            'rho_raw', 'rho_raw_null', 'p_rho_raw', 'diff_raw', 'diff_raw_null', 'p_diff_raw', ...
            'rho_partial', 'rho_partial_null', 'p_rho_partial', 'diff_partial', 'diff_partial_null', 'p_diff_partial'})]; %#ok<AGROW>
        fprintf('Session %2d %s: %d V1 x %d HPC, rho partial %.3f (null %.3f, p %.3f), diff %.4f\n', ...
            nsession, wn, numel(iv), numel(ih), W.partial.rho, W.partial.rho_null_mean, W.partial.p_rho, W.partial.diff);
    end
    S = [S R]; %#ok<AGROW>
end
save(fullfile(out_dir, 'coupling_POST.mat'), 'S', 'P', '-v7.3');
writetable(T, fullfile(out_dir, 'session_stats.csv'));
fprintf('Saved %s\n', out_dir);
end

function N = select_by_rate(N, rate_range)
% Keep the single units whose mean RUN1 rate is within rate_range (inclusive)
keep = N.rate_RUN(:) >= rate_range(1) & N.rate_RUN(:) <= rate_range(2);
if all(keep), return; end
N.cluster_id = N.cluster_id(keep); N.rate_RUN = N.rate_RUN(keep);
N.region = N.region(keep); N.z_track_diff = N.z_track_diff(keep);
for f = fieldnames(N.counts)'
    N.counts.(f{1}) = N.counts.(f{1})(:, keep);
end
end

function st = preference_stats(r, zv, zh, nPerm)
% Spearman rho of pair coupling vs z_V1 * z_HC, same - opposite preference mean, and
% their nulls from permuting preferences within each region
ok = isfinite(r);
[st.rho, st.diff] = pref_stat(r, ok, zv, zh);
st.rho_null = NaN(nPerm, 1); st.diff_null = NaN(nPerm, 1);
for k = 1:nPerm
    [st.rho_null(k), st.diff_null(k)] = pref_stat(r, ok, zv(randperm(numel(zv))), zh(randperm(numel(zh))));
end
st.rho_null_mean = mean(st.rho_null, 'omitnan');
st.diff_null_mean = mean(st.diff_null, 'omitnan');
st.p_rho = (1 + sum(st.rho_null >= st.rho)) / (1 + sum(isfinite(st.rho_null)));    % one-sided: rho > null
st.p_diff = (1 + sum(st.diff_null >= st.diff)) / (1 + sum(isfinite(st.diff_null)));
end

function [rho, d] = pref_stat(r, ok, zv, zh)
prod_z = zv(:) * zh(:)';
if sum(ok(:)) < 3
    rho = NaN; d = NaN; return
end
rho = corr(r(ok), prod_z(ok), 'type', 'Spearman');
d = mean(r(ok & prod_z > 0)) - mean(r(ok & prod_z < 0));
end

function r = corr_safe(x, y)
if std(x) == 0 || std(y) == 0
    r = NaN;
else
    r = corr(x, y);
end
end

function r = partial_corr(x, y, Z)
% Correlation of the residuals of x and y after least-squares regression on [1 Z]
% (constant covariates dropped)
Z = Z(:, std(Z, 0, 1) > 0);
X = [ones(size(x)) Z];
rx = x - X * (X \ x);
ry = y - X * (X \ y);
r = corr_safe(rx, ry);
end
