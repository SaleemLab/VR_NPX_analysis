function assemble_V1_split_half_PSTH_POST(in_dir)
% Merge run_V1_split_half_decoding_POST outputs across sessions into a ripple
% PSTH file for the non-duplicated POST ripple events, in the same format as
% KDE_reactivation_ripples_PSTH_POST_unique.mat (events x 200 bins, -1 to 1 s from
% onset, 10 ms step). Processing as the original PSTH (extract_KDE_reactivation_
% ripples_PSTH_PRE): SWS ripples, probe 1 then probe 2, sessions in order; first
% 200 bins of each [onset-1 onset+1] window; log odds z-scored per session and
% probe on finite values; +/-Inf capped by cap_log_odds_PSTH (existing nan_mask,
% pooled across sessions); then rows ripple_info.event_id.
% HPC fields are copied unchanged from the existing unique file (HPC not re-decoded).
%   in_dir : folder with POST_session_XX.mat (default C:\Data\PLS_split_half)
% Output: nature_neuroscience_revision_analysis\
%   KDE_reactivation_ripples_PSTH_POST_unique_V1_split_half.mat, struct
%   KDE_reactivation_ripples_PSTH with HPC_logodds_ripples, HPC_z_logodds_ripples,
%   V1A_*, V1B_* (split 1 halves) and V1full_* (if decoded) logodds and
%   z_logodds, nan_mask, event_id, split_info (per session).

if nargin < 1, in_dir = 'C:\Data\PLS_split_half'; end
if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');

load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');
nan_mask_full = KDE_reactivation_ripples_PSTH.nan_mask;
clear KDE_reactivation_ripples_PSTH
old = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
old = old.KDE_reactivation_ripples_PSTH;

nbins = 200;
nSessions = max(ripples_all(1).session_count);
files = arrayfun(@(s) fullfile(in_dir, sprintf('POST_session_%02d.mat', s)), 1:nSessions, 'UniformOutput', false);
have = cellfun(@(f) exist(f, 'file') > 0, files);
fprintf('%d of %d sessions decoded\n', sum(have), nSessions);
if ~all(have), warning('Missing sessions %s are NaN', mat2str(find(~have))); end
first = load(files{find(have, 1)}, 'out');
decs = {'A', 'B'};
if isfield(first.out, 'full'), decs = [{'full'}, decs]; end

% Rows of the full SWS PSTH: probe 1 sessions in order, then probe 2
rows = struct(); split_info = struct([]);
for d = decs, for m = ["logodds", "z"], rows.(d{1}).(m) = cell(2, 1); end, end
for nprobe = 1:2
    for nsession = 1:nSessions
        in_session = find(ripples_all(nprobe).session_count == nsession);
        [~, ripple_idx_session] = intersect(in_session, find(ripples_all(nprobe).SWS_index == 1));
        if ~have(nsession)
            for d = decs, for m = ["logodds", "z"]
                rows.(d{1}).(m){nprobe} = [rows.(d{1}).(m){nprobe}; NaN(numel(ripple_idx_session), nbins)];
            end, end
            continue
        end
        S = load(files{nsession}, 'out');
        if nprobe == 1
            N = S.out.neurons;
            A = ismember(N.cluster_id, S.out.split(1).A.cluster_id); B = ismember(N.cluster_id, S.out.split(1).B.cluster_id);
            split_info(nsession).session = nsession; %#ok<AGROW>
            split_info(nsession).n = [sum(A) sum(B)];
            split_info(nsession).n_T1_pref = [sum(A & N.preferred_track == 1) sum(B & N.preferred_track == 1)];
            split_info(nsession).n_T2_pref = [sum(A & N.preferred_track == 2) sum(B & N.preferred_track == 2)];
            split_info(nsession).median_abs_SI = [median(abs(N.SI(A))) median(abs(N.SI(B)))];
            split_info(nsession).median_sleep_rate = [median(N.rate_sleep(A)) median(N.rate_sleep(B))];
            split_info(nsession).cv_hit_rate = [mean(S.out.split(1).A.dec.cv_hit_rate) mean(S.out.split(1).B.dec.cv_hit_rate)];
        end
        for d = decs
            switch d{1}
                case 'full', R = S.out.full.R(nprobe);
                otherwise, R = S.out.split(1).(d{1}).R(nprobe);
            end
            [L, Z] = event_rows(R, ripple_idx_session, nbins);
            rows.(d{1}).logodds{nprobe} = [rows.(d{1}).logodds{nprobe}; L];
            rows.(d{1}).z{nprobe} = [rows.(d{1}).z{nprobe}; Z];
        end
    end
end

PSTH = struct();
PSTH.HPC_logodds_ripples = old.HPC_logodds_ripples;
PSTH.HPC_z_logodds_ripples = old.HPC_z_logodds_ripples;
for d = decs
    for m = ["logodds", "z"]
        full = [rows.(d{1}).(m){1}; rows.(d{1}).(m){2}];
        assert(size(full, 1) == size(nan_mask_full, 1), 'Row count %d differs from existing PSTH %d', size(full, 1), size(nan_mask_full, 1));
        C = cap_log_odds_PSTH(struct('V1_logodds_ripples', full, 'nan_mask', nan_mask_full));
        if m == "logodds", f = sprintf('V1%s_logodds_ripples', d{1}); else, f = sprintf('V1%s_z_logodds_ripples', d{1}); end
        PSTH.(f) = C.V1_logodds_ripples(ripple_info.event_id, :);
    end
end
PSTH.nan_mask = old.nan_mask;
PSTH.event_id = ripple_info.event_id;
PSTH.sessions_decoded = find(have);
PSTH.split_info = split_info;
assert(isequaln(PSTH.nan_mask, nan_mask_full(ripple_info.event_id, :)), 'nan_mask mismatch with existing unique file');

% Checks: event means (-0.2-0 and 0-0.2 s) of z log odds, decoded sessions only
bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
wins = {bc > -0.2 & bc < 0, bc > 0 & bc < 0.2}; wn = {'pre', 'post'};
in_dec = ismember(ripple_info.session_id(:), find(have));
pairs = {'V1_z_logodds_ripples', 'V1full_z_logodds_ripples'; 'V1full_z_logodds_ripples', 'V1A_z_logodds_ripples'; ...
    'V1full_z_logodds_ripples', 'V1B_z_logodds_ripples'; 'V1A_z_logodds_ripples', 'V1B_z_logodds_ripples'};
src = PSTH; src.V1_z_logodds_ripples = old.V1_z_logodds_ripples;
for k = 1:size(pairs, 1)
    if ~isfield(src, pairs{k, 1}) || ~isfield(src, pairs{k, 2}), continue; end
    for w = 1:2
        x = mean(src.(pairs{k, 1})(:, wins{w}) + PSTH.nan_mask(:, wins{w}), 2, 'omitnan');
        y = mean(src.(pairs{k, 2})(:, wins{w}) + PSTH.nan_mask(:, wins{w}), 2, 'omitnan');
        ok = isfinite(x) & isfinite(y) & in_dec;
        fprintf('%s vs %s (%s, %d events): r = %.3f, sign agreement %.1f%%\n', pairs{k, 1}, pairs{k, 2}, wn{w}, ...
            sum(ok), corr(x(ok), y(ok)), 100 * mean(sign(x(ok)) == sign(y(ok))));
    end
end

KDE_reactivation_ripples_PSTH = PSTH;
out_file = fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique_V1_split_half.mat');
save(out_file, 'KDE_reactivation_ripples_PSTH', '-v7.3');
fprintf('Saved %s\n', out_file);
end

function [L, Z] = event_rows(R, ripple_idx_session, nbins)
% Raw and z log odds (z per probe on finite values) for the given events, first nbins bins
lo = log(R.T1_probability(:) ./ R.T2_probability(:))';
fin = isfinite(lo);
z = (lo - mean(lo(fin))) / std(lo(fin));
L = NaN(numel(ripple_idx_session), nbins); Z = L;
for i = 1:numel(ripple_idx_session)
    idx = find(R.event_id == ripple_idx_session(i));
    idx = idx(1:min(nbins, numel(idx)));
    L(i, 1:numel(idx)) = lo(idx);
    Z(i, 1:numel(idx)) = z(idx);
end
end
