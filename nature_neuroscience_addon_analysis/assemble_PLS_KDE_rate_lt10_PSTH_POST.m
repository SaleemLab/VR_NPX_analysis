function assemble_PLS_KDE_rate_lt10_PSTH_POST(in_dir)
% Merge run_PLS_KDE_rate_lt10_POST outputs (HPC and V1 decoders with neurons < 10 Hz)
% into a ripple PSTH file for the non-duplicated POST ripple events, in the same format
% as KDE_reactivation_ripples_PSTH_POST_unique.mat (events x 200 bins, -1 to 1 s from
% onset, 10 ms step). Processing as assemble_V1_split_half_PSTH_POST: SWS ripples,
% probe 1 then probe 2, sessions in order; first 200 bins of each [onset-1 onset+1]
% window; log odds z-scored per session and probe on finite values; +/-Inf capped by
% cap_log_odds_PSTH (existing nan_mask, pooled across sessions); then rows
% ripple_info.event_id.
%   in_dir : folder with POST_session_XX_<region>.mat (default
%            nature_neuroscience_revision_analysis\PLS_KDE_rate_lt10)
% Output: nature_neuroscience_revision_analysis\
%   KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat (struct
%   KDE_reactivation_ripples_PSTH: HPC_logodds_ripples, HPC_z_logodds_ripples,
%   V1_logodds_ripples, V1_z_logodds_ripples, nan_mask, event_id, neuron counts)

if exist('D:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    analysis_folder = 'P:\corticohippocampal_replay';
end
rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
if nargin < 1, in_dir = fullfile(rev_dir, 'PLS_KDE_rate_lt10'); end

load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH.mat'), 'KDE_reactivation_ripples_PSTH');
nan_mask_full = KDE_reactivation_ripples_PSTH.nan_mask;
clear KDE_reactivation_ripples_PSTH
old = load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
old = old.KDE_reactivation_ripples_PSTH;

nbins = 200;
regions = {'HPC', 'V1'};
nSessions = max(ripples_all(1).session_count);
PSTH = struct();
neurons = struct();
for r = 1:numel(regions)
    files = arrayfun(@(s) fullfile(in_dir, sprintf('POST_session_%02d_%s.mat', s, regions{r})), 1:nSessions, 'UniformOutput', false);
    have = cellfun(@(f) exist(f, 'file') > 0, files);
    fprintf('%s: %d of %d sessions decoded\n', regions{r}, sum(have), nSessions);
    assert(all(have), 'Missing %s sessions %s', regions{r}, mat2str(find(~have)));
    rows = struct('logodds', {cell(2, 1)}, 'z', {cell(2, 1)});
    for nprobe = 1:2
        for nsession = 1:nSessions
            in_session = find(ripples_all(nprobe).session_count == nsession);
            [~, ripple_idx_session] = intersect(in_session, find(ripples_all(nprobe).SWS_index == 1));
            S = load(files{nsession}, 'out');
            if nprobe == 1
                neurons.(regions{r})(nsession, :) = [numel(S.out.neurons.keep) sum(S.out.neurons.keep)];
            end
            [L, Z] = event_rows(S.out.R(nprobe), ripple_idx_session, nbins);
            rows.logodds{nprobe} = [rows.logodds{nprobe}; L];
            rows.z{nprobe} = [rows.z{nprobe}; Z];
        end
    end
    for m = ["logodds", "z"]
        full = [rows.(m){1}; rows.(m){2}];
        assert(size(full, 1) == size(nan_mask_full, 1), 'Row count %d differs from existing PSTH %d', size(full, 1), size(nan_mask_full, 1));
        C = cap_log_odds_PSTH(struct('x_logodds', full, 'nan_mask', nan_mask_full));
        if m == "logodds", f = sprintf('%s_logodds_ripples', regions{r}); else, f = sprintf('%s_z_logodds_ripples', regions{r}); end
        PSTH.(f) = C.x_logodds(ripple_info.event_id, :);
    end
end
PSTH.nan_mask = old.nan_mask;
PSTH.event_id = ripple_info.event_id;
PSTH.neurons = neurons;            % per session [spatial neurons, kept (< 10 Hz)]
assert(isequaln(PSTH.nan_mask, nan_mask_full(ripple_info.event_id, :)), 'nan_mask mismatch with existing unique file');

% Check against the original (all neurons) decoders: event means of z log odds
bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
wins = {bc > 0 & bc < 0.1, bc > -0.2 & bc < 0, bc > 0 & bc < 0.2}; wn = {'0-0.1 s', '-0.2-0 s', '0-0.2 s'};
for r = 1:numel(regions)
    f = sprintf('%s_z_logodds_ripples', regions{r});
    for w = 1:numel(wins)
        x = mean(old.(f)(:, wins{w}) + PSTH.nan_mask(:, wins{w}), 2, 'omitnan');
        y = mean(PSTH.(f)(:, wins{w}) + PSTH.nan_mask(:, wins{w}), 2, 'omitnan');
        ok = isfinite(x) & isfinite(y);
        fprintf('%s original vs < 10 Hz (%s, %d events): r = %.3f, sign agreement %.1f%%\n', regions{r}, wn{w}, sum(ok), ...
            corr(x(ok), y(ok)), 100 * mean(sign(x(ok)) == sign(y(ok))));
    end
end

KDE_reactivation_ripples_PSTH = PSTH;
out_file = fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique_rate_lt10.mat');
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
