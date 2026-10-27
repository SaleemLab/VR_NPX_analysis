function [E, TC] = ar_epoch_table(sessions)
% Epoch-level track decoding around every lap transition (lap k-1 -> lap k),
% from the event-locked 100 ms bins of ar_bin_session (edges at the event):
%   end_prev   last 1 s of lap k-1 (track still on)          lap-end bins -10..-1
%   grey_early first 1 s of grey after lap k-1               lap-end bins 0..9
%   grey_late  last 1 s of grey before lap k                 onset bins -10..-1
%   onset      0-200 ms of lap k                             onset bins 0..1
%   onset_100  0-100 ms;  onset_200  100-200 ms              onset bin 0;  1
%   onset_late 200-1000 ms                                   onset bins 2..9
%   lap_rest   rest of lap k from 1 s (whole-session grid bins)
% Track log odds per bin (LO = log P(T1)/P(T2)) come from
% results\decoding\session_XX.mat; zLO = LO z-scored over the session's grid bins
% (per region and decoder) and clipped at their 0.5 / 99.5 percentiles.
% Decoders: bayes_multinomial, bayes_poisson, pls (block-held-out lobo3), their
% _lobo5 versions, and bayes_multinomial_all, pls_all (trained on all laps;
% valid for grey epochs only).
% Output (results\data\epoch_table.mat):
%   E  : one row per session x transition x region x decoder x epoch
%   TC : event-locked zLO per transition (100 ms bins, P.align_window around the
%        lap end ('off') and the onset ('on'); t_off / t_on = bin centres), with
%        a transition table and per-bin speed and region spike counts

P = ar_setup();
if nargin < 1, sessions = 1:22; end
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
decoders = {'bayes_multinomial', 'bayes_poisson', 'pls', 'bayes_multinomial_lobo5', 'bayes_poisson_lobo5', ...
    'pls_lobo5', 'bayes_multinomial_all', 'pls_all'};
w = P.bin_width;
rel = (round(P.align_window(1) / w):round(P.align_window(2) / w) - 1);
nRel = numel(rel);
% epoch: name, kind (1 lap end, 2 onset, 0 grid), rel range, lap the bins must lie in (-1: k-1, 0: grey, 1: k)
ep = {'end_prev', 1, -10:-1, -1; 'grey_early', 1, 0:9, 0; 'grey_late', 2, -10:-1, 0; 'onset', 2, 0:1, 1; ...
    'onset_100', 2, 0, 1; 'onset_200', 2, 1, 1; 'onset_late', 2, 2:9, 1; 'lap_rest', 0, [], 1};
nPos = numel(P.pos_edges) - 1;
TC.t_on = (rel + 0.5) * w;
TC.t_off = (rel + 0.5) * w;
TC.decoders = decoders; TC.regions = P.regions;

rows = cell(0, 1);
tc_trans = cell(0, 1); tc_on = cell(0, 1); tc_off = cell(0, 1);
for s = sessions
    f = fullfile(P.results, 'decoding', sprintf('session_%02d.mat', s));
    if ~exist(f, 'file'), fprintf('Session %d missing\n', s); continue; end
    S = load(f, 'out'); out = S.out;
    T = L(s).laps; K = L(s).licks; G = out.G;
    trans = find(T.valid_transition);
    nTr = numel(trans);

    D = table();
    D.session = repmat(s, nTr, 1);
    D.animal = repmat(L(s).animal_index, nTr, 1);
    D.lap = trans;
    D.track = T.track(trans);
    D.prev_track = T.prev_track(trans);
    D.is_switch = T.is_switch(trans);
    D.block = T.block(trans);
    D.lap_in_block = T.lap_in_block(trans);
    D.block_len = T.block_len(trans);
    D.run_len_prev = T.run_len_prev(trans);
    D.grey_dur = T.grey_before(trans);
    D.is_probe = T.is_probe(trans);
    D.prev_probe = T.is_probe(trans - 1);
    D.prev_lap_dur = T.t_off(trans - 1) - T.t_on(trans - 1);
    D.t_on = T.t_on(trans);
    tc_trans{end + 1} = D; %#ok<AGROW>

    % bin index per transition and relative bin, for each event kind
    idx = cell(1, 2);
    for kd = 1:2
        M = NaN(nTr, nRel);
        rows_k = find(G.kind == kd);
        [~, it] = ismember(G.trans(rows_k), trans);
        [~, ir] = ismember(G.rel(rows_k), rel);
        ok = it > 0 & ir > 0;
        M(sub2ind([nTr nRel], it(ok), ir(ok))) = rows_k(ok);
        idx{kd} = M;
    end
    % epoch bins per transition and lick intervals
    ebins = cell(nTr, size(ep, 1)); eiv = zeros(nTr, 2, size(ep, 1));
    for i = 1:nTr
        k = trans(i);
        for e = 1:size(ep, 1)
            if ep{e, 2} == 0
                b = find(G.kind == 0 & G.lap == k & G.t_start >= T.t_on(k) + 1);
                eiv(i, :, e) = [T.t_on(k) + 1, T.t_off(k)];
            else
                b = idx{ep{e, 2}}(i, ismember(rel, ep{e, 3}));
                b = b(~isnan(b));
                need = [k - 1, 0, k]; need = need(ep{e, 4} + 2);
                b = b(G.lap(b) == need);
                ev = [T.t_off(k - 1), T.t_on(k)]; ev = ev(ep{e, 2});
                eiv(i, :, e) = ev + [min(ep{e, 3}), max(ep{e, 3}) + 1] * w;
            end
            ebins{i, e} = b;
        end
    end

    for r = 1:numel(out.region)
        R = out.region(r);
        if isempty(R.bayes_multinomial.LO), continue; end
        nspk = G.n_spikes.(R.name);
        for d = 1:numel(decoders)
            [LO, post] = get_LO(R, decoders{d});
            zLO = normalise(LO, G.kind == 0);
            tc_on{end + 1} = tc_entry(s, R.name, decoders{d}, zLO, idx{2}, G.speed, nspk); %#ok<AGROW>
            tc_off{end + 1} = tc_entry(s, R.name, decoders{d}, zLO, idx{1}, G.speed, nspk); %#ok<AGROW>
            for e = 1:size(ep, 1)
                Mx = NaN(nTr, 12);
                for i = 1:nTr
                    b = ebins{i, e};
                    if isempty(b) || all(isnan(LO(b))), continue; end
                    lo = LO(b);
                    iv = eiv(i, :, e);
                    Mx(i, 1) = numel(b);
                    Mx(i, 2) = mean(lo, 'omitnan');
                    Mx(i, 3) = mean(zLO(b), 'omitnan');
                    lo_nz = lo(isfinite(lo) & lo ~= 0);
                    Mx(i, 4) = mean(lo_nz > 0);
                    Mx(i, 5) = mean(G.speed(b), 'omitnan');
                    Mx(i, 6) = sum(nspk(b)) / numel(b);
                    Mx(i, 7) = sum(K.t >= iv(1) & K.t < iv(2) & K.side == 1);
                    Mx(i, 8) = sum(K.t >= iv(1) & K.t < iv(2) & K.side == 2);
                    if ~isempty(post) && ep{e, 4} ~= 0
                        tr = T.track(trans(i) - (ep{e, 4} == -1));
                        pb = discretize(G.pos(b), P.pos_edges);
                        pm = post((tr - 1) * nPos + (1:nPos), b);
                        Mx(i, 9) = mean(sum(pm, 1), 'omitnan');
                        [~, mp] = max(pm, [], 1);
                        ok = ~isnan(pb) & G.speed(b) > P.bayes_speed;
                        Mx(i, 10) = median(abs(mp(ok)' - pb(ok))) * 10;
                        Mx(i, 11) = mean(G.pos(b), 'omitnan');
                    end
                    Mx(i, 12) = mean((lo_nz > 0) == (T.track(trans(i)) == 1));
                end
                X = D;
                X.region = repmat(string(R.name), nTr, 1);
                X.n_cells = repmat(R.n_cells, nTr, 1);
                X.decoder = repmat(string(decoders{d}), nTr, 1);
                X.epoch = repmat(string(ep{e, 1}), nTr, 1);
                X = [X array2table(Mx, 'VariableNames', {'n_win', 'LO', 'zLO', 'frac_T1', 'speed', ...
                    'spikes_per_win', 'licks_L', 'licks_R', 'p_true_track', 'pos_err', 'pos', 'frac_upcoming'})]; %#ok<AGROW>
                rows{end + 1} = X; %#ok<AGROW>
            end
        end
    end
    fprintf('Session %d: %d transitions\n', s, nTr);
end
E = vertcat(rows{:});
TC.trans = vertcat(tc_trans{:});
TC.on = [tc_on{:}];
TC.off = [tc_off{:}];
save(fullfile(P.results, 'data', 'epoch_table.mat'), 'E', 'TC', '-v7.3');
end

function c = tc_entry(s, region, decoder, zLO, M, speed, nspk)
% event-locked matrices (transitions x bins); NaN where a bin is missing
z = NaN(size(M), 'single'); sp = z; n = z;
ok = ~isnan(M);
z(ok) = zLO(M(ok)); sp(ok) = speed(M(ok)); n(ok) = nspk(M(ok));
c = struct('session', s, 'region', region, 'decoder', decoder, 'z', z, 'speed', sp, 'spikes', n);
end

function [LO, post] = get_LO(R, name)
% decoder names: <base>[_lobo5 | _all], base = bayes_multinomial, bayes_poisson, pls
post = [];
parts = regexp(name, '^(bayes_multinomial|bayes_poisson|pls)_?(lobo5|all)?$', 'tokens', 'once');
base = parts{1};
switch parts{2}
    case ''
        LO = double(R.(base).LO);
        if startsWith(base, 'bayes'), post = R.(base).posterior; end
    case 'lobo5'
        LO = double(R.(base).LO5);
        if startsWith(base, 'bayes'), post = R.(base).posterior5; end
    case 'all'
        LO = double(R.(base).LO_all);
end
end

function z = normalise(LO, ref)
% z-score with the mean / SD of the reference bins (session grid), clip at the
% 0.5 / 99.5 percentiles of the reference bins
ok = isfinite(LO) & ref;
if ~any(ok), z = NaN(size(LO)); return; end
mu = mean(LO(ok)); sd = std(LO(ok));
z = (LO - mu) / sd;
lim = prctile(z(ok), [0.5 99.5]);
z = min(max(z, lim(1)), lim(2));
end
