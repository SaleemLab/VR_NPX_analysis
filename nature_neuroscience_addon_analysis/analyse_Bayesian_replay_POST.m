function analyse_Bayesian_replay_POST(folder)
% folder: results folder inside nature_neuroscience_revision_analysis (default
% 'Bayesian_replay_decoding'; e.g. 'Bayesian_replay_decoding_rate_lt10')
% Follow-up of run_Bayesian_replay_decoding_POST (multinomial, 20 ms, POST ripples):
% 1. Weighted correlation pre vs post window, HPC and V1:
%    all scored events (both tracks averaged), HPC-significant events (post window, p <
%    0.05 on either track; HPC on its replayed track, V1 on the same and the other
%    track). HPC selected on post significance is biased towards post > pre, so the
%    mirror selection (HPC significant in the pre window) is shown alongside.
% 2. Posterior distribution during ripples: mean posterior over the 28 states (bins with
%    spikes) for pre / post windows and significant / non-significant events (each
%    region by its own post-window significance), per-bin maximum posterior and MAP
%    state histograms.
% 3. Cell firing during ripples by RUN place-field peak: number of cells, ripple-window
%    firing rate and ripple gain (ripple rate / whole-sleep rate) per peak state.
% Output: Bayesian_replay_decoding\ (wc_pre_vs_post, posterior_distribution,
% cell_ripple_firing_by_place_field) and Bayesian_replay_decoding\analysis\replay_analysis.mat

if nargin < 1, folder = 'Bayesian_replay_decoding'; end
rev = get_revision_folder();
in_dir = fullfile(rev, folder);
out_dir = fullfile(in_dir, 'analysis');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

files = dir(fullfile(in_dir, 'M*', 'session_*.mat'));
nS = numel(files);
regs = {'HPC', 'V1'}; nStates = 28;
A = struct(); EV = struct(); C = struct();
for r = 1:2
    EV.(regs{r}) = table();
    C.(regs{r}) = table();
    for f = {'post_all', 'pre_all', 'post_sig', 'post_nonsig'}
        A.(regs{r}).mean_post.(f{1}) = NaN(nS, nStates);
    end
    A.(regs{r}).maxpost = struct('post_sig', [], 'post_nonsig', [], 'pre_all', []);
    A.(regs{r}).map = struct('post_sig', zeros(1, nStates), 'post_nonsig', zeros(1, nStates));
end
subj = strings(nS, 1); sess = zeros(nS, 1);
for k = 1:nS
    S = load(fullfile(files(k).folder, files(k).name), 'out'); out = S.out;
    subj(k) = out.subject; sess(k) = out.session_index;
    fprintf('Session %d %s %s\n', out.session_index, out.subject, out.session);
    P = out.params; ob = P.onset_bin; nEv = numel(out.events.onset); nb = out.events.n_post_bins;
    H = out.region(1).post;
    [hsig, htrack] = significance(H);
    [hsig_pre, htrack_pre] = significance(out.region(1).pre);
    for r = 1:2
        R = out.region(r);
        [rsig, ~] = significance(R.post);
        % 1. weighted correlation pre vs post
        both = R.post.scored & R.pre.scored;
        Tk = table(repmat(subj(k), nEv, 1), repmat(sess(k), nEv, 1), (1:nEv)', both, hsig, htrack, hsig_pre, htrack_pre, ...
            'VariableNames', {'subject', 'session_index', 'event', 'both_scored', 'HPC_sig', 'HPC_track', 'HPC_sig_pre', 'HPC_track_pre'});
        Tk.awc_post_mean = mean(abs(R.post.wc), 2); Tk.awc_pre_mean = mean(abs(R.pre.wc), 2);
        Tk.z_post_mean = mean(R.post.z, 2); Tk.z_pre_mean = mean(R.pre.z, 2);
        same = sub2ind([nEv 2], (1:nEv)', htrack); other = sub2ind([nEv 2], (1:nEv)', 3 - htrack);
        same_pre = sub2ind([nEv 2], (1:nEv)', htrack_pre);
        Tk.awc_post_same = abs(R.post.wc(same)); Tk.awc_pre_same = abs(R.pre.wc(same));
        Tk.awc_post_other = abs(R.post.wc(other)); Tk.awc_pre_other = abs(R.pre.wc(other));
        Tk.z_post_same = R.post.z(same); Tk.z_pre_same = R.pre.z(same);
        Tk.z_post_other = R.post.z(other); Tk.z_pre_other = R.pre.z(other);
        Tk.awc_post_samepre = abs(R.post.wc(same_pre)); Tk.awc_pre_samepre = abs(R.pre.wc(same_pre));
        Tk.z_post_samepre = R.post.z(same_pre); Tk.z_pre_samepre = R.pre.z(same_pre);
        EV.(regs{r}) = [EV.(regs{r}); Tk];
        % 2. posterior distribution (bins with spikes in the window)
        win = @(e, w) ob + (w == 1) * (0:nb(e) - 1) + (w == 2) * (-nb(e):-1);
        groups = {R.post.scored, R.pre.scored, R.post.scored & rsig, R.post.scored & ~rsig};
        gname = {'post_all', 'pre_all', 'post_sig', 'post_nonsig'}; gw = [1 2 1 1];
        for g = 1:4
            acc = zeros(nStates, 1); n = 0; mx = []; mp = zeros(1, nStates);
            for e = find(groups{g})'
                Pw = double(R.posterior(:, win(e, gw(g)), e));
                Pw = Pw(:, all(isfinite(Pw), 1));
                acc = acc + sum(Pw, 2); n = n + size(Pw, 2);
                if g ~= 1
                    [m, i] = max(Pw, [], 1);
                    mx = [mx m]; %#ok<AGROW>
                    mp = mp + accumarray(i(:), 1, [nStates 1])';
                end
            end
            A.(regs{r}).mean_post.(gname{g})(k, :) = acc' / n;
            if g == 2, A.(regs{r}).maxpost.pre_all = [A.(regs{r}).maxpost.pre_all mx]; end
            if g == 3, A.(regs{r}).maxpost.post_sig = [A.(regs{r}).maxpost.post_sig mx]; A.(regs{r}).map.post_sig = A.(regs{r}).map.post_sig + mp; end
            if g == 4, A.(regs{r}).maxpost.post_nonsig = [A.(regs{r}).maxpost.post_nonsig mx]; A.(regs{r}).map.post_nonsig = A.(regs{r}).map.post_nonsig + mp; end
        end
    end
    % 3. cell firing during ripples
    info = experiment_info(out.session_index);
    is_post = strcmp(info.StimulusName, 'SleepChronic');
    options = info.session(is_post).probe(1);
    Sl = load(fullfile(options.ANALYSIS_DATAPATH, '..', 'session_clusters_Sleep.mat'), 'session_clusters');
    sc = Sl.session_clusters;
    keep = ismember(sc.cluster_id, out.cell_id);
    t = vertcat(sc.spike_times{keep}); id = vertcat(sc.spike_id{keep});
    [~, ci] = ismember(id, out.cell_id); [t, o] = sort(t); ci = ci(o);
    nC = numel(out.cell_id);
    T_sleep = t(end) - t(1);
    dur = nb * P.time_bin;
    lo = first_ge(t, out.events.onset); hi = first_ge(t, out.events.onset + dur) - 1;
    idx = cell2mat(arrayfun(@(e) (lo(e):hi(e))', (1:nEv)', 'UniformOutput', false));
    ripple_count = accumarray(ci(idx), 1, [nC 1]);
    part = zeros(nC, 1);
    for e = 1:nEv, u = unique(ci(lo(e):hi(e))); part(u) = part(u) + 1; end
    [peak_rate, peak_state] = max(out.rate_maps, [], 2);
    Ck = table(repmat(subj(k), nC, 1), repmat(sess(k), nC, 1), out.region_labels(:), peak_state, peak_rate, ...
        mean(out.rate_maps, 2), ripple_count / sum(dur), accumarray(ci, 1, [nC 1]) / T_sleep, part / nEv, ...
        'VariableNames', {'subject', 'session_index', 'region', 'peak_state', 'RUN_peak_rate', 'RUN_mean_rate', ...
        'ripple_rate', 'sleep_rate', 'participation'});
    used = true(nC, 1);
    if isfield(out, 'decoder_cell'), used = out.decoder_cell(:); end   % cells used by the decoder
    for r = 1:2
        in = contains(Ck.region, regs{r}) & used;
        C.(regs{r}) = [C.(regs{r}); Ck(in, :)];
    end
end
save(fullfile(out_dir, 'replay_analysis.mat'), 'A', 'EV', 'C', 'subj', 'sess', '-v7.3');
plot_wc_pre_post(EV, subj, sess, in_dir);
plot_posterior_distribution(A, subj, sess, in_dir);
plot_cell_firing(C, A, in_dir);
end

function [sig, track] = significance(X)
% Post-window significance on either track and the replayed track (lower p; ties: larger z)
[pm, track] = min(X.p, [], 2);
tie = X.p(:, 1) == X.p(:, 2);
[~, zt] = max(X.z, [], 2);
track(tie) = zt(tie);
sig = X.scored & pm < 0.05;
end

function plot_wc_pre_post(EV, subj, sess, out_dir)
animals = unique(subj); cmap = lines(numel(animals));
fig = figure('Name', 'wc pre vs post', 'Position', [30 30 1800 760], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
regs = {'HPC', 'V1'};
for r = 1:2
    E = EV.(regs{r});
    sets = {};
    sets{1} = struct('sel', E.both_scored, 'cols', {{'awc_pre_mean', 'awc_post_mean'}}, 'lab', {{'pre', 'post'}}, ...
        'yl', '|weighted corr|, mean of tracks', 'ttl', 'all scored events');
    sets{2} = struct('sel', E.both_scored, 'cols', {{'z_pre_mean', 'z_post_mean'}}, 'lab', {{'pre', 'post'}}, ...
        'yl', 'shuffle z, mean of tracks', 'ttl', 'all scored events');
    if r == 1
        sets{3} = struct('sel', E.both_scored & E.HPC_sig, 'cols', {{'awc_pre_same', 'awc_post_same'}}, 'lab', {{'pre', 'post'}}, ...
            'yl', '|wc|, HPC replayed track', 'ttl', 'HPC sig in POST window (selection bias: post > pre)');
        sets{4} = struct('sel', E.both_scored & E.HPC_sig_pre, 'cols', {{'awc_pre_samepre', 'awc_post_samepre'}}, 'lab', {{'pre', 'post'}}, ...
            'yl', '|wc|, track significant in pre', 'ttl', 'mirror: HPC sig in PRE window');
    else
        sets{3} = struct('sel', E.both_scored & E.HPC_sig, 'cols', {{'awc_pre_same', 'awc_post_same', 'awc_pre_other', 'awc_post_other'}}, ...
            'lab', {{'same pre', 'same post', 'other pre', 'other post'}}, 'yl', 'V1 |wc|', 'ttl', 'HPC-significant events');
        sets{4} = struct('sel', E.both_scored & E.HPC_sig, 'cols', {{'z_pre_same', 'z_post_same', 'z_pre_other', 'z_post_other'}}, ...
            'lab', {{'same pre', 'same post', 'other pre', 'other post'}}, 'yl', 'V1 shuffle z', 'ttl', 'HPC-significant events');
    end
    for q = 1:4
        s = sets{q};
        nexttile; hold on;
        Y = NaN(numel(subj), numel(s.cols));
        for k = 1:numel(subj)
            in = s.sel & E.session_index == sess(k);
            for c = 1:numel(s.cols), Y(k, c) = mean(E.(s.cols{c})(in), 'omitnan'); end
        end
        [~, ia] = ismember(subj, animals);
        for c = 1:2:numel(s.cols)
            for k = 1:numel(subj), plot([c c + 1], Y(k, [c c + 1]), '-', 'Color', [cmap(ia(k), :) 0.4]); end
        end
        for c = 1:numel(s.cols)
            scatter(c * ones(numel(subj), 1), Y(:, c), 16, cmap(ia, :), 'filled');
            plot(c + [-0.25 0.25], mean(Y(:, c), 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
        end
        ev_pre = E.(s.cols{1})(s.sel); ev_post = E.(s.cols{2})(s.sel);
        ttl = {sprintf('%s %s (n = %d events)', regs{r}, s.ttl, sum(s.sel)), ...
            sprintf('post-pre: sessions p = %.2g, events p = %.2g', signrank(Y(:, 2), Y(:, 1)), signrank(ev_post, ev_pre))};
        if numel(s.cols) == 4
            ttl{3} = sprintf('(same-other) post-pre: sessions p = %.2g', signrank((Y(:, 2) - Y(:, 1)) - (Y(:, 4) - Y(:, 3))));
        end
        title(ttl, 'FontWeight', 'normal', 'FontSize', 8, 'Interpreter', 'none');
        set(gca, 'XTick', 1:numel(s.cols), 'XTickLabel', s.lab, 'TickDir', 'out'); xlim([0.5 numel(s.cols) + 0.5]);
        ylabel(s.yl);
    end
end
[~, lab] = fileparts(out_dir);
title(tl, sprintf('%s: multinomial 20 ms weighted correlation, pre (same n bins before onset) vs post (onset to max(0.1 s, offset)); dots = session means', lab), 'Interpreter', 'none');
save_figure(fig, fullfile(out_dir, 'wc_pre_vs_post'));
close(fig);
end

function plot_posterior_distribution(A, subj, sess, out_dir)
regs = {'HPC', 'V1'};
fig = figure('Name', 'posterior distribution', 'Position', [30 30 1900 760], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
g = {'post_all', 'pre_all', 'post_sig', 'post_nonsig'};
gl = {'post, all', 'pre, all', 'post, significant', 'post, non-significant'};
gc = [0 0 0; 0.6 0.6 0.6; 0.85 0.2 0.2; 0.2 0.4 0.85];
x = 1:28;
for r = 1:2
    nexttile; hold on;
    for q = 1:4
        Y = A.(regs{r}).mean_post.(g{q});
        m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') / sqrt(sum(all(isfinite(Y), 2)));
        fill([x fliplr(x)], [m + se fliplr(m - se)], gc(q, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
        plot(x, m, '-', 'Color', gc(q, :), 'LineWidth', 1.5);
    end
    yline(1 / 28, 'k:'); xline(14.5, 'k-');
    set(gca, 'XTick', [1 7.5 14.5 21.5 28], 'XTickLabel', {'0', '70', '140|0', '70', '140'}, 'TickDir', 'out');
    ylabel('Mean posterior'); xlabel('Decoded position (cm), Track 1 | Track 2');
    legend(gl, 'Box', 'off', 'Location', 'north', 'FontSize', 7);
    title({sprintf('%s: mean posterior during ripples', regs{r}), '(session mean +/- sem; own-region significance)'}, 'FontWeight', 'normal', 'FontSize', 8);
    nexttile;
    Y = A.(regs{r}).mean_post.post_all;
    imagesc(Y * 28); colorbar; colormap(gca, parula); clim([0 3]);
    set(gca, 'YTick', 1:numel(subj), 'YTickLabel', compose('%s s%d', subj, sess), 'FontSize', 6, ...
        'XTick', [1 7.5 14.5 21.5 28], 'XTickLabel', {'0', '70', '140|0', '70', '140'});
    xline(14.5, 'w-', 'LineWidth', 1.5);
    title({sprintf('%s: post-window mean posterior x 28', regs{r}), '(1 = uniform), per session'}, 'FontWeight', 'normal', 'FontSize', 8);
    nexttile; hold on;
    edges = linspace(0, 1, 41);
    mp = A.(regs{r}).maxpost;
    histogram(mp.post_sig, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', gc(3, :), 'LineWidth', 1.5);
    histogram(mp.post_nonsig, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', gc(4, :), 'LineWidth', 1.5);
    histogram(mp.pre_all, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', gc(2, :), 'LineWidth', 1.5);
    xlabel('Max posterior in a 20 ms bin'); ylabel('Fraction of bins');
    legend({sprintf('post sig (median %.2f)', median(mp.post_sig)), sprintf('post non-sig (median %.2f)', median(mp.post_nonsig)), ...
        sprintf('pre all (median %.2f)', median(mp.pre_all))}, 'Box', 'off', 'FontSize', 7);
    title(sprintf('%s: posterior peakedness (1/28 = flat)', regs{r}), 'FontWeight', 'normal', 'FontSize', 8);
    set(gca, 'TickDir', 'out');
    nexttile; hold on;
    m = A.(regs{r}).map;
    b = bar(x, [m.post_sig / sum(m.post_sig); m.post_nonsig / sum(m.post_nonsig)]', 1, 'grouped', 'EdgeColor', 'none');
    b(1).FaceColor = gc(3, :); b(2).FaceColor = gc(4, :);
    yline(1 / 28, 'k:'); xline(14.5, 'k-');
    set(gca, 'XTick', [1 7.5 14.5 21.5 28], 'XTickLabel', {'0', '70', '140|0', '70', '140'}, 'TickDir', 'out');
    ylabel('Fraction of bins'); legend({'post sig', 'post non-sig'}, 'Box', 'off', 'FontSize', 7);
    title(sprintf('%s: MAP decoded state (post window bins)', regs{r}), 'FontWeight', 'normal', 'FontSize', 8);
end
[~, lab] = fileparts(out_dir);
title(tl, sprintf('%s: decoded posterior during POST ripples (multinomial, 20 ms; bins with spikes)', lab), 'Interpreter', 'none');
save_figure(fig, fullfile(out_dir, 'posterior_distribution'));
close(fig);
end

function plot_cell_firing(C, A, out_dir)
regs = {'HPC', 'V1'};
fig = figure('Name', 'cell firing', 'Position', [30 30 1900 760], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
x = 1:28;
xt = {'XTick', [1 7.5 14.5 21.5 28], 'XTickLabel', {'0', '70', '140|0', '70', '140'}, 'TickDir', 'out'};
for r = 1:2
    T = C.(regs{r});
    ok = T.sleep_rate > 0;
    nexttile;
    histogram(T.peak_state, 0.5:1:28.5, 'Normalization', 'probability', 'FaceColor', [0.5 0.5 0.5]);
    yline(1 / 28, 'k:'); xline(14.5, 'k-'); set(gca, xt{:});
    ylabel('Fraction of cells'); title(sprintf('%s: RUN place-field peak (n = %d cells)', regs{r}, height(T)), 'FontWeight', 'normal', 'FontSize', 8);
    vals = {T.ripple_rate, T.ripple_rate ./ T.sleep_rate, T.participation};
    yl = {'Ripple-window rate (Hz)', 'Ripple gain (ripple / whole-sleep rate)', 'Ripple participation (fraction)'};
    for q = 1:3
        nexttile; hold on;
        v = vals{q};
        m = accumarray(T.peak_state(ok), v(ok), [28 1], @(z) mean(z, 'omitnan'), NaN);
        se = accumarray(T.peak_state(ok), v(ok), [28 1], @(z) std(z, 'omitnan') / sqrt(numel(z)), NaN);
        errorbar(x, m, se, 'k.-', 'CapSize', 0);
        xline(14.5, 'k-'); set(gca, xt{:});
        mean_post = mean(A.(regs{r}).mean_post.post_all, 1, 'omitnan')';
        rr = corr(m, mean_post, 'rows', 'complete');
        ylabel(yl{q}); xlabel('RUN peak position (cm), Track 1 | Track 2');
        title({sprintf('%s by place-field peak', regs{r}), sprintf('r with decoded mean posterior = %.2f', rr)}, 'FontWeight', 'normal', 'FontSize', 8);
    end
end
[~, lab] = fileparts(out_dir);
title(tl, sprintf('%s: do cells with particular place-field positions dominate ripple firing? (decoder cells, post window; pooled over sessions, mean +/- sem)', lab), 'Interpreter', 'none');
save_figure(fig, fullfile(out_dir, 'cell_ripple_firing_by_place_field'));
close(fig);
end

function idx = first_ge(t, x)
idx = zeros(size(x));
n = numel(t);
for k = 1:numel(x)
    lo = 1; hi = n + 1;
    while lo < hi
        mid = floor((lo + hi) / 2);
        if t(mid) < x(k), lo = mid + 1; else, hi = mid; end
    end
    idx(k) = lo;
end
end

function save_figure(fig, save_path)
for attempt = 1:10
    try
        savefig(fig, [save_path '.fig']);
        exportgraphics(fig, [save_path '.png'], 'Resolution', 150);
        exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
        return
    catch err
        if attempt == 10, rethrow(err); end
        pause(3);
    end
end
end

function folder = get_revision_folder()
if exist('D:\corticohippocampal_replay', 'dir') > 0
    folder = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    folder = 'P:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
end
end
