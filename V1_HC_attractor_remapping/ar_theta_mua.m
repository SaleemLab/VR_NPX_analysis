function T = ar_theta_mua()
% Is HPC multi-unit activity (all HPC spatial cells pooled) theta-paced while
% the animal runs? No LFP needed (server down); LFP will validate it later.
% Per session: MUA in 5 ms bins; 2 s epochs taken from stretches where the
% condition holds throughout:
%   run>1 / run>5 : speed > 1 or > 5 cm/s, split into on-track and grey screen;
%                   run5_nolick: speed > 5 and no lick within 0.5 s
%   still         : speed < 1 cm/s, with or without a lick within 0.5 s (licking
%                   is 6-9 Hz in mice and could pace MUA or add artefacts)
% Power spectrum (Hann window, mean over epochs) of the real MUA and of the same
% spikes jittered uniformly by +-200 ms (keeps rate changes slower than ~2.5 Hz,
% removes theta). theta_ratio = mean power 6-10 Hz real / jittered; theta_peak =
% that ratio at 6-9 Hz over its mean at 3-5 and 11-14 Hz (a peak, not a
% broadband rise); peak_freq = frequency of the largest ratio within 5-12 Hz.
% Also: fraction of time at the start of the track (0-1 s after onset) and in the
% last 1 s of grey with speed > 1 and > 5 cm/s.
% Output: results\theta_mua\theta_mua.mat (T table, spectra), figure
P = ar_setup();
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'theta_mua');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
dt = 0.005; win = 2; nw = round(win / dt);
f = (0:nw / 2)' / win;
conds = {'run1_track', 'run1_grey', 'run5_track', 'run5_grey', 'run5_nolick', 'still_lick', 'still_nolick'};
rows = {}; SP = struct();
rng(7);
for s = 1:numel(D)
    Ds = D(s); Tl = L(s).laps;
    cells = find(contains(Ds.region, 'HPC'));
    st = Ds.spike_times(ismember(Ds.spike_id, Ds.cell_id(cells)));
    edges = (Ds.tvec(1):dt:Ds.tvec(end))';
    tc = edges(1:end - 1) + dt / 2;
    mua = histcounts(st, edges)';
    stj = st + (rand(size(st)) * 2 - 1) * 0.2;
    muaj = histcounts(stj, edges)';
    spd = interp1(Ds.tvec, Ds.speed, tc, 'linear', 0);
    inlap = false(size(tc));
    for k = 1:height(Tl), inlap(tc >= Tl.t_on(k) & tc < Tl.t_off(k)) = true; end
    % licking (6-9 Hz in mice) could pace MUA or add artefacts: bins within
    % 0.5 s of a lick
    lt = L(s).licks.t;
    near_lick = false(size(tc));
    if ~isempty(lt)
        lb = discretize(lt, edges);
        lb = lb(~isnan(lb));
        near_lick(lb) = true;
        near_lick = movmax(near_lick, 2 * round(0.5 / dt) + 1);
    end
    masks = {spd > 1 & inlap, spd > 1 & ~inlap, spd > 5 & inlap, spd > 5 & ~inlap, spd > 5 & ~near_lick, ...
        spd < 1 & near_lick, spd < 1 & ~near_lick};
    for c = 1:numel(conds)
        [pr, n_ep] = epoch_power(mua, masks{c}, nw);
        pj = epoch_power(muaj, masks{c}, nw);
        SP(s).(conds{c}) = struct('real', pr, 'jit', pj, 'n_epochs', n_ep);
        th = f >= 6 & f <= 10; band = f >= 5 & f <= 12;
        ratio = pr ./ pj;
        if n_ep >= 5
            tr = mean(pr(th)) / mean(pj(th));
            fb = f(band); [~, im] = max(ratio(band)); pf = fb(im);
            % peak prominence: a broadband rise (e.g. population bursts at rest)
            % raises the 6-10 Hz ratio without a theta peak
            tp = mean(ratio(f >= 6 & f <= 9)) / mean(ratio((f >= 3 & f <= 5) | (f >= 11 & f <= 14)));
        else
            tr = NaN; pf = NaN; tp = NaN;
        end
        rows(end + 1, :) = {s, L(s).animal_index, string(Ds.subject), string(conds{c}), numel(cells), n_ep, n_ep * win, ...
            tr, pf, tp, NaN, NaN, NaN, NaN}; %#ok<AGROW>
    end
    % running around track onset (stored on the session's last row)
    on = Tl.t_on(Tl.valid_transition);
    f_on = @(a, b, v) mean(arrayfun(@(t) mean(spd(tc >= t + a & tc < t + b) > v), on));
    run_on = [f_on(0, 1, 1), f_on(0, 1, 5), f_on(-1, 0, 1), f_on(-1, 0, 5)];
    rows(end, 11:14) = num2cell(run_on);
    r = rows(end - numel(conds) + 1:end, :);
    fprintf('Session %2d %s: %3d HPC cells | theta peak run>5 track %.2f, grey %.2f, still %.2f | peak run>5 %.1f Hz | onset 0-1 s running >1 %.2f >5 %.2f; last grey s >1 %.2f >5 %.2f\n', ...
        s, Ds.subject, numel(cells), r{3, 10}, r{4, 10}, r{7, 10}, r{3, 9}, run_on);
end
T = cell2table(rows, 'VariableNames', {'session', 'animal', 'subject', 'condition', 'n_cells', 'n_epochs', 'seconds', ...
    'theta_ratio', 'peak_freq', 'theta_peak', 'run1_onset', 'run5_onset', 'run1_late_grey', 'run5_late_grey'});
save(fullfile(out_dir, 'theta_mua.mat'), 'T', 'SP', 'f', 'conds');
writetable(T, fullfile(out_dir, 'theta_mua.csv'));
plot_theta(T, SP, f, conds, P, out_dir);
end

function [p, n] = epoch_power(x, mask, nw)
% mean power spectrum over non-overlapping epochs of nw samples where mask is
% true throughout; x mean-subtracted per epoch, Hann window
d = diff([0; mask(:); 0]);
st = find(d == 1); en = find(d == -1) - 1;
h = hann(nw);
acc = zeros(nw / 2 + 1, 1); n = 0;
for i = 1:numel(st)
    for a = st(i):nw:en(i) - nw + 1
        seg = x(a:a + nw - 1); seg = (seg - mean(seg)) .* h;
        F = fft(seg);
        acc = acc + abs(F(1:nw / 2 + 1)).^2;
        n = n + 1;
    end
end
p = acc / max(n, 1);
end

function plot_theta(T, SP, f, conds, P, out_dir)
acol = lines(6);
fig = figure('Position', [50 50 1500 420], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'loose');
nexttile(tl); hold on;
cc = [0 0 0; 0.55 0.55 0.55; 0.1 0.45 0.7; 0.45 0.7 0.9; 0.2 0.6 0.3; 0.85 0.3 0.3; 0.95 0.6 0.5];
keep = f > 1 & f <= 20;
for c = 1:numel(conds)
    R = NaN(numel(SP), sum(keep));
    for s = 1:numel(SP)
        x = SP(s).(conds{c});
        if x.n_epochs >= 5, R(s, :) = (x.real(keep) ./ x.jit(keep))'; end
    end
    m = mean(R, 1, 'omitnan'); se = std(R, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(R), 1));
    fill([f(keep); flipud(f(keep))]', [m + se fliplr(m - se)], cc(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
    plot(f(keep), m, 'Color', cc(c, :), 'LineWidth', 1.5);
end
xline(6, ':'); xline(10, ':'); yline(1, 'k:');
xlabel('Frequency (Hz)'); ylabel('MUA power, real / spike-jittered');
legend(strrep(conds, '_', ' '), 'Box', 'off');
title({'A  HPC MUA spectrum relative to jittered spikes', 'mean +- SEM over sessions'}, 'FontWeight', 'normal');
nexttile(tl); hold on;
for c = 1:numel(conds)
    x = T(T.condition == conds{c}, :);
    scatter(c + 0.1 * randn(height(x), 1), x.theta_peak, 20, acol(x.animal, :), 'filled', 'MarkerFaceAlpha', 0.7);
    plot(c + [-0.3 0.3], median(x.theta_peak, 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
end
yline(1, 'k:');
set(gca, 'XTick', 1:numel(conds), 'XTickLabel', strrep(conds, '_', ' ')); xtickangle(30);
ylabel('Theta peak (6-9 Hz over 3-5 & 11-14 Hz)');
title({'B  Theta peak prominence per session', '(colour = animal, bar = median; 1 = no peak)'}, 'FontWeight', 'normal');
nexttile(tl); hold on;
x = T(T.condition == string(conds{end}), :);
v = [x.run1_onset x.run5_onset x.run1_late_grey x.run5_late_grey];
for k = 1:4
    scatter(k + 0.1 * randn(height(x), 1), v(:, k), 20, acol(x.animal, :), 'filled', 'MarkerFaceAlpha', 0.7);
    plot(k + [-0.3 0.3], median(v(:, k), 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
end
set(gca, 'XTick', 1:4, 'XTickLabel', {'onset 0-1 s >1', 'onset 0-1 s >5', 'last grey s >1', 'last grey s >5'});
xtickangle(20); ylim([0 1]);
ylabel('Fraction of time running');
title({'C  Running around track onset', '(cm/s thresholds)'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, 'theta_mua'));
close(fig);
end
