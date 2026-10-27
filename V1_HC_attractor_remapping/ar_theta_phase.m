function T = ar_theta_phase(varargin)
% Theta phase from HPC population spiking, following Vollan et al. 2025 (Nature,
% doi 10.1038/s41586-024-08527-1; population method used there when only
% hippocampal probes were available):
%   spike counts of every HPC unit in 10 ms bins, band-passed 5-10 Hz (2nd-order
%   Butterworth, zero-phase), PCA with units as variables and bins as
%   observations; phase = angle of the projection on the PC1-PC2 plane; the
%   direction is set so that phase advances in time; phase 0 = the phase of
%   minimal population firing.
% Here only the HPC spatial cells of RUN1_decoding_input.mat are available (no
% interneurons, which Vollan et al. included); PCA is fitted on running bins
% (speed > P.theta_speed) and applied to the whole session.
% Validation per session: variance in PC1-2, instantaneous frequency while
% running, population firing modulation by phase, the fraction of cells
% phase-locked (Rayleigh p < 0.001 and first-half / second-half tuning r > 0, as
% in Vollan et al. with p < 0.01), and the consistency of the phase difference
% with the Hilbert phase of the band-passed summed MUA (circ_corr_hilbert column:
% |mean exp(i * difference)|, 1 = constant offset).
% Theta cycles: consecutive phase-0 crossings; kept if 100-200 ms long (5-10 Hz)
% and the animal runs (> P.theta_speed) throughout.
% Output: results\theta_phase\session_XX.mat (t, phase, cycles), summary csv and
% figure.
p = inputParser;
addParameter(p, 'speed', 1);
addParameter(p, 'sessions', 1:22);
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'theta_phase');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
dt = 0.01; fs = 1 / dt;
[bb, aa] = butter(2, [5 10] / (fs / 2), 'bandpass');
rows = {}; phase_hist = {};
for s = O.sessions
    Ds = D(s); Tl = L(s).laps;
    cells = find(contains(Ds.region, 'HPC'));
    edges = (Ds.tvec(1):dt:Ds.tvec(end))';
    t = edges(1:end - 1) + dt / 2;
    [~, ci] = ismember(Ds.spike_id, Ds.cell_id(cells));
    ok = ci > 0;
    b = discretize(Ds.spike_times(ok), edges);
    cio = ci(ok);
    k = ~isnan(b);
    C = accumarray([b(k) cio(k)], 1, [numel(t) numel(cells)]);
    F = filtfilt(bb, aa, C);
    spd = interp1(Ds.tvec, Ds.speed, t, 'linear', 0);
    run = spd > O.speed;
    mu = mean(F(run, :), 1);
    [coef, ~, ~, ~, expl] = pca(F(run, :));
    Z = (F - mu) * coef(:, 1:2);
    z = Z(:, 1) + 1i * Z(:, 2);
    % rotation of PC1 + i PC2 (amplitude-weighted mean step): sets the direction
    % so that phase advances, and gives the PCA-plane frequency
    rr = run(1:end - 1) & run(2:end);
    zz = z(2:end) .* conj(z(1:end - 1));
    omega = angle(mean(zz(rr)));
    if omega < 0, z = conj(z); omega = -omega; end
    freq_pca = omega / (2 * pi * dt);
    % ring index: PC2 should follow PC1 a quarter cycle later
    q = max(1, round(1 / (4 * max(freq_pca, 5)) / dt));
    a = real(z); bq = imag(z);
    ring = corr(a(find(run(1:end - q))), bq(find(run(1:end - q)) + q));
    % PCA-plane phase if the plane rotates at a theta frequency, else the Hilbert
    % phase of PC1 (a weighted population signal: one-dimensional oscillation)
    use_pca = freq_pca >= 5 && freq_pca <= 10 && abs(ring) >= 0.3;
    if use_pca
        ph = angle(z);
    else
        ph = angle(hilbert(real(z)));
    end
    dph = angle(exp(1i * diff(ph)));
    % zero = phase of minimal population firing (running bins)
    pb = discretize(ph, linspace(-pi, pi, 37));
    pop = sum(C, 2);
    rate_by_phase = accumarray(pb(run), pop(run), [36 1], @mean, NaN) / dt;
    centres = linspace(-pi, pi, 37); centres = centres(1:end - 1) + pi / 36;
    [~, imin] = min(rate_by_phase);
    ph = angle(exp(1i * (ph - centres(imin))));
    % validation: frequency = mean phase advance while running
    freq = mean(dph(rr)) / (2 * pi * dt);
    rate_by_phase = circshift(rate_by_phase, -(imin - 1));
    mod_depth = (max(rate_by_phase) - min(rate_by_phase)) / (max(rate_by_phase) + min(rate_by_phase));
    % agreement with the Hilbert phase of the band-passed summed MUA: consistency
    % of the phase difference (1 = constant offset). A circular correlation
    % coefficient is not defined for uniformly distributed phases.
    hil = angle(hilbert(filtfilt(bb, aa, pop)));
    ccorr = abs(mean(exp(1i * (ph(run) - hil(run)))));
    [frac_locked, mvl] = cell_locking(C, ph, run);
    % cycles: phase-0 crossings (wrap from + to -)
    wrap = find(ph(1:end - 1) > 2 & ph(2:end) < -2) + 1;
    cs = t(wrap(1:end - 1)) - dt / 2; ce = t(wrap(2:end)) - dt / 2;
    dur = ce - cs;
    run_cycle = arrayfun(@(a, e) all(run(t >= a & t < e)), cs, ce);
    keep = dur >= 0.1 & dur <= 0.2 & run_cycle;
    cycles = [cs(keep) ce(keep)];
    phase = single(ph);
    method = "pca_plane"; if ~use_pca, method = "hilbert_pc1"; end
    save(fullfile(out_dir, sprintf('session_%02d.mat', s)), 't', 'phase', 'cycles', 'freq', 'expl', 'rate_by_phase', ...
        'method', 'ring', 'freq_pca', '-v7.3');
    rows(end + 1, :) = {s, L(s).animal_index, string(Ds.subject), numel(cells), sum(expl(1:2)), freq_pca, ring, method, ...
        freq, mod_depth, ccorr, frac_locked, median(mvl, 'omitnan'), size(cycles, 1), median(diff(cycles, 1, 2)) * 1000, ...
        sum(diff(cycles, 1, 2)) / sum(run) / dt}; %#ok<AGROW>
    phase_hist{end + 1} = rate_by_phase'; %#ok<AGROW>
    fprintf('Session %2d %s: %3d cells, PC1+2 %.1f%% var, plane %.2f Hz ring %.2f -> %s | %.2f Hz, firing mod %.2f, circ corr MUA Hilbert %.2f, %.0f%% cells locked (MVL %.2f), %d cycles (median %.0f ms, %.0f%% of running)\n', ...
        s, Ds.subject, numel(cells), sum(expl(1:2)), freq_pca, ring, method, freq, mod_depth, ccorr, 100 * frac_locked, ...
        median(mvl, 'omitnan'), size(cycles, 1), median(diff(cycles, 1, 2)) * 1000, 100 * rows{end, 16});
end
T = cell2table(rows, 'VariableNames', {'session', 'animal', 'subject', 'n_cells', 'pc12_var_pct', 'plane_freq_hz', ...
    'ring_index', 'method', 'freq_hz', 'firing_mod_depth', 'circ_corr_hilbert', 'frac_cells_locked', 'median_mvl', ...
    'n_cycles', 'cycle_ms', 'frac_running_in_cycles'});
writetable(T, fullfile(out_dir, sprintf('theta_phase_summary_speed%d.csv', O.speed)));
plot_phase(T, phase_hist, out_dir, O.speed);
end

function [frac, mvl] = cell_locking(C, ph, run)
% per cell: Rayleigh test of spike phases (running bins) and first / second
% half tuning-curve correlation
n = size(C, 2); mvl = NaN(n, 1); locked = false(n, 1);
idx = find(run); half = idx(round(end / 2));
edges = linspace(-pi, pi, 19);
for c = 1:n
    w = C(run, c); p = ph(run);
    if sum(w) < 50, continue; end
    z = sum(w .* exp(1i * p));
    R = abs(z) / sum(w); mvl(c) = R;
    nn = sum(w);
    pr = exp(sqrt(1 + 4 * nn + 4 * (nn^2 - (R * nn)^2)) - (1 + 2 * nn));   % Rayleigh p (Zar)
    r1 = run; r1(half:end) = false; r2 = run; r2(1:half - 1) = false;
    t1 = accumarray(discretize(ph(r1), edges), C(r1, c), [18 1]);
    t2 = accumarray(discretize(ph(r2), edges), C(r2, c), [18 1]);
    rr = corr(t1, t2);
    locked(c) = pr < 0.001 && rr > 0 && ~isnan(rr);
end
frac = mean(locked(~isnan(mvl)));
end

function plot_phase(T, phase_hist, out_dir, speed)
acol = lines(6);
fig = figure('Position', [50 50 1500 400], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'loose');
nexttile(tl); hold on;
x = linspace(0, 360, 37); x = x(1:end - 1) + 5;
for i = 1:numel(phase_hist)
    v = phase_hist{i}; plot([x x + 360], [v v] / mean(v), 'Color', [acol(T.animal(i), :) 0.6]);
end
xlabel('Theta phase (deg; 0 = minimal firing)'); ylabel('HPC population rate (normalised)'); xlim([0 720]);
title({'A  Population firing by MUA-PCA theta phase', 'one line per session, two cycles'}, 'FontWeight', 'normal');
nexttile(tl); hold on;
scatter(T.freq_hz, T.frac_cells_locked, 30, acol(T.animal, :), 'filled');
xlabel('Theta frequency while running (Hz)'); ylabel('Fraction of HPC cells phase-locked');
title({'B  Frequency and cell locking per session', 'Rayleigh p < 0.001, stable halves'}, 'FontWeight', 'normal');
nexttile(tl); hold on;
scatter(T.firing_mod_depth, T.circ_corr_hilbert, 30, acol(T.animal, :), 'filled');
xlabel('Population firing modulation depth'); ylabel('Phase-difference consistency with summed-MUA phase');
title({'C  Agreement with a simpler phase estimate', 'colour = animal'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, sprintf('theta_phase_speed%d', speed)));
close(fig);
end
