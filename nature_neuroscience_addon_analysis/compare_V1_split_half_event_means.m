function T = compare_V1_split_half_event_means(session_file)
% Per non-duplicated POST ripple event of one session: mean V1 log odds and z log
% odds over -0.2-0 s and 0-0.2 s from ripple onset, for the original decoder
% (revision-folder KDE_reactivation_ripples_PSTH_POST_unique.mat) and for the
% full / half A / half B decoders of run_V1_split_half_decoding_POST.
% New decoders are processed as the original PSTH: rows matched via ripples_all
% (SWS ripples, probe 1 then probe 2), first 200 bins of the [onset-1 onset+1]
% window, z-scored per probe on finite values, +/-Inf capped at the 0.5/99.5
% prctile of the session's finite values (original: capped across sessions).
% Window means omit NaN; nan_mask is not applied (same bins for all decoders).
% Saves <session>_event_means.csv and _event_means.png/.fig next to session_file.

rev_dir = 'D:\corticohippocampal_replay\nature_neuroscience_revision_analysis';
load(session_file, 'out');
load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
nsession = out.session_count;

% Full SWS PSTH row -> probe and ripple index within the session (as the PSTH)
probe = []; session = []; idx_in_session = [];
for nprobe = 1:2
    for s = 1:max(ripples_all(nprobe).session_count)
        in_session = find(ripples_all(nprobe).session_count == s);
        [~, ia] = intersect(in_session, find(ripples_all(nprobe).SWS_index == 1));
        probe = [probe; nprobe * ones(numel(ia), 1)]; %#ok<AGROW>
        session = [session; s * ones(numel(ia), 1)]; %#ok<AGROW>
        idx_in_session = [idx_in_session; ia(:)]; %#ok<AGROW>
    end
end
rows_u = find(session(ripple_info.event_id) == nsession);      % rows in unique PSTH
full_row = ripple_info.event_id(rows_u);
ev_probe = probe(full_row); ev_idx = idx_in_session(full_row);
assert(all(ripple_info.session_id(rows_u) == nsession | ripple_info.session_count(rows_u) == nsession), ...
    'Session label mismatch between ripple_info and ripples_all');

bc = (-1:0.01:1); bc = bc(1:end-1) + 0.005;
wins = struct('name', {'pre', 'post'}, 'bins', {bc > -0.2 & bc < 0, bc > 0 & bc < 0.2});

dec_names = {'Original', 'Full', 'HalfA', 'HalfB'};
V = struct();
V.Original.logodds = KDE_reactivation_ripples_PSTH.V1_logodds_ripples(rows_u, :);
V.Original.z = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples(rows_u, :);
src = {out.full, out.split(1).A, out.split(1).B};
for d = 1:3
    [V.(dec_names{d+1}).logodds, V.(dec_names{d+1}).z] = event_matrix(src{d}.R, ev_probe, ev_idx);
end

T = table(rows_u(:), ev_probe, ev_idx, 'VariableNames', {'unique_row', 'probe', 'ripple_idx_in_session'});
for d = 1:numel(dec_names)
    for w = 1:numel(wins)
        for m = ["logodds", "z"]
            T.(sprintf('%s_%s_%s', dec_names{d}, m, wins(w).name)) = mean(V.(dec_names{d}).(m)(:, wins(w).bins), 2, 'omitnan');
        end
    end
end
[p, f] = fileparts(session_file);
writetable(T, fullfile(p, [f '_event_means.csv']));

% Scatter: original vs full, full vs half A, full vs half B, half A vs half B
pairs = {'Original', 'Full'; 'Full', 'HalfA'; 'Full', 'HalfB'; 'HalfA', 'HalfB'};
fig = figure('Name', [f ' event means'], 'Position', [100 100 1500 720]);
tiledlayout(2, size(pairs, 1), 'TileSpacing', 'compact');
fprintf('%s: %d non-duplicated events\n', f, height(T));
for w = 1:numel(wins)
    for k = 1:size(pairs, 1)
        x = T.(sprintf('%s_logodds_%s', pairs{k, 1}, wins(w).name));
        y = T.(sprintf('%s_logodds_%s', pairs{k, 2}, wins(w).name));
        xz = T.(sprintf('%s_z_%s', pairs{k, 1}, wins(w).name));
        yz = T.(sprintf('%s_z_%s', pairs{k, 2}, wins(w).name));
        ok = isfinite(x) & isfinite(y); okz = isfinite(xz) & isfinite(yz);
        r = corr(x(ok), y(ok)); rz = corr(xz(okz), yz(okz));
        agree = mean(sign(xz(okz)) == sign(yz(okz)));
        fprintf('  %s: %s vs %s  r = %.3f (z r = %.3f), sign agreement %.1f%%\n', wins(w).name, pairs{k, 1}, pairs{k, 2}, r, rz, 100 * agree);
        nexttile; hold on;
        scatter(x(ok), y(ok), 6, 'k', 'filled', 'MarkerFaceAlpha', 0.2);
        lim = prctile([x(ok); y(ok)], [0.5 99.5]);
        plot(lim, lim, 'r--'); xlim(lim); ylim(lim); axis square;
        xlabel(sprintf('%s mean log odds', pairs{k, 1})); ylabel(sprintf('%s mean log odds', pairs{k, 2}));
        title(sprintf('%s window: r = %.3f (z %.3f), sign %.0f%%', wins(w).name, r, rz, 100 * agree), 'FontWeight', 'normal');
        set(gca, 'TickDir', 'out', 'Box', 'off');
    end
end
sgtitle(sprintf('%s: mean V1 log odds per non-duplicated ripple (pre -0.2-0 s, post 0-0.2 s)', f), 'Interpreter', 'none');
savefig(fig, fullfile(p, [f '_event_means.fig']));
exportgraphics(fig, fullfile(p, [f '_event_means.png']), 'Resolution', 200);
exportgraphics(fig, fullfile(p, [f '_event_means.pdf']), 'ContentType', 'vector');
print(fig, fullfile(p, [f '_event_means.svg']), '-dsvg', '-vector');
end

function [L, Z] = event_matrix(R, ev_probe, ev_idx)
% events x 200 raw and z log odds (z per probe, finite values; +/-Inf capped)
nbins = 200;
L = NaN(numel(ev_idx), nbins); Z = L;
for nprobe = 1:numel(R)
    lo = log(R(nprobe).T1_probability(:) ./ R(nprobe).T2_probability(:))';
    fin = isfinite(lo);
    cap = prctile(lo(fin), [0.5 99.5]);
    z = (lo - mean(lo(fin))) / std(lo(fin));
    zcap = (cap - mean(lo(fin))) / std(lo(fin));
    lo(lo == Inf) = cap(2); lo(lo == -Inf) = cap(1);
    z(z == Inf) = zcap(2); z(z == -Inf) = zcap(1);
    for i = find(ev_probe(:)' == nprobe)
        b = find(R(nprobe).event_id == ev_idx(i));
        b = b(1:min(nbins, numel(b)));
        L(i, 1:numel(b)) = lo(b);
        Z(i, 1:numel(b)) = z(b);
    end
end
end
