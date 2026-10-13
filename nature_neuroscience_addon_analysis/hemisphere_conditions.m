function C = hemisphere_conditions(type, ripple_info)
% Hemisphere-matched event conditions for V1->HPC reactivation coherence AUC,
% as in POST (corticohippocampal_reactivation_temporal_KDE_diff_V1_detection.m):
% Track 1 (V1 bias >= th) uses hemisphere 2, Track 2 (V1 bias <= -th) uses
% hemisphere 1.
% type 'spindle' : within-session spindle amplitude percentile quartiles
% type 'SO'      : SO phase peak [-pi/2, pi/2] vs trough
% type 'SO_trough_bilateral' : matched-hemisphere SO trough with the other
%                  hemisphere in peak (unilateral) vs also in trough (bilateral)
% Returns struct array with fields label, c1 (Track 1), c2 (Track 2).

switch type
    case 'spindle'
        pct = ripple_info.spindle_percentile;
        edges = [0 25 50 75 100];
        labels = {'Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'};
        for q = 1:4
            if q == 1
                in_q = pct >= edges(q) & pct <= edges(q+1);
            else
                in_q = pct > edges(q) & pct <= edges(q+1);
            end
            C(q).label = labels{q}; %#ok<AGROW>
            C(q).c1 = in_q(:, 2);
            C(q).c2 = in_q(:, 1);
        end
    case 'SO'
        ph = ripple_info.SO_phase;
        peak = ph >= -pi/2 & ph <= pi/2;
        trough = (ph >= -pi & ph < -pi/2) | (ph > pi/2 & ph <= pi);
        C(1).label = 'SO Peak';   C(1).c1 = peak(:, 2);   C(1).c2 = peak(:, 1);
        C(2).label = 'SO Trough'; C(2).c1 = trough(:, 2); C(2).c2 = trough(:, 1);
    case 'SO_trough_bilateral'
        % Matched (dominant) hemisphere in SO trough; other hemisphere in peak
        % (unilateral trough) or also in trough (bilateral trough)
        ph = ripple_info.SO_phase;
        peak = ph >= -pi/2 & ph <= pi/2;
        trough = (ph >= -pi & ph < -pi/2) | (ph > pi/2 & ph <= pi);
        C(1).label = 'Unilateral Trough';
        C(1).c1 = trough(:, 2) & peak(:, 1);   C(1).c2 = trough(:, 1) & peak(:, 2);
        C(2).label = 'Bilateral Trough';
        C(2).c1 = trough(:, 2) & trough(:, 1); C(2).c2 = trough(:, 1) & trough(:, 2);
    otherwise
        error('Unknown condition type %s', type);
end
end
