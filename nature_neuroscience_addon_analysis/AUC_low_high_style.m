function S = AUC_low_high_style()
% Titles, condition labels and colours of the low vs high AUC bar plots
% (plot_AUC_low_vs_high and the control analyses), one field per analysis.
% colors: row 1 = first (Low / Peak) condition, row 2 = second (High / Trough).
% SO_trough_uni_bi: dominant hemisphere in SO trough, non-dominant in peak
% (unilateral) vs trough (bilateral); SO_peak_uni_bi: dominant in peak,
% non-dominant in peak (bilateral) vs trough (unilateral).
    green = [76 187 135] / 255; purple = [140 105 185] / 255;
    pale = @(c) 0.35 * c + 0.65 * [1 1 1];   % pale version for the Low condition
    S.ripple_power = struct('title', 'Ripple power (HPC->V1)', 'labels', {{'Low', 'High'}}, ...
        'colors', [pale(green); green]);
    S.spindle = struct('title', 'Spindle (V1->HPC)', 'labels', {{'Low', 'High'}}, ...
        'colors', [pale(purple); purple]);
    S.SO_phase = struct('title', 'SO phase (V1->HPC)', 'labels', {{'Peak', 'Trough'}}, ...
        'colors', [80 175 220; 230 100 165] / 255);
    S.SO_trough_uni_bi = struct('title', {{'Dominant SO trough;', 'non-dominant SO (V1->HPC)'}}, ...
        'labels', {{'Peak', 'Trough'}}, 'colors', [245 190 215; 200 40 120] / 255);
    S.SO_peak_uni_bi = struct('title', {{'Dominant SO peak;', 'non-dominant SO (V1->HPC)'}}, ...
        'labels', {{'Peak', 'Trough'}}, 'colors', [190 225 245; 40 120 190] / 255);
end
