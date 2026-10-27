%% V1-HC attractor remapping, step 5: theta check and report
% Theta from HPC MUA while running (> 1 and > 5 cm/s), on track and in the grey,
% against a spike-jitter null (no LFP: server down; LFP will validate later).
% Then the report figures and the self-contained HTML report.
% Output: results\theta_mua\, results\report\, results\V1_HC_remapping_report.html
P = ar_setup();

%% Theta from HPC MUA: power, then phase and cycles (Vollan et al. 2025 method,
% cycles from one minimal-firing point to the next as in Jezek et al. 2011)
ar_theta_mua();
ar_theta_phase('speed', 5);
ar_theta_phase('speed', 1);

%% Theta-cycle flicker around track onset (needs grey_coupling for the conflict model)
ar_theta_flicker('speed', 5);

%% Report
ar_report_figures();
ar_build_report();
