%% V1-HC attractor remapping, step 4: trial-to-trial dynamics
% Does the late-grey bias (last 1 s before onset) predict the onset bias (0-200,
% 0-100, 100-200 ms), beyond spikes, speed, licks and block (ar_trial_prediction),
% also on raw log odds; and how do HPC and V1 co-fluctuate in the early / late
% grey, and does that coupling strengthen the prediction (ar_grey_coupling)?
% Output: results\trial_prediction\, results\trial_prediction_rawLO\,
%         results\grey_coupling\
P = ar_setup();

%% Pre-onset -> onset prediction
ar_trial_prediction();
ar_trial_prediction('measure', 'LO');

%% Does track onset reset the ongoing state? (continuity across onset vs within grey)
ar_onset_reset();

%% Which cells carry it? (own-cell vs co-fluctuation decomposition, cell properties,
% decoder ablation against random drops)
ar_carryover_cells();

%% HPC-V1 grey-screen coupling
ar_grey_coupling();

%% HPC-V1 lead / lag inside the grey only (from 0.3 s after lap end to onset)
ar_grey_lag();
