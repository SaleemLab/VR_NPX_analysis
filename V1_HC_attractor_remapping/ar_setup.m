function P = ar_setup()
% Paths and shared parameters of the V1-HC attractor remapping project
% (track representation around the grey-screen inter-lap interval, block
% expectation). Adds this code folder to the path.
% Data:  D:\corticohippocampal_replay\V1_HC_attractor_remapping
%   inputs\   copies of the RUN1 spikes, behaviour and decoding outputs from the
%             revision folder (never edited)
%   results\  everything written by this project

addpath(genpath(fileparts(mfilename('fullpath'))));
if exist('D:\corticohippocampal_replay', 'dir') > 0
    P.root = 'D:\corticohippocampal_replay\V1_HC_attractor_remapping';
elseif exist('P:\corticohippocampal_replay', 'dir') > 0
    P.root = 'P:\corticohippocampal_replay\V1_HC_attractor_remapping';
end
P.inputs = fullfile(P.root, 'inputs');
P.results = fullfile(P.root, 'results');

P.subjects = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
P.V1_exclude_subjects = {'M24018'};   % V1_L has 3-15 cells: excluded from main V1 analyses
P.regions = {'HPC', 'V1', 'V1_L', 'V1_R'};
P.min_cells = 5;                      % regions with fewer cells are not decoded

% Decoding bins: 100 ms, non-overlapping; a grid over the whole RUN1 session plus
% bins locked to every lap end and track onset (edges at the event) over
% P.align_window s around it (ar_bin_session.m)
P.bin_width = 0.1;
P.align_window = [-1.5 1.5];
% Bayesian: 10 cm bins (14 per track, 28 joint states), running samples > 5 cm/s
P.pos_edges = 0:10:140;
P.bayes_speed = 5;
P.rate_floor = 0.1;
% PLS-KDE: 5 cm position bins, bins with speed >= 1 cm/s (as the paper decoder)
P.pls_speed = 1;
P.pls_pos_edges = 0:5:140;
P.pls_nShuffle = 500;

% Grey screen: a transition (lap k-1 -> lap k) is used if the grey interval is
% within these limits (s); shorter gaps are logging artefacts
P.grey_limits = [3.5 8];
end
