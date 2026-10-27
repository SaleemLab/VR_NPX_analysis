function [counts_name, out_name, cells_txt, rate_range] = pairwise_cell_set_folders(gain_control, cell_set, sleep)
% Folder names (under pairwise_V1_HC_coupling) of the pairwise coupling analysis:
%   cell_set 'lt10'       : counts\           -> gain_<gain_control>\             (cells < 10 Hz)
%   cell_set 'all_cells'  : counts_all_rates\ -> gain_<gain_control>_all_cells\   (all spatial cells)
%   cell_set 'rate_0p5_5' : counts_all_rates\ -> gain_<gain_control>_rate_0p5_5\  (spatial cells with
%                           mean RUN1 rate 0.5-5 Hz, selected in pairwise_V1_HC_coupling_POST)
%   sleep 'PRE' appends _PRE to both folders (default 'POST', no suffix)
% rate_range: [min max] mean RUN1 rate (Hz) applied when loading the counts
if nargin < 2 || isempty(cell_set), cell_set = 'lt10'; end
if nargin < 3 || isempty(sleep), sleep = 'POST'; end
rate_range = [0 Inf];
switch cell_set
    case 'lt10'
        counts_name = 'counts'; out_name = ['gain_' gain_control]; cells_txt = 'cells < 10 Hz';
    case 'all_cells'
        counts_name = 'counts_all_rates'; out_name = ['gain_' gain_control '_all_cells']; cells_txt = 'all spatial cells';
    case 'rate_0p5_5'
        counts_name = 'counts_all_rates'; out_name = ['gain_' gain_control '_rate_0p5_5'];
        cells_txt = 'spatial cells 0.5-5 Hz'; rate_range = [0.5 5];
    otherwise
        error('cell_set must be ''lt10'', ''all_cells'' or ''rate_0p5_5''');
end
switch sleep
    case 'POST'
    case 'PRE'
        counts_name = [counts_name '_PRE']; out_name = [out_name '_PRE']; cells_txt = [cells_txt ', PRE sleep'];
    otherwise
        error('sleep must be ''POST'' or ''PRE''');
end
end
