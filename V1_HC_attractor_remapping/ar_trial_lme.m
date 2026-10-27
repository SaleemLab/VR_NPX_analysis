function [lme, st] = ar_trial_lme(T, term, extra, response)
% LME of a transition-level response (default 'onset', the 0-200 ms zLO) on a
% predictor term plus the standard covariates (ar_trial_table columns):
%   lap_in_block, spikes and speed of the response epoch and of the late grey
%   (main effects and signed by the upcoming / previous track: log odds grow with
%   spike count in the decoded track's direction), licks (L - R) in both epochs,
%   track pair, random intercepts for animal, session and session:block.
% extra: further terms (variables or interactions 'a:b'). Rows with missing
% values in any used variable are dropped and unused pair levels removed.
% st: estimate, CI, p and n of the term.
if nargin < 3, extra = {}; end
if nargin < 4, response = 'onset'; end
sfx = ['_' response];
T.resp = T.(response);
T.spk_r = T.(['spk' sfx]); T.spd_r = T.(['spd' sfx]); T.lick_r = T.(['lick' sfx]);
T.spk_r_x = T.spk_r .* T.sgn_next; T.spd_r_x = T.spd_r .* T.sgn_next;
T.spk_gl_x = T.spk_grey_late .* T.sgn_prev; T.spd_gl_x = T.spd_grey_late .* T.sgn_prev;
extra_vars = unique(split(strjoin(string(extra), ':'), ':'))';
extra_vars = cellstr(extra_vars(extra_vars ~= ""));
vars = [{'resp', 'spk_r', 'spd_r', 'lick_r', 'spk_grey_late', 'spd_grey_late', 'lick_grey_late', 'lap_in_block'}, ...
    split(term, ':')', extra_vars];
vars = unique(vars);
T = T(all(isfinite(T{:, vars}), 2), :);
T.pair = removecats(T.pair);
f = ['resp ~ ' term ' + lap_in_block + spk_r + spk_r_x + spd_r + spd_r_x + spk_grey_late + spk_gl_x' ...
    ' + spd_grey_late + spd_gl_x + lick_r + lick_grey_late'];
if numel(categories(T.pair)) > 1, f = [f ' + pair']; end
for k = 1:numel(extra), f = [f ' + ' char(extra{k})]; end %#ok<AGROW>
re = ' + (1|sess_c) + (1|blk_c)';
if numel(unique(T.animal)) > 1, re = [' + (1|animal_c)' re]; end
lme = fitlme(T, [f re]);
st = ar_lme_coef(lme, term);
end
