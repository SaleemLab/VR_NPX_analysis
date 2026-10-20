function dec = PLS_KDE_train(tvec, position, speed, lap_times, track_ID, spikes, varargin)
% Train the PLS + KDE track decoder on RUN (minimal version of the training part
% of population_analysis/PLS_KDE_Reactivation_decoding.m, same steps and seeds).
%   1. Bin RUN spikes (0.1 s, speed >= 1), z-score per neuron, subtract the
%      per-bin mean across neurons.
%   2. 10-fold CV of PLS (3 components) + ridge logistic classifier, and the same
%      with neuron-ID shuffled training data, to find position bins where both
%      tracks are classified above the shuffled 97.5th prctile ("good" bins).
%   3. PLS on z-scored good bins -> orthonormal projection W.
%   4. KDE bandwidth (0.1:0.1:2) chosen by 10-fold CV AUC of T1 vs T2 bias.
% Inputs as PLS_KDE_Reactivation_decoding:
%   tvec, position, speed : RUN behaviour vectors
%   lap_times             : [start end] per lap;  track_ID : track per lap
%   spikes                : [time id], ids 1..nNeurons (already restricted to
%                           the neurons to use, e.g. one half)
% Options: 'nShuffle' (1000), 'nBootBW' (1000), 'nNeurons' (max id)
% Output dec: W, X_T1_proj, X_T2_proj (training projections for KDE), bandwidth,
%   good_position, n_good_bins, cv_hit_rate, bandwidth_AUC, nNeurons, and for
%   plotting pos_bins, acc_T1, acc_T2 (folds x positions), acc_shuffle_mean,
%   acc_shuffle_CI, T1_proj_position, T2_proj_position.

p = inputParser;
addParameter(p, 'nShuffle', 1000, @isnumeric);
addParameter(p, 'nBootBW', 1000, @isnumeric);
addParameter(p, 'nNeurons', max(spikes(:, 2)), @isnumeric);
parse(p, varargin{:});
nShuffle = p.Results.nShuffle;
nBootBW = p.Results.nBootBW;
nNeurons = p.Results.nNeurons;
bin_size = 0.1;
numFolds = 10;

%% 1. RUN spike counts
t_edges = interp1(tvec, tvec, tvec(1):bin_size:tvec(end));
speed_bins = interp1(tvec, speed, tvec(1):bin_size:tvec(end));
t_edges(speed_bins < 1) = NaN;
bins = [t_edges' t_edges' + bin_size];

pos = interp1(tvec, position, t_edges);
pos(isnan(pos)) = 0;
pos = discretize(pos, 28) * 5;

T1_bins = InIntervals(t_edges', lap_times(track_ID == 1, :));
T2_bins = InIntervals(t_edges', lap_times(track_ID == 2, :));
[~, ~, ~, ~, ~, n] = ActivityTemplates(spikes, 'bins', bins); % z-scored counts
n(:, end+1:nNeurons) = 0;   % neurons silent in RUN (columns missing from max id)
n(~isfinite(n)) = 0;        % zero-variance neurons contribute nothing
n_residual = n - mean(n, 2);

X = [n_residual(T1_bins, :); n_residual(T2_bins, :)];
Y = [ones(sum(T1_bins), 1); 2 * ones(sum(T2_bins), 1)];
position_Z = [pos(T1_bins) pos(T2_bins)];

%% 2. Cross-validated track classification: real vs neuron-ID shuffled
rng(200)
cv = cvpartition(Y, 'KFold', numFolds);
correct1 = cell(1, numFolds); incorrect1 = correct1; correct2 = correct1; incorrect2 = correct1;
hit_rate = NaN(1, numFolds);
for fold = 1:numFolds
    tr = training(cv, fold); te = test(cv, fold);
    [pred, hit_rate(fold)] = cv_classify(X(tr, :), Y(tr), X(te, :), Y(te));
    [correct1{fold}, incorrect1{fold}, correct2{fold}, incorrect2{fold}] = split_positions(position_Z(te), Y(te), pred);
end

sh_correct = cell(1, nShuffle); sh_incorrect = sh_correct;
parfor nshuffle = 1:nShuffle
    s = RandStream('mrg32k3a', 'Seed', nshuffle + 1000000);
    X_shuffled = X(:, randperm(s, size(X, 2)));
    c = []; ic = [];
    for fold = 1:numFolds
        tr = training(cv, fold); te = test(cv, fold);
        pred = cv_classify(X_shuffled(tr, :), Y(tr), X(te, :), Y(te));
        [c1, i1, c2, i2] = split_positions(position_Z(te), Y(te), pred);
        c = [c c1 c2]; ic = [ic i1 i2];
    end
    sh_correct{nshuffle} = c; sh_incorrect{nshuffle} = ic;
end

% Proportion correct per position bin (histcounts with nbins, as original)
x = unique(pos);
nb = length(x);
dist1 = NaN(numFolds, nb); dist2 = dist1;
for fold = 1:numFolds
    dist1(fold, :) = histcounts(correct1{fold}, nb) ./ (histcounts(incorrect1{fold}, nb) + histcounts(correct1{fold}, nb));
    dist2(fold, :) = histcounts(correct2{fold}, nb) ./ (histcounts(incorrect2{fold}, nb) + histcounts(correct2{fold}, nb));
end
dist_sh = NaN(nShuffle, nb);
for nshuffle = 1:nShuffle
    dist_sh(nshuffle, :) = histcounts(sh_correct{nshuffle}, nb) ./ ...
        (histcounts(sh_correct{nshuffle}, nb) + histcounts(sh_incorrect{nshuffle}, nb));
end
CI_U = prctile(dist_sh, 97.5, 1);
good_position = x(mean(dist1) > CI_U & mean(dist2) > CI_U);
good_idx = ismember(position_Z, good_position);

%% 3. PLS on good bins
X = zscore(X(good_idx, :));
X(~isfinite(X)) = 0;
Y = Y(good_idx);
XL = plsregress(X, Y, 3);
[W, ~] = qr(XL, 0);

%% 4. KDE bandwidth by cross-validated AUC
rng(200)
cv = cvpartition(Y, 'KFold', numFolds);
bw_list = 0.1:0.1:2;
AUC = NaN(numel(bw_list), nBootBW);
for i = 1:numel(bw_list)
    pdf1 = []; pdf2 = []; Y_CV = [];
    for fold = 1:numFolds
        tr = training(cv, fold); te = test(cv, fold);
        pdf1 = [pdf1; mvksdensity(X(Y == 1 & tr, :) * W, X(te, :) * W, 'Bandwidth', bw_list(i))]; %#ok<AGROW>
        pdf2 = [pdf2; mvksdensity(X(Y == 2 & tr, :) * W, X(te, :) * W, 'Bandwidth', bw_list(i))]; %#ok<AGROW>
        Y_CV = [Y_CV; Y(Y == 1 & te); Y(Y == 2 & te)]; %#ok<AGROW>
    end
    RUN_bias = pdf1 ./ (pdf1 + pdf2);
    a = NaN(1, nBootBW);
    parfor nboot = 1:nBootBW
        s = RandStream('mrg32k3a', 'Seed', nboot);
        idx = datasample(s, 1:length(RUN_bias), length(RUN_bias));
        [~, ~, ~, A] = perfcurve(Y_CV(idx), RUN_bias(idx), 1, 'XVals', 0:0.05:1, 'NBoot', 1, 'BootType', 'per');
        a(nboot) = A(1);
    end
    AUC(i, :) = a;
end
med_AUC = median(AUC, 2);
[~, k] = max(med_AUC);
if k == 1 % smallest bandwidth can be too local
    med_AUC(1) = NaN;
    [~, k] = max(med_AUC);
end

dec.W = W;
dec.X_T1_proj = X(Y == 1, :) * W;
dec.X_T2_proj = X(Y == 2, :) * W;
dec.bandwidth = bw_list(k);
dec.bandwidth_AUC = median(AUC(k, :));
dec.good_position = good_position;
dec.n_good_bins = sum(good_idx);
dec.cv_hit_rate = hit_rate;
dec.nNeurons = nNeurons;
% For plotting: per-position classification accuracy (folds x positions) and
% neuron-ID shuffled accuracy, and the positions of the PLS training bins
dec.pos_bins = x;
dec.acc_T1 = dist1;
dec.acc_T2 = dist2;
dec.acc_shuffle_mean = mean(dist_sh, 1, 'omitnan');
dec.acc_shuffle_CI = prctile(dist_sh, [2.5 97.5], 1);
dec.T1_proj_position = position_Z(good_idx & [true(sum(T1_bins), 1); false(sum(T2_bins), 1)]');
dec.T2_proj_position = position_Z(good_idx & [false(sum(T1_bins), 1); true(sum(T2_bins), 1)]');
end

function [pred, hit] = cv_classify(Xtr, Ytr, Xte, Yte)
% PLS (3 comp) on training bins, ridge logistic on PLS scores, predict test bins
XL = plsregress(Xtr, Ytr, 3);
[W, ~] = qr(XL, 0);
mdl = fitclinear(Xtr * W, Ytr, 'Learner', 'logistic', 'Regularization', 'ridge');
pred = predict(mdl, Xte * W);
hit = mean(pred == Yte);
end

function [c1, i1, c2, i2] = split_positions(pos, Y, pred)
c1 = pos(Y == pred & Y == 1); c2 = pos(Y == pred & Y == 2);
i1 = pos(Y ~= pred & Y == 1); i2 = pos(Y ~= pred & Y == 2);
end
