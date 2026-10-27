function setup = ar_pls_setup(X, Y, pos, nShuffle)
% Session-level choices of the PLS-KDE track decoder, following steps 2 and 4 of
% nature_neuroscience_addon_analysis/PLS_KDE_train.m (same seeds and models):
%   - "good" positions: 5 cm position bins where 10-fold cross-validated PLS
%     (3 components) + ridge logistic classification of track beats the
%     97.5th percentile of neuron-ID shuffled training data on both tracks
%   - KDE bandwidth (0.1:0.1:2) with the highest 10-fold cross-validated AUC of
%     T1 vs T2 bias on the good bins (smallest bandwidth not allowed)
% These two choices use all laps; the PLS projections themselves are refit on
% held-out blocks (ar_decode_session). Fields good_position (5 cm bin indices),
% n_good_bins, acc (positions x track CV accuracy), acc_shuffle_CI, bandwidth
% (NaN if fewer than 100 good bins: no PLS decoder for the region). Differences from PLS_KDE_train: per-position
% accuracy is counted per position bin (not histcounts over the data range), and
% the bandwidth AUC is the AUC itself (not the median of 1000 bootstraps).
%   X   : bins x neurons, residual z-scored counts (in-lap, speed >= 1, all laps)
%   Y   : track per bin (1/2);  pos : position bin index per bin (1..28)

numFolds = 10;
nPos = max(pos);
rng(200)
cv = cvpartition(Y, 'KFold', numFolds);
acc = NaN(numFolds, nPos, 2);
for fold = 1:numFolds
    tr = training(cv, fold); te = test(cv, fold);
    pred = cv_classify(X(tr, :), Y(tr), X(te, :));
    acc(fold, :, :) = position_accuracy(pos(te), Y(te), pred, nPos);
end
acc_sh = NaN(nShuffle, nPos);
parfor nshuffle = 1:nShuffle
    s = RandStream('mrg32k3a', 'Seed', nshuffle + 1000000);
    Xs = X(:, randperm(s, size(X, 2)));
    correct = zeros(nPos, 1); total = zeros(nPos, 1);
    for fold = 1:numFolds
        tr = training(cv, fold); te = test(cv, fold);
        pred = cv_classify(Xs(tr, :), Y(tr), X(te, :));
        correct = correct + accumarray(pos(te), pred == Y(te), [nPos 1]);
        total = total + accumarray(pos(te), 1, [nPos 1]);
    end
    acc_sh(nshuffle, :) = (correct ./ total)';
end
CI_U = prctile(acc_sh, 97.5, 1);
mean_acc = squeeze(mean(acc, 1, 'omitnan'));
good = find(mean_acc(:, 1)' > CI_U & mean_acc(:, 2)' > CI_U);

% Bandwidth by cross-validated AUC on the good bins
in = ismember(pos, good);
setup.acc = mean_acc;
setup.acc_shuffle_CI = prctile(acc_sh, [2.5 97.5], 1);
setup.good_position = good;
setup.n_good_bins = sum(in);
if sum(in) < 100 || numel(unique(Y(in))) < 2
    % no track-discriminative positions: no PLS decoder for this region
    setup.bandwidth = NaN; setup.bandwidth_list = []; setup.bandwidth_AUC = [];
    return
end
Xg = zscore(X(in, :)); Xg(~isfinite(Xg)) = 0;
Yg = Y(in);
XL = plsregress(Xg, Yg, 3);
[W, ~] = qr(XL, 0);
rng(200)
cv = cvpartition(Yg, 'KFold', numFolds);
bw_list = 0.1:0.1:2;
AUC = NaN(numel(bw_list), 1);
for i = 1:numel(bw_list)
    lo = NaN(numel(Yg), 1);
    for fold = 1:numFolds
        tr = training(cv, fold); te = test(cv, fold);
        lo(te) = ar_log_kde(Xg(Yg == 1 & tr, :) * W, Xg(te, :) * W, bw_list(i)) - ...
            ar_log_kde(Xg(Yg == 2 & tr, :) * W, Xg(te, :) * W, bw_list(i));
    end
    AUC(i) = rank_auc(lo(Yg == 1), lo(Yg == 2));
end
A = AUC; A(1) = NaN;
[~, k] = max(A);

setup.bandwidth = bw_list(k);
setup.bandwidth_list = bw_list;
setup.bandwidth_AUC = AUC;
end

function pred = cv_classify(Xtr, Ytr, Xte)
XL = plsregress(Xtr, Ytr, 3);
[W, ~] = qr(XL, 0);
mdl = fitclinear(Xtr * W, Ytr, 'Learner', 'logistic', 'Regularization', 'ridge');
pred = predict(mdl, Xte * W);
end

function a = position_accuracy(pos, Y, pred, nPos)
a = NaN(nPos, 2);
for t = 1:2
    sel = Y == t;
    a(:, t) = accumarray(pos(sel), pred(sel) == t, [nPos 1], @mean, NaN);
end
end

function a = rank_auc(x1, x2)
% P(x1 > x2) + 0.5 P(tie) (Mann-Whitney)
r = tiedrank([x1(:); x2(:)]);
n1 = numel(x1); n2 = numel(x2);
a = (sum(r(1:n1)) - n1 * (n1 + 1) / 2) / (n1 * n2);
end
