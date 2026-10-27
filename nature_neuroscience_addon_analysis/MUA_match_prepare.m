function P = MUA_match_prepare(ia, ib, sess, X, nQ, min_events)
% Within-session MUA equalisation of event sets ia and ib (row indices of X; one
% column of X per matched MUA measure). In each session every measure is split
% into nQ quantile bins, edges from ia and ib of that session pooled, giving nQ^d
% joint bins; a joint bin with na events of ia and nb of ib keeps k = min(na, nb)
% of each, drawn at random by MUA_match_draw. The two matched sets therefore have
% the same number of events in every joint bin (same joint distribution at bin
% resolution); events are not paired one to one. Sessions with sum(k) < min_events
% are dropped.
% P.ea / P.eb : events of ia / ib in each kept joint bin
% P.k         : events drawn per set from each joint bin
% P.session   : session of each joint bin
    P.ea = {}; P.eb = {}; P.k = []; P.session = [];
    for s = unique(sess([ia; ib]))'
        a = ia(sess(ia) == s); b = ib(sess(ib) == s);
        if isempty(a) || isempty(b), continue; end
        pooled = [a; b];
        ca = joint_bin(X, pooled, a, nQ);
        cb = joint_bin(X, pooled, b, nQ);
        ea = {}; eb = {}; k = [];
        for c = unique([ca; cb])'
            n = min(sum(ca == c), sum(cb == c));
            if n > 0
                ea{end+1} = a(ca == c); eb{end+1} = b(cb == c); k(end+1) = n; %#ok<AGROW>
            end
        end
        if sum(k) >= min_events
            P.ea = [P.ea ea]; P.eb = [P.eb eb]; P.k = [P.k k];
            P.session = [P.session s * ones(1, numel(k))];
        end
    end
end

function cid = joint_bin(X, pooled, idx, nQ)
    cid = ones(numel(idx), 1);
    for j = 1:size(X, 2)
        cid = (cid - 1) * nQ + quantile_bin(X(pooled, j), X(idx, j), nQ);
    end
end

function q = quantile_bin(ref, x, nQ)
    edges = unique(prctile(ref, linspace(0, 100, nQ + 1)));
    edges(1) = -Inf; edges(end) = Inf;
    q = discretize(x, edges);
end
