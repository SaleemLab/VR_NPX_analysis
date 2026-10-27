function idx = MUA_match_draw(P, s)
% One random MUA-matched draw from MUA_match_prepare, with RandStream s:
% idx{1} / idx{2} = P.k(j) events drawn without replacement from each joint bin j
% of the first / second event set
    sa = cell(numel(P.k), 1); sb = sa;
    for j = 1:numel(P.k)
        sa{j} = P.ea{j}(randperm(s, numel(P.ea{j}), P.k(j)));
        sb{j} = P.eb{j}(randperm(s, numel(P.eb{j}), P.k(j)));
    end
    idx = {vertcat(sa{:}), vertcat(sb{:})};
end
