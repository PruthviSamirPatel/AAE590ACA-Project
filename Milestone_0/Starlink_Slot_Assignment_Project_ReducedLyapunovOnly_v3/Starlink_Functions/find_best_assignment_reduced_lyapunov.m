function [best, allCases] = find_best_assignment_reduced_lyapunov(DV, cfg)
%FIND_BEST_ASSIGNMENT_REDUCED_LYAPUNOV Exhaustive assignment search.

    n = numel(cfg.sats);
    P = perms(1:n);

    rawTotal = zeros(size(P,1),1);
    marginTotal = zeros(size(P,1),1);
    permStr = strings(size(P,1),1);

    for k = 1:size(P,1)
        c = 0;
        for i = 1:n
            c = c + DV.total_mps(i, P(k,i));
        end
        rawTotal(k) = c;
        marginTotal(k) = (1 + cfg.assignment.dvMarginFraction)*c;
        permStr(k) = sprintf('%d ', P(k,:));
        permStr(k) = strtrim(permStr(k));
    end

    [bestObj, idx] = min(marginTotal);
    bestP = P(idx,:);

    [sortedObj, order] = sort(marginTotal);
    sortedRaw = rawTotal(order);
    sortedPerm = permStr(order);

    allCases = table((1:size(P,1))', sortedPerm, sortedRaw, sortedObj, ...
        'VariableNames', {'Rank','Permutation','RawDV_mps','DV_WithMargin_mps'});

    best = struct();
    best.permutation = bestP;
    best.rawDeltaV_mps = rawTotal(idx);
    best.objectiveDeltaV_mps = bestObj;
    best.assignment = repmat(struct('satIdx',[],'slotIdx',[],'deltaV_mps',[], ...
        'deltaV_withMargin_mps',[]), n, 1);

    for i = 1:n
        j = bestP(i);
        best.assignment(i).satIdx = i;
        best.assignment(i).slotIdx = j;
        best.assignment(i).deltaV_mps = DV.total_mps(i,j);
        best.assignment(i).deltaV_withMargin_mps = ...
            (1 + cfg.assignment.dvMarginFraction)*DV.total_mps(i,j);
    end
end
