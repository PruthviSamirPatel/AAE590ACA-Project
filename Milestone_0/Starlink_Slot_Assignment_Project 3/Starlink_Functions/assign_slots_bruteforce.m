function [assignment, C, pairEst] = assign_slots_bruteforce(sats, slots, cfg)
%ASSIGN_SLOTS_BRUTEFORCE Assign N satellites to N slots by exhaustive search.
%
% This avoids relying on matchpairs or other toolbox-specific functions.
% For N=5, there are only 120 permutations.
%
% MATLAB note:
%   Do not initialize pairEst with repmat(struct(),n,n). That creates a
%   no-field struct array, and assigning est, which has fields, causes:
%       "Subscripted assignment between dissimilar structures."
%   Initialize from the first actual estimate instead.

    n = numel(sats);
    if numel(slots) ~= n
        error('assign_slots_bruteforce requires equal number of satellites and slots.');
    end

    C = zeros(n,n);

    % Initialize the struct array using a real estimate so all fields match.
    firstEst = estimate_pmp_pair(sats(1), slots(1), cfg);
    pairEst = repmat(firstEst, n, n);

    for i = 1:n
        for j = 1:n
            if i == 1 && j == 1
                est = firstEst;
            else
                est = estimate_pmp_pair(sats(i), slots(j), cfg);
            end

            pairEst(i,j) = est;

            tf_days = est.tf/86400;
            mErr_deg = abs(rad2deg(est.finalMError));
            raanErr_deg = abs(rad2deg(est.finalRAANError));

            C(i,j) = cfg.assignment.timeWeight*tf_days + ...
                     cfg.assignment.meanAnomalyWeight_daysPerDeg*mErr_deg + ...
                     cfg.assignment.raanErrorWeight_daysPerDeg*raanErr_deg;
        end
    end

    P = perms(1:n);
    bestCost = inf;
    bestP = [];

    for k = 1:size(P,1)
        c = 0;
        for i = 1:n
            c = c + C(i, P(k,i));
        end

        if c < bestCost
            bestCost = c;
            bestP = P(k,:);
        end
    end

    assignment = repmat(struct('satIdx',[],'slotIdx',[],'cost',[]), n, 1);
    for i = 1:n
        assignment(i).satIdx = i;
        assignment(i).slotIdx = bestP(i);
        assignment(i).cost = C(i,bestP(i));
    end
end
