% function result = computeHist(state, N, current, P, edges)
%     tmp = randsample(N, current(state), true, P(state,:))';
%     result = histcounts(tmp, edges);
% end

function result = computeHist(state, N, currentstate, pstate, edges)
    tmp = randsample(N, currentstate, true, pstate)';
    result = histcounts(tmp, edges);
end