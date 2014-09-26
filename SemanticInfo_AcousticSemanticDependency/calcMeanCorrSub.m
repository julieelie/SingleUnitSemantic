function meanVal = calcMeanCorrSub(corrMat, ind)

% Calculates the mean of the offdiagonals for the matrix corrMat but
% only for the ind.  
% Make a vector of values

nInd = length(ind);
if nInd <= 1
    meanVal = 0;
    return;
end

corrVals = zeros(1, nInd*(nInd-1));

iCorr = 1;
for i=1:nInd
    for j = 1:nInd
        if i == j
            continue;
        end
        corrVals(iCorr) = corrMat(ind(i),ind(j));
        iCorr = iCorr +1;
    end
end

meanVal = mean(corrVals);
return

