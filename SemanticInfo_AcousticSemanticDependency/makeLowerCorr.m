function [indRand, changeSem, changeStim, meanCorrOut, meanCorrAll, sdCorrAll] = makeLowerCorr(corrMat, stimType, nRand, perChange)

% Makes an array of random index that preserve the mean correlation in
% corrMat obtained for each group in stimInd.
% indRand has size nRand x length(corrMat)
% perChange is a range [low high] of percent change for the number of
% row/columns to perturb



% Make space for output variables
indRand = zeros(nRand,length(corrMat));
meanCorrOut = zeros(1, nRand);
changeSem = zeros(1, nRand);
changeStim = zeros(1, nRand);

% Calculate the mean correlations in each group, the grand mean and the sd
% mean
uniqueStimType = unique(stimType);
nStimType = length(uniqueStimType);

meanCorr = zeros(1, nStimType);
indStim = cell(1, nStimType);

for i=1:nStimType
    indStim{i} = find(strcmp(stimType, uniqueStimType(i)));   
    meanCorr(i) = calcMeanCorrSub(corrMat, indStim{i});
end

meanCorrAll = mean(meanCorr);
sdCorrAll = std(meanCorr);

% Number to perturb
nStims = length(corrMat);


for ir=1:nRand   
    perChangeVal = perChange(1) + rand(1)*(perChange(2) - perChange(1));   
    nChange = fix(nStims*perChangeVal/100);
    changeStim(ir) = nChange/nStims;      % this is the change val as a
%     function of permutations...

% Permute nchange
    indChange = randperm(nStims, nChange);
    indRand(ir,:) = 1:nStims;
    indRand(ir,indChange) = indChange(randperm(nChange));
  % Calculate the actual number of stims that changed category 
    nChangeCat = 0;
    for is=1:nStims
        if ~strcmp(stimType(indRand(ir, is)), stimType(is))
            nChangeCat = nChangeCat + 1;
        end
    end
    changeSem(ir) = nChangeCat/nStims;    
           
    for i=1:nStimType
        meanCorr(i) = calcMeanCorrSub(corrMat(indRand(ir,:),indRand(ir,:)), indStim{i});
    end
    meanCorrOut(ir) = mean(meanCorr);
end

return



    
    

