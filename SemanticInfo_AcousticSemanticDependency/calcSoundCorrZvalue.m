function [errVal, zVal, pVal, MSE, slope] = calcSoundCorrZvalue(Corr, VocTypeSelSound, Conf, VocTypeSelAll)


% Debug fig flag - set to 1 to see results
debugFig = 0;

% Some error checking
nSounds = length(VocTypeSelSound);
nTotal = length(VocTypeSelAll);
if ( sum(strcmp(VocTypeSelSound, VocTypeSelAll(1:nSounds))) ~= nSounds)
    fprintf(1, 'ERROR: Call Type missmatch between acoustic and neural data\n');
end
    

% Get 500 random permutations by permuations rangin between 5 and 25 percent
nRand = 500;
[indRand, meanAcousCorr, meanActual, sdActual] = makeLowerCorr(Corr, VocTypeSelSound, nRand, [5 25]);

% Make cell array of index into stims
uniqueStimType = unique(VocTypeSelAll);
nStimType = length(uniqueStimType);
indStim = cell(1, nStimType);
for i=1:nStimType
    indStim{i} = find(strcmp(VocTypeSelAll, uniqueStimType(i)));   
end

% Calculate information for random matrix
mi_tot = zeros(2, nRand);
mi_tot2 = zeros(2, nRand);
mi_diag = zeros(2, nRand);
mi_error = zeros(2, nRand);
mi_diag_uni = zeros(2, nRand);
mi_all_error_uni = zeros(2, nRand);
mi_diag_uni_cat = zeros(2, nRand);
mi_real_error_uni = zeros(2, nRand);

for i = 1:nRand
    indRandFull = [indRand(i, :) nSounds+1:nTotal];
    [mi_tot(1,i), mi_tot2(1,i), mi_diag(1,i), mi_error(1,i), mi_diag_uni(1,i), ...
        mi_all_error_uni(1,i), mi_diag_uni_cat(1,i), mi_real_error_uni(1,i)] = info_matrix_perarea(Conf(indRandFull,indRandFull), indStim, 0);

end

% Calculate information for real matrix
[MI_tot, MI_tot2, MI_diag, MI_error, MI_diag_uni, ...
        MI_all_error_uni, MI_diag_uni_cat, MI_real_error_uni] = info_matrix_perarea(Conf, indStim, 0);

% Perform a linear model between the within category acoustical correlation
% and the information

mdl = LinearModel.fit(meanAcousCorr, mi_diag_uni_cat(1,:));
beta = mdl.Coefficients.Estimate;

% Calculate a prediction error
errVal = (MI_diag_uni_cat - mdl.predict(meanActual));
zVal = errVal/sqrt(mdl.MSE);
pVal = 2*(1 - normcdf(abs(zVal))); % two tail t-test
MSE = mdl.MSE;
slope = beta(2);

if (debugFig)
    figure(3);
    hold off;
    plot(meanAcousCorr, mi_diag_uni_cat(1,:), '+' );
    hold on;
    plot(meanActual, MI_diag_uni_cat, 'r+', 'MarkerSize', 14);
       
    x(1) = min(min(meanAcousCorr), meanActual);
    x(2) = max(max(meanAcousCorr), meanActual);
    
    y = beta(1) + beta(2).*x;
    plot(x,y,'k');
    
    % This also plots everything
    mdl.plot
    title(sprintf('Error %f p=%f', errVal, pVal));
    xlabel('Acoustic Correlation');
    ylabel('Category Information');
end

return






