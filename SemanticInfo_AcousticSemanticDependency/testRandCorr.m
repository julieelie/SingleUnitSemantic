function testRandCorr(Site)
% Test permuations of stimuli to form categories that have high acoustical
% similarities

% /auto/k6/julie/matfile/ConfMat/ConfMat_Site1_L1500R1500 goes with 
% statBlaBro09xxF_CorrFirstVoc_Site1

% Travaille peut-être avec GreBlu9508M ou YelBlu6903F car pour ces deux individus j'ai une liste de cellules bien typiques (linéaires vs non-linéaires).
% Par exemple: YelBlu6903F site2
% electrode L1000R900 e3 s0 ss1 fait plutôt partie des non-linéaires il me semble tandis que L1000R900 e28 s0 ss1 serait franchement linéaire (matrice de confusion diagonale)
% 
% Sinon GreBlu9508M (mais là mes deux exemples sont sur 2 sites donc matrices de corrélation de sons différents)
% Site 3 L1250R1650 e13 s0 ss1 est ma meilleure matrice diagonale
% Site2 L1100R1450 e 21 s0 ss1 est ma cellule qui ne répond pas du tout aux beggings donc plutôt non-linéaire. 
%
% trunk/src/discrimination/info_matrix_per_area

% First the correlations of the vocalizations
%load /Users/elie/Documents/MATLAB/data/matfile/StimCorrelationMatrices/statBlaBro09xxF_CorrFirstVoc_Site2.mat
%load statYelBlu6903F_CorrFirstVoc_Site2
load /Users/elie/Documents/MATLAB/data/matfile/StimCorrelationMatrices/statGreBlu9508M_CorrFirstVoc_Site2.mat


VocTypeSelSound = VocTypeSel;

% Now the neural data
% load 'Conf Data/ConfMat_Site2_L1100R1450_e21_s0_ss1'
[Pathsite, SiteCM, Siteext]=fileparts(Site);
Sitename = fullfile('/Users/elie/Documents/MATLAB/data/matfile/ConfMat',[SiteCM Siteext])
load(Sitename)

VocTypeSelAll = VocTypeSel;
clear VocTypeSel;

% Some error checking
nSounds = length(VocTypeSelSound);
nTotal = length(VocTypeSelAll);
if ( sum(strcmp(VocTypeSelSound, VocTypeSelAll(1:nSounds))) ~= nSounds)
    fprintf(1, 'ERROR: Call Type missmatch between acoustic and neural data\n');
end

% Open a pool of workers
%matlabpool open local 4
    

% Get 200 random permutations with rel high mean correlations
nRand = 500;
[indRand, changeVal, changeStim, meanAcousCorr, meanActual, sdActual] = makeLowerCorr(CORR, VocTypeSelSound, nRand, [0 100]);

% Max information of real matrix to get correct time window
indMax = find(mi_confusionCT == max(mi_confusionCT),1);

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
        mi_all_error_uni(1,i), mi_diag_uni_cat(1,i), mi_real_error_uni(1,i)] = info_matrix_perarea(confusionMatrix{indMax}(indRandFull,indRandFull), indStim, 0);

end

% Calculate information for real matrix

[MI_tot, MI_tot2, MI_diag, MI_error, MI_diag_uni, ...
        MI_all_error_uni, MI_diag_uni_cat, MI_real_error_uni] = info_matrix_perarea(confusionMatrix{indMax}, indStim, 0);
    
 % Plot histogram of the random values and actual values
 hist(mi_diag_uni_cat(1,:))
 L=line([MI_diag_uni_cat, MI_diag_uni_cat], [0,max(hist(mi_diag_uni_cat(1,:)))]);
 set(L, 'Color', [1,0,0])
 
 % Check if the MI calculated on totally random matrices with Background
 % category intact really have lower values than the random sound matrices

for i = 1:nRand
    indRandFull = [randperm(nSounds) nSounds+1:nTotal];
    [mi_tot(2,i), mi_tot2(2,i), mi_diag(2,i), mi_error(2,i), mi_diag_uni(2,i), ...
        mi_all_error_uni(2,i), mi_diag_uni_cat(2,i), mi_real_error_uni(2,i)] = info_matrix_perarea(confusionMatrix{indMax}(indRandFull,indRandFull), indStim, 0);
end

figure();
hold off;
hist(mi_diag_uni_cat',30)
legend('RandSound', 'Rand')
L=line([MI_diag_uni_cat, MI_diag_uni_cat], [0,max(max(hist(mi_diag_uni_cat')))]);
T=text(MI_diag_uni_cat, 27, 'Observed MI');
set(L, 'Color', [0,1,0]);
set(T, 'Color', [0,1,0]);
MeanRandSound = mean(mi_diag_uni_cat(1,:));
MeanRand = mean(mi_diag_uni_cat(2,:));

L2=line([MeanRandSound, MeanRandSound], [0,max(max(hist(mi_diag_uni_cat')))]);
T2=text(MeanRandSound, 27,'Mean MI RandSound');
set(L2, 'Color', [0.2,0.2,1]);
set(T2, 'Color', [0.2,0.2,1]);

L3=line([MeanRand, MeanRand], [0,max(max(hist(mi_diag_uni_cat')))]);
T3=text(MeanRand, 29,'Mean MI Rand');
set(L3, 'Color', [1,0.2,0.2]);
set(T3, 'Color', [1,0.2,0.2]);

%figure(3);
subplot(1,3,1);
hold off;
cubehelix_niceplot(1-meanAcousCorr,mi_diag_uni_cat(1,:),changeStim.*100,2)
%plot(1-meanAcousCorr, mi_diag_uni_cat(1,:), '+' );
hold on;
plot(1-meanActual, MI_diag_uni_cat, 'ro', 'MarkerSize', 14, 'MarkerFaceColor','r');

% Perform a linear model
mdlA = LinearModel.fit([1-meanAcousCorr; (1-meanAcousCorr).^2]', mi_diag_uni_cat(1,:));
x_min = min(min(1-meanAcousCorr), 1-meanActual);
x_max = max(max(1-meanAcousCorr), 1-meanActual);
beta = mdlA.Coefficients.Estimate;
x = x_min:(x_max-x_min)/100:x_max;
y = beta(1) + beta(2).*x + beta(3).*x.^2;
plot(x,y,'r');


% Calculate a prediction error
errVal = (MI_diag_uni_cat - mdlA.predict([1-meanActual; 1-meanActual.^2]'));
zval = errVal/sqrt(mdlA.MSE);
pval = 2*(1 - normcdf(abs(zval))); % two tail t-test

title(sprintf('A R2 Adj %.3f', mdlA.Rsquared.Adjusted));
xlabel('Acoustic Correlation Disruption');
ylabel('Semantic Information Disruption');

%figure(4)
subplot(1,3,2);
hold off;
xVal = changeVal;%here it was previously 1-changeVal
cubehelix_niceplot(xVal,mi_diag_uni_cat(1,:),changeStim.*100,2)
%plot(xVal, mi_diag_uni_cat(1,:), '+' );
hold on;
plot(0 , MI_diag_uni_cat, 'ro', 'MarkerSize', 14,'MarkerFaceColor','r');

% Perform a linear model
mdlS = LinearModel.fit([xVal; xVal.^2]', mi_diag_uni_cat(1,:));
x_min = 0;
x_max = 1;
x = 0:0.01:1;
beta = mdlS.Coefficients.Estimate;
y = beta(1) + beta(2).*x + beta(3).*x.^2;
plot(x,y,'r');

% This also plots everything
%mdl2.plot

% Calculate a prediction error
errVal = (MI_diag_uni_cat - mdlS.predict([1 ; 1]'));
zval = errVal/sqrt(mdlS.MSE);
pval = 2*(1 - normcdf(abs(zval))); % two tail t-test

title(sprintf('S R2 Adj %.3f', mdlS.Rsquared.Adjusted));
xlabel('Category Disruption');
ylabel('Semantic Information');


% Combined model
mdlAS = LinearModel.fit([(1-meanAcousCorr); (1-meanAcousCorr).^2; xVal; xVal.^2]', mi_diag_uni_cat(1,:));

% Model comparison A+S vs A
FvalAS_A = ((mdlA.SSE - mdlAS.SSE)/(mdlAS.NumCoefficients-mdlA.NumCoefficients))/(mdlAS.SSE/mdlAS.DFE);
pAS_A = 1 - fcdf(FvalAS_A,  mdlAS.NumCoefficients-mdlA.NumCoefficients, mdlAS.DFE);

% Model comparison A+S vs S
FvalAS_S = ((mdlS.SSE - mdlAS.SSE)/(mdlAS.NumCoefficients-mdlS.NumCoefficients))/(mdlAS.SSE/mdlAS.DFE);
pAS_S = 1 - fcdf(FvalAS_S,  mdlAS.NumCoefficients-mdlS.NumCoefficients, mdlAS.DFE);

%figure(5)
ss=subplot(1,3,3);
x_min = min(min(1-meanAcousCorr), 1-meanActual);
x_max = max(max(1-meanAcousCorr), 1-meanActual);
beta = mdlAS.Coefficients.Estimate;
x = x_min:(x_max-x_min)/100:x_max;
y=0:0.01:1;
z=beta(1) + beta(2).*x + beta(3).*x.^2 + beta(4).*y + beta(5).*y.^2;
C = [ones(numel(xVal),1) ; 2];
c = C(:);
S = repmat(40,length(xVal),1);
s=S(:);
cubehelix_niceplot_3D(1-meanAcousCorr,xVal,mi_diag_uni_cat(1,:),changeStim*100,2)
%scatter3(1-meanAcousCorr, xVal, mi_diag_uni_cat(1,:), s, '+');
hold on
%[xs,ys,zs]=sphere(30);
%sph=surf(xs*0.02+(1-meanActual), ys*0.02+0, zs*0.02+MI_diag_uni_cat,'EdgeColor','none','LineStyle','none','FaceColor','r', 'FaceLighting', 'phong','AmbientStrength', 0.5);   
%scatter3(1-meanActual, 0 , MI_diag_uni_cat, 500,'r','+');
plot3(1-meanActual,0,MI_diag_uni_cat,'ko','MarkerFaceColor',[1 0 0], 'MarkerSize',14)
hold on
plot3(x,y,z,'r', 'LineWidth',1)
set(gca,'YGrid','on','XGrid','on','ZGrid','on')
%mdlAS.plot;
title(sprintf('AS R2 Adj %.3f AS-S: F=%.2f p=%.3f AS-A: F=%.2f p=%.3f', mdlAS.Rsquared.Adjusted,FvalAS_S,pAS_S,FvalAS_A,pAS_A ));
xlabel('Acoustic Correlation Disruption')
ylabel('Category Disruption');
zlabel('Semantic Information');

% matlabpool close local


