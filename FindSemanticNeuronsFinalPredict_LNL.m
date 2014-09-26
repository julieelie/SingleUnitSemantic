%% This script is not used anymore, it was a tentative to estimate how linear or non-linear units were.
% Directory to save data
resultsDirectory='/auto/fdata/fet/julie/semanticProject';
addpath(resultsDirectory);
addpath(genpath('/auto/k2/share/matlab/matlab82/toolbox/stats'));  % need these - not sure why it is not automatic
rmpath('/auto/k2/share/matlab/matlab82/toolbox/stats/eml');  % but not these???


% Retrieve correlation matrices
cd /auto/k6/julie/matfile/StimCorrelationMatrices
StimCorrDir=pwd;
CorrMatfiles=dir(StimCorrDir);
NCMfiles = length(CorrMatfiles);

% Retrieve some output variables
% All the output from my spider code
load('/auto/k8/julie/SemanticAnalysis_ShuffPredict.mat');
% List_matfilepath is the cell containing the path to each ConfMat file of
% each unit
LM = length(List_matfilepath);
NbCell_estimate = LM;

% Start to loop through files and extract information

semIerror = zeros(1, LM);
semZscore = zeros(1, LM);
semPvalue = zeros(1, LM);
semMSE = zeros(1, LM);
semSlope = zeros(1, LM);

for hh = 601:LM
    % Load file containing confusion matrix
    ConfMAT = load(List_matfilepath{hh});
    
    % Find best confusion matrix
    % Max information of real matrix to get correct time window
    indMax = find(ConfMAT.mi_confusionCT == max(ConfMAT.mi_confusionCT),1);
    Conf = ConfMAT.confusionMatrix{indMax};
    VocTypeSelAll = ConfMAT.VocTypeSel;
    
    % find the sound corrrelation matrix corresponding to that unit
    [Idir, Matfile, ext] = fileparts(List_matfilepath{hh});
    found = 0;
    cmf = 0;
    while found==0
        cmf=cmf + 1;
        if strcmp(sprintf('stat%s_CorrFirstVoc_%s.mat', ConfMAT.subject, Matfile(9:13)), CorrMatfiles(cmf).name)
            CORRMAT=load(fullfile(StimCorrDir, CorrMatfiles(cmf).name));
            found=1;
        end
    end
    if cmf==NCMfiles && found==0
        fprintf(1, 'WARNING: the code could not find the correlation matrix corresponding to file\n%s\n',List_matfilepath{hh});
    end
    
    % Extract relevant parameters 
    Corr = CORRMAT.CORR;
    VocTypeSelSound = CORRMAT.VocTypeSel;
    
    % Calculate zscore values
    [semIerror(hh), semZscore(hh), semPvalue(hh), semMSE(hh), semSlope(hh)] = calcSoundCorrZvalue(Corr, VocTypeSelSound, Conf, VocTypeSelAll);
    
    if mod(hh, 10) == 0
        fprintf(1, 'Finished %d/%d neurons.  Mean Z score %.2f\n', hh, LM, mean(semZscore(1:hh)));
    end
    
end

% Save data
resultsDirectory='/auto/fdata/fet/julie/semanticProject';
save(fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict_LNL_601_end.mat'), 'List_matfilepath', 'semIerror', 'semZscore', 'semPvalue', 'semMSE', 'semSlope');
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict_LNL.mat'))

