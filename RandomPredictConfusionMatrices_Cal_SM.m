function [calfilename] = RandomPredictConfusionMatrices_Cal_SM(MatfilePath)
Bootmax=100;%Number of random matrices to compute per cell
MAT=load(MatfilePath);   % Ici MAT a plein de trucs qui ne vont pas etre les mêmes que ton mat mais c'est ici que tu charges ton fichiers de data obtenu par le code précédent
    
%% Retrieve the best Individual confusion matrix of that cell
% tu as besoin de faire cette étape (retrouver la matrice avec la fenêtre que tu as choisi) uniquement pour les matrices calculées
% à 2m puisque pour toutes les autres distances tu ne t'es concentrée que
% sur une seule fenêtre (fenêtre de calcul de la matrice de confusion).
% Ensuite tous le calcul dans ce code concerne le bootstrap pour une
% matrice de confusion par cellule, toi tu as plusieurs distances + 2
% catégories de cris (syn vs calls) donc il faut que tu fasse un loop pour
% faire tourner ce code sur toutes les matrices

Winsize=MAT.winSize;    
MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
optWin_CT=Winsize(MAXWinCT);
VocConfMat = MAT.confusionMatrix{find(Winsize==optWin_CT)}; % ici il faut que VocConfMat soit ta matrice MatConfCall
% convert the joint probability matrix to an event matrix
SF_ConfMat = VocConfMat.*MAT.neventMatrix(find(Winsize==optWin_CT)); % ici je convertis la matrice de proba jointe en matrice de contingence, tu n'as pas besoin de faire ça et peux travailler directement sur la matrice de probabilités jointes
% retrieve the vocalization type order
VocTypeSel = MAT.VocTypeSel; % ceci correspond pour toi à l'identité individuelle de chaque cri = identité individuelle de chaque cri
% retrieve the number of different categories
StimTypeCM=unique(VocTypeSel);
% Make sure stimType Background is at the end of the category list %% ceci
% n'est pas utile pour toi
BG_Ind = find(strcmp(StimTypeCM, 'BG'));
StimTypeCM = [StimTypeCM(1:(BG_Ind-1)); StimTypeCM((BG_Ind + 1):end); StimTypeCM(BG_Ind)];
NstimTypeCM=length(StimTypeCM);

%% find the vector of indices for spike trains obtained during sound (suppress those obtain during background) in the confusion matrix SF_ConfMat
%% Ceci n'est pas utile pour toi
IndBackground = find(strcmp(VocTypeSel, 'BG'));
IndVocOnly = setdiff(1:length(VocTypeSel), IndBackground);

%% Calculate bootmax random matrices with Background intact and bootmax matrices with Background shuffled
CM_IV_RandBG = cell(Bootmax, 1); % Pour mon code, je fais deux types de bootstrap pour chaque matrice. L'un où je ne shuffle que les vocalisations (variables avec BG à la fin) l'autre où je shuffle tout, même ma catégorie background (BG)
List_VocRandBG = cell(Bootmax, 1);
Nb_VocPerCatBG = cell(Bootmax, 1);
CM_IV_Rand = cell(Bootmax, 1);
List_VocRand = cell(Bootmax, 1);
Nb_VocPerCat = cell(Bootmax, 1);
for bb=1:Bootmax
    bb
    rng('shuffle'); % seeds the random number generator based on the current time
    
    % de la ligne 42 à 60, je fais le shuffle en gardant la catégorie BG
    % intacte
    VocTypeSel_randBG=[VocTypeSel(IndVocOnly(randperm(length(IndVocOnly)))); VocTypeSel(IndBackground)];
    RR=0;
    VocRandBG = zeros(1,length(VocTypeSel));
    NVPC = zeros(NstimTypeCM,1);
    for vtR=1:NstimTypeCM
        stR=StimTypeCM(vtR);
        selectorR=strcmp(VocTypeSel_randBG, stR);
        selectorIndR=find(selectorR);
        VocRandBG(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
        RR=RR+length(selectorIndR);
        NVPC(vtR) = length(selectorIndR);
    end
    confusion_matrix_vocalizationsRandBG = SF_ConfMat(:, VocRandBG);
    confusion_matrix_vocalizationsRandBG = confusion_matrix_vocalizationsRandBG./sum(sum(confusion_matrix_vocalizationsRandBG));
    CM_IV_RandBG{bb} = confusion_matrix_vocalizationsRandBG;
    List_VocRandBG{bb} = VocRandBG;
    Nb_VocPerCatBG{bb} = NVPC;
    
    % à partir de là commence le shuffle qui t'intéresse
    VocTypeSel_rand=VocTypeSel(randperm(length(VocTypeSel)));
    RR=0;
    VocRand = zeros(1,length(VocTypeSel));
    NVPC = zeros(NstimTypeCM,1);
    for vtR=1:NstimTypeCM
        stR=StimTypeCM(vtR);
        selectorR=strcmp(VocTypeSel_rand, stR);
        selectorIndR=find(selectorR);
        VocRand(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
        RR=RR+length(selectorIndR);
        NVPC(vtR) = length(selectorIndR);
    end
    confusion_matrix_vocalizationsRand = SF_ConfMat(:, VocRand);
    confusion_matrix_vocalizationsRand = confusion_matrix_vocalizationsRand./sum(sum(confusion_matrix_vocalizationsRand)); % cette ligne n'est pas nécessaire si tu travaille avec la matrice de proba jointes
    CM_IV_Rand{bb} = confusion_matrix_vocalizationsRand; % cette cellule sauve toutes les matrices random
    % calcules la matrice ConfMatBirds
    % sauve cette matrice ConfMatbirds
    % Calcul le MI de ConfMatBirds en appliquant info_matrix
    % Sauve cette valeur de MI dans un vecteur
    List_VocRand{bb} = VocRand;
    Nb_VocPerCat{bb} = NVPC;
end
RandMat.CM_IV_RandBG = CM_IV_RandBG;
RandMat.List_VocRandBG = List_VocRandBG;
RandMat.Nb_VocPerCatBG = Nb_VocPerCatBG;
RandMat.CM_IV_Rand = CM_IV_Rand;
RandMat.List_VocRand = List_VocRand;
RandMat.Nb_VocPerCat = Nb_VocPerCat;
RandMat.subject = MAT.subject;
RandMat.originalfile = MatfilePath;
[Path, Matfile] = fileparts(MatfilePath);

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',MAT.subject,['RandPMat_' Matfile(9:end) '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',['RandPMat_' Matfile(9:end) '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',MAT.subject,['RandPMat_' Matfile(9:end) '.mat']);
end

save(calfilename, '-struct', 'RandMat');
fprintf(1,'done making calculus on %s\nData save under %s\n', MatfilePath, calfilename);
end