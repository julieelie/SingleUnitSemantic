%% This code computes random matrices from confusion matrices but keeping
% the background category intact as it is the case for random
% sound-constrained matrices. I run this code both on my computer and on
% strfinator on 02/06/2014
% Number of random matrices to compute per cell
% Bootmax = 100;
% %cd /auto/k6/julie/matfile
% %resultsDirectory='/auto/k8/julie';
% resultsDirectory='/Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat';
% addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
% %addpath('/auto/k1/queued');
% 
% cd /Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat
% %cd /auto/k6/julie/matfile/ConfMat
% input_dir=pwd;
% 
% Idir=pwd;
% matfiles=dir(fullfile(Idir,'ConfMat*.mat'));
% lm=length(matfiles);
% SS_Ind=zeros(lm,1);
% for ff=1:lm
%     if ~isempty(strfind(matfiles(ff).name, 'ss'))
%         SS_Ind(ff)=1;
%     end
% end
% Indices=find(SS_Ind);
% LM=length(Indices);
% 
% for hh=1:LM
%     Matfile = matfiles(Indices(hh)).name
%     MatfilePath=fullfile(Idir, Matfile);
%     MAT=load(MatfilePath);
%     
%     %% Retrieve the best Individual confusion matrix of that cell
%     Winsize=MAT.winSize;
%     MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
%     optWin_CT=Winsize(MAXWinCT);
%     VocConfMat = MAT.confusionMatrix{find(Winsize==optWin_CT)};
%     % convert the joint probability matrix to an event matrix
%     SF_ConfMat = VocConfMat.*MAT.neventMatrix(find(Winsize==optWin_CT));
%     % retrieve the vocalization type order
%     VocTypeSel = MAT.VocTypeSel;
%     % retrieve the number of different categories
%     StimTypeCM=unique(VocTypeSel);
%     % Make sure stimType Background is at the end of the category list
%     BG_Ind = find(strcmp(StimTypeCM, 'BG'));
%     StimTypeCM = [StimTypeCM(1:(BG_Ind-1)); StimTypeCM((BG_Ind + 1):end); StimTypeCM(BG_Ind)];
%     NstimTypeCM=length(StimTypeCM);
%     
% %% find the vector of indices for spike trains obtained during sound (suppress those obtain during background) in the confusion matrix SF_ConfMat
%     IndBackground = find(strcmp(VocTypeSel, 'BG'));
%     IndVocOnly = setdiff(1:length(VocTypeSel), IndBackground);
%     CM_IV_Rand = cell(Bootmax, 1);
%     List_VocRand = cell(Bootmax, 1);
%     Nb_VocPerCat = cell(Bootmax, 1);
%     for bb=1:Bootmax
%         bb
%         rng('shuffle'); % seeds the random number generator based on the current time
%         VocTypeSel_rand=[VocTypeSel(IndVocOnly(randperm(length(IndVocOnly)))); VocTypeSel(IndBackground)];
%         %if RandChoice==2 || RandChoice==0
%         RR=0;
%         VocRand = zeros(1,length(VocTypeSel));
%         NVPC = zeros(NstimTypeCM,1);
%         for vtR=1:NstimTypeCM
%             stR=StimTypeCM(vtR);
%             selectorR=strcmp(VocTypeSel_rand, stR);
%             selectorIndR=find(selectorR);
%             VocRand(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
%             RR=RR+length(selectorIndR);
%             NVPC(vtR) = length(selectorIndR);
%         end
%         confusion_matrix_vocalizationsRand = SF_ConfMat(VocRand, VocRand);
%         confusion_matrix_vocalizationsRand = confusion_matrix_vocalizationsRand./sum(sum(confusion_matrix_vocalizationsRand));
%         CM_IV_Rand{bb} = confusion_matrix_vocalizationsRand;
%         List_VocRand{bb} = VocRand;
%         Nb_VocPerCat{bb} = NVPC;
%     end
%     RandMat.CM_IV_Rand = CM_IV_Rand;
%     RandMat.List_VocRand = List_VocRand;
%     RandMat.Nb_VocPerCat = Nb_VocPerCat;
%     RandMat.subject = MAT.subject;
%     RandMat.originalfile = MatfilePath;
%     calfilename=fullfile(resultsDirectory,'RandMat',['RandMat_BGintact_' Matfile(9:end)]);
%     save(calfilename, '-struct', 'RandMat');
% end

%% Now calculate values of MI for these random matrices
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath('/auto/k1/queued');

% retrieve data confusion matrices
%cd /Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat
%resultsDirectory='/Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat';
resultsDirectory='/auto/k8/julie';
Res=load(fullfile(resultsDirectory, 'Results131130.mat'), 'List_matfilepath');
Idir=pwd;
LM=length(Res.List_matfilepath); % here I'm using the same order as the one used to calculate random matrices, spike shape and anatomy list by FindSemanticNeurons.m

% retrieve random matrices
%cd /Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat/RandMat
cd /auto/k8/julie/RandMat
Rdir=pwd;
Rmatfiles=dir(fullfile(Rdir,'RandMat_BGintact*.mat'));
LR=length(Rmatfiles);
if LR~=LM
    fprintf('Problem we don''t have the same number of random files (%d) and data files(%d)', LR, LM)
end

AvMI_tot_Rand = zeros(LR,1); % first column is for rand sound matrices, second column for rand matrices
SDMI_tot_Rand = zeros(LR,1);
AvMI_diag_uni_Rand=zeros(LM,1);
SDMI_diag_uni_Rand=zeros(LM,1);
AVMI_diag_uni_cat_Rand=zeros(LM,1);
SDMI_diag_uni_cat_Rand=zeros(LM,1);
UnitNames = cell(LM,1);

ii=0;
for hh=1:LM
    hh
   [MATpath, Matfile]=fileparts(Res.List_matfilepath{hh});
    for rr = 1:LR
        RMatfile = Rmatfiles(rr).name;
        if strcmp(RMatfile(18:(end-4)), Matfile(9:end))
            break
        end
    end
    RMatfilePath=fullfile(Rdir, RMatfile);
    %load the random matrices file
    RMAT = load(RMatfilePath);
    if ~strcmp(RMAT.subject, 'WhiBlu5396M')% No randSOund matrices for WhiBlu5396M
        ii=ii+1;
        % keep track of the cell name in the outputs
        UnitNames{ii}=Matfile(9:end);
        %% Now calculate values of MI in different areas of random and random sound matrices
        LRM=length(RMAT.CM_IV_Rand);%number of random matrices
        MI_tot_rand=zeros(LRM,1);
        MI_diag_uni_rand=zeros(LRM,1);
        MI_diag_uni_cat_rand=zeros(LRM,1);
        for rm =1:LRM
            rm
            Rmat = RMAT.CM_IV_Rand{rm};
            % construct the cell array of indices for each random category
            List_indices= RMAT.List_VocRand{rm};
            Nb_VocPerCat = RMAT.Nb_VocPerCat{rm};
            NstimTypeRM = length(Nb_VocPerCat);
            cat_rand = cell(NstimTypeRM,1);
            nni = 0;
            for cc = 1:NstimTypeRM
                    cat_rand{cc} = (nni+1):(Nb_VocPerCat(cc)+nni);
                    nni = nni+Nb_VocPerCat(cc);
            end
            [ mi_tot, mi_tot2, mi_diag, mi_error, mi_diag_uni, mi_all_error_uni, mi_diag_uni_cat, mi_real_error_uni]=info_matrix_perarea(Rmat, cat_rand);
            MI_tot_rand(rm)=mi_tot;
            MI_diag_uni_rand(rm)=mi_diag_uni;
            MI_diag_uni_cat_rand(rm)=mi_diag_uni_cat; 
        end
        AvMI_tot_Rand(ii) = mean(MI_tot_rand);
        SDMI_tot_Rand(ii) = std(MI_tot_rand);
        AvMI_diag_uni_Rand(ii)=mean(MI_diag_uni_rand);
        SDMI_diag_uni_Rand(ii)=std(MI_diag_uni_rand);
        AVMI_diag_uni_cat_Rand(ii)=mean(MI_diag_uni_cat_rand);
        SDMI_diag_uni_cat_Rand(ii)=std(MI_diag_uni_cat_rand);
    end                  
end
MI.AvMI_tot_Rand = AvMI_tot_Rand(1:ii, :);
MI.SDMI_tot_Rand = SDMI_tot_Rand(1:ii, :);
MI.AvMI_diag_uni_Rand=AvMI_diag_uni_Rand(1:ii, :);
MI.SDMI_diag_uni_Rand=SDMI_diag_uni_Rand(1:ii, :);
MI.AVMI_diag_uni_cat_Rand=AVMI_diag_uni_cat_Rand(1:ii, :);
MI.SDMI_diag_uni_cat_Rand=SDMI_diag_uni_cat_Rand(1:ii, :);
MI.UnitNames=UnitNames(1:ii);
MIfilename='MI_analysis_RandBG_140207.mat';
save(fullfile(resultsDirectory, MIfilename), '-struct', 'MI');
% February 06 2014. I'm running the MI analysis of the Random matrices...
% with background category intact to test if the index of selectivity is
% higher for these matrices compare to totally random matrices. I'm running
% this both on my computer and on strfinator.
% February 07 2014. I found a bugg in my calculation of cat_rand so I'm
% running the MI calculations again on strfinator