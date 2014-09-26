%% First calculate the index of selectivity on the matrix diagonal of the...
...semantic group matrix for all the units. This code loops through all...
    ... confusion matrices to calculate MI confusion of CT matrices calculated...
    ...from individual call matrices without the diagonal. It also investigates...
    ...selectivity and invariance by looking at tghe proportion of information...
    ...in the diagonal vs outside of the diagonal of the individual call matrices
    ...This code is also going to calculate MI conf of rand matrices calculated...
    ...from individual call matrices without diagonals.
addpath(genpath('Users/elie/Documents/MATLAB/tlab/trunk/src'))
cd /Volumes/FREDERIC/ConfMat
%cd /Users/elie/Documents/MATLAB/data
%cd /auto/k6/julie/matfile/ConfMat
ConfMatPath = pwd;

cd /auto/k8/julie/RandMat
randdir=pwd;
randmatfiles = dir(fullfile(randdir, '*.mat'));
LRM = length(randmatfiles);

cd /auto/k6/julie/matfile/StimCorrelationMatrices
Ldir=pwd;
CorrMatfiles=dir(Ldir);
NCMfiles = length(CorrMatfiles);

cd /auto/k8/julie
input_dir=pwd;
load('Results131130.mat');
LM = length(List_matfilepath); % number of units
IS = zeros(LM,10);
DI = zeros(LM,9);
DI_tot = zeros(LM,9);
DI_error = zeros(LM,9);
II = zeros(LM,9);
II_uni = zeros(LM,9);
II_uni_DI = zeros(LM,9);
VT = cell(LM,1);
VS = zeros(LM,9);
NV_cat = zeros(LM,9); % number of stims in each call category for each unit
MIBCT_OptCT_woDiag = zeros(LM,1);
MIBCT_OptCT_woDiag2 = zeros(LM,1);
MIB_OptCT_woDiag = zeros(LM,1);

% Vectors that will containes values of MI for rand matrices calculated
% without diag in the individual call matrix
MeanMI_rand_wodiag = zeros(LM,1);
SDMI_rand_wodiag = zeros(LM,1);
MeanMI_randSound_wodiag = zeros(LM,1);
SDMI_randSound_wodiag = zeros(LM,1);

for cm=1:LM
    [MATpath, MATname]=fileparts(List_matfilepath{cm});
    MAT=load(fullfile(ConfMatPath, MATname));
    %MAT=load(List_matfilepath{cm});
    optWin_CT(cm);
    Winsize=MAT.winSize;
    Iwin=find(Winsize==optWin_CT(cm));
    MATRIX=MAT.confusionMatrix{Iwin};
    MATRIXwoDiag=MATRIX;
    n=size(MATRIX,1);
    MATRIXwoDiag(1:n+1:n*n) = 0;
    
    % calculating the corresponding Call type matrix with the call type in
    % the same order
    UVT=unique(MAT.VocTypeSel);
    BG = find(strcmp(UVT, 'BG'));
    UVT = [UVT(1:BG-1) ; UVT(BG+1:end) ; UVT(BG)];%keep background category at the end as for the original matrix
    if cm>1
        testVT = sum(strcmp(UVT2,UVT));
        if testVT~=length(UVT)
            printf('Problem in the order of the vocalization type in the output check VT');
        end
    end
    UVT2=UVT;
    NUVT=length(UVT2);
    confusion_matrix_CallType = zeros(NUVT, NUVT);
    confusion_matrix_CallType_woDiag = zeros(NUVT, NUVT);
    for vtR=1:NUVT
        stR=UVT2(vtR);
        selectorR=strcmp(MAT.VocTypeSel, stR);
        NV_cat(cm,vtR) = sum(selectorR);
        selectorIndR=find(selectorR);
        for vtC = 1:NUVT
            stC=UVT2(vtC);
            selectorC=strcmp(MAT.VocTypeSel, stC);
            selectorIndC=find(selectorC);
            confusion_matrix_CallType(vtR,vtC)=sum(sum(MATRIX(selectorIndR, selectorIndC)));
            confusion_matrix_CallType_woDiag(vtR,vtC)=sum(sum(MATRIXwoDiag(selectorIndR, selectorIndC)));
        end
    end
    Nevent = MAT.neventMatrix(Iwin);
    confusion_matrix_CallType_woDiag2 = confusion_matrix_CallType_woDiag .* Nevent;
    Nevent_local = sum(sum(confusion_matrix_CallType_woDiag2));
    confusion_matrix_CallType_woDiag2 = confusion_matrix_CallType_woDiag2./Nevent_local;
    MIBCT_OptCT_woDiag2(cm) = info_matrix(confusion_matrix_CallType_woDiag2);
    MIBCT_OptCT_woDiag(cm) = info_matrix(confusion_matrix_CallType_woDiag);
    
    MATRIXwoDiag2 = MATRIXwoDiag .* Nevent;
    Nevent_local2 = sum(sum(MATRIXwoDiag2));
    MATRIXwoDiag2 =  MATRIXwoDiag2 ./ Nevent_local2;
    MIB_OptCT_woDiag(cm) = info_matrix(MATRIXwoDiag2);
    
    ProbaCat = sum(confusion_matrix_CallType, 2);
    repProbaCat = repmat(ProbaCat, 1, size(confusion_matrix_CallType,2));
    confusion_matrix_CallType_cond = confusion_matrix_CallType ./ repProbaCat;
    
    % isolating the diagonal vectors
    Vect_class_cond = diag(confusion_matrix_CallType_cond,0);
    DI(cm, :) = Vect_class_cond(1:end-1); %only keep the vocalization categories
    Vect_class_tot = diag(confusion_matrix_CallType,0);
    DI_tot(cm, :) = Vect_class_tot(1:end-1);
    Vect_class_error = diag(confusion_matrix_CallType_woDiag,0);
    DI_error(cm, :) = Vect_class_error(1:end-1);
    II(cm, :) = (2*(DI_error(cm, :)) ./ DI_tot(cm, :) - 1) .* DI(cm, :);%here multiply by the percentage of classification correct from the semantic matrix with conditional probability.
    II_uni(cm, :) = DI_error(cm, :)./(DI_tot(cm, :).*(1-1./NV_cat(cm, 1:(end-1))));
    II_uni_DI(cm, :) = II_uni(cm, :) .* DI(cm, :);
    
    % calculating IS index of selectivity on the matrix diagonal
    for vt=1:NUVT
        Vect_class_temp = Vect_class_cond;
        Vect_class_temp(vt)=[];
        IS(cm,vt)=log2(Vect_class_cond(vt)/mean(Vect_class_temp));
    end
    
    %Keeping track of call type
    VT{cm}=UVT2;
    
    % Calculating the number of categories for which p>0.1*ll
    for ll=1:9
        VS(cm,ll)=sum(Vect_class_cond(1:end-1)>0.1*ll);
    end
    
%     %% Find the random matrices corresponding to that matfile and calculate the MI confusion of the Rand CT matrices after removing the diagonal in the IV Matrix
%     % find the corresponding randfile
%     for lrm=1:LRM
%         RandMatF = randmatfiles(lrm).name;
%         if strcmp(MATname(9:end), RandMatF(9:end))
%             RandMatfile = RandMatF;
%         end
%     end
%     % find the corresponding correlation sound matrix
%     for cmf=1:NCMfiles
%         if strcmp(sprintf('stat%s_CorrFirstVoc_%s.mat', MAT.subject, MATname(9:13)), CorrMatfiles(cmf).name)
%             CORR=load(fullfile(Ldir, CorrMatfiles(cmf).name));
%         end
%     end
%     
%      % retrieve the number of vocaliztaions in each vocalization category
%     Nvocs = CORR.NVT;
%     Nvocs = [Nvocs n-sum(Nvocs)]; % with n the size of the Individual vocalization matrix
%     
%     % loop through Rand file matrices
%     RandMAT = load(RandMatfile);
%     CM_IV_Rand = RandMAT.CM_IV_Rand;
%     CM_IV_RandSound = RandMAT.CM_IV_RandSound;
%     Nrm = length(CM_IV_Rand);
%     Nrms = length(CM_IV_RandSound);
%     if Nrms ~= Nrm
%         printf('not the same number of random and random sound constrained matrices FIX IT!!!');
%     end
%     
%     MI_CM_Rand_local = zeros(Nrm,1);
%     MI_CM_RandS_local = zeros(Nrm,1);
%     
%     for rm = 1:Nrm
%         % multiply by number of events to retrieve a contingent table
%         CMIV_rand = CM_IV_Rand{rm} .* Nevent;
%         CMIV_rands = CM_IV_RandSound{rm} .* Nevent;
%         
%         % calculate Matrices without the diagonal
%         MATRand_woDiag=CMIV_rand;
%         nrand=size(CMIV_rand,1);
%         if nrand~=n
%             printf('problem!! The random matrix has a different size than the original matrix of individual vocalizations/nThis is also a problem for the vector Nvocs/n')
%         end
%         MATRand_woDiag(1:n+1:n*n) = 0;
%         MATRandS_woDiag=CMIV_rands;
%         nrands=size(CMIV_rands,1);
%         if nrand~=n
%             printf('problem!! The random Sound matrix has a different size than the original matrix of individual vocalizations/nThis is also a problem for the vector Nvocs/n')
%         end
%         MATRandS_woDiag(1:n+1:n*n) = 0;
%         
%         % Compile Matrices
%         NUVT_here = length(Nvocs);
%         RR = 0;
%         if NUVT_here ~= NUVT
%             printf('Warning: the individual call type matrices (Original, rand and randsound) do not have the same number of categories');
%         end
%         confusion_matrix_Rand = zeros(NUVT, NUVT);
%         confusion_matrix_RandS = zeros(NUVT, NUVT);
%         for vtR=1:NUVT
%             selectorIndR = RR+1 : RR+Nvocs(vtR);
%             CC=0;
%             for vtC = 1:NUVT
%                 selectorIndC= CC+1 : CC+Nvocs(vtC);
%                 confusion_matrix_Rand(vtR,vtC)=sum(sum(MATRand_woDiag(selectorIndR, selectorIndC)));
%                 confusion_matrix_RandS(vtR,vtC)=sum(sum(MATRandS_woDiag(selectorIndR, selectorIndC)));
%                 CC = CC + Nvocs(vtC);
%             end
%             RR=RR+Nvocs(vtR);
%         end
%         
%         % Divide by the new number of events of the matrix wo diag
%         confusion_matrix_Rand = confusion_matrix_Rand ./ sum(sum(confusion_matrix_Rand));
%         confusion_matrix_RandS = confusion_matrix_RandS ./ sum(sum(confusion_matrix_RandS));
%        %Calculate MI and store temporarily values
%        MI_CM_Rand_local(rm) = info_matrix(confusion_matrix_Rand);
%        MI_CM_RandS_local(rm) = info_matrix(confusion_matrix_RandS);
%     end
%     % Calculate mean and SD values of MI for that cell
%     MeanMI_rand_wodiag(cm) = mean(MI_CM_Rand_local);
%     SDMI_rand_wodiag(cm) = sd(MI_CM_Rand_local);
%     MeanMI_randSound_wodiag(cm) = sd(MI_CM_RandS_local);
%     SDMI_randSound_wodiag(cm) = sd(MI_CM_RandS_local);
end
 
save('IndexSelectivity_confmat.mat', 'IS', 'VT', 'DI', 'DI_tot', 'DI_error', 'VS', 'MIBCT_OptCT_woDiag',  'MIBCT_OptCT_woDiag2', 'MIB_OptCT_woDiag','MeanMI_rand_wodiag','SDMI_rand_wodiag','MeanMI_randSound_wodiag','SDMI_randSound_wodiag', 'II', 'NV_cat', 'II_uni', 'II_uni_DI')
MIBCT_OptCT_woDiag_pSpike = MIBCT_OptCT_woDiag ./ MeanSR;
MIB_OptCT_woDiag_pSpike = MIB_OptCT_woDiag ./ MeanSR;


%% Plot the histograme of the number of vocalizations that are correctly classified above a varying threshold from 0.1 to 0.9
figure(2)
HIS = zeros(9,10);
for ll=1:9
    %subplot(3,3,ll)
    HIS(ll,:)=hist(VS(SemanticNeurons,ll),10);
    %title(sprintf('p>%.1f', ll*0.1));
end
bar3(HIS)
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', 0.1*(1:9));
ylabel('Probability of correct classification')

set(gca(), 'Xtick', 1:10);
set(gca(), 'XTickLabel', 0:9);
xlabel('Number of semantic categories')
zlabel('Number of units')

%% Plot for the number of cells that correctly classify the different categories of vocalization with various threshold values
figure(3)
NumCell = zeros(9,size(DI,2));
for ll=1:9
    %subplot(3,3,ll)
    NumCell(ll,:)=sum(DI>0.1*ll);
    %fig=bar(percCell);
    %title(sprintf('%% of discriminating cells with p>%.1f', ll*0.1))
    %for xx=1:9
    %    text(xx, -10*max(percCell)/100,UVT(xx))
    %end
end
bar3(NumCell)
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', 0.1*(1:9));
ylabel('Probability of correct classification')
set(gca(), 'Xtick', 1:9);
set(gca(), 'XTickLabel', VT{1}(1:end-1));
xlabel('Number of semantic categories')
zlabel('Number of units')

%% Plot values of MI for the individual call matrix and semantic matrix both without diagonals
figure(4)
plot(MIB_OptCT_woDiag_pSpike, MIBCT_OptCT_woDiag_pSpike,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of individual call matrix wo Diag (bits/spike)')
ylabel('Mutual Information of semantic group matrix wo Diag (bits/spike)')
axis([0 3 0 0.14])

figure(5)
plot(MIBCT_OptCT_pSpike, MIBCT_OptCT_woDiag_pSpike,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of semantic group matrix (bits/spike)')
ylabel('Mutual Information of semantic group matrix wo Diag (bits/spike)')
axis([0 0.25 0 0.25])

%% Investigate invariance of cells
% Invariance for category i is defined as the ratio of the sum of the proba off diagonal
% in the sub region of the category i in the individual call matrix divided
% by the same value in the subregion with a uni distribution of proba, in
% other words II_uni = DI_error/(DI_tot*(NV_cat^2 -NV_cat)/NV_cat^2

figure(6)
HIS = zeros(9,21);
for ll=1:9
    HIS(ll,:)=hist(II(SemanticNeurons,ll),-1:0.1:1);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', -1:0.1:1);
xlabel('Invariance Index')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')

figure(7)
subplot(1,2,1)
vectorbins= 0:0.1:1.3;
HIS = zeros(9,length(vectorbins));
for ll=1:9
    HIS(ll,:)=hist(II_uni(SemanticNeurons,ll),vectorbins);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', vectorbins);
xlabel('Invariance Index Uni')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')
title('Semantic Units')

subplot(1,2,2)
for ll=1:9
    HIS(ll,:)=hist(II_uni(NonSemanticNeurons,ll),vectorbins);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', vectorbins);
xlabel('Invariance Index Uni')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')
title('Non-Semantic Units')

figure(8)
subplot(1,2,1)
vectorbins=0:0.1:1;
HIS = zeros(9,length(vectorbins));
for ll=1:9
    HIS(ll,:)=hist(DI(SemanticNeurons,ll),vectorbins);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', vectorbins);
xlabel('Probability of correct classification')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')
title('Semantic Units')

subplot(1,2,2)
for ll=1:9
    HIS(ll,:)=hist(DI(NonSemanticNeurons,ll),vectorbins);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', vectorbins);
xlabel('Probability of correct classification')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')
title('Non-Semantic Units')


figure(9)
subplot(1,2,1)
vectorbins = 0:0.1:1;
HIS = zeros(9,length(vectorbins));
for ll=1:9
    HIS(ll,:)=hist(II_uni_DI(SemanticNeurons,ll),vectorbins);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', vectorbins);
xlabel('II*PCC')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')
title('Semantic Units')

subplot(1,2,2)
for ll=1:9
    HIS(ll,:)=hist(II_uni_DI(NonSemanticNeurons,ll),vectorbins);
end
bar3(HIS)
set(gca(), 'Xtick', 1:21);
set(gca(), 'XTickLabel', vectorbins);
xlabel('II*PCC')
set(gca(), 'Ytick', 1:9);
set(gca(), 'YTickLabel', VT{1}(1:end-1));
ylabel('Semantic categories')
zlabel('Number of units')
title('Non-Semantic Units')





%% Plot the values of IS
 %% mean LRI  and mean SSI for each call category
 % select semantic units
 semanticUnits=pv_MIBCT_OptCTSound<0.01;
 IS_sem= IS(find(semanticUnits),:);
 IS_nonsem = IS(find(~semanticUnits),:);
 figure(1)
 UVT = VT{1};
 NUVT = length(VT{1});
 Mean_IS=NaN(NUVT,1);
 Std_IS=NaN(NUVT,1);
 Nunits = size(IS_sem,1);
 for cc=1:NUVT
     IS_finite=find(isfinite(IS_sem(:,cc)));
     if length(IS_finite)~=Nunits
         fprintf('%d infinite index of selectivity values calculated for call category %s\n', Nunits-length(IS_finite), UVT{cc});
     end
     Mean_IS(cc)=nanmean(IS_sem(IS_finite,cc));
     Std_IS(cc)=nanstd(IS_sem(IS_finite,cc));
 end
 bar(Mean_IS);
 hold on
 errorbar(Mean_IS, Std_IS, 'b.');
 set(gca(), 'Xtick', 1:NUVT);
 set(gca(), 'XTickLabel', UVT);
 
 %% representing, for each call category the percentage of multiunits that have a IS above a certain threshold

 figure(2)
 sp =0;
 PercSem_Sel=NaN(NUVT,1);
 PercSem_nonSel=NaN(NUVT,1);
 PercnonSem_Sel=NaN(NUVT,1);
 PercnonSem_nonSel=NaN(NUVT,1);
     for cc=1:NUVT
        UnitsIdx=isfinite(IS_sem(:,cc));
        NbUnits_local=sum(UnitsIdx);
        UnitsIdx=find(UnitsIdx);
        IS_local=IS_sem(UnitsIdx,cc);
        UnitsIdxNS=isfinite(IS_nonsem(:,cc));
        NbUnits_localNS=sum(UnitsIdxNS);
        UnitsIdxNS=find(UnitsIdxNS);
        IS_localNS=IS_nonsem(UnitsIdxNS,cc);
        PercSem_Sel(cc)=100*sum(IS_local>=1)/NbUnits_local;
        PercSem_nonSel(cc)=100*sum(IS_local<1)/NbUnits_local;
        PercnonSem_Sel(cc)=100*sum(IS_localNS>=1)/NbUnits_localNS;
        PercnonSem_nonSel(cc)=100*sum(IS_localNS<1)/NbUnits_localNS;
     end
     subplot(2,2,1)
     bar(PercSem_Sel);
     title('Percentage of selective semantic units (IS>=1)');
     set(gca(), 'Xtick', 1:NUVT);
     set(gca(), 'XTickLabel', UVT);
     subplot(2,2,2)
     bar(PercSem_nonSel);
     title('Percentage of non selective semantic units (IS<1)');
     set(gca(), 'Xtick', 1:NUVT);
     set(gca(), 'XTickLabel', UVT);
     subplot(2,2,3)
      bar(PercnonSem_Sel);
     title('Percentage of selective non-semantic units (IS>=1)');
     set(gca(), 'Xtick', 1:NUVT);
     set(gca(), 'XTickLabel', UVT);
     subplot(2,2,4)
      bar(PercnonSem_nonSel);
     title('Percentage of non-selective non-semantic units (IS<1)');
     set(gca(), 'Xtick', 1:NUVT);
     set(gca(), 'XTickLabel', UVT);
     
%% plotting the percentage of semantic selective, semantic non selective and non semantic cells in each zone
% find selective and non-selectivecells
Sel_units=sum(IS>1,2)>0;
nonSel_units=~Sel_units;
Sem_sel_units = semanticUnits .* Sel_units;
Sem_nonsel_units = semanticUnits .* nonSel_units;
nonSem_units = ~semanticUnits;

ZO=unique(ZONES(:,1));
NZO = length(ZO);
pp=1;
PercCell = zeros(NZO, 3);
for ii = 1:NZO
    ss=ZO(ii);
    Local_indices=ZONES(:,1)==ss;
    PercCell(ii,1) = sum(nonSem_units .* Local_indices)/sum(Local_indices);
    PercCell(ii,2) = sum(Sem_nonsel_units .* Local_indices)/sum(Local_indices);
    PercCell(ii,3) = sum(Sem_sel_units .* Local_indices)/sum(Local_indices);
end
figure(4)
bar(PercCell, 'stacked')
colormap(cool)
set(gca(), 'Xtick', 1:NZO);
set(gca(), 'XTickLabel', ['Unknown' ZONES_List]);
legend('NonSemantic', 'Semantic non Selective', 'Semantic Selective');

   