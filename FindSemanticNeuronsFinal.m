%% This script gather data from all single units and is followed by FindSemanticNeuronsFinalPredict.m
resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')
NbCell_estimate = 1500;
FigFlag =0;

%% Retrieve spike shape info
%cd /Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/SpikeShape
cd /auto/fhome/julie/matlab/tlab/src/h5analysis/Julie_neuralcode/SpikeShape
Ldir=pwd;
SpikeShape=load(strcat(Ldir,'/spikeShapeResults.mat'),'allNames', 'indG5', 'spikeType');

%% Retrieve Histology info
cd /auto/fhome/julie/Documents/Histo
% cd /Volumes/FREDERIC/HistoJulie/Stacks
AnaLdir=pwd;
AnatTxt=dir(fullfile(AnaLdir, 'List_h5files*Histo.txt'));
Anat_root = cell(NbCell_estimate, 5);
nbTxt = length(AnatTxt);
Lanat = 0;
for tt=1:nbTxt
    textfile = AnatTxt(tt).name;
    fid = fopen(fullfile(AnaLdir, textfile));
    dataAnat = textscan(fid, '%s\t%s\t%f\t%f\t%s');
    anatnb = size(dataAnat{1}, 1);
    for ll = 1: anatnb
        Lanat=Lanat+1;
        for col=1:5
            Anat_root{Lanat,col}=dataAnat{col}(ll);
        end
    end
end
Anat_root = Anat_root(1:Lanat,:);


%% Setup some output variables
% cell containing the path of each unit
List_matfilepath = cell(NbCell_estimate, 1);

% cell containing the anatomy info of each unit
List_anat = cell(NbCell_estimate, 5);

% vectors containing optimal windows for each unit
optWin=zeros(NbCell_estimate,1); % based on the sound file matrix
optWin_BGwo=zeros(NbCell_estimate,1); % based on the sound file matrix without BG
optWin_CT= zeros(NbCell_estimate,1); % based on the calls type matrix

% vectors containing mutual information for each unit
MI_confB=zeros(NbCell_estimate,9);
MI_confBBGwo=zeros(NbCell_estimate,9);
MI_confBCT=zeros(NbCell_estimate,9);
MI_confBCTBGwo=zeros(NbCell_estimate,9);

% vector containing the mutual information of the best matrix for each unit
MIB_OptCT = zeros(NbCell_estimate,1);
MIBCT_OptCT = zeros(NbCell_estimate,1);
MIB_Opt = zeros(NbCell_estimate,1);
MIBCT_Opt = zeros(NbCell_estimate,1);

% vector containing the size of the confusion matrices for each unit
ConfSize = zeros(NbCell_estimate,2);
ConfSizeCT = zeros(NbCell_estimate,2);

% vector containing the diagonal of the best call type confusion matrix
DiagCMCT = zeros(NbCell_estimate,10); %only take into account call categories (no mln) and no whines since data for only 1 bird for whines + Bg category

% Vectors containing the mutual information per area for all units
MI_tot=zeros(NbCell_estimate,1);
MI_diag=zeros(NbCell_estimate,1);
MI_error=zeros(NbCell_estimate,1);
MI_diag_uni=zeros(NbCell_estimate,1);
MI_all_error_uni=zeros(NbCell_estimate,1);
MI_diag_uni_cat=zeros(NbCell_estimate,1);
MI_real_error_uni=zeros(NbCell_estimate,1);
MI_diag_uni_cat_maxInv = zeros(NbCell_estimate,1);
MI_diag_uni_cat_minInv = zeros(NbCell_estimate,1);

AvMI_tot_Rand = zeros(NbCell_estimate,2); % first column is for rand matrices, second column for rand matrices with Background sections not shuffled
SDMI_tot_Rand = zeros(NbCell_estimate,2);
AvMI_diag_uni_Rand=zeros(NbCell_estimate,2);
SDMI_diag_uni_Rand=zeros(NbCell_estimate,2);
AVMI_diag_uni_cat_Rand=zeros(NbCell_estimate,2);
SDMI_diag_uni_cat_Rand=zeros(NbCell_estimate,2);
MI_uni_diag_cat_Rand = zeros(NbCell_estimate, 100);
MI_uni_diag_cat_RandBG = zeros(NbCell_estimate, 100);


PCC_cat=cell(NbCell_estimate,1);
IS3 = cell(NbCell_estimate,1);
IS = zeros(NbCell_estimate,1);
IS2 = IS;
II2 = zeros(NbCell_estimate,1);
VS = zeros(NbCell_estimate,9);

% % Vectors for Global R2A Models resluts 
% Acoustic = zeros(NbCell_estimate,1);
% Semantic = zeros(NbCell_estimate,1);
% AcSem = zeros(NbCell_estimate,1);
% PvalueAcSem = Acoustic;
% PvalueSemAc = Acoustic;
% 
% %Matrices for R2A values
% RAcoustic = zeros(NbCell_estimate,140);
% RSemantic = RAcoustic;
% RAcSem = RAcoustic;
% 
% % Vectors for sum of squares of residuals (error) of each model and F
% % statistic
% SSAcoustic = zeros(NbCell_estimate,1);
% SSSemantic = zeros(NbCell_estimate,1);
% SSAcSem = zeros(NbCell_estimate,1);
% SStot=zeros(NbCell_estimate,1);
% FAcSem = zeros(NbCell_estimate,1);
% FSemAc = zeros(NbCell_estimate,1);
% 
% %MAtrices for nb of stim and number of cat on which predicting models are
% %calculated
% RNb_Stim = RAcoustic;
% RNb_Cat = RAcoustic;
% RNb_PC = RAcoustic;

% Vector of the mean spiking rate for each unit calculated from WholeVoc files on first extracts
MeanSR = zeros(NbCell_estimate,1);

% Vector of the spike shape and sort quality
Spike_shape = zeros(NbCell_estimate,1);

% Vector of the code for the subject
SUBJ = zeros(NbCell_estimate,2);
SUBJECTS = {'YelBlu6903F' 'BlaBro09xxF' 'LblBlu2028M' 'GreBlu9508M' 'WhiBlu5396M' 'WhiWhi4522M'};

%Vector of the code for the anatomical zones
ZONES = zeros(NbCell_estimate,2);
ZONES_List = {'CMM', 'CML', 'L1', 'L2A', 'L2B', 'L3', 'L', 'NCM', 'HP'};

%% Start to loop through individuals and files and extract information
cd /auto/k6/julie/matfile
input_dir=pwd;
Subjects = dir(input_dir);
ii=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        
        sprintf('Harvesting data of %s\n', Indiv)
        % retrieve Confusion matrices files
        matfiles=dir(fullfile(Idir,'ConfMat*.mat'));
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        Indices=find(SS_Ind);
        LM=length(Indices);
        
        %retrieve Random confusion matrices files
        Randmatfiles=dir(fullfile(Idir, 'RandMat*.mat'));
        rm=length(Randmatfiles);
        SS_Ind=zeros(rm,1);
        for ff=1:rm
            if ~isempty(strfind(Randmatfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        IndicesRM=find(SS_Ind);
        RM=length(IndicesRM);
        
        % retrieve WholeVoc files to calculate Mean rate
        WVmatfiles = dir(fullfile(Idir, 'WholeVoc*.mat'));
        wm=length(WVmatfiles);
        SS_Ind=zeros(wm,1);
        for ff=1:wm
            if ~isempty(strfind(WVmatfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        IndicesWV=find(SS_Ind);
        WV=length(IndicesWV);
        
        % retrieve Models files to extract Global R2A
        MDmatfiles = dir(fullfile(Idir, 'Models*.mat'));
        md = length(MDmatfiles);
        SS_Ind=zeros(md,1);
        for ff=1:md
            if ~isempty(strfind(MDmatfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        IndicesMD=find(SS_Ind);
        MD=length(IndicesMD);
        
        % check that we have the same number of files for all of them
        if LM~=RM
            sprintf('WARNING: the nb of observed matrices (%d) is different from the nb of randfiles (%d)\n', LM, RM)
            break
        end
        if LM~=WV
            sprintf('WARNING: the nb of observed matrices (%d) is different from the nb of WholeVocfiles (%d)\n', LM, WV)
            break
        end
        if LM~=MD
            sprintf('WARNING: the nb of observed matrices (%d) is different from the nb of Modelsfiles (%d)\n', LM, MD)
            break
        end
        
  %% loop through files and calculate all variables
        for hh=1:LM
            Matfile = matfiles(Indices(hh)).name;
            sprintf('Loading %s\n', Matfile)
            MatfilePath=fullfile(Idir, Matfile);
            MAT=load(MatfilePath);
            
            ii=ii+1;
            %% keep track of the site ID
            List_matfilepath{ii}=MatfilePath;
            
            %% Code the subject sex and identity
            sprintf('Subject ID\n')
            if strcmp(MAT.subject, 'BlaBro09xxF')
                SUBJ(ii,:)=[1 0];
            elseif strcmp(MAT.subject, 'LblBlu2028M')
                SUBJ(ii,:)=[2 1];
            elseif strcmp(MAT.subject, 'GreBlu9508M')
                SUBJ(ii,:)=[3 1];
            elseif strcmp(MAT.subject, 'WhiBlu5396M')
                SUBJ(ii,:)=[4 1];
            elseif strcmp(MAT.subject, 'WhiWhi4522M')
                SUBJ(ii,:)=[5 1];
            end
            
            %% Construct anatomy vectors
            sprintf('Anatomy\n')
            for ff = 1:size(Anat_root, 1)
                TempID = cell2mat(Anat_root{ff, 1});
                if strcmp(Matfile(9:end-8), TempID(1:end-3))
                    List_anat(ii, :) = Anat_root(ff, :);
                end
            end

            if strcmp(List_anat{ii,2}, 'CMM')
                ZONES(ii,1)=1;
            elseif strcmp(List_anat{ii,2}, 'CML')
                ZONES(ii,1)=2;
            elseif strcmp(List_anat{ii,2}, 'L1')
                ZONES(ii,1)=3;
            elseif strcmp(List_anat{ii,2}, 'L2A')
                ZONES(ii,1)=4;
            elseif strcmp(List_anat{ii,2}, 'L2B')
                ZONES(ii,1)=5;
            elseif strcmp(List_anat{ii,2}, 'L3')
                ZONES(ii,1)=6;
            elseif strcmp(List_anat{ii,2}, 'L')
                ZONES(ii,1)=7;
            elseif strcmp(List_anat{ii,2}, 'NCM')
                ZONES(ii,1)=8;
            elseif strcmp(List_anat{ii,2}, 'HP')
                ZONES(ii,1)=9;
            end
            if strcmp(List_anat{ii,5}, 'L')%lef or right hemisphere
                ZONES(ii,2)=1;
            end
            
            
            %% Collect info about the observed matrices
            sprintf('Observed Confusion Matrices\n')
            Winsize=MAT.winSize;
            MAXWin=find(MAT.mi_confusionB==max(MAT.mi_confusionB));
            MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
            MAXWin_BGwo=find(MAT.mi_confusionB_BGwo==max(MAT.mi_confusionB_BGwo));
            optWin(ii)=Winsize(MAXWin);
            optWin_BGwo(ii)=Winsize(MAXWin_BGwo);
            optWin_CT(ii)=Winsize(MAXWinCT);
            MI_confB(ii,:)=MAT.mi_confusionB;
            MI_confBBGwo(ii,:)=MAT.mi_confusionB_BGwo;
            MI_confBCT(ii,:)=MAT.mi_confusionCT;
            MI_confBCTBGwo(ii,:)=MAT.mi_confusionCT_BGwo;
            ConfSize(ii,:) = size(MAT.confusionMatrix{1});
            ConfSizeCT(ii,:) = size(MAT.confusionMatrixCT{1});
            MIB_OptCT(ii) = MI_confB(ii,find(Winsize==optWin_CT(ii)));
            MIB_Opt(ii) = MI_confB(ii,find(Winsize==optWin(ii)));
            MIBCT_OptCT(ii) = MI_confBCT(ii,find(Winsize==optWin_CT(ii)));
            MIBCT_Opt(ii) = MI_confBCT(ii,find(Winsize==optWin(ii)));

            %% Collect CT Matrices Diagonals
            DiagCMCT(ii,:)=diag(MAT.confusionMatrixCT{find(Winsize==optWin_CT(ii))});

            %% Calculate values of MI per area for the observed individual
            % vocalization matrices
            sprintf('Calculate values of MI per area in observed Matrix\n')
            Mat = MAT.confusionMatrix{find(Winsize==optWin_CT(ii))}; % we are choosing the window size that gives the best value of MI conf in the CT matrix
            VocTypeConfMat = MAT.confusionMatrixCT{find(Winsize==optWin_CT(ii))};
            SF_ConfMat = Mat.*MAT.neventMatrix(find(Winsize==optWin_CT(ii)));
            
            if FigFlag>0
                figure(1)
                subplot(1,2,1)
                imagesc(Mat)
                colorbar;
                xlabel('Predicted vocalization');
                ylabel('Actual vocalization');
                title(sprintf('Single vocalizations\n'));
                subplot(1,2,2)
                imagesc(VocTypeConfMat)
                colorbar;
                xlabel('Predicted vocalization');
                ylabel('Actual vocalization');
                title(sprintf('Voc Type \n'));
            end
            
            % construct the cell array of indices for each call category
            VocTypeSel = MAT.VocTypeSel;
            StimTypeCM=unique(VocTypeSel);
            IBG = find(strcmp(StimTypeCM, 'BG'));
            StimTypeCM = [StimTypeCM(1:(IBG-1)); StimTypeCM((IBG+1):end); StimTypeCM(IBG)];
            NstimTypeCM=length(StimTypeCM);
            cat = cell(NstimTypeCM,1);
            for cc = 1:NstimTypeCM
                cat{cc}=find(strcmp(StimTypeCM(cc), VocTypeSel));
            end
        
            %Calculate the values of MI in the different areas of the IV
            %confusion matrix
            [ mi_tot, mi_tot2, mi_diag, mi_error, mi_diag_uni, mi_all_error_uni, mi_diag_uni_cat, mi_real_error_uni, mi_diag_uni_cat_maxInv, mi_diag_uni_cat_minInv]=info_matrix_perarea_ext(Mat, cat, FigFlag);
            MI_tot(ii)=mi_tot;
            MI_diag(ii)=mi_diag;
            MI_error(ii)=mi_error;
            MI_diag_uni(ii)=mi_diag_uni;
            MI_all_error_uni(ii)=mi_all_error_uni;
            MI_diag_uni_cat(ii)=mi_diag_uni_cat;
            MI_real_error_uni(ii)=mi_real_error_uni;
            MI_diag_uni_cat_maxInv(ii) = mi_diag_uni_cat_maxInv;
            MI_diag_uni_cat_minInv(ii) = mi_diag_uni_cat_minInv;
        
            %% Calculate the alternative Index of Invariance II2
            % For all the category zones calculate observed entropy then max entropy if the joint probability was smoothed within each row 
            sprintf('Index of Invariance\n')
            Hobs_Inv = 0;
            Hmax_Inv = 0;
            nstim = size(Mat, 1);
            Ncat = length(cat);
            for aa = 1:nstim
                ProbaCat = 0;
                for nc = 1:Ncat
                        icat = cat{nc};
                        if ~isempty(intersect(aa, icat))
                            break
                        end
                end
                nstim_local = length(icat);
                for bb = 1:nstim_local;
                    Hobs_Inv = Hobs_Inv + Mat(aa,bb) * log2(Mat(aa,bb));
                    ProbaCat = ProbaCat + Mat(aa,bb);
                end
                Hmax_Inv = Hmax_Inv + ProbaCat*log2(ProbaCat/nstim_local);
            end
            II2(ii) = Hobs_Inv/Hmax_Inv;

            %% Calculate an index of selectivity for the different call categories for that unit
            sprintf('Indices of Selectivity\n')
            % First calculate the CT Matrix
            confusion_matrix_CallType = zeros(NstimTypeCM, NstimTypeCM);
            for vtR=1:NstimTypeCM
                stR=StimTypeCM(vtR);
                selectorR=strcmp(MAT.VocTypeSel, stR);
                selectorIndR=find(selectorR);
                for vtC = 1:NstimTypeCM
                    stC=StimTypeCM(vtC);
                    selectorC=strcmp(MAT.VocTypeSel, stC);
                    selectorIndC=find(selectorC);
                    confusion_matrix_CallType(vtR,vtC)=sum(sum(Mat(selectorIndR, selectorIndC)));
                end
            end

            %Convert joint probabilities to conditional probailities
            ProbaCat = sum(confusion_matrix_CallType, 2);
            repProbaCat = repmat(ProbaCat, 1, size(confusion_matrix_CallType,2));
            confusion_matrix_CallType_cond = confusion_matrix_CallType ./ repProbaCat;

            % isolate the diagonal vectors containing the percentage of correct
            % classification of each call category
            Vect_class_cond = diag(confusion_matrix_CallType_cond,0);
            PCC_cat{ii} = Vect_class_cond; %here we still have the background

            % calculating IS index of selectivity based on the entropy of the
            % diagonal of the CT matrix compare to the entropy of a uni diag
            % CT matrix without taking the background into account
            Vect_class = diag(confusion_matrix_CallType,0);
            Vect_class2 = Vect_class./sum(Vect_class);
            Hobs = 0;
            Hmax = sum(Vect_class(1:(end-1)))*log2(sum(Vect_class(1:(end-1)))/(NstimTypeCM-1));
            Hobs2 = 0;
            Hmax2 = -log2(NstimTypeCM-1);
            zero_ind = Vect_class == 0;
            Vect_class_for_entropy = Vect_class;
            Vect_class_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
            Vect_class_for_entropy2 = Vect_class2;
            Vect_class_for_entropy2(zero_ind) = 1; 
            for nt = 1:(NstimTypeCM-1)
                Hobs = Hobs + Vect_class_for_entropy(nt)*log2(Vect_class_for_entropy(nt));
                Hobs2 = Hobs2 + Vect_class_for_entropy2(nt)*log2(Vect_class_for_entropy2(nt));
            end
            IS(ii) = 1 - Hobs/Hmax;
            IS2(ii) = 1- Hobs2/Hmax2;% Normalized probabilities so to calculate more accurately the entropy in the diagonal because entropy calculation is complete when sum(proba)=1


            % calculating IS3 index of selectivity3 on the matrix diagonal
            % without taking the background
            IS3_temp=zeros((NstimTypeCM-1),1);
            for vt=1:(NstimTypeCM-1)
                Vect_class_temp = Vect_class_cond(1:(end-1));
                Vect_class_temp(vt)=[];
                IS3_temp(vt)=log2(Vect_class_cond(vt)/mean(Vect_class_temp));
            end
            IS3{ii}=IS3_temp;

            % Calculating the number of categories for which p>0.1*ll
            for ll=1:9
                VS(ii,ll)=sum(Vect_class_cond(1:end-1)>0.1*ll);
            end
        

 %% Upload the random matrices and calculate MI per area
            sprintf('Upload Random matrices and caluclate MI per area\n')
            for rr=1:RM
                RandMatfile = Randmatfiles(IndicesRM(rr)).name;
                if strcmp(Matfile(9:end), RandMatfile(9:end))
                    RMAT = load(fullfile(Idir, RandMatfile));
                end
            end
            rt = 1;
            while rt<=2
                if rt==1
                    LRM=length(RMAT.CM_IV_Rand);%number of random matrices
                elseif rt==2
                    LRM=length(RMAT.CM_IV_RandBG);%number of random matrices BG fixed
                end
                MI_tot_rand=zeros(LRM,1);
                MI_diag_uni_rand=zeros(LRM,1);
                MI_diag_uni_cat_rand=zeros(LRM,1);
                for rm =1:LRM
                    if rt==1
                        Rmat = RMAT.CM_IV_Rand{rm};
                    elseif rt==2
                        Rmat = RMAT.CM_IV_RandBG{rm};
                    end
                    
                    % construct the cell array of indices for each random category
                    if rt==1 % random matrices
                        Nb_VocPerCat = RMAT.Nb_VocPerCat{rm};
                    elseif rt==2 % random BG Matrices
                        Nb_VocPerCat = RMAT.Nb_VocPerCatBG{rm};
                    end
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
                AvMI_tot_Rand(ii,rt) = mean(MI_tot_rand);
                SDMI_tot_Rand(ii,rt) = std(MI_tot_rand);
                fprintf('the MI_confusion of the IV Matrix is %f for that cell\nThe average MI of random matrices is %f+/-%f\n', MI_tot(ii), AvMI_tot_Rand(ii),SDMI_tot_Rand(ii))
            %pause
                AvMI_diag_uni_Rand(ii,rt)=mean(MI_diag_uni_rand);
                SDMI_diag_uni_Rand(ii,rt)=std(MI_diag_uni_rand);
                AVMI_diag_uni_cat_Rand(ii,rt)=mean(MI_diag_uni_cat_rand);
                SDMI_diag_uni_cat_Rand(ii,rt)=std(MI_diag_uni_cat_rand);
                if rt==1
                    MI_uni_diag_cat_Rand(ii,:) = MI_diag_uni_cat_rand;
                elseif rt==2
                    MI_uni_diag_cat_RandBG(ii,:) = MI_diag_uni_cat_rand;
                end
                rt = rt +1;
            end
            
            
%% find the mean spike rate of the unit
            sprintf('Mean spike rate\n')
            for ww=1:WV
                WVMatfile = WVmatfiles(IndicesWV(ww)).name;
                if strcmp(Matfile(9:end), WVMatfile(10:end))
                    WMAT = load(fullfile(Idir, WVMatfile));
                end
            end
            FirstWholeVoc = find(WMAT.Voc_orders==1);
            MeanSR(ii)=nanmean(cell2mat(WMAT.MeanRate(FirstWholeVoc)));
            
%% find the Global R2A from the models for that unit MDmatfiles
%THIS STEP IS NOW DONE BY FINDSEMANTICNEURONSFINAL_ENCODINGM.M BESIDES THE
%F CALCULATION HERE IS FALSE
%             sprintf('Global R2A from models\n')
%             for mm=1:MD
%                 MDMatfile = MDmatfiles(IndicesMD(mm)).name;
%                 if strcmp(Matfile(9:end), MDMatfile(8:end))
%                     MMAT = load(fullfile(Idir, MDMatfile));
%                 end
%             end
%             
%             Acoustic(ii)=MMAT.GlobalR2A.Acoustic;
%             Semantic(ii)=MMAT.GlobalR2A.Semantic;
%             AcSem(ii)=MMAT.GlobalR2A.AcSem;
%             if ii==1
%                 RAcoustic = RAcoustic(:,1:length(MMAT.Wins));
%                 RSemantic = RSemantic(:,1:length(MMAT.Wins));
%                 RAcSem = RAcSem(:, 1:length(MMAT.Wins));
%                 Wins=MMAT.Wins;
%                 RNb_Stim = RNb_Stim(:,1:length(MMAT.Wins));
%                 RNb_Cat = RNb_Cat(:,1:length(MMAT.Wins));
%                 RNb_PC = RNb_PC(:,1:length(MMAT.Wins));
%             end
%             RAcoustic(ii,:) = MMAT.R2A.Acoustic;
%             RSemantic(ii,:) = MMAT.R2A.Semantic;
%             RAcSem(ii,:) = MMAT.R2A.AcSem;
%             
%             % extract the p-value of bootsrap differences Global R2A
%             PvalueAcSem(ii) =  MMAT.GlobalR2A.PvalAcSem;
%             PvalueSemAc(ii) = MMAT.GlobalR2A.PvalSemAc;
%             
%             % extract the Sum of squares Error (residuals) and total for
%             % the 3 models
%             SSAcoustic(ii) = nansum(MMAT.SSres.Acoustic);
%             SSSemantic(ii) = nansum(MMAT.SSres.Semantic);
%             SSAcSem(ii) = nansum(MMAT.SSres.AcSem);
%             SStot(ii) = nansum(MMAT.SStot);
%             
%             % extract the number of categories and nb of stims on wich the
%             % model are calculated for each time point
%             for vv=1:length(MMAT.Wins)
%                 RNb_Cat(ii,vv)=length(unique(MMAT.voc{vv}));
%                 RNb_Stim(ii,vv)=length(MMAT.voc{vv});
%                 RNb_PC(ii,vv)= MMAT.Best_nbPC(vv);
%             end
%             
%             % Calculate F statistic for SS
%             NT = sum(RNb_Stim(ii,:)); % # observations = sum of the number fo stims over all the time points
%             NDL = max(RNb_Cat(ii,:)); % # parameters semantic models = max number of categories used in the semantic model over all time points
%             NPC = nansum(RNb_PC(ii,:)); % # parameters acoustic models = sum of all differents PC used over time
%             FAcSem(ii) = (SSAcoustic(ii) - SSAcSem(ii))*(NT - NPC - NDL)/(NDL*SSAcSem(ii));
%             FSemAc(ii) = (SSSemantic(ii) - SSAcSem(ii))*(NT - NPC - NDL)/(NPC*SSAcSem(ii));
            
%% find the spike shape of the unit
            sprintf('Spike shape\n')
            uu=1;
            while strcmp(Matfile(9:end-4), SpikeShape.allNames{uu}(31:end-3))==0
                uu=uu+1;
            end
            if ~isempty(find(SpikeShape.indG5==uu))
                sshape=SpikeShape.spikeType(find(SpikeShape.indG5==uu));
                if strcmp(sshape, 'unclassified')
                    Spike_shape(ii)=1;
                elseif strcmp(sshape, 'large')
                    Spike_shape(ii)=2;
                elseif strcmp(sshape, 'narrow')
                    Spike_shape(ii)=3;
                elseif strcmp(sshape, 'wide')
                    Spike_shape(ii)=4;
                end
            end
        end
    end
end

%% Compile results into structures
sprintf('Compile results\n')
List_matfilepath=List_matfilepath(1:ii);
List_anat=List_anat(1:ii,:);
SUBJ = SUBJ(1:ii,:);
ZONES = ZONES(1:ii,:);
Spike_shape = Spike_shape(1:ii);
MeanSR = MeanSR(1:ii);


Conf.optWin=optWin(1:ii);
Conf.optWin_CT=optWin_CT(1:ii);
Conf.optWin_BGwo=optWin_BGwo(1:ii);

Conf.MI_confB=MI_confB(1:ii,:);
Conf.MI_confBBGwo=MI_confBBGwo(1:ii,:);
Conf.MI_confBCT=MI_confBCT(1:ii,:);
Conf.MI_confBCTBGwo=MI_confBCTBGwo(1:ii,:);
Conf.ConfSize = ConfSize(1:ii, :);
Conf.ConfSizeCT = ConfSizeCT(1:ii, :);
Conf.MIB_OptCT = MIB_OptCT(1:ii);
Conf.MIB_Opt = MIB_Opt(1:ii);
Conf.MIBCT_OptCT = MIBCT_OptCT(1:ii);
Conf.MIBCT_Opt = MIBCT_Opt(1:ii);
Conf.DiagCMCT=DiagCMCT(1:ii,:);

MI_perArea.MI_tot=MI_tot(1:ii);
MI_perArea.MI_diag=MI_diag(1:ii);
MI_perArea.MI_error=MI_error(1:ii);
MI_perArea.MI_diag_uni=MI_diag_uni(1:ii);
MI_perArea.MI_all_error_uni=MI_all_error_uni(1:ii);
MI_perArea.MI_diag_uni_cat=MI_diag_uni_cat(1:ii);
MI_perArea.MI_real_error_uni=MI_real_error_uni(1:ii);
MI_perArea.MI_diag_uni_cat_maxInv = MI_diag_uni_cat_maxInv(1:ii);
MI_perArea.MI_diag_uni_cat_minInv = MI_diag_uni_cat_minInv(1:ii);

MI_perArea.AvMI_tot_Rand = AvMI_tot_Rand(1:ii,:); % first column is for rand matrices, second column for rand matrices with Background sections not shuffled
MI_perArea.SDMI_tot_Rand = SDMI_tot_Rand(1:ii, :);
MI_perArea.AvMI_diag_uni_Rand=AvMI_diag_uni_Rand(1:ii, :);
MI_perArea.SDMI_diag_uni_Rand=SDMI_diag_uni_Rand(1:ii,:);
MI_perArea.AVMI_diag_uni_cat_Rand=AVMI_diag_uni_cat_Rand(1:ii,:);
MI_perArea.SDMI_diag_uni_cat_Rand=SDMI_diag_uni_cat_Rand(1:ii,:);
MI_perArea.MI_uni_diag_cat_Rand = MI_uni_diag_cat_Rand(1:ii,:);
MI_perArea.MI_uni_diag_cat_RandBG = MI_uni_diag_cat_RandBG(1:ii,:);

Selectivity.PCC_cat = PCC_cat(1:ii);
Selectivity.DiagEntropy = IS(1:ii);
Selectivity.DiagEntropyNorm = IS2(1:ii);
Selectivity.DiagLRI = IS3(1:ii);
Selectivity.NbCatHighProba = VS(1:ii,:);
Invariance = II2(1:ii);

SemanticIndex.Observed = MI_perArea.MI_diag_uni_cat ./ MI_perArea.MI_tot;
SemanticIndex.Rand = MI_perArea.AVMI_diag_uni_cat_Rand(:,1) ./ MI_perArea.AvMI_tot_Rand(:,1);
SemanticIndex.RandBG = MI_perArea.AVMI_diag_uni_cat_Rand(:,2) ./ MI_perArea.AvMI_tot_Rand(:,2);
SemanticIndex.Pvalue = normpdf((MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_Rand(:,1)) ./ MI_perArea.SDMI_diag_uni_cat_Rand(:,1));
SemanticIndex.PvalueBG = normpdf((MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_Rand(:,2)) ./ MI_perArea.SDMI_diag_uni_cat_Rand(:,2));
SemanticIndex.NonLinear = (MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_Rand(:,1)) ./ MI_perArea.MI_tot;
SemanticIndex.NonLinearBG = (MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_Rand(:,2)) ./ MI_perArea.MI_tot;

% GlobalR2A.Acoustic = Acoustic(1:ii);
% GlobalR2A.Semantic = Semantic(1:ii);
% GlobalR2A.AcSem = AcSem(1:ii);
% GlobalR2A.PvalueAcSem = PvalueAcSem(1:ii);
% GlobalR2A.PvalueSemAc = PvalueSemAc(1:ii);
% R2A.Acoustic = RAcoustic(1:ii,:);
% R2A.Semantic = RSemantic(1:ii,:);
% R2A.AcSem = RAcSem(1:ii,:);
% R2A.Wins = Wins;
% R2A.RNb_Cat = RNb_Cat(1:ii,:);
% R2A.RNb_Stim = RNb_Stim(1:ii,:);
% R2A.RNb_PC = RNb_PC(1:ii,:);
% R2A.SSAcoustic = SSAcoustic(1:ii);
% R2A.SSSemantic = SSSemantic(1:ii);
% R2A.SSAcSem = SSAcSem(1:ii);
% R2A.SStot = SStot(1:ii);
% R2A.FAcSem = FAcSem(1:ii);
% R2A.FSemAc = FSemAc(1:ii);




%% Save data
save(fullfile(resultsDirectory, 'SemanticAnalysis.mat'))
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis.mat'))



