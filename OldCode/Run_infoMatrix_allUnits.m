FigFlag=1; %to see the IV uniform confusion matrices onwhich the calculus of MI are made
RandSwitch = 0; %Set to 1 if you want the calculations to run on the random matrices too

%addpath('/auto/k1/queued');
%addpath(genpath('/auto/fhome/julie/matlab/tlab'));
%addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%resultsDirectory='/auto/k8/julie';
 resultsDirectory='/Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat';
Res=load(fullfile(resultsDirectory, 'Results131130.mat'), 'List_matfilepath');


% retrieve data confusion matrices
 cd /Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat
%cd /auto/k6/julie/matfile/ConfMat
Idir=pwd;
LM=length(Res.List_matfilepath); % here I'm using the same order as the one used to calculate random matrices, spike shape and anatomy list by FindSemanticNeurons.m

% retrieve random matrices
 cd /Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat/RandMat
%cd /auto/k8/julie/RandMat
Rdir=pwd;
Rmatfiles=dir(fullfile(Rdir,'RandMat*.mat'));
LR=length(Rmatfiles);
if LR~=LM
    fprintf('Problem we don''t have the same number of random files and data files')
end

%retrieve correlation matrices
 cd /Users/elie/Documents/MATLAB/data/matfile/StimCorrelationMatrices
%cd /auto/k6/julie/matfile/StimCorrelationMatrices
Ldir=pwd;
CorrMatfiles=dir(fullfile(Ldir, 'stat*.mat'));
NCMfiles = length(CorrMatfiles);

MI_tot=zeros(LM,1);
MI_diag=zeros(LM,1);
MI_error=zeros(LM,1);
MI_diag_uni=zeros(LM,1);
MI_all_error_uni=zeros(LM,1);
MI_diag_uni_cat=zeros(LM,1);
MI_real_error_uni=zeros(LM,1);
MI_diag_uni_cat_maxInv = zeros(LM,1);
MI_diag_uni_cat_minInv = zeros(LM,1);
if RandSwitch==1
    AvMI_tot_Rand = zeros(LR,2); % first column is for rand sound matrices, second column for rand matrices
    SDMI_tot_Rand = zeros(LR,2);
    AvMI_diag_uni_Rand=zeros(LM,2);
    SDMI_diag_uni_Rand=zeros(LM,2);
    AVMI_diag_uni_cat_Rand=zeros(LM,2);
    SDMI_diag_uni_cat_Rand=zeros(LM,2);
    MI_uni_diag_cat_RandSound = zeros(LM, 100);
    MI_uni_diag_cat_Rand = zeros(LM, 100);
end

PCC_cat=cell(LM,1);
IS3 = cell(LM,1);
IS = zeros(LM,1);
IS2 = IS;
II2 = zeros(LM,1);
VS = zeros(LM,9);

UnitNames = cell(LM,1);

ii=0;
for hh=1:LM
    hh
   [MATpath, Matfile]=fileparts(Res.List_matfilepath{hh});
    for rr = 1:LR
        RMatfile = Rmatfiles(rr).name;
        if strcmp(RMatfile(9:(end-4)), Matfile(9:end))
            break
        end
    end
    MatfilePath=fullfile(Idir, Matfile);
    RMatfilePath=fullfile(Rdir, RMatfile);
    MAT=load(MatfilePath);
    if ~strcmp(MAT.subject, 'WhiBlu5396M')% No randSOund matrices for WhiBlu5396M
        ii=ii+1;
        % keep track of the cell name in the outputs
        UnitNames{ii}=Matfile(9:end);
        %load the random matrices file
        RMAT = load(RMatfilePath);
        %find the sound correlation matrices for that cell
        for cmf=1:NCMfiles
            if strcmp(sprintf('stat%s_CorrFirstVoc_%s.mat', MAT.subject, Matfile(9:13)), CorrMatfiles(cmf).name)
                CORR=load(fullfile(Ldir, CorrMatfiles(cmf).name));
                break
            end
        end
        % find the best window to construct the IV confusion matrix of that cell 
        MI_ConfCT = max(MAT.mi_confusionCT);
        if FigFlag==1
            if zs_MI_diag_uni_cat(ii)>=2
                fprintf('%d: SEMANTIC Cell! the MI_confusion of the CT Matrix is %f\n', ii, MI_ConfCT);
            else
                fprintf('%d: the MI_confusion of the CT Matrix is %f\n',ii, MI_ConfCT);
            end
        else
            fprintf('%d: the MI_confusion of the CT Matrix is %f\n',ii, MI_ConfCT);    
        end
        MAXWinCT=find(MAT.mi_confusionCT==MI_ConfCT);% we are choosing the window size that gives the best value of MI conf in the CT matrix
        Mat = MAT.confusionMatrix{MAXWinCT};

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
        % For each category calculate observed entropy then max 
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
        IS2(ii) = 1- Hobs2/Hmax2;
            
        
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
        
        %% Plot the results of calculation for that cell to review calculations
        if FigFlag==1
            %calculate index of invariance for that unit
            InvInd = 1-(mi_diag_uni_cat-mi_diag_uni_cat_maxInv)/(mi_diag_uni_cat_minInv-mi_diag_uni_cat_maxInv);
            figure(4)
            GRAD2=cubehelix(ceil(max(logpv)), 0.5, -1.5, 2, 0.150, [0,1]);
            for jj=1:length(MI.MI_tot)
                if logpv(jj)<2
                    plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6,0.6,0.6]);
                else
                    plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(ceil(logpv(jj)),:));
                end
                    hold on
            end
            plot(mi_tot, mi_diag_uni_cat/mi_tot,'ko', 'MarkerFaceColor','k');
            text(2, 0.95, sprintf('Index of Invariance: %f', InvInd));
            text(2, 0.9, sprintf('Index of Selectivity: %f or %f (norm)', IS(ii), IS2(ii)));
            
            hold off
            xlabel('Mutual Information of the call matrix (bits)')
            ylabel('Semantic Index')
            title(sprintf('Observed Matrices\nz-axis: -log10(pvalue) on Semantic Index with RandSound Matrices'));
            colorbar
            colormap(GRAD2)
            axis([0 6 0 1])
            pause
        end
        
        %% Now calculate values of MI in different areas of random and random sound matrices
        if RandSwitch==1
            rt = 1;
            while rt<=2
                if rt==1
                    LRM=length(RMAT.CM_IV_RandSound);%number of random matrices
                elseif rt==2
                    LRM=length(RMAT.CM_IV_Rand);%number of random matrices
                end
                MI_tot_rand=zeros(LRM,1);
                MI_diag_uni_rand=zeros(LRM,1);
                MI_diag_uni_cat_rand=zeros(LRM,1);
                for rm =1:LRM
                    rm
                    if rt==1
                        Rmat = RMAT.CM_IV_RandSound{rm};
                    elseif rt==2
                        Rmat = RMAT.CM_IV_Rand{rm};
                    end
                    % construct the cell array of indices for each random category
                    if rt==1
                        Nvoc = sum(CORR.NVT); % number of vocalization in the matrix
                        StimTypeRM=CORR.UVT;
                        NstimTypeRM=length(StimTypeRM); % there should be NstimTypeCM-1 categories
                        if sum(strcmp(StimTypeRM, StimTypeCM(1:NstimTypeRM)))~=NstimTypeRM
                            fprintf('WARNING StimTypes are not in the same order in random sound and real matrices')
                        end
                        List_indices= RMAT.List_VocRandSound{rm};
                        cat_rand = cell(NstimTypeCM,1);
                        for cc = 1:NstimTypeRM
                            cat_rand{cc} = sum(CORR.NVT(1:(cc-1)))+1 : sum(CORR.NVT(1:cc));
                        end
                        cat_rand{NstimTypeCM} = (Nvoc+1):size(Rmat,1);%these are indices of the background category
                    elseif rt==2
                        List_indices= RMAT.List_VocRand{rm};
                        cat_rand = cell(NstimTypeCM,1);
                        CC=0;
                        for cc = 1:NstimTypeCM
                            st = StimTypeCM(cc);
                            CC_local = sum(strcmp(VocTypeSel,st));
                            cat_rand{cc}= (CC+1):(CC+CC_local);
                            CC=CC + CC_local;
                        end
                    end      
                    [ mi_tot, mi_tot2, mi_diag, mi_error, mi_diag_uni, mi_all_error_uni, mi_diag_uni_cat, mi_real_error_uni]=info_matrix_perarea(Rmat, cat_rand);
                    MI_tot_rand(rm)=mi_tot;
                    MI_diag_uni_rand(rm)=mi_diag_uni;
                    MI_diag_uni_cat_rand(rm)=mi_diag_uni_cat;    
                end
                AvMI_tot_Rand(ii,rt) = mean(MI_tot_rand);
                SDMI_tot_Rand(ii,rt) = std(MI_tot_rand);
                fprintf('the MI_confusion of the IV Matrix is %f for that cell\nThe average MI of random matrices is %f+/-%f\n', MI_tot(ii), AvMI_tot_Rand(ii),SDMI_tot_Rand(ii));
            %pause
                AvMI_diag_uni_Rand(ii,rt)=mean(MI_diag_uni_rand);
                SDMI_diag_uni_Rand(ii,rt)=std(MI_diag_uni_rand);
                AVMI_diag_uni_cat_Rand(ii,rt)=mean(MI_diag_uni_cat_rand);
                SDMI_diag_uni_cat_Rand(ii,rt)=std(MI_diag_uni_cat_rand);
                if rt==1
                    MI_uni_diag_cat_RandSound(ii,:) = MI_diag_uni_cat_rand;
                elseif rt==2
                    MI_uni_diag_cat_Rand(ii,:) = MI_diag_uni_cat_rand;
                end
                rt = rt +1;
            end
        end
    end
    
end

MI.MI_tot=MI_tot(1:ii);
MI.MI_diag=MI_diag(1:ii);
MI.MI_error=MI_error(1:ii);
MI.MI_diag_uni=MI_diag_uni(1:ii);
MI.MI_all_error_uni=MI_all_error_uni(1:ii);
MI.MI_diag_uni_cat=MI_diag_uni_cat(1:ii);
MI.MI_real_error_uni=MI_real_error_uni(1:ii);
MI.MI_diag_uni_cat_maxInv=MI_diag_uni_cat_maxInv(1:ii);
MI.MI_diag_uni_cat_minInv=MI_diag_uni_cat_minInv(1:ii);
if RandSwitch==1
    MI.AvMI_tot_Rand = AvMI_tot_Rand(1:ii, :);
    MI.SDMI_tot_Rand = SDMI_tot_Rand(1:ii, :);
    MI.AvMI_diag_uni_Rand=AvMI_diag_uni_Rand(1:ii, :);
    MI.SDMI_diag_uni_Rand=SDMI_diag_uni_Rand(1:ii, :);
    MI.AVMI_diag_uni_cat_Rand=AVMI_diag_uni_cat_Rand(1:ii, :);
    MI.SDMI_diag_uni_cat_Rand=SDMI_diag_uni_cat_Rand(1:ii, :);
    MI.MI_uni_diag_cat_RandSound = MI_uni_diag_cat_RandSound(1:ii, :);
    MI.MI_uni_diag_cat_Rand = MI_uni_diag_cat_Rand(1:ii,:);
end
MI.UnitNames=UnitNames(1:ii);
MI.PCC_cat=PCC_cat(1:ii);
MI.IS=IS(1:ii);
MI.IS3=IS3(1:ii);
MI.IS2=IS2(1:ii);
MI.VS=VS(1:ii,:);
MI.II2 = II2(1:ii);
MIfilename='MI_analysis_140207.mat';
save(fullfile(resultsDirectory, MIfilename), '-struct', 'MI');
%January 31 2014 I'm adding the selectivity part to that code and run
%all the code again under strfinator. I'm investigating selectivity using the data
%stored under IndexSelectivity_confmat.mat in cd /Frederic/ConfMat
%(hardrive of Frederic). Found out I did not save correctly IS3 so run the
%code again on 0205 2014. Found out that I did a mistaike in the category
%definition for random matrices so correct that and run the code again on
%strfinator on 02072014

%% Load the results
%MI=load('/Users/elie/Documents/MATLAB/data/MI_analysis.mat');
MI=load(fullfile(resultsDirectory,'MI_analysis_140207.mat'));%These are the MI values obtained with the code above.
MIBG = load(fullfile(resultsDirectory, 'MI_analysis_RandBG_140207.mat'));%these are the values obtained with random matrices with Background category intact (random_Matrices_Cal.m)
%% Bunch of figures to visualize and analyse the data 
figure(1)
 subplot(1,2,1)
 plot(MI.MI_tot, MI.MI_diag_uni,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Mutual Information of diagonal in the uniform call matrix (bits)')
hold on
MAX1 = max(max(MI.MI_all_error_uni), max(MI.MI_diag_uni));
plot([0,MAX1], [0,MAX1], 'r-');
hold off
 subplot(1,2,2)
 plot(MI.MI_tot, MI.MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Mutual Information of diagonal and category in the uniform call matrix (bits)')
hold on
MAX2=max(max(MI.MI_real_error_uni), max(MI.MI_diag_uni_cat));
plot([0,MAX2], [0,MAX2], 'r-');
hold on
plot([0,0], [0,max(MI.MI_diag_uni_cat)], 'b-');
hold off

%MI calculated on random and observed IV matrices should be the same
figure(2)
plot(MI.MI_tot, MI.AvMI_tot_Rand, 'ko')



%% Now calculate the z-score for the semantic units with only rand matrices (no sound constrains). Semantic units are...
...defined as units that have more info in the diagonal and category than expected by chance
zs_MI_diag_uni_cat_Rand = (MI.MI_diag_uni_cat - MI.AVMI_diag_uni_cat_Rand(:,2))./MI.SDMI_diag_uni_cat_Rand(:,2);
logpvRand = -log10(normpdf(zs_MI_diag_uni_cat_Rand));

zs_MI_diag_uni_cat_RandBG = (MI.MI_diag_uni_cat - MIBG.AVMI_diag_uni_cat_Rand)./MIBG.SDMI_diag_uni_cat_Rand;
logpvRandBG = -log10(normpdf(zs_MI_diag_uni_cat_RandBG));

addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')
%Plot the values of index of semantic of random and random sound
%constrained matrices
figure(1)
subplot(1,3,1)
plot(MI.AVMI_diag_uni_cat_Rand(:,2), MIBG.AVMI_diag_uni_cat_Rand,'ko')
hold on
plot(MI.AVMI_diag_uni_cat_Rand(WeirdCells,2), MIBG.AVMI_diag_uni_cat_Rand(WeirdCells),'bo')
xlabel('Mutual Information in diag and cat of the Random Marices (bits)')
ylabel('Mutual Information in diag and cat of the Random Marices BG intact (bits)')
hold on
plot([0, max(MI.AVMI_diag_uni_cat_Rand(:,2))], [0, max(MIBG.AVMI_diag_uni_cat_Rand)], 'r')
subplot(1,3,2)
plot(MIBG.AVMI_diag_uni_cat_Rand, MI.AVMI_diag_uni_cat_Rand(:,1),'ko')
hold on
plot(MIBG.AVMI_diag_uni_cat_Rand(WeirdCells), MI.AVMI_diag_uni_cat_Rand(WeirdCells,1),'bo')
xlabel('Mutual Information in diag and cat of the Random Marices BG intact (bits)')
ylabel('Mutual Information in diag and cat of the Random Marices Sound constrained (bits)')
hold on
plot([0, max(MIBG.AVMI_diag_uni_cat_Rand)], [0, max(MI.AVMI_diag_uni_cat_Rand(:,1))], 'r')
subplot(1,3,3)
plot(MI.AVMI_diag_uni_cat_Rand(:,1), MI.MI_diag_uni_cat,'ko')
hold on
plot(MI.AVMI_diag_uni_cat_Rand(WeirdCells,1), MI.MI_diag_uni_cat(WeirdCells),'bo')
xlabel('Mutual Information in diag and cat of the Random Marices Sound constrained (bits)')
ylabel('Mutual Information in diag and cat of the Random Marices observed (bits)')
hold on
plot([0, max(MI.AVMI_diag_uni_cat_Rand(:,1))], [0, max(MI.MI_diag_uni_cat)], 'r')
% there is a problem with the sound constrained matrices, there are not
% different enough from random BG matrices. So we are chossing to only use
% Random and Random BG matrices.


figure(3)
subplot(1,2,1)
%GRAD2=cubehelix(ceil(max(logpv)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor','k');
        hold on
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\n'));
axis([0 6 0 1])

subplot(1,2,2)
for jj=1:length(MI.MI_tot)
    plot(MI.MI_tot(jj), MIBG.AVMI_diag_uni_cat_Rand(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor', 'r');
    hold on
    plot(MI.MI_tot(jj), MI.AVMI_diag_uni_cat_Rand(jj,2)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor', 'y');
    hold on
    plot(MI.MI_tot(jj), MI.AVMI_diag_uni_cat_Rand(jj,1)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor', 'c');
    hold on
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Semantic Index')
legend('Random BG intact', 'Rand', 'Sound Constrained', 'Location', 'NorthWest');
title(sprintf('Random Matrices\n'));
axis([0 6 0 1])

figure(4)
GRAD2=cubehelix(ceil(max(logpvRandBG)), 0.5, -1.5, 2, 0.150, [0,1]);
subplot(1,2,1)
for jj=1:length(MI.MI_tot)
    if logpvRand(jj)>=2 && zs_MI_diag_uni_cat_Rand(jj)>0 && logpvRandBG(jj)>=2 && zs_MI_diag_uni_cat_RandBG(jj)>0
        plot(MI.MI_tot(jj), (MI.MI_diag_uni_cat(jj)-MIBG.AVMI_diag_uni_cat_Rand(jj))/MI.MI_tot(jj),'ko', 'MarkerFaceColor',GRAD2(ceil(logpvRandBG(jj)),:));
        hold on
    else
        plot(MI.MI_tot(jj), (MI.MI_diag_uni_cat(jj)-MIBG.AVMI_diag_uni_cat_Rand(jj))/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6 0.6 0.6]);
    end
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the Individual Vocalization matrix (bits)')
ylabel('Non-linear Semantic Index')
title(sprintf('Observed Matrices\n'));

subplot(1,2,2)
for jj=1:length(MI.MI_tot)
    if logpvRand(jj)>=2 && zs_MI_diag_uni_cat_Rand(jj)>0 && logpvRandBG(jj)>=2 && zs_MI_diag_uni_cat_RandBG(jj)>0
        plot(MI.MI_tot(jj), (MI.MI_diag_uni_cat(jj))/MI.MI_tot(jj),'ko', 'MarkerFaceColor',GRAD2(ceil(logpvRandBG(jj)),:));
        hold on
    else
        plot(MI.MI_tot(jj), (MI.MI_diag_uni_cat(jj))/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6 0.6 0.6]);
    end
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the Individual Vocalization matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\n'));

axis([0 6 0 1])

%% Now calculate the z-score for the semantic units. Semantic units are...
...defined as units that have more info in the diagonal and category than expected by chance
zs_MI_diag_uni_cat = (MI.MI_diag_uni_cat - MI.AVMI_diag_uni_cat_Rand(:,1))./MI.SDMI_diag_uni_cat_Rand(:,1);
logpv = -log10(normpdf(zs_MI_diag_uni_cat));
addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')
% we should do a t-test instead of a z-test (use ttest.m on distrib of unique(MI_Rand_diag_cat_uni) because the number of different possible sound groups is limited values minus observed MI_diag_cat_uni 
% add the pvalue info to the previous figure 3. Here it is calculated as
% pv2
pv2 = zeros(length(MI.MI_diag_uni_cat),1);
for mmi = 1:length(MI.MI_diag_uni_cat)
    [H,pv2(mmi)] = ttest(unique(MI.MI_uni_diag_cat_RandSound(mmi,:))-MI.MI_diag_uni_cat(mmi));
end
logpv2 = -log10(pv2);

figure(3)
subplot(1,3,3)
GRAD2=cubehelix(ceil(max(logpv2)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    if logpv2(jj)<2
        plot(MI.MI_tot(jj), (MI.MI_diag_uni_cat(jj)-MI.AVMI_diag_uni_cat_Rand(jj,1))/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6,0.6,0.6]);
    else
        plot(MI.MI_tot(jj), (MI.MI_diag_uni_cat(jj)-MI.AVMI_diag_uni_cat_Rand(jj,1))/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(ceil(logpv2(jj)),:));
    end
        hold on
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\nz-axis: -log10(pvalue) on Semantic Index with RandSound Matrices'));
colorbar
colormap(GRAD2)

figure(4)
subplot(1,2,1)
GRAD2=cubehelix(ceil(max(logpv)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    if logpv(jj)<2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6,0.6,0.6]);
    else
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(ceil(logpv(jj)),:));
    end
        hold on
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\nz-axis: -log10(pvalue) on Semantic Index with RandSound Matrices'));
colorbar
colormap(GRAD2)
axis([0 6 0 1])

subplot(1,2,2)
GRAD2=cubehelix(ceil(max(logpv)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    if logpv(jj)>=2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(ceil(logpv(jj)),:));
    end
        hold on
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('Observed Matrices\nz-axis: log10(p-value) on Semantic Index with RandSound Matrices'));
colorbar
colormap(GRAD2)
axis([0 6 0 1])

%% Explore the significance of the 3 randomizations
figure(5)
subplot(1,3,1)
GRAD1=cubehelix(ceil(max(logpvRand)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    if logpvRand(jj)<2 || zs_MI_diag_uni_cat_Rand(jj)<=0
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6,0.6,0.6]);
    elseif logpvRand(jj)>=2 && zs_MI_diag_uni_cat_Rand(jj)>0
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD1(ceil(logpvRand(jj)),:));
    end
        hold on
end
xlabel('Mutual Information of the Individual vocalization matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\nz-axis: -log10(pvalue) on Semantic Index with Rand Matrices'));
colorbar
colormap(GRAD1)
axis([0 6 0 1])
hold off

subplot(1,3,2)
GRAD2=cubehelix(ceil(max(logpvRandBG)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    if logpvRandBG(jj)<2 || zs_MI_diag_uni_cat_RandBG(jj)<=0
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6,0.6,0.6]);
    elseif logpvRandBG(jj)>=2 && zs_MI_diag_uni_cat_RandBG(jj)>0
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(ceil(logpvRandBG(jj)),:));
    end
        hold on
end
xlabel('Mutual Information of the Individual vocalization matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\nz-axis: -log10(pvalue) on Semantic Index with RandBG Matrices'));
colorbar
colormap(GRAD2)
axis([0 6 0 1])
hold off

subplot(1,3,3)
GRAD3=cubehelix(ceil(max(logpv)), 0.5, -1.5, 2, 0.150, [0,1]);
for jj=1:length(MI.MI_tot)
    if logpv(jj)<2 || zs_MI_diag_uni_cat(jj)<=0
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor',[0.6,0.6,0.6]);
    elseif logpv(jj)>=2 && zs_MI_diag_uni_cat(jj)>0
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD3(ceil(logpv(jj)),:));
    end
        hold on
end
xlabel('Mutual Information of the Individual vocalization matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Observed Matrices\nz-axis: -log10(pvalue) on Semantic Index with RandSound Matrices'));
colorbar
colormap(GRAD3)
axis([0 6 0 1])
hold off

% Find the number of cell that have significantly more info than random
% matrices
Sign_Rand = (logpvRand>=2) .* (zs_MI_diag_uni_cat_Rand>0);
Sign_RandBG = (logpvRandBG>=2) .* (zs_MI_diag_uni_cat_RandBG>0);
% identify the number of cells that become significant with the random BG
% constant
WeirdCells = find((~Sign_Rand).* Sign_RandBG);

%% Investigating the weird cells
figure(8)
subplot(1,3,1)
plot(MI.AVMI_diag_uni_cat_Rand(:,2), MIBG.AVMI_diag_uni_cat_Rand,'ko')
hold on
plot(MI.AVMI_diag_uni_cat_Rand(WeirdCells,2), MIBG.AVMI_diag_uni_cat_Rand(WeirdCells),'bo')
xlabel('Mutual Information in diag and cat of the Random Marices (bits)')
ylabel('Mutual Information in diag and cat of the Random Marices BG intact (bits)')
hold on
plot([0, max(MI.AVMI_diag_uni_cat_Rand(:,2))], [0, max(MIBG.AVMI_diag_uni_cat_Rand)], 'r')
subplot(1,3,3)
plot(MIBG.AVMI_diag_uni_cat_Rand, MI.MI_diag_uni_cat,'ko')
hold on
plot(MIBG.AVMI_diag_uni_cat_Rand(WeirdCells), MI.MI_diag_uni_cat(WeirdCells),'bo')
xlabel('Mutual Information in diag and cat of the Random Marices BG intact (bits)')
ylabel('Mutual Information in diag and cat of the observed Matrices (bits)')
hold on
plot([0, max(MIBG.AVMI_diag_uni_cat_Rand)], [0, max(MI.MI_diag_uni_cat)], 'r')
subplot(1,3,2)
plot(MI.AVMI_diag_uni_cat_Rand(:,2), MI.MI_diag_uni_cat,'ko')
hold on
plot(MI.AVMI_diag_uni_cat_Rand(WeirdCells,2), MI.MI_diag_uni_cat(WeirdCells),'bo')
xlabel('Mutual Information in diag and cat of the Random Marices (bits)')
ylabel('Mutual Information in diag and cat of the Matrices observed (bits)')
hold on
plot([0, max(MI.AVMI_diag_uni_cat_Rand(:,1))], [0, max(MI.MI_diag_uni_cat)], 'r')

% Find the hh indices of these weirdCells (position in
%%Res.List_matfilepath)
hhWeirdCells = zeros(size(WeirdCells));
resultsDirectory='/Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat';
Res=load(fullfile(resultsDirectory, 'Results131130.mat'), 'List_matfilepath');
LM = length(Res.List_matfilepath);
ww=0;
for hh = 1:LM
    for ii = 1:length(hhWeirdCells)
        [path, hhUnit] = fileparts(Res.List_matfilepath{hh});
        if strcmp(hhUnit(9:end), MI.UnitNames(WeirdCells(ii)))
            hhUnit(9:end)
            MI.UnitNames(WeirdCells(ii))
            ww=ww+1
            hh
            hhWeirdCells(ww) = hh;
            break
        end
    end
end
 
% Who's that?
cd /Users/elie/Documents/SauvegardeDisqueFREDERIC/ConfMat
Idir=pwd;
for ff =1:length(hhWeirdCells)
    [MATpath, Matfile]=fileparts(Res.List_matfilepath{hhWeirdCells(ff)});
    MatfilePath=fullfile(Idir, Matfile);
    MAT=load(MatfilePath);
    MAT.subject
end

%% Now look at the invariance of these semantic cells
% Invariance is defined as the proportion of information in the category
% area compare to the information in the diagonal + category area.

% Figure 5 with all cells on the left and only significant cells on the
% right
figure(5)
subplot(1,2,1)
Invariance = (MI.MI_diag_uni_cat - MI.MI_diag_uni)./MI.MI_diag_uni_cat;
GRAD2=cubehelix(ceil(max(Invariance*10000)), 0.5, -1.5, 2, 0.500, [0,1]);

%image(im)
for jj=1:length(MI.MI_tot)
    if logpv(jj)<2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor','k');
    else
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(Invariance(jj)*10000),:));
    end
    hold on
end
 %plot(MI_tot, MI_diag_uni,'ko','MarkerEdgeColor','k',...
 %   'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('All Units\nz-axis: Invariance non-significant cells in black'));
colorbar('YTickLabel', (0.1:0.1:0.9))
colormap(GRAD2)
axis([0 6 0 1])
hold off

subplot(1,2,2)
for jj=1:length(MI.MI_tot)
    if logpv(jj)>=2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(Invariance(jj)*10000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('Semantic Units\nz-axis: Invariance Proportion of MI in the categories'));
colorbar('YTickLabel', (0.1:0.1:0.9))
colormap(GRAD2)
axis([0 6 0 1])


% Figure 6 with random matrices on the left and observed matrices on the
% right
figure(6)
subplot(1,2,1)
InvarianceRand = (MI.AVMI_diag_uni_cat_Rand - MI.AvMI_diag_uni_Rand)./MI.AVMI_diag_uni_cat_Rand;
GRAD2=cubehelix(ceil(max(InvarianceRand*10000)), 0.5, -1.5, 2, 0.500, [0,1]);
%image(im)
for jj=1:length(MI.MI_tot)
    plot(MI.MI_tot(jj), MI.AVMI_diag_uni_cat_Rand(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(InvarianceRand(jj)*10000),:));
    hold on
end
 %plot(MI_tot, MI_diag_uni,'ko','MarkerEdgeColor','k',...
 %   'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('Random-constrained Matrices\nz-axis: Invariance Proportion of MI in the categories*10000 compare to diag+cat'));
colorbar
colormap(GRAD2)
hold off
%axis([0 6 0 0.9])
% hold on
% MAX1 = max(max(MI_tot), max(MI_diag_uni));
% plot([0,MAX1], [0,MAX1], 'r-');
% hold off

subplot(1,2,2)
InvarianceAll = (MI_diag_uni_cat - MI_diag_uni)./MI_tot;
InvarianceSub = (MI.MI_diag_uni_cat - MI.MI_diag_uni)./MI.MI_diag_uni_cat;
GRAD2=cubehelix(ceil(max(InvarianceSub*10000)), 0.5, -1.5, 2, 0.500, [0,1]);
for jj=1:length(MI.MI_tot)
    plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(InvarianceSub(jj)*10000),:));
    hold on
end
%plot(MI_tot, MI_diag_uni_cat,'ko','MarkerEdgeColor','k',...
%    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('Observed Matrices\nz-axis: Invariance Proportion of MI in the categories*10000'));
colorbar
colormap(GRAD2)


%% Now explore the selectivity index defined as 1-HdiagCTMatrix/HdiagUniCTMatrix
figure(7)
subplot(1,2,1)
GRAD2=cubehelix(ceil(max(MI.IS*10000)), 0.5, -1.5, 2, 0.500, [0,1]);

%image(im)
for jj=1:length(MI.MI_tot)
    if logpv(jj)<2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor','k');
    else
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(MI.IS(jj)*10000),:));
    end
    hold on
end
 %plot(MI_tot, MI_diag_uni,'ko','MarkerEdgeColor','k',...
 %   'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('All Units\nz-axis: Selectivity non-significant cells in black'));
colorbar('YTickLabel', (0.05:0.05:0.3))
colormap(GRAD2)
axis([0 6 0 1])
hold off

subplot(1,2,2)
for jj=1:length(MI.MI_tot)
    if logpv(jj)>=2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(MI.IS(jj)*10000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('Semantic Units\nz-axis: Selectivity 1-Proportion of entropy'));
colorbar('YTickLabel', (0.05:0.05:0.3))
colormap(GRAD2)
axis([0 6 0 1])

%% Now explore the evolution of spike rate in the new MI space

%% Now look at the real error of classification of these semantic cells
% Real error is defined as the proportion of information off category
% area and off diagonal compare to the information in the diagonal + category area.

% Figure 7 with all cells on the left and only significant cells on the
% right
figure(7)
subplot(1,2,1)
RealError = (MI.MI_real_error_uni)./MI.MI_tot;
GRAD2=cubehelix(ceil(max(RealError*10000)), 0.5, -1.5, 2, 0.500, [0,1]);

%image(im)
for jj=1:length(MI.MI_tot)
    if logpv(jj)<2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj),'ko', 'MarkerFaceColor','k');
    else
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(RealError(jj)*10000),:));
    end
    hold on
end
 %plot(MI_tot, MI_diag_uni,'ko','MarkerEdgeColor','k',...
 %   'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('All Units\nz-axis: RealError non-significant cells in black'));
colorbar('YTickLabel', (0.1:0.1:0.8))
colormap(GRAD2)
axis([0 6 0 1])
hold off

subplot(1,2,2)
for jj=1:length(MI.MI_tot)
    if logpv(jj)>=2
        plot(MI.MI_tot(jj), MI.MI_diag_uni_cat(jj)/MI.MI_tot(jj), 'ko', 'MarkerFaceColor',GRAD2(round(RealError(jj)*10000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the call matrix (bits)')
ylabel('Proportion of MI in diagonal and category')
title(sprintf('Semantic Units\nz-axis: RealError Proportion of MI off categories and off diag'));
colorbar('YTickLabel', (0.1:0.1:0.8))
colormap(GRAD2)
axis([0 6 0 1])

%% Investigate selectivity of these cells
%using the data previously calculated by IndexSelectivity_confmat.m. Note
%that the present code is now (31 january 2014) also calculating the same
%stuff
cd /Volumes/Frederic/ConfMat/
IndexSel = load('IndexSelectivity_confmat.mat');
% check for cells name concordance between MI_analysis.mat file and
% IndexSelectivity_confmat.mat
% the order of the file in IndexSel is in ListMatfilepath under results
load('Results131130.mat', List_matfilepath);

