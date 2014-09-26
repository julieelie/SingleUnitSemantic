cd /auto/k6/julie/matfile
resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')

cd /auto/k8/julie/RandMat
randdir=pwd;
randmatfiles = dir(fullfile(randdir, '*.mat'));
LRM = length(randmatfiles);

%cd /Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/SpikeShape
cd /auto/fhome/julie/matlab/tlab/src/h5analysis/Julie_neuralcode/SpikeShape
Ldir=pwd;
SpikeShape=load(strcat(Ldir,'/spikeShapeResults.mat'),'allNames', 'indG5', 'spikeType');
%cd /Users/elie/Documents/MATLAB/data/matfile
cd /auto/k6/julie/matfile
Ldir=pwd;
load(strcat(Ldir,'/AllUnits_LRI_DP.mat'), 'MAvgRcallcat', 'Unitnames');

%cd /Users/elie/Documents/MATLAB/data/matfile/StimCorrelationMatrices
cd /auto/k6/julie/matfile/StimCorrelationMatrices
Ldir=pwd;
CorrMatfiles=dir(Ldir);
NCMfiles = length(CorrMatfiles);

cd /auto/fhome/julie/Documents/Histo
% cd /Volumes/FREDERIC/HistoJulie/Stacks
AnaLdir=pwd;
AnatTxt=dir(fullfile(AnaLdir, 'List_h5files*Histo.txt'));


%cd /Volumes/FREDERIC/ConfMat
cd /auto/k6/julie/matfile/ConfMat
input_dir=pwd;

ii=0;

Idir=pwd;
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

% Retrieve anatomical data
Anat_root = cell(LM, 5);
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

% choose your kind of randomization, sound oriented or totally random
RandChoice = 0;%input('what kind of randomization?\n1:sound oriented\n2:totally random\n');

%this is to plot the histogram of correlation values of vocalizations selected in each category
plotfig = 0;                                 

% cell containing the path of each unit
List_matfilepath = cell(LM, 1);

% cell containing the anatomy info of each unit
List_anat = cell(LM, 5);

% vectors containing optimal windows for each unit
optWin=zeros(LM,1); % based on the sound file matrix
optWin_BGwo=zeros(LM,1); % based on the sound file matrix without BG
optWin_CT= zeros(LM,1); % based on the calls type matrix

% vectors containing mutual information for each unit
MI_confB=zeros(LM,9);
MI_confBBGwo=zeros(LM,9);
MI_confBCT=zeros(LM,9);
MI_confBCTBGwo=zeros(LM,9);

% vector containing the mutual information of the best matrix for each unit
MIB_OptCT = zeros(LM,1);
MIBCT_OptCT = zeros(LM,1);
MIB_Opt = zeros(LM,1);
MIBCT_Opt = zeros(LM,1);

% vector containing the size of the confusion matrices for each unit
ConfSize = zeros(LM,2);
ConfSizeCT = zeros(LM,2);

% vector containing the diagonal of the best call type confusion matrix
DiagCMCT = zeros(LM,10); %only take into account call categories (no mln) and no whines since data for only 1 bird for whines + Bg category

% Number of time you repeat the bootstrap
Bootmax=100;

%vector containing the mean and sd mutual information expected for randomly
%shrinked call type matrices
MeanMI_rand = zeros(LM,1);
SDMI_rand = zeros(LM,1);
MeanMI_randSound = zeros(LM,1);
SDMI_randSound = zeros(LM,1);
MIVal_randSound = zeros(LM,Bootmax);



% Vector of the mean spiking rate for each unit calculated on the mean spike rate per
% call category
MeanSR = zeros(LM,1);

% Vector of the spike shape and sort quality
Spike_shape = zeros(LM,1);

% Vector of the code for the subject
SUBJ = zeros(LM,2);
SUBJECTS = {'YelBlu6903F' 'BlaBro09xxF' 'LblBlu2028M' 'GreBlu9508M' 'WhiBlu5396M' 'WhiWhi4522M'};

%Vector of the code for the anatomical zones
ZONES = zeros(LM,2);
ZONES_List = {'CMM', 'CML', 'L1', 'L2A', 'L2B', 'L3', 'L', 'NCM', 'HP'};

for hh=1:LM
    Matfile = matfiles(Indices(hh)).name
    fileexist = 0;
    for lrm=1:LRM
        RandMatF = randmatfiles(lrm).name;
        if strcmp(Matfile(9:end), RandMatF(9:end))
            fileexist = 1;
            RandMatfile = RandMatF;
        end
    end
    if fileexist==1 %the following lines correct an error that was made in the previous code to be removed!!!!
        RandMatRes = load(fullfile(randdir, RandMatfile));
        
        GoodCM_CT_Rand = RandMatRes.CM_CT_RandSound;
        GoodCM_CT_Sound = RandMatRes.CM_CT_Rand;
        RandMatRes.CM_CT_Rand = GoodCM_CT_Rand;
        RandMatRes.CM_CT_RandSound = GoodCM_CT_Sound;

         calfilename=fullfile('/auto','k8','julie','RandMat',['RandMat_' Matfile(9:end)]);
         save(calfilename, '-struct', 'RandMatRes');
    else
            
        MatfilePath=fullfile(Idir, Matfile);
        MAT=load(MatfilePath);
        %if ~strcmp(MAT.subject, 'WhiBlu5396M')
        ii=ii+1;
        
        List_matfilepath{ii}=MatfilePath;
        
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
        if strcmp(List_anat{ii,5}, 'L')
            ZONES(ii,2)=1;
        end
        
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
        
        
        
        %Collect CT Matrices Diagonals
        DiagCMCT(ii,:)=diag(MAT.confusionMatrixCT{find(Winsize==optWin_CT(ii))});
        
        % find mean expected values of MI_conf random
        VocConfMat = MAT.confusionMatrix{find(Winsize==optWin_CT(ii))};
        VocTypeConfMat = MAT.confusionMatrixCT{find(Winsize==optWin_CT(ii))};
        SF_ConfMat = VocConfMat.*MAT.neventMatrix(find(Winsize==optWin_CT(ii)));
        VocTypeSel = MAT.VocTypeSel;
        StimTypeCM=unique(VocTypeSel);
        NstimTypeCM=length(StimTypeCM);
        mi_confCT_localSound=zeros(Bootmax, 1);
        mi_confCT_localRand=zeros(Bootmax, 1);
        
        if plotfig>0
            figure(1)
            subplot(1,2,1)
            imagesc(VocConfMat)
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
        %% find the vector of indices for spike trains obtained during sound (suppress those obtain during background) in the confusion matrix SF_ConfMat
        IndVocOnly = setdiff(1:length(VocTypeSel), find(strcmp(VocTypeSel, 'BG')));
        CM_CT_Rand =  cell(Bootmax, 1);
        CM_CT_RandSound =  cell(Bootmax, 1);
        CM_IV_Rand = cell(Bootmax, 1);
        CM_IV_RandSound = cell(Bootmax, 1);
        List_VocRand = cell(Bootmax, 1);
        List_VocRandSound=cell(Bootmax, 1);
        Mean_Corr_UVT=cell(Bootmax, 1);
        Median_Corr_UVT=cell(Bootmax, 1);
        Std_Corr_UVT=cell(Bootmax, 1);
        Min_Corr_UVT=cell(Bootmax, 1);
        Max_Corr_UVT=cell(Bootmax, 1);
        for bb=1:Bootmax
            confusion_matrix_CallTypeSound = zeros(NstimTypeCM, NstimTypeCM);
            confusion_matrix_CallTypeRand = zeros(NstimTypeCM, NstimTypeCM);
            
            rng('shuffle'); % seeds the random number generator based on the current time
            VocTypeSel_rand=VocTypeSel(randperm(length(VocTypeSel)));
            %if RandChoice==2 || RandChoice==0
            RR=0;
            VocRand = zeros(1,length(VocTypeSel));
            for vtR=1:NstimTypeCM
                stR=StimTypeCM(vtR);
                selectorR=strcmp(VocTypeSel_rand, stR);
                selectorIndR=find(selectorR);
                VocRand(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
                RR=RR+length(selectorIndR);
                
                for vtC = 1:NstimTypeCM
                    stC=StimTypeCM(vtC);
                    selectorC=strcmp(VocTypeSel_rand, stC);
                    selectorIndC=find(selectorC);
                    
                    confusion_matrix_CallTypeRand(vtR,vtC)=sum(sum(SF_ConfMat(selectorIndR, selectorIndC)));
                    
                    
                    
                end
                
            end
            confusion_matrix_vocalizationsRand = SF_ConfMat(VocRand, VocRand);
            
            %elseif RandChoice==1  || RandChoice==0
            if SUBJ(ii,1)~=4 %no pseudo-random matrices for WhiBlu5396M because the code can't find any good combnaisons of sounds
                RR=0;
                CC=0;
                VocRandSound = zeros(1,length(VocTypeSel));
                for cmf=1:NCMfiles
                    if strcmp(sprintf('stat%s_CorrFirstVoc_%s.mat', MAT.subject, Matfile(9:13)), CorrMatfiles(cmf).name)
                        CORR=load(fullfile(Ldir, CorrMatfiles(cmf).name));
                    end
                end
                CORR_Mat=CORR.CORR;
                Ncat=length(CORR.NVT);
                if Ncat+1~=NstimTypeCM
                    sprintf('WARNING: the number of categories in the single vocalization matrix is %d\nthe number of sound categories in the sound correlation matrix is %d\n it should be %d\n', NstimTypeCM, Ncat, NstimTypeCM-1);
                    return
                end
                Nvoc=size(SF_ConfMat,1);
                Ncorr=size(CORR_Mat,1);
                allcatok=1;
                while allcatok
                    fprintf(1, 'Finding random groups of vocalizations\n');
                    All_Indices = nan(Ncorr,1);
                    Mean_UVT=zeros(1,Ncat);
                    Median_UVT=zeros(1,Ncat);
                    Std_UVT=zeros(1,Ncat);
                    Min_UVT=zeros(1,Ncat);
                    Max_UVT=zeros(1,Ncat);
                    for ct = 1:Ncat
                        findallindices=1;
                        nbtrials=0;
                        while findallindices
                            nbtrials=nbtrials+1;
                            if nbtrials==101
                                % don't allow more than 100 trials
                                fprintf(1, '100 trials reached\ncannot constitute a group for category #%d\nstarting all over again\n', ct);
                                break
                            end
                            selectorIndC = nan(CORR.NVT(ct),1);
                            % randomly choose a seed for each category that is not in the
                            % indices already choosen
                            DispoVoc=setdiff(1:Ncorr,All_Indices);
                            randVoc=randperm(length(DispoVoc));
                            Seed = DispoVoc(randVoc(1));
                            selectorIndC(1) = Seed;
                            
                            % find the min and max values of correlation for that
                            % category
                            min_corr = CORR.Min_UVT(ct);
                            max_corr = CORR.Max_UVT(ct);
                            median_corr = CORR.Median_UVT(ct);
                            
                            % find all the vocalization that have correlation values for
                            % that seed in the interval for that category and randomly
                            % choose one, then find another one that have correlation
                            % values in the interval with the two previous one and so on
                            % until you reach the total number of stims for that category
                            for nvt = 1:(CORR.NVT(ct)-1) % -1 because we have already choosen 1
                                nind = sum(~isnan(selectorIndC));
                                PossCorr = [];
                                for NI = 1:nind
                                    % select all the vocalization indices that have
                                    % correlation values in the limits for that
                                    % vocalization NI
                                    PossCorr_local = find(CORR_Mat(selectorIndC(NI),:)>min_corr & CORR_Mat(selectorIndC(NI), :)<max_corr);
                                    % Make sure that these indices have not been
                                    % already choosen for that category
                                    PossCorr_local = setdiff(PossCorr_local, selectorIndC);
                                    % Make sure that these indices have not been
                                    % already choosen for previous
                                    % categories
                                    PossCorr_local = intersect(PossCorr_local, DispoVoc);
                                    
                                    % find the common possible indices between the
                                    % previous calculated ones and that one
                                    if isempty(PossCorr)
                                        PossCorr = PossCorr_local;
                                    else
                                        PossCorr = intersect(PossCorr, PossCorr_local);
                                    end
                                end
                                
                                if isempty(PossCorr)
                                    fprintf(1,'Not a single vocalization\nwithin the boundaries of correlation\nstarting category #%d with a new seed\n', ct);
                                    break
                                else
                                    % choose randomly one of the indices that works
                                    randCorr = randperm(length(PossCorr));
                                    selectorIndC(nvt+1) = PossCorr(randCorr(1));
                                end
                            end
                            if sum(isnan(selectorIndC))==0
                                findallindices=0;
                                fprintf(1, 'all the indices found for category #%d\n', ct);
                                
                                NInd = length(selectorIndC);
                                Corr_values=nan((NInd*NInd - NInd)/2, 1);
                                cor=0;
                                for rr=1:NInd-1
                                    IIr = selectorIndC(rr);
                                    for cc = rr+1:NInd
                                        IIc = selectorIndC(cc);
                                        cor=cor+1;
                                        Corr_values(cor)=CORR_Mat(IIr,IIc);
                                    end
                                end
                                Mean_UVT(ct)=nanmean(Corr_values);
                                Median_UVT(ct)=nanmedian(Corr_values);
                                Std_UVT(ct)=nanstd(Corr_values);
                                if ~isempty(Corr_values) % Corr_values is empty when there is only one stim in the vocalization category
                                    Min_UVT(ct)=min(Corr_values);
                                    Max_UVT(ct)=max(Corr_values);
                                end
                                
                                if plotfig>1 && ~isempty(PossCorr) %this is to plot the histogram of correlation values of vocalizations selected in each category
                                    
                                    figure(3)
                                    hist(Corr_values)
                                    oldaxis=axis;
                                    axis([0 1 0 oldaxis(4)])
                                    hold on
                                    vline(Mean_UVT(ct), 'r')
                                    vline(Median_UVT(ct), 'k')
                                    vline(Min_UVT(ct), 'g')
                                    vline(Max_UVT(ct), 'g')
                                    title(sprintf('Correlation values between voc of cat %s \n',ct));
                                    hold off
                                    pause
                                end
                            end
                        end
                        if nbtrials==101
                            % don't allow more than 100 trials for the same
                            % category start again from fresh
                            break
                        end
                        % record the indices selected for that category
                        All_Indices(sum(CORR.NVT(1:(ct-1)))+1 : sum(CORR.NVT(1:ct))) = selectorIndC;
                    end
                    if sum(isnan(All_Indices))==0
                        allcatok=0;
                        fprintf(1, 'all random indices found for matrix %d\n', bb);
                    end
                end
                
                % Now that we have all the indices for the different
                % categories, compile the matrix accordingly
                
                % first do the equivalence between indices from the
                % correlation matrix with indices from the confusion matrix
                % for each element of each category
                if  sum(strcmp(MAT.TDTwav(1:Ncorr), CORR.TDTwav))==Ncorr
                    fprintf('ça va d''aller!\n')
                else
                    fprintf('problem of correspondence of indices write the looop to check!! ln 245\n')
                    return
                end
                
                BackIndices = setdiff(1:Nvoc, All_Indices);
                for vtR=1:(NstimTypeCM)
                    if vtR<NstimTypeCM
                        selectorIndR= All_Indices(sum(CORR.NVT(1:(vtR-1)))+1 : sum(CORR.NVT(1:vtR)));
                    else
                        selectorIndR= BackIndices;
                    end
                    VocRandSound(RR+1:RR+length(selectorIndR)) = selectorIndR;
                    RR=RR+length(selectorIndR);
                    for vtC = 1:(NstimTypeCM)
                        if vtC<NstimTypeCM
                            selectorIndC=All_Indices(sum(CORR.NVT(1:(vtC-1)))+1 : sum(CORR.NVT(1:vtC)));
                        else
                            selectorIndC=BackIndices;
                        end
                        confusion_matrix_CallTypeSound(vtR,vtC)=sum(sum(SF_ConfMat(selectorIndR, selectorIndC)));
                        
                        
                    end
                end
                confusion_matrix_vocalizationsSound = SF_ConfMat(VocRandSound, VocRandSound);
                confusion_matrix_CallTypeSound = confusion_matrix_CallTypeSound./sum(sum(confusion_matrix_CallTypeSound));
                confusion_matrix_vocalizationsSound = confusion_matrix_vocalizationsSound./sum(sum(confusion_matrix_vocalizationsSound));
                mi_confCT_localSound(bb)=info_matrix(confusion_matrix_CallTypeSound);
                Mean_Corr_UVT{bb}=Mean_UVT;
                Median_Corr_UVT{bb}=Median_UVT;
                Std_Corr_UVT{bb}=Std_UVT;
                Min_Corr_UVT{bb}=Min_UVT;
                Max_Corr_UVT{bb}=Max_UVT;
                CM_CT_RandSound{bb} =  confusion_matrix_CallTypeSound;
                CM_IV_RandSound{bb} = confusion_matrix_vocalizationsSound;
                List_VocRandSound{bb}=VocRandSound;
            end
            
            confusion_matrix_CallTypeRand = confusion_matrix_CallTypeRand./sum(sum(confusion_matrix_CallTypeRand));
            
            confusion_matrix_vocalizationsRand = confusion_matrix_vocalizationsRand./sum(sum(confusion_matrix_vocalizationsRand));
            mi_confCT_localRand(bb)=info_matrix(confusion_matrix_CallTypeRand);
            CM_CT_Rand{bb} =  confusion_matrix_CallTypeRand;
            CM_IV_Rand{bb} = confusion_matrix_vocalizationsRand;
            List_VocRand{bb} = VocRand;
            
            
            if plotfig>0
                figure(2)
                subplot(2,2,1)
                imagesc(confusion_matrix_vocalizationsSound)
                colorbar;
                xlabel('Predicted vocalization');
                ylabel('Actual vocalization');
                title(sprintf('Single vocalizations Sound'));
                subplot(2,2,2)
                imagesc(confusion_matrix_vocalizationsRand)
                colorbar;
                xlabel('Predicted vocalization');
                ylabel('Actual vocalization');
                title(sprintf('Single vocalizations Rand'));
                subplot(2,2,3)
                imagesc(confusion_matrix_CallTypeSound)
                colorbar;
                xlabel('Predicted vocalization');
                ylabel('Actual vocalization');
                title(sprintf('Vocalization Categories miconf=%0.2f Sound', mi_confCT_localSound(bb)));
                subplot(2,2,4)
                imagesc(confusion_matrix_CallTypeRand)
                colorbar;
                xlabel('Predicted vocalization');
                ylabel('Actual vocalization');
                title(sprintf('Vocalization Categories miconf=%0.2f Rand', mi_confCT_localRand(bb)));
                hold off
                pause(0.5)
            end
            
        end
        MeanMI_rand(ii)=mean(mi_confCT_localRand);
        SDMI_rand(ii)=std(mi_confCT_localRand);
        MeanMI_randSound(ii)=mean(mi_confCT_localSound);
        SDMI_randSound(ii)=std(mi_confCT_localSound);
        MIVal_randSound(ii,:)=mi_confCT_localSound;
        
        % storing the random matrices
        RandMat.subject = MAT.subject;
        RandMat.originalfile = MatfilePath;
        RandMat.CM_CT_Rand = CM_CT_Rand;
        RandMat.CM_CT_RandSound = CM_CT_RandSound;
        RandMat.CM_IV_Rand = CM_IV_Rand;
        RandMat.CM_IV_RandSound = CM_IV_RandSound;
        RandMat.List_VocRand = List_VocRand;
        RandMat.List_VocRandSound = List_VocRandSound;
        RandMat.Mean_Corr_UVT=Mean_Corr_UVT;
        RandMat.Median_Corr_UVT=Median_Corr_UVT;
        RandMat.Std_Corr_UVT=Std_Corr_UVT;
        RandMat.Min_Corr_UVT=Min_Corr_UVT;
        RandMat.Max_Corr_UVT=Max_Corr_UVT;
        
        calfilename=fullfile('/auto','k8','julie','RandMat',['RandMat_' Matfile(9:end)]);
        save(calfilename, '-struct', 'RandMat');
        
        % find the mean spike rate of the unit
        uu=1;
        while strcmp(Matfile(9:end), Unitnames{uu}(7:end))==0
            uu=uu+1;
        end
        MeanSR(ii)=nanmean(MAvgRcallcat(uu,:));
        
        % find the spike shape of the unit
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
        
% hist(MIVal_randSound(ii,:));
% vline(MIBCT_OptCT(ii));
% zscore = (MIBCT_OptCT(ii) - mean(MIVal_randSound(ii,:)))./std(MIVal_randSound(ii,:))
% pvalue = log10(normpdf((MIBCT_OptCT(ii) - mean(MIVal_randSound(ii,:)))./std(MIVal_randSound(ii,:))))
% title(sprintf('log10 pvalue=%d, z-score=%f\n', pvalue, zscore));

end
List_matfilepath=List_matfilepath(1:ii);
List_anat=List_anat(1:ii,:);
SUBJ = SUBJ(1:ii,:);
ZONES = ZONES(1:ii,:);

MeanMI_rand=MeanMI_rand(1:ii);
optWin=optWin(1:ii);
optWin_CT=optWin_CT(1:ii);
optWin_BGwo=optWin_BGwo(1:ii);

MI_confB=MI_confB(1:ii);
MI_confBBGwo=MI_confBBGwo(1:ii);
MI_confBCT=MI_confBCT(1:ii);
MI_confBCTBGwo=MI_confBCTBGwo(1:ii);
ConfSize = ConfSize(1:ii);
ConfSizeCT = ConfSizeCT(1:ii);
MIB_OptCT = MIB_OptCT(1:ii);
MIB_Opt = MIB_Opt(1:ii);
MIBCT_OptCT = MIBCT_OptCT(1:ii);
MIBCT_Opt = MIBCT_Opt(1:ii);
DiagCMCT=DiagCMCT(1:ii);
MeanMI_rand=MeanMI_rand(1:ii);
SDMI_rand=SDMI_rand(1:ii);
MeanMI_randSound=MeanMI_randSound(1:ii);
SDMI_randSound=SDMI_randSound(1:ii);
MeanSR=MeanSR(1:ii);
Spike_shape = Spike_shape(1:ii);
MIVal_randSound = MIVal_randSound(1:ii, :);

% retrieve 

%% calculate the zscore of the actual value of MI confusion in the CT matrix
% compare to what is expected from a random shrinkage of the sound file
% matrix
zs_MIBCT_OptCT = (MIBCT_OptCT - MeanMI_rand)./SDMI_rand;
zs_MIBCT_OptCTSound = (MIBCT_OptCT - MeanMI_randSound)./SDMI_randSound;
zs_MIBCT_OptCTmix = (MIBCT_OptCT - MeanMI_randSound)./SDMI_rand;
pv_MIBCT_OptCTSound = normpdf(zs_MIBCT_OptCTSound);
logpv_MIBCT_OptCTSound = -log10(pv_MIBCT_OptCTSound);
pv_MIBCT_OptCTmix = normpdf(zs_MIBCT_OptCTmix);
logpv_MIBCT_OptCTmix = -log10(pv_MIBCT_OptCTmix);

MIB_OptCT_pSpike = MIB_OptCT ./ MeanSR;
MIBCT_OptCT_pSpike = MIBCT_OptCT ./ MeanSR;
MeanMI_rand_pSpike = MeanMI_rand./MeanSR;
MeanMI_randSound_pSpike= MeanMI_randSound./MeanSR;
SDMI_randSound_pSpike = SDMI_randSound./MeanSR;
SDMI_rand_pSpike = SDMI_rand./MeanSR;
SemanticNeurons = find(logpv_MIBCT_OptCTSound>2);
NonSemanticNeurons = find(logpv_MIBCT_OptCTSound<=2);



save(fullfile(resultsDirectory, 'Results131130.mat'))