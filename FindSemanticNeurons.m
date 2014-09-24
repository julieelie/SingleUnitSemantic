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
%     fileexist = 0;
%     for lrm=1:LRM
%         RandMatF = randmatfiles(lrm).name;
%         if strcmp(Matfile(9:end), RandMatF(9:end))
%             fileexist = 1;
%             RandMatfile = RandMatF;
%         end
%     end
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

%% How to best calculate the z-score: Are SD of rand matrices larger than...
...sd of randsound matrices? No!!!
figure(5)
subplot(1,3,1)
hist(SDMI_rand)
title('SD values of MI on Random matrices')
subplot(1,3,2)
hist(SDMI_randSound)
title('SD values of MI on sound-constrained matrices')
subplot(1,3,3)
hist(SDMI_rand-SDMI_randSound)
title('SD(MI randMat) - SD(MI sound-constrained mat)')


%% Plot an exemple of correlation values between sounds
% add some empty rows and lines between vocalization categories
MATRIX=CORR.CORR;
UVT = unique(CORR.VocTypeSel);
NUVT=length(UVT);
MATRIX2 = ones(size(MATRIX,1)+NUVT-1,size(MATRIX,2)+NUVT-1)*(max(max(MATRIX)))*2;
%first find min and max indices of each vocalization types
MIN=nan(NUVT,1);
MAX=nan(NUVT,1);
for tt = 1:NUVT
    vt=UVT{tt};
    MAX(tt)=max(find(strcmp(CORR.VocTypeSel, vt)));
    MIN(tt)=min(find(strcmp(CORR.VocTypeSel, vt)));
end

% then loop through rows and column to fill in matrices
for rr = 1:NUVT
    Irow2 = (MIN(rr):MAX(rr))+rr-1;
    Irow = MIN(rr):MAX(rr);
    for cc = 1:NUVT
        Icol2 = (MIN(cc):MAX(cc))+cc-1;
        Icol = MIN(cc):MAX(cc);
        MATRIX2(Irow2, Icol2)=MATRIX(Irow, Icol);
    end
end

Clim=[0 (max(max(MATRIX)))*2];
createcorrelationmatrix2(MATRIX2,MIN,MAX,UVT, Clim) % sound file matrix with white barres
axis([0.5 size(MATRIX,1)+0.5 0.5 size(MATRIX,1)+0.5])


%% Exploring the effect of window size in the calculation of confusion matrices on MI-confusion
% Plot the histograms of optimum window
figure(1)
subplot(2,1,1)
hist(optWin, 600)
axis([0 610 0 500])
hold on
xlabel('WinSize of confusion matrix')
ylabel('nb cells')
title(sprintf('Best window based on MI sound file matrix nb of cells:%d', ii));
hold off
subplot(2,1,2)
hist(optWin_CT, 600)
axis([0 610 0 500])
hold on
xlabel('WinSize of confusion matrix')
ylabel('nb cells')
title(sprintf('Best Window based on MI Call Type Matrix nb of cells: %d',ii))
hold off

WinSize=unique(optWin);
 
%Plot the values of normalized mi_conf for all cells depending on winsize
figure(2)
CO ={'b' ; 'g' ; 'r' ; 'c' ; 'm' ; 'y' ; 'k'}; 
for jj = 1:ii
    subplot(2,1,1)
    rn =randperm(7);
    plot(MI_confB(jj,:)/max(MI_confB(jj,:)), CO{rn(1)});
    set(gca,'XTickLabel',WinSize);
    axis([1 9 0 1.1])
    hold on
    subplot(2,1,2)
    plot(MI_confBCT(jj,:)/max(MI_confBCT(jj,:)), CO{rn(1)});
    axis([1 9 0 1.1])
    set(gca,'XTickLabel',WinSize);
    hold on
    %pause
end

hold off

% Plot the values of mi_conf calculated on Call type matrix for all cells depending on winsize by chunk of
% 30 units
figure(3)
CO ={'b' ; 'g' ; 'r' ; 'c' ; 'm' ; 'y' ; 'k'}; 
UU = randperm(ii);
for vv=1:floor(ii/30)
for kk = 1:30
    jj=UU(kk+(vv-1)*30);
    subplot(2,1,1)
    rn =randperm(7);
    plot(MI_confBCT(jj,:), CO{rn(1)});
    set(gca,'XTickLabel',WinSize);
    
    hold on
    subplot(2,1,2)
    plot(MI_confBCTBGwo(jj,:), CO{rn(1)});
    
    set(gca,'XTickLabel',WinSize);
    hold on
    %pause
end
pause
subplot(2,1,1)
hold off;
subplot(2,1,2)
hold off
end

% Plot the values of mi_conf calculated on sound file matrix for all cells depending on winsize by chunk of
% 30 units
figure(4)
CO ={'b' ; 'g' ; 'r' ; 'c' ; 'm' ; 'y' ; 'k'}; 
UU = randperm(ii);
for vv=1:floor(ii/30)
for kk = 1:30
    jj=UU(kk+(vv-1)*30);
    subplot(2,1,1)
    rn =randperm(7);
    plot(MI_confB(jj,:), CO{rn(1)});
    set(gca,'XTickLabel',WinSize);
    
    hold on
    subplot(2,1,2)
    plot(MI_confBBGwo(jj,:), CO{rn(1)});
    
    set(gca,'XTickLabel',WinSize);
    hold on
    %pause
end
pause
subplot(2,1,1)
hold off;
subplot(2,1,2)
hold off
end

% Plots mean of mi confusion depending on win size for the 4 matrices
% (sound file with or without background, call type with or without
% background)
figure(5)
subplot(4,1,1)
plot(mean(MI_confB));
subplot(4,1,2)
plot(mean(MI_confBCT));
subplot(4,1,3)
plot(mean(MI_confBBGwo));
subplot(4,1,4)
plot(mean(MI_confBCTBGwo));

% plot the log2 of Mi confusion 
figure(6)
plot(log2(MI_confB'))
hold on
set(gca,'XTickLabel',WinSize);
hold off

%% choisir la fenêtre optimale mi_confBCT calculer %collapse entre mi_confB
% histogram of optimal window with mi_confCT
figure(7)
hist(optWin_CT, 600)
axis([0 610 0 400])
hold on
xlabel('WinSize of confusion matrix')
ylabel('nb cells')
title(sprintf('MiConfB of cal type matrix nb of cells:%d', ii));

% collapse
% We are looking for semantic neurons so let's choose the confusion matrix
% window that give the highest value of mutual information in the call type
% matrix

figure(8)
subplot(2,2,1)
hist(MIB_OptCT)
title('MI sound file matrix with optimal window of call type matrix')
subplot(2,2,2)
hist(MIBCT_OptCT)
title('MI call type matrix with optimal window of call type matrix')
subplot(2,2,3)
hist(MIB_Opt)
title('MI sound file matrix with optimal window of sound file matrix')
subplot(2,2,4)
hist(MIBCT_Opt)
title('MI call type matrix with optimal window of sound file matrix')
% Note that even if I was choosing the best window for the sound file
% matrix, that might not change drastically the results, the histogram look
% really similar

% semantic neurons = neurons with high values of information in the sound file matrix and neurons
% that does not collapse when going down from sound file matrix to call
% type matrix
HighMIB = find(MIB_OptCT>=1);

%% Identify the collapsing neurons and plot them under the hist of MIB_Opt
% To compare values of MI I need to zscore the values of MI

% Calculate the expected minimum decrease of MI due to the change of matrix
% size only for each unit
MIB_max = zeros(ii,1);
MIBCT_max = zeros(ii,1);

for uu=1:ii
    MIB_max(uu) = info_matrix(diag(1/ConfSize(uu,1)*ones(ConfSize(uu,1),1)));
    MIBCT_max(uu) = info_matrix(diag(1/ConfSizeCT(uu,1)*ones(ConfSizeCT(uu,1),1)));
end


%Collapse = (MIB_OptCT - MIBCT_OptCT)./MIB_OptCT;
Collapse = (MIB_OptCT - MIBCT_OptCT)./(MIB_max-MIBCT_max);

%% Plot the MI values obtained for all units and the one expected for random...
... arrangement of sound files

figure(9)
subplot(1,2,1)
plot(MIB_OptCT, MIBCT_OptCT,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold on
plot(MIB_OptCT, MeanMI_rand, 'r+')
hold on
plot(MIB_OptCT, MeanMI_randSound, 'b+')
hold off
xlabel('Mutual Information Sound File matrix (bits)')
ylabel('Mutual Information Call Type matrix (bits)')
legend('Call Type grouping', 'Random grouping', 'Random grouping controlled for acoustic features correlation', 'Location', 'NorthWest');
legend boxoff

subplot(1,2,2)
plot(MIB_OptCT_pSpike, MIBCT_OptCT_pSpike,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold on
plot(MIB_OptCT_pSpike, MeanMI_rand_pSpike, 'r+')
hold on
plot(MIB_OptCT_pSpike, MeanMI_randSound_pSpike, 'b+')
hold off
xlabel('Mutual Information Sound File matrix (bits)')
ylabel('Mutual Information Call Type matrix (bits)')
legend('Call Type grouping', 'Random grouping', 'Random grouping controlled for acoustic features correlation', 'Location', 'NorthWest');
legend boxoff
axis([0 3 0 0.14])

figure(1)
subplot(1,2,1)
GRAD2=cubehelix(ceil(max(MeanSR)), 0.5, -1.5, 2, 0.425, [0,1]);
%image(im)
for jj=1:length(MIB_OptCT)
    jj
    if MeanSR(jj)<1
        plot(MIB_OptCT_pSpike(jj), MIBCT_OptCT_pSpike(jj), 'ko', 'MarkerFaceColor','k');
    else
        plot(MIB_OptCT_pSpike(jj), MIBCT_OptCT_pSpike(jj), 'ko', 'MarkerFaceColor',GRAD2(round(MeanSR(jj)),:));
    end
        hold on
end
hold off
colorbar
colormap(GRAD2)
axis([0 3 0 0.14])
xlabel('Mutual Information Sound File matrix (bits/spike)')
ylabel('Mutual Information Call Type matrix (bits/spike)')
title('Spike rate color coded, black dots with SR<1')
subplot(1,2,2)
GRAD2=cubehelix(ceil(max(MeanSR)), 0.5, -1.5, 2, 0.425, [0,1]);
%image(im)
for jj=1:length(MIB_OptCT)
    jj
    if MeanSR(jj)<1
        plot(MIB_OptCT(jj), MIBCT_OptCT(jj), 'ko', 'MarkerFaceColor','k');
    else
        plot(MIB_OptCT(jj), MIBCT_OptCT(jj), 'ko', 'MarkerFaceColor',GRAD2(round(MeanSR(jj)),:));
    end
        hold on
end
hold off
colorbar
colormap(GRAD2)
xlabel('Mutual Information of individual call matrix (bits)')
ylabel('Mutual Information Call Type matrix (bits)')
title('Spike rate color coded, black dots with SR<1')


figure(2)
plot(MIB_OptCT_pSpike, MIBCT_OptCT_pSpike,'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
xlabel('Mutual Information of individual call matrix (bits/spike)')
ylabel('Mutual Information of semantic group matrix (bits/spike)')
axis([0 3 0 0.14])

figure(3)
plot(MIB_OptCT_pSpike, MeanMI_rand_pSpike, 'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.9 0.7 0.20], 'MarkerSize', 6)
hold on
plot(MIB_OptCT_pSpike, MeanMI_randSound_pSpike, 'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.2 0.7 0.9], 'MarkerSize', 6)
axis([0 3 0 0.14])
legend('Random groups', 'Sound-constrained groups', 'Location', 'NorthWest');
legend boxoff
xlabel('Mutual Information of individual call matrix (bits/spike)')
ylabel('Mutual Information of semantic group matrix (bits/spike)')
[FORand, GRand]=fit(MIB_OptCT_pSpike, MeanMI_rand_pSpike, 'poly1');
hold on
PH=refline(FORand.p1, FORand.p2)
set(PH, 'Color', [0.9 0.7 0.20])
hold on
[FORandS, GRandS]=fit(MIB_OptCT_pSpike, MeanMI_randSound_pSpike, 'poly1');
PH=refline(FORandS.p1, FORandS.p2)
set(PH, 'Color', [0.2 0.7 0.9])

figure(4)
plot(MIB_OptCT_pSpike, MeanMI_randSound_pSpike-MeanMI_rand_pSpike, 'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.9 0.7 0.20], 'MarkerSize', 6)
xlabel('Mutual Information of individual call matrix (bits/spike)')
ylabel('Difference of Mutual Information between control matrices SC-Rand (bits/spike)')
axis([0 3 min(MeanMI_randSound_pSpike-MeanMI_rand_pSpike) max(MeanMI_randSound_pSpike-MeanMI_rand_pSpike)])
% proportion of unit that show higher mean values of sound constrained MI
% compare to random MI
sum((MeanMI_randSound_pSpike-MeanMI_rand_pSpike)>=0)./length(MeanMI_randSound_pSpike-MeanMI_rand_pSpike)

figure(5)
plot(MIB_OptCT_pSpike, SDMI_randSound_pSpike-SDMI_rand_pSpike, 'ko','MarkerEdgeColor','k',...
    'MarkerFaceColor', [0.9 0.7 0.20], 'MarkerSize', 6)
xlabel('Mutual Information of individual call matrix (bits/spike)')
ylabel('Difference of SD Mutual Information between control matrices SC-Rand (bits/spike)')
axis([0 3 min(SDMI_randSound_pSpike-SDMI_rand_pSpike) max(SDMI_randSound_pSpike-SDMI_rand_pSpike)])
% proportion of unit that show higher mean values of sound constrained SD MI
% compare to random MI
sum((SDMI_randSound_pSpike-SDMI_rand_pSpike)>=0)./length(SDMI_randSound_pSpike-SDMI_rand_pSpike)



figure(6)
GRAD2=cubehelix(ceil(max(logpv_MIBCT_OptCTSound)), 0.5, -1.5, 2, 0.150, [0,1]);
%image(im)
for jj=1:length(MIB_OptCT)
    jj
    if logpv_MIBCT_OptCTSound(jj)<2
        plot(MIB_OptCT_pSpike(jj), MIBCT_OptCT_pSpike(jj), 'ko', 'MarkerFaceColor','k');
    else
        plot(MIB_OptCT_pSpike(jj), MIBCT_OptCT_pSpike(jj), 'ko', 'MarkerFaceColor',GRAD2(round(logpv_MIBCT_OptCTSound(jj)),:));
    end
        hold on
end
axis([0 3 0 0.14])
PH=refline(FORand.p1, FORand.p2)
set(PH, 'Color', [0.9 0.7 0.20])
hold on
PH=refline(FORandS.p1, FORandS.p2)
set(PH, 'Color', [0.2 0.7 0.9])
hold off
colorbar
colormap(GRAD2)
xlabel('Mutual Information of individual call matrix (bits/spike)')
ylabel('Mutual Information of semantic group matrix (bits/spike)')


figure(7)
SemanticNeurons = find(logpv_MIBCT_OptCTSound>2);
NonSemanticNeurons = find(logpv_MIBCT_OptCTSound<=2);
plot(zs_MIBCT_OptCTSound(NonSemanticNeurons),MIBCT_OptCT_pSpike(NonSemanticNeurons),'ko', 'MarkerFaceColor','k');
hold on
plot(zs_MIBCT_OptCTSound(SemanticNeurons),MIBCT_OptCT_pSpike(SemanticNeurons),'ko', 'MarkerFaceColor','r');
legend('Non-semantic units', 'Semantic units', 'Location', 'NorthWest');
legend boxoff
xlabel('Z-scored value of MI of semantic group matrix (bits/spike)')
ylabel('Mutual Information of semantic group matrix (bits/spike)')

figure(10)
subplot(1,2,1)
hist(zs_MIBCT_OptCT, 70);
title('z-score Rand of the mutual information calculated on the call type confusion matrix');
vline(2,'g')
vline(5)
vline(8)
vline(20)
sum(zs_MIBCT_OptCT>=3)
sum(zs_MIBCT_OptCT>=5)
subplot(1,2,2)
hist(zs_MIBCT_OptCTSound, 70);
title('z-score RandSound of the mutual information calculated on the call type confusion matrix');
vline(2,'g')
vline(4)
vline(9.5)
vline(14)

%% find significant semantic units
IndSemU = find(zs_MIBCT_OptCT>=2);
IndSemU2_8 = intersect(IndSemU, find(zs_MIBCT_OptCT<8));
IndSemU8_20 = intersect(find(zs_MIBCT_OptCT>=8), find(zs_MIBCT_OptCT<20));
IndSemU20_40 = intersect(find(zs_MIBCT_OptCT>=20), find(zs_MIBCT_OptCT<40));
IndSemU40 = find(zs_MIBCT_OptCT>=40);
IndNonSemU = find(zs_MIBCT_OptCT<2);
IndSemUSound = find(zs_MIBCT_OptCTSound>=2);
IndSemU2_4Sound = intersect(IndSemUSound, find(zs_MIBCT_OptCTSound<4));
IndSemU4_95Sound = intersect(find(zs_MIBCT_OptCTSound>=4), find(zs_MIBCT_OptCTSound<9.5));
IndSemU95_14Sound = intersect(find(zs_MIBCT_OptCTSound>=9.5), find(zs_MIBCT_OptCTSound<14));
IndSemU14Sound = find(zs_MIBCT_OptCTSound>=14);
IndNonSemUSound = find(zs_MIBCT_OptCTSound<2);
figure(11)
subplot(1,2,1)
plot(MIB_OptCT(IndNonSemU), MIBCT_OptCT(IndNonSemU), 'k+')
hold on
plot(MIB_OptCT(IndSemU2_8), MIBCT_OptCT(IndSemU2_8), 'b+')
hold on
plot(MIB_OptCT(IndSemU8_20), MIBCT_OptCT(IndSemU8_20), 'c+')
hold on
plot(MIB_OptCT(IndSemU20_40), MIBCT_OptCT(IndSemU20_40), 'g+')
hold on
plot(MIB_OptCT(IndSemU40), MIBCT_OptCT(IndSemU40), 'r+')
hold off
title('Random z-score')

xlabel('Sound File matrix Mutual Information')
ylabel('Call Type matrix Mutual Information')
legend('z-score<2', '2<z-score<8', '8<z-score<20', '20<z-score<40','40<z-score', 'Location', 'NorthWest');

subplot(1,2,2)
figure(12)
plot(MIB_OptCT(IndNonSemUSound), MIBCT_OptCT(IndNonSemUSound), 'k+')
hold on
plot(MIB_OptCT(IndSemU2_4Sound), MIBCT_OptCT(IndSemU2_4Sound), 'b+')
hold on
plot(MIB_OptCT(IndSemU4_95Sound), MIBCT_OptCT(IndSemU4_95Sound), 'c+')
hold on
plot(MIB_OptCT(IndSemU95_14Sound), MIBCT_OptCT(IndSemU95_14Sound), 'g+')
hold on
plot(MIB_OptCT(IndSemU14Sound), MIBCT_OptCT(IndSemU14Sound), 'r+')
hold off
title('Sound z-score')

xlabel('Sound File matrix Mutual Information')
ylabel('Call Type matrix Mutual Information')
legend('z-score<2', '2<z-score<4', '4<z-score<9.5', '9.5<z-score<14','14<z-score', 'Location', 'NorthWest');



%% Super MI Cloud with z-score color graded
figure(13)
%[grad,im]=colorGradient([0 1 0],[1 0 0],ceil((max(zs_MIBCT_OptCT)-min(zs_MIBCT_OptCT))*100+1));
GRAD = cubehelix((max(zs_MIBCT_OptCT)-min(zs_MIBCT_OptCT))*100+2,0, 0.460, 3, 1.281, [0.061 1]);
GRAD2=cubehelix((max(zs_MIBCT_OptCT)-min(zs_MIBCT_OptCT))*100+2);
%image(im)
for jj=1:length(MIB_OptCT)
    jj
    plot(MIB_OptCT(jj), MIBCT_OptCT(jj), 'ko', 'MarkerFaceColor',GRAD2(round((zs_MIBCT_OptCT(jj)-min(zs_MIBCT_OptCT))*100)+1,:));
    hold on
end
hold off
colorbar
colormap(GRAD2)



figure(12)
subplot(1,2,2)
hist(zs_MIBCT_OptCT);
figure(10)
subplot(2,2,1)
hist(zs_Collapse) 
subplot(2,2,2)
hist(zs_Collapse(HighMIB)) 
subplot(2,2,3)
hist(Collapse) 
subplot(2,2,4)
hist(Collapse(HighMIB))
%axis([0.2 1 0 800])
% decide of a threshold of collapse to distinguish semantic neurons from
% non-semantic ones???!!!!

% to estimate the drop of info based on shrinkage claculate the mean info
% expected when bootstratping on

%% Figure of the confusion matrix for the different cells in each class of z-score
%MinUnits=IndNonSemU;
%MinUnits=IndSemU5_8;
%MinUnits=IndSemU8_20;
MinUnits=IndSemU20;

for mm=1:length(MinUnits)
MinUnit = MinUnits(mm);
MatfilePath=fullfile(Idir, matfiles(Indices(MinUnit)).name);
MAT=load(MatfilePath);
optWin_CT(MinUnit)

Iwin=find(Winsize==optWin_CT(MinUnit));

figure(15)
imagesc(MAT.confusionMatrix{Iwin});
colorbar;
xlabel('Model Vocalization');
ylabel('Actual Vocalization');
title(sprintf('All vocalizations\n Winsize=%d miconf=%0.2f PCC=%0.2f z_score=%0.2f',optWin_CT(MinUnit), MAT.mi_confusionB(Iwin), MAT.percorrectB(Iwin), zs_MIBCT_OptCT(MinUnit)));
UVT=unique(MAT.VocTypeSel);
TICK=zeros(length(UVT), 1);
TICKLABEL=cell(size(TICK));
for tt = 1:length(UVT)
    vt=UVT{tt};
    MAX=max(find(strcmp(MAT.VocTypeSel, vt)));
    MIN=min(find(strcmp(MAT.VocTypeSel, vt)));
    vline(MAX)
    hline(MAX)
    TICK(tt)=(MAX+MIN)/2;
    TICKLABEL{tt}=vt;
    text(TICK(tt), 145, vt)
end


figure(16);
imagesc(MAT.confusionMatrixCT{Iwin});
colorbar;
xlabel('Model vocalization');
ylabel('Actual vocalization');
title(sprintf('Vocalization Categories p(x,y)\n Winsize=%d miconf=%0.2f PCC=%0.2f', optWin_CT(MinUnit), MAT.mi_confusionCT(Iwin), MAT.percorrectCT(Iwin)));
set(gca(), 'Ytick', 1:length(MAT.StimType));
set(gca(), 'YTickLabel', MAT.StimType);
set(gca(), 'Xtick', 1:length(MAT.StimType));
set(gca(), 'XTickLabel', MAT.StimType);
pause
end

%% Figures of the confusion matrix of 
%the cell with the highest value of mutual information on the sound file matrix
MaxUnit = find(MIB_OptCT==max(MIB_OptCT))
FiguresConfusionMatrices(List_matfilepath, MaxUnit, optWin_CT);

%the cell with the mean value of mutual information on the sound file matrix
MaxUnit = intersect(find(MIB_OptCT>(mean(MIB_OptCT)+1.5)), find(MIB_OptCT<(mean(MIB_OptCT)+1.51)))
FiguresConfusionMatrices(List_matfilepath, MaxUnit(1), optWin_CT);

%the cell with the highest value of mutual information on the call type matrix
MaxUnit = find(MIB_OptCT==max(MIB_OptCT))
FiguresConfusionMatrices(List_matfilepath, MaxUnit, optWin_CT);

%the cell with the highest value of mutual information/spike on the call type matrix
MaxUnit = find(MIBCT_OptCT_pSpike==max(MIBCT_OptCT_pSpike))
FiguresConfusionMatrices(List_matfilepath, MaxUnit, optWin_CT);

% the cell with the highest value of zscored mutual information random
% sound groups
MaxUnit = find(zs_MIBCT_OptCT==max(zs_MIBCT_OptCT))
FiguresConfusionMatrices(List_matfilepath, MaxUnit, optWin_CT);

% the cell with the highest value of z-scored MI acoutically supervised
% sound groups
MaxUnit = find(zs_MIBCT_OptCTSound==max(zs_MIBCT_OptCTSound))
FiguresConfusionMatrices(List_matfilepath, MaxUnit, optWin_CT);

%% figure of the confusion matrix of the cell with the magenta MI confusion
MinUnits = intersect(intersect(find(MIBCT_OptCT>0.1), find(MIBCT_OptCT<0.2)), intersect(find(MIB_OptCT>2.8), find(MIB_OptCT<3)));
for mm=1:length(MinUnits)
MinUnit = MinUnits(mm);
MatfilePath=fullfile(Idir, matfiles(Indices(MinUnit)).name);
MAT=load(MatfilePath);
optWin_CT(MinUnit)
Iwin=find(Winsize==optWin_CT(MinUnit));

figure(15)
imagesc(MAT.confusionMatrix{Iwin});
colorbar;
xlabel('Model Vocalization');
ylabel('Actual Vocalization');
title(sprintf('All vocalizations\n Winsize=%d miconf=%0.2f PCC=%0.2f',optWin_CT(MinUnit), MAT.mi_confusionB(Iwin), MAT.percorrectB(Iwin)));
UVT=unique(MAT.VocTypeSel);
TICK=zeros(length(UVT), 1);
TICKLABEL=cell(size(TICK));
for tt = 1:length(UVT)
    vt=UVT{tt};
    MAX=max(find(strcmp(MAT.VocTypeSel, vt)));
    MIN=min(find(strcmp(MAT.VocTypeSel, vt)));
    vline(MAX)
    hline(MAX)
    TICK(tt)=(MAX+MIN)/2;
    TICKLABEL{tt}=vt;
    text(TICK(tt), 145, vt)
end


figure(16);
imagesc(MAT.confusionMatrixCT{Iwin});
colorbar;
xlabel('Model vocalization');
ylabel('Actual vocalization');
title(sprintf('Vocalization Categories p(x,y)\n Winsize=%d miconf=%0.2f PCC=%0.2f', optWin_CT(MinUnit), MAT.mi_confusionCT(Iwin), MAT.percorrectCT(Iwin)));
set(gca(), 'Ytick', 1:length(MAT.StimType));
set(gca(), 'YTickLabel', MAT.StimType);
set(gca(), 'Xtick', 1:length(MAT.StimType));
set(gca(), 'XTickLabel', MAT.StimType);
pause
end



%% figure of the confusion matrix of the cell with the low MI confusion
MinUnits = intersect(intersect(find(MIBCT_OptCT>0.08), find(MIBCT_OptCT<0.15)), intersect(find(MIB_OptCT>0.5), find(MIB_OptCT<1)));
MinUnits = intersect(find(MIB_OptCT>2), find(MIB_OptCT<3));
MinIndex = find(MIBCT_OptCT(MinUnits)==min(MIBCT_OptCT(MinUnits)));
MinUnit = MinUnits(MinIndex);

%for mm=1:length(MinUnits)
%MinUnit = MinUnits(mm);
MatfilePath=fullfile(Idir, matfiles(Indices(MinUnit)).name);
MAT=load(MatfilePath);
optWin_CT(MinUnit)
Iwin=find(Winsize==optWin_CT(MinUnit));

figure(15)
imagesc(MAT.confusionMatrix{Iwin});
colorbar;
xlabel('Model Vocalization');
ylabel('Actual Vocalization');
title(sprintf('All vocalizations\n Winsize=%d miconf=%0.2f PCC=%0.2f',optWin_CT(MinUnit), MAT.mi_confusionB(Iwin), MAT.percorrectB(Iwin)));
UVT=unique(MAT.VocTypeSel);
TICK=zeros(length(UVT), 1);
TICKLABEL=cell(size(TICK));
for tt = 1:length(UVT)
    vt=UVT{tt};
    MAX=max(find(strcmp(MAT.VocTypeSel, vt)));
    MIN=min(find(strcmp(MAT.VocTypeSel, vt)));
    vline(MAX)
    hline(MAX)
    TICK(tt)=(MAX+MIN)/2;
    TICKLABEL{tt}=vt;
    text(TICK(tt), 145, vt)
end


figure(16);
imagesc(MAT.confusionMatrixCT{Iwin});
colorbar;
xlabel('Model vocalization');
ylabel('Actual vocalization');
title(sprintf('Vocalization Categories p(x,y)\n Winsize=%d miconf=%0.2f PCC=%0.2f', optWin_CT(MinUnit), MAT.mi_confusionCT(Iwin), MAT.percorrectCT(Iwin)));
set(gca(), 'Ytick', 1:length(MAT.StimType));
set(gca(), 'YTickLabel', MAT.StimType);
set(gca(), 'Xtick', 1:length(MAT.StimType));
set(gca(), 'XTickLabel', MAT.StimType);
%pause
%end

%% Figure of the cell cloud with size of confusion window color coded
Winsize
IndWinU2_10 = find(optWin_CT<=10);
IndWinU20_50 = intersect(find(optWin_CT>10), find(optWin_CT<100));
IndWinU100_300 = intersect(find(optWin_CT>=100), find(optWin_CT<600));
IndWinU600 = find(optWin_CT == 600);

figure(17)
plot(MIB_OptCT(IndWinU2_10), MIBCT_OptCT(IndWinU2_10), 'k+')
hold on
plot(MIB_OptCT(IndWinU20_50), MIBCT_OptCT(IndWinU20_50), 'c+')
hold on
plot(MIB_OptCT(IndWinU100_300), MIBCT_OptCT(IndWinU100_300), 'g+')
hold on
plot(MIB_OptCT(IndWinU600), MIBCT_OptCT(IndWinU600), 'r+')
hold off
xlabel('Mutual Information Sound file marix')
ylabel('Mutual information Call Type matrix')


%% figure of the cells cloud with type of spike color coded
figure(18)
IndSS_Mult=find(Spike_shape==0);
IndSS_Unclass=find(Spike_shape==1);
IndSS_Narrow=find(Spike_shape==3);
IndSS_Large=find(Spike_shape==2);
IndSS_Wide=find(Spike_shape==4);

plot(MIB_OptCT(IndSS_Mult), MIBCT_OptCT(IndSS_Mult), 'ko')
hold on
plot(MIB_OptCT(IndSS_Unclass), MIBCT_OptCT(IndSS_Unclass), 'k+')
hold on
plot(MIB_OptCT(IndSS_Narrow), MIBCT_OptCT(IndSS_Narrow), 'c+')
hold on
plot(MIB_OptCT(IndSS_Large), MIBCT_OptCT(IndSS_Large), 'g+')
hold on
plot(MIB_OptCT(IndSS_Wide), MIBCT_OptCT(IndSS_Wide), 'r+')
hold off
xlabel('Mutual Information Sound file marix')
ylabel('Mutual information Call Type matrix')
legend('MultiU', 'Unclassified', 'Narrow', 'Large', 'Wide', 'Location','NorthWest');

%% figure of the cells cloud with the mean spike rate color coded or 3D
figure(19)
subplot(1,2,1)
plot3(MIB_OptCT, MIBCT_OptCT,MeanSR, 'k+');
subplot(1,2,2)
hist(MeanSR, 30)

figure(20)
IndSRU0_10 = find(MeanSR<=10);
IndSRU10_50 = intersect(find(MeanSR>10), find(MeanSR<=50));
IndSRU50_100 = intersect(find(MeanSR>=50), find(MeanSR<100));
IndSRU100 = find(MeanSR> 100);
plot(MIB_OptCT(IndSRU0_10), MIBCT_OptCT(IndSRU0_10), 'k+')
hold on
plot(MIB_OptCT(IndSRU10_50), MIBCT_OptCT(IndSRU10_50), 'c+')
hold on
plot(MIB_OptCT(IndSRU50_100), MIBCT_OptCT(IndSRU50_100), 'g+')
hold on
plot(MIB_OptCT(IndSRU100), MIBCT_OptCT(IndSRU100), 'r+')
hold off
xlabel('Mutual Information Sound file marix (bits/s)')
ylabel('Mutual information Call Type matrix (bits/s)')
legend('MeanSR<10', '10<MeanSR<50', '50<MeanSR<100', 'MeanSR>100', 'Location', 'NorthWest');


figure(21)
plot(MIB_OptCT_pSpike(IndSRU0_10), MIBCT_OptCT_pSpike(IndSRU0_10), 'k+')
hold on
plot(MIB_OptCT_pSpike(IndSRU10_50), MIBCT_OptCT_pSpike(IndSRU10_50), 'c+')
hold on
plot(MIB_OptCT_pSpike(IndSRU50_100), MIBCT_OptCT_pSpike(IndSRU50_100), 'g+')
hold on
plot(MIB_OptCT_pSpike(IndSRU100), MIBCT_OptCT_pSpike(IndSRU100), 'r+')
hold off
xlabel('Mutual Information Sound file marix (bits/spike)')
ylabel('Mutual information Call Type matrix (bits/spike)')
legend('MeanSR<10', '10<MeanSR<50', '50<MeanSR<100', 'MeanSR>100', 'Location', 'NorthEast');
axis([0 3 0 0.15])

%% figure of the cloud in bits per spike with z-score value color coded
figure(30)
GRAD = cubehelix((max(zs_MIBCT_OptCT)-min(zs_MIBCT_OptCT))*100+2,0, 0.460, 3, 1.281, [0.061 1]);
GRAD2=cubehelix((max(zs_MIBCT_OptCT)-min(zs_MIBCT_OptCT))*100+2);
%image(im)
for jj=1:length(MIB_OptCT_pSpike)
    jj
    plot(MIB_OptCT_pSpike(jj), MIBCT_OptCT_pSpike(jj), 'ko', 'MarkerFaceColor',GRAD2(round((zs_MIBCT_OptCT(jj)-min(zs_MIBCT_OptCT))*100)+1,:));
    hold on
end
axis([0 3 0 0.14])
hold off
colorbar
colormap(GRAD2)

%% figure of the cells cloud in bits/spike with type of spike color coded
figure(23)
plot(MIB_OptCT_pSpike(IndSS_Mult), MIBCT_OptCT_pSpike(IndSS_Mult), 'ko')
hold on
plot(MIB_OptCT_pSpike(IndSS_Unclass), MIBCT_OptCT_pSpike(IndSS_Unclass), 'k+')
hold on
plot(MIB_OptCT_pSpike(IndSS_Narrow), MIBCT_OptCT_pSpike(IndSS_Narrow), 'c+')
hold on
plot(MIB_OptCT_pSpike(IndSS_Large), MIBCT_OptCT_pSpike(IndSS_Large), 'g+')
hold on
plot(MIB_OptCT_pSpike(IndSS_Wide), MIBCT_OptCT_pSpike(IndSS_Wide), 'r+')
hold off
xlabel('Mutual Information Sound file marix')
ylabel('Mutual information Call Type matrix')
legend('MultiU', 'Unclassified', 'Narrow', 'Large', 'Wide', 'Location','NorthEast');
axis([0 3 0 0.15])

%% figure of the cells cloud with the subject color coded
figure(24)
SUBJS = unique(SUBJ(:,1));
NSUBJ=length(SUBJS);
COLOR = get(gca,'ColorOrder');
for ii = 1:NSUBJ
    ss=SUBJS(ii)
    plot(MIB_OptCT_pSpike(find(SUBJ(:,1)==ss)), MIBCT_OptCT_pSpike(find(SUBJ(:,1)==ss)),'+', 'Color', COLOR(ii,:))
    hold on
end
legend(SUBJECTS(SUBJS+1), 'Location', 'NorthWest')
xlabel('Mutual Information Sound file marix')
ylabel('Mutual information Call Type matrix')
hold off

%% figure of the cells cloud with the anatomical zone color coded
figure(25)
ZO=unique(ZONES(:,1));
NZO = length(ZO);
GRAD2=cubehelix(ceil(max(logpv_MIBCT_OptCTSound)), 0.5, -1.5, 2, 0.150, [0,1]);
for ii = 1:NZO
    subplot(2,5,ii)
    ss=ZO(ii)
    Local_indices=find(ZONES(:,1)==ss);
    for jj = 1:length(Local_indices)
        if logpv_MIBCT_OptCTSound(Local_indices(jj))<2
            plot(MIB_OptCT_pSpike(Local_indices(jj)), MIBCT_OptCT_pSpike(Local_indices(jj)), 'ko', 'MarkerFaceColor','k');
        else
            plot(MIB_OptCT_pSpike(Local_indices(jj)), MIBCT_OptCT_pSpike(Local_indices(jj)), 'ko', 'MarkerFaceColor',GRAD2(round(logpv_MIBCT_OptCTSound(Local_indices(jj))),:))
        end
        hold on
    end
    axis([0 3 0 0.14])
    axis([0 3 0 0.14])
    PH=refline(FORand.p1, FORand.p2)
    set(PH, 'Color', [0.9 0.7 0.20])
    hold on
    PH=refline(FORandS.p1, FORandS.p2)
    set(PH, 'Color', [0.2 0.7 0.9])
    if ss~=0
        title(ZONES_List(ss))
    else
        title('Unknown')
    end
    hold on
    xlabel('MI of individual calls matrix')
    ylabel('MI of semantic groups matrix')
    hold off
    pause 
end
%legend(['Unknown', ZONES_List], 'Location', 'NorthWest')

figure(26)
ZO=unique(ZONES(:,1));
NZO = length(ZO);
GRAD2=cubehelix(ceil(max(logpv_MIBCT_OptCTSound)), 0.5, -1.5, 2, 0.150, [0,1]);
pp=1;
for ii = 1:NZO
    subplot(2,3,pp)
    ss=ZO(ii)
    Local_indices=find(ZONES(:,1)==ss);
    for jj = 1:length(Local_indices)
        if logpv_MIBCT_OptCTSound(Local_indices(jj))<2
            plot(MIB_OptCT_pSpike(Local_indices(jj)), MIBCT_OptCT_pSpike(Local_indices(jj)), 'ko', 'MarkerFaceColor','k');
        else
            plot(MIB_OptCT_pSpike(Local_indices(jj)), MIBCT_OptCT_pSpike(Local_indices(jj)), 'ko', 'MarkerFaceColor',GRAD2(round(logpv_MIBCT_OptCTSound(Local_indices(jj))),:))
        end
        hold on
    end
        axis([0 3 0 0.14])
        if ss~=0
            title(ZONES_List(ss))
        else
            title('Unknown')
        end
        hold on
         xlabel('MI of individual calls matrix')
        ylabel('MI of semantic groups matrix')
       
        pause
        if ii<4
            hold off
            pp=pp+1;
            axis([0 3 0 0.14])
            PH=refline(FORand.p1, FORand.p2)
            set(PH, 'Color', [0.9 0.7 0.20])
            hold on
            PH=refline(FORandS.p1, FORandS.p2)
            set(PH, 'Color', [0.2 0.7 0.9])
        elseif ii>7
            hold off
            pp=pp+1;
            axis([0 3 0 0.14])
            PH=refline(FORand.p1, FORand.p2)
            set(PH, 'Color', [0.9 0.7 0.20])
            hold on
            PH=refline(FORandS.p1, FORandS.p2)
            set(PH, 'Color', [0.2 0.7 0.9])
        end
    
end
hold off

%% Figure per zone for MI semantic and semantic z-score
figure(28)
ZO=unique(ZONES(:,1));
NZO = length(ZO);
pp=1;
for ii = 1:NZO
    subplot(2,3,pp)
    ss=ZO(ii)
    Local_indices=find(ZONES(:,1)==ss);
    Local_indicesSem  = intersect(Local_indices, SemanticNeurons);
    Local_indicesnonSem  = intersect(Local_indices, NonSemanticNeurons);
    plot(zs_MIBCT_OptCTSound(Local_indicesnonSem),MIBCT_OptCT_pSpike(Local_indicesnonSem),'ko', 'MarkerFaceColor','k');
    hold on
    plot(zs_MIBCT_OptCTSound(Local_indicesSem),MIBCT_OptCT_pSpike(Local_indicesSem),'ko', 'MarkerFaceColor','r');
    hold on
    legend('Non-semantic units', 'Semantic units', 'Location', 'NorthWest');
    legend boxoff
    xlabel('Z-score(MI) of semantic matrix (bits/spike)')
    ylabel('MI of semantic matrix (bits/spike)')
    if ss~=0
            title(ZONES_List(ss))
    else
            title('Unknown')
    end
    if ii<4
            hold off
            pp=pp+1;
            axis([-5 25 0 0.14])
    elseif ii>7
            hold off
            pp=pp+1;
            axis([-5 25 0 0.14])
    end
end
%% Figure per zone with spike shape
figure(27)
ZO=unique(ZONES(:,1));
NZO = length(ZO);
COLOR1 = cubehelix(length(unique(Spike_shape)),0, 0.460, 3, 1.281, [0.061 1]);
COLOR2 = cubehelix(length(unique(Spike_shape)));
pp=1;
for ii = 1:NZO
    subplot(2,3,pp)
    ss=ZO(ii)
    Local_indices=find(ZONES(:,1)==ss);
    for jj = 1:length(Local_indices)
        %subplot(2,5,ii)
        jj
        plot(MIB_OptCT_pSpike(Local_indices(jj)), MIBCT_OptCT_pSpike(Local_indices(jj)), 'ko', 'MarkerFaceColor',COLOR1(Spike_shape(Local_indices(jj))+1,:))
        hold on
    end
    legend('MultiU', 'Unclassified', 'Narrow', 'Large', 'Wide', 'Location','NorthEast');
    axis([0 3 0 0.14])
        if ss~=0
            title(ZONES_List(ss))
        else
            title('Unknown')
        end
        hold on
        xlabel('Mutual Information Sound file matrix')
        ylabel('Mutual Information Call Type matrix')
       
        pause
        if ii<4
            hold off
            pp=pp+1;
        elseif ii>7
            hold off
            pp=pp+1;
        end
    
end
hold off


%% %% Kmean on diagonal of matrices: finding population of neurons with the same confusion matrix profil
 fprintf(1, 'Calculate PC on CT Matrices Diagonals\n');
x=DiagCMCT; %only take into account call categories (no mln) and no whines since data for only 1 bird for whines + Bg category
[COEFF,SCORE,latent,tsquare]=princomp(x,'econ');

% find optimal number of PC
figure(11)
 plot(cumsum(latent/sum(latent)))

 % scatter plot of the data in the PC space
 figure(12)
 plot3(SCORE(:,1), SCORE(:,2), SCORE(:,3), '+')
 
 % I choose to work with the first 4 PC
 X = SCORE(:,1:4);
 Nunits=size(DiagCMCT, 1);
 Ncluster=200;
 %empty cluster when reaching k=289 below
 IDXTot=zeros(Nunits, Ncluster);
 sumdTot = cell(Ncluster,1);
 sumdTotRand = cell(Ncluster,1);
 % also generate random values to test the relationship between the
 % increasing cluster number and the diminution of distance
 MU=ones(Nunits, size(X,2));
 SIGMA=ones(Nunits, size(X,2));
 MeanX=mean(X);
 STDX=std(X);
 for ss = 1:size(X,2)
    MU(:,ss)= ones(Nunits,1)*MeanX(ss);
    SIGMA(:,ss)= ones(Nunits,1)*STDX(ss);
 end
    
 XRand = normrnd(MU, SIGMA);
 
 %Calculate K-means with increasing value of k (nb of groups)
 for k=1:Ncluster
     k
     [IDX,C,sumd,D] = kmeans(X,k);
     [IDXRand,CRand,sumdRand,DRand] = kmeans(XRand,k);
     IDXTot(:,k)=IDX;
     sumdTot{k}=sumd;
     sumdTotRand{k}=sumdRand;
 end
 
 % plot the sum of the the within-cluster sums of point-to-centroid
 % distances (quality of the k-mean) against increasing number of clusters
 SumD=zeros(1,Ncluster);
 SumDRand=zeros(1,Ncluster);
for k=1:Ncluster
    SumD(k) = sum(sumdTot{k});
    SumDRand(k) = sum(sumdTotRand{k});
end
figure(13)
 plot(SumD);
 hold on
 plot(SumDRand, 'r');
MI=min(SumD);
NClusters=find(SumD==MI);
vline(NClusters)
hold off
axis([0 60 0 10])

NCchoosen=[3 4 5 6 7 8 9 10];
COLOR = {'r+' 'b+' 'g+' 'c+' 'm+' 'kd' 'rd' 'k+' 'bd' 'md'};

figure(15)
plot(SumD-SumDRand)
vline(4)
axis([0 60 -1 0.4])

NCC = length(NCchoosen);
figure(14)
for ncc=1:NCC
    ncluster = NCchoosen(ncc);
    IDXKmean=IDXTot(:,ncluster);
    subplot(2,4,ncc)

    for ii=1:ncluster
        plot3(SCORE(IDXKmean==ii,1), SCORE(IDXKmean==ii,2), SCORE(IDXKmean==ii,3), COLOR{ii});
        hold on
    end
    hold off
end


