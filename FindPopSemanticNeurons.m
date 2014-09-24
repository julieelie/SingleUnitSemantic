%cd /auto/k6/julie/matfile
% resultsDirectory='/auto/k6/julie/matfile';
% addpath('/auto/k1/queued');
cd /Volumes/FREDERIC/ConfMat

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

% vectors containing optimal windows for each unit
optWin=zeros(LM,1); % based on the sound file matrix
optWin_BGwo=zeros(LM,1); % based on the sound file matrix without BG
optWin_CT= zeros(LM,1); % based on the calls type matrix

% vectors containing mutual information for each unit
MI_confB=zeros(LM,9);
MI_confBBGwo=zeros(LM,9);
MI_confBCT=zeros(LM,9);
MI_confBCTBGwo=zeros(LM,9);

% vector containing the size of the confusion matrices for each unit
ConfSize = zeros(LM,2);
ConfSizeCT = zeros(LM,2);
DiagCMCT = zeros(LM,10); %only take into account call categories (no mln) and no whines since data for only 1 bird for whines + Bg category

        for hh=1:LM
            
            MatfilePath=fullfile(Idir, matfiles(Indices(hh)).name)
            ii=ii+1;
            MAT=load(MatfilePath);
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
            
            %Collect CT Matrices Diagonals
            DiagCMCT(ii,:)=diag(MAT.confusionMatrixCT{find(Winsize==optWin_CT(ii))});

            
            
        end
%     end
% end
optWin=optWin(1:ii);
optWin_CT=optWin_CT(1:ii);
optWin_BGwo=optWin_BGwo(1:ii);


%% Plot the histograms of optimum window
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
%% Plot the values of normalized mi_conf for all cells depending on winsize
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

MIB_OptCT = zeros(ii,1);
MIBCT_OptCT = zeros(ii,1);
MIB_Opt = zeros(ii,1);
MIBCT_Opt = zeros(ii,1);
for uu = 1:ii
    MIB_OptCT(uu) = MI_confB(uu,find(WinSize==optWin_CT(uu)));
    MIB_Opt(uu) = MI_confB(uu,find(WinSize==optWin(uu)));
    MIBCT_OptCT(uu) = MI_confBCT(uu,find(WinSize==optWin_CT(uu)));
    MIBCT_Opt(uu) = MI_confBCT(uu,find(WinSize==optWin(uu)));
end
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

zs_MIB_OptCT = (MIB_OptCT-mean(MIB_OptCT))./std(MIB_OptCT);
zs_MIBCT_OptCT = (MIBCT_OptCT-mean(MIBCT_OptCT))./std(MIBCT_OptCT);
zs_Collapse = (zs_MIB_OptCT - zs_MIBCT_OptCT);
%Collapse = (MIB_OptCT - MIBCT_OptCT)./MIB_OptCT;
Collapse = (MIB_OptCT - MIBCT_OptCT)./(MIB_max-MIBCT_max);
figure(9)
subplot(1,2,1)
hist(zs_MIB_OptCT);
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
     [IDXRand,CRand,sumdRand,DRand] = kmeans(X,k);
     IDXTot(:,k)=IDX;
     sumdTot{k}=sumd;
 end
 
 % plot the sum of the the within-cluster sums of point-to-centroid
 % distances (quality of the k-mean) against increasing number of clusters
 SumD=zeros(1,Ncluster);
for k=1:Ncluster
    SumD(k) = sum(sumdTot{k});
end
figure(13)
 plot(SumD);
MI=min(SumD);
NClusters=find(SumD==MI);
vline(NClusters)

NCchoosen=[3 4 5 6 7 8 9 10];
COLOR = {'r+' 'b+' 'g+' 'c+' 'm+' 'kd' 'rd' 'k+' 'bd' 'md'};

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


