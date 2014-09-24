% This script calculate for each unit and each category the invariance =
% the entropy of each bloc in the bloc diagonal

resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')


%% Retrieve some output variables
% cell containing the path of each unit
load('/auto/k8/julie/SemanticAnalysis_ShuffPredict_PCCAbCh.mat');
LM = length(List_matfilepath);
NbCell_estimate = LM;

%% Get ready some output variables
% Matrix containing the invariance value for each cell and each category
H_invariance = nan(NbCell_estimate,length(StimTypeCM));
H_invariancemax = H_invariance;
H_invariancemin = H_invariance;

%% Start to loop through files and extract information
ii = 0;

for hh = 1:LM
    [Idir, Matfile, ext] = fileparts(List_matfilepath{hh});
    
    % retrieve Confusion matrix files 
    MAT = load(List_matfilepath{hh});
    sprintf(1,'loaded %s\n', Matfile);
    ii=ii+1;
    Winsize=MAT.winSize;
    MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
    Mat = MAT.confusionMatrix{MAXWinCT};% we are choosing the window size that gives the best value of MI conf in the CT matrix
    
    
    % construct the cell array of indices for each call category
    VocTypeSel = MAT.VocTypeSel;
    StimTypeCM=unique(VocTypeSel);
    IBG = find(strcmp(StimTypeCM, 'BG'));
    StimTypeCM = [StimTypeCM(1:(IBG-1)); StimTypeCM((IBG+1):end); StimTypeCM(IBG)];
    NstimTypeCM=length(StimTypeCM);
    cat = cell(NstimTypeCM,1);
    for cc = 1:NstimTypeCM
        cat{cc}=find(strcmp(StimTypeCM(cc), VocTypeSel));
        Nb_stim = length(cat{cc});
        Psum=zeros(Nb_stim,1);
        for rr=1:Nb_stim
            rro = cat{cc}(rr);
            for co = 1:Nb_stim
                cco = cat{cc}(co);
                Psum(rr)=Psum(rr) + Mat(rro,cco);
            end
        end
        Mat_for_entropy = reshape(Mat(cat{cc}, cat{cc}),Nb_stim.*Nb_stim,1);
        Mat_for_entropy = Mat_for_entropy./sum(Psum);%some values are going to be NAN but we don't worry too much since these are non discriminant units
        MZeroIndices=find(Mat_for_entropy==0);
        Mat_for_entropy(MZeroIndices)=1;%to make sure entropy is zero when p=0 because 1.*log2(1)=0
        Psum_for_entropy = Psum./sum(Psum); % This is to calculate the minimum entropy possible given de correct classification values
        PZeroIndices=find(Psum_for_entropy==0);
        Psum_for_entropy(PZeroIndices)=1;
        H_invariance(ii,cc) = sum(-Mat_for_entropy.*log2(Mat_for_entropy));
        H_invariancemax(ii,cc) = -log2(1./(Nb_stim.*Nb_stim));
        H_invariancemin(ii,cc) = sum(-Psum_for_entropy.*log2(Psum_for_entropy));
    end
end

%% Compile results into structures
sprintf('Compile results\n')
Invariance.GlobalInvariance = Invariance;
Invariance.PerCat = H_invariance(1:ii,:);
Invariance.PerCatmax = H_invariancemax(1:ii,:);
Invariance.PerCatmin = H_invariancemin(1:ii,:);


%% Save data
save(fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance.mat'))
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance.mat'))



