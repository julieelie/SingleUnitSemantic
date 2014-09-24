% This script extract the number of vocalizations per category for each
% unit, calculate the mean spike rate per category for each unit, retrieve
% the number of trials for each vocalization for each unit

resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')


%% Retrieve some output variables
% cell containing the path of each unit
load('/auto/k8/julie/SemanticAnalysis_ShuffPredict.mat');
LM = length(List_matfilepath);
NbCell_estimate = LM;

%% Get ready some output variables
% Vectors containing the number of vocalizations per category, the total number of vocalizations and the proba of each category for all units
NVocPerCat = zeros(NbCell_estimate,length(StimTypeCM));
NVoc = zeros(NbCell_estimate,1);
ProbaVocPerCat = NVocPerCat;
VocCatNames = cell(NbCell_estimate,1);

%Vectors containing the mean spike rate per category per unit
SpikeRate_Cat = zeros(NbCell_estimate,length(StimTypeCM)-1); %We don't need the spike rate during background that's why: StimTypeCM-1

%Matrix containing the number of trial per vocalization per unit
NbTrial_voc = zeros(NbCell_estimate,length(StimTypeCM));

%Cell arrays containing the mean spike rate and the TDT name for each vocalization for each
%unit
SpikeRate_Voc = cell(NbCell_estimate,1);
TDT_names = cell(NbCell_estimate,1);

%% Start to loop through files and extract information
ii = 0;

for hh = 1:LM
    [Idir, Matfile, ext] = fileparts(List_matfilepath{hh});
    
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
    
    % retrieve FirstVoc files to extract the number
    % of trials for all stims within same category
    FVmatfiles = dir(fullfile(Idir, 'FirstVoc*.mat'));
   fm=length(FVmatfiles);
    FSS_Ind=zeros(fm,1);
    for ff=1:fm
        if ~isempty(strfind(FVmatfiles(ff).name, 'ss'))
            FSS_Ind(ff)=1;
        end
    end
    IndicesFV=find(FSS_Ind);
    FV=length(IndicesFV);
    
    
    % load the matfile 
    MAT = load(List_matfilepath{hh});
    fprintf('Loaded file\n%s\n',List_matfilepath{hh} );
    ii=ii+1;
    
    NVoc(ii) = size(MAT.VocTypeSel,1);
    
    StimTypeCM=unique(MAT.VocTypeSel);
    IBG = find(strcmp(StimTypeCM, 'BG'));
    StimTypeCM = [StimTypeCM(1:(IBG-1)); StimTypeCM((IBG+1):end); StimTypeCM(IBG)];
    NstimTypeCM=length(StimTypeCM);
    VocCatNames{ii} = StimTypeCM;
    
    for cc=1:NstimTypeCM
        cat = StimTypeCM(cc);
        NVocPerCat(ii,cc) = sum(strcmp(MAT.VocTypeSel, cat));
        ProbaVocPerCat(ii,cc) = NVocPerCat(ii,cc)./NVoc(ii);
    end
    
    %% find the mean spike rate of the unit per call category
    sprintf('Category Mean spike rate\n')
    for ww=1:WV
        WVMatfile = WVmatfiles(IndicesWV(ww)).name;
        if strcmp(Matfile(9:end), WVMatfile(10:end-4))
            WMAT = load(fullfile(Idir, WVMatfile));
            break
        end
    end
    FirstWholeVoc = find(WMAT.Voc_orders==1);
    VocTypeFirstSect = WMAT.VocType(FirstWholeVoc);
    MeanSR_sections=cell2mat(WMAT.MeanRate(FirstWholeVoc));
    Icat_all = cell(NstimTypeCM-1, 1);
    for cc=1:(NstimTypeCM-1)%There is no value of mean spike rate for backgroung here
        cat = StimTypeCM(cc);
        Icat_all{cc} = find(strcmp(VocTypeFirstSect,cat));
        SpikeRate_Cat(ii,cc) = nanmean(MeanSR_sections(Icat_all{cc}));
    end
    
    %% find the mean spike rate for each first vocalization
    sprintf('Individual Voc Mean spike rate\n')
    SpikeRate_Voc{ii} = MeanSR_sections;
    TDT_names{ii} = WMAT.TDT_wavfiles(FirstWholeVoc);

    %% Find the sum of trials for all vocalizations per category
    sprintf('Nb Trials for all vocalizations within each category\n')
    for fw=1:FV
        FVMatfile = FVmatfiles(IndicesFV(fw)).name;
        if strcmp(Matfile(9:end), FVMatfile(10:end-4))
            FMAT = load(fullfile(Idir, FVMatfile));
            break
        end
    end
    % retrieve values for non background sounds
    for cc=1:(NstimTypeCM-1)
        cat = StimTypeCM(cc);
        Icat_local = find(strcmp(FMAT.VocType,cat));
        Nsections=length(Icat_local);
        for ss = 1:Nsections
            NbTrial_voc(ii,cc) = NbTrial_voc(ii,cc) + length(FMAT.Trials{Icat_local(ss)});
        end
    end
    %retrieve values for background sounds
    BGInd = find(strcmp(MAT.VocTypeSel, 'BG'));
    TDTnames = MAT.TDTwav(BGInd);
    for bg=1:length(TDTnames)
        TDTname = TDTnames(bg);
        ITDTname = find(strcmp(FMAT.TDT_wavfiles, TDTname));
        NbTrial_voc(ii,NstimTypeCM) = NbTrial_voc(ii,NstimTypeCM) + length(FMAT.Trials_BG{ITDTname});
    end
end

%% Compile results into structures
sprintf('Compile results\n')
Selectivity.NVocPerCat = NVocPerCat(1:ii,:);
Selectivity.ProbaVocPerCat = ProbaVocPerCat(1:ii,:);
Selectivity.NVoc = NVoc(1:ii);
Selectivity.VocCatNames = VocCatNames(1:ii);
Selectivity.NbTrial_voc = NbTrial_voc(1:ii,:);
MeanSR.MeanSR=MeanSR;
MeanSR.SpikeRate_Cat = SpikeRate_Cat(1:ii,:);
MeanSR.SpikeRate_Voc = SpikeRate_Voc(1:ii);
MeanSR.TDT_names = TDT_names(1:ii);


%% Save data
save(fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict_PCCAbCh.mat'))
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict_PCCAbCh.mat'))



