
%load('/Users/elie/Documents/MATLAB/data/matfile/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance.mat', 'MI_perArea', 'SemanticIndex', 'R2A', 'GlobalR2A', 'Selectivity', 'Invariance', 'Conf', 'List_matfilepath', 'List_anat', 'ZONES', 'ZONES_List', 'SUBJ', 'Spike_shape', 'MeanSR', 'StimTypeCM','SUBJECTS')
%load('/Users/elie/Documents/MATLAB/data/matfile/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance2.mat', 'MI_perArea', 'SemanticIndex', 'R2A', 'GlobalR2A', 'Selectivity', 'Invariance', 'Conf', 'List_matfilepath', 'List_anat', 'ZONES', 'ZONES_List', 'SUBJ', 'Spike_shape', 'MeanSR', 'StimTypeCM','SUBJECTS');
%load('/Users/elie/Documents/MATLAB/data/matfile/SemanticAnalysis_ShuffPredict_PCCAbCh_Invarianceonly.mat', 'Invariance');
load('/Users/elie/Documents/CODE/data/matfile/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance3.mat', 'MI_perArea', 'SemanticIndex', 'R2A', 'GlobalR2A', 'Selectivity', 'Invariance', 'Conf', 'List_matfilepath', 'List_anat', 'ZONES', 'ZONES_List', 'SUBJ', 'Spike_shape', 'MeanSR', 'StimTypeCM','SUBJECTS');
%FE=load('/Users/elie/Documents/MATLAB/data/matfile/SemanticAnalysis_ShuffPredict_LNL.mat'); %NO USE ANYMORE: results of the analysis of Linearity that predict the expected value of MI_cat based on the linear model that predict values of MI_cat for different level of perturbation of the acoustic correlation of vocalizations within the confusion matrix
AS = load('/Users/elie/Documents/CODE/data/matfile/SemanticAnalysis_ShuffPredict_ASModels.mat'); %Results of the acoustic and semantic perturbations of confusion matrices
PCC_DFA=[52.32 90.96 42.12 41.5 41.04 67.64 82.92 68.97 58.64]./100;%according to email from frederic 03172014
NU = length(MI_perArea.MI_tot);
MU=find(Spike_shape==0);
SU=find(Spike_shape>0);


% First check that we have the same order of data
% if sum(strcmp(List_matfilepath, FE.List_matfilepath))==NU
%     fprintf('Tout est bon chez dupont\n')
% else
%     fprintf('Panique!!! Data are not in the same order!')
% end

if sum(strcmp(List_matfilepath, AS.List_matfilepath))==NU
    fprintf('Tout est bon chez dupont\n')
else
    fprintf('Panique!!! Data are not in the same order!')
end

%% Calculate the invariance index and find the max value of invariance for each unit
Ncat = length(StimTypeCM);
II = Invariance.PerCat ./Invariance.PerCatmax;
II2 = (Invariance.PerCat - Invariance.PerCatmin) ./(Invariance.PerCatmax-Invariance.PerCatmin);
II_max = max(II2(:,1:Ncat-1),[],2);
II_mean = mean(II2(:,1:Ncat-1),2);


%% get the coordinates of the cells and sites ready
SUBJ_Site = nan(NU,1);
LR = nan(NU,1);
DV = nan(NU,1);
RC = nan(NU,1);
LR_theo = LR;
DV_theo = DV;
RC_theo = RC;
LR_Site = nan(NU,1);
DV_Site = nan(NU,1);
RC_Site = nan(NU,1);
LR_theo_Site = LR;
DV_theo_Site = DV;
RC_theo_Site = RC;
Units_PerSite = cell(NU,1);
RC_theo_surf = [flip((1:8).*250) (1:8).*250 (1:8).*250 flip((1:8).*250)];% Rostro-caudal position of the 32 electrodes at the surface for all birds except LblBlu2028M (SUBJ=2)
RC_theo_surf2= [(1:8).*250 flip((1:8).*250) flip((1:8).*250) (1:8).*250];% Rostro-caudal position of the 32 electrodes at the surface for LblBlu2028M
LR_theo_surf = [repmat(-700,1,8) repmat(-1200,1,8) repmat(250,1,8) repmat(750,1,8)];%Medio-lateral position of the 32 electrodes at the surface for all birds except LblBlu2028M (SUBJ=2)
LR_theo_surf2 = [repmat(250,1,8) repmat(750,1,8) repmat(-700,1,8) repmat(-1200,1,8)];%Medio-lateral position of the 32 electrodes at the surface for LblBlu2028M
ii=0;
for uu=1:NU
    % Fisrt retrieve theoric values
    % what is the electrode number?
    i1=strfind(List_matfilepath{uu}, '_e');
    i2=strfind(List_matfilepath{uu}, '_s');
    i2=i2(1);
    E=str2num(List_matfilepath{uu}((i1+2):(i2-1)));
    % what is the deepness? find the position indicating deepeness in the
    % name of site
    Lpos=strfind(List_matfilepath{uu},'L');
    if length(Lpos)>1 %that's to take care of LblBlu2028 that also has an L in the name of bird
        Lpos = Lpos(end);
    end
    Rpos=strfind(List_matfilepath{uu},'R');
    if isempty(Rpos)
        Rpos = i1;% That's for WhiBlu5396M that does not have any R
    end
    endR=strfind(List_matfilepath{uu},'_e');
    endR=endR-1;
    %retrieve theoretical positions
    if SUBJ(uu,1)==2 %(LblBlu2028M)
        if LR_theo_surf2(E)<0 %left hemisphere
            DV_local=List_matfilepath{uu};
            DV_local=str2num(DV_local(Lpos+1:Rpos-1));
            DV_local2 = DV_local ./ cosd(15);
            LR_theo(uu) = LR_theo_surf2(E) + DV_local .* tand(15);
            RC_theo(uu) = RC_theo_surf2(E);
        else %right hemisphere
            DV_local=List_matfilepath{uu};
            DV_local=str2num(DV_local(Rpos+1:endR));
            DV_local2 = DV_local ./ cosd(17);
            LR_theo(uu) = LR_theo_surf2(E);
            RC_theo(uu) = RC_theo_surf2(E) + DV_local .* tand(17);
        end
    elseif SUBJ(uu,1)==4 %WhiBlu5396M values of motor were noted as entries for site names
        if LR_theo_surf(E)<0 %left hemisphere
            DV_local=List_matfilepath{uu};
            DV_local=str2num(DV_local(Lpos+1:Rpos-1));
            DV_local2 = DV_local ./(2.* cosd(15));
            LR_theo(uu) = LR_theo_surf(E) + DV_local .* tand(15)./2;
            RC_theo(uu) = RC_theo_surf(E);
        else %right hemisphere
            DV_local=List_matfilepath{uu};
            DV_local=str2num(DV_local(Rpos+1:endR));
            DV_local2 = DV_local ./ (2.*cosd(17));
            LR_theo(uu) = LR_theo_surf(E);
            RC_theo(uu) = RC_theo_surf(E) + DV_local .* tand(17)./2;
        end
    else
        if LR_theo_surf(E)<0 %left hemisphere
            DV_local=List_matfilepath{uu};
            DV_local=str2num(DV_local(Lpos+1:Rpos-1));
            DV_local2 = DV_local ./ cosd(15);
            LR_theo(uu) = LR_theo_surf(E) + DV_local .* tand(15);
            RC_theo(uu) = RC_theo_surf(E);
        else %right hemisphere
            DV_local=List_matfilepath{uu};
            DV_local=str2num(DV_local(Rpos+1:endR));
            DV_local2 = DV_local ./ cosd(17);
            LR_theo(uu) = LR_theo_surf(E);
            RC_theo(uu) = RC_theo_surf(E) + DV_local .* tand(17);
        end
    end
    if uu==1
        ii=ii+1;
        Units_PerSite{ii} = uu;
        DV_theo(uu)=-DV_local2;
        LR_theo_Site(ii) = LR_theo(uu);
        DV_theo_Site(ii) = DV_theo(uu);
        RC_theo_Site(ii) = RC_theo(uu);
        SUBJ_Site(ii)=SUBJ(uu);
    else
        if LR_theo(uu-1)==LR_theo(uu) && RC_theo(uu-1)==RC_theo(uu) %cells recorded at the excat same spot
            DV_theo(uu)=DV_theo(uu-1)-50;%place the sphere of that one a little lower so as to see all of them
            Units_PerSite{ii} = [Units_PerSite{ii} uu]; % Add that unit ID number to the list of cell for that electrode-site
        else %this is the position of a new site
            ii=ii+1;
            Units_PerSite{ii} = uu;
            DV_theo(uu)=-DV_local2;
            LR_theo_Site(ii) = LR_theo(uu);
            DV_theo_Site(ii) = DV_theo(uu);
            RC_theo_Site(ii) = RC_theo(uu);
            SUBJ_Site(ii)=SUBJ(uu);
        end
    end
    
    %only retrieve exact xyz positions of cells when we really know all of the
    %coordinates
    NAN=sum(isnan(cell2mat(List_anat(uu,3:4))));
    if NAN==0 && isempty(List_anat{uu,1})==0
        if SUBJ(uu,1)==3 %(GreBlu9508M) error in GreBlu9508M thickness of slices estimated to be 4x
            RC(uu)=4*1000*List_anat{uu,4}; % values are put in micrometers here
        else
            RC(uu)=1000*List_anat{uu,4}; % values are put in micrometers here
        end
        Lpos=cell2mat(strfind(List_anat{uu,1},'L'));
        Rpos=cell2mat(strfind(List_anat{uu,1},'R'));
        endR=cell2mat(strfind(List_anat{uu,1},'e'));
        endR=endR(2)-2;
        if strcmp(List_anat{uu,5}, 'L')
            LR(uu)=-1000*List_anat{uu,3};%if in left hemisphere put a negative sign on the medio-lateral value so left will be on the negative part of the 3D graph
            DV_local=cell2mat(List_anat{uu,1});
            DV_local=DV_local(Lpos+1:Rpos-1);
        else
            LR(uu)=1000*List_anat{uu,3};
            DV_local=cell2mat(List_anat{uu,1});
            DV_local=DV_local(Rpos+1:endR);
        end
        if uu==1
           DV(uu)=-str2num(DV_local);
           DV_Site(ii) = DV(uu);
           RC_Site(ii) = RC(uu);
           LR_Site(ii) = LR(uu);
        else
            if LR(uu-1)==LR(uu) && RC(uu-1)==RC(uu) %cells recorded at the excat same spot
                DV(uu)=DV(uu-1)-50;%place the sphere of that one a little lower so as to see all of them
            else %new electrode-site
                DV(uu)=-str2num(DV_local);
                DV_Site(ii) = DV(uu);
                RC_Site(ii) = RC(uu);
                LR_Site(ii) = LR(uu);
            end
        end
    end
end
DV_Site = DV_Site(1:ii);
RC_Site = RC_Site(1:ii);
LR_Site = LR_Site(1:ii);
LR_theo_Site = LR_theo_Site(1:ii);
DV_theo_Site = DV_theo_Site(1:ii);
RC_theo_Site = RC_theo_Site(1:ii);
Units_PerSite = Units_PerSite(1:ii);
SUBJ_Site = SUBJ_Site(1:ii);
NS = ii;%total number of electrode-sites

%% try to find a significant threshold of semantic info
sum(MI_perArea.AVMI_diag_uni_cat_RandP(:,1)>0.543)
sum(MI_perArea.AVMI_diag_uni_cat_RandP(:,2)>0.863) %for these matrices the response to silence where not shuffled
% if MI_cat>0.863 then the proba that MI_cat is randomly getting this value
% is below 0.01 according to the distribution of AVMI_diag_uni_cat_RandP BG intact So
% cells with values of MI_cat above 0.863 are semantic cells (another
% criteria)
SemCell = find(MI_perArea.MI_diag_uni_cat>0.863);
NonSemCell = setdiff(1:NU,SemCell);
%semantic cells based on p-value calculated per cell
SemCellPV = find(SemanticIndex.PvalueBGP<0.01);
NonSemCellPV = setdiff(1:NU,SemCellPV);
SemIndex = MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_RandP(:,2);
SemSU = intersect(SU, SemCellPV);
NonSemSU=intersect(SU,NonSemCellPV);
SemMU=intersect(MU,SemCellPV);

%% identifying semantic cells from non semantic ones
ID_Sem_NonSem = ones(NU);
ID_Sem_NonSem(SemCellPV) = 2;

%% significance of the semantic index
logpvRand = -log10(SemanticIndex.Pvalue);
logpvRandBG = -log10(SemanticIndex.PvalueBG);
logpvRandP = -log10(SemanticIndex.PvalueP);
logpvRandBGP = -log10(SemanticIndex.PvalueBGP);
DiffSem = MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_Rand(:,2);
DiffSemP = MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_RandP(:,2);

%% Explore linear and non-linear cells with acoustic/correlation disruption
% Semantic cells are
%AS.pAS_A<0.05
% Semantic and linear cells (gradient colored cells in figure12)
SemLN=find((AS.pAS_A<0.05).*(AS.pAS_S<0.05));
% Semantic and non-linear cells (Red cells in figure 12)
SemNLN=find((AS.pAS_A<0.05).*(AS.pAS_S>=0.05));%largest vector
% Non Semantic Cells
%AS.pAS_A>0.05
% Non Semantic acoustic cells (White cells in figure 12)
NonSemAc=find((AS.pAS_A>=0.05).*(AS.pAS_S<0.05));
% Non Semantic poorly acoustic cells (grey cells in figure 12)
NonSemNA=find((AS.pAS_A>=0.05).*(AS.pAS_S>=0.05));
% Create a vector that indicate the category of cell
LN_identity = zeros(NU,1);
LN_identity(SemLN)=1;
LN_identity(SemNLN)=2;
LN_identity(NonSemAc)=3;
LN_identity(NonSemNA)=4;

%% identifying $ categpry of cells (+S+A, +
ID_AS = ones(NU);%(1=+A+S units)
ID_AS(SemNLN) = 2;%(2=a+S units)
ID_AS(NonSemAc) = 3;%(3=+As units)
ID_AS(NonSemNA) = 4;%(4=as units)

%% Investigating the effect of anatomy
% Calculate the percentage of semantic units, the mean values of MI_tot and MI_diag_cat in each area
UZ = unique(ZONES(:,1));
NZ = length(UZ);
PercSSUCells_perarea = nan(NZ,3);
PercSSCells_perarea = nan(NZ,3);
TotCells_perarea = nan(2,NZ);
TotSSCells_perarea = nan(2,NZ);
MI_perzone = nan(NZ,2); % first column is mean, second is sem
MI_HSCells_perzone=MI_perzone;
MI_NonHSCells_perzone=MI_perzone;
MI_SSCells_perzone=MI_perzone;
MI_NonSSCells_perzone=MI_perzone;
SemanticIndex_perzone = MI_perzone;
SemanticIndex_HSCells_perzone=MI_perzone;
SemanticIndex_NonHSCells_perzone=MI_perzone;
SemanticIndex_SSCells_perzone=MI_perzone;
SemanticIndex_NonSSCells_perzone=MI_perzone;
SemIndex_perzone = MI_perzone;
SemIndex_HSCells_perzone=MI_perzone;
SemIndex_NonHSCells_perzone=MI_perzone;
SemIndex_SSCells_perzone=MI_perzone;
SemIndex_NonSSCells_perzone=MI_perzone;
UnitsAllZonesSS=cell(1,NZ);
UnitsAllZonesSSU=cell(1,NZ);


for zz = 1:NZ
    TotCells_perarea(1,zz) = sum(ZONES(:,1) == UZ(zz));
    TotCells_perarea(2,zz) = sum(ZONES(SU,1) == UZ(zz));
    Unitinzone = find(ZONES(:,1) == UZ(zz));
    UnitinzoneSSU = intersect(Unitinzone, SemSU);
    UnitinzoneNonSSU = intersect(Unitinzone, NonSemSU);
    UnitinzoneSS = intersect(Unitinzone, SemCellPV);
    TotSSCells_perarea(1,zz) = length(UnitinzoneSS);
    TotSSCells_perarea(2,zz) = length(UnitinzoneSSU);
    UnitsAllZonesSS{zz} = UnitinzoneSS;
    UnitsAllZonesSSU{zz} = UnitinzoneSSU;
    UnitinzoneNonSS = intersect(Unitinzone, NonSemCellPV);
    PercSSUCells_perarea(zz,1) = length(UnitinzoneSSU)/TotCells_perarea(2,zz);
    PercSSCells_perarea(zz,1) = length(UnitinzoneSS)/TotCells_perarea(1,zz);
    PercSSUCells_perarea(zz,2) = length(intersect(UnitinzoneSSU,SemLN))/length(UnitinzoneSSU);
    PercSSCells_perarea(zz,2) = length(intersect(UnitinzoneSS,SemLN))/length(UnitinzoneSS);
    PercSSUCells_perarea(zz,3) = length(intersect(UnitinzoneSSU,SemNLN))/length(UnitinzoneSSU);
    PercSSCells_perarea(zz,3) = length(intersect(UnitinzoneSS,SemNLN))/length(UnitinzoneSS);
    MI_perzone(zz,1) = mean(MI_perArea.MI_tot(Unitinzone));
    MI_perzone(zz,2) = std(MI_perArea.MI_tot(Unitinzone))/sqrt(length(Unitinzone));
    MI_HSCells_perzone(zz,1)=mean(MI_perArea.MI_tot(UnitinzoneSSU));
    MI_SSCells_perzone(zz,1)=mean(MI_perArea.MI_tot(UnitinzoneSS));
    MI_HSCells_perzone(zz,2)=std(MI_perArea.MI_tot(UnitinzoneSSU))/sqrt(length(UnitinzoneSSU));
    MI_SSCells_perzone(zz,2)=std(MI_perArea.MI_tot(UnitinzoneSS))/sqrt(length(UnitinzoneSS));
    MI_NonHSCells_perzone(zz,1)=mean(MI_perArea.MI_tot(UnitinzoneNonSSU));
    MI_NonHSCells_perzone(zz,2)=std(MI_perArea.MI_tot(UnitinzoneNonSSU))/sqrt(length(UnitinzoneNonSSU));
    MI_NonSSCells_perzone(zz,1)=mean(MI_perArea.MI_tot(UnitinzoneNonSS));
    MI_NonSSCells_perzone(zz,2)=std(MI_perArea.MI_tot(UnitinzoneNonSS))/sqrt(length(UnitinzoneNonSS));
    SemanticIndex_perzone(zz,1) = mean(SemanticIndex.Observed(Unitinzone));
    SemanticIndex_perzone(zz,2) = std(SemanticIndex.Observed(Unitinzone))/sqrt(length(Unitinzone));
    SemanticIndex_HSCells_perzone(zz,1)=mean(SemanticIndex.Observed(UnitinzoneSSU));
    SemanticIndex_HSCells_perzone(zz,2)=std(SemanticIndex.Observed(UnitinzoneSSU))/sqrt(length(UnitinzoneSSU));
    SemanticIndex_NonHSCells_perzone(zz,1)=mean(SemanticIndex.Observed(UnitinzoneNonSSU));
    SemanticIndex_NonHSCells_perzone(zz,2)=std(SemanticIndex.Observed(UnitinzoneNonSSU))/sqrt(length(UnitinzoneNonSSU));
    SemanticIndex_SSCells_perzone(zz,1)=mean(SemanticIndex.Observed(UnitinzoneSS));
    SemanticIndex_SSCells_perzone(zz,2)=std(SemanticIndex.Observed(UnitinzoneSS))/sqrt(length(UnitinzoneSS));
    SemanticIndex_NonSSCells_perzone(zz,1)=mean(SemanticIndex.Observed(UnitinzoneNonSS));
    SemanticIndex_NonSSCells_perzone(zz,2)=std(SemanticIndex.Observed(UnitinzoneNonSS))/sqrt(length(UnitinzoneNonSS));
    SemIndex_perzone(zz,1) = mean(SemIndex(Unitinzone));
    SemIndex_perzone(zz,2) = std(SemIndex(Unitinzone))/sqrt(length(Unitinzone));
    SemIndex_HSCells_perzone(zz,1)=mean(SemIndex(UnitinzoneSSU));
    SemIndex_HSCells_perzone(zz,2)=std(SemIndex(UnitinzoneSSU))/sqrt(length(UnitinzoneSSU));
    SemIndex_NonHSCells_perzone(zz,1)=mean(SemIndex(UnitinzoneNonSSU));
    SemIndex_NonHSCells_perzone(zz,2)=std(SemIndex(UnitinzoneNonSSU))/sqrt(length(UnitinzoneNonSSU));
    SemIndex_SSCells_perzone(zz,1)=mean(SemIndex(UnitinzoneSS));
    SemIndex_SSCells_perzone(zz,2)=std(SemIndex(UnitinzoneSS))/sqrt(length(UnitinzoneSS));
    SemIndex_NonSSCells_perzone(zz,1)=mean(SemIndex(UnitinzoneNonSS));
    SemIndex_NonSSCells_perzone(zz,2)=std(SemIndex(UnitinzoneNonSS))/sqrt(length(UnitinzoneNonSS));
end

%units index classified as "L1, L2, L3, CMM, CML, NCM, unknown"
Units7ZonesSS = cell(2,7);
Units7ZonesSS{1,1}=UnitsAllZonesSS{4};
Units7ZonesSS{1,2}=[UnitsAllZonesSS{5};UnitsAllZonesSS{6}];
Units7ZonesSS{1,3}=UnitsAllZonesSS{7};
Units7ZonesSS{1,4}=UnitsAllZonesSS{2};
Units7ZonesSS{1,5}=UnitsAllZonesSS{3};
Units7ZonesSS{1,6}=UnitsAllZonesSS{9};
Units7ZonesSS{1,7}=[UnitsAllZonesSS{1};UnitsAllZonesSS{8};UnitsAllZonesSS{10}];
Units7ZonesSS{2,1}=UnitsAllZonesSSU{4};
Units7ZonesSS{2,2}=[UnitsAllZonesSSU{5};UnitsAllZonesSSU{6}];
Units7ZonesSS{2,3}=UnitsAllZonesSSU{7};
Units7ZonesSS{2,4}=UnitsAllZonesSSU{2};
Units7ZonesSS{2,5}=UnitsAllZonesSSU{3};
Units7ZonesSS{2,6}=UnitsAllZonesSSU{9};
Units7ZonesSS{2,7}=[UnitsAllZonesSSU{1};UnitsAllZonesSSU{8};UnitsAllZonesSSU{10}];

Units2ZonesSS = cell(2,2);
Units2ZonesSS{1,1} = [UnitsAllZonesSS{4};UnitsAllZonesSS{5};UnitsAllZonesSS{6};UnitsAllZonesSS{7}];
Units2ZonesSS{1,2} = [UnitsAllZonesSS{2};UnitsAllZonesSS{3};UnitsAllZonesSS{9}];
Units2ZonesSS{2,1} = [UnitsAllZonesSSU{4};UnitsAllZonesSSU{5};UnitsAllZonesSSU{6};UnitsAllZonesSSU{7}];
Units2ZonesSS{2,2} = [UnitsAllZonesSSU{2};UnitsAllZonesSSU{3};UnitsAllZonesSSU{9}];

TotCells_per7area = zeros(2,7);
TotCells_per7area(:,1)=TotCells_perarea(:,4);
TotCells_per7area(:,2)=TotCells_perarea(:,5) + TotCells_perarea(:,6);
TotCells_per7area(:,3)=TotCells_perarea(:,7);
TotCells_per7area(:,4)=TotCells_perarea(:,2);
TotCells_per7area(:,5)=TotCells_perarea(:,3);
TotCells_per7area(:,6)=TotCells_perarea(:,9);
TotCells_per7area(:,7)=TotCells_perarea(:,1) + TotCells_perarea(:,8) + TotCells_perarea(:,10);

TotSSCells_per7area = zeros(2,7);
TotSSCells_per7area(:,1)=TotSSCells_perarea(:,4);
TotSSCells_per7area(:,2)=TotSSCells_perarea(:,5) + TotSSCells_perarea(:,6);
TotSSCells_per7area(:,3)=TotSSCells_perarea(:,7);
TotSSCells_per7area(:,4)=TotSSCells_perarea(:,2);
TotSSCells_per7area(:,5)=TotSSCells_perarea(:,3);
TotSSCells_per7area(:,6)=TotSSCells_perarea(:,9);
TotSSCells_per7area(:,7)=TotSSCells_perarea(:,1) + TotSSCells_perarea(:,8) + TotSSCells_perarea(:,10);

ZONES7=zeros(NU,1);
ZONES7(find(ZONES==3))=1;
ZONES7([find(ZONES==4);find(ZONES==5)])=2;
ZONES7(find(ZONES==6))=3;
ZONES7(find(ZONES==1))=4;
ZONES7(find(ZONES==2))=5;
ZONES7(find(ZONES==8))=6;
ZONES7([find(ZONES==9);find(ZONES==7);find(ZONES==0)])=7;


ZONES2=zeros(NU,1);
ZONES2([find(ZONES==3);find(ZONES==4);find(ZONES==5);find(ZONES==6);find(ZONES==7)])=1;
ZONES2([find(ZONES==1);find(ZONES==2);find(ZONES==8)])=2;


PercSSCells_per7area = zeros(7,5);
for zz=1:length(Units7ZonesSS)
    PercSSCells_per7area(zz,1) = size(Units7ZonesSS{1,zz},1)./TotCells_per7area(1,zz);
    PercSSCells_per7area(zz,2) = length(intersect(Units7ZonesSS{1,zz},SemLN))./TotCells_per7area(1,zz);
    PercSSCells_per7area(zz,3) = length(intersect(Units7ZonesSS{1,zz},SemNLN))./TotCells_per7area(1,zz);
    PercSSCells_per7area(zz,4) = length(intersect(Units7ZonesSS{1,zz},NonSemAc))./TotCells_per7area(1,zz);
    PercSSCells_per7area(zz,5) = length(intersect(Units7ZonesSS{1,zz},NonSemNA))./TotCells_per7area(1,zz);
end

PercSSUCells_per7area = zeros(7,5);
for zz=1:length(Units7ZonesSS)
    PercSSUCells_per7area(zz,1) = size(Units7ZonesSS{2,zz},1)./TotCells_per7area(2,zz);
    PercSSUCells_per7area(zz,2) = length(intersect(Units7ZonesSS{2,zz},SemLN))./TotCells_per7area(2,zz);
    PercSSUCells_per7area(zz,3) = length(intersect(Units7ZonesSS{2,zz},SemNLN))./TotCells_per7area(2,zz);
    PercSSUCells_per7area(zz,4) = length(intersect(Units7ZonesSS{2,zz},NonSemAc))./TotCells_per7area(2,zz);
    PercSSUCells_per7area(zz,5) = length(intersect(Units7ZonesSS{2,zz},NonSemNA))./TotCells_per7area(2,zz);
end

List_zones={'Unknown',ZONES_List{1:end}};
List_zones7={'L1','L2','L3','CMM','CML','NCM','Unknown'};

% find the cells for which we know the position
Loc = cell(3,1);
Loc{1} = find(ZONES(:,1)>0);%all cells
Loc{2} = intersect(Loc{1},SemCell);% HS cells
Loc{3} = intersect(Loc{1},SemCellPV);% SS cells
CTs = {'HS cells', 'SS cells'};

%% Find the max selectivity value for each unit
LRI = cell2mat(Selectivity.DiagLRI')';
[LRI_max, LRI_max_CT] = max(LRI');
InfInd = find(LRI_max==Inf);
maxLRI_max = max(LRI_max([1:(InfInd-1),(InfInd+1):end]));
LRI_max(InfInd)=maxLRI_max+1;

%% Calculate the proportion of semantic units that are selective (LRI>1.6) in each zone for each category
LRI = cell2mat(Selectivity.DiagLRI')';
PercSelSSCells_per7area = zeros(7,Ncat);
PercSelSSUCells_per7area = zeros(7,Ncat);
for zz=1:length(Units7ZonesSS)
    PercSelSSCells_per7area(zz,1) = sum(sum(LRI(Units7ZonesSS{1,zz},:)>1.6,2)>0)./length(Units7ZonesSS{1,zz});%Propotion of selective cells among the semantic cells in that area
    PercSelSSUCells_per7area(zz,1) = sum(sum(LRI(Units7ZonesSS{2,zz},:)>1.6,2)>0)./length(Units7ZonesSS{2,zz});%Propotion of selective cells among the semantic single units in that area
    for ct=1:(Ncat-1)
        PercSelSSCells_per7area(zz,ct+1) = sum(LRI(Units7ZonesSS{1,zz},ct)>1.6)./length(Units7ZonesSS{1,zz});%Propotion of cell selective for calltype ct among the semantic cells in that area
        PercSelSSUCells_per7area(zz,ct+1) = sum(LRI(Units7ZonesSS{2,zz},ct)>1.6)./length(Units7ZonesSS{2,zz});%Propotion of cell selective for calltype ct among the semantic single units in that area
    end
end

%% Indentify the Invariance value of the category that is best discriminated by each unit
II_PCCmax=zeros(NU,1);
for ii=1:NU
    II_PCCmax(ii) = II2(ii,LRI_max_CT(ii));
end

%% Calculate the proportion of semantic and selective (LRI>1.6) units that are more invariant (>median) in each zone for each category
PercInvSelSSCells_per7area = zeros(7,Ncat);
PercInvSelSSUCells_per7area = zeros(7,Ncat);
Inv_thresh = nanmedian(II_PCCmax(SemCellPV));
Inv_threshSU = nanmedian(II_PCCmax(SemSU));
for zz=1:length(Units7ZonesSS)
    %indices of semantic and selective units in that zz area
    I_SemSel = intersect(Units7ZonesSS{1,zz},find(LRI_max>1.6));
    I_SemSelSU = intersect(Units7ZonesSS{2,zz},find(LRI_max>1.6));
    %indices of semantic and selective and more invariant units in that
    %area
    I_SemSelInv = intersect(I_SemSel,find(II_PCCmax>Inv_thresh));
    I_SemSelInvSU = intersect(I_SemSelSU,find(II_PCCmax>Inv_threshSU));
    PercInvSelSSCells_per7area(zz,1) = length(I_SemSelInv)./length(I_SemSel);
    PercInvSelSSUCells_per7area(zz,1) = length(I_SemSelInvSU)./length(I_SemSelSU);
    for ct=1:(Ncat-1)
        PercInvSelSSCells_per7area(zz,ct+1) = sum(LRI_max_CT(I_SemSelInv)==ct)./length(I_SemSel);%Propotion of semantic and selective cell invariant for calltype ct among the semantic and selective cells in that area
        PercInvSelSSUCells_per7area(zz,ct+1) = sum(LRI_max_CT(I_SemSelInvSU)==ct)./length(I_SemSelSU);%Propotion of semantic and selective single units invariant for calltype ct among the semantic and selective single units in that area
    end
end
%% Discrimination performance compare to DFA: investigating PCC 
% Calculate mean PCC for each cell over call categories
PCC = cell2mat(Selectivity.PCC_cat');
PCC_inv = PCC';
PCC_mean = mean(PCC_inv(:,1:9),2);
PCC_se = std(PCC(1:9,:))./sqrt(size(PCC_inv,2));
PCC_se = PCC_se';

%For highly semantic cells HS
PCC_HScell = cell2mat(Selectivity.PCC_cat(SemCell)');
PCC_HScell_inv = PCC_HScell';
PCC_HScell_mean = mean(PCC_HScell_inv,1);
PCC_HScell_mean_DFA = PCC_HScell_mean(1:(end-1))./PCC_DFA;
PCC_HScell_se = std(PCC_HScell_inv)./sqrt(size(PCC_HScell_inv,1));
PCC_HScell_inv_DFA = zeros(length(SemCell), length(PCC_DFA));
for ii=1:size(PCC_HScell_inv,1)
    PCC_HScell_inv_DFA(ii,:)= PCC_HScell_inv(ii,1:(end-1)) ./ PCC_DFA;
end
PCC_HScell_se_DFA = std(PCC_HScell_inv_DFA)./sqrt(size(PCC_HScell_inv_DFA,1));

%For significantly semantic cells
PCC_SScell = cell2mat(Selectivity.PCC_cat(SemCellPV)');
PCC_SSUcell = cell2mat(Selectivity.PCC_cat(SemSU)');
PCC_SScell_inv = PCC_SScell(1:9,:)';
PCC_SSUcell_inv = PCC_SSUcell(1:9,:)';
PCC_SScell_mean = mean(PCC_SScell_inv,1);
PCC_SSUcell_mean = mean(PCC_SSUcell_inv,1);
[PCC_SScell_max, PCC_SScell_max_ind] = max(PCC_SScell_inv);
[PCC_SSUcell_max, PCC_SSUcell_max_ind] = max(PCC_SSUcell_inv);
%PCC_SScell_mean_DFA = PCC_SScell_mean(1:(end-1))./PCC_DFA;%old single value of DFA
PCC_SScell_se = std(PCC_SScell_inv)./sqrt(size(PCC_SScell_inv,1));
PCC_SSUcell_se = std(PCC_SSUcell_inv)./sqrt(size(PCC_SSUcell_inv,1));
% PCC_SScell_inv_DFA = zeros(length(SemCellPV), length(PCC_DFA));
% for ii=1:size(PCC_SScell_inv,1)
%     PCC_SScell_inv_DFA(ii,:)= PCC_SScell_inv(ii,1:(end-1)) ./ PCC_DFA;
% end
% PCC_SScell_se_DFA = std(PCC_SScell_inv_DFA)./sqrt(size(PCC_SScell_inv_DFA,1));

DFA=load('/Users/elie/Documents/CODE/data/matfile/DFAAcoustic.mat');
PCC_DFA_mean = mean(DFA.PCC_cat);
PCC_DFA_se = std(DFA.PCC_cat)./sqrt(size(DFA.PCC_cat,1));

%% Identify cells that have significant PCC compare to chance level for each call category
NbTrial_voc = Selectivity.NbTrial_voc; 
p=Selectivity.ProbaVocPerCat(:,:);
PCC_effectif = PCC_inv(:,:).*NbTrial_voc;
y = 1-binocdf(PCC_effectif, NbTrial_voc, p);
Signif_PCC = cell(Ncat-1,1);
PCC_signifSScell_mean=nan(2,Ncat-1);
PCC_signifSScell_se = PCC_signifSScell_mean;

for cc=1:Ncat-1
    Signif_PCC{cc} = find(y(:,cc)<0.01);
    PCC_signifSScell_mean(1,cc) = mean(PCC_inv(intersect(SemCellPV, Signif_PCC{cc}),cc));
    PCC_signifSScell_se(1,cc) = std(PCC_inv(intersect(SemCellPV, Signif_PCC{cc}),cc))./sqrt(length(intersect(SemCellPV, Signif_PCC{cc})));
    PCC_signifSScell_mean(2,cc) = mean(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc));
    PCC_signifSScell_se(2,cc) = std(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc))./sqrt(length(intersect(SemSU, Signif_PCC{cc})));
end
%% Identify cells that have significant PCC compare to chance level for at least one call category
Signif_PCC_Cell = find(sum(y(:,(1:(Ncat-1)))<0.01,2)~=0);
%% Calculate the number of category significantly discriminated by each cell
Signif_PCC_Cellnb = sum(y(:,(1:(Ncat-1)))<0.01,2);
mean(Signif_PCC_Cellnb(SemSU))
std(Signif_PCC_Cellnb(SemSU))
min(Signif_PCC_Cellnb(SemSU))
max(Signif_PCC_Cellnb(SemSU))


%% Calculate average invariance of each unit for the categories for which it is significantly discriminant
II_Mean_disccat = nan(NU,1);
for ii=1:NU
    Signif_cat = find(y(ii,1:9)<0.01);
    II_Mean_disccat(ii)=nanmean(II2(ii,Signif_cat));
end

%% calculate mean average PCC and selectivity over categories per zone
NZ7=size(Units7ZonesSS(2,:),2);
PCC_AV = mean(PCC(1:(end-1),:),1);
[PCC_max, PCC_max_CT] = max(PCC(1:(end-1),:));
PCC_AVstd = std(PCC(1:(end-1),:));
Mean_PCC_AV_perzone = nan(NZ7,1);
Mean_PCC_AV_perzoneSU = nan(NZ7,1);
Mean_GlobalSel_AV_perzoneSU = nan(NZ7,1);
SE_GlobalSel_AV_perzoneSU = nan(NZ7,1);
Mean_PCC_max_SU_perzone = nan(NZ7,1);
Mean_PCC_max_SS_perzone = nan(NZ7,1);
SE_PCC_AV_perzone = nan(NZ7,1);
SE_PCC_AV_perzoneSU = nan(NZ7,1);
SE_PCC_max_HS_perzone = nan(NZ7,1);
SE_PCC_max_SS_perzone = nan(NZ7,1);
Mean_PCC_HS_AV_perzone = nan(NZ7,1);
SE_PCC_HS_AV_perzone = nan(NZ7,1);
Mean_PCC_SS_AV_perzone = nan(NZ7,1);
SE_PCC_SS_AV_perzone = nan(NZ7,1);
Mean_II_Meandiscat = nan(NZ7,1);
SE_II_Meandiscat = nan(NZ7,1);
for zz = 1:NZ7
    Unitinzone = find(ZONES(:,1) == UZ(zz));
    UnitinzoneSSU = Units7ZonesSS{2,zz};
    UnitinzoneSS = intersect(Unitinzone, SemCellPV);
    Mean_PCC_AV_perzone(zz)=mean(PCC_AV(Unitinzone));
    Mean_PCC_AV_perzoneSU(zz)=mean(PCC_AV(UnitinzoneSSU));
    Mean_GlobalSel_AV_perzoneSU(zz)=mean(Selectivity.DiagEntropyNorm(UnitinzoneSSU));
    SE_GlobalSel_AV_perzoneSU(zz)=std(Selectivity.DiagEntropyNorm(UnitinzoneSSU))./sqrt(length(UnitinzoneSSU));
    Mean_II_Meandiscat(zz) = nanmean(II_Mean_disccat(UnitinzoneSSU));
    SE_II_Meandiscat(zz) = nanstd(II_Mean_disccat(UnitinzoneSSU))./sqrt(length(II_Mean_disccat(UnitinzoneSSU)));
    Mean_PCC_max_SU_perzone(zz)=mean(PCC_max(UnitinzoneSSU));
    Mean_PCC_max_SS_perzone(zz)=mean(PCC_max(UnitinzoneSS));
    SE_PCC_AV_perzone(zz) = std(PCC_AV(Unitinzone))./sqrt(length(Unitinzone));
    SE_PCC_AV_perzoneSU(zz) = std(PCC_AV(UnitinzoneSSU))./sqrt(length(UnitinzoneSSU));
    SE_PCC_max_HS_perzone(zz)=std(PCC_max(UnitinzoneSSU))./sqrt(length(UnitinzoneSSU));
    SE_PCC_max_SS_perzone(zz)=std(PCC_max(UnitinzoneSS))./sqrt(length(UnitinzoneSS));
    Mean_PCC_HS_AV_perzone(zz)=mean(PCC_AV(UnitinzoneSSU));
    SE_PCC_HS_AV_perzone(zz) = std(PCC_AV(UnitinzoneSS))./sqrt(length(UnitinzoneSS));
     Mean_PCC_SS_AV_perzone(zz)=mean(PCC_AV(UnitinzoneSS));
    SE_PCC_SS_AV_perzone(zz) = std(PCC_AV(UnitinzoneSS))./sqrt(length(UnitinzoneSS));
end


%% Calculate PCC significance compare to chance level
% identify the values of PCC that are above chance for each call category
% using y = binocdf(x,N,p,'upper') with p=Selectivity.ProbaVocPerCat and
% x=PCC values and N = nb voc * nb trials=Selectivity.NbTrial_voc.
% Signif_PCC contains the indices of units that show significant PCC for
% each category
VerifCallTypeNU=0;
for ii=1:NU
    VocCatName=Selectivity.VocCatNames{ii};
    VerifCallType=0;
    for jj=1:length(VocCatName)
        VerifCallType = strcmp(VocCatName{jj},StimTypeCM(jj)) + VerifCallType;
    end
    if VerifCallType==length(StimTypeCM) && length(VocCatName)==length(StimTypeCM)
        VerifCallTypeNU = VerifCallTypeNU + 1;
    end
end
if VerifCallTypeNU==NU
    fprintf('Tout est bon chez Dupont\n');
else
    fprintf('ALERT:Problem with call type correspondence\n');
end
 


%% Calculate the proportion of semantic units that are signifiacntly discriminant in each zone for each category
PercDiscSSCells_per7area = zeros(6,Ncat-1);
PercDiscSSUCells_per7area = zeros(6,Ncat-1);
for zz=1:length(Units7ZonesSS)
    for ct=1:(Ncat-1)
        PercDiscSSCells_per7area(zz,ct) = length(intersect(Signif_PCC{ct},Units7ZonesSS{1,zz}))./length(Units7ZonesSS{1,zz});%Propotion of cell discriminant for calltype ct among the semantic cells in that area
        PercDiscSSUCells_per7area(zz,ct) = length(intersect(Signif_PCC{ct},Units7ZonesSS{2,zz}))./length(Units7ZonesSS{2,zz});%Propotion of cell selective for calltype ct among the semantic single units in that area
    end
end

%% Calculate the proportion of semantic units that are selective (LRI>1.6) in each zone for each category
LRI = cell2mat(Selectivity.DiagLRI')';
PercSelDiscSSCells_per7area = zeros(7,Ncat-1);
PercSelDiscSSUCells_per7area = zeros(7,Ncat-1);
for zz=1:length(Units7ZonesSS)
    for ct=1:(Ncat-1)
        PercSelDiscSSCells_per7area(zz,ct) = sum(LRI(intersect(Signif_PCC{ct},Units7ZonesSS{1,zz}),ct)>1.6)./length(intersect(Signif_PCC{ct},Units7ZonesSS{1,zz}));%Propotion of cell selective for calltype ct among the semantic discriminative cells in that area
        PercSelDiscSSUCells_per7area(zz,ct) = sum(LRI(intersect(Signif_PCC{ct},Units7ZonesSS{2,zz}),ct)>1.6)./length(intersect(Signif_PCC{ct},Units7ZonesSS{2,zz}));%Propotion of cell selective for calltype ct among the semantic discriminative single units in that area
    end
end

%% Calculate the proportion of semantic units that are selective (LRI>1.6) in each zone for each category
PercInvDiscSSCells_per7area = zeros(7,Ncat-1);
PercInvDiscSSUCells_per7area = zeros(7,Ncat-1);
II_SSU_Signif_PCC=[];
II_SS_Signif_PCC=[];
for ct=1:(Ncat-1)
    TC=intersect(Signif_PCC{ct},SemCellPV);
    TCU=intersect(Signif_PCC{ct},SemSU);
    II_SSU_Signif_PCC=[II_SSU_Signif_PCC;II(TCU,ct)];
    II_SS_Signif_PCC=[II_SS_Signif_PCC;II(TC,ct)];
end
Inv_thresh=nanmedian(II_SS_Signif_PCC);
Inv_threshSU=nanmedian(II_SSU_Signif_PCC);

for zz=1:length(Units7ZonesSS)
    for ct=1:(Ncat-1)
        TC=intersect(Signif_PCC{ct},Units7ZonesSS{1,zz});
        TCU=intersect(Signif_PCC{ct},Units7ZonesSS{2,zz});
        PercInvDiscSSCells_per7area(zz,ct) = length(intersect(TC,find(II(:,ct)>Inv_thresh)))./length(TC);%Propotion of cell selective for calltype ct among the semantic discriminative cells in that area
        PercInvDiscSSUCells_per7area(zz,ct) = length(intersect(TCU,find(II(:,ct)>Inv_threshSU)))./length(TCU);%Propotion of cell selective for calltype ct among the semantic discriminative single units in that area
    end
end

%% Stuffs about the temporal vs rate code: size of confusion matrix windows
hist(Conf.optWin_CT,600)
mean(Conf.optWin_CT)
std(Conf.optWin_CT)
MI_CT_gain=zeros(NU,1);
MI_gain=zeros(NU,1);
for uu=1:NU
    if max(Conf.MI_confBCT(uu,:))==Conf.MIBCT_OptCT(uu)
        MI_CT_gain(uu) = max(Conf.MI_confBCT(uu,:))-min(Conf.MI_confBCT(uu,:));
        MI_gain(uu) = max(Conf.MI_confB(uu,:))-min(Conf.MI_confB(uu,:));
    else
        frpintf(1,'Pb of MI choice for %d unit', uu)
    end
end
mean(MI_gain)
std(MI_gain)
mean(MI_CT_gain)
std(MI_CT_gain)
max(max(Conf.MI_confB,[],2))
max(Conf.MIBCT_OptCT)
max(Conf.MIB_OptCT)
quantile(Conf.optWin_CT, [0.25 0.5 0.75])

%% Calculate the proportion of stims that were cut or kept intact for each units
Perc_FullStim = NB.Nb_FullStim./(NB.Nb_CutStim+NB.Nb_FullStim);
mean(Perc_FullStim)
min(Perc_FullStim)



%% PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT

%% Histogram of semantic information for observed and random matrices for all cells

%%%%%%%%%%FIG3b SEMANTICITY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,X]=hist([MI_perArea.MI_diag_uni_cat MI_perArea.AVMI_diag_uni_cat_RandP], 40);
figure(2)
ss=subplot(1,1,1)
hist([MI_perArea.MI_diag_uni_cat(SemCellPV) [MI_perArea.MI_diag_uni_cat(NonSemCellPV); nan(length(SemCellPV)-length(NonSemCellPV),1)]],X)
h = findobj(gca,'Type','patch');
set(h(2),'FaceColor',[0.5 0 0.5],'EdgeColor',[0.5 0 0.5])
set(h(1),'FaceColor',[1 0.59 1],'EdgeColor',[1 0.59 1])
axis([0 3.7 0 300])
legend('Significantly Semantic Cells (SS cells)', 'Non-SS cells')
legend BOXOFF
%line([0.863 0.863], [0 300], 'Color', [0 .3 .6])
%text(0.88, 250, sprintf('Semantic information\nthreshold for highly\nsemantic cells'), 'Color', [0 .3 .6])
set(get(ss, 'YLabel'), 'String', '# Units');
set(get(ss, 'XLabel'), 'String', 'Semantic Information for all cells');
text(1.5, 150, sprintf('%d out of %d units are significantly semantic (SS cells)', length(SemCellPV),NU));
text(1.5, 130, sprintf('%d out of %d units are highly semantic (HS cells)', sum(MI_perArea.MI_diag_uni_cat>=0.863),NU));

figure(2)
ss=subplot(1,1,1)
hist([[MI_perArea.MI_diag_uni_cat(SemSU);nan(length(NonSemSU)-length(SemSU),1)] MI_perArea.MI_diag_uni_cat(NonSemSU)],X)
h = findobj(gca,'Type','patch');
set(h(2),'FaceColor',[0.5 0 0.5],'EdgeColor',[0.5 0 0.5])
set(h(1),'FaceColor',[1 0.59 1],'EdgeColor',[1 0.59 1])
axis([0 3.7 0 300])
legend('Significantly Semantic Cells (SS single units)', 'Non-SS cells')
legend BOXOFF
%line([0.863 0.863], [0 300], 'Color', [0 .3 .6])
%text(0.88, 250, sprintf('Semantic information\nthreshold for highly\nsemantic cells'), 'Color', [0 .3 .6])
set(get(ss, 'YLabel'), 'String', '# Units');
set(get(ss, 'XLabel'), 'String', 'Semantic Information for all cells');
text(1.5, 150, sprintf('%d out of %d single units are significantly semantic (SS cells)', length(SemSU),length(SU)));
text(1.5, 130, sprintf('%d out of %d single units are highly semantic (HS cells)', sum(MI_perArea.MI_diag_uni_cat>=0.863),length(SU)));

% h = findobj(gca,'Type','patch');
% set(h(1),'FaceColor',[0 .3 .6],'EdgeColor',[0 .3 .6])
% set(h(2),'FaceColor',[0 .8 .5],'EdgeColor',[0 .8 .5])
% line([0.543 0.543], [0 600], 'Color', [0 .8 .5])
% text(0.55, 500, 'Acoustic information threshold', 'Color', [0 .8 .5])
% line([0.863 0.863], [0 600], 'Color', [0 .3 .6])
% text(0.87, 500, sprintf('Semantic\ninformation\nthreshold'), 'Color', [0 .3 .6])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%FIG3a SEMANTICITY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find an example SS cell to show the distribution of random values of
% MI_cat
ExSSCell = find(MI_perArea.MI_diag_uni_cat>2);
ExSSCell_SU = intersect(SemSU,ExSSCell);
figure()
hist(MI_perArea.MI_uni_diag_cat_RandBGP(ExSSCell(1),:), X)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 0.5 0.5],'EdgeColor',[0 0.5 0.5])
axis([0 3.7 0 600])
ylabel('Number of shuffled matrices');
xlabel('Semantic Information');
line([MI_perArea.MI_diag_uni_cat(ExSSCell(1)) MI_perArea.MI_diag_uni_cat(ExSSCell(1))], [0 600], 'Color', [0.5 0 .5])
text(2.3, 300, sprintf('Semantic information\nvalue observed\nfor that semantic cell'), 'Color', [0.5 0 0.5])
figure()
hist(MI_perArea.MI_uni_diag_cat_RandBGP(ExSSCell_SU(1),:), X)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 0.5 0.5],'EdgeColor',[0 0.5 0.5])
axis([0 3.7 0 600])
ylabel('Number of shuffled matrices');
xlabel('Semantic Information');
line([MI_perArea.MI_diag_uni_cat(ExSSCell_SU(1)) MI_perArea.MI_diag_uni_cat(ExSSCell_SU(1))], [0 600], 'Color', [0.5 0 .5])
text(2.3, 300, sprintf('Semantic information\nvalue observed\nfor that semantic cell'), 'Color', [0.5 0 0.5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate correlation value between PCC and Semantic Information
[Rho, Pval]=corr(MI_perArea.MI_diag_uni_cat(SemSU), PCC_mean(SemSU),'Spearman');
%% Explore non-Linearity of the cells with the correlation between values...
...of MI_Cat and % of disruption of the semantic categories of the confusion...
    ...matrices or % acoustic correlation within disrupted categories

%%%%%%%FIGURE 3 SEMANTICITY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% last histogram of different cell categories for SemCellPV only
ss=subplot(3,2,5)
Vec1 = [MI_perArea.MI_diag_uni_cat(intersect(SemLN,SemCellPV));nan(length(intersect(SemCellPV,SemNLN))-length(intersect(SemCellPV,SemLN)),1)];
Vec2 = MI_perArea.MI_diag_uni_cat(intersect(SemCellPV,SemNLN));
Vec3 = [MI_perArea.MI_diag_uni_cat(intersect(SemCellPV,NonSemAc));nan(length(intersect(SemCellPV,SemNLN))-length(intersect(SemCellPV,NonSemAc)),1)];
Vec4 = [MI_perArea.MI_diag_uni_cat(intersect(SemCellPV,NonSemNA));nan(length(intersect(SemCellPV,SemNLN))-length(intersect(SemCellPV,NonSemNA)),1)];
hist([Vec1 Vec2 Vec3 Vec4],X)
%h = findobj(gca,'Type','patch');
%set(h(1),'FaceColor',[0.5 0 0.5],'EdgeColor',[0.5 0 0.5])
%set(h(2),'FaceColor',[1 0 1],'EdgeColor',[1 0 1])
axis([0 3.7 0 110])
legend('+A+S','A+S','+AS','AS')
legend BOXOFF
%line([0.863 0.863], [0 300], 'Color', [0 .3 .6])
%text(0.88, 80, sprintf('Semantic\ninformation\nthreshold\nfor highly\nsemantic cells'), 'Color', [0 .3 .6])
set(get(ss, 'YLabel'), 'String', '# Units');
set(get(ss, 'XLabel'), 'String', 'Semantic Information of significantly semantic cells (SS cells)');
text(1.5, 50, sprintf('%d out of %d units are +A+S', length(intersect(SemCellPV,SemLN)),length(SemCellPV)));
text(1.5, 40, sprintf('%d out of %d units are A+S', length(intersect(SemCellPV,SemNLN)),length(SemCellPV)));
text(1.5, 30, sprintf('%d out of %d units are +AS', length(intersect(SemCellPV,NonSemAc)),length(SemCellPV)));
text(1.5, 20, sprintf('%d out of %d units are AS', length(intersect(SemCellPV,NonSemNA)),length(SemCellPV)));

figure()
ss=subplot(1,1,1);
Vec1 = [MI_perArea.MI_diag_uni_cat(intersect(SemLN,SemSU));nan(length(intersect(SemSU,SemNLN))-length(intersect(SemSU,SemLN)),1)];
Vec2 = MI_perArea.MI_diag_uni_cat(intersect(SemSU,SemNLN));
Vec3 = [MI_perArea.MI_diag_uni_cat(intersect(SemSU,NonSemAc));nan(length(intersect(SemSU,SemNLN))-length(intersect(SemSU,NonSemAc)),1)];
Vec4 = [MI_perArea.MI_diag_uni_cat(intersect(SemSU,NonSemNA));nan(length(intersect(SemSU,SemNLN))-length(intersect(SemSU,NonSemNA)),1)];
hist([Vec1 Vec2 Vec3 Vec4],X)
%h = findobj(gca,'Type','patch');
%set(h(1),'FaceColor',[0.5 0 0.5],'EdgeColor',[0.5 0 0.5])
%set(h(2),'FaceColor',[1 0 1],'EdgeColor',[1 0 1])
axis([0 3.7 0 80])
legend('+A+S','A+S','+AS','AS')
legend BOXOFF
%line([0.863 0.863], [0 300], 'Color', [0 .3 .6])
%text(0.88, 80, sprintf('Semantic\ninformation\nthreshold\nfor highly\nsemantic cells'), 'Color', [0 .3 .6])
set(get(ss, 'YLabel'), 'String', '# Single Units');
set(get(ss, 'XLabel'), 'String', 'Semantic Information of significantly semantic cells (SS cells)');
text(1.5, 50, sprintf('%d out of %d units are +A+S', length(intersect(SemSU,SemLN)),length(SemSU)));
text(1.5, 40, sprintf('%d out of %d units are A+S', length(intersect(SemSU,SemNLN)),length(SemSU)));
text(1.5, 30, sprintf('%d out of %d units are +AS', length(intersect(SemSU,NonSemAc)),length(SemSU)));
text(1.5, 20, sprintf('%d out of %d units are AS', length(intersect(SemSU,NonSemNA)),length(SemSU)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%FIG3 SEMANTICITY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the fitting curves for the perturbation experiment REMEMBER you have
% to change the path for the correct correlation matrix of sound in
% testRandCorr depending on the unit choosen to work with
testRandCorr(List_matfilepath{ExSSCell(1)}); %a BlaBro09xxF site2 unit
testRandCorr(List_matfilepath{ExSSCell_SU(1)});%a GreBlu9508M Site2 Unit
AS.pAS_S(ExSSCell_SU(1))
AS.FvalAS_S(ExSSCell_SU(1))
AS.r2adjA(ExSSCell_SU(1))
AS.r2adjS(ExSSCell_SU(1))
AS.r2adjAS(ExSSCell_SU(1))

% Get the bar graph of the number of cells per category (possible insert of 3D) 
figure()
bar([length(intersect(SemCellPV,SemLN)) length(intersect(SemCellPV,SemNLN)) length(intersect(SemCellPV,NonSemAc)) length(intersect(SemCellPV,NonSemNA))], 'k');
%ylim([0 1])
set(gca, 'XTickLabel',{'+A+S' 'A+S' '+AS' 'AS'});
ylabel('Number of significantly semantic units');
[P,Chi2stat,Df]=chi2test([length(intersect(SemCellPV,SemLN)) length(intersect(SemCellPV,SemNLN)) length(intersect(SemCellPV,NonSemAc)) length(intersect(SemCellPV,NonSemNA))])%significant
figure()
bar([length(intersect(SemSU,SemLN)) length(intersect(SemSU,SemNLN)) length(intersect(SemSU,NonSemAc)) length(intersect(SemSU,NonSemNA))], 'k');
%ylim([0 1])
set(gca, 'XTickLabel',{'+A+S' 'A+S' '+AS' 'AS'});
ylabel('Number of significantly semantic single units');
[P,Chi2stat,Df]=chi2test([length(intersect(SemSU,SemLN)) length(intersect(SemSU,SemNLN)) length(intersect(SemSU,NonSemAc)) length(intersect(SemSU,NonSemNA))])%significant

% Version 2 of the previous figure with percentage instead of actual
% numbers
figure()
bar([length(intersect(SemCellPV,SemLN))./length(SemCellPV) length(intersect(SemCellPV,SemNLN))./length(SemCellPV) length(intersect(SemCellPV,NonSemAc))./length(SemCellPV) length(intersect(SemCellPV,NonSemNA))./length(SemCellPV)], 'k');
set(gca, 'XTickLabel',{'+A+S' 'A+S' '+AS' 'AS'});
ylabel('Proportion of significantly semantic units');
figure()
bar([length(intersect(SemSU,SemLN))./length(SemSU) length(intersect(SemSU,SemNLN))./length(SemSU) length(intersect(SemSU,NonSemAc))./length(SemSU) length(intersect(SemSU,NonSemNA))./length(SemSU)], 'k');
set(gca, 'XTickLabel',{'+A+S' 'A+S' '+AS' 'AS'});
ylabel('Proportion of significantly semantic single units');

% average R2adj values for +S cells
mean(AS.r2adjAS([intersect(SemSU,SemLN) ;intersect(SemSU,SemNLN)]))
std(AS.r2adjAS([intersect(SemSU,SemLN) ;intersect(SemSU,SemNLN)]))

% Get the bar graph of the mean and sd values of R2adj AS per cell category
% (insert for 3D)
MEAN_R2adj = [mean(AS.r2adjAS(intersect(SemCellPV,SemLN))) mean(AS.r2adjAS(intersect(SemCellPV,SemNLN))) mean(AS.r2adjAS(intersect(SemCellPV,NonSemAc))) mean(AS.r2adjAS(intersect(SemCellPV,NonSemNA)))];
SD_R2adj = [std(AS.r2adjAS(intersect(SemCellPV,SemLN))) std(AS.r2adjAS(intersect(SemCellPV,SemNLN))) std(AS.r2adjAS(intersect(SemCellPV,NonSemAc))) std(AS.r2adjAS(intersect(SemCellPV,NonSemNA)))];
figure()
ss=subplot(1,1,1);
bar(MEAN_R2adj, 'k')
hold on
errorb(MEAN_R2adj,SD_R2adj);
set(gca, 'XTickLabel',{'AS' 'aS' 'As' 'as'});
ylabel('Adjusted R2')
title(ss,'all semantic units')
hold off
MEAN_R2adj = [mean(AS.r2adjAS(intersect(SemSU,SemLN))) mean(AS.r2adjAS(intersect(SemSU,SemNLN))) mean(AS.r2adjAS(intersect(SemSU,NonSemAc))) mean(AS.r2adjAS(intersect(SemSU,NonSemNA)))];
SD_R2adj = [std(AS.r2adjAS(intersect(SemSU,SemLN))) std(AS.r2adjAS(intersect(SemSU,SemNLN))) std(AS.r2adjAS(intersect(SemSU,NonSemAc))) std(AS.r2adjAS(intersect(SemSU,NonSemNA)))];
figure()
ss=subplot(1,1,1);
bar(MEAN_R2adj, 'k')
hold on
errorb(MEAN_R2adj,SD_R2adj);
set(gca, 'XTickLabel',{'AS' 'aS' 'As' 'as'});
ylabel('Adjusted R2')
title(ss,'Semantic single units')
hold off
[p,table,stats]=anova1(AS.r2adjAS(SemCellPV), LN_identity(SemCellPV))
c=multcompare(stats,'ctype','scheffe', 'alpha',0.01)
[p,table,stats]=anova1(AS.r2adjAS(SemSU), LN_identity(SemSU))
c=multcompare(stats,'ctype','scheffe','alpha',0.01)

% Get the bar graph of the number of semantic and non semantic cells
% (insert for 3B)
figure()
bar([length(NonSemCellPV) length(SemCellPV)], 'k');
%ylim([0 1])
set(gca, 'XTickLabel',{'Non Significantly Semantic' 'Significantly Semantic'});
ylabel('Number of units')
figure()
bar([length(NonSemSU) length(SemSU)], 'k');
%ylim([0 1])
set(gca, 'XTickLabel',{'Non Significantly Semantic' 'Significantly Semantic'});
ylabel('Number of single units')
% Get the bar graph of the proportion of semantic and non semantic cells
% (insert for 3B)
figure()
bar([length(NonSemCellPV)./NU length(SemCellPV)./NU], 'k');
%ylim([0 1])
set(gca, 'XTickLabel',{'Non Significantly Semantic' 'Significantly Semantic'});
ylabel('Proportion of units')
figure()
bar([length(NonSemSU)./length(SU) length(SemSU)./length(SU)], 'k');
%ylim([0 1])
set(gca, 'XTickLabel',{'Non Significantly Semantic' 'Significantly Semantic'});
ylabel('Proportion of single units')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%% Where are the semantic cells? Investigate the anatomy of MI confusion values

%%%%%%%%%%%%%%%%%%%%%%%%%%%Fig4d Anatomy1
%%%%%%%%%%%%%%%%%%%%%%%%%%%paper%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,2,2)
figure(18)
NZ7 = size(PercSSCells_per7area,1);
[BBaxes, BBbar, BBline]=plotyy(1:NZ7, PercSSCells_per7area(:,2:end), 1:NZ7, TotCells_per7area(1,:), 'bar', 'plot');
for bb=1:4
    set(BBbar(bb), 'BarLayout', 'stacked')
end
set(BBline,'LineWidth',2,'Color',[1,0,0],'Marker','o');
%set(BBbar, 'FaceColor', 'k')
LL=line([0,8], [length(SemCellPV)/NU length(SemCellPV)/NU]);
set(LL, 'Color', [0 0 0], 'LineStyle','--');
xlim([0 8])
ylabel(BBaxes(1),'Proportion of Semantic units');
ylabel(BBaxes(2), 'Number of units recorded');
xlabel('anatomical zones');
set(BBaxes(1),'XTickLabel',List_zones7);
set(BBaxes(2),'XTickLabel',{});
%T=text(11,0.53, sprintf('Total percentage\nof SS units\n'))
T=text(4.5,0.55, sprintf('Total percentage of Semantic units\n'))
set(T, 'Color', [0 0 0])
% MPSC=mean(PercSSCells_per7area(:,1));
% MPSCse = std(PercSSCells_per7area(:,1))/sqrt(length(PercSSCells_per7area));
% LL=line([0,8], [MPSC MPSC]);
% set(LL, 'Color', [1 0.5 0]);
% LL=line([0,8], [MPSC-MPSCse MPSC-MPSCse]);
% set(LL, 'Color', [1 0.5 0], 'LineStyle', '--');
% LL=line([0,8], [MPSC+MPSCse MPSC+MPSCse]);
% set(LL, 'Color', [1 0.5 0], 'LineStyle', '--')
% T=text(6,MPSC-0.02, sprintf('Average percentage\nof Semantic units\n'))
% set(T, 'Color', [1 0.5 0])
figure(18)
NZ7 = size(PercSSUCells_per7area,1);
[BBaxes, BBbar, BBline]=plotyy(1:NZ7, PercSSUCells_per7area(:,2:end), 1:NZ7, TotCells_per7area(2,:), 'bar', 'plot');
for bb=1:4
    set(BBbar(bb), 'BarLayout', 'stacked')
end
set(BBline,'LineWidth',2,'Color',[1,0,0],'Marker','o');
%set(BBbar, 'FaceColor', 'k')
LL=line([0,8], [length(SemSU)/length(SU) length(SemSU)/length(SU)]);
set(LL, 'Color', [0 0 0], 'LineStyle','--');
xlim([0 8])
ylabel(BBaxes(1),'Proportion of Semantic single units');
ylabel(BBaxes(2), 'Number of single units recorded');
xlabel('anatomical zones');
set(BBaxes(1),'XTickLabel',List_zones7);
set(BBaxes(2),'XTickLabel',{});
%T=text(11,0.53, sprintf('Total percentage\nof SS units\n'))
T=text(4.5,0.55, sprintf('Total percentage of Semantic single units\n'))
set(T, 'Color', [0 0 0])
ylim(BBaxes(1), [0 1])
ylim(BBaxes(2), [0 500])

% Final figure with error bars (confidence interval instead of nb of units
% with red line)
[phat,pci] = binofit(PercSSUCells_per7area(:,1)'.*TotCells_per7area(2,:),TotCells_per7area(2,:),0.05);
figure(18)
NZ7 = size(PercSSUCells_per7area,1);
H=bar(1:NZ7, PercSSUCells_per7area(:,2:end), 'stacked');
LL=line([0,8], [length(SemSU)/length(SU) length(SemSU)/length(SU)]);
set(LL, 'Color', [0 0 0], 'LineStyle','--');
ylabel('Proportion of Semantic single units');
xlabel('anatomical zones');
set(gca,'XTickLabel',List_zones7);
%T=text(11,0.53, sprintf('Total percentage\nof SS units\n'))
T=text(4.5,0.55, sprintf('Total percentage of Semantic single units\n'))
set(T, 'Color', [0 0 0])
ylim([0 1])
%plot the confidence intervals
hold on
for ee=1:size(pci,1)
    plot([ee ee],[pci(ee,1) pci(ee,2)],'color','k');
    hold on
    plot([ee-0.1 ee+0.1],[pci(ee,1) pci(ee,1)],'color','k');
    hold on
    plot([ee-0.1 ee+0.1],[pci(ee,2) pci(ee,2)],'color','k');
end
hold off

% Test that the number of semantic units are different between regions
[P,Chi2stat,Df]=chi2test([PercSSCells_per7area(1:6,1).*TotCells_per7area(1,1:6)' TotCells_per7area(1,1:6)'-PercSSCells_per7area(1:6,1).*TotCells_per7area(1,1:6)'])%significant
[P,Chi2stat,Df]=chi2test([PercSSUCells_per7area(1:6,1).*TotCells_per7area(2,1:6)' TotCells_per7area(2,1:6)'-PercSSUCells_per7area(1:6,1).*TotCells_per7area(2,1:6)'])%significant

% Test that the proportion of CellType is different between regions
[P,Chi2stat,Df]=chi2test(PercSSCells_per7area(1:6,2:end).*repmat(TotCells_per7area(1,1:6)',1,4))%NOT significant
[P,Chi2stat,Df]=chi2test(PercSSUCells_per7area(1:6,2:end).*repmat(TotCells_per7area(2,1:6)',1,4))%NOT significant
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Selectivity and comparison with Acoustic DFA Investigate the PCC of semantic cells
%%%%%%%%Figure4 Selectivity Comparison with DFA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bargraph of the mean PCC values of significantly discriminant semantic
% cells
figure()
errorb([PCC_DFA_mean; 100.*PCC_signifSScell_mean(1,:)]', [PCC_DFA_se; 100.*PCC_signifSScell_se(1,:)]')
set(gca, 'XTickLabel', StimTypeCM);
ylabel('PCC')
legend('DFA', 'Semantic Units')

figure()
errorb([PCC_DFA_mean; 100.*PCC_signifSScell_mean(2,:)]', [PCC_DFA_se; 100.*PCC_signifSScell_se(2,:)]')
set(gca, 'XTickLabel', StimTypeCM);
ylabel('PCC')
legend('DFA', 'Semantic Single Units', 'Location','NorthWest')

% Compute an anova2 on the PCC (effect of DFA vs Units and Call Type)
ANOVA_PCC_SS = nan(size(DFA.PCC_cat,1).*size(DFA.PCC_cat,2)+size(PCC_SScell_inv,1).*size(PCC_SScell_inv,2),1);
ANOVA_group_SS1 = cell(size(ANOVA_PCC_SS));
ANOVA_group_SS2 = ANOVA_group_SS1;
aa=0;
for vv=1:size(DFA.PCC_cat,1)
    DFAi=DFA.PCC_cat(vv,:);
    VocType=DFA.VocTypes{vv};
    for cc=1:length(DFAi)
        aa=aa+1;
        ANOVA_PCC_SS(aa)=DFAi(cc);
        ANOVA_group_SS1{aa}=VocType{cc};
        ANOVA_group_SS2{aa}='DFA';
    end
end
ANOVA_PCC_SSU = ANOVA_PCC_SS;
ANOVA_group_SSU1 = ANOVA_group_SS1;
ANOVA_group_SSU2 = ANOVA_group_SS2;
bb=aa;

for cc=1:(length(StimTypeCM)-1)
    PCC_cc=PCC_inv(intersect(SemCellPV, Signif_PCC{cc}),cc);
    for uu=1:length(PCC_cc)
        aa=aa+1;
        ANOVA_PCC_SS(aa)=100.*PCC_cc(uu);
        ANOVA_group_SS1{aa} = StimTypeCM{cc};
        ANOVA_group_SS2{aa}='Unit';
    end
end
for cc=1:(length(StimTypeCM)-1)
    PCC_cc=PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc);
    for uu=1:length(PCC_cc)
        bb=bb+1;
        ANOVA_PCC_SSU(bb)=100.*PCC_cc(uu);
        ANOVA_group_SSU1{bb} = StimTypeCM{cc};
        ANOVA_group_SSU2{bb}='Unit';
    end
end
 ANOVA_PCC_SSU = ANOVA_PCC_SSU(1:bb);
 ANOVA_group_SSU1 = ANOVA_group_SSU1(1:bb);
 ANOVA_group_SSU2 = ANOVA_group_SSU2(1:bb);
 ANOVA_PCC_SS = ANOVA_PCC_SS(1:aa);
 ANOVA_group_SS1 = ANOVA_group_SS1(1:aa);
 ANOVA_group_SS2 = ANOVA_group_SS2(1:aa);
    
[p,table,stats]=anovan(ANOVA_PCC_SS, {ANOVA_group_SS1 ANOVA_group_SS2}, 'varnames',{'CallType';'DFAvsUNITS'},'model','interaction')    
multcompare(stats, 'dimension',[1 2],'alpha',0.01)
[p,table,stats]=anovan(ANOVA_PCC_SSU, {ANOVA_group_SSU1 ANOVA_group_SSU2}, 'varnames',{'CallType';'DFAvsUNITS'},'model','interaction')    
multcompare(stats, 'dimension',[1 2],'alpha',0.01)
%any correlation between DFA and unit results??
corrcoef(PCC_DFA_mean,100.*PCC_signifSScell_mean(2,:))
[rho,pval]=corr(PCC_DFA_mean',100.*PCC_signifSScell_mean(2,:)','type','Spearman')

% Analyze Cell PCC with a Kruskall Wallis%%% Effect of call type
KW_dataSU=nan(length(SemSU).*9,1);
KW_groupSU=cell(size(KW_dataSU));
KW_dataDFA=nan(size(DFA.PCC_cat,1).*9,1);
KW_groupDFA=cell(size(KW_dataDFA));
bb=0;
aa=0;
for cc=1:9
    PCC_cc_SU=PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc);
    
    for uu=1:length(intersect(SemSU, Signif_PCC{cc}))
        bb=bb+1;
        KW_dataSU(bb)=PCC_cc_SU(uu);
        KW_groupSU{bb}=StimTypeCM{cc};
    end
    for dd=1:size(DFA.PCC_cat,1)
        aa=aa+1;
        KW_dataDFA(aa)=DFA.PCC_cat(dd,cc);
        KW_groupDFA{aa}=StimTypeCM{cc};
    end
    % Compare units and DFA performance with a Mann whitney test for each
    % category
    %[p,h,stats]=ranksum(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc), DFA.PCC_cat(:,cc))
end
KW_dataSU=KW_dataSU(1:bb);
KW_groupSU=KW_groupSU(1:bb);
KW_dataDFA=KW_dataDFA(1:aa);
KW_groupDFA=KW_groupDFA(1:aa);
[p,tbl,stats]=kruskalwallis(KW_dataSU,KW_groupSU)%significant
c=multcompare(stats, 'alpha',0.01)
[p,tbl,stats]=kruskalwallis(KW_dataDFA,KW_groupDFA)%significant
c=multcompare(stats, 'alpha',0.01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Look at values of PCC above chance per category per region




% represent in selectivity and discrimination space only the cells that are
% significantly discriminating + insert of the number of significantly
% discriminating vs not discriminating cells
% Plot for each type call the PCC values per category vs the LRI for that category for HS cells with regions color coded


%%%%%%%FIGURE4 Remember to answer 0 to cubehelix and just manually change
%%%%%%%the spike rate multiplier by dividing by 10^3 (color scale)
PCC_DFA_mean_div=PCC_DFA_mean./100;
PCC_DFA_se_div=PCC_DFA_se./100;
GRAD1=cubehelix(ceil(max(max(MeanSR.SpikeRate_Cat(SemCellPV,:)))), 0.5, -1.1, 1.5, 0.5, [1,0]);
figure(56)
for cc=1:(Ncat-1)
    subplot(3,3,cc)
    PCC_x=PCC_inv(intersect(SemCellPV, Signif_PCC{cc}),cc);
    LRI_y = LRI(intersect(SemCellPV, Signif_PCC{cc}),cc);
    MeanSR_z = MeanSR.SpikeRate_Cat(intersect(SemCellPV, Signif_PCC{cc}),cc);
    NU_local = length(PCC_x);
    for jj=1:NU_local
        plot(PCC_x(jj), LRI_y(jj),'ko', 'MarkerFaceColor',GRAD1(ceil(MeanSR_z(jj)),:));
        hold on
    end
    xlabel('PCC')
    ylabel('LSI')
    title(sprintf('%s',StimTypeCM{cc}));
    %colorbar('YTickLabel', [10^(-20) 10^-40 10^-60 10^-80 10^-100 10^-120 10^-140 10^-160]);
    colormap(GRAD1)
    axis([0 1 -3 4])
    LL2=line([0 1], [1 1]);
    set(LL2, 'Color', 'r','LineWidth', 1,'LineStyle','--')
    LL2=line([0 1], [-1 -1]);
    set(LL2, 'Color', 'r', 'LineWidth', 1,'LineStyle','--')
    LL2=line([0 1], [0 0]);
    set(LL2, 'Color', 'r','LineWidth', 1)
    LL1=line([PCC_DFA_mean_div(cc) PCC_DFA_mean_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1)
    LL1=line([PCC_DFA_mean_div(cc)+PCC_DFA_se_div(cc) PCC_DFA_mean_div(cc)+PCC_DFA_se_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1,'LineStyle','--')
    LL1=line([PCC_DFA_mean_div(cc)-PCC_DFA_se_div(cc) PCC_DFA_mean_div(cc)-PCC_DFA_se_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1,'LineStyle','--')
    hold off
end
colorbar()

GRAD1=cubehelix(ceil(max(max(II(SemCellPV,:).*1000))), 0.5, -1.6, 1.5, 1.4, [1,0]);
figure(56)%bis version with Invariance index as z axis
for cc=1:(Ncat-1)
    subplot(3,3,cc)
    PCC_x=PCC_inv(intersect(SemCellPV, Signif_PCC{cc}),cc);
    LRI_y = LRI(intersect(SemCellPV, Signif_PCC{cc}),cc);
    II_z = II(intersect(SemCellPV, Signif_PCC{cc}),cc).*1000;
    ZZ=find(II_z==0);
    II_z(ZZ)=1;%make sure that we don't have zero indices and apply the lowest value
    NU_local = length(PCC_x);
    for jj=1:NU_local
        if ~isnan(II_z(jj))
            plot(PCC_x(jj), LRI_y(jj),'ko', 'MarkerFaceColor',GRAD1(ceil(II_z(jj)),:));
            hold on
        end
    end
    xlabel('PCC')
    ylabel('LSI')
    title(sprintf('%s',StimTypeCM{cc}));
    %colorbar('YTickLabel', [10^(-20) 10^-40 10^-60 10^-80 10^-100 10^-120 10^-140 10^-160]);
    colormap(GRAD1)
    axis([0 1 -3 4])
    LL2=line([0 1], [1 1]);
    set(LL2, 'Color', 'r','LineWidth', 1,'LineStyle','--')
    LL2=line([0 1], [-1 -1]);
    set(LL2, 'Color', 'r', 'LineWidth', 1,'LineStyle','--')
    LL2=line([0 1], [0 0]);
    set(LL2, 'Color', 'r','LineWidth', 1)
    LL1=line([PCC_DFA_mean_div(cc) PCC_DFA_mean_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1)
    LL1=line([PCC_DFA_mean_div(cc)+PCC_DFA_se_div(cc) PCC_DFA_mean_div(cc)+PCC_DFA_se_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1,'LineStyle','--')
    LL1=line([PCC_DFA_mean_div(cc)-PCC_DFA_se_div(cc) PCC_DFA_mean_div(cc)-PCC_DFA_se_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1,'LineStyle','--')
    hold off
end
colorbar()

% GRAD1=cubehelix(ceil(max(max(II2(SemCellPV,:).*1000))), 0.5, -1.6, 1.5, 1.4, [1,0]);
% GRAD1=cubehelix(ceil(max(max(II2(SemCellPV,:).*1000))), 1.8, -0.4, 2.5, 1.4, [1,0]);
GRAD1=cubehelix(ceil(max(max(II2(SemSU,:).*1000))), 1.7, -0.4, 2.5, 1, [1,0]);
II2_SemSUDisc=[];
PCC_SemSUDisc=[];
LRI_SemSUDisc=[];
figure(57)%bis version with Invariance index as z axis
for cc=1:(Ncat-1)
    subplot(3,3,cc)
    PCC_x=PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc);
    LRI_y = LRI(intersect(SemSU, Signif_PCC{cc}),cc);
    II_z = II2(intersect(SemSU, Signif_PCC{cc}),cc);
    II2_SemSUDisc = [II2_SemSUDisc ;II_z];
    PCC_SemSUDisc=[PCC_SemSUDisc; PCC_x];
    LRI_SemSUDisc=[LRI_SemSUDisc; LRI_y];
    max(II_z)
    II_bzeros = find(II_z<0);% Minmum value happen to be -1e-16 for some categories
    II_z(II_bzeros)=0;
    II_z = II_z.*1000;
    ZZ=find(II_z==0);
    II_z(ZZ)=1;%make sure that we don't have zero indices and apply the lowest value
    NU_local = length(PCC_x);
    for jj=1:NU_local
        if ~isnan(II_z(jj))
            plot(PCC_x(jj), LRI_y(jj),'ko', 'MarkerFaceColor',GRAD1(ceil(II_z(jj)),:));
            hold on
        end
    end
    xlabel('PCC')
    ylabel('LSI')
    title(sprintf('%s',StimTypeCM{cc}));
    %colorbar('YTickLabel', [10^(-20) 10^-40 10^-60 10^-80 10^-100 10^-120 10^-140 10^-160]);
    %colormap(GRAD1)
    axis([0 1 -3 4])
    LL2=line([0 1], [1.75 1.75]);
    set(LL2, 'Color', 'r','LineWidth', 1,'LineStyle','--')
    LL2=line([0 1], [-1.75 -1.75]);
    set(LL2, 'Color', 'r', 'LineWidth', 1,'LineStyle','--')
    LL2=line([0 1], [0 0]);
    set(LL2, 'Color', 'r','LineWidth', 1)
    LL1=line([PCC_DFA_mean_div(cc) PCC_DFA_mean_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1)
    LL1=line([PCC_DFA_mean_div(cc)+PCC_DFA_se_div(cc) PCC_DFA_mean_div(cc)+PCC_DFA_se_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1,'LineStyle','--')
    LL1=line([PCC_DFA_mean_div(cc)-PCC_DFA_se_div(cc) PCC_DFA_mean_div(cc)-PCC_DFA_se_div(cc)], [-3 4]);
    set(LL1,'Color', 'b', 'LineWidth', 1,'LineStyle','--')
    hold on
end
colorbar()
colormap(GRAD1)

quantile(PCC_SemSUDisc,0.95)
quantile(LRI_SemSUDisc,0.95)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%Fig Selectivity,Discrimination Invariance%%%%%%%%%%%%%
% Bar graph du nombre de units avec un PCC significatif pour chaque
% categorie de cri
figure(57)
for cc=1:(Ncat-1)
    ss=subplot(3,3,cc);
    NbSS_PCCsignif = length(intersect(SemCellPV, Signif_PCC{cc}));
    NbSS_PCCnonsignif = length(SemCellPV) - NbSS_PCCsignif;
    bar([NbSS_PCCnonsignif NbSS_PCCsignif], 'k')
    title(ss,sprintf('%s',StimTypeCM{cc}))
    set(ss, 'YLim', [0 600]);
end

NbSS_PCCsignif=nan(Ncat-1,1);
NbSS_PCCnonsignif=nan(Ncat-1,1);
figure(57)
for cc=1:(Ncat-1)
    ss=subplot(3,3,cc);
    NbSS_PCCsignif(cc) = length(intersect(SemSU, Signif_PCC{cc}));
    NbSS_PCCnonsignif(cc) = length(SemSU) - NbSS_PCCsignif(cc);
    bar([NbSS_PCCnonsignif(cc) NbSS_PCCsignif(cc)], 'k')
    title(ss,sprintf('%s',StimTypeCM{cc}))
    set(ss, 'YLim', [0 300]);
end
% Test that the number of discriminative units are different between call
% types
[P,Chi2stat,Df]=chi2test([NbSS_PCCsignif;NbSS_PCCnonsignif])
% Test that the proportion of discriminative units is correlated with the
% performance of these units
corrcoef(NbSS_PCCsignif./length(SemSU),100.*PCC_signifSScell_mean(2,:))
[rho,pval]=corr(NbSS_PCCsignif./length(SemSU),100.*PCC_signifSScell_mean(2,:)', 'type', 'Kendall')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find a threshold to consider units as selective. A selective units is a unit that would have no more than one category above the selectivity threshold
figure()
for ii=0:0.1:3
    plot(ii,sum(sum(LRI(SemCellPV,:)>ii,2)>1),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','b')
    hold on
    plot(ii,sum(sum(LRI(SemCellPV,:)>ii,2)>0),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','r')
    hold on
end
%According to this figure, LRI=1.6 or 2 might be better threshold to consider
%cells as selective compare to LRI=1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%selectivity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
for ii=0:0.05:3
    plot(ii,sum(sum(LRI(SemSU,:)>ii,2)>1)/length(SemSU),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','b')
    hold on
    plot(ii,sum(sum(LRI(SemSU,:)>ii,2)>0)/length(SemSU),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','r')
    hold on
end
legend('selective for more than 1 category','selective for at least 1 category')
xlabel('Sel')
ylabel('Proportion of Semantic Units')
sum(sum(LRI(SemSU,:)>1.75,2)>1)
length(SemSU)*0.05
%According to this figure, LRI=1.75 (5%) or 2.1(1%) might be better threshold to consider
%single units as selective compare to LRI=1.
Sel_Thresh=1.75;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Figure4 selectivity%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find for each call category the number of semantic units that are selective to
% that category (higher than threshold value) among the significantly
% discriminating units
% Construct the vector of selectivity
SelCat_aboveThresh=zeros(NU, Ncat-1);
for uu=1:NU
    SelCat_aboveThresh(uu,:)=(LRI(uu,:)>Sel_Thresh);
end
Units_SelThresh_perCatSU=zeros(2,Ncat-1);
Nb_SignifSUnits=zeros(Ncat-1,1);
for cc=1:Ncat-1
    Units_SelThresh_perCatSU(1,cc)= sum((y(SemSU,cc)<0.01).*SelCat_aboveThresh(SemSU,cc));
    Units_SelThresh_perCatSU(2,cc)=Units_SelThresh_perCatSU(1,cc)./sum((y(SemSU,cc)<0.01));
    Nb_SignifSUnits(cc)=sum((y(SemSU,cc)<0.01));
end

figure()
bar(Units_SelThresh_perCatSU(2,:),'k')
set(gca,'XTickLabel',StimTypeCM)
ylabel('Proportion of selective single units among discriminant ones');

%plot the confidence intervals
[phat,pci] = binofit(Units_SelThresh_perCatSU(1,:),Nb_SignifSUnits',0.05);
hold on
for ee=1:size(pci,1)
    plot([ee ee],[pci(ee,1) pci(ee,2)],'color','k');
    hold on
    plot([ee-0.1 ee+0.1],[pci(ee,1) pci(ee,1)],'color','k');
    hold on
    plot([ee-0.1 ee+0.1],[pci(ee,2) pci(ee,2)],'color','k');
end
hold off
ylim([0 0.4])


% Test if the number of units selective for the different call type is
% different
[p,Chi2stat,Df]=chi2test([Units_SelThresh_perCatSU(1,:);Nb_SignifSUnits'-Units_SelThresh_perCatSU(1,:)])%highly significant
% Any correlation between the performance of DFA and the proportion of
% significantly discriminant units?
corrcoef(PCC_DFA_mean,Nb_SignifSUnits./length(SemSU))
[rho,pval]=corr(PCC_DFA_mean',Nb_SignifSUnits./length(SemSU)')
corrcoef(PCC_DFA_mean,Units_SelThresh_perCatSU(2,:))
[rho,pval]=corr(PCC_DFA_mean',Units_SelThresh_perCatSU(2,:)')

% Find for each call category the number of semantic units that are selective to
% that category (highest value of PCC similar to highest value of LRI)
% Construct the vector of selectivity
SelCat_PCCmax=zeros(NU, Ncat-1);
for uu=1:NU
    SelCat_PCCmax(uu,:)=(PCC_inv(uu,1:Ncat-1)==max(PCC_inv(uu,1:Ncat-1)));
end
Units_Sel_perCat=zeros(3,Ncat-1);
Nb_SignifUnits=zeros(Ncat-1,1);
for cc=1:Ncat-1
    Units_Sel_perCat(1,cc)= sum((y(SemCellPV,cc)<0.01).*SelCat_PCCmax(SemCellPV,cc));
    Units_Sel_perCat(2,cc)=Units_Sel_perCat(1,cc)./size(SemCellPV,1);
    Units_Sel_perCat(3,cc)=Units_Sel_perCat(1,cc)./sum((y(SemCellPV,cc)<0.01));
    Nb_SignifUnits(cc)=sum((y(SemCellPV,cc)<0.01));
end
Units_Sel_perCatSU=zeros(3,Ncat-1);
for cc=1:Ncat-1
    Units_Sel_perCatSU(1,cc)= sum((y(SemSU,cc)<0.01).*SelCat_PCCmax(SemSU,cc));
    Units_Sel_perCatSU(2,cc)=Units_Sel_perCatSU(1,cc)./size(SemSU,1);
    Units_Sel_perCatSU(3,cc)=Units_Sel_perCatSU(1,cc)./sum((y(SemSU,cc)<0.01));
end
figure()
bar(Units_Sel_perCat(2,:),'k')
set(gca,'XTickLabel',StimTypeCM)
ylabel('Proportion of selective cells'); 

figure()
bar(Units_Sel_perCatSU(2,:),'k')
set(gca,'XTickLabel',StimTypeCM)
ylabel('Proportion of selective single units');

% Test if the number of units selective for the different call type is
% different
[p,Chi2stat,Df]=chi2test([Units_Sel_perCat(1,:);Nb_SignifUnits'-Units_Sel_perCat(1,:)])%highly significant
[p,Chi2stat,Df]=chi2test([Units_Sel_perCatSU(1,:);Nb_SignifSUnits'-Units_Sel_perCatSU(1,:)])%highly significant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Kruskall wallis on selectivity values effect of call type

KW_dataSU=nan(length(SemSU).*9,1);
KW_groupSU=cell(size(KW_dataSU));
LRI_mean = nan(1,9);
LRI_se = nan(1,9);
bb=0;
for cc=1:9
    LRI_cc_SU=LRI(intersect(SemSU, Signif_PCC{cc}),cc);
    LRI_mean(cc)=mean(LRI_cc_SU(isfinite(LRI_cc_SU)));
    LRI_se(cc)=std(LRI_cc_SU(isfinite(LRI_cc_SU)))./sqrt(sum(isfinite(LRI_cc_SU)));
    for uu=1:length(intersect(SemSU, Signif_PCC{cc}))
        bb=bb+1;
        KW_dataSU(bb)=LRI_cc_SU(uu);
        KW_groupSU{bb}=StimTypeCM{cc};
    end
end
KW_dataSU=KW_dataSU(1:bb);
KW_groupSU=KW_groupSU(1:bb);
[p,tbl,stats]=kruskalwallis(KW_dataSU,KW_groupSU)%significant
c=multcompare(stats, 'alpha',0.01)

%%%%%%%%%%%%%%%%%%%%%%%Fig Discrimination, Selectivity,
%%%%%%%%%%%%%%%%%%%%%%%Invariance%%%%%%%%%%%%%%
figure()
bar(LRI_mean,'k')
hold on
errorb(LRI_mean, LRI_se)
ylabel('Index of Selectivity Single SS units significantly discriminant')
set(gca,'XTickLabel',StimTypeCM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the mean +/- sd of invariance values for semantic cells
II_95 = quantile(reshape(II2, size(II2,1).*size(II2,2),1),0.95);
II_95 = quantile(II2_SemSUDisc,0.95)
II_mean2 = nan(2,9);
II_extremes = nan(2,9);
II_PercUnitsAb75 = nan(1,9);
II_se2 = II_mean2;
KW_data=nan(length(SemCellPV).*9,1);
KW_group=cell(size(KW_data));
KW_dataSU=nan(length(SemCellPV).*9,1);
KW_groupSU=cell(size(KW_data));
aa=0;
bb=0;
for cc=1:9
    II_mean2(1,cc)=nanmean(II2(intersect(SemCellPV, Signif_PCC{cc}),cc));
    II_mean2(2,cc)=nanmean(II2(intersect(SemSU, Signif_PCC{cc}),cc));
    II_se2(1,cc) = nanstd(II2(intersect(SemCellPV, Signif_PCC{cc}),cc))./sqrt(sum(~isnan(intersect(SemCellPV, Signif_PCC{cc}))));
    II_se2(2,cc) = nanstd(II2(intersect(SemSU, Signif_PCC{cc}),cc))./sqrt(sum(~isnan(intersect(SemSU, Signif_PCC{cc}))));
    II_PercUnitsAb75(cc) = sum(II2(intersect(SemCellPV, Signif_PCC{cc}))>0.8);%./length(intersect(SemCellPV, Signif_PCC{cc}));
    II_extremes(1,cc) = nanmin(II2(intersect(SemCellPV, Signif_PCC{cc}),cc));
    II_extremes(2,cc) = nanmax(II2(intersect(SemCellPV, Signif_PCC{cc}),cc));
    II2_cc=II2(intersect(SemCellPV, Signif_PCC{cc}),cc);
    II2_cc_SU=II2(intersect(SemSU, Signif_PCC{cc}),cc);
    for uu=1:length(intersect(SemCellPV, Signif_PCC{cc}))
        aa=aa+1;
        KW_data(aa)=II2_cc(uu);
        KW_group{aa}=StimTypeCM{cc};
    end
    for uu=1:length(intersect(SemSU, Signif_PCC{cc}))
        bb=bb+1;
        KW_dataSU(bb)=II2_cc_SU(uu);
        KW_groupSU{bb}=StimTypeCM{cc};
    end
end
KW_data=KW_data(1:aa);
KW_group=KW_group(1:aa);
KW_dataSU=KW_dataSU(1:bb);
KW_groupSU=KW_groupSU(1:bb);


%%%%%%%%%%%%%%%%%%%%%%%Fig Discrimination, Selectivity,
%%%%%%%%%%%%%%%%%%%%%%%Invariance%%%%%%%%%%%%%%
figure()
bar(II_mean2(1,:),'k')
hold on
errorb(II_mean2(1,:), II_se2(1,:))
ylim([0.4 0.7])
ylabel('Index of Invariance All SS units significantly Discriminant')

figure()
bar(II_mean2(2,:),'k')
hold on
errorb(II_mean2(2,:), II_se2(2,:))
ylim([0 0.4])
ylabel('Index of Invariance Single SS units significantly discriminant')
set(gca,'XTickLabel',StimTypeCM)

%test for difference of invariance between call types
[p,tbl,stats]=kruskalwallis(KW_data,KW_group)%significant
c=multcompare(stats, 'alpha',0.01)
[p,tbl,stats]=kruskalwallis(KW_dataSU,KW_groupSU)%significant
c=multcompare(stats, 'alpha',0.01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Corelation between PCC, LRI and Invariance of Semantic cells for
%%%%%%%% each call category
for cc=1:9
    StimTypeCM(cc)
    corrcoef(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc),LRI(intersect(SemSU, Signif_PCC{cc}),cc))
    [rho,pval]=corr(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc),LRI(intersect(SemSU, Signif_PCC{cc}),cc),'type','Spearman')
    if cc==4
        Indices_Local=intersect(SemSU, Signif_PCC{cc});%For WhiBlu5396M units, distress call invariance does not mean much since there is only one distress call in the matrix, these lines get rid of these units
        Indices_Local = intersect(Indices_Local, find(~isnan(II2(intersect(SemSU, Signif_PCC{cc}),cc))));
        corrcoef(PCC_inv(Indices_Local,cc),II2(Indices_Local,cc))
        [rho,pval]=corr(PCC_inv(Indices_Local,cc),II2(Indices_Local,cc),'type','Spearman')
        corrcoef(LRI(Indices_Local,cc),II2(Indices_Local,cc))
        [rho,pval]=corr(LRI(Indices_Local,cc),II2(Indices_Local,cc),'type','Spearman')
    else
        corrcoef(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc),II2(intersect(SemSU, Signif_PCC{cc}),cc))
        [rho,pval]=corr(PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc),II2(intersect(SemSU, Signif_PCC{cc}),cc),'type','Spearman')
        corrcoef(LRI(intersect(SemSU, Signif_PCC{cc}),cc),II2(intersect(SemSU, Signif_PCC{cc}),cc))
        [rho,pval]=corr(LRI(intersect(SemSU, Signif_PCC{cc}),cc),II2(intersect(SemSU, Signif_PCC{cc}),cc),'type','Spearman')
    end
    pause
end
PCC_effectif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




PCC_cat_zoneHS = nan(NZ-1,Ncat);
PCC_cat_zoneHSse = nan(NZ-1,Ncat);
PCC_cat_zoneSS = nan(NZ-1,Ncat);
PCC_cat_zoneSSse = nan(NZ-1,Ncat);
PCC_cat_zone = nan(NZ-1,Ncat);
PCC_cat_zonese = nan(NZ-1,Ncat);
for zz=1:(NZ-1)%only look at real region ("0" is unknown)
    Unitinzone = find(ZONES(:,1) == zz);
    UnitinzoneSSU = intersect(Unitinzone, SemCell);
    UnitinzoneSS = intersect(Unitinzone, SemCellPV);
    PCC_cat_zoneHS(zz,:) = mean(PCC_inv(UnitinzoneSSU,:),1);
    PCC_cat_zoneHSse(zz,:) = std(PCC_inv(UnitinzoneSSU,:))./sqrt(length(UnitinzoneSSU));
    PCC_cat_zoneSS(zz,:) = mean(PCC_inv(UnitinzoneSS,:),1);
    PCC_cat_zoneSSse(zz,:) = std(PCC_inv(UnitinzoneSS,:))./sqrt(length(UnitinzoneSS));
    PCC_cat_zone(zz,:) = mean(PCC_inv(Unitinzone,:),1);
    PCC_cat_zonese(zz,:) = std(PCC_inv(Unitinzone,:))./sqrt(length(Unitinzone));
end



% Plot for each type call the PCC values per category vs the LRI for that category for HS cells with regions color coded
LRI = cell2mat(Selectivity.DiagLRI')';
LRI_max = max(LRI');
InfInd = find(LRI_max==Inf);
maxLRI_max = max(LRI_max([1:(InfInd-1),(InfInd+1):end]));
LRI_max(InfInd)=maxLRI_max+1;


% Plot for semantic cells index of category discrimination (mean PCC or max
% PCC) vs index of selectivity (diag entropy or max LRI) and z axis regions
% and then for each cell the max PCC vs the LRI for the category that is
% best classified.
%find indices of category best classified for each cell
LRI_maxPCCcat = nan(NU,1);
for uu=1:NU
    Ind = find(PCC(:,uu)==max(PCC(1:(end-1),uu)));
    LRI_maxPCCcat(uu) = LRI(uu,Ind(1));
end
InfInd=find(LRI_maxPCCcat==Inf);
maxLRI_maxPCCcat = max(LRI_maxPCCcat([1:(InfInd-1),(InfInd+1):end]));
LRI_maxPCCcat(InfInd)=maxLRI_maxPCCcat+1;


% Plot bar graphs of selectivity per region for semantic cells
Mean_LRI_max_perzone=nan(NZ,2);
SE_LRI_max_perzone = nan(NZ,2);
Mean_DiagEnt_perzone = nan(NZ,2);
SE_DiagEnt_perzone = nan(NZ,2);

for zz = 1:NZ
    Unitinzone = find(ZONES(:,1) == UZ(zz));
    UnitinzoneSSU = intersect(Unitinzone, SemCell);
    UnitinzoneSS = intersect(Unitinzone, SemCellPV);
    Mean_LRI_max_perzone(zz,1) = mean(LRI_max(UnitinzoneSSU));
    SE_LRI_max_perzone(zz,1) = std(LRI_max(UnitinzoneSSU))./length(UnitinzoneSSU);
    Mean_DiagEnt_perzone(zz,1) = mean(Selectivity.DiagEntropyNorm(UnitinzoneSSU));
    SE_DiagEnt_perzone(zz,1) = std(Selectivity.DiagEntropyNorm(UnitinzoneSSU))./length(UnitinzoneSSU);
    Mean_LRI_max_perzone(zz,2) = mean(LRI_max(UnitinzoneSS));
    SE_LRI_max_perzone(zz,2) = std(LRI_max(UnitinzoneSS))./length(UnitinzoneSS);
    Mean_DiagEnt_perzone(zz,2) = mean(Selectivity.DiagEntropyNorm(UnitinzoneSS));
    SE_DiagEnt_perzone(zz,2) = std(Selectivity.DiagEntropyNorm(UnitinzoneSS))./length(UnitinzoneSS);
end

% Look at LRI per region per call category

LRI_cat_zoneHS = nan(NZ-1,(Ncat-1));%Here LRI is calculated only for vocalization cat (not for BG)
LRI_cat_zoneHSse = nan(NZ-1,(Ncat-1));
LRI_cat_zoneSS = nan(NZ-1,Ncat-1);
LRI_cat_zoneSSse = nan(NZ-1,Ncat-1);

for zz=1:(NZ-1)%only look at real region ("0" is unknown)
    Unitinzone = find(ZONES(:,1) == zz);
    UnitinzoneSSU = intersect(Unitinzone, SemCell);
    UnitinzoneSS = intersect(Unitinzone, SemCellPV);
    LRI_cat_zoneHS(zz,:) = mean(LRI(UnitinzoneSSU,:),1);
    LRI_cat_zoneHSse(zz,:) = std(LRI(UnitinzoneSSU,:))./sqrt(length(UnitinzoneSSU));
    LRI_cat_zoneSS(zz,:) = mean(LRI(Unitinzone,:),1);
    LRI_cat_zoneSSse(zz,:) = std(LRI(Unitinzone,:))./sqrt(length(Unitinzone));
end

% Look at the selectivity measure for the different significantly semantic
% cell type: linear/non linear
 %first replot figure 38 now looking at the cell linearity
 %first need to construct a vector of cell types. Code this as 1:linear
 %semantic cells, 2:non-linear semantic cells, 3:linear non semantic cells,
 %4:non-linear non semantic cells
 % Semantic and linear cells (gradient colored cells in figure12)
SemLN=find((AS.pAS_A<0.05).*(AS.pAS_S<0.05));
% Semantic and non-linear cells (Red cells in figure 12)
SemNLN=find((AS.pAS_A<0.05).*(AS.pAS_S>=0.05));%largest vector
% Non Semantic Cells
%AS.pAS_A>0.05
% Non Semantic acoustic cells (White cells in figure 12)
NonSemAc=find((AS.pAS_A>=0.05).*(AS.pAS_S<0.05));
% Non Semantic poorly acoustic cells (grey cells in figure 12)
NonSemNA=find((AS.pAS_A>=0.05).*(AS.pAS_S>=0.05));
Linearity = (AS.pAS_A<0.05).*(AS.pAS_S<0.05) + (AS.pAS_A<0.05).*(AS.pAS_S>=0.05).*2 + (AS.pAS_A>=0.05).*(AS.pAS_S<0.05).*3 + (AS.pAS_A>=0.05).*(AS.pAS_S>=0.05).*4;


% Plot bar graphs of selectivity, average discrimination (AV_PCC) and max discrimination level (max_PCC) per type of cells for semantic cells
Mean_LRI_max_perCT=nan(4,2);
SE_LRI_max_perCT = nan(4,2);
Mean_DiagEnt_perCT = nan(4,2);
SE_DiagEnt_perCT = nan(4,2);
Mean_PCC_max_perCT=nan(4,2);
SE_PCC_max_perCT = nan(4,2);
Mean_PCC_av_perCT=nan(4,2);
SE_PCC_av_perCT = nan(4,2);

for zz = 1:4
    UnitinCT = find(Linearity == zz);
    UnitinCTHS = intersect(UnitinCT, SemCell);
    UnitinCTSS = intersect(UnitinCT, SemCellPV);
    Mean_LRI_max_perCT(zz,1) = mean(LRI_max(UnitinCTHS));
    SE_LRI_max_perCT(zz,1) = std(LRI_max(UnitinCTHS))./length(UnitinCTHS);
    Mean_DiagEnt_perCT(zz,1) = mean(Selectivity.DiagEntropyNorm(UnitinCTHS));
    SE_DiagEnt_perCT(zz,1) = std(Selectivity.DiagEntropyNorm(UnitinCTHS))./length(UnitinCTHS);
    Mean_LRI_max_perCT(zz,2) = mean(LRI_max(UnitinCTSS));
    SE_LRI_max_perCT(zz,2) = std(LRI_max(UnitinCTSS))./length(UnitinCTSS);
    Mean_DiagEnt_perCT(zz,2) = mean(Selectivity.DiagEntropyNorm(UnitinCTSS));
    SE_DiagEnt_perCT(zz,2) = std(Selectivity.DiagEntropyNorm(UnitinCTSS))./length(UnitinCTSS);
    Mean_PCC_max_perCT(zz,1)=mean(PCC_max(UnitinCTHS));
    SE_PCC_max_perCT (zz,1)= std(PCC_max(UnitinCTHS))./length(UnitinCTHS);
    Mean_PCC_av_perCT(zz,1)=mean(PCC_AV(UnitinCTHS));
    SE_PCC_av_perCT (zz,1)= std(PCC_AV(UnitinCTHS))./length(UnitinCTHS);
    Mean_PCC_max_perCT(zz,2)=mean(PCC_max(UnitinCTSS));
    SE_PCC_max_perCT(zz,2) = std(PCC_max(UnitinCTSS))./length(UnitinCTHS);
    Mean_PCC_av_perCT(zz,2)=mean(PCC_AV(UnitinCTSS));
    SE_PCC_av_perCT(zz,2) = std(PCC_AV(UnitinCTSS))./length(UnitinCTSS);
end 

% Plot in Discrimination and selectivity space the name of subject
figure(50)
for cc=1:(Ncat-1)
    ss=subplot(3,3,cc)
    gg=gscatter(PCC_inv(:,cc), LRI(:,cc), SUBJ(:,1), 'mrgckb', '.....',[20 20 20 20 20 20])
    title(ss,sprintf('LRI (y) vs PCC (x) %s',StimTypeCM{cc}))
    axis([0 1 -6 4])
    LL1=line([1/9 1/9], [-6 4])
    LL2=line([0 1], [1 1])
    set(LL2, 'Color', 'r','LineWidth', 2)
    LL2=line([0 1], [-1 -1])
    set(LL2, 'Color', 'r', 'LineWidth', 2)
    LL2=line([0 1], [0 0])
    set(LL2, 'Color', 'r')
end

% Plot in Discrimination and selectivity space the sex of subject
figure(51)
for cc=1:(Ncat-1)
    ss=subplot(3,3,cc)
    gg=gscatter(PCC_inv(:,cc), LRI(:,cc), SUBJ(:,2), 'mb', '..',[20 20])
    title(ss,sprintf('LRI (y) vs PCC (x) %s',StimTypeCM{cc}))
    axis([0 1 -6 4])
    LL1=line([1/9 1/9], [-6 4])
    LL2=line([0 1], [1 1])
    set(LL2, 'Color', 'r','LineWidth', 2)
    LL2=line([0 1], [-1 -1])
    set(LL2, 'Color', 'r', 'LineWidth', 2)
    LL2=line([0 1], [0 0])
    set(LL2, 'Color', 'r')
end

% Plot in Discrimination and selectivity space the spike shape
figure(52)

for cc=1:(Ncat-1)
    ss=subplot(3,3,cc)
    %gg=gscatter(PCC_inv(:,cc), LRI(:,cc), Spike_shape, 'kbgcm', '.....',[20 20 20 20 20])
    gg=gscatter(PCC_inv(MU,cc), LRI(MU,cc), Spike_shape(MU), 'kbgcm', '.....',[20 20 20 20 20])
    title(ss,sprintf('LRI (y) vs PCC (x) %s',StimTypeCM{cc}))
    axis([0 1 -6 4])
    LL1=line([1/9 1/9], [-6 4])
    LL2=line([0 1], [1 1])
    set(LL2, 'Color', 'r','LineWidth', 2)
    LL2=line([0 1], [-1 -1])
    set(LL2, 'Color', 'r', 'LineWidth', 2)
    LL2=line([0 1], [0 0])
    set(LL2, 'Color', 'r')
end
%% Explore the effect of single vs multiunits on all the variables
% Effect on semantic info
[H,P] = ttest2(SemanticIndex.Observed(MU), SemanticIndex.Observed(SU))
ranksum(SemanticIndex.Observed(MU), SemanticIndex.Observed(SU))%significant test on all units
fprintf(1,'Semantic info of all %d Multi-units:%f+/-%f\n',length(MU),mean(SemanticIndex.Observed(MU)), std(SemanticIndex.Observed(MU)));
fprintf(1,'Semantic info of all %d single-units:%f+/-%f\n',length(SU),mean(SemanticIndex.Observed(SU)), std(SemanticIndex.Observed(SU)));
ranksum(SemanticIndex.Observed(intersect(MU,SemCellPV)), SemanticIndex.Observed(intersect(SU,SemCellPV)))%significant test on semantic units
fprintf(1,'Semantic info of the %d semantic Multi-units:%f+/-%f\n',length(intersect(MU,SemCellPV)),mean(SemanticIndex.Observed(intersect(MU,SemCellPV))), std(SemanticIndex.Observed(intersect(MU,SemCellPV))));
fprintf(1,'Semantic info of the %d semantic single-units:%f+/-%f\n',length(intersect(SU,SemCellPV)),mean(SemanticIndex.Observed(intersect(SU,SemCellPV))), std(SemanticIndex.Observed(intersect(SU,SemCellPV))));


% Effect on PCC_mean
ranksum(PCC_mean(SemMU), PCC_mean(SemSU))
fprintf(1,'Mean discrimination (PCC_mean) of the %d semantic Multi-units:%f+/-%f\n',length(SemMU),mean(PCC_mean(SemMU)), std(PCC_mean(SemMU)));
fprintf(1,'Mean discrimination (PCC_mean) of the %d semantic single-units:%f+/-%f\n',length(SemSU),mean(PCC_mean(SemSU)), std(PCC_mean(SemSU)));

% Effect on PCC, LRI Spike rate and invariance for each call category
for cc=1:(Ncat-1)
    PCC_su=PCC_inv(intersect(SemSU, Signif_PCC{cc}),cc);
    PCC_mu=PCC_inv(intersect(SemMU, Signif_PCC{cc}),cc);
    ranksum(PCC_mu, PCC_su)
    fprintf(1,'Discrimination (PCC) for %s of the %d semantic Multi-units:%f+/-%f\n',StimTypeCM{cc},length(PCC_mu),mean(PCC_mu), std(PCC_mu));
    fprintf(1,'Discrimination (PCC) for %s of the %d semantic single-units:%f+/-%f\n',StimTypeCM{cc},length(PCC_su),mean(PCC_su), std(PCC_su));
    LRI_su=LRI(intersect(SemSU, Signif_PCC{cc}),cc);
    LRI_mu=LRI(intersect(SemMU, Signif_PCC{cc}),cc);
    ranksum(LRI_mu, LRI_su)
    fprintf(1,'Selectivity (LRI) for %s of the %d semantic Multi-units:%f+/-%f\n',StimTypeCM{cc},length(LRI_mu),mean(LRI_mu), std(LRI_mu));
    fprintf(1,'Selectivity (LRI) for %s of the %d semantic single-units:%f+/-%f\n',StimTypeCM{cc},length(LRI_su),mean(LRI_su), std(LRI_su));
    SR_su=MeanSR.SpikeRate_Cat(intersect(SemSU, Signif_PCC{cc}),cc);
    SR_mu=MeanSR.SpikeRate_Cat(intersect(SemMU, Signif_PCC{cc}),cc);
    ranksum(SR_mu, SR_su)
    fprintf(1,'Spike Rate for %s of the %d semantic Multi-units:%f+/-%f\n',StimTypeCM{cc},length(SR_mu),mean(SR_mu), std(SR_mu));
    fprintf(1,'Spike Rate for %s of the %d semantic single-units:%f+/-%f\n',StimTypeCM{cc},length(SR_su),mean(SR_su), std(SR_su));
    II_su=II(intersect(SemSU, Signif_PCC{cc}),cc);
    II_mu=II(intersect(SemMU, Signif_PCC{cc}),cc);
    ranksum(II_mu, II_su)
    fprintf(1,'Invariance for %s of the %d semantic Multi-units:%f+/-%f\n',StimTypeCM{cc},length(II_mu),mean(II_mu), std(II_mu));
    fprintf(1,'Invariance for %s of the %d semantic single-units:%f+/-%f\n',StimTypeCM{cc},length(II_su),mean(II_su), std(II_su));
end

%% LOGISTIC ANALYSIS Find cells that are highly discriminative (PCC>0.5) and highly selective
% (LSI>3)
PCC_units = PCC_inv(:,1:9);
Cell4Logistic = cell(9,1);
Ncat = size(PCC_units,2);
for cc=1:(Ncat)
    Cell4Logistic_local = intersect(find(PCC_units(:,cc)>0.5), find(LRI(:,cc)>3));
    if ~isempty(Cell4Logistic_local)
        fprintf('There are %d interesting cells for %s\n', length(Cell4Logistic_local), StimTypeCM{cc});
        %check that these cells are in SemCellPV and Signif_PCC
        Cell4Logistic_local2 = intersect(intersect(SemCellPV, Signif_PCC{cc}),Cell4Logistic_local);
        if ~isempty(Cell4Logistic_local2);
            fprintf('%d out of %d are significantly discriminant semantic units\n', length(Cell4Logistic_local2), length(Cell4Logistic_local));
            Cell4Logistic{cc}=List_matfilepath(Cell4Logistic_local2);
        else
               fprintf('%d out of %d are significantly discriminant semantic units\n', length(Cell4Logistic_local2), length(Cell4Logistic_local));
        end
    else
        fprintf('There are NO interesting cells for %s\n', StimTypeCM{cc});
    end
end
%MeanSpikeRate.Values =  MeanSR.SpikeRate_Voc;
%MeanSpikeRate.TDT_StimNames = MeanSR.TDT_names;
%save('Data4Logistic.mat', 'Cell4Logistic','PCC_units', 'LRI', 'SemCellPV', 'Signif_PCC', 'List_matfilepath', 'StimTypeCM', 'MeanSpikeRate')


%% 3D BRAIN PLOTS
% Plot  a 3D representation of semantic information for all cells for each individual

Subjects=unique(SUBJ);
for subj=1:length(Subjects)
    subject=Subjects(subj)
    xyzCells=intersect(find(~isnan(DV_theo)), find(SUBJ==subject));
    if ~isempty(xyzCells)
        PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),0,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
%         if subject==3
            PlotSphere(LR(xyzCells), 4*RC(xyzCells), DV(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),0,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
% 
%         else
%             PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MeanSR(xyzCells),0);
%         end
        set(get(gca,'Title'), 'String', sprintf('%s Semantic Information', SUBJECTS{subject+1}));
    end
end

% Plot all the cells in an hypothetic brain
xyzCells = find(~isnan(DV_theo));
PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),1,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),0,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),1,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
xyzCells = find(~isnan(DV));
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),0,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),0,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),1,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_diag_uni_cat(xyzCells),1,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),PCC_mean(xyzCells),1,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),PCC_mean(xyzCells),1,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),2.^(LRI_max(xyzCells)),1,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),LRI_max(xyzCells),1,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),LRI_max(xyzCells),0,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),LRI_max(xyzCells),0,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fig4a 4b 4c Anatomy1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyzCells = find(~isnan(DV));%all cells for which we have the anatomical positions
PlotAnat2D(LR(xyzCells), RC(xyzCells), DV(xyzCells),PCC_mean(xyzCells),0,1,ZONES7(xyzCells),List_zones7,0,1);
PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),PCC_mean(xyzCells),0,1,ID_Sem_NonSem(xyzCells),{'NonSemantic' 'Semantic'},0,1);
SemAnat = intersect(SemCellPV, xyzCells);%all cells for which we have the anatomical positions + are significantly semantic
PlotSphere(LR(SemAnat), RC(SemAnat), DV(SemAnat),PCC_mean(SemAnat),0,1,ID_AS(SemAnat),{'+A+S' 'a+S' '+As' 'as'},0,1);

xyzCellsSU = intersect(SU,find(~isnan(DV)));%all single units for which we have the anatomical positions
PlotAnat2D(LR(xyzCellsSU), RC(xyzCellsSU), DV(xyzCellsSU),PCC_mean(xyzCellsSU),0,1,ZONES7(xyzCellsSU),List_zones7,0,1);
PlotSphere(LR(xyzCellsSU), RC(xyzCellsSU), DV(xyzCellsSU),PCC_mean(xyzCellsSU),0,1,ZONES7(xyzCellsSU),List_zones7,0,1);
PlotAnat2D(LR(xyzCellsSU), RC(xyzCellsSU), DV(xyzCellsSU),PCC_mean(xyzCellsSU),0,1,ID_Sem_NonSem(xyzCellsSU),{'NonSemantic' 'Semantic'},0,1);
PlotSphere(LR(xyzCellsSU), RC(xyzCellsSU), DV(xyzCellsSU),PCC_mean(xyzCellsSU),0,1,ID_Sem_NonSem(xyzCellsSU),{'NonSemantic' 'Semantic'},0,1);
SemAnatSU = intersect(SemSU, xyzCellsSU);%all single units for which we have the anatomical positions + are significantly semantic
PlotSphere(LR(SemAnatSU), RC(SemAnatSU), DV(SemAnatSU),PCC_mean(SemAnatSU),0,1,ID_AS(SemAnatSU),{'+A+S' 'a+S' '+As' 'as'},0,1);
PlotAnat2D(LR(SemAnatSU), RC(SemAnatSU), DV(SemAnatSU),PCC_mean(SemAnatSU),0,1,ID_AS(SemAnatSU),{'+A+S' 'a+S' '+As' 'as'},0,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



SemAnatDisc = intersect(SemAnat, Signif_PCC_Cell);% all cells that are anatomically localized, significantly semantic and discriminate at least one category above chance (all cells except one...)
PlotSphere(LR(SemAnatDisc), RC(SemAnatDisc), DV(SemAnatDisc),PCC_mean(SemAnatDisc),1,1,ZONES7(SemAnatDisc),List_zones7,0,1);
PlotSphere(LR(SemAnatDisc), RC(SemAnatDisc), DV(SemAnatDisc),LRI_max(SemAnatDisc),1,1,ZONES7(SemAnatDisc),List_zones7,0,1);

SemAnatDiscSU = intersect(SemAnatSU, Signif_PCC_Cell);% single cells that are anatomically localized, significantly semantic and discriminate at least one category above chance (all cells except one...)
PlotSphere(LR(SemAnatDiscSU), RC(SemAnatDiscSU), DV(SemAnatDiscSU),PCC_mean(SemAnatDiscSU),1,1,ZONES7(SemAnatDiscSU),List_zones7,0,1);
PlotSphere(LR(SemAnatDiscSU), RC(SemAnatDiscSU), DV(SemAnatDiscSU),LRI_max(SemAnatDiscSU),1,1,ZONES7(SemAnatDiscSU),List_zones7,0,1);

SemAnatDiscwoUn = setdiff(SemAnatDisc, find(ZONES7==7));
SemAnatDiscwoUnSU = setdiff(SemAnatDiscSU, find(ZONES7==7));
SemAnatwoUnSU = setdiff(SemAnatSU, find(ZONES7==7));
x=[PCC_mean(SemAnatDiscwoUn) LRI_max(SemAnatDiscwoUn)' II_PCCmax(SemAnatDiscwoUn)];
xSU=[PCC_mean(SemAnatDiscwoUnSU) LRI_max(SemAnatDiscwoUnSU)' II_PCCmax(SemAnatDiscwoUnSU)];
x2=PCC_mean(SemAnatDiscwoUn);
[d,p,stats]=manova1(x,ZONES7(SemAnatDiscwoUn));% significant at 2 dims
[d,p,stats]=manova1(xSU,ZONES7(SemAnatDiscwoUnSU));% significant at 1 dim
%[d,p,stats]=manova1(x,ZONES2(SemAnatDisc));

% Calculate and plot group stats
[mean_groups std_groups sem_groups]= grpstats(stats.canon(:,1:2),ZONES7(SemAnatDiscwoUn), {'mean', 'std', 'sem'} ); % groups along first two principal components

 % Plot groups along first 2 canonicals
figure()
hold on;
plot(mean_groups(:,1),mean_groups(:,2), 'k+');
hold on;
for ia = 1:length(List_zones7(1:end-1))
    plot( mean_groups(ia,1) + (sem_groups(ia,1).*cos([0:pi/10:2*pi])), mean_groups(ia,2) + (sem_groups(ia,2).*sin([0:pi/10:2*pi])),'k');
    text(mean_groups(ia,1)+0.1,mean_groups(ia,2)-.1, List_zones7(ia));
end
hold off;



[p,tbl,stats]=kruskalwallis(PCC_mean(SemAnatDiscwoUn),ZONES7(SemAnatDiscwoUn))%significant
c=multcompare(stats, 'alpha',0.01)%some signif diff
[p,tbl,stats]=kruskalwallis(LRI_max(SemAnatDiscwoUn),ZONES7(SemAnatDiscwoUn))%significant
c=multcompare(stats, 'alpha',0.01)%some signif diff

[p,tbl,stats]=kruskalwallis(PCC_mean(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%significant
c=multcompare(stats, 'alpha',0.01)%NS
[P,H,STATS]=ranksum(PCC_mean(intersect(SemAnatDiscwoUnSU, find(ZONES2==1))), PCC_mean(intersect(SemAnatDiscwoUnSU, find(ZONES2==2))))
[p,tbl,stats]=kruskalwallis(PCC_mean(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%significant
c=multcompare(stats, 'alpha',0.01)%Significant

[p,tbl,stats]=kruskalwallis(PCC_max(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%significant
c=multcompare(stats, 'alpha',0.01)%some significant
[p,tbl,stats]=kruskalwallis(PCC_max(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%significant
c=multcompare(stats, 'alpha',0.01)%Significant

[p,tbl,stats]=kruskalwallis(LRI_max(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%NS
c=multcompare(stats, 'alpha',0.01)%NS
[p,tbl,stats]=kruskalwallis(LRI_max(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%NS
c=multcompare(stats, 'alpha',0.01)%NS

[p,tbl,stats]=kruskalwallis(Selectivity.DiagEntropyNorm(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%NS
c=multcompare(stats, 'alpha',0.01)%NS
[p,tbl,stats]=kruskalwallis(Selectivity.DiagEntropyNorm(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%NS
c=multcompare(stats, 'alpha',0.01)%NS

[p,tbl,stats]=kruskalwallis(II_max(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%Significant
c=multcompare(stats, 'alpha',0.01)%NS
[p,tbl,stats]=kruskalwallis(II_max(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%Significant
c=multcompare(stats, 'alpha',0.01)%significant

[p,tbl,stats]=kruskalwallis(II_mean(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%Significant
c=multcompare(stats, 'alpha',0.01)%NS
[p,tbl,stats]=kruskalwallis(II_mean(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%Significant
c=multcompare(stats, 'alpha',0.01)%significant

[p,tbl,stats]=kruskalwallis(II_Mean_disccat(SemAnatDiscwoUnSU),ZONES7(SemAnatDiscwoUnSU))%NS
c=multcompare(stats, 'alpha',0.01)%NS
[P,H,STATS]=ranksum(II_Mean_disccat(intersect(SemAnatDiscwoUnSU, find(ZONES2==1))), II_Mean_disccat(intersect(SemAnatDiscwoUnSU, find(ZONES2==2))))%NS
[p,tbl,stats]=kruskalwallis(II_Mean_disccat(SemAnatDiscwoUnSU),ZONES2(SemAnatDiscwoUnSU))%NS
c=multcompare(stats, 'alpha',0.01)%NS



% Run anovan to test both call type and zones on the different variables
PCC_anova = reshape(PCC_inv(SemAnatDiscwoUnSU,1:9),9.*length(SemAnatDiscwoUnSU),1);
LRI_anova = reshape(LRI(SemAnatDiscwoUnSU,1:9),9.*length(SemAnatDiscwoUnSU),1);
InfInd = find(LRI_anova==-Inf);
minLRI_anova = min(LRI_anova(setdiff(1:length(LRI_anova),InfInd)));
LRI_anova(InfInd)=minLRI_anova-1;
II_anova = reshape(II(SemAnatDiscwoUnSU,1:9),9.*length(SemAnatDiscwoUnSU),1);
ZONES7_anova = repmat(ZONES7(SemAnatDiscwoUnSU),9,1);
CallType_anova = reshape(repmat(StimTypeCM(1:9)',length(SemAnatDiscwoUnSU),1),9.*length(SemAnatDiscwoUnSU),1);
[p,tbl,stats]=anovan(PCC_anova, {ZONES7_anova, CallType_anova}, 'model', 'full')%each var and interaction significant 
c=multcompare(stats,'dim',[1 2])%some significant

% [p,tbl,stats]=anovan(LRI_anova, {ZONES7_anova, CallType_anova}, 'model', 'full')%calltype and interaction significant
% c=multcompare(stats,'dim',[1 2])
% 
% [p,tbl,stats]=anovan(II_anova, {ZONES7_anova, CallType_anova}, 'model', 'full')%all var and interaction significant
% c=multcompare(stats,'dim',[1 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fig AnatCallType Matrices of PCC,LRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%and II mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%values%%%%%%%%%%%%%%%%%%%%%
% Calculate and look effect of zones and calltype on LRI and II with anovan for discriminant units only
LRI_anova_disc=cell(Ncat-1,1);
II_anova_disc=LRI_anova_disc;
PCC_anova_disc=LRI_anova_disc;
ZONES7_anova_disc=LRI_anova_disc;
ZONES7_KW_perCT=LRI_anova_disc;
ZONES2_KW_perCT=LRI_anova_disc;
CallType_anova_disc=LRI_anova_disc;
PCC_KW_perCT = LRI_anova_disc;
for cc=1:(Ncat-1)
    LRI_anova_disc{cc}=LRI(intersect(SemAnatDiscwoUnSU,Signif_PCC{cc}),cc);
    ZONES7_anova_disc{cc}=ZONES7(intersect(SemAnatDiscwoUnSU,Signif_PCC{cc}));
    ZONES7_KW_perCT{cc}=ZONES7(SemAnatwoUnSU);
    ZONES2_KW_perCT{cc}=ZONES2(SemAnatwoUnSU);
    II_anova_disc{cc}=II2(intersect(SemAnatDiscwoUnSU,Signif_PCC{cc}),cc);
    PCC_anova_disc{cc}=PCC_inv(intersect(SemAnatDiscwoUnSU,Signif_PCC{cc}),cc);
    PCC_KW_perCT{cc}=PCC_inv(SemAnatwoUnSU,cc);
    CallType_anova_disc{cc} = repmat(cc,length(intersect(SemAnatDiscwoUnSU,Signif_PCC{cc})),1);
%     fprintf(1,'%s: %s',StimTypeCM{cc},'Discrimination PCC all SU')
%     [p,tbl,stats]=kruskalwallis(PCC_KW_perCT{cc},ZONES7_KW_perCT{cc})% NS
%     c=multcompare(stats, 'alpha',0.01)%NS
%     [P,H,STATS]=ranksum(PCC_KW_perCT{cc}(find(ZONES2_KW_perCT{cc}==1)),PCC_KW_perCT{cc}(find(ZONES2_KW_perCT{cc}==2)))
%     pause
%     fprintf(1,'%s: %s',StimTypeCM{cc},'Discrimination PCC')
%     [p,tbl,stats]=kruskalwallis(PCC_anova_disc{cc},ZONES7_anova_disc{cc})% NS
%     c=multcompare(stats, 'alpha',0.01)%NS
%     pause
%     fprintf(1,'%s: %s',StimTypeCM{cc},'Selectivity LRI')
%     [p,tbl,stats]=kruskalwallis(LRI_anova_disc{cc},ZONES7_anova_disc{cc})% only significant for DC
%     c=multcompare(stats, 'alpha',0.01)%for DC CMM>L1
%     pause
    fprintf(1,'%s: %s',StimTypeCM{cc},'Invariance II')
    [p,tbl,stats]=kruskalwallis(II_anova_disc{cc},ZONES7_anova_disc{cc})%NS
    c=multcompare(stats, 'alpha',0.01)%NS
    pause
end
LRI_anova_disc=cell2mat(LRI_anova_disc);
II_anova_disc=cell2mat(II_anova_disc);
ZONES7_anova_disc=cell2mat(ZONES7_anova_disc);
CallType_anova_disc=cell2mat(CallType_anova_disc);
[p,tbl,stats]=anovan(LRI_anova_disc, {ZONES7_anova_disc, CallType_anova_disc}, 'model', 'full')%calltype and zones significant, tendancy for interaction
c=multcompare(stats,'dim',[1 2])
[p,tbl,stats]=anovan(II_anova_disc, {ZONES7_anova_disc, CallType_anova_disc}, 'model', 'full')%call type significant, tendancy for zones and interaction NS
c=multcompare(stats,'dim',[1 2])




% Loop through regions and plot the average values of PCC, LRI and II per
% call type
ZZ=unique(ZONES7(SemAnatDiscwoUnSU));
PCC_CallType_Zones_mean=nan((Ncat-1),length(ZZ));
PCC_CallType_Zones_sd=PCC_CallType_Zones_mean;
PCC_CallType_Zones_mean2=PCC_CallType_Zones_mean;
PCC_CallType_Zones_sd2=PCC_CallType_Zones_mean;
LRI_CallType_Zones_mean=PCC_CallType_Zones_mean;
LRI_CallType_Zones_sd=PCC_CallType_Zones_mean;
II_CallType_Zones_mean=PCC_CallType_Zones_mean;
II_CallType_Zones_sd=PCC_CallType_Zones_mean;
Min_LRI=min(min(LRI.*(isfinite(LRI))));
for zz=1:length(ZZ)
    zzi = ZZ(zz);
    Uzones = intersect(SemAnatDiscwoUnSU,find(ZONES7==zzi));
    for cc=1:(Ncat-1)
        PCC_CallType_Zones_mean(cc,zzi)=nanmean(PCC_inv(Uzones,cc));
        PCC_CallType_Zones_sd(cc,zzi)=nanstd(PCC_inv(Uzones,cc));
        UzonesDisc=intersect(Uzones, Signif_PCC{cc});% This is to investigate only the selectivity and invariance of discriminant cells
        PCC_CallType_Zones_mean2(cc,zzi)=nanmean(PCC_inv(UzonesDisc,cc));
        PCC_CallType_Zones_sd2(cc,zzi)=nanstd(PCC_inv(UzonesDisc,cc));
        LRI_local=LRI(UzonesDisc,cc);
        LRI_local(find(LRI_local==-Inf)) = Min_LRI;
        LRI_CallType_Zones_mean(cc,zzi)=nanmean(LRI_local);
        LRI_CallType_Zones_sd(cc,zzi)=nanstd(LRI_local);
        II_CallType_Zones_mean(cc,zzi)=nanmean(II2(UzonesDisc,cc));
        II_CallType_Zones_sd(cc,zzi)=nanstd(II2(UzonesDisc,cc));
    end
end

figure()
ss=subplot(1,3,1)
[B,IPCC]=sort(mean(PCC_CallType_Zones_mean,2),'descend');
imagesc(PCC_CallType_Zones_mean(IPCC,:))
colorbar
colormap('Hot')
Hotmap=colormap;
colormap(flip(Hotmap))
set(gca,'XTickLabel', List_zones7(1:6),'XTick',1:6)
set(gca,'YTickLabel', StimTypeCM(IPCC))
title('Discrimination (PCC) Semantic single units')
ss=subplot(1,3,2)
[B,ILRI]=sort(mean(LRI_CallType_Zones_mean,2),'descend');
imagesc(LRI_CallType_Zones_mean(ILRI,:))
colorbar
colormap('Hot')
Hotmap=colormap;
colormap(flip(Hotmap))
set(gca,'XTickLabel', List_zones7(1:6),'XTick',1:6)
set(gca,'YTickLabel', StimTypeCM(ILRI))
title('Selectivity (LRI) Semantic discriminant single units')
ss=subplot(1,3,3)
[B,III]=sort(mean(II_CallType_Zones_mean,2),'descend');
imagesc(II_CallType_Zones_mean(III,:))
colorbar
colormap('Hot')
Hotmap=colormap;
colormap(flip(Hotmap))
set(gca,'XTickLabel', List_zones7(1:6),'XTick',1:6)
set(gca,'YTickLabel', StimTypeCM(III))
title('Invariance (II) Semantic discriminant single units')

% Plot average value of PCC_mean fro each region
figure(13)
plot(1:6,Mean_PCC_AV_perzoneSU(1:6),'ks','MarkerFaceColor','k','MarkerSize',14);
hold on
errorb(Mean_PCC_AV_perzoneSU(1:6),SE_PCC_AV_perzoneSU(1:6));
xlim([0.5 6.5])
ylim([0.1 0.25])
set(gca,'XTick',1:6)
set(gca,'XTickLabel',List_zones7(1:6))
ylabel('Average PCC')
title('Semantic Single Units')

% Plot average value of Global selectivity for each zone
figure(11)
plot(1:6,Mean_GlobalSel_AV_perzoneSU(1:6),'ks','MarkerFaceColor','k','MarkerSize',14);
hold on
errorb(Mean_GlobalSel_AV_perzoneSU(1:6),SE_GlobalSel_AV_perzoneSU(1:6));
xlim([0.5 6.5])
ylim([0.3 0.4])
set(gca,'XTick',1:6)
set(gca,'XTickLabel',List_zones7(1:6))
ylabel('Average GlobalSelectivity')
title('Semantic Single Units')

% Plot average value of Average invariance for each zone
figure()
plot(1:6,Mean_II_Meandiscat(1:6),'ks','MarkerFaceColor','k','MarkerSize',14);
hold on
errorb(Mean_II_Meandiscat(1:6),SE_II_Meandiscat(1:6));
xlim([0.5 6.5])
%ylim([0.3 0.4])
set(gca,'XTick',1:6)
set(gca,'XTickLabel',List_zones7(1:6))
ylabel('Average Invariance Index')
title('Semantic Single Units')

% side view of 3D plot with Mean PCC for semantic cells
PlotAnat2D(LR(SemAnatSU), RC(SemAnatSU), DV(SemAnatSU),PCC_mean(SemAnatSU),1,1,ZONES7(SemAnatSU),List_zones7,0,1,16);
set(gca,'ZLim',[-3000 0])
%side view of 3D plot with mean(II) of discriminated categories for semantic single cells
IINan = find(isnan(II_Mean_disccat));
SemAnatSU_NN=setdiff(SemAnatSU,IINan);
PlotAnat2D(LR(SemAnatSU_NN), RC(SemAnatSU_NN), DV(SemAnatSU_NN),II_Mean_disccat(SemAnatSU_NN,1),1,1,ZONES7(SemAnatSU_NN),List_zones7,0,1,16);
set(gca,'ZLim',[-3000 0])
%%%%%%%%%%% side view of 3D plot with only selective units in color and the type of call they
% are selective to. Non-selective units are grey colored. Size of ball =
% global selectivity
LRI_sel_CT=PCC_max_CT;
LRI_sel_CT(find(LRI_max<=Sel_Thresh))=10;
PlotAnat2D(LR(SemAnatSU), RC(SemAnatSU), DV(SemAnatSU),Selectivity.DiagEntropyNorm(SemAnatSU),1,1, LRI_sel_CT(SemAnatSU),{StimTypeCM{1:end-1} 'Non-selective'},0,1,10);
set(gca,'ZLim',[-3000 0])


% side view of 3D plot with LRI for semantic cells discriminating DC
SemAnatSU_DiscDC = intersect(SemAnatSU,Signif_PCC{3});
NANInd = find(isnan(LRI(SemAnatSU_DiscDC,3)));
SemAnatSU_DiscDC = setdiff(SemAnatSU_DiscDC,NANInd);
InfInd = ~isfinite(LRI(SemAnatSU_DiscDC,3));
NonInfInd = isfinite(LRI(SemAnatSU_DiscDC,3));
LRI_local = LRI(SemAnatSU_DiscDC,3);
MaxLRI=max(LRI_local(NonInfInd));
LRI_local(InfInd)=2*MaxLRI;
PlotAnat2D(LR(SemAnatSU_DiscDC), RC(SemAnatSU_DiscDC), DV(SemAnatSU_DiscDC),LRI_local,1,1,ZONES7(SemAnatSU_DiscDC),List_zones7,0,1,16);
set(gca,'ZLim',[-3000 0])
%side view of 3D plot with II for semantic cells discriminating Agg
SemAnatSU_DiscAgg = intersect(SemAnatSU,Signif_PCC{1});
PlotAnat2D(LR(SemAnatSU_DiscAgg), RC(SemAnatSU_DiscAgg), DV(SemAnatSU_DiscAgg),II2(SemAnatSU_DiscAgg,1),1,1,ZONES7(SemAnatSU_DiscAgg),List_zones7,0,1,20);
set(gca,'ZLim',[-3000 0])
set(gca,'YLim',[-1500 1500])
% side view of 3D plot with only selective units and the type of call they
% are selective to
SemAnatDiscSelSU=intersect(SemAnatDiscSU, find(LRI_max>Sel_Thresh));
PlotAnat2D(LR(SemAnatDiscSelSU), RC(SemAnatDiscSelSU), DV(SemAnatDiscSelSU),PCC_max(SemAnatDiscSelSU),1,1,PCC_max_CT(SemAnatDiscSelSU),StimTypeCM(1:end-1),0,1);
set(gca,'ZLim',[-3000 0])

% side view of 3D plot with only selective units in color and the type of call they
% are selective to. Non-selective units are grey colored. Size of ball =
% PCC
PlotAnat2D(LR(SemAnatSU), RC(SemAnatSU), DV(SemAnatSU),PCC_max(SemAnatSU),1,1, LRI_sel_CT(SemAnatSU),{StimTypeCM{1:end-1} 'Non-selective'},0,1,16);
set(gca,'ZLim',[-3000 0])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%Proportion of discriminant, selective and invariant single units per
%%%%%%%%%%%%%%%%%%%%%%%region%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
subplot(1,3,1)
imagesc(PercDiscSSUCells_per7area(1:6,:)')
colorbar
set(gca,'XTickLabel', List_zones7(1:6),'XTick',1:6)
set(gca,'YTickLabel', StimTypeCM(1:9))
xlabel('anatomical zones');
ylabel('Call Type')
title('Proportion of semantic units Discriminant')
subplot(1,3,2)
imagesc(PercSelDiscSSUCells_per7area(1:6,:)')
colorbar
set(gca,'XTickLabel', List_zones7(1:6),'XTick',1:6)
set(gca,'YTickLabel', StimTypeCM(1:9))
xlabel('anatomical zones');
ylabel('Call Type')
title('Proportion of discriminant semantic units Selective (>1.6)')
subplot(1,3,3)
imagesc(PercInvDiscSSUCells_per7area(1:6,:)')
colorbar
set(gca,'XTickLabel', List_zones7(1:6),'XTick',1:6)
set(gca,'YTickLabel', StimTypeCM(1:9))
xlabel('anatomical zones');
ylabel('Call Type')
title('Proportion of discriminant semantic units Invariant (>0.6)')


figure()
GRAD2=[1 1 0; 1 0 0; 1 0.6 0; 0 0.5 1; 0 1 1; 0 1 0 ; 0.8 0.8 0.8];
for zz=1:max(ZONES7)
    Ind = intersect(SemAnatDisc, find(ZONES7==zz));
    plot(PCC_AV(Ind), LRI_max(Ind),'MarkerFaceColor', GRAD2(zz,:), 'LineStyle', 'none','Marker','o','MarkerSize',10,'MarkerEdgeColor',GRAD2(zz,:));
    axis([0 0.7 0 4])
    hold on
    pause()
end
legend(List_zones7)

ss=subplot(1,1,1);
gg=gscatter(PCC_AV(SemAnatDisc), LRI_max(SemAnatDisc), ZONES7(SemAnatDisc), 'yrmbcgk', '.......',[20 20 20 20 20 20 20]);
set( get(ss,'XLabel'), 'String', 'Discrimination index: Average PCC over categories' );
set( get(ss,'YLabel'), 'String', 'Selectivity Index: max LRI on diag PCC' );
title(ss, sprintf('%d SS cells', length(SemAnatDisc)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fig6d 6e Anatomy2 plot selective cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%per call cat%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SemAnatDiscSel=intersect(SemAnatDisc, find(LRI_max>1.6));
PlotSphere(LR(SemAnatDiscSel), RC(SemAnatDiscSel), DV(SemAnatDiscSel),PCC_max(SemAnatDiscSel),1,1,PCC_max_CT(SemAnatDiscSel),StimTypeCM(1:end-1),0,1);
figure()
bar(PercSelSSCells_per7area(:,2:end),'stacked')
legend(StimTypeCM(1:end-1), 'Location','NorthWest')
xlabel('anatomical zones');
ylabel('Proportion of selective units')
set(gca, 'XTickLabel',List_zones7);

SemAnatDiscSelSU=intersect(SemAnatDiscSU, find(LRI_max>1.6));
PlotSphere(LR(SemAnatDiscSelSU), RC(SemAnatDiscSelSU), DV(SemAnatDiscSelSU),PCC_max(SemAnatDiscSelSU),1,1,PCC_max_CT(SemAnatDiscSelSU),StimTypeCM(1:end-1),0,1);
figure()
bar(PercSelSSUCells_per7area(:,2:end),'stacked')
legend(StimTypeCM(1:end-1), 'Location','NorthWest')
xlabel('anatomical zones');
ylabel('Proportion of selective units')
set(gca, 'XTickLabel',List_zones7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fig6f 6g Anatomy2 plot slective cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%per call cat%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II_opt = (II_PCCmax-min(II_PCCmax(SemAnatDiscSel))+0.1)./(max(II_PCCmax(SemAnatDiscSel))-min(II_PCCmax(SemAnatDiscSel)));
figure()
hist(II_PCCmax(SemAnatDiscSel))
vline(median(II_PCCmax(SemAnatDiscSel)))
vline(median(II_PCCmax(SemAnatDisc)), 'k')
PlotSphere(LR(SemAnatDiscSel), RC(SemAnatDiscSel), DV(SemAnatDiscSel),II_opt(SemAnatDiscSel),1,1,PCC_max_CT(SemAnatDiscSel),StimTypeCM(1:end-1),0,1);
figure()
subplot(1,3,1)
plot(PCC_max(SemAnatDiscSel), II_PCCmax(SemAnatDiscSel),'LineStyle','none','Marker','o')
LinearModel.fit(PCC_max(SemAnatDiscSel),II_PCCmax(SemAnatDiscSel))
subplot(1,3,2)
plot(PCC_mean(SemAnatDiscSel), II_PCCmax(SemAnatDiscSel),'LineStyle','none','Marker','o')
LinearModel.fit(PCC_mean(SemAnatDiscSel),II_PCCmax(SemAnatDiscSel))
subplot(1,3,3)
plot(PCC_mean(SemAnatDiscSel), LRI_max(SemAnatDiscSel),'LineStyle','none','Marker','o')
LinearModel.fit(PCC_mean(SemAnatDiscSel),LRI_max(SemAnatDiscSel))

figure()
bar(PercInvSelSSCells_per7area(:,2:end),'stacked')
legend(StimTypeCM(1:end-1), 'Location','NorthWest')
xlabel('anatomical zones');
ylabel('Proportion of invariant units')
set(gca, 'XTickLabel',List_zones7);
ylim([0 1])

PlotSphere(LR(SemAnatDisc), RC(SemAnatDisc), DV(SemAnatDisc),II_max(SemAnatDisc),1,1,ZONES7(SemAnatDisc),List_zones7,0,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SemAnat = intersect(SemCellPV, xyzCells);
PlotSphere(LR(SemAnat), RC(SemAnat), DV(SemAnat),PCC_mean(SemAnat),1,1,ZONES(SemAnat)+1,['Unknown', ZONES_List],1);
PlotSphere(LR(SemAnat), RC(SemAnat), DV(SemAnat),2.^(LRI_max(SemAnat)),1,1,ZONES(SemAnat)+1,['Unknown', ZONES_List],1);

for subj=1:length(Subjects)
    subject=Subjects(subj)
    xyzCells=intersect(find(~isnan(DV_theo)), find(SUBJ==subject));
    if ~isempty(xyzCells)
        %PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_tot(xyzCells));
        PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),MeanSR(xyzCells),0);
        set(get(gca,'Title'), 'String', sprintf('Theoretical plot subject %s', SUBJECTS{subject+1}));
    end
end

% Now plot for each cell in each subject the PCC value of the call type
% that is best categorized with a color code for call type
for subj=1:length(Subjects)
    subject=Subjects(subj)
    xyzCells=intersect(find(~isnan(DV_theo)), find(SUBJ==subject));
    if ~isempty(xyzCells)
        %PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_tot(xyzCells));
        PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),PCC_max(xyzCells),1,1,PCC_max_CT(xyzCells),StimTypeCM(1:end-1),1);
        set(get(gca,'Title'), 'String', sprintf('max(PCC) Theoretical plot subject %s', SUBJECTS{subject+1}));
    end
end

% Now plot for each cell in each subject the LRI value of the call type
% that is best categorized with a color code for call type
for subj=1:length(Subjects)
    subject=Subjects(subj)
    xyzCells=intersect(find(~isnan(DV_theo)), find(SUBJ==subject));
    if ~isempty(xyzCells)
        %PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_tot(xyzCells));
        PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),LRI_maxPCCcat(xyzCells),1,1,PCC_max_CT(xyzCells),StimTypeCM(1:end-1),round(max(LRI_maxPCCcat)));
        set(get(gca,'Title'), 'String', sprintf('max(LRI) Theoretical plot subject %s', SUBJECTS{subject+1}));
    end
end

% Now plot for each cell in each subject the PCC value of the call type
% that is best categorized with a color code for anatomical region
for subj=1:length(Subjects)
    subject=Subjects(subj)
    xyzCells=intersect(find(~isnan(DV)), find(SUBJ==subject));
    if ~isempty(xyzCells)
        %PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_tot(xyzCells));
        PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),PCC_max(xyzCells),1,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],1);
        set(get(gca,'Title'), 'String', sprintf('max(PCC) Theoretical plot subject %s', SUBJECTS{subject+1}));
    end
end

% Now plot for each cell in each subject the LRI value of the call type
% that is best categorized with a color code for anatomical region
for subj=1:length(Subjects)
    subject=Subjects(subj)
    xyzCells=intersect(find(~isnan(DV_theo)), find(SUBJ==subject));
    if ~isempty(xyzCells)
        %PlotSphere(LR(xyzCells), RC(xyzCells), DV(xyzCells),MI_perArea.MI_tot(xyzCells));
        PlotSphere(LR_theo(xyzCells), RC_theo(xyzCells), DV_theo(xyzCells),LRI_maxPCCcat(xyzCells),1,1,ZONES(xyzCells)+1,['Unknown', ZONES_List],round(max(LRI_maxPCCcat)));
        set(get(gca,'Title'), 'String', sprintf('max(LRI) Theoretical plot subject %s', SUBJECTS{subject+1}));
    end
end

% Calculate the mean PCC_mean, mean maxLRI and mean invariance per electrode site
PCC_mean_Site = nan(1:NS);
LRI_max_Site = nan(1:NS);
II_m
for ss=1:NS
end
    

%% Plot the cells in MI space with spike rate on zaxes

figure(4)
GRAD2=cubehelix(ceil(max(MeanSR*1000)), 0.5, -1.1, 1.5, 0.5, [1,0]);
subplot(2,2,1)
for jj=1:NU
    plot(MI_perArea.MI_tot(jj), SemanticIndex.NonLinearBG(jj),'ko', 'MarkerFaceColor',GRAD2(ceil(MeanSR(jj)*1000),:));
    hold on
end
xlabel('Mutual Information of the Individual Vocalization matrix (bits)')
ylabel('Non-linear Semantic Index')
title(sprintf('All cells\nSpike Rate Hz'));
colorbar('YTickLabel', (50:50:250))
colormap(GRAD2)
axis([0 6 -0.2 0.5])
hold off

subplot(2,2,2)
for jj=1:NU
    if logpvRand(jj)>=2 && SemanticIndex.NonLinear(jj)>0 && logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.NonLinearBG(jj),'ko', 'MarkerFaceColor',GRAD2(ceil(MeanSR(jj)*1000),:));
        hold on
    end
end
xlabel('Mutual Information of the Individual Vocalization matrix (bits)')
ylabel('Non-linear Semantic Index')
title(sprintf('Semantic Cells\nSpike Rate Hz'));
colorbar('YTickLabel', (50:50:250))
colormap(GRAD2)
axis([0 6 -0.2 0.5])
hold off

subplot(2,2,3)
for jj=1:NU
    plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj),'ko', 'MarkerFaceColor',GRAD2(ceil(MeanSR(jj)*1000),:));
    hold on
end
xlabel('Mutual Information of the Individual Vocalization matrix (bits)')
ylabel('Semantic Index')
title(sprintf('All Cells\nSpike Rate Hz'));
colorbar('YTickLabel', (50:50:250))
colormap(GRAD2)
axis([0 6 0 1.4])
hold off

subplot(2,2,4)
for jj=1:NU
    if logpvRand(jj)>=2 && SemanticIndex.NonLinear(jj)>0 && logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj),'ko', 'MarkerFaceColor',GRAD2(ceil(MeanSR(jj)*1000),:));
        hold on
    end
end
xlabel('Mutual Information of the Individual Vocalization matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nSpike Rate Hz'));
colorbar('YTickLabel', (50:50:250))
colormap(GRAD2)
axis([0 6 0 1.4])
hold off

%% Plot the cells in MI space with two indices of selectivity
% Extract max value of LRI on PCC for each cell
MaxLRI = nan(NU,1);
for jj=1:NU
    MaxLRI(jj) = max(Selectivity.DiagLRI{jj});
    if isinf(MaxLRI(jj)) %get rid of the only infinite value. Note that this cell is not a semantic one anyways
        MaxLRI(jj) = 6;
    end
end

figure(5)
subplot(1,2,1)
GRAD3=cubehelix(ceil(max(Selectivity.DiagEntropyNorm*10000)), 0.5, -0.8, 2, 1, [1,0]);
for jj=1:NU
    if logpvRand(jj)>=2 && SemanticIndex.NonLinear(jj)>0 && logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.NonLinearBG(jj), 'ko', 'MarkerFaceColor',GRAD3(round(Selectivity.DiagEntropyNorm(jj)*10000),:));
    end
    hold on
end
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Non-linear Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Selectivity 1-Proportion entropy in categories\n'));
colorbar()
colormap(GRAD3)
axis([0 6 -0.2 0.5])
hold off

subplot(1,2,2)
for jj=1:NU
    if logpvRand(jj)>=2 && SemanticIndex.NonLinear(jj)>0 && logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj), 'ko', 'MarkerFaceColor',GRAD3(round(Selectivity.DiagEntropyNorm(jj)*10000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Selectivity 1-Proportion of entropy'));
colorbar()
colormap(GRAD3)
axis([0 6 0 1.4])


figure(6)
subplot(1,2,1)
GRAD4=cubehelix(ceil(max(MaxLRI*10000)), 0.5, -1, 2, 0.6, [1,0]);
for jj=1:NU
    if logpvRand(jj)>=2 && SemanticIndex.NonLinear(jj)>0 && logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.NonLinear(jj), 'ko', 'MarkerFaceColor',GRAD4(round(MaxLRI(jj)*10000),:));
    end
    hold on
end
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Selectivity LRI on PCC'));
colorbar()
colormap(GRAD4)
axis([0 6 -0.2 0.5])
hold off

subplot(1,2,2)
for jj=1:NU
    if logpvRand(jj)>=2 && SemanticIndex.NonLinear(jj)>0 && logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj), 'ko', 'MarkerFaceColor',GRAD4(round(MaxLRI(jj)*10000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Units\nz-axis: Selectivity LRI on PCC'));
colorbar()
colormap(GRAD4)
axis([0 6 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables for Neelabh Kumar (09/12/2022)

NK_Data.PCC = PCC_inv(:,1:9);
NK_Data.Sel = LRI(:,1:9);
NK_Data.Inv = II(:,1:9);
NK_Data.StimType = StimTypeCM(1:9);
NK_Data.SemUnitIndices = SemCellPV;
NK_Data.UnitName = cell(NU,1);
for ss=1:NU
    [~,UnitName,~]=fileparts(List_matfilepath{ss});
    NK_Data.UnitName{ss} = UnitName(9:end);
end

save('/Users/elie/Documents/CODE/data/matfile/DataElieTheunissen_EJN2015_4Neelabh.mat', '-struct', 'NK_Data')
