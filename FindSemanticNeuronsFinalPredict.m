%% This script is following FindSemanticNeuronsFinal and calculate...
    ...the Semantic Information for random matrices where only the columns...
        ... have been shuffled to estimate the significance of SI for each unit
resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')

FigFlag =0;


%% Retrieve some output variables
% cell containing the path of each unit
load('/auto/k8/julie/SemanticAnalysis.mat');
LM = length(List_matfilepath);
NbCell_estimate = LM;

% Vectors containing the mutual information per area for Random matrices of all units

AvMI_tot_RandP = zeros(NbCell_estimate,2); % first column is for rand matrices, second column for rand matrices with Background sections not shuffled
SDMI_tot_RandP = zeros(NbCell_estimate,2);
AvMI_diag_uni_RandP=zeros(NbCell_estimate,2);
SDMI_diag_uni_RandP=zeros(NbCell_estimate,2);
AVMI_diag_uni_cat_RandP=zeros(NbCell_estimate,2);
SDMI_diag_uni_cat_RandP=zeros(NbCell_estimate,2);
MI_uni_diag_cat_RandP = zeros(NbCell_estimate, 1000);
MI_uni_diag_cat_RandBGP = zeros(NbCell_estimate, 1000);


%% Start to loop through files and extract information
cd /auto/k6/julie/matfile
input_dir=pwd;
ii = 0;

for hh = 1:LM
    [Idir, Matfile, ext] = fileparts(List_matfilepath{hh});
    
    %retrieve Random confusion matrices files
    Randmatfiles=dir(fullfile(Idir, 'RandPMat*.mat'));
    rm=length(Randmatfiles);
    SS_Ind=zeros(rm,1);
    for ff=1:rm
        if ~isempty(strfind(Randmatfiles(ff).name, 'ss'))
            SS_Ind(ff)=1;
        end
    end
    IndicesRM=find(SS_Ind);
    RM=length(IndicesRM);
     
    ii=ii+1;
            
        

 %% Upload the random matrices and calculate MI per area
    fprintf('%d/%d Upload Random matrices and caluclate MI per area\n', hh, LM)
    for rr=1:RM
        RandMatfile = Randmatfiles(IndicesRM(rr)).name;
        if strcmp(Matfile(9:end), RandMatfile(10:end-4))
            RMAT = load(fullfile(Idir, RandMatfile));
            break
        end
    end
    rt = 2;%I set it to 2 on 14/07/29 because we are not using values of random matrices without BG fixed
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
            if rem(rm,100)==0
                rm
            end
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
        AvMI_tot_RandP(ii,rt) = mean(MI_tot_rand);
        SDMI_tot_RandP(ii,rt) = std(MI_tot_rand);
        fprintf('the MI_confusion of the IV Matrix is %f for that cell\nThe average MI of random matrices is %f+/-%f\n', MI_tot(ii), AvMI_tot_RandP(ii,rt),SDMI_tot_RandP(ii,rt))
    %pause
        AvMI_diag_uni_RandP(ii,rt)=mean(MI_diag_uni_rand);
        SDMI_diag_uni_RandP(ii,rt)=std(MI_diag_uni_rand);
        AVMI_diag_uni_cat_RandP(ii,rt)=mean(MI_diag_uni_cat_rand);
        SDMI_diag_uni_cat_RandP(ii,rt)=std(MI_diag_uni_cat_rand);
        if rt==1
            MI_uni_diag_cat_RandP(ii,:) = MI_diag_uni_cat_rand;
        elseif rt==2
            MI_uni_diag_cat_RandBGP(ii,:) = MI_diag_uni_cat_rand;
        end
        rt = rt +1;
    end
            
            

end

%% Compile results into structures
sprintf('Compile results\n')
MI_perArea.AvMI_tot_RandP = AvMI_tot_RandP(1:ii,:); % first column is for rand matrices, second column for rand matrices with Background sections not shuffled
MI_perArea.SDMI_tot_RandP = SDMI_tot_RandP(1:ii, :);
MI_perArea.AvMI_diag_uni_RandP=AvMI_diag_uni_RandP(1:ii, :);
MI_perArea.SDMI_diag_uni_RandP=SDMI_diag_uni_RandP(1:ii,:);
MI_perArea.AVMI_diag_uni_cat_RandP=AVMI_diag_uni_cat_RandP(1:ii,:);
MI_perArea.SDMI_diag_uni_cat_RandP=SDMI_diag_uni_cat_RandP(1:ii,:);
MI_perArea.MI_uni_diag_cat_RandP = MI_uni_diag_cat_RandP(1:ii,:);
MI_perArea.MI_uni_diag_cat_RandBGP = MI_uni_diag_cat_RandBGP(1:ii,:);


SemanticIndex.Observed = MI_perArea.MI_diag_uni_cat ./ MI_perArea.MI_tot;
SemanticIndex.RandP = MI_perArea.AVMI_diag_uni_cat_RandP(:,1) ./ MI_perArea.AvMI_tot_RandP(:,1);
SemanticIndex.RandBGP = MI_perArea.AVMI_diag_uni_cat_RandP(:,2) ./ MI_perArea.AvMI_tot_RandP(:,2);
SemanticIndex.PvalueP = normpdf((MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_RandP(:,1)) ./ MI_perArea.SDMI_diag_uni_cat_RandP(:,1));
SemanticIndex.PvalueBGP = normpdf((MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_RandP(:,2)) ./ MI_perArea.SDMI_diag_uni_cat_RandP(:,2));
SemanticIndex.NonLinearP = (MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_RandP(:,1)) ./ MI_perArea.MI_tot;
SemanticIndex.NonLinearBGP = (MI_perArea.MI_diag_uni_cat - MI_perArea.AVMI_diag_uni_cat_RandP(:,2)) ./ MI_perArea.MI_tot;

%% Save data
save(fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict.mat'))
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis_ShuffPredict.mat'))



