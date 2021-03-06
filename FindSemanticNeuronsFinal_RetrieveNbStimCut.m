
%% This script comes after FindSemanticNeuronsFinalPredict_PCCaboveChance_Invariance.m
resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')

%% Retrieve some output variables
% cell containing the path of each unit
load('/auto/k8/julie/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance.mat');
LM = length(List_matfilepath);
NbCell_estimate = LM;

% Vectors containing the nb of cut stimuli or intact stimuli of all units
Nb_CutStim = zeros(NbCell_estimate,1);
Nb_FullStim = zeros(NbCell_estimate,1);

%% Start to loop through files and extract information
cd /auto/k6/julie/matfile
input_dir=pwd;

for hh = 1:LM
    CM=load(List_matfilepath{hh});
    Nb_CutStim(hh)=CM.NbCutStim;
    Nb_FullStim(hh)=CM.NbFullStim;
    fprintf(1,'%d/%d\n',hh,LM)
end

%% Compile results into structures
sprintf('Compile results\n')
Conf.NbCutStim=Nb_CutStim;
Conf.Nb_FullStim=Nb_FullStim;

%% Save data
save(fullfile(resultsDirectory, 'SemanticAnalysis_LastONE.mat'),'Nb_CutStim','Nb_FullStim')
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis_LastONE.mat'))



