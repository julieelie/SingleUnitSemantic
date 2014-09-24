%% Find cells that are highly discriminative (PCC>0.5) and highly selective (LRI>3)
% This file goes along with Data4Logistic.mat
% PCC_units is the value of PCC for each unit (line) and each call category
% (column)
% LRI is the index of selectivity for each unit (line) and each call category
% (column)
% Note that PCC_units and LRI containes the values for all 1401 single
% units
% SemCellPV are the indices of significantly semantic cells
% Signif_PCC is a cell array that contains for each call category the indices of
% units that have a significant values of PCC compare to chance level
% List_matfilepath contains the exact name of each unit
% StimTypeCM contains the names of the call categories (legend of the column of LRI, PCC_units)
% MeanSpikeRate is a structure that contains the Values of spike rate for
% each stim and each unit, and the TDT names of each stim for each unit.

load('/Users/frederictheunissen/Documents/Data/Julie/Semantic Project/Data4Logistic.mat');

Cell4LogisticName = cell(9,1);
Cell4LogisticId = cell(9,1);
Ncat = size(PCC_units,2);
for cc=1:(Ncat)
    Cell4Logistic_local = intersect(find(PCC_units(:,cc)>0.5), find(LRI(:,cc)>3));
    if ~isempty(Cell4Logistic_local)
        fprintf('There are %d interesting cells for %s\n', length(Cell4Logistic_local), StimTypeCM{cc});
        %check that these cells are in SemCellPV and Signif_PCC
        Cell4Logistic_local2 = intersect(intersect(SemCellPV, Signif_PCC{cc}),Cell4Logistic_local);
        if ~isempty(Cell4Logistic_local2);
            fprintf('%d out of %d are significantly discriminant semantic units\n', length(Cell4Logistic_local2), length(Cell4Logistic_local));
            Cell4LogisticName{cc}=List_matfilepath(Cell4Logistic_local2);
            Cell4LogisticId{cc} = Cell4Logistic_local2;
        else
               fprintf('%d out of %d are significantly discriminant semantic units\n', length(Cell4Logistic_local2), length(Cell4Logistic_local));
        end
    else
        fprintf('There are NO interesting cells for %s\n', StimTypeCM{cc});
    end
end

%% Now that we have list of interesting cells - let's do some plots

% Start with one for a Distance Call
cellId = Cell4LogisticId{3}(1);
matFileName = Cell4LogisticName{3}{1};
spikeRates = MeanSpikeRate.Values{cellId};

% Load acoustical analysis file
load('/Users/frederictheunissen/Documents/Data/Julie/Acoustical Analysis/matfile/BlaBro09xxF/DFAvocCuts_Site3_L2500R2300.mat');

% Some checking
% Are we dealing with the same vocalization type?
if ( ~strcmp(vocTypes{3}, StimTypeCM{3}) )
    fprintf(1, 'Error: mismatch in voc type between Acoustical Analysis and Neural Data\n');
end


% Do we have the same stim in the two mat files? No we certainly don't
% exaclty nStimsRate should always be bigger than nStimsCut.
nStimsCut = length(stimNameCuts);
nStimsRate = length(spikeRates);
indFix = zeros(1, nStimsRate);   % because this is smaller....

for is=1:nStimsCut
    for is2=1:nStimsRate      
        if (strcmp(MeanSpikeRate.TDT_StimNames{cellId}{is2}, stimNameCuts{is}) )
            indFix(is2) = is;
        end
    end
end

countMiss = 0;
for is=1:nStimsCut
    if (isempty(find(is == indFix)))
        countMiss = countMiss + 1;
        fprintf(1,'(%d) Missing %s stim in Mean Spike Rate:\n\t%s\n', countMiss, vocTypeCuts{is}, stimNameCuts{is});
    end
end

countMiss = 0;
indFixRate = zeros(1,nStimsRate);
for is2=1:nStimsRate
    if (indFix(is2) == 0 )
      countMiss = countMiss + 1;
      fprintf(1,'(%d) Missing stim in Stim Name:\n\t%s\n', countMiss, MeanSpikeRate.TDT_StimNames{cellId}{is2}); 
    else
        indFixRate(is2) = is2;
    end
end


% Plotting...
indFix2 = indFix;
indFix2(indFix2 == 0) = [];
indFixRate2 = indFixRate;
indFixRate2(indFixRate2 == 0) = [];

figure(1);
for ip=1:length(indFix2)
    if (strcmp(vocTypeCuts{indFix2(ip)}, 'DC'))
        h1= plot(PCCoeffLR(indFix2(ip),3), spikeRates(indFixRate2(ip)), 'r+');
    else
        h2 = plot(PCCoeffLR(indFix2(ip),3), spikeRates(indFixRate2(ip)), 'k+');
    end
    
    hold on;
end
hold off;
xlabel('Sound Coefficient on Logistic Function');
ylabel('Mean Rate');
legend([h1 h2], 'DC', 'Other');


