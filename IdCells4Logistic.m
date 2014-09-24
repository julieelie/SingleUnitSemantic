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