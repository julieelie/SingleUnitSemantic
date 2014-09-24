function FiguresConfusionMatrices(MatfilePath)
%List_matfilepath, MaxUnit
%MatfilePath=List_matfilepath{MaxUnit};
MAT=load(MatfilePath);
Iwin=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));% we are choosing the window size that gives the best value of MI conf in the CT matrix


% add some empty rows and lines between vocalization categories
MATRIX=MAT.confusionMatrix{Iwin};
UVT = unique(MAT.VocTypeSel);
BG = find(strcmp(UVT, 'BG'));
UVT = [UVT(1:BG-1) ; UVT(BG+1:end) ; UVT(BG)];%keep background category at the end as for the original matrix
NUVT=length(UVT);
MATRIX2 = ones(size(MATRIX,1)+NUVT-1,size(MATRIX,2)+NUVT-1)*0.02;
%first find min and max indices of each vocalization types
MIN=nan(NUVT,1);
MAX=nan(NUVT,1);
for tt = 1:NUVT
    vt=UVT{tt};
    MAX(tt)=max(find(strcmp(MAT.VocTypeSel, vt)));
    MIN(tt)=min(find(strcmp(MAT.VocTypeSel, vt)));
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

%Clim=[0 (max(max(MATRIX)))*2];
Clim=[0 0.02];
createconfusionmatrix2(MATRIX2,MIN,MAX,UVT, Clim) % sound file matrix with white barres

createconfusionmatrix(MATRIX, Clim, MIN,MAX,UVT) % sound file matrix without the white barres

% recalculating the corresponding Call type matrix with the call type in
% the same order
UVT=unique(MAT.VocTypeSel);
BG = find(strcmp(UVT, 'BG'));
UVT = [UVT(1:BG-1) ; UVT(BG+1:end) ; UVT(BG)];%keep background category at the end as for the original matrix
NUVT=length(UVT);
confusion_matrix_CallType = zeros(NUVT, NUVT);
for vtR=1:NUVT
    stR=UVT(vtR);
    selectorR=strcmp(MAT.VocTypeSel, stR);
    selectorIndR=find(selectorR);
    for vtC = 1:NUVT
        stC=UVT(vtC);
        selectorC=strcmp(MAT.VocTypeSel, stC);
        selectorIndC=find(selectorC);
        confusion_matrix_CallType(vtR,vtC)=sum(sum(MATRIX(selectorIndR, selectorIndC)));
    end
end
ProbaCat = sum(confusion_matrix_CallType, 2);
repProbaCat = repmat(ProbaCat, 1, size(confusion_matrix_CallType,2));
confusion_matrix_CallType_cond = confusion_matrix_CallType ./ repProbaCat;

Clim=[0 1];
createconfusionmatrixCT(confusion_matrix_CallType_cond, UVT, Clim);

% recalculating a random Call type matrix with the call type in
% the same order
UVT=unique(MAT.VocTypeSel);
BG = find(strcmp(UVT, 'BG'));
UVT = [UVT(1:BG-1) ; UVT(BG+1:end) ; UVT(BG)];%keep background category at the end as for the original matrix
NUVT=length(UVT);
confusion_matrix_CallType = zeros(NUVT, NUVT);
rng('shuffle'); % seeds the random number generator based on the current time
VocTypeSel_rand=MAT.VocTypeSel(randperm(length(MAT.VocTypeSel)));
for vtR=1:NUVT
    stR=UVT(vtR);
    selectorR=strcmp(VocTypeSel_rand, stR);
    selectorIndR=find(selectorR);
    for vtC = 1:NUVT
        stC=UVT(vtC);
        selectorC=strcmp(VocTypeSel_rand, stC);
        selectorIndC=find(selectorC);
        confusion_matrix_CallType(vtR,vtC)=sum(sum(MATRIX(selectorIndR, selectorIndC)));
    end
end
ProbaCat = sum(confusion_matrix_CallType, 2);
repProbaCat = repmat(ProbaCat, 1, size(confusion_matrix_CallType,2));
confusion_matrix_CallType_cond = confusion_matrix_CallType ./ repProbaCat;
Clim=[0 1];
UVTrand={'Cat1' 'Cat2' 'Cat3' 'Cat4' 'Cat5' 'Cat6' 'Cat7' 'Cat8' 'Cat9' 'Cat10'};
createconfusionmatrixCT(confusion_matrix_CallType_cond, UVTrand, Clim);

end