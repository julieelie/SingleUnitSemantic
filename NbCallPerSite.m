function NbCallPerSite
%% complete the statictics of NeuroVocalizationBank.mat created by hand that gives the number of vocalizations for each call type and for each individual in the whole neuro bank
%here we add the average number of vocalizations per recording site per
%category

cd /auto/k6/julie/matfile
input_dir = pwd;
Subjects = dir(input_dir);
StimTypeCM={'Ag','Be','DC','Di','Ne','LT','Te','Th','song'};
NStimperSite=nan(25,9);
NbTrialsperSite=nan(25,1);
Site=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        allFiles = dir(fullfile(input_dir, Indiv,'ConfMat*.mat'));
        NF = length(allFiles);
        for nsite=1:8
            % find a matfile for each site
            Matsite='';
            for nf = 1:NF
                File_n=allFiles(nf).name;
                if str2num(File_n(13))==nsite
                    Matsite = allFiles(nf).name;
                    break
                end
            end
            if isempty(Matsite)
                fprintf(1, 'No file could be find for %s Site %d\n', Indiv, nsite);
            else
                Site=Site+1;
                % Retrieve the number of stims in each category
                MAT = load(fullfile(input_dir, Indiv,Matsite));
                InterestVocType=zeros(size(MAT.VocTypeSel));
                for cc=1:length(StimTypeCM)
                    ct = StimTypeCM(cc);
                    NStimperSite(Site,cc)=sum(strcmp(ct,MAT.VocTypeSel));
                    InterestVocType = InterestVocType + strcmp(ct,MAT.VocTypeSel);
                end
                % Retrieave the number of trials for each stim
                NbTrialsperstim=sum(MAT.confusionMatrix{1} .*MAT.neventMatrix(1),2);
                NbTrialsperSite(Site)=mean(NbTrialsperstim(find(InterestVocType)));
                
                
            end
        end
    end
end
save('/auto/k6/julie/matfile/NeuroVocalizationBankNstimperSite.mat', 'NStimperSite','StimTypeCM','NbTrialsperSite');


        
                    