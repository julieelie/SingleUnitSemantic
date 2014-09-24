resultsDirectory='/auto/k8/julie';
%Create some output variables
NS=30;%estimation of the number of sites
Opt_PC=nan(NS,2);
PCC_DFA = nan(NS,9);
Voc_DFA = cell(NS,1);
Site_names = cell(NS,1);
%Run through individuals and DFA files
cd /auto/fdata/fet/julie/Acoustical' Analysis'/
input_dir=pwd;
Subjects = dir(input_dir);
ii=0;

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        DFAfiles = dir(fullfile(Idir,'DFA*.mat'));
        Ndfa = length(DFAfiles);
        for dd=1:Ndfa
            DFA = load(fullfile(Idir,DFAfiles(dd).name));
            fprintf('Loaded %s\n',DFAfiles(dd).name)
            ii=ii+1;
            Site_names{ii}=DFAfiles(dd).name;

            % find the optimum number of PCs for each site (max value of DFA)
            Nsteps = length(DFA.PCC_info);
            PCC_Local=nan(Nsteps,2);
            for cc=1:Nsteps
                PCC_Local(cc,1)=DFA.PCC_info(cc).PCC_Total_DFA;
                PCC_Local(cc,2)=DFA.PCC_info(cc).PCC_M_DFA;
            end
            MaxPCC = max(PCC_Local);    
            Opt_PC(ii,1)=DFA.PCC_info(find(PCC_Local(:,1)==MaxPCC(1))).nb;
            Opt_PC(ii,2)=DFA.PCC_info(find(PCC_Local(:,2)==MaxPCC(2))).nb;
            if Opt_PC(ii,1)~=Opt_PC(ii,2)
                fprintf('Optimal nb of PC is different between Total and M calculation of PCC\nTotal=%d\nM=%d\n',Opt_PC(ii,1),Opt_PC(ii,2));
            end
            
            %Extract PCC of each category each site
            PCC_DFA(ii,:)=diag(DFA.PCC_info(find(PCC_Local(:,2)==MaxPCC(2))).PCC_group_DFA);
            Voc_DFA{ii} = DFA.vocTypes;
        end
    end
end
DFA.PCC_cat=PCC_DFA(1:ii,:);
DFA.VocTypes = Voc_DFA(1:ii);
DFA.Opt_PCs=Opt_PC(1:ii,:);
DFA.Site_names=Site_names(1:ii);
save(fullfile(resultsDirectory,'DFAAcoustic.mat'), '-struct', 'DFA');

%%%Few figures
DFA=load('/Users/elie/Documents/MATLAB/data/matfile/DFAAcoustic.mat');
figure()
subplot(1,2,1)
hist(DFA.Opt_PCs(:,1))
subplot(1,2,2)
hist(DFA.Opt_PCs(:,2))


