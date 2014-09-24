cd /auto/k6/julie/matfile
input_dir = pwd;
Subjects = dir(input_dir);
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        allFiles = dir(fullfile(input_dir, Indiv,'FirstVoc*.mat'));
        NF = length(allFiles);
        for nsite=1:8
            % find a matfile for each site
            Matsite='';
            for nf = 1:NF
                File_n=allFiles(nf).name;
                if str2num(File_n(14))==nsite
                    Matsite = allFiles(nf).name;
                    break
                end
            end
            if isempty(Matsite)
                fprintf(1, 'No file could be find for %s Site %d\n', Indiv, nsite);
            else
                % Retrieve the Matfile that contains the spectrograms
                MAT = load(fullfile(input_dir, Indiv,Matsite));
                
                
                %%%%YOUR CODE HERE
                
                
                % save the output matrix
                Res.CORR=CORR;
                Res.VocTypeSel=VocTypeSel;
                Res.TDTwav=TDTwav;
                Res.StimIndices= StimIndices;
                if ismac()
                    [status username] = system('who am i');
                    if strcmp(strtok(username), 'frederictheunissen')
                        if strncmp('/auto/fdata/solveig',stim_name, 19)
                        elseif strncmp('/auto/fdata/julie',stim_name, 17)
                            filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['FirstVoc_' Res.Site '.mat']);
                        end
                    elseif strcmp(strtok(username), 'elie')
                        filename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',Res.subject,['FirstVoc_' Res.Site '.mat']);
                    end
                else
                    filename=fullfile('/auto','k6','julie','matfile',Indiv,sprintf('%s_CorrFirstVoc_Site%d.mat',Indiv,nsite));
                end
                save(filename, '-struct', 'Res');
                fprintf('saved data under: %s\n', filename);
            end
        end
    end
end
