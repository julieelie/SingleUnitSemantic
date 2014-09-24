function [calfilename]=nb_cutStim_cal(MatfilePath)
%% This code run through FirstVoc files and calculate the number of vocalization/pst that were truncated or kept intact. the output values are added to the confusion matrix file.
Res=load(MatfilePath, 'VocType', 'Trials', 'Trials_BG','TDT_wavfiles', 'SectionLength', 'Site', 'subject','Section_cat');
%% Restrict analysis to stim that are not whine or mlnoise
RestrictAna=intersect(find(~strcmp(Res.VocType, 'Wh')),find(~strcmp(Res.VocType, 'mlnoise')));% localize the stims we don't want to study here and get rid of them
%% Retrieve the number of cut and non-cut stims (psths)
NbCutStim=sum(strcmp(Res.Section_cat(RestrictAna),'cut'));
NbFullStim=sum(strcmp(Res.Section_cat(RestrictAna),'full'));
%% Store Values
if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['ConfMat_' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',['ConfMat_' Res.Site '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',Res.subject,['ConfMat_' Res.Site '.mat']);
end

save(calfilename, 'NbCutStim','NbFullStim','-append');
fprintf(1,'done making calculus on %s\nData save under %s\n',MatfilePath, calfilename);
clear Res 
end 
