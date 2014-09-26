function plot_Subject_AvgSR_perunit(Indiv)

cd /Users/elie/Documents/MATLAB/data/calmatfile
input_dir=pwd;
Idir = fullfile(input_dir, Indiv);
calfiles=dir(fullfile(Idir,'Cal*_nc.mat'));
LC=length(calfiles);
for cc=1:LC
    calfile=fullfile(Idir, calfiles(cc).name);
    fprintf('Loading Matfile of %s\n', calfiles(cc).name);
    Cal=load(calfile);
    figure(7)
    H=bar(Cal.AvgRcallcat);
    hold on
    h2=errorbar(Cal.AvgRcallcat, Cal.StdRcallcat, 'b.');
    set(gca(), 'Xtick', 1:length(Cal.StimType));
    set(gca(), 'XTickLabel', Cal.StimType);
    pause
end