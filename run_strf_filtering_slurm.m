jackKnives = [5 8 9 10];
directory = '07032013T1342_ieee_sentences_noisexbabble_0.0SNR';
threshold = 0.5;

cmd = 'ap; strf_filtering_optimize_wavsqerror(''%s'', %d, %2.1f);';
resultsDirectory = '/auto/k8/tlee/noise_filtering/results';

for ii = 1: length(jackKnives);
    jobParams = struct;
    jobParams.partition = 'all';
    jobParams.cpus = 4;
    jobParams.memory = 7000;
    jobParams.out = fullfile(resultsDirectory, directory, sprintf('slurm_opt_%d.txt', jackKnives(ii)));
    jobParams.err = jobParams.out;
    icmd = sprintf(cmd, directory, jackKnives(ii), threshold); 
    fprintf('%d. Calling slurm_sbatch with command %s\n',ii, icmd);
    slurm_sbatch(icmd,jobParams);
end
