function test_gamma_one()
here=fileparts(mfilename('fullpath')); addpath(here);
addpath('C:\Users\Adam\Documents\GitHub\ZelanoLabScripts');
addpath('C:\Users\Adam\Documents\GitHub\Superlets\matlab-pure');
base='C:\Users\Adam\AppData\Local\Temp\claude\C--Users-Adam-Documents-GitHub-ZelanoLabScripts-plotsForBruce-2\2dd116fc-2bdb-4993-a330-8afe33295cdb\scratchpad\testfinals\';

tests = { 'cueTask',    [base '250623_Dupi_NMH_KS_1_cueTaskPreproc.mat'],       'macBP1';
          'audiobook',  [base '250623_Dupi_NMH_KS_1_breathingPreProc.mat'],     'macBP2';
          'focusedBreathing', [base '250623_Dupi_NMH_KS_1_breathingPreProc.mat'],'macBP2'};
for i=1:size(tests,1)
    tr=tests{i,1}; fp=tests{i,2}; ch=tests{i,3};
    fprintf('\n===== %s (%s) =====\n', tr, ch);
    t0=tic;
    out = extract_gamma_session(fp, tr, ch, struct());
    fprintf('  time=%.1fs  nBreaths=%d  perBreath=[%d x %d]\n', toc(t0), out.nBreaths, height(out.perBreath), width(out.perBreath));
    if height(out.perBreath)>0
        T=out.perBreath;
        fprintf('  vars(1:14): %s\n', strjoin(T.Properties.VariableNames(1:min(14,width(T))), ', '));
        showcols = {'peakZ','peakLatMs','peakFreq','anyBurst','burstLatMs','timeAboveMs','nBursts','dutyCycle','chirpSlope','freqSpan','apExp','gammaBumpDb','gammaPeakPresent','w2_rfreqPW','w2_rpowZ','p1_rfreqPW'};
        for c=1:numel(showcols)
            if any(strcmp(T.Properties.VariableNames,showcols{c}))
                v=T.(showcols{c});
                fprintf('   %-16s mean=%.3g  median=%.3g  nNaN=%d\n', showcols{c}, mean(v,'omitnan'), median(v,'omitnan'), sum(isnan(v)));
            end
        end
        fprintf('   goodBreath frac=%.2f ; mean nBreaths landmark returnCrossMs=%.0f inhalePeakMs=%.0f\n', ...
            mean(T.goodBreath,'omitnan'), mean(T.returnCrossMs,'omitnan'), mean(T.inhalePeakMs,'omitnan'));
    end
end
fprintf('\nDONE test_gamma_one\n');
end
