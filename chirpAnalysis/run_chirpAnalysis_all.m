function run_chirpAnalysis_all(stage, sessFilter, opts)
% RUN_CHIRPANALYSIS_ALL  Master wrapper for the chirp (one-vs-two oscillator) analysis.
% =========================================================================
%   Stages, in dependency order (each also runnable standalone):
%     (0) 'phase0'   -- chirp_phase0_harness: synthetic sweep that GATES interpretation.
%                       Heavy; run on the lab desktop. Writes <outDir>/phase0/.
%     (1) 'pipeline' -- run_chirp_analysis_pipeline: per-session battery over FRESH E: cue
%                       finals (D5). Writes per-trial CSVs + <id>_chirp.mat + figures.
%     (2) 'group'    -- chirp_group_adjudicate: per-session verdict vectors + group tests.
%     (3) 'report'   -- knit chirpAnalysisReport.Rmd -> HTML.
%
%   USAGE
%     run_chirpAnalysis_all('all')                 % phase0 + pipeline + group + report
%     run_chirpAnalysis_all('pipeline', {'<id>'})  % one session
%     run_chirpAnalysis_all('phase0', [], struct('quick',true))   % small validation sweep
%
%   Outputs go to CHIRP_OUTDIR (env) or <labPaths.figPath>/groupStatFigs/chirpAnalysis.
% =========================================================================
    if nargin<1||isempty(stage), stage='all'; end
    if nargin<2, sessFilter=[]; end
    if nargin<3, opts=struct(); end
    setup_chirpAnalysis_paths();
    C = chirp_config();
    outDir = C.outDir; if isempty(outDir), try, L=labPaths(); outDir=fullfile(L.figPath,'groupStatFigs','chirpAnalysis'); catch, outDir=fullfile(pwd,'chirpOut'); end, end
    if ~isfolder(outDir), mkdir(outDir); end
    do = @(s) any(strcmp(stage,{s,'all'}));

    if do('phase0')
        fprintf('\n########## (0) Phase-0 harness ##########\n');
        sw = struct('outDir', fullfile(outDir,'phase0'));
        if isfield(opts,'tiny') && opts.tiny
            sw.overlapGrid=500; sw.DfGrid=7; sw.ampGrid={[1 1]}; sw.nTrials=5;
            sw.families={'singleLin','dualIndep'}; sw.windowSources={'pipeline'};
        elseif isfield(opts,'quick') && opts.quick
            sw.overlapGrid=[200 500 1000 2000]; sw.DfGrid=[5 10]; sw.ampGrid={[1 1]}; sw.nTrials=20;
            sw.families={'singleLin','singleNL','dualIndep','dualLocked'};
            C.ridge.nSurr=40; C.ridge.nNullTrials=2; C.surr.n=100; C.phase.nPerm=500;  % lighter for the sweep
        end
        chirp_phase0_harness(sw, C);
    end
    if do('pipeline')
        fprintf('\n########## (1) per-session pipeline ##########\n');
        run_chirp_analysis_pipeline(sessFilter, opts);
    end
    if do('group')
        fprintf('\n########## (2) group adjudication ##########\n');
        try, chirp_group_adjudicate(outDir, C); catch ME, fprintf('group skipped: %s\n', ME.message); end
    end
    if do('report')
        fprintf('\n########## (3) knit report ##########\n');
        knitReport(C, outDir);
    end
    fprintf('\nrun_chirpAnalysis_all(%s) done.\n', stage);
end

function knitReport(C, outDir)
    rmd = fullfile(fileparts(mfilename('fullpath')), 'chirpAnalysisReport.Rmd');
    cands = {'C:\Program Files\R\R-4.4.1\bin\Rscript.exe','C:\Program Files\R\R-4.4.0\bin\Rscript.exe','Rscript'};
    rs=''; for k=1:numel(cands), if strcmp(cands{k},'Rscript')||isfile(cands{k}), rs=cands{k}; break; end, end
    if isempty(rs), fprintf('  Rscript not found; knit by hand:\n  rmarkdown::render(''%s'', params=list(outDir=''%s''))\n', strrep(rmd,'\','/'), strrep(outDir,'\','/')); return; end
    expr = sprintf('rmarkdown::render(''%s'', output_file=''chirpAnalysisReport.html'', output_dir=''%s'', params=list(outDir=''%s''), quiet=TRUE)', ...
        strrep(rmd,'\','/'), strrep(outDir,'\','/'), strrep(outDir,'\','/'));
    cmd = ['"' rs '" -e "' expr '"'];
    st = system(cmd);
    if st==0, fprintf('  report -> %s\n', fullfile(outDir,'chirpAnalysisReport.html')); else, fprintf('  knit returned %d; run by hand.\n', st); end
end
