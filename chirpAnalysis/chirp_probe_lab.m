% chirp_probe_lab -- one-time lab environment + data-contract probe for the chirp analysis.
%   Run on the lab desktop: matlab -batch "run('E:\chirp_tmp\chirp_probe_lab.m')"
%   Prints toolbox licenses (decides TFR/ridge/chirplet engine), labPaths roots, and the
%   data contract of one cue final (labels, macBP set, behDat, finalOnset, fs).
%   Loads a STALE final only to inspect STRUCTURE (unaffected by the bestMac/noise rerun).

fprintf('=== MATLAB %s ===\n', version);
tb = {'Wavelet_Toolbox','Signal_Toolbox','Statistics_Toolbox','Curve_Fitting_Toolbox', ...
      'Optimization_Toolbox','Distrib_Computing_Toolbox'};
for i = 1:numel(tb), fprintf('  license %-24s = %d\n', tb{i}, license('test', tb{i})); end
fprintf('  exist tfridge=%d wsst=%d fsst=%d pspectrum=%d emd=%d\n', ...
    exist('tfridge'), exist('wsst'), exist('fsst'), exist('pspectrum'), exist('emd'));

repo = 'G:\My Drive\GitHub\ZelanoLabScripts';
addpath(repo); addpath(fullfile(repo,'cueAnalysis'));
fprintf('\n=== labPaths ===\n');
try
    L = labPaths();
    fns = fieldnames(L);
    for i = 1:numel(fns)
        v = L.(fns{i});
        if ischar(v) || isstring(v), fprintf('  L.%-14s = %s\n', fns{i}, char(v)); end
    end
catch ME, fprintf('  labPaths FAILED: %s\n', ME.message); end

% inspect one final's structure (STALE file is fine for structure)
cands = { ...
  'E:\Lab_Common\Dupi\250623_Dupi_NMH_KS_1\preProc\250623_Dupi_NMH_KS_1_cueTaskPreproc.mat', ...
  'E:\Lab_Common\Dupi\250623_Dupi_NMH_KS_2\preProc\250623_Dupi_NMH_KS_2_cueTaskPreproc.mat'};
fp = ''; for i = 1:numel(cands), if isfile(cands{i}), fp = cands{i}; break; end, end
fprintf('\n=== probe final: %s ===\n', fp);
if ~isempty(fp)
    try
        s = load(fp); fn = fieldnames(s); od = s.(fn{1}); clear s;
        fprintf('  top var = %s\n', fn{1});
        labs  = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
        isMac = cellfun(@(x) contains(x,'macBP'), labs);
        isRsp = cellfun(@(x) contains(x,'rsp'),   labs);
        fprintf('  nChan=%d  fs=%g  task=%s  type=%s\n', numel(labs), od.fs, ...
            char(string(od.task)), char(string(getfieldsafe(od,'type'))));
        fprintf('  nMacBP=%d  labels: %s\n', sum(isMac), strjoin(labs(isMac), ', '));
        fprintf('  rsp idx=%d flip=%d (rsp chans: %s)\n', getfieldsafe(od,'rspIDX'), ...
            getfieldsafe(od,'rspFlip'), strjoin(labs(isRsp), ', '));
        if isfield(od,'bestMac'), fprintf('  bestMac = %s\n', char(string(od.bestMac)));
        else, fprintf('  bestMac: ABSENT (stale/pre-rerun)\n'); end
        fprintf('  data size = [%d x %d]  (%.1f s)\n', size(od.data,1), size(od.data,2), size(od.data,2)/od.fs);
        if isfield(od,'behDat') && istable(od.behDat)
            bd = od.behDat;
            fprintf('  behDat [%d x %d] vars: %s\n', height(bd), width(bd), ...
                strjoin(bd.Properties.VariableNames, ','));
            if ismember('finalOnset', bd.Properties.VariableNames)
                fo = double(bd.finalOnset);
                fprintf('  finalOnset: n=%d  range=[%d %d] samp  (%.0f..%.0f s)\n', ...
                    numel(fo), min(fo), max(fo), min(fo)/od.fs, max(fo)/od.fs);
                d = diff(sort(fo(isfinite(fo))));
                fprintf('  median inter-onset gap = %.0f ms (min %.0f, p10 %.0f)\n', ...
                    median(d)/od.fs*1000, min(d)/od.fs*1000, prctile(d,10)/od.fs*1000);
            end
            if ismember('sniffLabel', bd.Properties.VariableNames)
                sl = string(bd.sniffLabel); u = unique(sl);
                fprintf('  sniffLabel values: %s\n', strjoin(cellstr(u)', ','));
            end
        end
        if isfield(od,'TTL'), fprintf('  TTL class=%s\n', class(od.TTL)); end
    catch ME, fprintf('  LOAD/INSPECT FAILED: %s\n', ME.message); end
end
fprintf('PROBE_DONE\n');

function v = getfieldsafe(s, f)
    if isfield(s, f), v = s.(f); else, v = NaN; end
end
