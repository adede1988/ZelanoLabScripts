function o15_calibrate_noise()
% O15_CALIBRATE_NOISE  Calibrate the RELATIVE sharp-deflection threshold K for the
%   O15 dataset (the cue analysis's relative rule, ported -- mirrors
%   cue_calibrate_noise but the rejection unit here is the SNIFF, not the trial).
%
%   For each fresh session's bestMac: d = 10 ms max-min over the whole recording;
%   robust-z it: zd = (d - median(d)) / (1.4826*MAD(d)). Per SNIFF (epoch
%   [-1.75,5.75]s) take max(zd). Sweep K and report the dataset-wide fraction of
%   sniffs that would be rejected (max zd > K), the worst single-session fraction,
%   how many sessions would be fully excluded, and a per-sniff-type breakdown
%   (start/free/confirm). Pick K for >0%, <=~20% overall loss with no session at
%   100%, then set it in o15_noise_K.m and do the rerun.
%
%   Writes the raw per-sniff max-z to groupStatFigs\_calib_sniffMaxZ_O15.csv.
%   Read-only w.r.t. the finals; does NOT change any results by itself.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'cueAnalysis')); addpath(fullfile(repo,'O15Analysis'));
    groupDir = getenv('O15_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(labPaths().figPath, 'groupStatFigs'); end

    epWin = [-1.75 5.75]; win = 5;   % 10 ms window @ 500 Hz (matches cue_noise_trials)
    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);

    rows = {};   % sessID, group, sniffType, sniff#, sniffMaxZ
    for i = 1:height(T)
        id = char(T.sessID(i));
        try
            Sv = load(char(T.path(i))); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
        catch, continue; end
        if ~isfield(od,'bestMac') || isempty(od.bestMac), clear od; continue; end
        labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
        ci = find(strcmp(od.bestMac, labs), 1); if isempty(ci), clear od; continue; end
        x = double(od.data(ci,:)); fs = od.fs;
        d = movmax(x,win) - movmin(x,win);
        med = median(d,'omitnan'); sg = 1.4826*median(abs(d-med),'omitnan'); if sg<=0, sg = eps; end
        zd = (d - med)/sg;
        s0 = round(epWin(1)*fs); s1 = round(epWin(2)*fs); Tt = numel(x);
        bd = od.behDat; sl = string(bd.sniffLabel); fo = round(bd.finalOnset);
        for t = 1:numel(fo)
            a = fo(t)+s0; b = fo(t)+s1; if a<1 || b>Tt, continue; end
            rows(end+1,:) = {id, char(T.group(i)), char(sl(t)), t, max(zd(a:b))}; %#ok<AGROW>
        end
        fprintf('%-26s med(d)=%.1f MAD-sd=%.1f  sniffs=%d\n', id, med, sg, numel(fo));
        clear od x d zd;
    end

    if isempty(rows), fprintf('No sessions with bestMac found.\n'); return; end
    C = cell2table(rows, 'VariableNames', {'sessID','group','sniffType','sniff','sniffMaxZ'});
    writetable(C, fullfile(groupDir, '_calib_sniffMaxZ_O15.csv'));

    z = C.sniffMaxZ; sess = string(C.sessID); usess = unique(sess); typ = string(C.sniffType);
    fprintf('\n=== %d sniffs over %d sessions ===\n', numel(z), numel(usess));
    fprintf('sniffMaxZ quantiles: p50=%.1f p80=%.1f p90=%.1f p95=%.1f p99=%.1f max=%.1f\n', ...
        prctile(z,50),prctile(z,80),prctile(z,90),prctile(z,95),prctile(z,99),max(z));
    fprintf('\n  K   overall%%rej   worstSess%%   #sess@100%%   start%%   free%%   confirm%%\n');
    for K = [4 5 6 7 8 9 10 12 15 20]
        rej = z > K;
        perSess = arrayfun(@(s) mean(rej(sess==s))*100, usess);
        rs = pct(rej, typ=="start"); rf = pct(rej, typ=="free"); rc = pct(rej, typ=="confirm");
        fprintf('  %-4g  %8.1f     %8.1f    %6d     %6.1f  %6.1f   %6.1f\n', ...
            K, mean(rej)*100, max(perSess), sum(perSess>=100), rs, rf, rc);
    end
    fprintf(['\nPick K for >0%%, <=~20%% overall loss and no session at 100%%,\n' ...
             'then set it in o15_noise_K.m and rerun (run_o15_ztfr -> run_o15_gamma_epochs\n' ...
             '-> run_o15_task_group -> knit). Current o15_noise_K = %g.\n'], o15_noise_K());
end

function p = pct(rej, mask)
    if any(mask), p = mean(rej(mask))*100; else, p = NaN; end
end
