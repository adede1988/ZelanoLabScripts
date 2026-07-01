function thresh_calibrate_noise()
% THRESH_CALIBRATE_NOISE  Calibrate the RELATIVE sharp-deflection threshold.
%   For each session's bestMac: d = 10ms max-min; robust z zd=(d-median)/(1.4826*MAD)
%   over the whole recording; per trial (epoch [-1.75,5.75]s about TTL.start) take
%   max(zd). Sweep K and report the dataset-wide fraction of trials that would be
%   rejected (max zd>K), plus the worst single-session fraction, so we can pick K
%   for >0%, <=20% loss. Writes _calib_trialMaxZ.csv under the thresh groupDir.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'threshAnalysis'));
    L = labPaths();
    groupDir = getenv('THRESH_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end

    epWin = [-1.75 5.75]; win = 5;   % 10ms window @500Hz
    T = thresh_session_table(false); T = T(T.onDisk, :);

    rows = {};   % sessID, group, trial, trialMaxZ
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
        ts = round(od.TTL.start);
        for t = 1:numel(ts)
            a = ts(t)+s0; b = ts(t)+s1; if a<1 || b>Tt, continue; end
            rows(end+1,:) = {id, char(T.group(i)), t, max(zd(a:b))}; %#ok<AGROW>
        end
        fprintf('%-26s med(d)=%.1f MAD-sd=%.1f  trials=%d\n', id, med, sg, numel(ts));
        clear od x d zd;
    end

    C = cell2table(rows, 'VariableNames', {'sessID','group','trial','trialMaxZ'});
    writetable(C, fullfile(groupDir, '_calib_trialMaxZ.csv'));

    z = C.trialMaxZ; sess = string(C.sessID); usess = unique(sess);
    fprintf('\n=== %d trials over %d sessions ===\n', numel(z), numel(usess));
    fprintf('trialMaxZ quantiles: p50=%.1f p80=%.1f p90=%.1f p95=%.1f p99=%.1f max=%.1f\n', ...
        prctile(z,50),prctile(z,80),prctile(z,90),prctile(z,95),prctile(z,99),max(z));
    fprintf('\n  K   overall%%rej   worstSess%%   #sess@100%%\n');
    for K = [4 5 6 7 8 9 10 12 15 20]
        rej = z>K;
        perSess = arrayfun(@(s) mean(rej(sess==s))*100, usess);
        fprintf('  %-4g  %8.1f     %8.1f    %6d\n', K, mean(rej)*100, max(perSess), sum(perSess>=100));
    end
end
