function rsa_extract()
% Extract the minimum needed for the respiration->HRV transfer function analysis
% from each finished breathing final, one session per iteration with clears.

maxNumCompThreads(3);
P      = ohrv_config();
outDir = P.work;
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg  = ohrv_config();

root = cfg.dataRoot;
d = dir(fullfile(root, '*', 'preProc', '*breathingPreProc.mat'));
d = d(~contains({d.name}, 'separate'));
fprintf('found %d breathing finals\n', numel(d));

for ii = 1:numel(d)
    f = fullfile(d(ii).folder, d(ii).name);
    [~, base] = fileparts(d(ii).name);
    sessID = erase(base, '_breathingPreProc');
    outF = fullfile(outDir, [sessID '_slim.mat']);
    if exist(outF, 'file'), fprintf('[%2d/%d] %s  SKIP (done)\n', ii, numel(d), sessID); continue; end

    t0 = tic;
    try
        S = load(f);
        fn = fieldnames(S);
        pick = fn{1};
        for k = 1:numel(fn)
            if ismember(fn{k}, {'outDat','chanDat','out'}), pick = fn{k}; break; end
        end
        D = S.(pick);
        clear S

        if ~isfield(D, 'behDat') || ~isfield(D, 'bmFeatures')
            fprintf('[%2d/%d] %s  SKIP (unfinished: no behDat/bmFeatures)\n', ii, numel(d), sessID);
            clear D; continue;
        end

        b  = D.behDat;
        fs = double(D.fs);

        isRR  = cellfun(@(x) contains(x, 'RRint'), D.labels);
        if ~any(isRR)
            fprintf('[%2d/%d] %s  SKIP (no RRint channel)\n', ii, numel(d), sessID);
            clear D; continue;
        end
        rrint = double(D.data(find(isRR, 1), :));
        nSamp = numel(rrint);
        clear D

        vn = b.Properties.VariableNames;
        need = {'goodBreath','RR_max_min','maxRR','minRR','length','finalOnset', ...
                'bm_inhaleDurations','bm_inhaleVolumes','condition','task','noseMouth'};
        miss = need(~ismember(need, vn));
        if ~isempty(miss)
            fprintf('[%2d/%d] %s  SKIP (missing: %s)\n', ii, numel(d), sessID, strjoin(miss, ','));
            clear b rrint; continue;
        end

        T = table();
        T.goodBreath  = double(b.goodBreath);
        T.RR_max_min  = double(b.RR_max_min);
        T.maxRR       = double(b.maxRR);
        T.minRR       = double(b.minRR);
        T.len         = double(b.length);
        T.finalOnset  = double(b.finalOnset);
        T.inhDur      = double(b.bm_inhaleDurations);
        T.inhVol      = double(b.bm_inhaleVolumes);
        T.condition   = double(b.condition);
        T.task        = safeStr(b.task);
        T.noseMouth   = safeStr(b.noseMouth);

        save(outF, 'T', 'rrint', 'fs', 'nSamp', 'sessID', '-v7');
        fprintf('[%2d/%d] %s  ok  nBreaths=%d  nSamp=%d  (%.1f s)\n', ...
            ii, numel(d), sessID, height(T), nSamp, toc(t0));
        clear T b rrint
    catch ME
        fprintf('[%2d/%d] %s  ERROR: %s\n', ii, numel(d), sessID, ME.message);
    end
end

fprintf('EXTRACT DONE\n');
end

function s = safeStr(v)
% Element-wise conversion that tolerates empties, nested cells and non-char
% entries (some sessions carry a stray non-string in behDat.task).
if isstring(v) || ischar(v) || iscategorical(v)
    s = string(v); s = s(:); return
end
n = numel(v); s = strings(n, 1);
for k = 1:n
    try
        e = v{k};
        while iscell(e) && ~isempty(e), e = e{1}; end
        if ischar(e) || isstring(e)
            s(k) = string(e);
        elseif isnumeric(e) && isscalar(e)
            s(k) = string(e);
        else
            s(k) = "";
        end
    catch
        s(k) = "";
    end
end
end
