function rsa_olf()
% Compute the three olfactory scores per session, for every session that has a
% breathing final. Formulas follow R_groupLevel\stanfordHelpers.R.

maxNumCompThreads(3);
P      = ohrv_config();
sp     = P.code;
outDir = P.work;
cfg  = ohrv_config();
root = cfg.dataRoot;

d = dir(fullfile(root, '*', 'preProc', '*breathingPreProc.mat'));
d = d(~contains({d.name}, 'separate'));
sess = erase({d.name}, '_breathingPreProc.mat');

% incremental: keep sessions already scored, only compute the new ones
csvF = fullfile(outDir, 'olfactory_scores.csv');
if exist(csvF, 'file')
    R = readtable(csvF, 'TextType', 'string');
    have = R.sessID;
else
    R = table(); have = strings(0, 1);
end

for ii = 1:numel(sess)
    sessID = sess{ii};
    pdir   = d(ii).folder;
    if ismember(string(sessID), have)
        fprintf('%-32s  cached\n', sessID); continue
    end
    row = table(string(sessID), NaN, NaN, NaN, NaN, NaN, NaN, ...
        'VariableNames', {'sessID','cue_d','cue_HR','cue_FA','thresh_low_cal','thresh_high_cal','O15_acc'});

    % ---- cue: d-prime with log-linear correction ----
    try
        b = loadBeh(pdir, [sessID '_cueTask*.mat']);
        if ~isempty(b) && all(ismember({'cue','odor','respString'}, b.Properties.VariableNames))
            b  = dedupeTrials(b);
            cu = double(b.cue); od = double(b.odor);
            rs = strtrim(string(b.respString));
            ok = isfinite(cu) & isfinite(od) & (rs == "Yes" | rs == "No");
            cu = cu(ok); od = od(ok); rs = rs(ok);
            hits = sum(cu == od & rs == "Yes");  miss = sum(cu == od & rs == "No");
            fas  = sum(cu ~= od & rs == "Yes");  crs  = sum(cu ~= od & rs == "No");
            nSig = hits + miss;  nNoi = fas + crs;
            if nSig > 0, row.cue_HR = hits / nSig; end
            if nNoi > 0, row.cue_FA = fas / nNoi;  end
            if nSig > 0 && nNoi > 0
                hrA = (hits + 0.5) / (nSig + 1);
                faA = (fas  + 0.5) / (nNoi + 1);
                row.cue_d = norminv(hrA) - norminv(faA);
            end
        end
    catch ME, fprintf('  %s cue: %s\n', sessID, ME.message); end

    % ---- thresh: mean intensity per odor level, air-calibrated ----
    try
        b = loadBeh(pdir, [sessID '_PEA_threshold*.mat']);
        if ~isempty(b) && all(ismember({'odor','intensity'}, b.Properties.VariableNames))
            b  = dedupeTrials(b);
            od = double(b.odor); it = double(b.intensity);
            m = @(v) mean(it(od == v & isfinite(it)), 'omitnan');
            none = m(1); low = m(2); high = m(3);
            if isfinite(none) && isfinite(low),  row.thresh_low_cal  = low  - none; end
            if isfinite(none) && isfinite(high), row.thresh_high_cal = high - none; end
        end
    catch ME, fprintf('  %s thresh: %s\n', sessID, ME.message); end

    % ---- O15: sum expScore over unique trials / 15 ----
    try
        b = loadBeh(pdir, [sessID '_O15*.mat']);
        if ~isempty(b) && all(ismember({'n','expScore'}, b.Properties.VariableNames))
            nn = double(b.n); es = double(b.expScore);
            ok = isfinite(nn); nn = nn(ok); es = es(ok);
            [un, ia] = unique(nn, 'stable');
            if ~isempty(un), row.O15_acc = sum(es(ia), 'omitnan') / 15; end
        end
    catch ME, fprintf('  %s O15: %s\n', sessID, ME.message); end

    fprintf('%-32s  d=%6.3f  low=%7.2f  high=%7.2f  O15=%5.3f\n', ...
        sessID, row.cue_d, row.thresh_low_cal, row.thresh_high_cal, row.O15_acc);
    R = [R; row]; %#ok<AGROW>
end

writetable(R, fullfile(outDir, 'olfactory_scores.csv'));
fprintf('OLF DONE  (%d sessions)\n', height(R));
end

function b = loadBeh(pdir, pat)
b = [];
g = dir(fullfile(pdir, pat));
g = g(~contains({g.name}, '.tmp'));
if isempty(g), return; end
S = load(fullfile(g(1).folder, g(1).name));
fn = fieldnames(S); pick = fn{1};
for k = 1:numel(fn)
    if ismember(fn{k}, {'outDat','chanDat','out'}), pick = fn{k}; break; end
end
D = S.(pick); clear S
if isfield(D, 'behDat'), b = D.behDat; end
end

function b = dedupeTrials(b)
% cue/thresh are one sniff per trial, but dedupe on n defensively
if ismember('n', b.Properties.VariableNames)
    [~, ia] = unique(double(b.n), 'stable');
    b = b(ia, :);
end
end
