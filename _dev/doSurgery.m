function doSurgery()
% Reliable anchor-based surgery for cueTask_makeOutDat.m and
% threshPreProc_makeOutDat.m: fix machine paths, delete the hard-coded array
% blocks, and insert the applyParams drop-in block before the main loop.
% Byte-level read/write preserves CRLF line endings.

    repo = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
    NL = sprintf('\r\n');

    % ---------------- cueTask_makeOutDat.m ----------------
    f = fullfile(repo, 'cueTask_makeOutDat.m');
    s = rd(f);
    s = strrep(s, 'G:\My Drive\GitHub\', 'C:\Users\Adam\Documents\GitHub\');
    s = strrep(s, 'C:\Users\dtf8829\Documents\eeglab2025.0.0', ...
                  'C:\Users\Adam\Documents\eeglab2026.0.0');
    s = delSpan(s, 'datPre = {', '%hard code flip');
    cueCfg = sprintf(['cfg        = applyParams(''cueTask'',''makeOutDat'');\r\n' ...
                      'sessionIDs = cfg.sessionIDs;\r\n' ...
                      'datPre     = cfg.datPre;\r\n' ...
                      'datPrei    = cfg.datPrei;\r\n' ...
                      'newSet     = cfg.newIDs;\r\n' ...
                      'rspIDX     = cfg.rspIDX;\r\n' ...
                      'rspFlip    = cfg.rspFlip;\r\n\r\n']);
    s = insertBefore1(s, 'for sessi = 1:length(sessionIDs)', cueCfg);
    wr(f, s);
    fprintf('cue done\n');

    % ---------------- threshPreProc_makeOutDat.m ----------------
    f = fullfile(repo, 'preproc', 'threshPreProc_makeOutDat.m');
    s = rd(f);
    s = strrep(s, 'G:\My Drive\GitHub\', 'C:\Users\Adam\Documents\GitHub\');
    s = strrep(s, 'C:\Users\dtf8829\Documents\eeglab2025.0.0', ...
                  'C:\Users\Adam\Documents\eeglab2026.0.0');
    % delete datPre block (keep behDatPath_newSet which follows it)
    s = delSpan(s, 'datPre    = {', '};');
    % delete datPrei..rspFlip (comment header through rspFlip line)
    s = delSpan(s, '% datPre index for each session', '%hard code flip');
    thrCfg = sprintf(['cfg        = applyParams(''threshTask'',''makeOutDat'');\r\n' ...
                      'sessionIDs = cfg.sessionIDs;\r\n' ...
                      'datPre     = cfg.datPre;\r\n' ...
                      'datPrei    = cfg.datPrei;\r\n' ...
                      'newSet     = cfg.newIDs;\r\n' ...
                      'rspIDX     = cfg.rspIDX;\r\n' ...
                      'rspFlip    = cfg.rspFlip;\r\n\r\n']);
    s = insertBefore1(s, 'for sessi = 1:numel(sessionIDs)', thrCfg);
    wr(f, s);
    fprintf('thresh done\n');
end

function s = rd(f)
    fid = fopen(f, 'r'); s = fread(fid, '*char')'; fclose(fid);
end
function wr(f, s)
    fid = fopen(f, 'w'); fwrite(fid, s); fclose(fid);
end

function raw = delSpan(raw, a, b)
    i1 = strfind(raw, a);
    assert(~isempty(i1), 'start anchor missing: %s', a);
    i1 = i1(1);
    bs = strfind(raw, b);
    bs = bs(bs >= i1);
    assert(~isempty(bs), 'end anchor missing after start: %s', b);
    i2 = bs(1) + length(b) - 1;
    j = i2;
    while j <= numel(raw) && raw(j) ~= char(10)
        j = j + 1;
    end
    if j <= numel(raw), j = j; end   % j now at the \n (or past end)
    raw = [raw(1:i1-1), raw(min(j+1,numel(raw)+1):end)];
end

function s = insertBefore1(s, anchor, ins)
    i = strfind(s, anchor);
    assert(~isempty(i), 'insert anchor missing: %s', anchor);
    i = i(1);
    s = [s(1:i-1), ins, s(i:end)];
end
