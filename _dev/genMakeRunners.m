function genMakeRunners()
% Generate _dev execution copies of the cue/thresh makeOutDat scripts:
%   * skip-check forced TRUE (always (re)process -> overwrite from raw)
%   * per-session body wrapped in try/catch (one failure won't abort the batch)
% The deliverable scripts stay untouched. Byte-level IO preserves CRLF.

    dev  = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\_dev';
    repo = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';

    % ---- cue ----
    s = rd(fullfile(repo, 'cueTask_makeOutDat.m'));
    s = strrep(s, 'if ~exist([datPre{datPrei(sessi)}', 'if true || exist([datPre{datPrei(sessi)}');
    s = wrapLoop(s, 'for sessi = 1:length(sessionIDs)', 'sessionIDs{sessi}');
    wr(fullfile(dev, 'run_cue_make.m'), s);
    fprintf('wrote run_cue_make.m\n');

    % ---- thresh ----
    s = rd(fullfile(repo, 'preproc', 'threshPreProc_makeOutDat.m'));
    s = strrep(s, 'if ~exist(fullfile(preProcDir', 'if true || exist(fullfile(preProcDir');
    s = wrapLoop(s, 'for sessi = 1:numel(sessionIDs)', 'sessID');
    wr(fullfile(dev, 'run_thresh_make.m'), s);
    fprintf('wrote run_thresh_make.m\n');
end

function s = wrapLoop(s, forLine, idExpr)
    assert(~isempty(strfind(s, forLine)), 'for-line missing: %s', forLine); %#ok<STREMP>
    % insert: try + optional RUNSUBSET filter (env var of space-separated indices)
    inject = [forLine sprintf(['\r\ntry\r\n' ...
        'if ~isempty(getenv(''RUNSUBSET'')) && ~ismember(sessi, str2num(getenv(''RUNSUBSET''))), continue; end'])];
    s = strrep(s, forLine, inject);
    t = deblank(s);
    assert(endsWith(t, 'end'), 'file does not end with end');
    pos = length(t) - 2;   % index of final 'e' in 'end'
    catchBlk = sprintf(['catch MErun\r\n' ...
        '    warning(''RUNFAIL %%s: %%s'', %s, MErun.message);\r\n' ...
        'end\r\n' ...
        'clear outDat dat dat1 dat2 comboDat photoDiode rawData behDat behDat1 behDat2 TTLs TTLs2 savedTTL data\r\n'], idExpr);
    s = [s(1:pos-1), catchBlk, s(pos:end)];
end

function s = rd(f)
    fid = fopen(f, 'r'); s = fread(fid, '*char')'; fclose(fid);
end
function wr(f, s)
    fid = fopen(f, 'w'); fwrite(fid, s); fclose(fid);
end
