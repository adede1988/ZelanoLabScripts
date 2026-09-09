function genReprocMake()
% Generate portable, batch-safe make-runners (raw -> intermediate) for the
% reprocessing campaign, from the deliverable *_makeOutDat scripts:
%   * skip-guard forced TRUE  -> always (re)process, overwrite from raw
%   * per-session try/catch    -> one failure won't abort the batch
%   * RUNSUBSET env filter      -> process only chosen session indices
% The deliverables use labPaths(), so the generated copies are portable.
% Output: _reprocess/reproc_{cue,thresh,breathing}_make.m  (CRLF preserved).

    repo   = fileparts(fileparts(mfilename('fullpath')));   % _reprocess/ is one level down
    outdir = fullfile(repo, '_reprocess');

    % ---- cue ----
    s = rd(fullfile(repo, 'cueTask_makeOutDat.m'));
    s = strrep(s, 'if ~exist([datPre{datPrei(sessi)}', 'if true || exist([datPre{datPrei(sessi)}');
    s = wrapLoop(s, 'for sessi = 1:length(sessionIDs)', 'sessionIDs{sessi}');
    wr(fullfile(outdir, 'reproc_cue_make.m'), s);
    fprintf('wrote reproc_cue_make.m\n');

    % ---- thresh ----
    s = rd(fullfile(repo, 'preproc', 'threshPreProc_makeOutDat.m'));
    s = strrep(s, 'if ~exist(fullfile(preProcDir', 'if true || exist(fullfile(preProcDir');
    s = wrapLoop(s, 'for sessi = 1:numel(sessionIDs)', 'sessID');
    wr(fullfile(outdir, 'reproc_thresh_make.m'), s);
    fprintf('wrote reproc_thresh_make.m\n');

    % ---- breathing (already has per-iteration try/catch; convert parfor->for,
    %      force guard, inject the RUNSUBSET filter) ----
    s = rd(fullfile(repo, 'breathingTask_makeOutDat.m'));
    s = replaceFirst(s, 'if ~exist([datPre{datPrei(sessi)}', 'if true || exist([datPre{datPrei(sessi)}');  % skip-guard only, not the mkdir guard
    fromP = 'parfor sessi = 1:length(sessionIDs)';
    toP   = ['for sessi = 1:length(sessionIDs)' char(13) char(10) ...
             'if ~isempty(getenv(''RUNSUBSET'')) && ~ismember(sessi, str2num(getenv(''RUNSUBSET''))), continue; end'];
    assert(contains(s, fromP), 'breathing parfor line missing');
    s = strrep(s, fromP, toP);
    wr(fullfile(outdir, 'reproc_breathing_make.m'), s);
    fprintf('wrote reproc_breathing_make.m\n');
end

function s = wrapLoop(s, forLine, idExpr)
    assert(contains(s, forLine), 'for-line missing: %s', forLine);
    inject = [forLine char(13) char(10) 'try' char(13) char(10) ...
        'if ~isempty(getenv(''RUNSUBSET'')) && ~ismember(sessi, str2num(getenv(''RUNSUBSET''))), continue; end'];
    s = strrep(s, forLine, inject);
    t = deblank(s);
    assert(endsWith(t, 'end'), 'file does not end with end');
    pos = length(t) - 2;   % index of final 'e' in 'end'
    catchBlk = sprintf(['catch MErun' char(13) char(10) ...
        '    warning(''RUNFAIL %%s: %%s'', %s, MErun.message);' char(13) char(10) ...
        'end' char(13) char(10) ...
        'clear outDat dat dat1 dat2 comboDat photoDiode rawData behDat behDat1 behDat2 TTLs TTLs2 savedTTL data' char(13) char(10)], idExpr);
    s = [s(1:pos-1), catchBlk, s(pos:end)];
end

function s = replaceFirst(s, from, to)
    i = strfind(s, from);
    assert(~isempty(i), 'pattern missing: %s', from);
    i = i(1);
    s = [s(1:i-1), to, s(i+length(from):end)];
end

function s = rd(f), fid = fopen(f, 'r'); s = fread(fid, '*char')'; fclose(fid); end
function wr(f, s), fid = fopen(f, 'w'); fwrite(fid, s); fclose(fid); end
