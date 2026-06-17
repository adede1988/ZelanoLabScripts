function syncAdmin()
% Safely add the parameter columns (28..44) from the repo's param-enriched
% dataTracking.xlsx into the lab's Admin master, writing ONLY those columns of
% Sheet1 (preserving every other column and sheet). Verifies row alignment by
% Subject ID before writing, and backs up the Admin file first.

    repoXlsx  = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\dataTracking.xlsx';
    adminXlsx = 'R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx';

    R = readcell(repoXlsx,  'Sheet', 'Sheet1');
    A = readcell(adminXlsx, 'Sheet', 'Sheet1');
    fprintf('repo  %dx%d\nadmin %dx%d\n', size(R,1), size(R,2), size(A,1), size(A,2));

    nr = min(size(R,1), size(A,1));
    % verify Subject ID alignment (col 1), rows 3..end (data)
    nMis = 0;
    for r = 3:nr
        a = cellStr(A{r,1}); b = cellStr(R{r,1});
        if ~strcmpi(strtrim(a), strtrim(b))
            nMis = nMis + 1;
            if nMis <= 5, fprintf('  row %d ID mismatch: admin="%s" repo="%s"\n', r, a, b); end
        end
    end
    if size(R,1) ~= size(A,1)
        fprintf('ROW COUNT DIFFERS -> aborting sync.\n'); return;
    end
    if nMis > 0
        fprintf('ABORT: %d Subject ID misalignments; will NOT write Admin.\n', nMis); return;
    end
    fprintf('Subject IDs aligned across all %d rows.\n', nr-2);

    % backup admin
    bak = strrep(adminXlsx, '.xlsx', '_backup_preParamSync.xlsx');
    if ~exist(bak,'file')
        copyfile(adminXlsx, bak);
        fprintf('backed up Admin -> %s\n', bak);
    else
        fprintf('backup already exists: %s\n', bak);
    end

    % write repo param columns (28..44), rows 2..end (header + data) to Admin AB2
    block = R(2:end, 28:44);
    writecell(block, adminXlsx, 'Sheet', 'Sheet1', 'Range', 'AB2', ...
              'AutoFitWidth', false);
    fprintf('wrote %dx%d param block to Admin Sheet1!AB2\n', size(block,1), size(block,2));

    % verify: applyParams should now work against the Admin default
    clear functions %#ok<CLFUNC>  % drop applyParams sheet cache
    cfg = applyParams('breathingTask','main', adminXlsx);
    P   = applyParams('O15', '250818_Dupi_NMH_JH_1', adminXlsx);
    fprintf('VERIFY admin: breathing sessions=%d ; JH_1 O15 respThresh=%g hasEEG=%d\n', ...
        numel(cfg.sessionIDs), P.respThresh, P.hasEEG);
    fprintf('ADMIN_SYNC_OK\n');
end

function s = cellStr(v)
    if isa(v,'missing'), s=''; return; end
    if ischar(v), s=v; return; end
    if isstring(v), s=char(v); return; end
    if isnumeric(v)||islogical(v), s=num2str(v); return; end
    s = char(string(v));
end
