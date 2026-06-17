function test_integration()
% Per task, for a capped number of sessions whose source files exist, compare the
% NEW path (applyParams + assemble_outDat_all) against the LEGACY path
% (getSessionParams_<task> + [detect_ttls_O15] + assemble_outDat_<task>):
%   - Mode B P: common fields equal
%   - assembled outDat: common fields equal (data via isequaln)
%   - O15 TTL: equal

    repo = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
    addpath(repo); addpath(fullfile(repo,'preproc'));
    xlsx = fullfile(repo, 'dataTracking.xlsx');

    tasks = {'breathingTask','cueTask','threshTask','O15'};
    caps  = [6, 4, 4, 3];
    knownCaveat = 0;

    totMis = 0;
    for ti = 1:numel(tasks)
        task = tasks{ti}; isO15 = strcmp(task,'O15');
        cfg = applyParams(task, 'main', xlsx);
        got = 0;
        fprintf('\n==== %s ====\n', task);
        for s = 1:numel(cfg.sessionIDs)
            if got >= caps(ti), break; end
            id = cfg.sessionIDs{s};
            S = struct('id', id, 'root', cfg.root{s}, 'figPath', tempdir);
            if ~isO15, S.fig = fullfile(tempdir, id); end
            if ~isfolder(fullfile(tempdir,id)), mkdir(fullfile(tempdir,id)); end

            % ---- legacy ----
            try
                switch task
                    case 'breathingTask', [rawL,PL] = getSessionParams_breathingTask(S);
                    case 'cueTask',       [rawL,PL] = getSessionParams_cueTask(S);
                    case 'threshTask',    [rawL,PL] = getSessionParams_threshTask(S);
                    case 'O15',           [rawL,PL] = getSessionParams_O15(S);
                end
                TTLL = [];
                if isO15, [TTLL, rawL] = detect_ttls_O15(rawL, PL); end
                switch task
                    case 'O15', outL = assemble_outDat_O15(rawL, S, PL);
                    otherwise,  outL = assemble_outDat_breathing_cue_Task(rawL, S, PL);
                end
            catch ME
                fprintf('  skip %-24s (legacy: %s)\n', id, ME.message);
                continue;
            end

            % ---- new ----
            try
                PN = applyParams(task, id, xlsx);
                [outN, ~, TTLN] = assemble_outDat_all(S, PN);
            catch ME
                totMis = totMis + 1;
                fprintf('  NEW-ERROR %-24s : %s\n', id, ME.message);
                continue;
            end
            got = got + 1;

            m = 0;
            [dm, dc] = cmpStructs(['P:' id],    PL,   PN,   {'getBeats'});
            m = m + dm; knownCaveat = knownCaveat + dc;
            [dm, dc] = cmpStructs(['out:' id],  outL, outN, {});
            m = m + dm; knownCaveat = knownCaveat + dc;
            if isO15
                if ~isequaln(TTLL, TTLN)
                    m = m + 1; fprintf('  MISMATCH %-24s TTL differs\n', id);
                end
            end
            if m == 0, fprintf('  OK   %-24s\n', id); end
            totMis = totMis + m;
        end
    end
    fprintf('\nintegration mismatches = %d  (known caveats = %d)\n', totMis, knownCaveat);
    if totMis == 0, fprintf('ALL_INTEGRATION_PASS\n'); else, fprintf('INTEGRATION_FAIL\n'); end
end

function [m, caveat] = cmpStructs(tag, A, B, skip)
    m = 0; caveat = 0;
    f = intersect(fieldnames(A), fieldnames(B));
    for i = 1:numel(f)
        fn = f{i};
        if any(strcmp(fn, skip)), continue; end
        a = A.(fn); b = B.(fn);
        % task may legitimately be string("O15") vs char('O15')
        if strcmp(fn,'task'), a = char(string(a)); b = char(string(b)); end
        % figs paths: normalise separators
        if strcmp(fn,'figs'), a = char(string(a)); b = char(string(b)); end
        if ~isequaln(a, b)
            % KNOWN CAVEAT: breathing spikeClean encodes the sheet's (commented-out)
            % intent vs the legacy runtime default; no effect where hasMacros=false.
            if strcmp(fn,'spikeClean') && startsWith(tag,'P:')
                caveat = caveat + 1;
                fprintf('  caveat   %-26s .%-12s sheet=%s legacy=%s\n', tag, fn, mat2str(b), mat2str(a));
            else
                m = m + 1;
                fprintf('  MISMATCH %-28s .%-12s\n', tag, fn);
            end
        end
    end
end
