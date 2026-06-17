function test_applyParams()
% Validate applyParams against the legacy getSessionParams_<task> loaders.
%   - Mode B P parity: every common field must match (per session, where the
%     source .mat exists so the legacy loader can run).
%   - getBeats cross-check: the sheet beatSpec must reproduce the exact legacy
%     getBeats_<name> assigned to that session in the breathing switch.
%   - Mode A cfg sanity.

    repo = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
    addpath(repo);
    xlsx = fullfile(repo, 'dataTracking.xlsx');   % the param-enriched copy

    tasks   = {'breathingTask','cueTask','threshTask','O15'};
    loaders = {@getSessionParams_breathingTask, @getSessionParams_cueTask, ...
               @getSessionParams_threshTask, @getSessionParams_O15};

    rng(11); ECGz = randn(3, 120000); ECGz(:,1:40:end) = ECGz(:,1:40:end)*4;

    totMis = 0; totCmp = 0; totSkip = 0;
    for ti = 1:numel(tasks)
        task = tasks{ti};
        cfg  = applyParams(task, 'main', xlsx);
        fprintf('\n==== %s : %d sessions ====\n', task, numel(cfg.sessionIDs));
        for s = 1:numel(cfg.sessionIDs)
            id = cfg.sessionIDs{s};
            S = struct('id', id, 'root', cfg.root{s}, ...
                       'figPath', tempdir, 'fig', fullfile(tempdir, id));
            % legacy loader (needs the source .mat) -------------------------
            try
                [~, Pleg] = loaders{ti}(S);
            catch ME
                totSkip = totSkip + 1;
                continue;   % .mat not present -> can't run oracle; skip silently
            end
            Pnew = applyParams(task, id, xlsx);
            totCmp = totCmp + 1;

            % compare common scalar/array fields ---------------------------
            f = intersect(fieldnames(Pleg), fieldnames(Pnew));
            for k = 1:numel(f)
                fn = f{k};
                if strcmp(fn, 'getBeats'), continue; end   % handled below
                a = Pleg.(fn); b = Pnew.(fn);
                if ~isequaln(a, b)
                    totMis = totMis + 1;
                    fprintf('  MISMATCH %-24s .%-14s legacy=%s new=%s\n', ...
                        id, fn, briefv(a), briefv(b));
                end
            end

            % getBeats cross-check (breathing) -----------------------------
            if isfield(Pleg, 'getBeats') && isfield(Pnew, 'getBeats')
                gl = Pleg.getBeats(ECGz, 200);
                gn = Pnew.getBeats(ECGz, 200);
                if ~isequal(gl, gn)
                    totMis = totMis + 1;
                    fprintf('  MISMATCH %-24s .getBeats legacy=%d new=%d beats\n', ...
                        id, numel(gl), numel(gn));
                end
            end
        end
    end

    fprintf('\n--- Mode A cfg sanity ---\n');
    for ti = 1:numel(tasks)
        cfg = applyParams(tasks{ti}, 'makeOutDat', xlsx);
        assert(numel(cfg.sessionIDs) == numel(unique(lower(cfg.sessionIDs))), ...
            'duplicate sessionIDs in %s', tasks{ti});
        assert(all(cfg.datPrei >= 1 & cfg.datPrei <= numel(cfg.datPre)), ...
            'datPrei out of range in %s', tasks{ti});
        assert(numel(cfg.datPre) >= 3, 'datPre missing canonical roots');
        fprintf('%-14s n=%-3d datPre=%d roots, newIDs=%d, datPrei in [%d %d]\n', ...
            tasks{ti}, numel(cfg.sessionIDs), numel(cfg.datPre), ...
            numel(cfg.newIDs), min(cfg.datPrei), max(cfg.datPrei));
    end

    fprintf('\napplyParams: compared=%d sessions, skipped(no .mat)=%d, mismatches=%d\n', ...
        totCmp, totSkip, totMis);
    if totMis == 0
        fprintf('ALL_APPLYPARAMS_PASS\n');
    else
        fprintf('APPLYPARAMS_FAIL count=%d\n', totMis);
    end
end

function s = briefv(v)
    if islogical(v), s = mat2str(v); return; end
    if isnumeric(v)
        if isempty(v), s = '[]'; elseif isscalar(v), s = num2str(v); else, s = mat2str(v); end
        return;
    end
    if ischar(v), s = ['''' v '''']; return; end
    if isstring(v), s = ['"' char(v) '"']; return; end
    s = ['<' class(v) '>'];
end
