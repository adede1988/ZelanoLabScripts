function L = labPaths()
%LABPATHS  Single source of truth for all machine-specific paths.
%
%   L = labPaths() auto-detects the current machine (by Windows USERNAME, with
%   COMPUTERNAME available as a tiebreaker) and returns a struct holding every
%   path the preprocessing pipeline needs. The 7 deliverable scripts
%   (*_makeOutDat.m, *_main.m) and the loaders (applyParams / writeParams /
%   writePreProcX) all read their paths from here, so NOTHING else hard-codes a
%   machine path.
%
%   To run on a NEW machine you do exactly one thing: add a `case` to the switch
%   below (or drop in a labPaths_local.m, see "Local override"). Unknown machines
%   error with a copy-pasteable template.
%
%   ---- Base fields (machine-specific; the only things you set per machine) ----
%     .codePre    GitHub repos root that contains ZelanoLabScripts/,
%                 slowBreathing/, closed-loop-respiration/  (TRAILING filesep)
%     .eeglab     eeglab folder to addpath(genpath(...))
%     .labCommon  R:\...\Lab_Common base for shared lab data  (TRAILING filesep)
%     .gdrive     Google-Drive "My Drive" base (TRAILING filesep), used only for
%                 the breathing target-trace files; '' if not present here.
%
%   ---- Derived fields (built from the bases; never set these per machine) -----
%     .repo .eegLocCsv .slowBreathing .closedLoopResp .procBehavior .figPath
%     .adminXlsx .rootDupi .rootOBE .rootEEG .behCue .behThresh .targTraceDir
%     .labCommonCanon  (the Lab_Common prefix as stored in dataTracking.xlsx;
%                       lets loaders rebase the sheet's absolute paths onto a
%                       machine whose labCommon differs)
%
%   ---- Local override (optional, for a machine you don't want to commit) ----
%     If a file labPaths_local.m exists on the path it is used instead of the
%     switch. It must return ONLY the four base fields above; labPaths fills in
%     the derived ones. labPaths_local.m is git-ignored.

    % ---- optional untracked local override ----
    if exist('labPaths_local', 'file') == 2
        L = deriveLabPaths(labPaths_local());
        return;
    end

    user = lower(strtrim(getenv('USERNAME')));
    host = lower(strtrim(getenv('COMPUTERNAME')));

    switch user
        case 'adam'            % Adam's home desktop (COMPUTERNAME: B95)
            L.codePre   = 'C:\Users\Adam\Documents\GitHub\';
            L.eeglab    = 'C:\Users\Adam\Documents\eeglab2026.0.0';
            L.labCommon = 'R:\Neurology\Zelano_Lab\Lab_Common\';
            L.gdrive    = 'G:\My Drive\';

        case 'dtf8829'         % Northwestern lab workstation
            L.codePre   = 'G:\My Drive\GitHub\';
            L.eeglab    = 'C:\Users\dtf8829\Documents\eeglab2025.0.0';
            L.labCommon = 'R:\Neurology\Zelano_Lab\Lab_Common\';
            L.gdrive    = 'G:\My Drive\';

        otherwise
            error('labPaths:unknownMachine', ...
                ['Unknown machine (USERNAME="%s", COMPUTERNAME="%s").\n' ...
                 'Add a case to labPaths.m (paths must end with a filesep where shown):\n\n' ...
                 '    case ''%s''\n' ...
                 '        L.codePre   = ''<...>\\GitHub\\'';     %% holds ZelanoLabScripts, slowBreathing, closed-loop-respiration\n' ...
                 '        L.eeglab    = ''<...>\\eeglab2026.0.0'';\n' ...
                 '        L.labCommon = ''R:\\Neurology\\Zelano_Lab\\Lab_Common\\'';\n' ...
                 '        L.gdrive    = ''G:\\My Drive\\'';        %% '''''''' if Google Drive is not mapped here\n'], ...
                 user, host, user);
    end

    L = deriveLabPaths(L);
end

% ============================ helpers ============================

function L = deriveLabPaths(L)
% Fill in everything that is a function of the four base paths. Keeping this in
% one place means a new machine only ever sets the bases.

    reqd = {'codePre', 'eeglab', 'labCommon', 'gdrive'};
    for i = 1:numel(reqd)
        if ~isfield(L, reqd{i})
            error('labPaths:missingBase', ...
                'Base field "%s" is missing (check your labPaths case or labPaths_local).', reqd{i});
        end
    end

    cp = L.codePre;
    lc = L.labCommon;

    % code side (all under the GitHub repos root)
    L.repo           = fullfile(cp, 'ZelanoLabScripts');
    L.eegLocCsv      = fullfile(L.repo, 'eegLocs_standard_coords.csv');
    L.slowBreathing  = fullfile(cp, 'slowBreathing');
    L.closedLoopResp = fullfile(cp, 'closed-loop-respiration');
    L.procBehavior   = fullfile(L.closedLoopResp, 'processedBehavior');

    % lab-common data side. Roots used in raw string concatenation downstream
    % (e.g. [root sessID '\preProc\...']) so they MUST keep a trailing filesep.
    L.figPath   = [lc 'Adam\Dupi_processing\'];
    if ~isfield(L, 'adminXlsx')   % labPaths_local may point the tracker elsewhere
        L.adminXlsx = fullfile(lc, 'Admin', 'dataTracking.xlsx');
    end
    L.rootDupi  = [lc 'Dupi\'];
    L.rootOBE   = [lc 'OBEControl\'];
    L.rootEEG   = [lc 'AllStudyData\EEGbreathing\'];
    L.behCue    = fullfile(lc, 'OBE_task_backup', 'tasks', 'olf_cuetask', 'results', 'olf_cuetask_results');
    L.behThresh = fullfile(lc, 'OBE_task_backup', 'tasks', 'pea_threshold', 'results');

    % Google Drive side (breathing target traces only)
    if isempty(L.gdrive)
        L.targTraceDir = '';
    else
        L.targTraceDir = fullfile(L.gdrive, 'cZelano', 'breathingDataFiles');
    end

    % The Lab_Common prefix exactly as the absolute paths are stored in
    % dataTracking.xlsx (datPre column). Loaders rebase those paths from this
    % canonical prefix onto L.labCommon so the sheet stays portable.
    L.labCommonCanon = 'R:\Neurology\Zelano_Lab\Lab_Common\';
end
