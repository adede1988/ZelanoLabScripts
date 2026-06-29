function o15_init_paths()
% O15_INIT_PATHS  Put EEGLAB + FieldTrip + repo (+ cueAnalysis, O15Analysis) on
%   the path for the O15 analysis. Mirrors cue_init_paths exactly (same EEGLAB-
%   then-FieldTrip ordering for the shared FOOOF/newtimef names); also adds
%   cueAnalysis so the O15 drivers can REUSE cue_fooof_macBP / cue_noise_trials /
%   cue_plot_* unchanged. Run once per MATLAB session.

    % self-bootstrap the repo onto the path so labPaths resolves regardless of
    % the current folder (this file lives in <repo>\O15Analysis).
    here = fileparts(mfilename('fullpath'));      % ...\ZelanoLabScripts\O15Analysis
    repo = fileparts(here);
    addpath(repo); addpath(fullfile(repo,'cueAnalysis')); addpath(here);

    L = labPaths();
    addpath(L.repo);
    addpath(fullfile(L.repo, 'cueAnalysis'));     % reuse cue_fooof_macBP, cue_noise_trials, cue_plot_*
    addpath(fullfile(L.repo, 'O15Analysis'));

    % --- EEGLAB (headless) ---
    if exist('eeglab','file') ~= 2
        addpath(L.eeglab);
    end
    eeglab nogui;

    % --- FieldTrip (machine-aware via labPaths; add AFTER eeglab so it wins shared names) ---
    ftDir = L.fieldtrip;
    if isempty(ftDir) || ~isfolder(ftDir)
        error('o15_init_paths:noFieldTrip', ...
            'FieldTrip not found at "%s" (set L.fieldtrip for this machine in labPaths.m)', ftDir);
    end
    addpath(ftDir);
    ft_defaults;
    bdir = fullfile(ftDir, 'external', 'brainstorm');   % brainstorm process_fooof
    if isfolder(bdir), addpath(bdir); end
    if exist('newtimef','file') ~= 2 || exist('eeglab2fieldtrip','file') ~= 2
        addpath(genpath(L.eeglab));
    end
end
