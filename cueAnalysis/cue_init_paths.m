function cue_init_paths()
% CUE_INIT_PATHS  Put EEGLAB + FieldTrip + repo on the path for the cue analysis.
%   Order matters: start EEGLAB (gives eeglab2fieldtrip, newtimef), then add the
%   real FieldTrip LAST so it shadows EEGLAB's bundled Fieldtrip-lite for the
%   handful of shared function names. Run once per MATLAB session.

    L = labPaths();
    addpath(L.repo);
    addpath(fullfile(L.repo, 'cueAnalysis'));

    % --- EEGLAB (headless) ---
    addpath(L.eeglab);   % unconditional: a stale/saved eeglab path entry in a detached
                         % (WMI) MATLAB makes an exist() guard skip this, then `eeglab`
                         % resolves to a broken reference and errors. Prepending always wins.
    eeglab nogui;

    % --- FieldTrip (machine-aware via labPaths; add AFTER eeglab so it wins shared names) ---
    ftDir = L.fieldtrip;
    if isempty(ftDir) || ~isfolder(ftDir)
        error('cue_init_paths:noFieldTrip', ...
            'FieldTrip not found at "%s" (set L.fieldtrip for this machine in labPaths.m)', ftDir);
    end
    addpath(ftDir);
    ft_defaults;
    % brainstorm FOOOF (process_fooof) is normally added on-demand by
    % ft_freqanalysis; add it now so it can be called directly too.
    bdir = fullfile(ftDir, 'external', 'brainstorm');
    if isfolder(bdir), addpath(bdir); end
    % ft_defaults can prepend external/ dirs that shadow EEGLAB; make sure the
    % EEGLAB functions we need are still reachable.
    if exist('newtimef','file') ~= 2 || exist('eeglab2fieldtrip','file') ~= 2
        addpath(genpath(L.eeglab));
    end
end
