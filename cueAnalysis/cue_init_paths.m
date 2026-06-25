function cue_init_paths()
% CUE_INIT_PATHS  Put EEGLAB + FieldTrip + repo on the path for the cue analysis.
%   Order matters: start EEGLAB (gives eeglab2fieldtrip, newtimef), then add the
%   real FieldTrip LAST so it shadows EEGLAB's bundled Fieldtrip-lite for the
%   handful of shared function names. Run once per MATLAB session.

    L = labPaths();
    addpath(L.repo);
    addpath(fullfile(L.repo, 'cueAnalysis'));

    % --- EEGLAB (headless) ---
    if exist('eeglab','file') ~= 2
        addpath(L.eeglab);
    end
    eeglab nogui;

    % --- FieldTrip (this machine; add AFTER eeglab so it wins shared names) ---
    ftDir = 'C:\Users\Adam\Documents\fieldtrip-20260518';
    if ~isfolder(ftDir)
        error('cue_init_paths:noFieldTrip', 'FieldTrip not found at %s', ftDir);
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
