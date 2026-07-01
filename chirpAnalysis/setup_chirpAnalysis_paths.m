function setup_chirpAnalysis_paths(verbose)
% SETUP_CHIRPANALYSIS_PATHS  Put the chirp-analysis code + its external tools on the path.
%   Idempotent. Resolves the GitHub root via labPaths (so it works on both the home
%   desktop C:\Users\Adam\Documents\GitHub\ and the lab desktop G:\My Drive\GitHub\).
%
%   Tools (cloned as siblings of ZelanoLabScripts under the GitHub root):
%     Superlets\matlab-pure   -> nfaslt (FASLT TFR; PURE matlab, no compiled dll)
%     mpact\mpc_analysis +
%     mpact\utility           -> MPACT adaptive chirplet transform (MPEM)
%     Synchrosqueezed-chirplet-transforms -> sqSTCT (Chen & Wu) [optional]
%   Ridge extraction uses MATLAB's licensed tfridge (no clone). MODA (Iatsenko ecurve)
%   is a fallback only and is NOT added by default to avoid name shadowing.

    if nargin < 1, verbose = false; end

    here = fileparts(mfilename('fullpath'));        % ...\ZelanoLabScripts\chirpAnalysis
    repo = fileparts(here);
    addpath(repo); addpath(here);

    % GitHub root (sibling of the repo) via labPaths, with a safe fallback
    codePre = '';
    try, L = labPaths(); codePre = char(L.codePre); catch, end
    if isempty(codePre) || ~isfolder(codePre), codePre = fileparts(repo); end

    tools = { ...
        fullfile(codePre,'Superlets','matlab-pure'), ...
        fullfile(codePre,'mpact','mpc_analysis'), ...
        fullfile(codePre,'mpact','utility'), ...
        fullfile(codePre,'Synchrosqueezed-chirplet-transforms') };
    for k = 1:numel(tools)
        if isfolder(tools{k}), addpath(tools{k});
        elseif verbose, fprintf('  [setup] MISSING tool dir: %s\n', tools{k}); end
    end

    % FieldTrip (for ft_preproc_bandpassfilter 'firws'); add + ft_defaults once
    try
        L = labPaths();
        if isfolder(L.fieldtrip) && isempty(which('ft_preproc_bandpassfilter'))
            addpath(L.fieldtrip); evalc('ft_defaults');
        end
    catch ME
        if verbose, fprintf('  [setup] fieldtrip not added: %s\n', ME.message); end
    end

    % sanity
    need = {'nfaslt','tfridge','hilbert','mp_act_signal'};
    for k = 1:numel(need)
        if isempty(which(need{k})) && verbose
            fprintf('  [setup] NOT on path: %s\n', need{k});
        end
    end
    if verbose, fprintf('  [setup] chirpAnalysis paths ready (codePre=%s)\n', codePre); end
end
