function run_o15_singletrials(sessFilter)
% RUN_O15_SINGLETRIALS  Generate the per-trial multichannel noise-diagnostic
%   plots (o15_plot_session_singletrials) for each fresh O15 session (or a
%   subset). STANDALONE diagnostic -- NOT part of run_O15Analysis_all; it reads
%   the finals and writes only <figs>\<id>\O15\singleTrials\trial_NN.png. Use it
%   to eyeball the RELATIVE noise detector (zd>K, K=o15_noise_K) against the raw
%   macBP signal across every trial before committing to a K and rerunning.
%
%   run_o15_singletrials()            % all fresh on-disk O15 sessions
%   run_o15_singletrials({'<id>'})    % one or a subset

    if nargin < 1, sessFilter = []; end
    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'cueAnalysis')); addpath(fullfile(repo,'O15Analysis'));

    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    for i = 1:height(T)
        id = char(T.sessID(i));
        fprintf('\n== %d/%d %s ==\n', i, height(T), id);
        try
            o15_plot_session_singletrials(id);
        catch ME
            fprintf('  FAILED: %s\n', ME.message);
        end
    end
    fprintf('\nO15 per-trial single-trial diagnostic pass done.\n');
end
