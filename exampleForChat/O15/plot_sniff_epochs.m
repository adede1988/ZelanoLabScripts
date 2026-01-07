function [epochs, t, hFig] = plot_sniff_epochs(outDat, R)
% plot_sniff_epochs  Epoch R.smoothR around finalOnset and plot all sniffs.
%
% [epochs, t, hFig] = plot_sniff_epochs(outDat, R)
%   outDat.behDat.finalOnset : sample indices (at outDat.fs)
%   R.smoothR                : 1 x T vector (same Fs as outDat.fs)
%
% Outputs:
%   epochs : nSniff x nSamp matrix (NaNs where an epoch would exceed bounds)
%   t      : 1 x nSamp time vector in seconds, from -1 to +2
%   hFig   : figure handle

    fs = outDat.fs;
    if ~isfield(outDat, 'behDat') || ~ismember('finalOnset', outDat.behDat.Properties.VariableNames)
        error('plot_sniff_epochs:NoFinalOnset', ...
              'outDat.behDat.finalOnset is required.');
    end

    % Gather onsets (clean)
    onsets = round(outDat.behDat.finalOnset(:));
    onsets = onsets(isfinite(onsets) & onsets > 0);

    % Epoch window in samples
    pre  = round(2 * fs);   % -1 s
    post = round(4 * fs);   % +2 s
    offs = (-pre:post);     % 1 x nSamp
    t    = offs / fs;       % time axis

    % Build index matrix for all sniffs (vectorized)
    idxMat = onsets + offs;                        % nSniff x nSamp
    T = numel(R.smoothR);
    valid = idxMat >= 1 & idxMat <= T;             % mask for in-bounds
    idxClamped = idxMat;
    idxClamped(~valid) = 1;                        % temporary valid index

    % Pull epochs and apply NaNs for out-of-bounds
    epochs = R.smoothR(idxClamped);                % nSniff x nSamp
    epochs(~valid) = NaN;

    % ----- Plot -----
    hFig = figure('Color','w'); hold on;
    plot(t, epochs', 'Color', [0.7 0.7 0.7]);      % all sniffs
    mu = nanmean(epochs, 1);
    plot(t, mu, 'k', 'LineWidth', 2);              % grand mean
    xline(0, '--');
    
    xlabel('Time (s)'); ylabel('Respiration (a.u., smoothed)');
    title(sprintf('Sniff-aligned respiration (%d epochs)', size(epochs,1)));
    box on; grid on;

    saveas(hFig,fullfile(outDat.figs, ...
                    ['allSniffs' '.jpg']));
end
