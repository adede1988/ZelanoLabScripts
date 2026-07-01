function thresh_plot_singletrial(NT, sessID, bestMac, sep, outPng)
% THRESH_PLOT_SINGLETRIAL  Stacked single-trial raw plot; NOISY trials drawn red.
%   NT = output of thresh_noise_trials. Traces offset `sep` uV apart; y-ticks are
%   original trial numbers; TTL.start at the blue line.

    okIdx = find(NT.ok); nT = numel(okIdx); times = NT.times;
    f = figure('visible','off','position',[0 0 1100 1000]); hold on;
    for r = 1:nT
        i = okIdx(r); y = NT.D(:,i) + (r-1)*sep;
        if NT.noisy(i), c = [0.85 0 0]; lw = 0.9; else, c = [0.10 0.10 0.10]; lw = 0.4; end
        plot(times, y, 'color', c, 'linewidth', lw);
    end
    xline(0, 'b', 'linewidth', 1);
    ridx = 1:5:nT; yt = (ridx-1)*sep; ytl = arrayfun(@(r) num2str(okIdx(r)), ridx, 'uni', 0);
    set(gca, 'ytick', yt, 'yticklabel', ytl);
    ylim([-sep, nT*sep]); xlim([times(1) times(end)]);
    nNoise = sum(NT.noisy(okIdx));
    xlabel('Time from TTL.start (ms)');
    ylabel(sprintf('trial   (offset %d uV apart; red = noise)', sep));
    title(sprintf('%s   %s   single-trial raw macBP (n=%d, %d noise: relative zd>%.1f)', ...
        sessID, bestMac, nT, nNoise, NT.K), 'interpreter','none');
    saveas(f, outPng); close(f);
end
