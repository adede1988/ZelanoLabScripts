function v2_figs_phase(pp, figDir, sessID, C)
% V2_FIGS_PHASE  Per-subject phase-continuity + burst figures (spec chirpAnalysis_3.md sec 4).
%   Phase: ONE set of axes, peak-centered (+-1000 ms from gammaPeakTime), peak (blue) vs ridge
%   (red) consistency traces per trial + thicker across-trial means. Burst: overlaid histograms
%   of peakBurstLength and ridgeBurstLength. Overwrites prior versions.

    if ~isfolder(figDir), mkdir(figDir); end
    fs = C.fs;
    [Mpk, relT] = v2_peakcenter(pp.peakPhaseConsistency,  pp.tMs, pp.gammaPeakTime, 1000, fs);
    [Mrg, ~   ] = v2_peakcenter(pp.ridgePhaseConsistency, pp.tMs, pp.gammaPeakTime, 1000, fs);
    good = find(~isnan(pp.gammaPeakTime));
    if isempty(good), return; end

    % --- phase continuity: single axes, peak-centered, peak vs ridge ---
    fh = figure('Visible','off','Position',[100 100 780 540]); hold on;
    for i = good', plot(relT, abs(Mpk(i,:)), '-', 'Color',[0 0.3 0.8 0.10]); end
    for i = good', plot(relT, abs(Mrg(i,:)), '-', 'Color',[0.85 0.2 0 0.10]); end
    hP = plot(relT, mean(abs(Mpk(good,:)),1,'omitnan'), '-', 'Color',[0 0.2 0.7], 'LineWidth',2.6);
    hR = plot(relT, mean(abs(Mrg(good,:)),1,'omitnan'), '-', 'Color',[0.75 0.1 0], 'LineWidth',2.6);
    yline(C.pp.threshNarrow,'k:'); yline(C.pp.threshWide,'k--');
    xlim([-1000 1000]); ylim([0 pi]);
    xlabel('time from gammaPeakTime (ms)'); ylabel('|phase diff| (rad)');
    legend([hP hR], {'peakPhaseConsistency (mean)','ridgePhaseConsistency (mean)'}, 'Location','north');
    title(sprintf('%s  phase continuity (peak-centered, n=%d)', strrep(sessID,'_','\_'), numel(good)),'Interpreter','tex');
    exportgraphics(fh, fullfile(figDir, sprintf('phaseConsistency_%s.png', sessID)), 'Resolution',120); close(fh);

    % --- burst length: peak vs ridge histograms ---
    bp = pp.peakBurstLength(~isnan(pp.peakBurstLength));
    br = pp.ridgeBurstLength(~isnan(pp.ridgeBurstLength));
    edges = 0:100:2600;
    fh = figure('Visible','off','Position',[100 100 640 460]); hold on;
    histogram(bp, edges, 'FaceColor',[0 0.3 0.8], 'FaceAlpha',0.5, 'EdgeColor','none');
    histogram(br, edges, 'FaceColor',[0.85 0.2 0], 'FaceAlpha',0.5, 'EdgeColor','none');
    xlabel('burst length (ms)'); ylabel('trials');
    legend({sprintf('peak (med %.0f)',median(bp)), sprintf('ridge (med %.0f)',median(br))});
    title(sprintf('%s  burst length: peak vs ridge', strrep(sessID,'_','\_')),'Interpreter','tex');
    exportgraphics(fh, fullfile(figDir, sprintf('burstLength_%s.png', sessID)), 'Resolution',120); close(fh);
end
