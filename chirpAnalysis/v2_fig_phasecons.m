function v2_fig_phasecons(Mpk, Mrg, gc, relT, figPath, titleStr, C)
% V2_FIG_PHASECONS  Peak-centered phase continuity on ONE set of axes: peak (blue) vs ridge (red)
%   |phase diff| per trial (light) + across-trial means (thick). gc = trial indices to include.
    d = fileparts(figPath); if ~isempty(d) && ~isfolder(d), mkdir(d); end
    gc = gc(:)'; if isempty(gc), return; end
    fh = figure('Visible','off','Position',[100 100 780 540]); hold on;
    for i = gc, plot(relT, abs(Mpk(i,:)), '-', 'Color',[0 0.3 0.8 0.10]); end
    for i = gc, plot(relT, abs(Mrg(i,:)), '-', 'Color',[0.85 0.2 0 0.10]); end
    hP = plot(relT, mean(abs(Mpk(gc,:)),1,'omitnan'), '-', 'Color',[0 0.2 0.7], 'LineWidth',2.6);
    hR = plot(relT, mean(abs(Mrg(gc,:)),1,'omitnan'), '-', 'Color',[0.75 0.1 0], 'LineWidth',2.6);
    yline(C.pp.threshNarrow,'k:'); yline(C.pp.threshWide,'k--'); xlim([-1000 1000]); ylim([0 pi]);
    xlabel('time from gammaPeakTime (ms)'); ylabel('|phase diff| (rad)');
    legend([hP hR], {'peakPhaseConsistency (mean)','ridgePhaseConsistency (mean)'}, 'Location','north');
    title(titleStr, 'Interpreter','tex');
    exportgraphics(fh, figPath, 'Resolution', 120); close(fh);
end
