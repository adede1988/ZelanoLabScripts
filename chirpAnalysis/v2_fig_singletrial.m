function v2_fig_singletrial(Z, freqs, tMs, pfi, sfi, figPath, titleStr, C)
% V2_FIG_SINGLETRIAL  One trial's z-TFR with primary (white) + secondary (magenta) ridge. clim [-10 20].
    d = fileparts(figPath); if ~isempty(d) && ~isfolder(d), mkdir(d); end
    if all(~isfinite(Z(:))), return; end
    fh = figure('Visible','off','Position',[100 100 800 470]);
    imagesc(tMs, freqs, Z); axis xy; hold on; colormap(parula); cb = colorbar; cb.Label.String = 'z-power';
    clim([-10 20]); ylim(C.faslt.range); xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');
    plot(tMs, pfi, '-',  'Color',[1 1 1], 'LineWidth',1.8);
    plot(tMs, sfi, '--', 'Color',[1 0 1], 'LineWidth',1.2);
    yline(C.ridge.band(1),'w:'); yline(C.ridge.band(2),'w:');
    legend({'primary','secondary'}, 'TextColor','w','Color','k','Location','northeast');
    title(titleStr, 'Interpreter','tex');
    exportgraphics(fh, figPath, 'Resolution', 110); close(fh);
end
