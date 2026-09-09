function v2_fig_meantfr(M, freqs, tMs, figPath, titleStr, C)
% V2_FIG_MEANTFR  Mean TFR (finalOnset-locked), clim [-10 20]. Used per subject/category & group.
    d = fileparts(figPath); if ~isempty(d) && ~isfolder(d), mkdir(d); end
    fh = figure('Visible','off','Position',[100 100 820 470]);
    imagesc(tMs, freqs, M); axis xy; colormap(parula); cb = colorbar; cb.Label.String = 'mean z-power';
    clim([-10 20]); ylim(C.faslt.range); xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');
    title(titleStr, 'Interpreter','tex');
    exportgraphics(fh, figPath, 'Resolution', 120); close(fh);
end
