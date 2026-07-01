function v2_fig_peaktfr(M, relFreq, relTimeMs, figPath, titleStr)
% V2_FIG_PEAKTFR  Save a peak-aligned mean TFR (freq centered on gammaPeakFrequency, time on
%   gammaPeakTime). M [nRelF x nPk], relFreq (Hz rel peak), relTimeMs (ms rel peak). clim [-10 20].
    d = fileparts(figPath); if ~isempty(d) && ~isfolder(d), mkdir(d); end
    fh = figure('Visible','off','Position',[100 100 720 470]);
    imagesc(relTimeMs, relFreq, M); axis xy; colormap(parula);
    cb = colorbar; cb.Label.String = 'mean z-power'; clim([-10 20]);
    xlabel('time from gammaPeakTime (ms)'); ylabel('frequency - gammaPeakFrequency (Hz)');
    yline(0,'w:'); xline(0,'w:');
    title(titleStr, 'Interpreter','tex');
    exportgraphics(fh, figPath, 'Resolution',120); close(fh);
end
