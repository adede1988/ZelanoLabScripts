function v2_fig_bursthist(peakBL, ridgeBL, figPath, titleStr)
% V2_FIG_BURSTHIST  Overlaid peak (blue) vs ridge (red) burst-length histograms.
    d = fileparts(figPath); if ~isempty(d) && ~isfolder(d), mkdir(d); end
    bp = peakBL(~isnan(peakBL)); br = ridgeBL(~isnan(ridgeBL)); edges = 0:100:2600;
    fh = figure('Visible','off','Position',[100 100 640 460]); hold on;
    if ~isempty(bp), histogram(bp, edges, 'FaceColor',[0 0.3 0.8], 'FaceAlpha',0.5, 'EdgeColor','none'); end
    if ~isempty(br), histogram(br, edges, 'FaceColor',[0.85 0.2 0], 'FaceAlpha',0.5, 'EdgeColor','none'); end
    xlabel('burst length (ms)'); ylabel('trials');
    legend({sprintf('peak (med %.0f)', medSafe(bp)), sprintf('ridge (med %.0f)', medSafe(br))});
    title(titleStr, 'Interpreter','tex');
    exportgraphics(fh, figPath, 'Resolution', 120); close(fh);
end
function v = medSafe(x), if isempty(x), v=NaN; else, v=median(x); end, end
