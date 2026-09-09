function chirp_fig_montage(tfr, R, ch, sessID, chi, figDir, C)
% CHIRP_FIG_MONTAGE  ~12 example single-trial TFRs + the trial-AVERAGED power TFR (spec 6.3).
%   Power-averaging preserves induced gamma (a time-domain average would destroy it).

    if ~isfolder(figDir), mkdir(figDir); end
    f = tfr.freqs; t = tfr.tMs;
    vi = find(tfr.valid); if isempty(vi), return; end
    pick = vi(round(linspace(1, numel(vi), min(12, numel(vi)))));

    fh = figure('Visible','off','Position',[60 60 1200 720]);
    tl = tiledlayout(fh, 4, 4, 'TileSpacing','compact', 'Padding','compact');
    for k = 1:numel(pick)
        nexttile; i = pick(k);
        imagesc(t, f, 10*log10(tfr.power(:,:,i)+eps)); axis xy; hold on;
        if ~isempty(R.trial(i).fhat), plot(t, R.trial(i).fhat, 'w-', 'LineWidth',1); end
        ylim(C.ridge.band); set(gca,'XTick',[]); title(sprintf('tr %d',i),'FontSize',8);
    end
    nexttile([1 2]);
    imagesc(t, f, 10*log10(tfr.meanPower+eps)); axis xy; hold on;
    yline(R.anchors.f_hi,'r:'); yline(R.anchors.f_lo,'r:');
    ylim(C.faslt.range); xlabel('ms'); ylabel('Hz'); title('trial-avg power TFR');
    title(tl, sprintf('%s  ch=%s  (n=%d valid)', strrep(sessID,'_','\_'), strrep(ch,'_','\_'), numel(vi)), ...
        'Interpreter','tex');

    fn = fullfile(figDir, sprintf('montage_%s_ch-%d.png', sessID, chi));
    exportgraphics(fh, fn, 'Resolution', 110);
    close(fh);
end
