function v2_figs_ridge(T, ridgeInfo, figDir, sessID, chi, C)
% V2_FIGS_RIDGE  Single-trial TF (z) + primary/secondary ridge overlays, and both ridge
%   ensembles over the mean TF (spec 4). Saves PNGs under figDir. chi = channel label (string).

    if ~isfolder(figDir), mkdir(figDir); end
    freqs = T.freqs; tMs = T.tMs; good = find(T.good);
    pf = ridgeInfo.primaryRidge.f; sf = ridgeInfo.secondaryRidge.f;
    chStr = char(string(chi));
    if isempty(good), return; end
    avg = mean(T.zTFR(:,:,good), 3, 'omitnan');

    % --- example single-trial TF + both ridges ---
    pick = unique(round(linspace(1, numel(good), min(C.nFigTrials, numel(good)))));
    for j = pick
        i = good(j);
        fh = figure('Visible','off','Position',[100 100 800 470]);
        imagesc(tMs, freqs, T.zTFR(:,:,i)); axis xy; hold on;
        colormap(parula); cb = colorbar; cb.Label.String = 'z-power'; clim([-10 20]);
        ylim(C.faslt.range); xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');
        plot(tMs, pf(i,:), '-',  'Color',[1 1 1],      'LineWidth',1.8);
        plot(tMs, sf(i,:), '--', 'Color',[1 0 1],      'LineWidth',1.2);
        yline(C.ridge.band(1),'w:'); yline(C.ridge.band(2),'w:');
        legend({'primary','secondary'}, 'TextColor','w','Color','k','Location','northeast');
        title(sprintf('%s  ch=%s  trial %d', strrep(sessID,'_','\_'), strrep(chStr,'_','\_'), i),'Interpreter','tex');
        fn = fullfile(figDir, sprintf('sub-%s_ch-%s_trial-%03d.png', sessID, chStr, i));
        exportgraphics(fh, fn, 'Resolution', 110); close(fh);
    end

    % --- subject mean TFR (clim [-10 20]) ---
    fh = figure('Visible','off','Position',[100 100 820 470]);
    imagesc(tMs, freqs, avg); axis xy; colormap(parula); cb = colorbar; cb.Label.String = 'mean z-power';
    clim([-10 20]); ylim(C.faslt.range); xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');
    title(sprintf('%s  mean TFR (n=%d)', strrep(sessID,'_','\_'), numel(good)),'Interpreter','tex');
    exportgraphics(fh, fullfile(figDir, sprintf('meanTFR_%s.png', sessID)), 'Resolution',120); close(fh);

    % --- all primary ridges over mean TF ---
    ensembleFig(tMs, freqs, avg, pf, good, [1 1 1 0.15], 'r', C, ...
        sprintf('%s  primary ridges (n=%d) over mean TF', strrep(sessID,'_','\_'), numel(good)), ...
        fullfile(figDir, sprintf('primaryRidges_%s.png', sessID)));

    % --- all secondary ridges over mean TF ---
    ensembleFig(tMs, freqs, avg, sf, good, [1 1 1 0.15], 'm', C, ...
        sprintf('%s  secondary ridges (n=%d) over mean TF', strrep(sessID,'_','\_'), numel(good)), ...
        fullfile(figDir, sprintf('secondaryRidges_%s.png', sessID)));
end

function ensembleFig(tMs, freqs, avg, rf, good, lineCol, meanCol, C, ttl, fn)
    fh = figure('Visible','off','Position',[100 100 840 480]);
    imagesc(tMs, freqs, avg); axis xy; hold on;
    colormap(parula); cb = colorbar; cb.Label.String = 'mean z-power'; clim([-10 20]);
    ylim(C.faslt.range); xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');
    for i = good', plot(tMs, rf(i,:), '-', 'Color', lineCol, 'LineWidth', 0.8); end
    plot(tMs, mean(rf(good,:), 1, 'omitnan'), '-', 'Color', meanCol, 'LineWidth', 2);
    yline(C.ridge.band(1),'w:'); yline(C.ridge.band(2),'w:');
    title(ttl, 'Interpreter','tex');
    exportgraphics(fh, fn, 'Resolution', 120); close(fh);
end
