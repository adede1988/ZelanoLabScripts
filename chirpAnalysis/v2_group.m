function v2_group(aggDir, groupFigDir, outDir, C)
% V2_GROUP  Group-level figures + R-stats CSVs (spec chirpAnalysis_3.md sec 4/5).
%   Reads per-subject aggregates (aggDir\*_agg.mat), groups by type (Dupi/OBE), and writes to
%   groupFigDir: group mean TFR, group peak-aligned mean TFR, group phase-continuity (subject
%   light lines + group-mean thick, peak vs ridge), group burst histogram. Writes to outDir:
%   phaseConsistency_forR.csv (group,subject,timeMs,peakCons,ridgeCons) and
%   burst_forR.csv (group,subject,medPeakBurst,medRidgeBurst) for the R permutation stats.

    if ~isfolder(groupFigDir), mkdir(groupFigDir); end
    d = dir(fullfile(aggDir,'*_agg.mat'));
    if isempty(d), fprintf('v2_group: no agg files in %s\n', aggDir); return; end
    A = [];
    for k = 1:numel(d), S = load(fullfile(d(k).folder,d(k).name),'agg'); A = [A; S.agg]; end %#ok<AGROW>
    groups = unique({A.group});
    freqs=A(1).freqs; tMs=A(1).tMs; relFreq=A(1).relFreq; relTimeMs=A(1).relTimeMs; phaseRelT=A(1).phaseRelT;

    for gi = 1:numel(groups)
        g = groups{gi}; sel = find(strcmp({A.group}, g)); ns = numel(sel);

        % --- group mean TFR (avg of subject means) ---
        TT = nan(numel(freqs), numel(tMs), ns);
        for j = 1:ns, TT(:,:,j) = A(sel(j)).meanTFR; end
        fh = figure('Visible','off','Position',[100 100 820 470]);
        imagesc(tMs, freqs, mean(TT,3,'omitnan')); axis xy; colormap(parula);
        cb = colorbar; cb.Label.String='mean z-power'; clim([-10 20]);
        ylim(C.faslt.range); xlabel('time from finalOnset (ms)'); ylabel('frequency (Hz)');
        title(sprintf('GROUP %s  mean TFR (n=%d subj)', g, ns));
        exportgraphics(fh, fullfile(groupFigDir,sprintf('groupMeanTFR_%s.png',g)),'Resolution',120); close(fh);

        % --- group peak-aligned mean TFR ---
        PPm = nan(numel(relFreq), numel(relTimeMs), ns);
        for j = 1:ns, PPm(:,:,j) = A(sel(j)).meanPeakTFR; end
        v2_fig_peaktfr(mean(PPm,3,'omitnan'), relFreq, relTimeMs, ...
            fullfile(groupFigDir,sprintf('groupPeakTFR_%s.png',g)), ...
            sprintf('GROUP %s  peak-aligned mean TFR (n=%d)', g, ns));

        % --- group phase continuity (subject light lines + group-mean thick; peak vs ridge) ---
        PC = nan(ns,numel(phaseRelT)); RC = nan(ns,numel(phaseRelT));
        for j = 1:ns, PC(j,:) = A(sel(j)).peakConsMean; RC(j,:) = A(sel(j)).ridgeConsMean; end
        fh = figure('Visible','off','Position',[100 100 780 540]); hold on;
        for j = 1:ns, plot(phaseRelT, PC(j,:), '-', 'Color',[0 0.3 0.8 0.25]); end
        for j = 1:ns, plot(phaseRelT, RC(j,:), '-', 'Color',[0.85 0.2 0 0.25]); end
        hP = plot(phaseRelT, mean(PC,1,'omitnan'), '-', 'Color',[0 0.2 0.7], 'LineWidth',3);
        hR = plot(phaseRelT, mean(RC,1,'omitnan'), '-', 'Color',[0.75 0.1 0], 'LineWidth',3);
        yline(C.pp.threshNarrow,'k:'); yline(C.pp.threshWide,'k--'); xlim([-1000 1000]); ylim([0 pi]);
        xlabel('time from gammaPeakTime (ms)'); ylabel('|phase diff| (rad)');
        legend([hP hR], {'peak (group mean)','ridge (group mean)'}, 'Location','north');
        title(sprintf('GROUP %s  phase continuity (n=%d subj)', g, ns));
        exportgraphics(fh, fullfile(groupFigDir,sprintf('groupPhaseConsistency_%s.png',g)),'Resolution',120); close(fh);

        % --- group burst histogram (all trials in group) ---
        allP=[]; allR=[];
        for j = 1:ns, allP=[allP; A(sel(j)).peakBurst(:)]; allR=[allR; A(sel(j)).ridgeBurst(:)]; end %#ok<AGROW>
        fh = figure('Visible','off','Position',[100 100 640 460]); hold on; edges = 0:100:2600;
        histogram(allP, edges, 'FaceColor',[0 0.3 0.8], 'FaceAlpha',0.5, 'EdgeColor','none');
        histogram(allR, edges, 'FaceColor',[0.85 0.2 0], 'FaceAlpha',0.5, 'EdgeColor','none');
        xlabel('burst length (ms)'); ylabel('trials');
        legend({sprintf('peak (med %.0f)',median(allP,'omitnan')), sprintf('ridge (med %.0f)',median(allR,'omitnan'))});
        title(sprintf('GROUP %s  burst length peak vs ridge', g));
        exportgraphics(fh, fullfile(groupFigDir,sprintf('groupBurst_%s.png',g)),'Resolution',120); close(fh);
    end

    % --- CSVs for R (across all subjects) ---
    nT = numel(phaseRelT); nSub = numel(A); tot = nSub*nT;
    Gc=strings(tot,1); Sc=strings(tot,1); Tc=zeros(tot,1); Pc=zeros(tot,1); Rc=zeros(tot,1); ptr=0;
    for j = 1:nSub
        r = ptr + (1:nT);
        Gc(r)=string(A(j).group); Sc(r)=string(A(j).sessID); Tc(r)=phaseRelT(:);
        Pc(r)=A(j).peakConsMean(:); Rc(r)=A(j).ridgeConsMean(:); ptr = ptr + nT;
    end
    writetable(table(Gc,Sc,Tc,Pc,Rc,'VariableNames',{'group','subject','timeMs','peakCons','ridgeCons'}), ...
        fullfile(outDir,'phaseConsistency_forR.csv'));
    bu = table(string({A.group}'), string({A.sessID}'), [A.medPeakBurst]', [A.medRidgeBurst]', ...
        'VariableNames',{'group','subject','medPeakBurst','medRidgeBurst'});
    writetable(bu, fullfile(outDir,'burst_forR.csv'));
    fprintf('v2_group: %d groups (%s); group figs -> %s; CSVs -> %s\n', ...
        numel(groups), strjoin(groups,','), groupFigDir, outDir);
end
