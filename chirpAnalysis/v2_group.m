function v2_group(aggDir, groupFigDir, outDir, C, task, cats)
% V2_GROUP  Group-level figures + R-stats CSVs, split by trial-type category (spec 4/5).
%   v2_group(aggDir, groupFigDir, outDir, C, task, cats)
%   Reads per-subject aggregates (each with agg.byCat over categories, agg.group = Dupi/OBE).
%   For every (category, group): group mean TFR, peak-aligned mean TFR, phase continuity
%   (subject light lines + group-mean thick), burst histogram -> groupFigDir, named
%   group<...>_<task>_<cat>_<group>.png. CSVs (per task) carry a category column for R.

    if ~isfolder(groupFigDir), mkdir(groupFigDir); end
    d = dir(fullfile(aggDir,'*_agg.mat'));
    if isempty(d), fprintf('v2_group(%s): no agg files\n', task); return; end
    A = [];
    for k = 1:numel(d), S = load(fullfile(d(k).folder,d(k).name),'agg'); A = [A; S.agg]; end %#ok<AGROW>
    freqs=A(1).freqs; tMs=A(1).tMs; relFreq=A(1).relFreq; relTimeMs=A(1).relTimeMs; phaseRelT=A(1).phaseRelT;

    % flatten to per (subject, category) records
    R = struct('group',{},'cat',{},'sessID',{},'meanTFR',{},'meanPeakTFR',{}, ...
        'peakConsMean',{},'ridgeConsMean',{},'medPeakBurst',{},'medRidgeBurst',{},'peakBurst',{},'ridgeBurst',{});
    for j = 1:numel(A)
        for c = 1:numel(A(j).byCat)
            bc = A(j).byCat(c);
            R(end+1) = struct('group',A(j).group,'cat',char(bc.cat),'sessID',A(j).sessID, ...
                'meanTFR',bc.meanTFR,'meanPeakTFR',bc.meanPeakTFR, ...
                'peakConsMean',bc.peakConsMean,'ridgeConsMean',bc.ridgeConsMean, ...
                'medPeakBurst',bc.medPeakBurst,'medRidgeBurst',bc.medRidgeBurst, ...
                'peakBurst',bc.peakBurst,'ridgeBurst',bc.ridgeBurst); %#ok<AGROW>
        end
    end
    groups = unique({R.group});

    for ci = 1:numel(cats)
        cat = cats{ci};
        for gi = 1:numel(groups)
            g = groups{gi};
            sel = find(strcmp({R.cat},cat) & strcmp({R.group},g)); ns = numel(sel);
            if ns < 1, continue; end
            tag = sprintf('%s_%s_%s', task, cat, g);

            TT = nan(numel(freqs),numel(tMs),ns);
            for j=1:ns, TT(:,:,j)=R(sel(j)).meanTFR; end
            v2_fig_meantfr(mean(TT,3,'omitnan'), freqs, tMs, ...
                fullfile(groupFigDir,sprintf('groupMeanTFR_%s.png',tag)), ...
                sprintf('GROUP %s  mean TFR (n=%d)', tag, ns), C);

            PP = nan(numel(relFreq),numel(relTimeMs),ns);
            for j=1:ns, PP(:,:,j)=R(sel(j)).meanPeakTFR; end
            v2_fig_peaktfr(mean(PP,3,'omitnan'), relFreq, relTimeMs, ...
                fullfile(groupFigDir,sprintf('groupPeakTFR_%s.png',tag)), ...
                sprintf('GROUP %s  peak-aligned mean TFR (n=%d)', tag, ns));

            PC=nan(ns,numel(phaseRelT)); RC=nan(ns,numel(phaseRelT));
            for j=1:ns, PC(j,:)=R(sel(j)).peakConsMean; RC(j,:)=R(sel(j)).ridgeConsMean; end
            fh=figure('Visible','off','Position',[100 100 780 540]); hold on;
            for j=1:ns, plot(phaseRelT,PC(j,:),'-','Color',[0 0.3 0.8 0.25]); end
            for j=1:ns, plot(phaseRelT,RC(j,:),'-','Color',[0.85 0.2 0 0.25]); end
            hP=plot(phaseRelT,mean(PC,1,'omitnan'),'-','Color',[0 0.2 0.7],'LineWidth',3);
            hR=plot(phaseRelT,mean(RC,1,'omitnan'),'-','Color',[0.75 0.1 0],'LineWidth',3);
            yline(C.pp.threshNarrow,'k:'); yline(C.pp.threshWide,'k--'); xlim([-1000 1000]); ylim([0 pi]);
            xlabel('time from gammaPeakTime (ms)'); ylabel('|phase diff| (rad)');
            legend([hP hR],{'peak (group mean)','ridge (group mean)'},'Location','north');
            title(sprintf('GROUP %s  phase continuity (n=%d)', tag, ns));
            exportgraphics(fh, fullfile(groupFigDir,sprintf('groupPhaseConsistency_%s.png',tag)),'Resolution',120); close(fh);

            allP=[]; allR=[];
            for j=1:ns, allP=[allP; R(sel(j)).peakBurst(:)]; allR=[allR; R(sel(j)).ridgeBurst(:)]; end %#ok<AGROW>
            v2_fig_bursthist(allP, allR, fullfile(groupFigDir,sprintf('groupBurst_%s.png',tag)), ...
                sprintf('GROUP %s  burst length peak vs ridge', tag));
        end
    end

    % --- CSVs for R (per task; carry category) ---
    nT = numel(phaseRelT); tot = numel(R)*nT;
    Gc=strings(tot,1); Cc=strings(tot,1); Sc=strings(tot,1); Tc=zeros(tot,1); Pc=zeros(tot,1); Rc=zeros(tot,1); ptr=0;
    for j=1:numel(R)
        r = ptr+(1:nT);
        Gc(r)=string(R(j).group); Cc(r)=string(R(j).cat); Sc(r)=string(R(j).sessID);
        Tc(r)=phaseRelT(:); Pc(r)=R(j).peakConsMean(:); Rc(r)=R(j).ridgeConsMean(:); ptr=ptr+nT;
    end
    writetable(table(Gc,Cc,Sc,Tc,Pc,Rc,'VariableNames',{'group','category','subject','timeMs','peakCons','ridgeCons'}), ...
        fullfile(outDir,sprintf('phaseConsistency_forR_%s.csv',task)));
    bu = table(string({R.group}'), string({R.cat}'), string({R.sessID}'), [R.medPeakBurst]', [R.medRidgeBurst]', ...
        'VariableNames',{'group','category','subject','medPeakBurst','medRidgeBurst'});
    writetable(bu, fullfile(outDir,sprintf('burst_forR_%s.csv',task)));
    fprintf('v2_group(%s): %d cats x %d groups; figs -> %s; CSVs -> %s\n', task, numel(cats), numel(groups), groupFigDir, outDir);
end
