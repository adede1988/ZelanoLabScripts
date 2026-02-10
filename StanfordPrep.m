datPre  = 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\';

sessionIDs = {'250818_Dupi_NMH_JH_1', ... 
           '250623_DUPI_NMH_KS_2',...
           '250623_Dupi_NMH_KS_1',...
           '250818_Dupi_NMH_JH_2',...
           '250811_Dupi_NMH_TPB_1',...
           '250811_Dupi_NMH_TB_2',...
           '250929_Dupi_NMH_GH_1',...
           '251002_Dupi_NMH_AB_1',...
           '251027_Dupi_NMH_DL_1', ...
           '250929_Dupi_NMH_GH_2',...
           '251002_Dupi_NMH_AB_2',...
           '251013_Dupi_NMH_JN_2',...
            '251030_Dupi_NMH_DB_1',...
            '251110_Dupi_NMH_PC_1'};

figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';

taskNames = {'breathingPreProc'};

% -----------------------------
% Parameters
% -----------------------------
plotWinMS   = [-1000 3000];
statWinMS   = [-250 1000];

summaryVars = table;

% -----------------------------
% Loop subjects
% -----------------------------
for s = 1:numel(sessionIDs)

    for t = 1:numel(taskNames)

        % (avoid variable carryover)
        clear bestChanIdx bestChanLabel macroIdx macroLabels labelsNorm psd fooofOut

        % -----------------------------
        % Load data
        % -----------------------------
        fname = fullfile(datPre, sessionIDs{s}, 'preProc', ...
            [sessionIDs{s} '_' taskNames{t} '.mat']);

        if ~exist(fname, 'file')
            warning('Missing file: %s', fname);
            continue
        end

        tmp = load(fname);
        try
            outDat = tmp.out;
        catch
            outDat = tmp.outDat;
        end

        % ensure output figs folder is defined
        if ~isfield(outDat, 'figs') || isempty(outDat.figs)
            outDat.figs = figPath;
        end

        fs = outDat.fs;
        labels = outDat.labels;
        data   = outDat.data;
        behDat = outDat.behDat;

        inhaleIdx = behDat.finalOnset;
        inhaleIdx = inhaleIdx(~isnan(inhaleIdx));

        % -----------------------------
        % Identify channels (robust labels)
        % -----------------------------
        labelsNorm = strings(size(labels));
        for iL = 1:numel(labels)
            x = labels{iL};

            if ischar(x)
                labelsNorm(iL) = string(x);

            elseif isstring(x)
                if ~isempty(x)
                    labelsNorm(iL) = x(1);
                else
                    labelsNorm(iL) = "";
                end

            elseif iscell(x)
                if ~isempty(x)
                    y = x{1};
                    if ischar(y) || isstring(y)
                        labelsNorm(iL) = string(y);
                    else
                        labelsNorm(iL) = "";
                    end
                else
                    labelsNorm(iL) = "";
                end
            else
                labelsNorm(iL) = "";
            end
        end

        badIdx   = find(labelsNorm == "badTS");
        blinkIdx = find(labelsNorm == "blinkIndicator");
        spikeIdx = find(labelsNorm == "spikeCleanVec");

        % -----------------------------
        % Build time axes
        % -----------------------------
        plotWinSamp = round(plotWinMS/1000 * fs);
        statWinSamp = round(statWinMS/1000 * fs); 

        tPlot = (plotWinSamp(1):plotWinSamp(2)) / fs * 1000;

        % -----------------------------
        % Extract inhale-locked snippets (QC vectors)
        % -----------------------------
        nEvt = numel(inhaleIdx);
        nT   = numel(plotWinSamp(1):plotWinSamp(2));

        badMat   = nan(nEvt, nT);
        blinkMat = nan(nEvt, nT);
        spikeMat = nan(nEvt, nT);

        for e = 1:nEvt
            idx0 = inhaleIdx(e);
            idx  = idx0 + (plotWinSamp(1):plotWinSamp(2));

            if min(idx) < 1 || max(idx) > size(data,2)
                continue
            end

            if ~isempty(badIdx)
                badMat(e,:) = data(badIdx(1), idx);
            end
            if ~isempty(blinkIdx)
                blinkMat(e,:) = data(blinkIdx(1), idx);
            end
            if ~isempty(spikeIdx)
                spikeMat(e,:) = data(spikeIdx(1), idx);
            end
        end

        % -----------------------------
        % Flag bad trials (as requested) + compute badTrialPrct
        % -----------------------------
        if ~isempty(badIdx)
            badTrials = sum(badMat, 2, 'omitnan') > (size(badMat,2)/5);
            badTrialPrct = sum(badTrials) / numel(badTrials);
        else
            % can't compute; keep all trials for downstream analyses
            badTrials = false(nEvt,1);
            badTrialPrct = NaN;
        end
        goodTrials = ~badTrials;

        % -----------------------------
        % Average time courses (ONLY GOOD TRIALS)
        % -----------------------------
        if any(goodTrials)
            badMean   = nanmean(badMat(goodTrials,:),   1);
            blinkMean = nanmean(blinkMat(goodTrials,:), 1);
            spikeMean = nanmean(spikeMat(goodTrials,:), 1);
        else
            badMean   = nan(1, nT);
            blinkMean = nan(1, nT);
            spikeMean = nan(1, nT);
        end

        % -----------------------------
        % Summary statistics window (ONLY GOOD TRIALS)
        % -----------------------------
        statIdx = find(tPlot >= statWinMS(1) & tPlot <= statWinMS(2));

        pctBad   = mean(badMean(statIdx),   'omitnan') * 100;
        pctBlink = mean(blinkMean(statIdx), 'omitnan') * 100;
        pctSpike = mean(1 - spikeMean(statIdx), 'omitnan') * 100; % spikeCleanVec: 0 = bad

        % -----------------------------
        % Store summary row (robust to evolving columns)
        % -----------------------------
        rowIdx = height(summaryVars) + 1;

        coreVars = {'sessionID','task','pctBadTS','pctBlink','pctSpikeNoise','badTrialPrct'};
        for iV = 1:numel(coreVars)
            v = coreVars{iV};
            if ~ismember(v, summaryVars.Properties.VariableNames)
                switch v
                    case {'sessionID','task'}
                        summaryVars.(v) = strings(height(summaryVars),1);
                    otherwise
                        summaryVars.(v) = nan(height(summaryVars),1);
                end
            end
        end

        summaryVars.sessionID(rowIdx)     = string(sessionIDs{s});
        summaryVars.task(rowIdx)          = string(taskNames{t});
        summaryVars.pctBadTS(rowIdx)      = pctBad;
        summaryVars.pctBlink(rowIdx)      = pctBlink;
        summaryVars.pctSpikeNoise(rowIdx) = pctSpike;
        summaryVars.badTrialPrct(rowIdx)  = badTrialPrct;

        % ==============================
        % Participant/session parsing + normalization
        % ==============================
        if ~ismember('session_num', summaryVars.Properties.VariableNames)
            summaryVars.session_num = nan(height(summaryVars),1);
        end
        if ~ismember('participant_ID', summaryVars.Properties.VariableNames)
            summaryVars.participant_ID = strings(height(summaryVars),1);
        end

        sid = string(sessionIDs{s});
        tok = regexp(sid, '^(.*)_([123])$', 'tokens', 'once');
        if ~isempty(tok)
            pid = string(tok{1});
            sn  = str2double(tok{2});
        else
            pid = sid;
            sn  = NaN;
        end

        pid = regexprep(pid, 'DUPI', 'Dupi');
        pid = regexprep(pid, 'dupi', 'Dupi');

        if pid == "250811_Dupi_TPB"
            pid = "250811_Dupi_NMH_TB";
        end

        summaryVars.session_num(rowIdx)    = sn;
        summaryVars.participant_ID(rowIdx) = pid;

        % -----------------------------
        % Plot QC figure (ONLY GOOD TRIALS)
        % -----------------------------
        figure('Color','w','Position',[200 200 800 500]); hold on
        plot(tPlot, badMean,      'r', 'LineWidth',1.5)
        plot(tPlot, blinkMean,    'b', 'LineWidth',1.5)
        plot(tPlot, 1-spikeMean,  'k', 'LineWidth',1.5)

        xline(0,'k--');
        xlabel('Time from inhale onset (ms)')
        ylabel('Probability')
        legend({'badTS','blink','spike noise'}, 'Location','NorthWest')
        title(sprintf('%s | %s (good trials only)', sessionIDs{s}, taskNames{t}), 'interpreter', 'none')
        ylim([0 1])

        figTitle = 'QC_sniffLocked';
        saveas(gcf, fullfile(outDat.figs, [figTitle '.jpg']));
        close

        % ==============================
        % Spectral analysis (intranasal macro channels) - GOOD TRIALS ONLY
        % ==============================
        specWins = struct( ...
            'base',  [-500 0], ...
            'early', [0 500], ...
            'late',  [500 1500]);

        gammaBand = [30 80];

        isMacro = startsWith(labelsNorm, "macBP");
        macroIdx = find(isMacro);
        macroLabels = cellstr(labelsNorm(macroIdx));

        goodInhaleIdx = inhaleIdx(goodTrials);

        if isempty(macroIdx)
            warning('No macBP channels found for %s %s', sessionIDs{s}, taskNames{t});
        elseif isempty(goodInhaleIdx)
            warning('No GOOD trials available for spectra in %s %s', sessionIDs{s}, taskNames{t});
        else

            nfft = 2^nextpow2(round(1.5 * fs));
            f = (0:nfft/2)' * fs / nfft;
            keepF = f >= 2 & f <= 150;

            psd = struct();
            fooofOut = struct();

            winNames = fieldnames(specWins);

            for w = 1:numel(winNames)
                wname = winNames{w};
                wms   = specWins.(wname);
                wsamp = round(wms/1000 * fs);

                winVec = (wsamp(1):wsamp(2)-1);
                winLen = numel(winVec);

                psd.(wname) = nan(numel(macroIdx), nnz(keepF));

                for c = 1:numel(macroIdx)

                    segs = nan(numel(goodInhaleIdx), winLen);
                    nGoodSeg = 0;

                    for e = 1:numel(goodInhaleIdx)
                        idx0 = goodInhaleIdx(e);
                        idx  = idx0 + winVec;

                        if min(idx) < 1 || max(idx) > size(data,2)
                            continue
                        end

                        x = data(macroIdx(c), idx);
                        x = detrend(x);
                        x = x(:)' .* hann(numel(x))';   % row

                        nGoodSeg = nGoodSeg + 1;
                        segs(nGoodSeg,:) = x;
                    end

                    if nGoodSeg == 0
                        continue
                    end

                    segs = segs(1:nGoodSeg,:);

                    X = fft(segs, nfft, 2);
                    P = mean(abs(X(:,1:nfft/2+1)).^2, 1, 'omitnan');

                    psd.(wname)(c,:) = P(keepF);
                end
            end

            % choose channel with strongest gamma peak (baseline)
            gammaF = (f(keepF) >= gammaBand(1)) & (f(keepF) <= gammaBand(2));
            baseBlock = psd.base(:,gammaF);

            if all(isnan(baseBlock(:)))
                warning('All baseline spectra are NaN for %s %s (good trials). Skipping FOOOF.', sessionIDs{s}, taskNames{t});
            else
                [~, bestChanIdx] = max(max(baseBlock, [], 2, 'omitnan'));
                bestChanLabel = macroLabels{bestChanIdx};

                % plot spectra
                figure('Color','w','Position',[100 100 1100 400])
                colors = lines(numel(macroIdx));

                for w = 1:numel(winNames)
                    wname = winNames{w};
                    subplot(1,3,w); hold on

                    for c = 1:numel(macroIdx)
                        if c == bestChanIdx
                            lw = 2.5; col = [0 0 0];
                        else
                            lw = 1;   col = colors(c,:);
                        end
                        plot(f(keepF), log10(psd.(wname)(c,:)), 'Color', col, 'LineWidth', lw)
                    end

                    title(wname)
                    xlabel('Frequency (Hz)')
                    ylabel('log10 Power')
                    xlim([2 150])

                    % fooof on selected channel
                    fooofOut.(wname) = fooof_basic( ...
                        f(keepF), psd.(wname)(bestChanIdx,:)', ...
                        'f_range', [2 150], ...
                        'aperiodic_mode', 'knee', ...
                        'max_peaks', 6);
                end

                sgtitle(sprintf('%s | %s | macro=%s (good trials only)', ...
                    sessionIDs{s}, taskNames{t}, bestChanLabel), 'interpreter', 'none')

                figTitle = 'macroSpectra';
                saveas(gcf, fullfile(outDat.figs, [figTitle '.jpg']));
                close

                % append FOOOF metrics
                summaryVars.bestMacroChannel(rowIdx) = string(bestChanLabel);

                for w = 1:numel(winNames)
                    wname = winNames{w};
                    F = fooofOut.(wname);

                    summaryVars.([wname '_offset'])(rowIdx)   = F.ap.offset;
                    summaryVars.([wname '_exponent'])(rowIdx) = F.ap.exponent;

                    nPeaks = size(F.peaks,1);
                    for p = 1:nPeaks
                        summaryVars.([wname '_peakFreq' num2str(p)])(rowIdx) = F.peaks(p,1);
                        summaryVars.([wname '_peakAmp'  num2str(p)])(rowIdx) = F.peaks(p,2);
                    end
                end
            end
        end

        % ==============================
        % Time-frequency (sniff-locked) plot over plotWinMS - GOOD TRIALS ONLY
        % ==============================
        if exist('macroIdx','var') && ~isempty(macroIdx) && ~isempty(goodInhaleIdx)

            if exist('bestChanIdx','var') && ~isempty(bestChanIdx) && bestChanIdx <= numel(macroIdx)
                tfChanIdx = macroIdx(bestChanIdx);
            else
                tfChanIdx = macroIdx(1);
            end

          % --- Time-frequency params ---
            tf_f_range = [2 150];
            tf_win_sec = 0.25;
            tf_ovlp    = 0.80;
            
            tf_win_samp = max(32, round(tf_win_sec * fs));
            tf_win_samp = tf_win_samp + mod(tf_win_samp,2);
            tf_noverlap = round(tf_ovlp * tf_win_samp);
            tf_nfft     = 2^nextpow2(tf_win_samp);
            
            % --- extend TF extraction by +/- 2 sec, then trim to plotWinMS for plotting ---
            padSec  = 2;
            padSamp = round(padSec * fs);
            
            plotWinVec   = (plotWinSamp(1):plotWinSamp(2));          % target plotting window (samples)
            extWinSamp   = [plotWinSamp(1)-padSamp, plotWinSamp(2)+padSamp];
            extWinVec    = (extWinSamp(1):extWinSamp(2));            % extended window for TF computation
            
            Psum   = [];
            nGoodTF = 0;
            Fkeep  = [];
            TmsExt = [];
            
            win = hann(tf_win_samp, 'periodic');  % STFT window
            allSnips = []; 
            for e = 1:numel(goodInhaleIdx)
                idx0 = goodInhaleIdx(e);
                idx  = idx0 + extWinVec;
            
                if min(idx) < 1 || max(idx) > size(data,2)
                    continue
                end
            
                x = data(tfChanIdx, idx);
                x = detrend(x(:));  % no global taper
                allSnips = [allSnips, x]; 
                % Use PSD output (scaled) from spectrogram
                [~,F,T,P] = spectrogram(x, win, tf_noverlap, tf_nfft, fs);  % P: freq x time
            
                if isempty(Fkeep)
                    Fkeep  = (F >= tf_f_range(1) & F <= tf_f_range(2));
                    % Time axis in ms relative to sniff onset:
                    % ext window starts at (plotWinMS(1) - padSec*1000)
                    extStartMS = plotWinMS(1) - padSec*1000;
                    TmsExt = (T * 1000) + extStartMS;  % ms
            
                    Psum = zeros(nnz(Fkeep), numel(T), 'double');
                end
            
                Psum = Psum + P(Fkeep,:);
                nGoodTF = nGoodTF + 1;
            end
            
            % ---- trim to plotWinMS for plotting ----
            if nGoodTF > 0
                PmeanExt = Psum ./ nGoodTF;
                Fplot = F(Fkeep);
            
                trimMask = (TmsExt >= plotWinMS(1)) & (TmsExt <= plotWinMS(2));
                Tms = TmsExt(trimMask);
                Pmean = PmeanExt(:, trimMask);
            
                % Optional baseline-normalize in dB using -500..0ms (within trimmed window)
                doBaselineDB = true;
                if doBaselineDB
                    baseMask = (Tms >= -500 & Tms <= 0);
                    if any(baseMask)
                        basePow = mean(Pmean(:, baseMask), 2);
                        Pplot = 10*log10(Pmean ./ max(basePow, eps));
                        cbarLabel = 'Power (dB vs -500..0 ms)';
                    else
                        Pplot = 10*log10(Pmean);
                        cbarLabel = 'Power (10*log10)';
                    end
                else
                    Pplot = 10*log10(Pmean);
                    cbarLabel = 'Power (10*log10)';
                end
            
                % ---- plot ----
                figure('Color','w','Position',[150 150 900 500]);
                imagesc(Tms, Fplot, Pplot);
                axis xy
                xline(0,'k--','LineWidth',1);
                xlabel('Time from inhale onset (ms)');
                ylabel('Frequency (Hz)');
                title(sprintf('%s | %s | TF: %s', sessionIDs{s}, taskNames{t}, labelsNorm(tfChanIdx)), 'interpreter','none');
                cb = colorbar; ylabel(cb, cbarLabel);
            
                % ---- save ----
                figTitle = 'sniff_TF';
                saveas(gcf, fullfile(outDat.figs, [figTitle '.jpg']));
                close
            end
        end

        % ==============================
        % CueTask performance: d'  (UNCHANGED)
        % ==============================
        if strcmp(taskNames{t}, 'cueTaskPreProc') && istable(behDat) && any(strcmp(behDat.Properties.VariableNames,'type'))

            if ~ismember('performance', summaryVars.Properties.VariableNames)
                summaryVars.performance = nan(height(summaryVars),1);
            end
            if ~ismember('performanceType', summaryVars.Properties.VariableNames)
                summaryVars.performanceType = strings(height(summaryVars),1);
            end

            tt = string(behDat.type);
            tt = tt(~ismissing(tt));

            nHit  = sum(tt == "hit");
            nMiss = sum(tt == "miss");
            nFA   = sum(tt == "fa");
            nCR   = sum(tt == "cr");

            nSignal = nHit + nMiss;
            nNoise  = nFA  + nCR;

            dprime = NaN;
            if nSignal > 0 && nNoise > 0
                Hc  = (nHit + 0.5) / (nSignal + 1);
                FAc = (nFA  + 0.5) / (nNoise  + 1);
                dprime = norminv(Hc) - norminv(FAc);
            end

            summaryVars.performance(rowIdx)     = dprime;
            summaryVars.performanceType(rowIdx) = "d_prime";
        end

        % ==============================
        % O15 performance: mean expScore across 15 targets (UNCHANGED)
        % ==============================
        if strcmp(taskNames{t}, 'O15preproc') && istable(behDat) && ...
                all(ismember({'target','expScore'}, behDat.Properties.VariableNames))

            if ~ismember('performance', summaryVars.Properties.VariableNames)
                summaryVars.performance = nan(height(summaryVars),1);
            end
            if ~ismember('performanceType', summaryVars.Properties.VariableNames)
                summaryVars.performanceType = strings(height(summaryVars),1);
            end

            T = behDat(:, {'target','expScore'});
            T = T(~isnan(T.expScore), :);

            if ~iscategorical(T.target)
                T.target = string(T.target);
            end

            G = findgroups(T.target);
            expPerTarget = splitapply(@(x) x(1), T.expScore, G);

            acc = mean(expPerTarget, 'omitnan');

            summaryVars.performance(rowIdx)     = acc;
            summaryVars.performanceType(rowIdx) = "accuracy";
        end

    end
end

% -----------------------------
% Save summary table
% -----------------------------
writetable(summaryVars, fullfile(figPath, 'QC_summary_metrics.csv'));
