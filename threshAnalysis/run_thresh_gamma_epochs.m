function run_thresh_gamma_epochs(sessFilter)
% RUN_THRESH_GAMMA_EPOCHS  Time-resolved gamma via FOOOF, split by ODOR CONDITION.
%   Per session (bestMac, finalOnset-locked, noise-clean trials) and per condition
%   c in {low, med, air}: Morlet power 2-58 Hz (100 linear freqs), averaged over
%   that condition's trials & 250 ms windows into 10 sequential power spectra
%   (-500..+2000 ms); FOOOF each; the periodic spectra are plotted (ochre->purple)
%   to <figs>/gammaTimeProgression_<cond>.png and the largest 30-58 Hz peak per
%   (condition x epoch) is written to threshTask_gammaEpochs.csv (with a `cond`
%   column).
%
%   Low-N guard: a condition with fewer than MINN clean trials is skipped and
%   logged (no noisy FOOOF is emitted).

    if nargin < 1, sessFilter = []; end
    thresh_init_paths(); L = labPaths();
    groupDir = getenv('THRESH_GROUPDIR');                    % per-job override (parallel runs)
    if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end
    epoCsv   = fullfile(groupDir, 'threshTask_gammaEpochs.csv');

    MINN    = 8;                            % min clean trials/condition to run a progression
    epWin   = [-1.75 5.75];                 % s, finalOnset-locked
    cycles  = [3 0.8];
    winEdges = -500:250:2000;               % ms -> 10 windows
    nWin    = numel(winEdges) - 1;
    gammaBand = [30 58];
    conds = {'low','med','air'};
    opt = build_fooof_opt([2 58]);          % fit 2-58 Hz with a knee (matches thresh_fooof_macBP)
    hasOpt = exist('fmincon','file') > 0;

    T = thresh_session_table(false); T = T(T.onDisk, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    hdr = {'sessID','subID','sessNum','type','group','cond','epoch','centerMs', ...
           'gammaPeakFreq','gammaPeakPower','gammaDetected'};
    newRows = {};
    doneIDs = strings(0,1);

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n== %d/%d %s ==\n', i, height(T), id);
        try
            Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
            if ~isfield(od,'bestMac') || isempty(od.bestMac), fprintf('  no bestMac -> skip\n'); clear od; continue; end
            labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
            ci = find(strcmp(od.bestMac, labs), 1); if isempty(ci), clear od; continue; end
            sig = double(od.data(ci, :)); fs = od.fs;

            % noise rejection on TTL.start; bucket clean finalOnset by odor condition
            ts = od.TTL.start;
            NT = thresh_noise_trials(sig, fs, ts);
            bd = od.behDat; typ = string(bd.type);
            foByCond = struct('low',[],'med',[],'air',[]);
            for j = 1:height(bd)
                k = bd.n(j);
                if ~(k >= 1 && k <= numel(ts) && k <= numel(NT.noisy)), continue; end
                cond = char(typ(j));
                if ~ismember(cond, conds), continue; end
                if ~NT.noisy(k) && isfinite(bd.finalOnset(j))
                    foByCond.(cond)(end+1,1) = bd.finalOnset(j);
                end
            end

            figDir = fullfile(L.figPath, id, 'threshTask'); if ~isfolder(figDir), mkdir(figDir); end
            anyDone = false;

            for cc = 1:numel(conds)
                cond = conds{cc}; foV = foByCond.(cond);
                if numel(foV) < MINN
                    fprintf('  %s: %d clean trials (< %d) -> skipped\n', cond, numel(foV), MINN);
                    continue;
                end

                % epoch + Morlet power (per trial) via newtimef
                s0 = round(epWin(1)*fs); s1 = round(epWin(2)*fs); nF = s1 - s0 + 1; Tn = numel(sig);
                ep = zeros(nF, numel(foV)); kept = 0;
                for e = foV(:)'
                    a = round(e)+s0; b = round(e)+s1;
                    if a < 1 || b > Tn, continue; end
                    seg = sig(a:b); if any(~isfinite(seg)), continue; end
                    kept = kept + 1; ep(:,kept) = seg(:);
                end
                ep = ep(:, 1:kept);
                if kept < MINN, fprintf('  %s: %d in-bounds (< %d) -> skipped\n', cond, kept, MINN); continue; end

                [~,~,~,times,fout,~,~,atf] = newtimef(ep, nF, epWin*1000, fs, cycles, ...
                    'freqs', [2 58], 'nfreqs', 100, 'freqscale', 'linear', ...
                    'baseline', NaN, 'plotersp','off','plotitc','off','verbose','off');
                P = abs(atf).^2;                 % [F x time x trial]
                meanP = mean(P, 3, 'omitnan');   % [F x time]
                fU = fout(:)';                   % newtimef freq grid

                % 10 windowed power spectra
                spec = zeros(nWin, numel(fU)); centerMs = zeros(nWin,1);
                for w = 1:nWin
                    tmask = times >= winEdges(w) & times < winEdges(w+1);
                    spec(w,:) = mean(meanP(:, tmask), 2, 'omitnan')';
                    centerMs(w) = mean(winEdges(w:w+1));
                end

                % FOOOF each epoch spectrum. process_fooof expects TF [nSig x 1 x nFreq].
                TF = zeros(nWin, 1, numel(fU)); TF(:,1,:) = spec;
                [~, fg] = process_fooof('FOOOF_matlab', TF, fU, opt, hasOpt);
                fMask = fU >= opt.freq_range(1) & fU <= opt.freq_range(2);
                fAxis = fU(fMask);
                specM = spec(:, fMask);

                gF = nan(nWin,1); gP = nan(nWin,1); det = zeros(nWin,1);
                apfit = zeros(nWin, numel(fAxis));
                for w = 1:nWin
                    apfit(w,:) = fg(w).ap_fit(:)';
                    pp = fg(w).peak_params;      % [center height bw]
                    if ~isempty(pp)
                        ing = pp(:,1) >= gammaBand(1) & pp(:,1) <= gammaBand(2);
                        if any(ing)
                            sub = pp(ing,:); [~,mx] = max(sub(:,2));
                            gF(w) = sub(mx,1); gP(w) = sub(mx,2); det(w) = 1;
                        end
                    end
                end

                % per-subject, per-condition gamma time-progression figure
                plot_gamma_progression(fAxis, specM, apfit, centerMs, gammaBand, id, od.bestMac, cond, ...
                    fullfile(figDir, ['gammaTimeProgression_' cond '.png']));

                for w = 1:nWin
                    newRows(end+1,:) = {id, char(T.subID(i)), T.sessNum(i), char(T.type(i)), ...
                        char(T.group(i)), cond, w, centerMs(w), gF(w), gP(w), det(w)}; %#ok<AGROW>
                end
                anyDone = true;
                fprintf('  %s: gamma epochs done: %d/%d epochs with peak (n=%d trials)\n', cond, sum(det), nWin, kept);
            end

            if anyDone, doneIDs(end+1) = string(id); end %#ok<AGROW>
            clear od sig;
        catch ME
            fprintf('  FAILED: %s\n', ME.message);
            if ~isempty(ME.stack), fprintf('  @ %s line %d\n', ME.stack(1).name, ME.stack(1).line); end
            try, clear od; catch, end
        end
    end

    % merge into the epoch CSV (replace rows for processed sessions)
    if isempty(newRows)
        fprintf('\nNo sessions produced gamma epochs; CSV unchanged.\n'); return;
    end
    NRT = cell2table(newRows, 'VariableNames', hdr);
    if isfile(epoCsv)
        E = readtable(epoCsv); E.sessID = string(E.sessID);
        E(ismember(E.sessID, doneIDs), :) = [];
        NRT2 = NRT; NRT2.sessID = string(NRT2.sessID);
        E = [E; NRT2];
        writetable(E, epoCsv);
    else
        writetable(NRT, epoCsv);
    end
    fprintf('\ngamma-epoch pass done -> %s\n', epoCsv);
end

% ======================================================================
function opt = build_fooof_opt(freqRange)
% mirror ft_freqanalysis's FOOOF opt construction (brainstorm process_fooof),
% overriding to match thresh_fooof_macBP's settings for comparability.
    ob = getfield(process_fooof('GetDescription'), 'options'); %#ok<GFLD>
    opt = struct();
    opt.freq_range          = freqRange;
    opt.peak_width_limits   = [1 12];
    opt.max_peaks           = 6;
    opt.min_peak_height     = 0;
    opt.aperiodic_mode      = 'knee';
    opt.peak_threshold      = 2;
    opt.return_spectrum     = 1;
    opt.border_threshold    = 1;
    opt.power_line          = 'inf';                 % 2-58 Hz: no line-noise notch
    opt.peak_type           = ob.peaktype.Value;     % 'gaussian'
    opt.proximity_threshold = ob.proxthresh.Value{1};
    opt.guess_weight        = ob.guessweight.Value;  % 'none'
    opt.thresh_after        = true;
    opt.sort_type           = ob.sorttype.Value;
    opt.sort_param          = ob.sortparam.Value;
    opt.sort_bands          = ob.sortbands.Value;
end

% ----------------------------------------------------------------------
function plot_gamma_progression(f, spec, apfit, centerMs, gammaBand, sessID, bestMac, cond, outPng)
    nW = size(spec,1);
    c1 = [0.85 0.60 0.10]; c2 = [0.40 0.00 0.55];   % ochre -> purple
    cols = [linspace(c1(1),c2(1),nW)', linspace(c1(2),c2(2),nW)', linspace(c1(3),c2(3),nW)'];
    fig = figure('visible','off','position',[0 0 860 540]); hold on;
    for w = 1:nW
        periodic = 10*log10(spec(w,:) ./ apfit(w,:));   % dB over aperiodic
        plot(f, periodic, 'color', cols(w,:), 'linewidth', 1.6, ...
            'DisplayName', sprintf('%+d ms', round(centerMs(w))));
    end
    yl = ylim;
    patch([gammaBand(1) gammaBand(2) gammaBand(2) gammaBand(1)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
    set(gca,'children',flipud(get(gca,'children')));
    xlim([f(1) f(end)]); ylim(yl);
    xlabel('Frequency (Hz)'); ylabel('Periodic power (dB over aperiodic)');
    legend('show','location','northeastoutside');
    title(sprintf('%s  %s  %s-odor gamma time-progression (finalOnset epochs)', sessID, bestMac, cond), 'interpreter','none');
    saveas(fig, outPng); close(fig);
end
