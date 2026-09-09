function run_o15_gamma_epochs(sessFilter)
% RUN_O15_GAMMA_EPOCHS  Time-resolved gamma via FOOOF for the O15 task, computed
%   SEPARATELY for each of the three sniff categories (start / free / confirm).
%   Per session x sniff type (bestMac, that-type-locked, noise-clean sniffs):
%   Morlet power 5-58 Hz (100 linear freqs), averaged over trials & 250 ms windows
%   into 10 sequential power spectra (-500..+2000 ms); FOOOF each; the periodic
%   spectra are plotted (ochre->purple) to <figs>/gammaTimeProgression_<type>.png
%   and the largest 30-58 Hz peak per epoch is written to O15Task_gammaEpochs.csv
%   (with a `locking` column = start/free/confirm). There are far more free sniffs
%   than start/confirm; that is expected and fine.

    if nargin < 1, sessFilter = []; end
    o15_init_paths(); L = labPaths();
    groupDir = getenv('O15_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(labPaths().figPath, 'groupStatFigs'); end
    epoCsv   = fullfile(groupDir, 'O15Task_gammaEpochs.csv');

    epWin   = [-1.75 5.75];
    cycles  = [3 0.8];
    winEdges = -500:250:2000;
    nWin    = numel(winEdges) - 1;
    gammaBand = [30 58];
    opt = build_fooof_opt([5 58]);
    hasOpt = exist('fmincon','file') > 0;
    types = {'start','free','confirm'};
    Knoise = o15_noise_K();   % relative sharp-deflection threshold (robust-SDs), single source

    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    hdr = {'sessID','subID','sessNum','type','group','locking','epoch','centerMs', ...
           'gammaPeakFreq','gammaPeakPower','gammaDetected'};
    newRows = {}; doneIDs = strings(0,1);

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n== %d/%d %s ==\n', i, height(T), id);
        try
            Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
            if ~isfield(od,'bestMac') || isempty(od.bestMac), fprintf('  no bestMac -> skip\n'); clear od; continue; end
            labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
            ci = find(strcmp(od.bestMac, labs), 1); if isempty(ci), clear od; continue; end
            sig = double(od.data(ci, :)); fs = od.fs; Tn = numel(sig);
            bd = od.behDat; sl = string(bd.sniffLabel); fo = bd.finalOnset;

            figDir = fullfile(L.figPath, id, 'O15'); if ~isfolder(figDir), mkdir(figDir); end
            legacy = fullfile(figDir, 'gammaTimeProgression.png');   % start-only file from the prior version
            if isfile(legacy), delete(legacy); end

            anyType = false;
            for tj = 1:numel(types)
                tp = types{tj};
                onsAll = fo(sl == tp);
                NT = cue_noise_trials(sig, fs, onsAll, [], Knoise);   % RELATIVE rule (zd>K)
                ons = onsAll(NT.ok & ~NT.noisy);

                s0 = round(epWin(1)*fs); s1 = round(epWin(2)*fs); nF = s1 - s0 + 1;
                ep = zeros(nF, numel(ons)); kept = 0;
                for e = ons(:)'
                    a = round(e)+s0; b = round(e)+s1;
                    if a < 1 || b > Tn, continue; end
                    seg = sig(a:b); if any(~isfinite(seg)), continue; end
                    kept = kept + 1; ep(:,kept) = seg(:);
                end
                ep = ep(:, 1:kept);
                if kept < 3, fprintf('  %-7s: too few clean sniffs (%d) -> skip\n', tp, kept); continue; end

                [~,~,~,times,fout,~,~,atf] = newtimef(ep, nF, epWin*1000, fs, cycles, ...
                    'freqs', [5 58], 'nfreqs', 100, 'freqscale', 'linear', ...
                    'baseline', NaN, 'plotersp','off','plotitc','off','verbose','off');
                P = abs(atf).^2; meanP = mean(P, 3, 'omitnan'); fU = fout(:)';

                spec = zeros(nWin, numel(fU)); centerMs = zeros(nWin,1);
                for w = 1:nWin
                    tmask = times >= winEdges(w) & times < winEdges(w+1);
                    spec(w,:) = mean(meanP(:, tmask), 2, 'omitnan')';
                    centerMs(w) = mean(winEdges(w:w+1));
                end

                TF = zeros(nWin, 1, numel(fU)); TF(:,1,:) = spec;
                [~, fg] = process_fooof('FOOOF_matlab', TF, fU, opt, hasOpt);
                fMask = fU >= opt.freq_range(1) & fU <= opt.freq_range(2);
                fAxis = fU(fMask); specM = spec(:, fMask);

                gF = nan(nWin,1); gP = nan(nWin,1); det = zeros(nWin,1);
                apfit = zeros(nWin, numel(fAxis));
                for w = 1:nWin
                    apfit(w,:) = fg(w).ap_fit(:)';
                    pp = fg(w).peak_params;
                    if ~isempty(pp)
                        ing = pp(:,1) >= gammaBand(1) & pp(:,1) <= gammaBand(2);
                        if any(ing)
                            sub = pp(ing,:); [~,mx] = max(sub(:,2));
                            gF(w) = sub(mx,1); gP(w) = sub(mx,2); det(w) = 1;
                        end
                    end
                end

                plot_gamma_progression(fAxis, specM, apfit, centerMs, gammaBand, id, od.bestMac, ...
                    tp, kept, fullfile(figDir, ['gammaTimeProgression_' tp '.png']));

                for w = 1:nWin
                    newRows(end+1,:) = {id, char(T.subID(i)), T.sessNum(i), char(T.type(i)), ...
                        char(T.group(i)), tp, w, centerMs(w), gF(w), gP(w), det(w)}; %#ok<AGROW>
                end
                anyType = true;
                fprintf('  %-7s: %d/%d epochs with peak (n=%d sniffs)\n', tp, sum(det), nWin, kept);
            end
            if anyType, doneIDs(end+1) = string(id); end %#ok<AGROW>
            clear od sig P meanP spec atf;
        catch ME
            fprintf('  FAILED: %s\n', ME.message); try, clear od; catch, end
        end
    end

    if isempty(newRows)
        fprintf('\nNo sessions produced gamma epochs; CSV unchanged.\n'); return;
    end
    NRT = cell2table(newRows, 'VariableNames', hdr); NRT.sessID = string(NRT.sessID);
    if isfile(epoCsv)
        E = readtable(epoCsv); E.sessID = string(E.sessID);
        % schema changed (added `locking`) or mismatched columns -> rebuild from this run
        if ~ismember('locking', E.Properties.VariableNames) || ...
           ~isequal(sort(E.Properties.VariableNames), sort(hdr))
            E = NRT;
        else
            E(ismember(E.sessID, doneIDs), :) = [];
            E = [E; NRT];
        end
        writetable(E, epoCsv);
    else
        writetable(NRT, epoCsv);
    end
    fprintf('\nO15 per-sniff-type gamma-epoch pass done -> %s\n', epoCsv);
end

% ======================================================================
function opt = build_fooof_opt(freqRange)
    ob = getfield(process_fooof('GetDescription'), 'options'); %#ok<GFLD>
    opt = struct();
    opt.freq_range          = freqRange;
    opt.peak_width_limits   = [1 12];
    opt.max_peaks           = 6;
    opt.min_peak_height     = 0;
    opt.aperiodic_mode      = 'fixed';
    opt.peak_threshold      = 2;
    opt.return_spectrum     = 1;
    opt.border_threshold    = 1;
    opt.power_line          = 'inf';
    opt.peak_type           = ob.peaktype.Value;
    opt.proximity_threshold = ob.proxthresh.Value{1};
    opt.guess_weight        = ob.guessweight.Value;
    opt.thresh_after        = true;
    opt.sort_type           = ob.sorttype.Value;
    opt.sort_param          = ob.sortparam.Value;
    opt.sort_bands          = ob.sortbands.Value;
end

% ----------------------------------------------------------------------
function plot_gamma_progression(f, spec, apfit, centerMs, gammaBand, sessID, bestMac, tp, nSniff, outPng)
    nW = size(spec,1);
    c1 = [0.85 0.60 0.10]; c2 = [0.40 0.00 0.55];
    cols = [linspace(c1(1),c2(1),nW)', linspace(c1(2),c2(2),nW)', linspace(c1(3),c2(3),nW)'];
    fig = figure('visible','off','position',[0 0 860 540]); hold on;
    for w = 1:nW
        periodic = 10*log10(spec(w,:) ./ apfit(w,:));
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
    title(sprintf('%s  %s  gamma time-progression (%s-sniff epochs, n=%d)', sessID, bestMac, tp, nSniff), ...
        'interpreter','none');
    saveas(fig, outPng); close(fig);
end
