
addpath('G:\My Drive\GitHub\ZelanoLabScripts')
addpath 'G:\My Drive\GitHub\slowBreathing'
eegLocs = readtable("G:\My Drive\GitHub\ZelanoLabScripts\eegLocs_standard_coords.csv"); 

figSaveDir = 'G:\My Drive\cZelano\FigsStanford';

curDat = load('R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\260105_OBE_NWU_ZF_1\raw\raw_stimFaces\raw_stimFaces.mat')
curDat = curDat.curDat; 

rspDat = load('R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\260105_OBE_NWU_ZF_1\rawBinaps\rawBinaps_stimFaces\rawBinaps_stimFaces.mat')
rspDat = rspDat.curDat; 


rspDat2 = load('R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\260105_OBE_NWU_ZF_1\raw\raw_focusedBreathing_multiSniff_echem\raw_focusedBreathing_multiSniff_echem.mat')
rspDat2 = rspDat2.curDat; 


outDat = struct; 
outDat.data = curDat.rawData.trial{1}; 

outDat.fs = curDat.rawData.fsample; 
outDat.data(isnan(outDat.data)) = 0; 
outDat.labels = curDat.outLabs;
outDat2 = downsample_data(outDat, 2000); 
test = curDat.rawData.trial{1}(34,:); 

stimRsp = struct; 
stimRsp.data = rspDat.rawData.trial{1}; 
stimRsp.fs = rspDat.rawData.fsample; 
stimRsp.data(isnan(stimRsp.data)) = 0; 
stimRsp = downsample_data(stimRsp, 2000); 

conResp = struct; 
conResp.data = rspDat2.rawData.trial{1}; 
conResp.fs = rspDat2.rawData.fsample; 
conResp.data(isnan(conResp.data)) = 0; 
conResp = downsample_data(conResp, 2000); 


% --- theme colors ---
bgCol    = [26 24 56]/255;      % #1A1838
labCol   = [255 234 177]/255;   % #FFEAB1
rspCol   = [223 230 218]/255;   % sage green
ttlCol   = [190 160 255]/255;   % light purple
ephysCol = [230 110 110]/255;   % softer red

tim = 1/2000 : 1/2000 : length(outDat.data(1,:)) * (1/2000);

%% =========================
% Figure 1
% =========================
fig = figure('Color', bgCol, 'Position', [80 80 1200 550], 'InvertHardcopy', 'off');
ax = axes('Parent', fig);
hold(ax, 'on');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax, 'off');

yyaxis left
plot(tim, outDat.data(44,:) + 2000, 'Color', ttlCol, 'LineWidth', 2.0)
ax.YAxis(1).Color = labCol;

xlabel('time (seconds)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);
ylabel('voltage (\muV)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

ylim([-350 400])
% xlim([20 40])

% relabel left-axis ticks to undo the /10-style display scaling convention
yyaxis left
yt = yticks;
yticklabels(string(round(yt * 10)))

yyaxis right
plot(tim(1:length(stimRsp.data(1,:))), stimRsp.data(1,:), ...
    'Color', rspCol, 'LineWidth', 3.0)
ax.YAxis(2).Color = labCol;

ylabel('respiration', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

lgd = legend({'Trigger Output From Stim Generator', 'Respiration'}, ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'TextColor', labCol, ...
    'Color', bgCol, ...
    'Box', 'off', ...
    'location', 'south');

% title('Respiration and channel 44 signal', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

set(gca, 'Layer', 'top');

outFile = fullfile(figSaveDir, 'stimWide.jpg');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);


%% =========================
% Figure 2
% =========================
fig = figure('Color', bgCol, 'Position', [100 100 1200 550], 'InvertHardcopy', 'off');
ax = axes('Parent', fig);
hold(ax, 'on');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax, 'off');

yyaxis left
plot(tim, outDat.data(44,:) + 2300, 'Color', ttlCol, 'LineWidth', 2.0)
hold on
plot(tim, outDat2.data(34,:) ./ 10, 'Color', ephysCol, 'LineWidth', 1, 'linestyle', '-')
ax.YAxis(1).Color = labCol;

xlabel('time (seconds)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);
ylabel('voltage (\muV)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

ylim([-300 400])
xlim([100 107])

% relabel left-axis ticks so displayed values correspond to 10x-scaled µV
yyaxis left
yt = yticks;
yticklabels(string(round(yt * 10)))

yyaxis right
plot(tim(1:length(stimRsp.data(1,:))), stimRsp.data(1,:), ...
    'Color', rspCol, 'LineWidth', 3.0)
ax.YAxis(2).Color = labCol;

ylabel('respiration', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

lgd = legend({'Trigger Output From Stim Generator', ...
              'Target OB Electrode', ...
              'Respiration'}, ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'TextColor', labCol, ...
    'Color', bgCol, ...
    'Box', 'off', ...
    'location', 'northeast');



% title('Respiration, channel 44 signal, and scaled ephys trace', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

set(gca, 'Layer', 'top');

outFile = fullfile(figSaveDir, 'stimZoom.jpg');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);






%stim Dat focused boundaries: 
oni = [33060*4 199799*4];
offi = [73437*4 239315*4];

rspDat = [stimRsp.data(1,oni(1):offi(1)) stimRsp.data(1,oni(2):offi(2))]; 
StimRespStats = breathTemplates4(rspDat, 2000);
    %col 1: onset Y value
    %col 2: onset tim
    %col 3: peak Y value
    %col 4: peak tim
    %col 5: end Y value
    %col 6: end tim
    %col 7: length (end tim - onset tim)
    %col 8: amp (peak Y - avg of two ends)
    %col 9: idx of peak in rspSig2
    %col10: exhale peak Y value
    %col11: exhale peak tim
    %col12: condition
    %col13: empty
    %col14: index

%control Dat focused boundaries: 

tim = 1/2000:1/2000:length(conResp.data(1,:))*(1/2000); 
figure; plot(tim, conResp.data(42,:))
rspDat2 = conResp.data(42,tim<300); 
ConRespStats = breathTemplates4(rspDat2, 2000);


figure; %plot of breath length: 
histogram(ConRespStats(:,7), [1:.33:8])
hold on 
histogram(StimRespStats(:,7), [1:.33:8])

figure; %plot of breath amplitude: 
histogram(ConRespStats(:,8), [0:.33:4])
hold on 
histogram(StimRespStats(:,8), [0:.33:4])



%time vector in seconds for stim data
stimTim = 1/2000:1/2000:length(rspDat)*(1/2000);
%onset in seconds: 
stimidx = StimRespStats(:,2);
%time vector in seconds for control data
conTim = 1/2000:1/2000:length(rspDat2)*(1/2000);
%onset in seconds: 
conidx = ConRespStats(:,2);


rspDat2 = smoothdata(rspDat2, 'gaussian', 200).*3; 
% ============================================
% Native-time breath overlays: -0.5 to +6 s
% with thick mean traces
% ============================================

% --- theme colors ---
bgCol    = [26 24 56]/255;      % #1A1838
labCol   = [255 234 177]/255;   % #FFEAB1
conCol   = [223 230 218]/255;   % sage green
stimCol  = [230 110 110]/255;   % soft light red

% pseudo-alpha by blending with background
alphaVal = 0.35;
conPlotCol  = alphaVal*conCol  + (1-alphaVal)*bgCol;
stimPlotCol = alphaVal*stimCol + (1-alphaVal)*bgCol;

fs = 2000;

% time vectors in seconds
stimTim = 1/fs : 1/fs : length(rspDat)  * (1/fs);
conTim  = 1/fs : 1/fs : length(rspDat2) * (1/fs);

% onset times in seconds
stimOnSec = StimRespStats(:,2);
conOnSec  = ConRespStats(:,2);

% onset sample indices
stimidx = round(stimOnSec * fs);
conidx  = round(conOnSec  * fs);

% keep valid
stimidx = stimidx(isfinite(stimidx) & stimidx >= 1 & stimidx <= numel(rspDat));
conidx  = conidx(isfinite(conidx)  & conidx  >= 1 & conidx  <= numel(rspDat2));

stimidx = unique(stimidx(:));
conidx  = unique(conidx(:));

% epoch: -0.5s to +6s
winSamp = round(-0.5*fs) : round(6*fs);
tRel    = winSamp / fs;
nT      = numel(winSamp);

% preallocate epoch matrices for means
stimMat = nan(numel(stimidx), nT);
conMat  = nan(numel(conidx),  nT);

% --- figure ---
fig = figure('Color', bgCol, 'Position', [80 80 600 650], 'InvertHardcopy', 'off');
ax = axes('Parent', fig);
hold(ax, 'on');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax, 'off');

% ---------- control breaths ----------
for k = 1:numel(conidx)
    idx0 = conidx(k);
    useIdx = idx0 + winSamp;

    if min(useIdx) < 1 || max(useIdx) > numel(rspDat2)
        continue
    end

    seg = double(rspDat2(useIdx));
    if any(~isfinite(seg))
        continue
    end

    conMat(k,:) = seg;
    plot(tRel, seg, 'Color', [conCol alphaVal], 'LineWidth', 1.5)
end

% ---------- stim breaths ----------
for k = 1:numel(stimidx)
    idx0 = stimidx(k);
    useIdx = idx0 + winSamp;

    if min(useIdx) < 1 || max(useIdx) > numel(rspDat)
        continue
    end

    seg = double(rspDat(useIdx));
    if any(~isfinite(seg))
        continue
    end

    stimMat(k,:) = seg;
    plot(tRel, seg, 'Color', [stimCol alphaVal], 'LineWidth', 1.5)
end

% ---------- thick mean traces ----------
conMean  = mean(conMat,  1, 'omitnan');
stimMean = mean(stimMat, 1, 'omitnan');

plot(tRel, conMean,  'Color', conCol,  'LineWidth', 6.0);
plot(tRel, stimMean, 'Color', stimCol, 'LineWidth', 6.0);

% onset reference line
xline(0, '--', 'Color', labCol, 'LineWidth', 2);

xlabel('time from breath onset (seconds)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

ylabel('respiration', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

% title('Breath-aligned respiration traces', ...
%     'FontSize', 20, ...
%     'FontWeight', 'bold', ...
%     'FontName', 'Dotum', ...
%     'Color', labCol);

xlim([-0.5 6])

% legend handles
h1 = plot(nan, nan, 'Color', conCol,  'LineWidth', 4);
h2 = plot(nan, nan, 'Color', stimCol, 'LineWidth', 4);

lgd = legend([h1 h2], {'Non-stim session', 'Stim session'}, ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'TextColor', labCol, ...
    'Color', bgCol, ...
    'Box', 'off', ...
    'Location', 'northeast');

set(gca, 'Layer', 'top');


outFile = fullfile(figSaveDir, 'stimBreathTraceEffect.jpg');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);




% --- theme colors ---
bgCol   = [26 24 56]/255;      % #1A1838
labCol  = [255 234 177]/255;   % #FFEAB1
conCol  = [223 230 218]/255;   % sage green
stimCol = [230 110 110]/255;   % soft light red

%% =========================
% Figure 1: breath length
% =========================
fig = figure('Color', bgCol, 'Position', [80 80 950 600], 'InvertHardcopy', 'off');
ax = axes('Parent', fig);
hold(ax, 'on');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax, 'off');

edges1 = 1:0.33:8;

histogram(ConRespStats(:,7), edges1, ...
    'FaceColor', conCol, ...
    'EdgeColor', conCol, ...
    'FaceAlpha', 0.65, ...
    'LineWidth', 1.5);

histogram(StimRespStats(:,7), edges1, ...
    'FaceColor', stimCol, ...
    'EdgeColor', stimCol, ...
    'FaceAlpha', 0.65, ...
    'LineWidth', 1.5);

xlabel('breath length (s)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

ylabel('count', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

title('Breath length distribution', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

lgd = legend({'Control', 'Stim'}, ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'TextColor', labCol, ...
    'Color', bgCol, ...
    'Box', 'off', ...
    'Location', 'northeast');

set(gca, 'Layer', 'top');

outFile = fullfile(figSaveDir, 'stimBreathLengthHist.jpg');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);

%% =========================
% Figure 2: breath amplitude
% =========================
fig = figure('Color', bgCol, 'Position', [100 100 950 600], 'InvertHardcopy', 'off');
ax = axes('Parent', fig);
hold(ax, 'on');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax, 'off');

edges2 = 0:0.33:4;

histogram(ConRespStats(:,8), edges2, ...
    'FaceColor', conCol, ...
    'EdgeColor', conCol, ...
    'FaceAlpha', 0.65, ...
    'LineWidth', 1.5);

histogram(StimRespStats(:,8), edges2, ...
    'FaceColor', stimCol, ...
    'EdgeColor', stimCol, ...
    'FaceAlpha', 0.65, ...
    'LineWidth', 1.5);

xlabel('breath amplitude (a.u.)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

ylabel('count', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

title('Breath amplitude distribution', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

lgd = legend({'Control', 'Stim'}, ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'TextColor', labCol, ...
    'Color', bgCol, ...
    'Box', 'off', ...
    'Location', 'northeast');

set(gca, 'Layer', 'top');


outFile = fullfile(figSaveDir, 'stimBreathAmpHist.jpg');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);




















figure; plot(stimTim(stimTim>20 & stimTim<70), rspDat(stimTim>20 & stimTim<70))
hold on 
plot(conTim(conTim>20 & conTim<70), smoothdata(rspDat2(conTim>20 & conTim<70), 'gaussian', 200).*3)




%% =========================
%  Detect stim events + plot ERPs across scalp channels
%  Assumes:
%    outDat.data      : [nChan x nTime]
%    outDat.labels    : channel labels in row order of outDat.data
%    eegLocs          : table with columns Label, X2D_right, Y2D_front
% ==========================

fs = 500;                 % Hz
ttlRow = 44;              % TTL signal row
eegRows = 1:32;           % scalp EEG rows
thr = -1400;              % detection threshold

tWin = [-1 4];            % seconds relative to stim
sampWin = round(tWin * fs);
winSamps = sampWin(1):sampWin(2);
t = winSamps / fs;

%% --- Convert labels into clean cell arrays of char/string ---
% outDat.labels may be char matrix, cellstr, or string array
if ischar(outDat.labels)
    outLabels = cellstr(outDat.labels);
elseif isstring(outDat.labels)
    outLabels = cellstr(outDat.labels);
else
    outLabels = outDat.labels;
end
outLabels = strtrim(outLabels);

% eegLocs.Label may also vary in type
if ischar(eegLocs.Label)
    locLabels = cellstr(eegLocs.Label);
elseif isstring(eegLocs.Label)
    locLabels = cellstr(eegLocs.Label);
elseif iscategorical(eegLocs.Label)
    locLabels = cellstr(string(eegLocs.Label));
else
    locLabels = eegLocs.Label;
end
locLabels = strtrim(locLabels);

%% --- Detect TTL crossings: above -> below -1400 ---
ttl = double(outDat.data(ttlRow,:));
stimIdx = find(ttl(1:end-1) > thr & ttl(2:end) <= thr) + 1;

fprintf('Detected %d stimulation events.\n', numel(stimIdx));

%% --- Keep only events with a full epoch available ---
nTime = size(outDat.data,2);
goodStim = stimIdx + sampWin(1) >= 1 & stimIdx + sampWin(2) <= nTime;
stimIdx = stimIdx(goodStim);

fprintf('Kept %d events with full [%0.1f %0.1f] s windows.\n', numel(stimIdx), tWin(1), tWin(2));

%% --- Extract EEG epochs: [channel x time x event] ---
nChan = numel(eegRows);
nEvt  = numel(stimIdx);
nSamp = numel(winSamps);

epochs = nan(nChan, nSamp, nEvt);

for e = 1:nEvt
    idx = stimIdx(e) + winSamps;
    epochs(:,:,e) = double(outDat.data(eegRows, idx));
end

%% --- Optional baseline correction using -200 to 0 ms ---
doBaseline = true;
if doBaseline
    baseIdx = t >= -0.2 & t <= 0;
    for e = 1:nEvt
        b = mean(epochs(:,baseIdx,e), 2, 'omitnan');
        epochs(:,:,e) = epochs(:,:,e) - b;
    end
end

%% --- Average across events to get ERP ---
erp = mean(epochs, 3, 'omitnan');   % [channel x time]

%% --- Match EEG rows (outDat 1:32) to eegLocs rows by label ---
eegLabels = outLabels(eegRows);

[isMatched, locIdx] = ismember(upper(eegLabels), upper(locLabels));

fprintf('Matched %d/%d EEG channels to eegLocs.\n', sum(isMatched), nChan);

% Keep only matched channels for plotting
erpPlot    = erp(isMatched, :);
plotLabels = eegLabels(isMatched);
plotLocIdx = locIdx(isMatched);

x = eegLocs.X2D_right(plotLocIdx);
y = eegLocs.Y2D_front(plotLocIdx);

%% --- Normalize scalp coordinates into figure positions ---
x = double(x);
y = double(y);

xNorm = (x - min(x)) / (max(x) - min(x) + eps);
yNorm = (y - min(y)) / (max(y) - min(y) + eps);

% axis box size in normalized figure coordinates
w = 0.12;
h = 0.10;

% margins
leftMargin   = 0.06;
rightMargin  = 0.06;
bottomMargin = 0.08;
topMargin    = 0.08;

xPos = leftMargin + xNorm * (1 - leftMargin - rightMargin - w);
yPos = bottomMargin + yNorm * (1 - bottomMargin - topMargin - h);

%% --- Common y-limits across channels ---
yl = [min(erpPlot(:), [], 'omitnan'), max(erpPlot(:), [], 'omitnan')];
if any(~isfinite(yl)) || yl(1) == yl(2)
    yl = [-1 1];
end

%% --- Plot scalp-layout ERPs ---
figure('Color', 'w', 'Name', 'Stimulus-locked ERPs');

for ch = 1:size(erpPlot,1)
    ax = axes('Position', [xPos(ch), yPos(ch), w, h]); %#ok<LAXES>
    plot(t, erpPlot(ch,:), 'k', 'LineWidth', 1);
    hold on
    xline(0, 'r--', 'LineWidth', 1);
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    xlim(tWin);
    ylim([-20 20]);
    title(plotLabels{ch}, 'FontSize', 8, 'Interpreter', 'none');
    set(ax, 'FontSize', 7, 'Box', 'off', 'XTick', [-1 0 2 4]);

    % Only keep full tick labels for a few channels to reduce clutter
    if ch ~= 1
        ax.XTickLabel = [];
        ax.YTickLabel = [];
    else
        xlabel('Time (s)');
        ylabel('\muV');
    end
end

sgtitle(sprintf('ERPs relative to stimulation (%d events)', nEvt));

%% --- Optional: also plot all ERPs in one standard stacked/tiled figure ---
figure('Color', 'w', 'Name', 'ERPs tiled by channel');
tl = tiledlayout(6,6, 'TileSpacing', 'compact', 'Padding', 'compact');

for ch = 1:nChan
    nexttile
    plot(t, erp(ch,:), 'k', 'LineWidth', 1); hold on
    xline(0, 'r--');
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    xlim(tWin);
    ylim(yl);
    title(eegLabels{ch}, 'FontSize', 8, 'Interpreter', 'none');
    set(gca, 'FontSize', 7, 'Box', 'off');
end
title(tl, sprintf('Stimulus-locked ERPs (%d events)', nEvt));