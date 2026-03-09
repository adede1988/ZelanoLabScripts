%% scratch figure generation code for methods showing: 
chanDat = load('R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\CHANDAT_processed\250904_OBE_NWU_TI_macro_breathingTask_48.mat');
rawDat = load('R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\CHANDAT\250904_OBE_NWU_TI_macro_breathingTask_48.mat');

figSaveDir = 'G:\My Drive\cZelano\FigsStanford';
chanDat = chanDat.chanDat; 
rawDat  = rawDat.chanDat; 



fb     = chanDat.pac.diag.fb; 
bOnset = chanDat.behDat.finalOnset; 

set(0, 'defaultfigurewindowstyle', 'docked')

figure; plot(movmean(fb, 10))

exOn = bOnset(82);




% prep the gamma data: 
gamMed = median(chanDat.fooof.gamma_peaks, 'all', 'omitnan')

fs = chanDat.fs;               % Hz
    
halfBW = 5;                    % +/- 5 Hz
bp   = double([gamMed-halfBW, gamMed+halfBW]);

xAmpSrc = rawDat.data; 

bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', bp(1), ...
    'HalfPowerFrequency2', bp(2), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');

xBP = filtfilt(bpFilt, xAmpSrc);
A   = abs(hilbert(xBP));

lpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency', min(20, fs/2-1), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
A_lp = filtfilt(lpFilt, A);


%prep the respiration data: 

rsp = rawDat.rsp; 
bwMin = 0.06; 
bwFrac = .3; 
f = mean(fb(82:92));


 bw = max(bwMin, bwFrac * f);
        fLo = f - bw/2;
        fHi = f + bw/2;
bpFilt = designfilt('bandpassiir', ...
            'FilterOrder', 8, ...
            'HalfPowerFrequency1', fLo, ...
            'HalfPowerFrequency2', fHi, ...
            'SampleRate', fs, ...
            'DesignMethod', 'butter');

ybp = filtfilt(bpFilt, rsp);

z   = hilbert(ybp);
U   = exp(1i * angle(z));

%% make a figure! 

tim = 1/500: 1/500: 50; 
idx = exOn:exOn+length(tim)-1; 
zoomTo = 4.0; 
zoomLen = 8.5; 
figure
subplot(2,3,[1 2])
plot(tim, rsp(idx))

hold on 
plot(tim, ybp(idx))

% ---------- unit-vector time series on right axis ----------


idx = idx(tim>zoomTo & tim<zoomLen+zoomTo);
tim = tim(tim>zoomTo & tim<zoomLen+zoomTo);
subplot 233
plot(tim, rsp(idx))

hold on 
plot(tim, ybp(idx))
yyaxis right
hold on

Uwin = U(idx);

% plot every 50th vector
step = 50;
sIdx = 1:step:length(tim);

tVec   = tim(sIdx);
y0Vec  = ybp(sIdx); 
theta  = angle(Uwin(sIdx));   % extract phase angle explicitly

% choose displayed vector length
L = 0.5;

% trig-based endpoint offsets
dx = L * cos(theta);
dy = L * sin(theta);

% draw baseline
plot([tim(1) tim(end)], [0 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)

% draw vectors as little line segments
for k = 1:numel(tVec)
    x0 = tVec(k);
    y0 = 0;

    x1 = x0 + dx(k);
    y1 = y0 + dy(k);

    plot([x0 x1], [y0 y1], 'Color', [0.1 0.45 0.95], 'LineWidth', 1.5)
    plot(x1, y1, '.', 'Color', [0.1 0.45 0.95], 'MarkerSize', 10)
end

rVal = max(tim) - min(tim);
ylim([-rVal/2 rVal/2])
ylabel('Phase unit vectors')





subplot 234
plot(tim, rsp(idx))

hold on 
plot(tim, ybp(idx))

plot(tim, movmean(A_lp(idx)*15-200,100))

yyaxis right
hold on

Uwin = U(idx);

% plot every 50th vector
step = 50;
sIdx = 1:step:length(tim);
yScaleVec = A_lp(idx); 
yScaleVec = movmean(yScaleVec, 100); 
tVec   = tim(sIdx);
yScaleVec  = yScaleVec(sIdx); 
yScaleVec = yScaleVec - min(yScaleVec); 
yScaleVec = yScaleVec ./ prctile(yScaleVec, 80); 
theta  = angle(Uwin(sIdx));   % extract phase angle explicitly

% choose displayed vector length
L = yScaleVec .*.5;

% trig-based endpoint offsets
dx = L .* cos(theta);
dy = L .* sin(theta);

% draw baseline
plot([tim(1) tim(end)], [0 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)

% draw vectors as little line segments
for k = 1:numel(tVec)
    x0 = tVec(k);
    y0 = 0;

    x1 = x0 + dx(k);
    y1 = y0 + dy(k);

    plot([x0 x1], [y0 y1], 'Color', [0.1 0.45 0.95], 'LineWidth', 1.5)
    plot(x1, y1, '.', 'Color', [0.1 0.45 0.95], 'MarkerSize', 10)
end

rVal = max(tim) - min(tim);
ylim([-rVal/2 rVal/2])
ylabel('Phase unit vectors')


% show it polar: 
% --- sample the unit vectors you want to display ---

axTmp = subplot(2,3,5);
pos = axTmp.Position;
delete(axTmp)

% now place polar axes in that same spot
pax = polaraxes('Position', pos);

step = 50;
sIdx = 1:step:length(Uwin);

Usub  = Uwin(sIdx);
theta = angle(Usub);              % radians
rFix  = 1;                        % all vectors have equal length


hold(pax, 'on')

for k = 1:numel(theta)
    polarplot(pax, [theta(k) theta(k)], [0 rFix], ...
        'LineWidth', 1.5, 'Color', [0.1 0.45 0.95]);
end

rlim(pax, [0 1.05])
rticks(pax, [0.5 1.0])
title('Unit vectors on polar axes')


% --- sample the unit vectors and matching scale values ---
axTmp = subplot(2,3,6);
pos = axTmp.Position;
delete(axTmp)

% now place polar axes in that same spot
pax = polaraxes('Position', pos);

step = 50;
sIdx = 1:step:length(Uwin);

Usub  = Uwin(sIdx);
theta = angle(Usub);
rVal  = yScaleVec;   % scaled length for each vector

hold(pax, 'on')

for k = 1:numel(theta)
    polarplot(pax, [theta(k) theta(k)], [0 rVal(k)], ...
        'LineWidth', 1.5, 'Color', [0.1 0.45 0.95]);
    polarplot(pax, theta(k), rVal(k), 'o', ...
        'MarkerFaceColor', [0.1 0.45 0.95], ...
        'MarkerEdgeColor', [0.1 0.45 0.95], ...
        'MarkerSize', 4);
end

rlim(pax, [0 max(rVal)*1.05])
title('Scaled vectors on polar axes')

%% =========================
%  Styled single-panel figures
%  =========================

% --- colors / theme ---
bgCol    = [26 24 56]/255;      % #1A1838
labCol   = [255 234 177]/255;   % #FFEAB1

rawCol   = [	223	230	218]/255;    % light sage green
bpCol    = [255  92  92]/255;   % warm coral-red
vecCol   = [190 160 255]/255;   % light purple
ampCol = [153 236 230]/255;   % light turquoise
baseCol  = [160 150 190]/255;   % muted lavender-gray
gridCol  = [110 105 145]/255;   % dim grid / polar grid

% --- hard-coded analysis setup (kept as-is) ---
tim = 1/500 : 1/500 : 50;
idx = exOn : exOn + length(tim) - 1;

zoomTo  = 4.0;
zoomLen = 8.5;
set(0, 'defaultfigurewindowstyle', 'normal')
% =========================
% FIGURE 1: full respiration + bandpassed respiration
% =========================
fig = figure('Color', bgCol, 'Position', [50 50 1100 550], 'InvertHardcopy', 'off');
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

plot(tim, rsp(idx), 'Color', rawCol, 'LineWidth', 2.0)
plot(tim, ybp(idx), 'Color', bpCol,  'LineWidth', 2.8)

xlabel('time (seconds)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);
ylabel('Respiration (L/min)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

% title('Respiration signal and respiration-frequency bandpassed signal', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

legend({'raw respiration','bandpassed respiration'}, ...
    'TextColor', labCol, 'Color', bgCol, 'EdgeColor', 'none', ...
    'FontName', 'Dotum', 'FontSize', 14, 'Location', 'northeast');

set(gca, 'Layer', 'top');

outFile = fullfile(figSaveDir, 'bandBassBreathExamp.png');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);

% =========================
% prepare zoomed interval
% =========================
idx = idx(tim > zoomTo & tim < zoomLen + zoomTo);
tim = tim(tim > zoomTo & tim < zoomLen + zoomTo);

Uwin = U(idx);

% =========================
% FIGURE 2: zoomed respiration + unit vectors over time
% =========================
fig = figure('Color', bgCol, 'Position', [80 80 900 550], 'InvertHardcopy', 'off');
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
plot(tim, rsp(idx), 'Color', rawCol, 'LineWidth', 1.8)
hold on
plot(tim, ybp(idx), 'Color', bpCol,  'LineWidth', 2.8)
ax.YColor = labCol;
ylabel('Respiration / bandpassed respiration', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

yyaxis right
hold on

% plot every 50th vector
step = 50;
sIdx = 1:step:length(tim);

tVec  = tim(sIdx);
theta = angle(Uwin(sIdx));

% choose displayed vector length
L = 0.8;

% trig-based endpoint offsets
dx = L * cos(theta);
dy = L * sin(theta);

% draw baseline
% plot([tim(1) tim(end)], [0 0], 'Color', baseCol, 'LineWidth', 2.0)

% draw vectors
for k = 1:numel(tVec)
    x0 = tVec(k);
    y0 = 0;
    x1 = x0 + dx(k);
    y1 = y0 + dy(k);

    plot([x0 x1], [y0 y1], '.','Color', vecCol, 'LineWidth', 2.2, 'linestyle','-')
    % plot(x1, y1, '.', 'Color', vecCol, 'MarkerSize', 16, ...
    %     'Marker','.')
end

rVal = max(tim) - min(tim);
ylim([-rVal/2 rVal/2])
ax.YColor = labCol;
ax.YAxis(2).Color = bgCol;
ax.YAxis(1).Color = bgCol;
yticks([])
% ylabel('Phase unit vectors', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

xlabel('time (seconds)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);
legend({'raw respiration','bandpassed respiration', 'Directional Unit Vectors'}, ...
    'TextColor', labCol, 'Color', bgCol, 'EdgeColor', 'none', ...
    'FontName', 'Dotum', 'FontSize', 14, 'Location', 'southeast');
% title('Instantaneous phase represented as equal-length unit vectors', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

set(gca, 'Layer', 'top');

outFile = fullfile(figSaveDir, 'breathZoom_unitVectors.png');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);

% =========================
% FIGURE 3: zoomed respiration + amplitude-scaled vectors over time
% =========================
fig = figure('Color', bgCol, 'Position', [110 110 900 550], 'InvertHardcopy', 'off');
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
plot(tim, rsp(idx), 'Color', rawCol, 'LineWidth', 1.8)
hold on
plot(tim, ybp(idx), 'Color', bpCol,  'LineWidth', 2.8)
plot(tim, movmean(A_lp(idx)*15-200,100), 'Color', ampCol, 'LineWidth', 2.6)
ax.YColor = labCol;
ylabel('Respiration / bandpassed / scaled amplitude trace', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

yyaxis right
hold on

step = 50;
sIdx = 1:step:length(tim);

yScaleVec = A_lp(idx);
yScaleVec = movmean(yScaleVec, 100);
tVec      = tim(sIdx);
yScaleVec = yScaleVec(sIdx);
yScaleVec = yScaleVec - min(yScaleVec);
yScaleVec = yScaleVec ./ prctile(yScaleVec, 80);
theta     = angle(Uwin(sIdx));

% choose displayed vector lengths
L = yScaleVec .* 1;

dx = L .* cos(theta);
dy = L .* sin(theta);

% plot([tim(1) tim(end)], [0 0], 'Color', baseCol, 'LineWidth', 2.0)

for k = 1:numel(tVec)
    x0 = tVec(k);
    y0 = 0;
    x1 = x0 + dx(k);
    y1 = y0 + dy(k);

    plot([x0 x1], [y0 y1],'.', 'Color', vecCol, 'LineWidth', 2.2,...
        'linestyle', '-')
    % plot(x1, y1, '.', 'Color', vecCol, 'MarkerSize', 16)
end

rVal = max(tim) - min(tim);
ylim([-rVal/2 rVal/2])
ax.YColor = labCol;
ax.YAxis(2).Color = bgCol;
ax.YAxis(1).Color = bgCol;
% ylabel('Amplitude-weighted phase vectors', ...
    % 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

xlabel('time (seconds)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

legend({'raw respiration','bandpassed respiration', ...
    '\gamma Amplitude', '\gamma Amplitude Scaled Vectors'}, ...
    'TextColor', labCol, 'Color', bgCol, 'EdgeColor', 'none', ...
    'FontName', 'Dotum', 'FontSize', 14, 'Location', 'southeast');

outFile = fullfile(figSaveDir, 'breathZoom_scaledVectors.png');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);

% title('Instantaneous phase vectors scaled by gamma amplitude', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

set(gca, 'Layer', 'top');

% =========================
% FIGURE 4: polar plot of equal-length unit vectors
% =========================
fig = figure('Color', bgCol, 'Position', [140 140 700 700], 'InvertHardcopy', 'off');
pax = polaraxes('Parent', fig);
hold(pax, 'on');

step = 50;
sIdx = 1:step:length(Uwin);

Usub  = Uwin(sIdx);
theta = angle(Usub);
rFix  = 1;

for k = 1:numel(theta)
    polarplot(pax, [theta(k) theta(k)], [0 rFix], ...
        'LineWidth', 2.0, 'Color', vecCol);
end

rlim(pax, [0 1.05])
rticks(pax, [0.5 1.0])

pax.Color          = bgCol;
pax.FontSize       = 16;
pax.FontWeight     = 'bold';
pax.FontName       = 'Dotum';
pax.ThetaColor     = labCol;
pax.RColor         = labCol;
pax.GridColor      = gridCol;
pax.MinorGridColor = gridCol;
pax.LineWidth      = 2.0;
outFile = fullfile(figSaveDir, 'polar_unitVectors.png');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);
% title(pax, 'Equal-length phase unit vectors on polar axes', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

% =========================
% FIGURE 5: polar plot of amplitude-scaled vectors
% =========================
fig = figure('Color', bgCol, 'Position', [170 170 700 700], 'InvertHardcopy', 'off');
pax = polaraxes('Parent', fig);
hold(pax, 'on');

step = 50;
sIdx = 1:step:length(Uwin);

Usub  = Uwin(sIdx);
theta = angle(Usub);
rVal  = yScaleVec;

for k = 1:numel(theta)
    polarplot(pax, [theta(k) theta(k)], [0 rVal(k)], ...
        'LineWidth', 2.0, 'Color', vecCol);
    polarplot(pax, theta(k), rVal(k), 'o', ...
        'MarkerFaceColor', vecCol, ...
        'MarkerEdgeColor', vecCol, ...
        'MarkerSize', 5);
end

rlim(pax, [0 max(rVal)*1.05])

pax.Color          = bgCol;
pax.FontSize       = 16;
pax.FontWeight     = 'bold';
pax.FontName       = 'Dotum';
pax.ThetaColor     = labCol;
pax.RColor         = labCol;
pax.GridColor      = gridCol;
pax.MinorGridColor = gridCol;
pax.LineWidth      = 2.0;
outFile = fullfile(figSaveDir, 'polar_scaledVectors.png');

exportgraphics(fig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);
% title(pax, 'Amplitude-scaled phase vectors on polar axes', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol); 