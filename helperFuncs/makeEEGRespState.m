function [fig, itpcout, powout] = makeEEGRespState(chanDat, taskVec, conds, opts)
%MAKEEEGTF  Plot mean time-frequency (powZ) heatmaps by condition.
%
%   fig = makeEEGTF(chanDat, taskVec, conds)
%   fig = makeEEGTF(chanDat, taskVec, conds, opts)
%
% INPUTS
%   chanDat.tf.powZ : [nBreath x nTime x nFreq] z-scored power TF
%   chanDat.tf.tVec : [1 x nTime] time vector (optional)
%   chanDat.tf.frex : [1 x nFreq] frequency vector (optional)
%   chanDat.use     : [nBreath x 1] logical/0-1
%   taskVec         : string array length nBreath with condition labels
%   conds           : string array, e.g. ["audio","focus","shadow"]
%   opts (optional struct):
%       .clim       (default [-5 5])
%       .colormap   (default parula)
%       .minN       (default 5) minimum breaths to plot a condition
%
% OUTPUT
%   fig : figure handle

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'clim')     || isempty(opts.clim),     opts.clim = [-5 5]; end
if ~isfield(opts,'colormap') || isempty(opts.colormap), opts.colormap = parula; end
if ~isfield(opts,'minN')     || isempty(opts.minN),     opts.minN = 5; end

powZ = chanDat.tf.powZ;  % [breath x time x freq]

[nBreath, nTime, nFreq] = size(powZ);

freqBreaks = [2, 8; ...
              8, 14; ...
              14,30];
freqLabs = {'theta','alpha','beta'};

% colors: forest green, burnt orange, dark mauve
frexCols = [0.13 0.55 0.13;
            0.80 0.33 0.00;
            0.55 0.41 0.53];

HRV_RMS  = chanDat.behDat.HRV_RMSSD30;
HRV_SDNN = chanDat.behDat.HRV_SDNN30;
HRV_RSA  = chanDat.behDat.HRV_RSAamp; 


useVec = logical(chanDat.use(:));

taskVec = string(taskVec(:));

conds = string(conds(:))';
blen = chanDat.behDat.length(:); 
bamp = chanDat.behDat.amp(:); 
tVec = 1:nTime;
epochNames = {'inhale rise','inhale fall','exhale rise','exhale fall','pause'};
phaseEdges = [10.5 20.5 30.5 40.5]; 

fVec = double(chanDat.tf.frex(:))';


fig = figure('Color','w', 'visible', true, 'position', [0,0,1200, 700]);
tlo = tiledlayout(fig, 3, 4, 'Padding','compact', 'TileSpacing','compact');


%% look at the HRV measures against each other

ax = nexttile([1 4]);


% --- subset to used breaths ---
mUse = useVec(:)~=0;
rms  = double(HRV_RMS(mUse));
sdnn = double(HRV_SDNN(mUse));
rsa  = double(HRV_RSA(mUse));
task = taskVec(mUse);

x = (1:numel(rms))';

% --- robust min-max scaling to [0,1] (handles NaNs / constants) ---
minmax01 = @(v) local_minmax01(v);
rms01  = minmax01(rms);
sdnn01 = minmax01(sdnn);
rsa01  = minmax01(rsa);

% --- overall state = mean of the 3 scaled measures ---
overall = mean([rms01 sdnn01 rsa01], 2, 'omitnan');

% --- 66th percentile line (across session / plotted points) ---
ovFinite = overall(isfinite(overall));
p66 = prctile(ovFinite, 66);

% --- separation points between task conditions ---
changeIdx = find(task(2:end) ~= task(1:end-1)) + 1;    % indices where a new task starts
xBound = changeIdx - 0.5;                              % draw between points

% --- plot ---
cla(ax); hold(ax,'on');

h1 = plot(ax, x, rms01,  'LineWidth', 1.25);
h2 = plot(ax, x, sdnn01, 'LineWidth', 1.25);
h3 = plot(ax, x, rsa01,  'LineWidth', 1.25);

hO = plot(ax, x, overall, 'k-', 'LineWidth', 3); % bold overall state

% horizontal threshold
% yline(ax, p66, 'k--', 'LineWidth', 1.75, 'Label', '66th pct', ...
%     'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
%     'Interpreter','none');

% vertical task boundaries
for ii = 1:numel(xBound)
    xline(ax, xBound(ii), ':', 'Color', [0 0 0]*0.35, 'LineWidth', 1.25);
end

% --- styling ---
grid(ax,'on'); ax.GridAlpha = 0.15;
box(ax,'off');
set(ax, 'LineWidth', 1.5, 'XColor',[0 0 0], 'YColor',[0 0 0], 'FontSize', 11, ...
    'TickDir','out', 'TickLength',[0.008 0.008]);

xlabel(ax, 'Breath index', 'Interpreter','none');
ylabel(ax, 'Min–max scaled HRV', 'Interpreter','none');

legend(ax, [h1 h2 h3 hO], {'RMSSD30','SDNN30','RSAamp','Overall state'}, ...
    'Location', 'eastoutside', 'Box','off', 'Interpreter','none',...
    'AutoUpdate', 'off');

ylim(ax, [-0.05 1.2]);
% --- task epoch labels (centered at top of plot) ---
yl = ylim(ax);
yText = yl(2) - 0.02*diff(yl);   % a little below the top edge

% task = taskVec(useVec) from your plotting code (string vector, length = numel(x))
chg = [true; task(2:end) ~= task(1:end-1); true];
ix  = find(chg);
starts = ix(1:end-1);
ends   = ix(2:end)-1;

for k = 1:numel(starts)
    lab = task(starts(k));
    if ~isfinite(starts(k)) || strlength(lab)==0, continue; end

    xc = mean(x([starts(k) ends(k)]));  % center of epoch in x-units
    text(ax, xc, yText, char(lab), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'FontSize', 12, ...
        'Interpreter','none', ...
        'Clipping','on', ...
        'BackgroundColor','w', ...
        'Margin', 2);
end



hold(ax,'off');




taskUse = string(taskVec(mUse));
yLen    = blen(mUse);


% colors
C = lines(numel(conds));

% ---------- Tile 1: overall vs breath length ----------
ax = nexttile;
cla(ax); hold(ax,'on');

h = gobjects(numel(conds),1);
for ii = 1:numel(conds)
    idx = (taskUse == conds(ii)) & isfinite(overall) & isfinite(yLen);
    h(ii) = scatter(ax, overall(idx), yLen(idx), 28, 'filled', ...
        'MarkerFaceColor', C(ii,:), 'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'none');
end

xlabel(ax, 'overall HRV score', 'Interpreter','none');
ylabel(ax, 'breath length (s)', 'Interpreter','none');
legend(ax, h, cellstr(conds), 'Location','best', 'Box','off', 'Interpreter','none');

grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
hold(ax,'off');

% ---------- Tile 2: RSA vs breath length ----------
ax = nexttile;
cla(ax); hold(ax,'on');

h = gobjects(numel(conds),1);
for ii = 1:numel(conds)
    idx = (taskUse == conds(ii)) & isfinite(rsa01) & isfinite(yLen);
    h(ii) = scatter(ax, rsa01(idx), yLen(idx), 28, 'filled', ...
        'MarkerFaceColor', C(ii,:), 'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'none');
end

xlabel(ax, 'RSA', 'Interpreter','none');
ylabel(ax, 'breath length (s)', 'Interpreter','none');
legend(ax, h, cellstr(conds), 'Location','best', 'Box','off', 'Interpreter','none');

grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
hold(ax,'off');



taskUse = string(taskVec(mUse));
yLen    = bamp(mUse);


% colors
C = lines(numel(conds));

% ---------- Tile 1: overall vs breath length ----------
ax = nexttile;
cla(ax); hold(ax,'on');

h = gobjects(numel(conds),1);
for ii = 1:numel(conds)
    idx = (taskUse == conds(ii)) & isfinite(overall) & isfinite(yLen);
    h(ii) = scatter(ax, overall(idx), yLen(idx), 28, 'filled', ...
        'MarkerFaceColor', C(ii,:), 'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'none');
end

xlabel(ax, 'overall HRV score', 'Interpreter','none');
ylabel(ax, 'breath amp (au)', 'Interpreter','none');
legend(ax, h, cellstr(conds), 'Location','best', 'Box','off', 'Interpreter','none');

grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
hold(ax,'off');

% ---------- Tile 2: RSA vs breath length ----------
ax = nexttile;
cla(ax); hold(ax,'on');

h = gobjects(numel(conds),1);
for ii = 1:numel(conds)
    idx = (taskUse == conds(ii)) & isfinite(rsa01) & isfinite(yLen);
    h(ii) = scatter(ax, rsa01(idx), yLen(idx), 28, 'filled', ...
        'MarkerFaceColor', C(ii,:), 'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'none');
end

xlabel(ax, 'RSA', 'Interpreter','none');
ylabel(ax, 'breath amp (au)', 'Interpreter','none');
legend(ax, h, cellstr(conds), 'Location','best', 'Box','off', 'Interpreter','none');

grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
hold(ax,'off');




ax = nexttile; 

pacFrex = chanDat.breathLock.frex;
fb      = 1./chanDat.behDat.length(useVec);
bLock   = chanDat.breathLock.ispcZ_blk(useVec,:); 
fbi     = arrayfun(@(x) find(min(abs(x - pacFrex(:))) == ...
                        abs(x - pacFrex(:))), fb, 'uniformoutput', false);  
fbi     = cell2mat(fbi);
idx     = 1:length(fbi); 
bLock   = arrayfun(@(x,y) bLock(x, y), idx(:), fbi(:)); 

scatter(ax, overall, bLock)
xlabel(ax, 'overall HRV score', 'Interpreter','none');
ylabel(ax, 'breath - EEG ISPC', 'Interpreter','none');
grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
title(ax, chanDat.labels{chanDat.chi})

ax = nexttile; 

pacFrex = chanDat.breathLock.frex;
fb      = 1./chanDat.behDat.length(useVec);
bLock   = chanDat.breathLock.powObs(useVec,:); 
fbi     = arrayfun(@(x) find(min(abs(x - pacFrex(:))) == ...
                        abs(x - pacFrex(:))), fb, 'uniformoutput', false);  
fbi     = cell2mat(fbi);
idx     = 1:length(fbi); 

bLock   = arrayfun(@(x,y) bLock(x, y) ./ median(bLock(x,pacFrex<1)), idx(:), fbi(:)); 

scatter(ax, overall, bLock)
xlabel(ax, 'overall HRV score', 'Interpreter','none');
ylabel(ax, 'EEG pow at rsp freq (au)', 'Interpreter','none');
grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
title(ax, chanDat.labels{chanDat.chi})

ax = nexttile; 

phase      = chanDat.tf.phase(useVec,:,:); 
frex       = chanDat.tf.frex; 
freqBreaks = [2, 8; ...
              8, 14; ...
              14,30];
freqLabs = {'theta','alpha','beta'};

% colors: forest green, burnt orange, dark mauve
frexCols = [0.13 0.55 0.13;
            0.80 0.33 0.00;
            0.55 0.41 0.53];
splitidx = overall > median(overall);

hold(ax, 'on')
hLeg = gobjects(3,1);
itpcout = nan(3, 2, 50); 

for fi = 1:3
    fMask = frex>=freqBreaks(fi,1) & frex<freqBreaks(fi,2);

    itpc1  = squeeze(abs(mean(exp(1i * phase(splitidx,  :, fMask)), 1, 'omitnan')));
    itpc1  = max(itpc1, [], 2);

    itpc0  = squeeze(abs(mean(exp(1i * phase(~splitidx, :, fMask)), 1, 'omitnan')));
    itpc0  = max(itpc0, [], 2);

    hLeg(fi) = plot(ax, movmean(itpc1,3), 'LineWidth', 1.75, 'Color', frexCols(fi,:));
               plot(ax, movmean(itpc0,3), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', frexCols(fi,:));
    itpcout(fi, 1, :) = itpc1; 
    itpcout(fi, 2, :) = itpc0;
end

legend(ax, hLeg, freqLabs, 'Location','best', 'Box','off')
xlabel('normalized time in breath')
ylabel('ITPC')
grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
title(ax, chanDat.labels{chanDat.chi})


ax = nexttile; 

powZ      = chanDat.tf.powZ(useVec,:,:); 
frex       = chanDat.tf.frex; 
freqBreaks = [2, 8; ...
              8, 14; ...
              14,30];
freqLabs = {'theta','alpha','beta'};

% colors: forest green, burnt orange, dark mauve
frexCols = [0.13 0.55 0.13;
            0.80 0.33 0.00;
            0.55 0.41 0.53];
splitidx = overall > median(overall);

hold(ax, 'on')
hLeg = gobjects(3,1);
powout = nan(3, 2, 50);

for fi = 1:3
    fMask = frex>=freqBreaks(fi,1) & frex<freqBreaks(fi,2);

    pow1  = squeeze(mean(powZ(splitidx,  :, fMask), 1, 'omitnan'));
    pow1  = max(pow1, [], 2);

    pow2  = squeeze(mean(powZ(~splitidx,  :, fMask), 1, 'omitnan'));
    pow2  = max(pow2, [], 2);

    hLeg(fi) = plot(ax, movmean(pow1,3), 'LineWidth', 1.75, 'Color', frexCols(fi,:));
               plot(ax, movmean(pow2,3), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', frexCols(fi,:));

    powout(fi, 1, :) = pow1; 
    powout(fi, 2, :) = pow2; 

end

legend(ax, hLeg, freqLabs, 'Location','best', 'Box','off')
xlabel('normalized time in breath')
ylabel('power (z-scored)')
grid(ax,'on'); ax.GridAlpha = 0.15;
set(ax,'LineWidth',1.5,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out','FontSize',11);
title(ax, chanDat.labels{chanDat.chi})










end







% ---------------- local helper ----------------
function y = local_minmax01(v)
    v = v(:);
    mn = min(v(isfinite(v)));
    mx = max(v(isfinite(v)));
    if isempty(mn) || isempty(mx) || mx <= mn
        y = nan(size(v));
        if any(isfinite(v)), y(isfinite(v)) = 0; end
        return
    end
    y = (v - mn) ./ (mx - mn);
end