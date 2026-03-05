function plotERP(varargin)
% plotERP(data1, data2, ..., eegLocs, channelInfo, typeSel, sessSel, [subSel], erpLabels)
%
% Optional subSel = participant/session ID (channelInfo col 7).
%
% Examples:
%   plotERP(data1, data2, eegLocs, channelInfo, "OBE", 1, {'peaks','controlPoints'})
%   plotERP(data1, data2, eegLocs, channelInfo, "Dupi", 1, 23, {'peaks','controlPoints'})
%   plotERP(data1, [],   eegLocs, channelInfo, [], [], {'peaks'})

% ---------------- locate eegLocs (table) and parse tail args robustly ----------------
iLoc = find(cellfun(@istable, varargin), 1, 'last');
if isempty(iLoc), error('Could not find eegLocs table in inputs.'); end

eegLocs     = varargin{iLoc};
channelInfo = varargin{iLoc+1};
typeSel     = varargin{iLoc+2};
sessSel     = varargin{iLoc+3};

if numel(varargin) == iLoc+4
    subSel    = [];
    erpLabels = varargin{iLoc+4};
elseif numel(varargin) == iLoc+5
    subSel    = varargin{iLoc+4};
    erpLabels = varargin{iLoc+5};
else
    error('Input parsing failed. Expected [data..., eegLocs, channelInfo, typeSel, sessSel, (subSel), erpLabels].');
end

dataArgs = varargin(1:iLoc-1);

% keep only non-empty ERP matrices
keep = ~cellfun(@isempty, dataArgs);
dataArgs = dataArgs(keep);
nERP = numel(dataArgs);

% labels
if isstring(erpLabels); erpLabels = cellstr(erpLabels); end
if ~iscell(erpLabels), error('erpLabels must be a cell array (or string array).'); end
if numel(erpLabels) ~= nERP
    error('erpLabels length (%d) must match number of non-empty data matrices (%d).', numel(erpLabels), nERP);
end

% ---------------- basic sizes ----------------
if ~iscell(channelInfo) || size(channelInfo,2) < 7
    error('channelInfo must be channels x 7 cell array.');
end
nChan = size(channelInfo,1);

nTime = [];
for k = 1:nERP
    D = dataArgs{k};
    if ~isnumeric(D) || ndims(D)~=2, error('Each data matrix must be numeric 2D: channels x time.'); end
    if size(D,1) ~= nChan, error('Data matrix %d has %d channels; expected %d.', k, size(D,1), nChan); end
    if isempty(nTime), nTime = size(D,2); end
    if size(D,2) ~= nTime, error('All data matrices must have same timepoints.'); end
end

% ---------------- eegLocs ----------------
locLabels = strtrim(string(eegLocs.Label(:)));
x = eegLocs.X2D_right(:);
y = eegLocs.Y2D_front(:);
nLoc = numel(locLabels);

xN = (x - min(x)) ./ max(eps, (max(x)-min(x)));
yN = (y - min(y)) ./ max(eps, (max(y)-min(y)));

% ---------------- channelInfo fields ----------------
chanLabels = strtrim(string(channelInfo(:,4)));
ptype      = lower(strtrim(string(channelInfo(:,6))));

% session number (col 5)
sessNum = nan(nChan,1);
for i = 1:nChan
    v = channelInfo{i,5};
    if isnumeric(v); sessNum(i) = double(v);
    else;            sessNum(i) = str2double(string(v));
    end
end

% subject/session group ID (col 7)
subID = nan(nChan,1);
for i = 1:nChan
    v = channelInfo{i,7};
    if isnumeric(v); subID(i) = double(v);
    else;            subID(i) = str2double(string(v));
    end
end

% ---------------- selection mask ----------------
mask = true(nChan,1);

if ~isempty(typeSel) && ~(isstring(typeSel) && strlength(typeSel)==0)
    t = lower(strtrim(string(typeSel)));
    if t == "dupi", t = "dupi"; end
    mask = mask & (ptype == t);
end

if ~isempty(sessSel)
    mask = mask & (sessNum == double(sessSel));
end

if ~isempty(subSel)
    mask = mask & (subID == double(subSel));
end

% ---------------- plotting params ----------------
tVec = -4000:2:4000;            % replace with your real time vector if desired
smoothWin = 25;            % 25-point moving average
vlineX    = 0;          % vertical reference line index

% make mini axes larger
scale = sqrt(32 / max(1,nLoc));
axW = min(0.18, max(0.05, 0.13 * scale));   % larger than before
axH = min(0.14, max(0.04, 0.10 * scale));   % larger than before
marg = 0.04;                                 % slightly smaller margin to fit

% ---------------- create clean figure (no background axes) ----------------
figure('Color','w');
clf;  % important: nukes any default axes

legH = gobjects(nERP,1);
axLegendParent = [];  % use a real mini-axes as legend parent
gotLegendHandles = false;

for li = 1:nLoc
    lab = locLabels(li);
    mChan = mask & (chanLabels == lab);

    left   = marg + xN(li)*(1-2*marg) - axW/2;
    bottom = marg + yN(li)*(1-2*marg) - axH/2;
    left   = min(max(left,   0.001), 0.999-axW);
    bottom = min(max(bottom, 0.001), 0.999-axH);

    ax = axes('Units','normalized', 'Position',[left bottom axW axH]); %#ok<LAXES>
    hold(ax,'on');

    if any(mChan)
        if isempty(axLegendParent), axLegendParent = ax; end

        for k = 1:nERP
            mu = mean(dataArgs{k}(mChan,:), 1, 'omitnan');

            % 25-pt moving average smoothing (fallback if movmean unavailable)
            mu = local_smooth1d(mu, smoothWin);
            
            h = plot(ax, tVec, mu, 'LineWidth', 1.5);

            if ~gotLegendHandles
                legH(k) = h;
            end
        end
        gotLegendHandles = gotLegendHandles || all(isgraphics(legH));

        % reference lines
        yline(ax, 0, ':');
        xline(ax, vlineX, 'k:');

        axis(ax,'tight');

        % remove boxes/frames/ticks
        box(ax,'off');
        set(ax, 'XTick', [], 'YTick', [], ...
                'XColor','none', 'YColor','none', ...
                'Color','none');   % transparent background
        ylim(ax, [-7 7])
        xlim(ax, [-1000 4000])
        title(ax, char(lab), 'Interpreter','none', 'FontSize', 7);

    else
        axis(ax,'off');
    end
end


%% ---- Reference axes with labeled units (ms, µV) ----
% Drop this in AFTER you’ve plotted all inset ERPs.

% 1) Find an example inset axis that has ERP lines on it
axAll = findall(gcf, 'Type','axes');
axEx  = [];

for a = axAll(:)'
    % pick an axis that is visible and contains at least one line
    if strcmp(get(a,'Visible'),'on') && ~isempty(findobj(a,'Type','line'))
        axEx = a;
        break
    end
end

if ~isempty(axEx) && isgraphics(axEx)
    xl = get(axEx,'XLim');
    yl = get(axEx,'YLim');

    % 2) Create a new small axis (bottom-left); adjust Position as desired
    axRef = axes('Units','normalized', 'Position',[0.04 0.02 axW axH], ...
                 'Color','none', 'Box','off', 'TickDir','out', ...
                 'LineWidth',1);

    xlim(axRef, xl);
    ylim(axRef, yl);

    % 3) Sensible ticks (min/0/max if 0 is inside range; else just min/max)
    xt = [xl(1) xl(2)];
    yt = [yl(1) yl(2)];
    if xl(1) < 0 && xl(2) > 0, xt = [xl(1) 0 xl(2)]; end
    if yl(1) < 0 && yl(2) > 0, yt = [yl(1) 0 yl(2)]; end
    set(axRef, 'XTick', xt, 'YTick', yt);

    % 4) Labels
    xlabel(axRef, 'Time (ms)');
    ylabel(axRef, 'Voltage (\muV)');

    % (Optional) draw baseline on the reference axis
    hold(axRef,'on');
    yline(axRef, 0, ':');
    % If your time-zero is at x=0 ms and you want it marked:
    xline(axRef, 0, ':');
end


% ---- legend attached to the reference axis ----
if gotLegendHandles && exist('axRef','var') && isgraphics(axRef)
    lg = legend(axRef, legH, erpLabels, 'Interpreter','none', ...
        'Box','off', 'Orientation','horizontal', 'Location','southoutside');
    lg.ItemTokenSize = [10 10];
end

% figure-level title without sgtitle (prevents full-figure hidden axes)
selTxt = "";
if ~isempty(typeSel), selTxt = selTxt + string(typeSel); end
if ~isempty(sessSel), selTxt = selTxt + " | sess " + string(sessSel); end
if ~isempty(subSel),  selTxt = selTxt + " | subID " + string(subSel); end
if strlength(selTxt)==0, selTxt = "all participants/sessions"; end

annotation(gcf,'textbox',[0 0.965 1 0.03], ...
    'String', char(selTxt), ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'Interpreter','none', ...
    'FontWeight','bold');

end

% -------- local helper: smooth row vector with moving average ----------
function y = local_smooth1d(x, win)
x = x(:)'; % row
try
    y = movmean(x, win, 2, 'Endpoints','shrink');
catch
    k = ones(1,win) / win;
    y = conv(x, k, 'same');
end
end