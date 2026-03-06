function hFig = plot_dumbbell_sess1_sess2(allSubIDs, allBehDat, xVarName)
% One-off helper for this pipeline.
% Dumbbell plot of Session 1 -> Session 2 changes:
%   Y-axis: HRV index (always)
%   X-axis: mean of xVarName within session (circular mean if values in [-pi,+pi])
%
% Only includes breaths where useVec==1 inside each allBehDat{row} table.
% Automatically includes only participants with BOTH session 1 and session 2.
% allSubIDs has redundant rows per subject/session; this function uses only 1 row per subject/session.
%
% Markers:
%   Session 1 = square
%   Session 2 = triangle
% Each participant has a unique color.

% ---- pull subject name + session number ----
nRow = size(allSubIDs,1);

subName = strtrim(string(allSubIDs(:,8)));

sessNum = nan(nRow,1);
for i = 1:nRow
    v = allSubIDs{i,5};
    if isnumeric(v); sessNum(i) = double(v);
    else;            sessNum(i) = str2double(string(v));
    end
end

% ---- select ONE representative row per (subject, session) ----
key = subName + "_sess" + string(sessNum);
[~, repRows] = unique(key, 'stable');   % first occurrence for each subject/session

repSub  = subName(repRows);
repSess = sessNum(repRows);

% ---- participants (by name) that have both sessions ----
uSub = unique(repSub, 'stable');

x1 = []; y1 = []; x2 = []; y2 = []; subKeep = strings(0,1);

for s = 1:numel(uSub)
    r1 = repRows(repSub==uSub(s) & repSess==1);
    r2 = repRows(repSub==uSub(s) & repSess==2);
    if isempty(r1) || isempty(r2)
        continue
    end
    r1 = r1(1); r2 = r2(1);

    % ---- session 1 ----
    T1 = allBehDat{r1};
    m1 = (T1.useVec == 1);

    xraw1 = double(T1.(xVarName)(m1));
    hrv1  = 1000*double(T1.HRV_RMSSD30(m1)) + 100*double(T1.HRV_SDNN30(m1)) + double(T1.HRV_RSAamp(m1));

    % ---- session 2 ----
    T2 = allBehDat{r2};
    m2 = (T2.useVec == 1);

    xraw2 = double(T2.(xVarName)(m2));
    hrv2  = 1000*double(T2.HRV_RMSSD30(m2)) + 100*double(T2.HRV_SDNN30(m2)) + double(T2.HRV_RSAamp(m2));

    % ---- clean finite ----
    xraw1 = xraw1(isfinite(xraw1));
    xraw2 = xraw2(isfinite(xraw2));

    hrv1  = hrv1(isfinite(hrv1));
    hrv2  = hrv2(isfinite(hrv2));

    if isempty(xraw1) || isempty(xraw2) || isempty(hrv1) || isempty(hrv2)
        continue
    end

    % ---- circular vs linear mean on X ----
    % (rule: if values lie in ~[-pi,+pi], treat as circular)
    isCirc1 = max(abs(xraw1)) <= (pi + 0.25);
    isCirc2 = max(abs(xraw2)) <= (pi + 0.25);

    if isCirc1
        xmu1 = angle(mean(exp(1i*xraw1), 'omitnan'));
    else
        xmu1 = mean(xraw1, 'omitnan');
    end

    if isCirc2
        xmu2 = angle(mean(exp(1i*xraw2), 'omitnan'));
    else
        xmu2 = mean(xraw2, 'omitnan');
    end

    % ---- HRV mean on Y ----
    ymu1 = mean(hrv1, 'omitnan');
    ymu2 = mean(hrv2, 'omitnan');

    x1(end+1,1) = xmu1; %#ok<AGROW>
    y1(end+1,1) = ymu1; %#ok<AGROW>
    x2(end+1,1) = xmu2; %#ok<AGROW>
    y2(end+1,1) = ymu2; %#ok<AGROW>
    subKeep(end+1,1) = uSub(s); %#ok<AGROW>
end

% ---- plot ----
hFig = figure('Color','w');
ax = axes('Parent', hFig); hold(ax,'on');

nP = numel(subKeep);
cols = lines(max(nP,1));

for i = 1:nP
    plot(ax, [x1(i) x2(i)], [y1(i) y2(i)], '-', 'Color', cols(i,:), 'LineWidth', 1.5);

 
end

legend(ax, cellstr(subKeep), 'Interpreter','none', 'Box','off', 'Location','eastoutside', 'AutoUpdate','off');


for i = 1:nP
  

    scatter(ax, x1(i), y1(i), 60, 's', ...
        'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', cols(i,:));

    scatter(ax, x2(i), y2(i), 70, '^', ...
        'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', cols(i,:));
end

xlabel(ax, xVarName, 'Interpreter','none');
ylabel(ax, 'HRV index = 1000*HRV\_RMSSD30 + 100*HRV\_SDNN30 + HRV\_RSAamp');
title(ax, sprintf('Session 1 \x2192 2 dumbbells: %s (X) vs HRV index (Y)', xVarName), 'Interpreter','none');

box(ax,'off');
grid(ax,'on');

% (optional) legend; comment out if too many participants
if nP > 0
    
end

end