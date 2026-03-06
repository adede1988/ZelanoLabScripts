function hFig = topo_fromBehVarStem(allSubIDs, allBehDat, idxVec, varStem)
% topo_fromBehVarStem(allSubIDs, allBehDat, idxVec, varStem)
% One-off helper for this pipeline.
%
% KEY FIX: aligns channels by matching eegLocs.Label <-> allSubIDs(:,4) WITHIN EACH SESSION,
% so channel ordering differences across subjects/sessions do not matter.
%
% Assumes eegLocs exists in base workspace and miniTopo(m,x,y) is on path.

% pull eegLocs from base workspace (pipeline assumption)
eegLocs = evalin('base','eegLocs');
eegLabels = strtrim(string(eegLocs.Label(:)));
x = eegLocs.X2D_right(:);
y = eegLocs.Y2D_front(:);
nLoc = numel(eegLabels);

% index mask
nRow = size(allSubIDs,1);
if islogical(idxVec)
    mUse = idxVec(:);
else
    mUse = false(nRow,1);
    mUse(idxVec(:)) = true;
end

% subject/session ID (col 7)
sid = nan(nRow,1);
for i = 1:nRow
    v = allSubIDs{i,7};
    if isnumeric(v); sid(i) = double(v);
    else;            sid(i) = str2double(string(v));
    end
end

uSID = unique(sid(mUse), 'stable');
uSID = uSID(isfinite(uSID));
nSess = numel(uSID);

% per-session topo vectors in eegLocs order
sessTopo = nan(nLoc, nSess);

chanLabAll = strtrim(string(allSubIDs(:,4)));

% ---- build each session topo in eegLocs label order ----
tol = 0.25; % radians tolerance for "is [-pi,+pi]?"
for s = 1:nSess
    mSess = mUse & (sid == uSID(s));
    rowsSess = find(mSess);

    % labels present for THIS session (may be in any order)
    chanLabSess = chanLabAll(rowsSess);

    % map eegLocs labels -> row index for this session
    [tf, posInSess] = ismember(eegLabels, chanLabSess);  % tf/li indicates if label exists
    % tf(li)==true => matching row is rowsSess(posInSess(li))

    for li = 1:nLoc
        if ~tf(li), continue; end

        r = rowsSess(posInSess(li));     % absolute row index for this session+electrode
        T = allBehDat{r};

        % variable suffix: last 2 chars of electrode label
        EE = char(eegLabels(li));
        if numel(EE) > 2
            EE = EE(end-1:end);
        end

        fullVar = char(string(varStem) + string(EE));

        if ismember(fullVar, T.Properties.VariableNames)
            vdat = T.(fullVar);

            if ismember('useVec', T.Properties.VariableNames)
                vdat = vdat(T.useVec == 1);
            end

            v = double(vdat(:));
            v = v(isfinite(v));
            if isempty(v), continue; end

            % detect circular radians in [-pi,+pi] and average appropriately
            isCirc = max(abs(v)) <= (pi + tol);
            if isCirc
                sessTopo(li,s) = angle(mean(exp(1i*v)));
            else
                sessTopo(li,s) = mean(v, 'omitnan');
            end
        end
    end
end

% group topo (average across sessions)
groupTopo = mean(sessTopo, 2, 'omitnan');

% ---- title: varStem + "Name(sessN), Name(sessN), ..." for included sessions ----
repRows = nan(nSess,1);
for s = 1:nSess
    repRows(s) = find(mUse & sid==uSID(s), 1, 'first');
end

nameStr = strings(nSess,1);
for s = 1:nSess
    nm = string(allSubIDs{repRows(s),8});
    ss = allSubIDs{repRows(s),5};
    if isnumeric(ss); ssn = double(ss); else; ssn = str2double(string(ss)); end
    nameStr(s) = strtrim(nm) + "(sess" + string(ssn) + ")";
end

titleStr = string(varStem) + " | " + strjoin(nameStr, ", ");

% ---- plot ----
hFig = figure('Color','w');
miniTopo(groupTopo, x, y);
title(char(titleStr), 'Interpreter','none');

end