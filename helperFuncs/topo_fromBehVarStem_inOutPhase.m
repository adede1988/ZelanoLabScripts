function hFig = topo_fromBehVarStem_inOutPhase(allSubIDs, allBehDat, idxVec, varStem, targetidx50)
% topo_fromBehVarStem_inOutPhase(allSubIDs, allBehDat, idxVec, varStem, targetidx50)
% One-off helper for this pipeline.
%
% Plots two topo panels:
%   1) trials IN phase with targetidx50 (within 5 circular positions on 1..50 scale)
%   2) trials OUT of phase with targetidx50
%
% KEY FIX retained from original:
% aligns channels by matching eegLocs.Label <-> allSubIDs(:,4) WITHIN EACH SESSION,
% so channel ordering differences across subjects/sessions do not matter.
%
% Assumes:
%   - eegLocs exists in base workspace
%   - miniTopo(m,x,y) is on path
%   - each T has field gamRspPACphase_idx50
%
% Circular phase-bin rule:
%   dist = min(abs(idx50-targetidx50), 50-abs(idx50-targetidx50))
%   inPhase = dist <= 5

% pull eegLocs from base workspace (pipeline assumption)
eegLocs   = evalin('base','eegLocs');
eegLabels = strtrim(string(eegLocs.Label(:)));
x         = eegLocs.X2D_right(:);
y         = eegLocs.Y2D_front(:);
nLoc      = numel(eegLabels);

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
    if isnumeric(v)
        sid(i) = double(v);
    else
        sid(i) = str2double(string(v));
    end
end

uSID = unique(sid(mUse), 'stable');
uSID = uSID(isfinite(uSID));
nSess = numel(uSID);

% per-session topo vectors in eegLocs order
sessTopo_in  = nan(nLoc, nSess);
sessTopo_out = nan(nLoc, nSess);

chanLabAll = strtrim(string(allSubIDs(:,4)));

tol = 0.25; % radians tolerance for "is [-pi,+pi]?"

% ---- build each session topo in eegLocs label order ----
for s = 1:nSess
    mSess = mUse & (sid == uSID(s));
    rowsSess = find(mSess);

    % labels present for THIS session (may be in any order)
    chanLabSess = chanLabAll(rowsSess);

    % map eegLocs labels -> row index for this session
    [tf, posInSess] = ismember(eegLabels, chanLabSess);

    for li = 1:nLoc
        if ~tf(li), continue; end

        r = rowsSess(posInSess(li));   % absolute row index for this session+electrode
        T = allBehDat{r};

        % variable suffix: last 2 chars of electrode label
        EE = char(eegLabels(li));
        if numel(EE) > 2
            EE = EE(end-1:end);
        end

        fullVar = char(string(varStem) + string(EE));

        if ~ismember(fullVar, T.Properties.VariableNames)
            continue
        end
        if ~ismember('gamRspPACphase_idx50', T.Properties.VariableNames)
            continue
        end

        vdat   = T.(fullVar);
        idx50  = T.gamRspPACphase_idx50;

        % apply useVec if present
        if ismember('useVec', T.Properties.VariableNames)
            uv = T.useVec == 1;
            vdat  = vdat(uv);
            idx50 = idx50(uv);
        end

        vdat  = double(vdat(:));
        idx50 = double(idx50(:));

        % keep aligned finite entries
        good = isfinite(vdat) & isfinite(idx50);
        vdat  = vdat(good);
        idx50 = idx50(good);

        if isempty(vdat), continue; end

        % wrap idx50 to 1..50 just in case
        idx50 = mod(idx50 - 1, 50) + 1;
        tgt   = mod(targetidx50 - 1, 50) + 1;

        % circular distance on 1..50 ring
        d = abs(idx50 - tgt);
        d = min(d, 50 - d);

        inPhase  = d <= 5;
        outPhase = d > 5;

        if any(inPhase)
            sessTopo_in(li,s) = meanMaybeCirc(vdat(inPhase), tol);
        end
        if any(outPhase)
            sessTopo_out(li,s) = meanMaybeCirc(vdat(outPhase), tol);
        end
    end
end

% ---- group topo (average across sessions) ----
groupTopo_in  = nan(nLoc,1);
groupTopo_out = nan(nLoc,1);

for li = 1:nLoc
    groupTopo_in(li)  = meanMaybeCirc(sessTopo_in(li,:),  tol);
    groupTopo_out(li) = meanMaybeCirc(sessTopo_out(li,:), tol);
end

% ---- title text ----
repRows = nan(nSess,1);
for s = 1:nSess
    repRows(s) = find(mUse & sid==uSID(s), 1, 'first');
end

nameStr = strings(nSess,1);
for s = 1:nSess
    nm = string(allSubIDs{repRows(s),8});
    ss = allSubIDs{repRows(s),5};
    if isnumeric(ss)
        ssn = double(ss);
    else
        ssn = str2double(string(ss));
    end
    nameStr(s) = strtrim(nm) + "(sess" + string(ssn) + ")";
end

titleStr = string(varStem) + " | targetidx50=" + string(targetidx50) + ...
    " | " + strjoin(nameStr, ", ");

% ---- use matched color scaling when variable is not circular radians ----
allVals = [groupTopo_in(:); groupTopo_out(:)];
allVals = allVals(isfinite(allVals));
useMatchedCLim = false;
clims = [];

if ~isempty(allVals)
    isCircAll = max(abs(allVals)) <= (pi + tol);
    if ~isCircAll
        useMatchedCLim = true;
        mx = max(abs(allVals));
        if isfinite(mx) && mx > 0
            clims = [-mx mx];
        end
    end
end

% ---- plot ----
hFig = figure('Color','w');

subplot(1,2,1)
miniTopo(groupTopo_in, x, y);
title(sprintf('In phase (within 5 of %d)', targetidx50), 'Interpreter','none');
if useMatchedCLim && ~isempty(clims)
    clim(clims)
end
colorbar
subplot(1,2,2)
miniTopo(groupTopo_out, x, y);
title(sprintf('Out of phase'), 'Interpreter','none');
if useMatchedCLim && ~isempty(clims)
    clim(clims)
end
colorbar
sgtitle(char(titleStr), 'Interpreter','none');

end


function m = meanMaybeCirc(v, tol)
% average helper:
%   - circular mean if values look like radians in [-pi,+pi]
%   - arithmetic mean otherwise

v = double(v(:));
v = v(isfinite(v));

if isempty(v)
    m = NaN;
    return
end

isCirc = max(abs(v)) <= (pi + tol);
if isCirc
    m = angle(mean(exp(1i*v)));
else
    m = mean(v, 'omitnan');
end
end