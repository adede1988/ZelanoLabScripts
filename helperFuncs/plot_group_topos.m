function plot_group_topos(allPowShuf, allSubIDs, eegLocs)
% plot_group_topos
% Creates 9 figures total:
%   (OBE, Dupi/DUPI sess1, Dupi/DUPI sess2) × (3 frequency bands)
% Each figure is a 5×3 grid: 5 time-epochs (rows) × 3 conditions (cols).
%
% Requires: miniTopo(m, x, y)

% -------------------- basic checks --------------------
assert(ndims(allPowShuf)==4, 'allPowShuf must be [nChan x nBand x nCond x nTime]');
[nChan, nBand, nCond, nTime] = size(allPowShuf);
assert(nBand==3 && nCond==3 && nTime==50, 'Expected [*,3,3,50] for allPowShuf');
assert(size(allSubIDs,1)==nChan, 'allSubIDs rows must match allPowShuf channels');

% -------------------- helper: robust omitnan mean --------------------
meanOmit = @(A,dim) mean(A, dim, 'omitnan');
try
    meanOmit(rand(3,3),2); 
catch
    % older MATLAB fallback
    meanOmit = @(A,dim) nanmean(A, dim); 
end

% -------------------- EEG locs & labels --------------------
eegLabels = strtrim(string(eegLocs.Label(:)));
x = eegLocs.X2D_right(:);
y = eegLocs.Y2D_front(:);
nLoc = numel(eegLabels);

% -------------------- parse metadata from allSubIDs --------------------
chanLabelsAll = strtrim(string(allSubIDs(:,4)));
typeAll       = strtrim(string(allSubIDs(:,6)));

% session number (col 5) and group integer (col 7) may be numeric-in-cell, string, etc.
sessNumAll  = nan(nChan,1);
groupIDAll  = nan(nChan,1);
for i = 1:nChan
    v5 = allSubIDs{i,5};
    v7 = allSubIDs{i,7};

    if isnumeric(v5); sessNumAll(i) = double(v5);
    else;             sessNumAll(i) = str2double(string(v5));
    end

    if isnumeric(v7); groupIDAll(i) = double(v7);
    else;             groupIDAll(i) = str2double(string(v7));
    end
end

% -------------------- define epochs --------------------
epRanges = {1:10, 11:20, 21:30, 31:40, 41:50};
epNames  = {'inhale rise','inhale fall','exhale rise','exhale fall','pause'};
nEpoch   = numel(epRanges);

% condition & band labels 
condNames = {'low HRV', 'high HRV', 'high - low'};
bandNames = {'theta', 'alpha', 'beta'}; 

% -------------------- participant/session groups via groupID (col 7) --------------------
uG = unique(groupIDAll(isfinite(groupIDAll)));
uG = uG(:);
nG = numel(uG);

gType = strings(nG,1);
gSess = nan(nG,1);
for g = 1:nG
    idx0 = find(groupIDAll==uG(g), 1, 'first');
    gType(g) = typeAll(idx0);
    gSess(g) = sessNumAll(idx0);
end

% normalize type for matching
gTypeLower = lower(strtrim(gType));

mOBE   = (gTypeLower == "obe");
mDupi1 = (gTypeLower == "dupi") & (gSess == 1);
mDupi2 = (gTypeLower == "dupi") & (gSess == 2);

groupDefs = { ...
    struct('name','OBE',        'mask', mOBE), ...
    struct('name','Dupi_sess1', 'mask', mDupi1), ...
    struct('name','Dupi_sess2', 'mask', mDupi2) ...
};

% -------------------- compute + plot --------------------
for gi = 1:numel(groupDefs)

    gname = groupDefs{gi}.name;
    gmask = groupDefs{gi}.mask;

    gIDs = uG(gmask);
    if isempty(gIDs)
        warning('No participant-sessions found for group "%s". Skipping.', gname);
        continue
    end

    nPart = numel(gIDs);

    % participant-level aligned data:
    %   [nLoc x nBand x nCond x nEpoch x nPart]
    partTopo = nan(nLoc, nBand, nCond, nEpoch, nPart);

    for p = 1:nPart
        rows = find(groupIDAll == gIDs(p));                 % channels for this participant-session
        partLabels = chanLabelsAll(rows);                   % labels in allPowShuf order

        % map participant channels -> eegLocs order
        [tf, locPos] = ismember(partLabels, eegLabels);

        % epoch means in participant channel order, then insert into aligned positions
        for e = 1:nEpoch
            tmp = meanOmit(allPowShuf(rows,:,:,epRanges{e}), 4); % [nChanPart x nBand x nCond]

            if any(tf)
                partTopo(locPos(tf),:,:,e,p) = tmp(tf,:,:);
            end
        end

        % (optional) warn on missing labels
        if any(~tf)
            % comment this out if it gets noisy
            warning('Group %s, participant-session %d: %d/%d channel labels not found in eegLocs.', ...
                gname, gIDs(p), sum(~tf), numel(tf));
        end
    end

    % group average across participant-sessions (omit NaNs from missing channels)
    groupTopo = meanOmit(partTopo, 5);  % [nLoc x nBand x nCond x nEpoch]

    % ---- plot 3 figures (one per band) for this group ----
    % for b = 1:nBand
    b = 1
        figure('Color','w', 'Name', sprintf('%s | %s', gname, bandNames{b}));
        t = tiledlayout(nEpoch, nCond, 'TileSpacing','compact', 'Padding','compact');
        title(t, sprintf('%s | %s', gname, bandNames{b}), 'Interpreter','none');

        for e = 1:nEpoch
            for c = 1:nCond
                nexttile;

                m = squeeze(groupTopo(:,b,c,e));   % [nLoc x 1]
               

                miniTopo(m, x, y);
                colorbar
                % clim([.4 .7])

                % top row: condition titles
                if e == 1
                    title(condNames{c}, 'Interpreter','none');
                end
                % left column: epoch labels (as y-labels so they sit outside the topo)
                if c == 1
                    ylabel(epNames{e}, 'Interpreter','none');
                   
                end
            end
        end
    % end
end
end