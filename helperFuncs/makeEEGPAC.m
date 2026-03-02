function fig = makeEEGPAC(chanDat, macChan, taskVec, conds, opts)


x = 5

useVec = chanDat.use == 1; 

fig = figure('Color','w', 'visible', true, 'position', [0,0,1200, 700]);
tlo = tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');
pacFrex = chanDat.pac.PACfrex;
ax = nexttile(tlo, 1);

cols = lines(numel(conds));
cla(ax); hold(ax,'on');
maxi = 0; 
maxVal = 0; 
for c = 1:numel(conds)
    m = useVec & (taskVec == conds(c));
    if nnz(m) > 5
        curpac = squeeze(mean(median(chanDat.pac.pac(m,1:50,:,10), 2, 'omitnan'), 1, 'omitnan'));
        [curMax, curi] = max(curpac); 
        if curMax > maxVal
            maxVal = curMax; 
            maxi = curi; 
        end
        plot(ax, pacFrex, curpac, 'Color', cols(c,:), 'LineWidth', 2);
    end
end
xline(pacFrex(maxi))
% ---- styling ----
title(ax, 'OB \gamma amp - scalp phase PAC', 'FontWeight','bold');
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'z-scored PAC');
legend(ax, cellstr(conds), 'Location','southeast');
grid(ax,'on');
ax.GridAlpha = 0.15;          % light grid
ax.LineWidth = 1.5;           % bolder axes
ax.FontSize  = 11;
box(ax,'off');

hold(ax,'off');

%% phase preference at peak freq: 
ax = nexttile(tlo,2);

% --- condition info (same as prior) ---
taskVec = string(chanDat.behDat.task);
conds   = ["audio","focus","shadow"];
cols    = lines(numel(conds));   % matches your scatter colors

% --- data ---
phi      = chanDat.pac.pac(:,1:50, maxi, 5);

nRow  = height(chanDat.behDat);

m0 = useVec; 


% --- binning (frequency polygon) ---
nBin  = 12;
edges = linspace(0, 2*pi, nBin+1);
cent  = (edges(1:end-1) + edges(2:end))/2;

phiW = mod(phi, 2*pi);

% --- plot ---

pos = ax.Position;
delete(ax);

ax = polaraxes('Position', pos);
hold(ax,'on');
conds = ["audio","focus","shadow"];

for c = 1:numel(conds)
    mc = m0 & (taskVec == conds(c));
    if nnz(mc) < 5, continue; end
    curVals = phiW(mc,:); 
    curVals = curVals(:); 
    cnt = histcounts(curVals, edges);
    r   = cnt / sum(cnt);                 % normalize to probability (area-free)
    th  = [cent, cent(1)+2*pi];           % close polygon
    rr  = [r,   r(1)];

    polarplot(ax, th, rr, 'LineWidth', 2.0, 'color', cols(c,:));
end
rlim([0 .2])
title(ax, ['PAC phase at ' num2str(round(pacFrex(maxi))) ' Hz']);

hold(ax,'off');


%% condition tf plots of PAC strength


% --- condition info (same as prior) ---
taskVec = string(chanDat.behDat.task);
conds   = ["audio","focus","shadow"];
cols    = lines(numel(conds));   % matches your scatter colors

for c = 1:numel(conds)
    ax = nexttile(tlo, c+3); 
    mc = m0 & (taskVec == conds(c));
    if nnz(mc) < 5, continue; end
    curVals = mean(chanDat.pac.pac(mc, 1:50, :, 10), 1, 'omitnan'); 
    imagesc(ax, squeeze(curVals)')
    axis(ax,'xy')
    yticks(ax, 5:5:50)
    yticklabels(ax, round(pacFrex(5:5:50)))
  
    xline([10, 20, 30, 40])
    title(ax, sprintf('%s (n=%d)', conds(c), nnz(mc)));
   
    xlabel(ax,'Time'); ylabel(ax,'Freq (Hz)');
    caxis(ax, [-1 1]);

    cb = colorbar(ax);
    cb.Label.String = 'PAC z';
end






%% what's the breath by breath change in coupling strength? 

ax = nexttile(tlo,3);

% --- condition info (same as prior) ---
taskVec = string(chanDat.behDat.task);
conds   = ["audio","focus","shadow"];
cols    = lines(numel(conds));   % matches your scatter colors

% --- data ---
pac      = chanDat.pac.pac(:,1:50, maxi, 10);
pac      = prctile(pac, 80, 2);
nRow  = height(chanDat.behDat);

m0 = useVec; 




% --- plot ---
hold(ax,'on');
conds = ["audio","focus","shadow"];

for c = 1:numel(conds)
    mc = m0 & (taskVec == conds(c));
    if nnz(mc) < 5, continue; end
    curVals = pac(mc); 
  

    plot(ax, 1:numel(curVals), movmean(curVals, 3), 'LineWidth', 2.0, 'color', cols(c,:));
end

title(ax, ['Single breath PAC strength at ' num2str(round(pacFrex(maxi))) ' Hz']);

hold(ax,'off');


% %% pac at gamma peak 
% ax = nexttile(tlo,3);
% 
% cols = lines(numel(conds));
% cla(ax); hold(ax,'on');
% 
% for c = 1:numel(conds)
%     m = useVec & (taskVec == conds(c));
%     if nnz(m) > 5
%         curpac = squeeze(mean(chanDat.pac.pac(m,51,:,10), 1, 'omitnan'));
%         plot(ax, pacFrex, curpac, 'Color', cols(c,:), 'LineWidth', 2);
%     end
% end
% 
% % ---- styling ----
% title(ax, 'PAC at \gamma peak', 'FontWeight','bold');
% xlabel(ax, 'Frequency (Hz)');
% ylabel(ax, 'z-scored PAC');
% 
% grid(ax,'on');
% ax.GridAlpha = 0.15;          % light grid
% ax.LineWidth = 1.5;           % bolder axes
% ax.FontSize  = 11;
% box(ax,'off');
% 
% hold(ax,'off');
% 
% 
% ax = nexttile(tlo,3);
% 
% cols = lines(numel(conds));
% cla(ax); hold(ax,'on');
% 
% for c = 1:numel(conds)
%     m = useVec & (taskVec == conds(c));
%     if nnz(m) > 5
%         curpac = squeeze(mean(chanDat.pac.pac(m,51,:,10), 1, 'omitnan'));
%         plot(ax, pacFrex, curpac, 'Color', cols(c,:), 'LineWidth', 2);
%     end
% end
% 
% % ---- styling ----
% title(ax, 'PAC at \gamma peak', 'FontWeight','bold');
% xlabel(ax, 'Frequency (Hz)');
% ylabel(ax, 'z-scored PAC');
% 
% grid(ax,'on');
% ax.GridAlpha = 0.15;          % light grid
% ax.LineWidth = 1.5;           % bolder axes
% ax.FontSize  = 11;
% box(ax,'off');
% 
% hold(ax,'off');












end