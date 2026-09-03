function assemble_spectrograms()
% ASSEMBLE_SPECTROGRAMS  Build the 5x5 group inhale-locked gamma spectrograms.
% Rows = {cueTask, threshTask, O15, audiobook, focusedBreathing}
% Cols = {control, S2/3 responder, S1 responder, S2/3 non-responder, S1 non-responder}
% Group map = mean over breaths of per-breath single-trial within-frequency z
%             (== myChanZscore, no baseline), combining all breaths in the cell.

here=fileparts(mfilename('fullpath')); proj=fileparts(here);
spdir=fullfile(proj,'out','gamma','spectro');
tdir =fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs');
if ~exist(fdir,'dir'), mkdir(fdir); end

resp = readtable(fullfile(tdir,'responder_table.csv'),'TextType','string');
classOf = containers.Map('KeyType','char','ValueType','char');
for i=1:height(resp), classOf(char(resp.participant(i))) = char(resp.class(i)); end

rows = {'cueTask','threshTask','O15','audiobook','focusedBreathing'};
cols = {'control','S2/3 responder','S1 responder','S2/3 non-responder','S1 non-responder'};
nR=numel(rows); nC=numel(cols);
Zsum=cell(nR,nC); Nsum=zeros(nR,nC); Ssum=zeros(nR,nC); tMs=[]; freqs=[];

files = dir(fullfile(spdir,'*.mat'));
for k=1:numel(files)
    S=load(fullfile(spdir,files(k).name));
    if ~isfield(S,'sumZwhole')||isempty(S.sumZwhole), continue; end
    m=S.meta; tr=char(m.taskRow);
    ri=find(strcmp(rows,tr),1); if isempty(ri), continue; end
    ci=assign_col(m, classOf);
    if isempty(ci), continue; end
    if isempty(tMs), tMs=m.tMs; freqs=m.freqs; end
    if isempty(Zsum{ri,ci}), Zsum{ri,ci}=zeros(size(S.sumZwhole)); end
    Zsum{ri,ci}=Zsum{ri,ci}+S.sumZwhole;
    Nsum(ri,ci)=Nsum(ri,ci)+m.nBreaths;
    Ssum(ri,ci)=Ssum(ri,ci)+1;
end

% build maps
maps=cell(nR,nC);
for ri=1:nR, for ci=1:nC
    if ~isempty(Zsum{ri,ci}) && Nsum(ri,ci)>0
        maps{ri,ci}=Zsum{ri,ci}/Nsum(ri,ci);
    end
end, end
save(fullfile(proj,'out','gamma','group_spectrograms.mat'),'maps','Nsum','Ssum','tMs','freqs','rows','cols','-v7');

% ---- plot ----
fig=figure('Position',[50 50 1700 1500],'Color','w','Visible','off');
tl=tiledlayout(nR,nC,'TileSpacing','compact','Padding','compact');
% common color scale from robust percentiles across all maps
allv=[]; for ri=1:nR, for ci=1:nC, if ~isempty(maps{ri,ci}), allv=[allv; maps{ri,ci}(:)]; end, end, end %#ok<AGROW>
cl = prctile(allv(isfinite(allv)),[2 98]); cl=max(abs(cl)); if isempty(cl)||cl==0, cl=1; end
for ri=1:nR
  for ci=1:nC
    nexttile;
    if ~isempty(maps{ri,ci})
        imagesc(tMs, freqs, maps{ri,ci}, [-cl cl]); axis xy; hold on;
        yline(25,'w:'); yline(58,'w:'); xline(0,'w-','LineWidth',1);
        title(sprintf('%s / %s\nn=%d br, %d sess', rows{ri}, cols{ci}, Nsum(ri,ci), Ssum(ri,ci)),'FontSize',7,'Interpreter','none');
    else
        axis off; title(sprintf('%s / %s\n(no data)', rows{ri}, cols{ci}),'FontSize',7,'Interpreter','none');
    end
    if ci==1, ylabel('Hz','FontSize',7); end
    if ri==nR, xlabel('ms','FontSize',7); end
    set(gca,'FontSize',6);
  end
end
try, colormap(turbo); catch, colormap(jet); end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='mean within-freq z (no baseline)';
title(tl,'Inhale/sniff-locked group gamma spectrograms (superlet; z within frequency, no baseline)','FontSize',11);
exportgraphics(fig, fullfile(fdir,'spectrograms_5x5.png'),'Resolution',140);
close(fig);
fprintf('wrote spectrograms_5x5.png ; cells with data: %d/%d\n', sum(Nsum(:)>0), nR*nC);
end

function ci = assign_col(m, classOf)
ci=[];
coh=''; if isfield(m,'cohort'), coh=char(m.cohort); end
part=''; if isfield(m,'participant'), part=char(m.participant); end
sn=1; if isfield(m,'sessNum'), sn=m.sessNum; end
if strcmpi(coh,'OBE'), ci=1; return; end          % control
if ~strcmpi(coh,'Dupi'), return; end
cl=''; if isKey(classOf, part), cl=classOf(part); end
isResp = strcmpi(cl,'responder');
isNon  = strcmpi(cl,'non-responder');
if ~isResp && ~isNon, return; end                  % unclassified -> excluded
if isResp
    if ismember(sn,[2 3]), ci=2; else, ci=3; end
else
    if ismember(sn,[2 3]), ci=4; else, ci=5; end
end
end
