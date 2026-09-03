function assemble_spectrograms()
% ASSEMBLE_SPECTROGRAMS  Build 3 sets of 5x5 group spectrograms (freq ranges
% low 2-25, mid 25-58, high 60-150 Hz), z-scored to the -500..-100 ms pre-inhale
% baseline, with the mean respiration waveform overlaid. One shared color scale
% per frequency-range set. Reads out/spectro2/*.mat + responder_table.csv.

here=fileparts(mfilename('fullpath')); proj=fileparts(here);
sdir=fullfile(proj,'out','spectro2'); tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs');
if ~exist(fdir,'dir'), mkdir(fdir); end

resp=readtable(fullfile(tdir,'responder_table.csv'),'TextType','string');
classOf=containers.Map('KeyType','char','ValueType','char');
for i=1:height(resp), classOf(char(resp.participant(i)))=char(resp.class(i)); end

rows={'cueTask','threshTask','O15','audiobook','focusedBreathing'};
cols={'control','S2/3 responder','S1 responder','S2/3 non-responder','S1 non-responder'};
sets={'low','mid','high'}; setLabel={'2-25 Hz','25-58 Hz','60-150 Hz'};
nR=numel(rows); nC=numel(cols);
CLIPZ=10;   % winsorize each session's z map to +-CLIPZ before averaging (suppress single-session broadband artifacts)

files=dir(fullfile(sdir,'*.mat'));
% Per the analysis convention: z-score WITHIN CHANNEL (per session) first, then
% AVERAGE the z-scored maps ACROSS SUBJECTS (sessions). So we build each session's
% within-channel baseline-z map from its sufficient stats and accumulate the sum.
Zsum=struct(); for s=1:3, Zsum.(sets{s})=cell(nR,nC); end
Ns=zeros(nR,nC); Ss=zeros(nR,nC); Rsp=cell(nR,nC); Nrsp=zeros(nR,nC); Subj=cell(nR,nC);
tMs=[]; freqs=struct();

for k=1:numel(files)
    S=load(fullfile(sdir,files(k).name)); if ~isfield(S,'out'), continue; end
    o=S.out; if o.nBreaths==0, continue; end
    ri=find(strcmp(rows,char(o.taskRow)),1); if isempty(ri), continue; end
    ci=assign_col(o, classOf); if isempty(ci), continue; end
    if isempty(tMs), tMs=o.tMs; for s=1:3, freqs.(sets{s})=o.(sets{s}).freqs; end, end
    nb = o.nBreaths * o.nBaseSamp;
    for s=1:3
        meanPow = o.(sets{s}).sumPow / o.nBreaths;           % session mean power
        baseMu  = o.(sets{s}).baseSum / nb;                   % baseline mean over (breaths x baseline-times)
        baseVar = o.(sets{s}).baseSumSq/nb - baseMu.^2; baseVar(baseVar<=0)=eps;
        % myChanZscore convention: baseline pooled over this session's breaths;
        % denominator = std(pooled baseline)/sqrt(nBreaths) (bootstrap-of-means / SEM).
        % Analytic equivalent of myChanZscore(power[time x breaths], baselinePeriod)
        % averaged over breaths -> avoids single-breath dominance.
        zsess   = (meanPow - baseMu) ./ sqrt(baseVar / o.nBreaths);
        zsess   = max(min(zsess, CLIPZ), -CLIPZ);            % winsorize: cap single-session broadband artifacts
        if isempty(Zsum.(sets{s}){ri,ci}), Zsum.(sets{s}){ri,ci}=zeros(size(zsess)); end
        Zsum.(sets{s}){ri,ci}=Zsum.(sets{s}){ri,ci}+zsess;
    end
    Ns(ri,ci)=Ns(ri,ci)+o.nBreaths; Ss(ri,ci)=Ss(ri,ci)+1;
    Subj{ri,ci}=[Subj{ri,ci}, {char(o.participant)}];
    if isempty(Rsp{ri,ci}), Rsp{ri,ci}=zeros(1,numel(tMs)); end
    Rsp{ri,ci}=Rsp{ri,ci}+o.sumRsp; Nrsp(ri,ci)=Nrsp(ri,ci)+o.nRsp;
end
Nsubj=zeros(nR,nC);
for ri=1:nR, for ci=1:nC, if ~isempty(Subj{ri,ci}), Nsubj(ri,ci)=numel(unique(Subj{ri,ci})); end, end, end

for s=1:3
    F=freqs.(sets{s});
    maps=cell(nR,nC); allv=[];
    for ri=1:nR, for ci=1:nC
        if ~isempty(Zsum.(sets{s}){ri,ci}) && Ss(ri,ci)>0
            maps{ri,ci}=Zsum.(sets{s}){ri,ci}/Ss(ri,ci);      % mean across sessions of within-channel z
            allv=[allv; maps{ri,ci}(:)]; %#ok<AGROW>
        end
    end, end
    cl=prctile(allv(isfinite(allv)),[2 98]); cl=max(abs(cl)); if isempty(cl)||cl==0, cl=1; end
    % global rsp range for overlay scaling (shared across the set)
    rg=[]; for ri=1:nR, for ci=1:nC, if Nrsp(ri,ci)>0, rg=[rg, Rsp{ri,ci}/Nrsp(ri,ci)]; end, end, end %#ok<AGROW>
    rmin=prctile(rg,1); rmax=prctile(rg,99); if rmax<=rmin, rmax=rmin+1; end

    fig=figure('Units','pixels','Position',[20 20 2600 2050],'Color','w','Visible','off');
    tl=tiledlayout(nR,nC,'TileSpacing','compact','Padding','compact');
    for ri=1:nR
      for ci=1:nC
        ax=nexttile;
        if ~isempty(maps{ri,ci})
            imagesc(tMs, F, maps{ri,ci}, [-cl cl]); axis xy; hold on;
            xline(0,'w-','LineWidth',1);
            if Nrsp(ri,ci)>0    % respiration overlay scaled to freq axis
                r=Rsp{ri,ci}/Nrsp(ri,ci);
                ry=F(1)+(r-rmin)/(rmax-rmin)*(F(end)-F(1));
                plot(tMs, ry, 'k-', 'LineWidth',1.6);
            end
            text(tMs(1)+150, F(end), sprintf(' %d subj, %d sess', Nsubj(ri,ci), Ss(ri,ci)), ...
                'Color','k','FontSize',18,'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','left');
        else
            imagesc(tMs, F, zeros(numel(F),numel(tMs)), [-cl cl]); axis xy;
            text(mean(tMs), mean(F), 'no data','Color','k','FontSize',13,'HorizontalAlignment','center');
        end
        if ri==1, title(cols{ci},'FontSize',20,'FontWeight','bold','Interpreter','none'); end
        if ci==1, ylabel(rows{ri},'FontSize',17,'FontWeight','bold','Interpreter','none');
        else, set(ax,'YTickLabel',[]); end
        if ri<nR, set(ax,'XTickLabel',[]); end
        set(ax,'FontSize',14);
      end
    end
    try, colormap(turbo); catch, colormap(jet); end
    cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='baseline-z (-500..-100 ms)'; cb.FontSize=15; cb.Label.FontSize=15;
    xlabel(tl,'Time from inhale/sniff onset (ms)','FontSize',17);
    title(tl, sprintf('Inhale/sniff-locked group spectrograms — %s (columns = group, rows = task; black = mean respiration)', setLabel{s}),'FontSize',19,'FontWeight','bold');
    exportgraphics(fig, fullfile(fdir, sprintf('spectrograms_%s.png', sets{s})),'Resolution',130);
    close(fig);
    fprintf('wrote spectrograms_%s.png (cells with data %d/%d, clim +-%.2f)\n', sets{s}, sum(Ns(:)>0), nR*nC, cl);
end
end

function ci=assign_col(o, classOf)
ci=[]; coh=char(o.cohort); part=char(o.participant); sn=o.sessNum;
if strcmpi(coh,'OBE'), ci=1; return; end
if ~strcmpi(coh,'Dupi'), return; end
cl=''; if isKey(classOf,part), cl=classOf(part); end
isR=strcmpi(cl,'responder'); isN=strcmpi(cl,'non-responder');
if ~isR && ~isN, return; end
if isR, if ismember(sn,[2 3]), ci=2; else, ci=3; end
else,   if ismember(sn,[2 3]), ci=4; else, ci=5; end, end
end
