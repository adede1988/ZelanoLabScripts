function cluster_perm_audio_focus()
% CLUSTER_PERM_AUDIO_FOCUS  Control-only comparison of audiobook vs
% focusedBreathing group spectrograms (gammaPowerZ ~ task) at every
% time x frequency pixel, with Maris-Oostenveld cluster-based permutation
% correction, for the low (2-25), mid (25-58) and high (60-150 Hz) bands.
%
% Design: controls = OBE only. gammaPowerZ is the per-session within-channel
% baseline-z (-500..-100 ms) map, identical to assemble_spectrograms. Sessions
% are aggregated to SUBJECT level within condition (subject AD's two sessions
% count once). Every control breathing session has BOTH conditions (after the
% focus-label fix: focus | naturalFocus | slowFocus), so this is a WITHIN-SUBJECT
% PAIRED contrast: for each subject with both conditions, D = z(audiobook) -
% z(focusedBreathing). Pixel statistic = one-sample t of D across subjects
% (df = nSubj-1). Clusters = 4-connected pixels with |t|>t_crit(p<.05,two-tailed);
% cluster mass = sum(t). Null = max |cluster mass| over all 2^nSubj exhaustive
% sign-flips of the paired differences (the exact paired permutation scheme).
% A cluster is significant at alpha=.05 if its |mass| exceeds the 95th percentile
% of that null (p = fraction of sign-flips with a larger max |mass|). Outputs 3
% difference-score spectrograms (mean D = audiobook - focus) with significant
% clusters outlined solid and 0.05<=p<0.10 "trend" clusters dashed, plus a stats CSV.
%
% If some subjects lack one condition, they are dropped from the paired test
% (reported). With nSubj<=~14, exhaustive sign-flip is used; else 2000 random.

here=fileparts(mfilename('fullpath')); proj=fileparts(here);
sdir=fullfile(proj,'out','spectro2');
fdir=fullfile(proj,'out','figs','scratch'); if ~exist(fdir,'dir'), mkdir(fdir); end
tdir=fullfile(proj,'out','tables'); if ~exist(tdir,'dir'), mkdir(tdir); end
sets={'low','mid','high'}; setLabel={'2-25 Hz','25-58 Hz','60-150 Hz'};
CLIPZ=10; ALPHA_PIX=0.05; ALPHA_CLU=0.05;

% ---- collect control (OBE) audiobook + focusedBreathing sessions ----
fa=dir(fullfile(sdir,'*OBE*__audiobook.mat'));
recs=struct('file',{},'subj',{},'cond',{});
for i=1:numel(fa)
    base=strrep(fa(i).name,'__audiobook.mat','');
    subj=subj_from_id(base);
    recs(end+1)=struct('file',fullfile(sdir,fa(i).name),'subj',subj,'cond','audiobook'); %#ok<AGROW>
    ffoc=fullfile(sdir,[base '__focusedBreathing.mat']);
    if exist(ffoc,'file'), recs(end+1)=struct('file',ffoc,'subj',subj,'cond','focusedBreathing'); end %#ok<AGROW>
end

% ---- load once, compute per-session z per band, drop nBreaths==0 ----
tMs=[]; freqs=struct();
Z=struct(); for s=1:3, Z.(sets{s})=struct('subj',{{}},'cond',{{}},'map',{{}}); end
for r=1:numel(recs)
    o=load(recs(r).file); o=o.out;
    if o.nBreaths==0, continue; end
    if isempty(tMs), tMs=o.tMs; for s=1:3, freqs.(sets{s})=o.(sets{s}).freqs; end, end
    nb=o.nBreaths*o.nBaseSamp;
    for s=1:3
        meanPow=o.(sets{s}).sumPow/o.nBreaths;
        baseMu =o.(sets{s}).baseSum/nb;
        baseVar=o.(sets{s}).baseSumSq/nb - baseMu.^2; baseVar(baseVar<=0)=eps;
        zsess=(meanPow-baseMu)./sqrt(baseVar/o.nBreaths);
        zsess=max(min(zsess,CLIPZ),-CLIPZ);
        Z.(sets{s}).subj{end+1}=recs(r).subj;
        Z.(sets{s}).cond{end+1}=recs(r).cond;
        Z.(sets{s}).map{end+1}=zsess;
    end
end

allStats=table();
for s=1:3
    band=sets{s}; F=freqs.(band);
    subj=Z.(band).subj; cond=Z.(band).cond; maps=Z.(band).map;
    % ---- subject-level maps per condition, then keep subjects with BOTH ----
    [A,subjA]=subj_mean(maps, subj, cond, 'audiobook');
    [Fo,subjF]=subj_mean(maps, subj, cond, 'focusedBreathing');
    paired=intersect(subjA, subjF, 'stable');
    D=cell(1,numel(paired));
    for i=1:numel(paired)
        D{i}=A{strcmp(subjA,paired{i})} - Fo{strcmp(subjF,paired{i})};   % audiobook - focus
    end
    nS=numel(paired);
    fprintf('[%s] paired within-subject n=%d (%s)\n', band, nS, strjoin(paired,','));
    if nS<3, fprintf('   too few paired subjects; skipping.\n'); continue; end
    Dstack=cat(3,D{:});                    % [freq x time x subj]
    df=nS-1; tcrit=tinv(1-ALPHA_PIX/2, df);

    % ---- observed one-sample t-map + clusters ----
    signs0=ones(1,nS);
    tObs=tmap_paired(Dstack, signs0);
    [cmObsPos, cLabPos]=clusters(tObs, tcrit, +1);
    [cmObsNeg, cLabNeg]=clusters(tObs,-tcrit, -1);

    % ---- exhaustive (or sampled) sign-flip permutation null ----
    if nS<=14
        P=2^nS; signMat=sign_combos(nS); nP=P;         % exact
    else
        nP=2000; signMat=(randi([0 1],nP,nS)*2-1); signMat(1,:)=1; %#ok<*RAND>
    end
    maxNull=zeros(nP,1);
    for p=1:nP
        tp=tmap_paired(Dstack, signMat(p,:));
        cp=clusters(tp, tcrit, +1); cn=clusters(tp,-tcrit,-1);
        m=0; if ~isempty(cp), m=max(m,max(abs(cp))); end
        if ~isempty(cn), m=max(m,max(abs(cn))); end
        maxNull(p)=m;
    end
    thr=prctile(maxNull, 100*(1-ALPHA_CLU));

    % ---- p per observed cluster; masks ----
    diffMap=mean(Dstack,3);                 % mean audiobook - focus
    sigMask=false(size(tObs)); trendMask=false(size(tObs));
    obs=[cmObsPos(:); cmObsNeg(:)];
    obsLab={}; for k=1:numel(cmObsPos), obsLab{end+1}=sprintf('pos%d',k); end %#ok<AGROW>
    for k=1:numel(cmObsNeg), obsLab{end+1}=sprintf('neg%d',k); end %#ok<AGROW>
    for k=1:numel(obs)
        pval=mean(maxNull >= abs(obs(k)));
        sig = pval < ALPHA_CLU;
        if k<=numel(cmObsPos), mask=(cLabPos==k); else mask=(cLabNeg==(k-numel(cmObsPos))); end
        if sig, sigMask=sigMask | mask; elseif pval<0.10, trendMask=trendMask | mask; end
        allStats=[allStats; table({band},obsLab(k),obs(k),pval,sig,sum(mask(:)), ...
            'VariableNames',{'band','cluster','mass','p_perm','sig','nPix'})]; %#ok<AGROW>
    end

    plot_diff(diffMap, tMs, F, sigMask, trendMask, setLabel{s}, nS, ...
        fullfile(fdir, sprintf('clusterperm_audio_vs_focus_%s.png', band)), thr, tcrit);
    nsig = nnz(allStats.sig(strcmp(allStats.band,band)));
    fprintf('   %d clusters; %d significant (p<%.2f); maxNull95=%.1f\n', numel(obs), nsig, ALPHA_CLU, thr);
end
writetable(allStats, fullfile(tdir,'clusterperm_audio_vs_focus_stats.csv'));
fprintf('wrote figs/scratch/clusterperm_audio_vs_focus_{low,mid,high}.png + stats CSV\n');
end

% ================= helpers =================
function s=subj_from_id(id)
p=strsplit(id,'_'); if numel(p)>=4, s=p{4}; else s=id; end
end

function [out,subjOut]=subj_mean(maps, subj, cond, want)
sel=strcmp(cond,want); us=unique(subj(sel),'stable'); out={}; subjOut={};
for i=1:numel(us)
    idx=find(sel & strcmp(subj,us{i}));
    m=maps{idx(1)}; for j=2:numel(idx), m=m+maps{idx(j)}; end
    out{end+1}=m/numel(idx); subjOut{end+1}=us{i}; %#ok<AGROW>
end
end

function t=tmap_paired(Dstack, signs)
% one-sample t across subjects of the sign-flipped paired differences
n=size(Dstack,3);
sg=reshape(signs(:),1,1,n);
X=Dstack.*sg;
m=mean(X,3); sd=std(X,0,3);
den=sd/sqrt(n); den(den<=0)=eps;
t=m./den;
end

function M=sign_combos(n)
% all 2^n +-1 sign vectors (row 1 = all +1)
P=2^n; M=ones(P,n);
for p=0:P-1
    b=bitget(p, 1:n); M(p+1,:)=1-2*b;   % 0->+1, 1->-1
end
M(1,:)=1;
end

function [mass, lab]=clusters(t, thr, sgn)
if sgn>0, bw=t>=thr; else bw=t<=thr; end
lab=bwlabel_local(bw);
K=max(lab(:)); mass=zeros(1,K);
for k=1:K, mass(k)=sum(t(lab==k)); end
end

function L=bwlabel_local(bw)
if ~isempty(which('bwlabel'))
    L=bwlabel(bw,4); return;
end
L=zeros(size(bw)); [R,C]=size(bw); cur=0;
for i=1:R, for j=1:C
    if bw(i,j) && L(i,j)==0
        cur=cur+1; stack=[i j];
        while ~isempty(stack)
            p=stack(end,:); stack(end,:)=[]; a=p(1); b=p(2);
            if a<1||a>R||b<1||b>C, continue; end
            if ~bw(a,b) || L(a,b)~=0, continue; end
            L(a,b)=cur;
            stack=[stack; a-1 b; a+1 b; a b-1; a b+1]; %#ok<AGROW>
        end
    end
end, end
end

function plot_diff(diffMap, tMs, F, sigMask, trendMask, bandLbl, nS, outpng, thr, tcrit)
fig=figure('Units','pixels','Position',[40 40 1500 620],'Color','w','Visible','off');
cl=prctile(abs(diffMap(:)),99); if ~isfinite(cl)||cl==0, cl=1; end
imagesc(tMs, F, diffMap, [-cl cl]); axis xy; hold on;
try, colormap(bluewhitered_local(256)); catch, colormap(parula); end
xline(0,'k-','LineWidth',1);
if any(trendMask(:))
    contour(tMs, F, double(trendMask), [0.5 0.5], '--', 'Color',[.35 .35 .35], 'LineWidth', 1.6);
end
if any(sigMask(:))
    contour(tMs, F, double(sigMask), [0.5 0.5], 'k-', 'LineWidth', 2.4);
end
xlabel('Time from inhale onset (ms)','FontSize',14);
ylabel('Frequency (Hz)','FontSize',14);
title(sprintf('Control: audiobook - focusedBreathing gammaPowerZ  (%s)   [paired, n=%d subj]', bandLbl, nS), ...
      'FontSize',15,'FontWeight','bold');
subtitle_local(sprintf('warm = higher during audiobook; solid = sig cluster (p<.05); dashed = trend (.05<=p<.10). paired sign-flip; |t|>%.2f forms clusters; max-null95 mass=%.0f', tcrit, thr));
cb=colorbar; cb.Label.String='\Deltaz (audiobook - focus)'; cb.FontSize=12;
set(gca,'FontSize',12);
exportgraphics(fig, outpng, 'Resolution',140); close(fig);
end

function subtitle_local(str)
try, subtitle(str,'FontSize',11,'Color',[.25 .25 .25]); catch
    ax=gca; text(ax,0.5,1.02,str,'Units','normalized','HorizontalAlignment','center','FontSize',10,'Color',[.25 .25 .25]);
end
end

function cmap=bluewhitered_local(m)
half=floor(m/2); up=linspace(0,1,half)';
top=[ones(half,1), flipud(up), flipud(up)];
bot=[up, up, ones(half,1)];
mid=[1 1 1]; cmap=[bot; mid; top]; cmap=cmap(1:m,:);
end
