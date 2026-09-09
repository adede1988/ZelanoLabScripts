function taskcmp_roi()
% TASKCMP_ROI  Lens F of the audiobook-vs-focusedBreathing sweep: a-priori
% time x frequency ROI paired contrasts on the control spectro2 z-maps. Unlike
% the whole-map cluster-permutation test (which controls FWER over ~60k pixels
% and can miss focal effects), this averages within a handful of hypothesis-driven
% ROIs and does a paired test across the 6 control subjects (higher power, few
% comparisons). gammaPowerZ = within-channel baseline-z (-500..-100 ms), same as
% assemble_spectrograms. Subjects with two sessions (AD) are pooled to one map.
here=fileparts(mfilename('fullpath')); proj=fileparts(here);
sdir=fullfile(proj,'out','spectro2'); tdir=fullfile(proj,'out','tables');
fdir=fullfile(proj,'out','figs','taskcmp'); if ~exist(fdir,'dir'), mkdir(fdir); end
CLIPZ=10; sets={'low','mid','high'};

fa=dir(fullfile(sdir,'*OBE*__audiobook.mat'));
recs=struct('file',{},'subj',{},'cond',{});
for i=1:numel(fa)
    base=strrep(fa(i).name,'__audiobook.mat',''); p=strsplit(base,'_'); subj=p{4};
    recs(end+1)=struct('file',fullfile(sdir,fa(i).name),'subj',subj,'cond','audiobook'); %#ok<AGROW>
    ff=fullfile(sdir,[base '__focusedBreathing.mat']);
    if exist(ff,'file'), recs(end+1)=struct('file',ff,'subj',subj,'cond','focus'); end %#ok<AGROW>
end
tMs=[]; F=struct(); Z=struct();
for s=1:3, Z.(sets{s})=struct('subj',{{}},'cond',{{}},'map',{{}}); end
for r=1:numel(recs)
    o=load(recs(r).file); o=o.out; if o.nBreaths==0, continue; end
    if isempty(tMs), tMs=o.tMs; for s=1:3, F.(sets{s})=o.(sets{s}).freqs; end, end
    nb=o.nBreaths*o.nBaseSamp;
    for s=1:3
        mp=o.(sets{s}).sumPow/o.nBreaths; bmu=o.(sets{s}).baseSum/nb;
        bv=o.(sets{s}).baseSumSq/nb-bmu.^2; bv(bv<=0)=eps;
        z=(mp-bmu)./sqrt(bv/o.nBreaths); z=max(min(z,CLIPZ),-CLIPZ);
        Z.(sets{s}).subj{end+1}=recs(r).subj; Z.(sets{s}).cond{end+1}=recs(r).cond; Z.(sets{s}).map{end+1}=z;
    end
end
% subject-level maps per condition (pool sessions)
subM=struct();
for s=1:3
    band=sets{s}; su=unique(Z.(band).subj,'stable'); subM.(band)=struct();
    for c={'audiobook','focus'}
        M=cell(1,numel(su));
        for i=1:numel(su)
            idx=find(strcmp(Z.(band).subj,su{i}) & strcmp(Z.(band).cond,c{1}));
            if isempty(idx), M{i}=[]; continue; end
            m=Z.(band).map{idx(1)}; for j=2:numel(idx), m=m+Z.(band).map{idx(j)}; end
            M{i}=m/numel(idx);
        end
        subM.(band).(c{1})=M; subM.(band).subj=su;
    end
end

% ---- a-priori ROIs: {name, band, tMs[lo hi], Hz[lo hi]} ----
R={ 'periOnset_gamma','mid',[0 300],[25 58];
    'early_gamma','mid',[0 500],[30 58];
    'mid_gamma','mid',[500 1500],[30 58];
    'late_gamma','mid',[1500 3000],[30 58];
    'low_early','low',[0 1000],[2 25];
    'low_late','low',[1000 3000],[2 25];
    'high_early','high',[0 500],[60 150];
    'high_mid','high',[500 1500],[60 150]};
rows=table();
for k=1:size(R,1)
    nm=R{k,1}; band=R{k,2}; tw=R{k,3}; fw=R{k,4};
    su=subM.(band).subj; A=subM.(band).audiobook; Fo=subM.(band).focus;
    tm=tMs>=tw(1)&tMs<=tw(2); fm=F.(band)>=fw(1)&F.(band)<=fw(2);
    dA=[]; dF=[]; keep={};
    for i=1:numel(su)
        if isempty(A{i})||isempty(Fo{i}), continue; end
        a=mean(A{i}(fm,tm),'all'); f=mean(Fo{i}(fm,tm),'all');
        dA(end+1)=a; dF(end+1)=f; keep{end+1}=su{i}; %#ok<AGROW>
    end
    d=dA-dF; n=numel(d);
    if n>=3
        dz=mean(d)/std(d); [~,pt]=ttest(dA,dF); pw=signrank(dA,dF);
    else, dz=NaN; pt=NaN; pw=NaN; end
    rows=[rows; table({nm},{band},tw(1),tw(2),fw(1),fw(2),n,mean(dA),mean(dF),mean(d),dz,pt,pw, ...
        'VariableNames',{'roi','band','t0','t1','f0','f1','n','mean_audiobook','mean_focus','mean_diff','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end
rows.t_p_fdr = bh_fdr(rows.t_p);
writetable(rows, fullfile(tdir,'taskcmp_F_roi.csv'));
disp(rows(:,{'roi','n','mean_diff','dz','t_p','t_p_fdr','wilcox_p'}));

% plot: ROI dz with significance
fig=figure('Position',[40 40 1100 520],'Color','w','Visible','off');
dz=rows.dz; [~,ord]=sort(dz); barh(dz(ord)); yticks(1:height(rows)); yticklabels(strrep(rows.roi(ord),'_','\_'));
hold on; xline(0,'k-'); sig=rows.t_p_fdr(ord)<0.05;
for i=find(sig)', text(dz(ord(i)),i,'  *','FontSize',16,'FontWeight','bold'); end
xlabel('Cohen dz (audiobook - focus)  [* = FDR<0.05]','FontSize',12); title('Lens F: a-priori TF-ROI paired contrasts (control, n=6)','FontSize',13,'FontWeight','bold');
set(gca,'FontSize',11); exportgraphics(fig,fullfile(fdir,'F_roi_dz.png'),'Resolution',140); close(fig);
fprintf('wrote taskcmp_F_roi.csv + F_roi_dz.png\n');
end

function q=bh_fdr(p)
p=p(:); n=numel(p); q=nan(n,1); ok=~isnan(p); pv=p(ok); [sp,ix]=sort(pv); m=numel(sp);
adj=sp.*m./(1:m)'; for i=m-1:-1:1, adj(i)=min(adj(i),adj(i+1)); end
qq=nan(m,1); qq(ix)=min(adj,1); q(ok)=qq;
end
