function taskcmp_bandz()
% TASKCMP_BANDZ  The exact metric the user described: mean baseline-z across the
% WHOLE 25-58 Hz band at each time point, then averaged over a post-onset time
% window. Paired audiobook-vs-focus across the 6 control subjects. Also prints the
% full band-z time course per condition to test the "audiobook lower early,
% higher/equal later" crossover. Uses the spectro2 mid-band z-maps (same
% within-channel baseline-z, -500..-100 ms, as assemble_spectrograms).
here=fileparts(mfilename('fullpath')); proj=fileparts(here);
sdir=fullfile(proj,'out','spectro2'); tdir=fullfile(proj,'out','tables');
fdir=fullfile(proj,'out','figs','taskcmp'); if ~exist(fdir,'dir'), mkdir(fdir); end
CLIPZ=10;

fa=dir(fullfile(sdir,'*OBE*__audiobook.mat')); recs=struct('file',{},'subj',{},'cond',{});
for i=1:numel(fa)
    base=strrep(fa(i).name,'__audiobook.mat',''); p=strsplit(base,'_'); subj=p{4};
    recs(end+1)=struct('file',fullfile(sdir,fa(i).name),'subj',subj,'cond','audiobook'); %#ok<AGROW>
    ff=fullfile(sdir,[base '__focusedBreathing.mat']);
    if exist(ff,'file'), recs(end+1)=struct('file',ff,'subj',subj,'cond','focus'); end %#ok<AGROW>
end
tMs=[]; Fmid=[]; SB=struct('subj',{{}},'cond',{{}},'bandz',{{}});
for r=1:numel(recs)
    o=load(recs(r).file); o=o.out; if o.nBreaths==0, continue; end
    if isempty(tMs), tMs=o.tMs; Fmid=o.mid.freqs; end
    nb=o.nBreaths*o.nBaseSamp;
    mp=o.mid.sumPow/o.nBreaths; bmu=o.mid.baseSum/nb; bv=o.mid.baseSumSq/nb-bmu.^2; bv(bv<=0)=eps;
    z=(mp-bmu)./sqrt(bv/o.nBreaths); z=max(min(z,CLIPZ),-CLIPZ);      % [freq x time]
    fm=Fmid>=25 & Fmid<=58;                                          % whole 25-58 Hz band
    SB.subj{end+1}=recs(r).subj; SB.cond{end+1}=recs(r).cond; SB.bandz{end+1}=mean(z(fm,:),1);  % 1 x time
end
% subject-level band-z time course per condition (pool sessions)
su=unique(SB.subj,'stable'); nS=numel(su); nT=numel(tMs);
A=nan(nS,nT); Fo=nan(nS,nT);
for i=1:nS
    ia=find(strcmp(SB.subj,su{i}) & strcmp(SB.cond,'audiobook'));
    ic=find(strcmp(SB.subj,su{i}) & strcmp(SB.cond,'focus'));
    if ~isempty(ia), A(i,:)=mean(cat(1,SB.bandz{ia}),1); end
    if ~isempty(ic), Fo(i,:)=mean(cat(1,SB.bandz{ic}),1); end
end
keep=all(isfinite(A),2)&all(isfinite(Fo),2); A=A(keep,:); Fo=Fo(keep,:); su=su(keep); nS=size(A,1);
fprintf('paired subjects n=%d (%s)\n', nS, strjoin(su,','));

% ---- windowed paired tests (THE requested metric = window 0-500) ----
W={ '0-500',[0 500]; '500-1000',[500 1000]; '1000-1500',[1000 1500];
    '1500-2000',[1500 2000]; '2000-3000',[2000 3000]; '0-2000',[0 2000] };
rows=table();
for k=1:size(W,1)
    tw=W{k,2}; tm=tMs>=tw(1)&tMs<=tw(2);
    a=mean(A(:,tm),2); f=mean(Fo(:,tm),2); d=a-f;
    dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f);
    rows=[rows; table(W(k,1),mean(a),mean(f),mean(d),dz,pt,pw, ...
        'VariableNames',{'window_ms','mean_audiobook','mean_focus','diff_A_minus_F','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end
writetable(rows, fullfile(tdir,'taskcmp_bandz.csv'));
disp(rows);

% ---- time-course plot: mean +- SEM across subjects, both conditions + difference ----
mA=mean(A,1); sA=std(A,0,1)/sqrt(nS); mF=mean(Fo,1); sF=std(Fo,0,1)/sqrt(nS);
fig=figure('Position',[40 40 1150 560],'Color','w','Visible','off');
subplot(1,2,1); hold on;
fill([tMs fliplr(tMs)],[mA+sA fliplr(mA-sA)],[0.2 0.4 0.8],'FaceAlpha',.18,'EdgeColor','none');
fill([tMs fliplr(tMs)],[mF+sF fliplr(mF-sF)],[0.85 0.3 0.1],'FaceAlpha',.18,'EdgeColor','none');
pA=plot(tMs,mA,'-','Color',[0.2 0.4 0.8],'LineWidth',2);
pF=plot(tMs,mF,'-','Color',[0.85 0.3 0.1],'LineWidth',2);
xline(0,'k-'); yline(0,'k:'); xlim([-1000 4000]);
xlabel('Time from inhale onset (ms)'); ylabel('mean 25-58 Hz baseline-z'); legend([pA pF],{'audiobook','focus'},'Location','northeast');
title(sprintf('Whole-band (25-58 Hz) gamma z time course (control, n=%d)',nS)); set(gca,'FontSize',11);
subplot(1,2,2); hold on;
d=A-Fo; md=mean(d,1); sd=std(d,0,1)/sqrt(nS);
fill([tMs fliplr(tMs)],[md+sd fliplr(md-sd)],[0.3 0.3 0.3],'FaceAlpha',.2,'EdgeColor','none');
plot(tMs,md,'k-','LineWidth',2); xline(0,'k-'); yline(0,'k:'); xlim([-1000 4000]);
xlabel('Time from inhale onset (ms)'); ylabel('\Deltaz (audiobook - focus)');
title('Paired difference (+ = audiobook higher)'); set(gca,'FontSize',11);
exportgraphics(fig, fullfile(fdir,'bandz_timecourse.png'),'Resolution',140); close(fig);
fprintf('wrote taskcmp_bandz.csv + bandz_timecourse.png\n');
end
