function scratch_theta_pac()
% SCRATCH: theta(4-8 Hz)-gamma(25-58 Hz) coupling + theta power, control breathing,
% audiobook vs focus. Per session -> subject (AD's 2 sessions averaged) -> paired (n=6).
% (1) ITPC of the 4-8 Hz neural phase across breaths at 100 respiratory-cycle-warped
%     alignments (onset->next-onset warped to 100 points; normalises longer focus breaths).
% (2) Tort MI PAC (4-8 Hz phase -> 25-58 Hz amplitude) within FIVE landmark breath periods:
%       inhale rise onset->inhale peak; inhale fall peak->downward cross back through onset
%       level; exhale rise cross->trough; exhale fall trough->85% recovery toward onset;
%       pause 85% recovery->next onset.
% (3) Theta power (|4-8 Hz|, baseline-z to -500..-100 ms) time course over first 1000 ms.
% (4) Mean theta-z within each of the five breath periods.
% Best macBP channel per session. Scratch only.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
idx=readtable(fullfile(proj,'out','tables','session_index.csv'),'TextType','string');
best=readtable(fullfile(proj,'out','tables','macbp_best.csv'),'TextType','string');
idx=idx(idx.cohort=="OBE" & idx.task=="breathingTask" & idx.onDisk==1, :);
fs=500; NPH=100; NBIN=18; NTS=round(1.0*fs)+1; periods={'inhale_rise','inhale_fall','exhale_rise','exhale_fall','pause'};
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus');

Sres=struct('subj',{},'cond',{},'itpc',{},'pac',{},'thTS',{},'thPer',{},'nbr',{});
for r=1:height(idx)
    id=char(idx.sessID(r)); fp=char(idx.finalPath(r));
    m=strcmp(best.sessID,id)&strcmp(best.task,"breathingTask"); if ~any(m), continue; end
    bl=char(string(best.bestLabel(find(m,1)))); pp=strsplit(id,'_'); subj=pp{4};
    fprintf('%s (%s) bl=%s ... ', id, subj, bl); t0=tic;
    od=[]; for a=1:3, try, S=load(fp); fn=fieldnames(S); od=S.(fn{1}); clear S; break; catch e, fprintf('[retry %d] ',a); pause(3); end, end
    if isempty(od), fprintf('LOAD FAIL\n'); continue; end
    labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1);
    if isempty(ci), fprintf('no bl\n'); continue; end
    sig=double(od.data(ci,:)); sig=fillmissing(sig,'linear'); sig=fillmissing(sig,'nearest'); N=numel(sig);
    isR=cellfun(@(x)~isempty(regexpi(x,'rsp','once')),labs); rspAll=od.data(isR,:);
    ri=1; if isfield(od,'rspIDX'),ri=double(od.rspIDX);end; rf=1; if isfield(od,'rspFlip'),rf=double(od.rspFlip);end
    if ri>size(rspAll,1),ri=1; end
    rsp=double(rspAll(ri,:))*rf; rsp=fillmissing(rsp,'linear'); rsp=fillmissing(rsp,'nearest');
    thA=hilbert(bp(sig,4,8,fs)); thPh=angle(thA); thAmp=abs(thA);      % theta analytic/phase/amp
    gAmp=abs(hilbert(bp(sig,25,58,fs)));                              % gamma amplitude
    bd=od.behDat; tc=strtrim(string(getcol(bd,'task')));
    onS=round(coerce(bd.finalOnset)); inMax=round(coerce(getnum(bd,'inMaxTim')));
    exMin=round(coerce(getnum(bd,'exMinTim'))); Yon=coerce(getnum(bd,'Yonset'));
    good=~isnan(onS)&onS>0; ons=sort(onS(good));
    for cnd=1:2
        sel = (cnd==1).*strcmpi(tc,'audio') + (cnd==2).*isFocus(tc); sel=logical(sel(:).') & good(:).';
        rows=find(sel);
        itpcSum=zeros(1,NPH); itpcN=0; tsSum=zeros(1,NTS); tsN=0; thPer=[];
        phBins=cell(1,5); ampBins=cell(1,5); for pI=1:5, phBins{pI}=[]; ampBins{pI}=[]; end
        for kk=rows
            o=onS(kk); ip=inMax(kk); tr=exMin(kk); y0=Yon(kk);
            no=next_onset(ons,o,N); if isnan(no)||no<=o, continue; end
            if isnan(ip)||ip<=o||ip>=no, continue; end
            rc=downcross(rsp,ip,y0,no); if isnan(rc)||isnan(tr)||tr<=rc||tr>=no, continue; end
            lvl=tr_level(rsp,tr,y0); r85=upcross(rsp,tr,lvl,no); if isnan(r85)||r85>=no, r85=no; end
            bnds={[o ip],[ip rc],[rc tr],[tr r85],[r85 no]};
            % ITPC (breath-warped)
            seg=thA(o:no); warped=interp1(1:numel(seg), seg, linspace(1,numel(seg),NPH), 'linear');
            itpcSum=itpcSum + exp(1i*angle(warped)); itpcN=itpcN+1;
            % PAC bins
            for pI=1:5, s=max(round(bnds{pI}(1)),1); e=min(round(bnds{pI}(2)),N);
                if e>s+2, phBins{pI}=[phBins{pI}, thPh(s:e)]; ampBins{pI}=[ampBins{pI}, gAmp(s:e)]; end, end
            % theta power baseline-z
            b0=o-round(0.5*fs); b1=o-round(0.1*fs);
            if b0>=1
                bmu=mean(thAmp(b0:b1)); bsd=std(thAmp(b0:b1)); if bsd<=0, bsd=eps; end
                e1=o+NTS-1; if e1<=N, tsSum=tsSum+(thAmp(o:e1)-bmu)/bsd; tsN=tsN+1; end
                tp=nan(1,5); for pI=1:5, s=max(round(bnds{pI}(1)),1); e=min(round(bnds{pI}(2)),N);
                    if e>s, tp(pI)=(mean(thAmp(s:e))-bmu)/bsd; end, end
                thPer=[thPer; tp]; %#ok<AGROW>
            end
        end
        if itpcN<5, fprintf('[cond%d thin] ',cnd); end
        pac=nan(1,5); for pI=1:5, pac(pI)=tortMI(phBins{pI}, ampBins{pI}, NBIN); end
        thPerM=mean(thPer,1,'omitnan'); if isempty(thPerM), thPerM=nan(1,5); end
        Sres(end+1)=struct('subj',subj,'cond',ternary(cnd==1,'audiobook','focus'), ...
            'itpc',abs(itpcSum)/max(itpcN,1),'pac',pac,'thTS',tsSum/max(tsN,1),'thPer',thPerM,'nbr',itpcN); %#ok<AGROW>
    end
    fprintf('(%.0fs)\n', toc(t0));
end
save(fullfile(proj,'out','tables','scratch_theta_pac.mat'),'Sres','periods','NPH','NTS','fs');
aggregate_and_plot(proj, Sres, periods, NPH, NTS, fs);
fprintf('DONE scratch_theta_pac (%d records)\n', numel(Sres));
end

% ================= aggregation + plots =================
function aggregate_and_plot(proj, Sres, periods, NPH, NTS, fs)
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
subs=unique({Sres.subj},'stable'); conds={'audiobook','focus'};
ITP=struct(); TS=struct(); for ci=1:2, ITP.(conds{ci})=nan(numel(subs),NPH); TS.(conds{ci})=nan(numel(subs),NTS); end
PACt=table(); THt=table();
for si=1:numel(subs), for ci=1:2
    m=strcmp({Sres.subj},subs{si}) & strcmp({Sres.cond},conds{ci}); if ~any(m), continue; end
    R=Sres(m);
    ITP.(conds{ci})(si,:)=mean(cell2mat({R.itpc}'),1,'omitnan');
    TS.(conds{ci})(si,:)=mean(cell2mat({R.thTS}'),1,'omitnan');
    pac=mean(cell2mat({R.pac}'),1,'omitnan'); thp=mean(cell2mat({R.thPer}'),1,'omitnan');
    PACt=[PACt; table(string(subs{si}),string(conds{ci}),pac(1),pac(2),pac(3),pac(4),pac(5),'VariableNames',[{'subject','condition'},periods])]; %#ok<AGROW>
    THt=[THt; table(string(subs{si}),string(conds{ci}),thp(1),thp(2),thp(3),thp(4),thp(5),'VariableNames',[{'subject','condition'},periods])]; %#ok<AGROW>
end, end
writetable(PACt, fullfile(tdir,'scratch_theta_pac_subject.csv'));
writetable(THt, fullfile(tdir,'scratch_thetaPower_period_subject.csv'));
St_pac=paired_stats(PACt, periods, subs); writetable(St_pac, fullfile(tdir,'scratch_theta_pac_stats.csv'));
St_th =paired_stats(THt , periods, subs); writetable(St_th , fullfile(tdir,'scratch_thetaPower_period_stats.csv'));
cat_disp('PAC (Tort MI) paired stats', St_pac); cat_disp('Theta-power per-period paired stats', St_th);

% ITPC plot (2 panels)
xph=linspace(0,100,NPH); ymax=max(0.02, max([ITP.audiobook(:); ITP.focus(:)],[],'omitnan')*1.1);
fig=figure('Position',[30 30 1250 500],'Color','w','Visible','off');
for ci=1:2, subplot(1,2,ci); hold on; M=ITP.(conds{ci}); cm=lines(size(M,1));
  for si=1:size(M,1), if all(isfinite(M(si,:))), plot(xph,M(si,:),'-','Color',[cm(si,:) 0.3]); end, end
  plot(xph,mean(M,1,'omitnan'),'k-','LineWidth',2.5); xlim([0 100]); ylim([0 ymax]);
  xlabel('% respiratory cycle (onset\rightarrownext onset)'); ylabel('4-8 Hz ITPC across breaths');
  title(sprintf('ITPC — %s',conds{ci}),'FontWeight','bold'); set(gca,'FontSize',10); end
sgtitle('Control: 4-8 Hz inter-breath phase coherence, breath-warped (faint = subject)','FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'theta_itpc.png'),'Resolution',140); close(fig);

% PAC + theta-period paired bar plots
paired_bar(PACt, periods, subs, 'theta(4-8)\rightarrowgamma(25-58) PAC (Tort MI) by breath period','Tort MI', fullfile(fdir,'theta_pac_periods.png'));
paired_bar(THt , periods, subs, 'Mean theta-z by breath period','theta power (baseline-z)', fullfile(fdir,'thetaPower_periods.png'));

% Theta power time course (0-1000 ms) + difference
tms=(0:NTS-1)/fs*1000; A=TS.audiobook; F=TS.focus; nS=size(A,1);
mA=mean(A,1,'omitnan'); sA=std(A,0,1,'omitnan')/sqrt(nS); mF=mean(F,1,'omitnan'); sF=std(F,0,1,'omitnan')/sqrt(nS);
fig=figure('Position',[30 30 1200 500],'Color','w','Visible','off');
subplot(1,2,1); hold on;
fill([tms fliplr(tms)],[mA+sA fliplr(mA-sA)],[0.2 0.4 0.8],'FaceAlpha',.18,'EdgeColor','none');
fill([tms fliplr(tms)],[mF+sF fliplr(mF-sF)],[0.85 0.3 0.1],'FaceAlpha',.18,'EdgeColor','none');
pa=plot(tms,mA,'-','Color',[0.2 0.4 0.8],'LineWidth',2); pf=plot(tms,mF,'-','Color',[0.85 0.3 0.1],'LineWidth',2);
yline(0,'k:'); xlabel('ms after inhale onset'); ylabel('theta power (baseline-z)'); legend([pa pf],{'audiobook','focus'});
title(sprintf('Theta power time course (n=%d)',nS),'FontWeight','bold'); set(gca,'FontSize',10);
subplot(1,2,2); hold on; d=A-F; md=mean(d,1,'omitnan'); sd=std(d,0,1,'omitnan')/sqrt(nS);
fill([tms fliplr(tms)],[md+sd fliplr(md-sd)],[.3 .3 .3],'FaceAlpha',.2,'EdgeColor','none');
plot(tms,md,'k-','LineWidth',2); yline(0,'k:'); xlabel('ms after inhale onset'); ylabel('\Delta theta-z (audiobook - focus)');
title('Paired difference (+ = audiobook higher)','FontWeight','bold'); set(gca,'FontSize',10);
exportgraphics(fig, fullfile(fdir,'thetaPower_timecourse.png'),'Resolution',140); close(fig);
fprintf('wrote theta_itpc.png, theta_pac_periods.png, thetaPower_periods.png, thetaPower_timecourse.png + CSVs\n');
end

function St=paired_stats(Tb, periods, subs)
St=table();
for pI=1:5
  A=nan(numel(subs),1); F=nan(numel(subs),1);
  for si=1:numel(subs)
    a=Tb.(periods{pI})(Tb.subject==subs{si}&Tb.condition=="audiobook");
    f=Tb.(periods{pI})(Tb.subject==subs{si}&Tb.condition=="focus");
    if ~isempty(a),A(si)=a(1);end; if ~isempty(f),F(si)=f(1);end
  end
  ok=isfinite(A)&isfinite(F); a=A(ok); f=F(ok); d=a-f; n=numel(d);
  if n>=3, dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f); else dz=NaN;pt=NaN;pw=NaN; end
  St=[St; table(periods(pI),n,mean(a),mean(f),mean(d),dz,pt,pw,'VariableNames',{'period','n','mean_audiobook','mean_focus','diff','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end
end
function paired_bar(Tb, periods, subs, ttl, ylab, outpng)
St=paired_stats(Tb, periods, subs); fig=figure('Position',[30 30 1150 540],'Color','w','Visible','off'); hold on;
x=1:5; w=0.35; bar(x-w/2, St.mean_audiobook, w,'FaceColor',[0.2 0.4 0.8]); bar(x+w/2, St.mean_focus, w,'FaceColor',[0.85 0.3 0.1]);
for pI=1:5, for si=1:numel(subs)
    a=Tb.(periods{pI})(Tb.subject==subs{si}&Tb.condition=="audiobook"); f=Tb.(periods{pI})(Tb.subject==subs{si}&Tb.condition=="focus");
    if ~isempty(a)&&~isempty(f), plot([pI-w/2 pI+w/2],[a(1) f(1)],'-','Color',[.5 .5 .5 .5]); end, end, end
for pI=1:5, if isfinite(St.t_p(pI))&&St.t_p(pI)<0.1, text(pI, max(St.mean_audiobook(pI),St.mean_focus(pI)), sprintf('p=%.2f',St.t_p(pI)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9); end, end
set(gca,'XTick',x,'XTickLabel',strrep(periods,'_','\_')); ylabel(ylab); legend({'audiobook','focus'},'Location','best');
title(['Control: ' ttl],'FontWeight','bold'); set(gca,'FontSize',11); exportgraphics(fig,outpng,'Resolution',140); close(fig);
end
function cat_disp(ttl, St), fprintf('\n=== %s ===\n',ttl); disp(St); end

% ================= signal helpers =================
function y=bp(x,lo,hi,fs), [b,a]=butter(4,[lo hi]/(fs/2),'bandpass'); y=filtfilt(b,a,double(x(:)')); end
function no=next_onset(ons,o,N), j=find(ons>o,1,'first'); if isempty(j), no=N; else, no=ons(j); end, end
function rc=downcross(rsp,ip,y0,no), rc=NaN; if isnan(ip)||isnan(y0), return; end
  for t=ip+1:min(no,numel(rsp)), if rsp(t)<=y0, rc=t; return; end, end, end
function lvl=tr_level(rsp,tr,y0), lvl=rsp(tr)+0.85*(y0-rsp(tr)); end
function uc=upcross(rsp,tr,lvl,no), uc=NaN; if isnan(tr), return; end
  for t=tr+1:min(no,numel(rsp)), if rsp(t)>=lvl, uc=t; return; end, end, end
function mi=tortMI(ph, amp, nbin), mi=NaN; ph=ph(:); amp=amp(:); ok=isfinite(ph)&isfinite(amp); ph=ph(ok); amp=amp(ok);
  if numel(ph)<nbin*10, return; end
  edges=linspace(-pi,pi,nbin+1); [~,bi]=histc(ph,edges); bi(bi==nbin+1)=nbin; %#ok<HISTC>
  ma=zeros(nbin,1); for b=1:nbin, v=amp(bi==b); if isempty(v), ma(b)=eps; else, ma(b)=mean(v); end, end
  p=ma/sum(ma); p(p<=0)=eps; H=-sum(p.*log(p)); mi=(log(nbin)-H)/log(nbin); end
function v=coerce(c), if isnumeric(c), v=double(c); return; end
  if iscell(c), v=nan(numel(c),1); for i=1:numel(c), try, e=c{i}; if isnumeric(e)&&~isempty(e), v(i)=double(e(1)); else, v(i)=double(string(e)); end, catch, v(i)=NaN; end, end, return; end
  try, v=double(string(c)); catch, v=nan(numel(c),1); end, end
function v=getnum(bd,nm), if any(strcmp(bd.Properties.VariableNames,nm)), v=coerce(bd.(nm)); else, v=nan(height(bd),1); end, end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
function y=ternary(c,a,b), if c, y=a; else, y=b; end, end
