function taskcmp_groups_gamma()
% Gamma whole-band (25-58 Hz) baseline-z time course, audiobook vs focus, PAIRED
% within subject, for Control (OBE, CP excluded), Dupi S1, and Dupi S2 (S3 also
% tabulated). Full epoch (-1000..4000 ms) so we can test the LATE window (3000-4000)
% where focus gamma may exceed audiobook. From the spectro2 maps (no re-extraction).
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
sdir=fullfile(proj,'out','spectro2'); tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
conds={'audiobook','focus'}; CLIPZ=10; EXCL_CTRL={'CP'};

fa=dir(fullfile(sdir,'*__audiobook.mat')); tG=[];
D=struct('grp',{},'subj',{},'cond',{},'tc',{});
for i=1:numel(fa)
    base=strrep(fa(i).name,'__audiobook.mat','');
    for c=1:2
        f=fullfile(sdir,[base '__' tern(c==1,'audiobook','focusedBreathing') '.mat']);
        if ~exist(f,'file'), continue; end; o=load(f); o=o.out; if o.nBreaths==0, continue; end
        if isempty(tG), tG=o.tMs; end
        coh=char(string(o.cohort)); sn=o.sessNum; part=char(string(o.participant));
        if strcmpi(coh,'OBE'), grp='Control'; if any(strcmpi(part,EXCL_CTRL)), continue; end
        elseif strcmpi(coh,'Dupi'), grp=sprintf('DupiS%d',sn); else, continue; end
        nb=o.nBreaths*o.nBaseSamp; mp=o.mid.sumPow/o.nBreaths; bmu=o.mid.baseSum/nb;
        bv=o.mid.baseSumSq/nb-bmu.^2; bv(bv<=0)=eps; z=(mp-bmu)./sqrt(bv/o.nBreaths); z=max(min(z,CLIPZ),-CLIPZ);
        fm=o.mid.freqs>=25&o.mid.freqs<=58;
        D(end+1)=struct('grp',grp,'subj',part,'cond',conds{c},'tc',mean(z(fm,:),1)); %#ok<AGROW>
    end
end
groups={'Control','DupiS1','DupiS2','DupiS3'};
wins={'0-500',[0 500];'500-1000',[500 1000];'1000-2000',[1000 2000];'2000-3000',[2000 3000];'3000-4000',[3000 4000]};

% subject-level time courses per group/cond
TC=struct();
for gi=1:numel(groups), g=groups{gi};
  for c=1:2, subs=unique({D(strcmp({D.grp},g)&strcmp({D.cond},conds{c})).subj},'stable');
    M=nan(numel(subs),numel(tG));
    for si=1:numel(subs), m=strcmp({D.grp},g)&strcmp({D.cond},conds{c})&strcmp({D.subj},subs{si});
      M(si,:)=mean(cell2mat({D(m).tc}'),1,'omitnan'); end
    TC.(g).(conds{c})=M; TC.(g).([conds{c} '_subs'])=subs;
  end
end

% windowed paired stats (focus - audiobook)
St=table();
for gi=1:numel(groups), g=groups{gi};
  as=TC.(g).audiobook_subs; fs2=TC.(g).focus_subs; subs=intersect(as,fs2,'stable');
  for wi=1:size(wins,1), tw=wins{wi,2}; tm=tG>=tw(1)&tG<=tw(2);
    a=nan(numel(subs),1); f=nan(numel(subs),1);
    for si=1:numel(subs)
      a(si)=mean(TC.(g).audiobook(strcmp(as,subs{si}),tm),'all');
      f(si)=mean(TC.(g).focus(strcmp(fs2,subs{si}),tm),'all');
    end
    d=f-a; n=numel(d); dz=mean(d)/std(d); [~,pt]=ttest(f,a); pw=signrank(f,a);
    St=[St; table({g},wins(wi,1),n,mean(a),mean(f),mean(d),dz,pt,pw, ...
      'VariableNames',{'group','window_ms','n','mean_audiobook','mean_focus','focus_minus_audiobook','dz','t_p','wilcox_p'})]; %#ok<AGROW>
  end
end
writetable(St, fullfile(tdir,'taskcmp_groups_gamma.csv'));
cat_show(St);

% ---- plot: rows = Control/DupiS1/DupiS2, cols = time course + difference ----
plotg={'Control','DupiS1','DupiS2'};
fig=figure('Position',[20 20 1350 1050],'Color','w','Visible','off');
allm=[]; for gi=1:3, g=plotg{gi}; allm=[allm; mean(TC.(g).audiobook,1,'omitnan')'; mean(TC.(g).focus,1,'omitnan')']; end %#ok<AGROW>
yl=[min(allm) max(allm)]+[-.1 .1];
for gi=1:3, g=plotg{gi}; nA=size(TC.(g).audiobook,1); nF=size(TC.(g).focus,1);
  mA=mean(TC.(g).audiobook,1,'omitnan'); sA=std(TC.(g).audiobook,0,1,'omitnan')/sqrt(nA);
  mF=mean(TC.(g).focus,1,'omitnan'); sF=std(TC.(g).focus,0,1,'omitnan')/sqrt(nF);
  subplot(3,2,(gi-1)*2+1); hold on;
  fill([tG fliplr(tG)],[mA+sA fliplr(mA-sA)],[.2 .4 .8],'FaceAlpha',.15,'EdgeColor','none');
  fill([tG fliplr(tG)],[mF+sF fliplr(mF-sF)],[.85 .3 .1],'FaceAlpha',.15,'EdgeColor','none');
  pa=plot(tG,mA,'-','Color',[.2 .4 .8],'LineWidth',2); pf=plot(tG,mF,'-','Color',[.85 .3 .1],'LineWidth',2);
  xline(0,'k-'); yline(0,'k:'); xlim([-1000 4000]); ylim(yl);
  patch([3000 4000 4000 3000],[yl(1) yl(1) yl(2) yl(2)],[.9 .9 .5],'FaceAlpha',.15,'EdgeColor','none');
  ylabel('25-58 Hz z'); title(sprintf('%s  (n=%d)  gamma band-z',g,min(nA,nF)),'FontWeight','bold');
  if gi==1, legend([pa pf],{'audiobook','focus'},'Location','northeast'); end
  if gi==3, xlabel('ms from inhale onset'); end; set(gca,'FontSize',10);
  subplot(3,2,(gi-1)*2+2); hold on; d=mF-mA;
  plot(tG,d,'k-','LineWidth',1.8); xline(0,'k-'); yline(0,'k:'); xlim([-1000 4000]);
  yl2=[min(d)-.1 max(d)+.1]; ylim(yl2);
  patch([3000 4000 4000 3000],[yl2(1) yl2(1) yl2(2) yl2(2)],[.9 .9 .5],'FaceAlpha',.2,'EdgeColor','none');
  sr=St(strcmp(St.group,g)&strcmp(St.window_ms,'3000-4000'),:);
  ylabel('\Delta z (focus - audiobook)'); title(sprintf('%s: diff (3-4s dz=%.2f, p=%.2f)',g,sr.dz,sr.t_p),'FontWeight','bold');
  if gi==3, xlabel('ms from inhale onset'); end; set(gca,'FontSize',10);
end
sgtitle('Gamma band-z: audiobook vs focus by group (yellow = 3000-4000 ms late window)','FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'groups_gamma_timecourse.png'),'Resolution',140); close(fig);
fprintf('wrote taskcmp_groups_gamma.csv, groups_gamma_timecourse.png\n');
end

function cat_show(St)
fprintf('\n=== Gamma band-z, focus - audiobook, paired by group/window ===\n');
for g=unique(St.group)', sub=St(strcmp(St.group,g{1}),:);
  fprintf('\n%s (n=%d):\n',g{1},sub.n(1)); disp(sub(:,{'window_ms','mean_audiobook','mean_focus','focus_minus_audiobook','dz','t_p','wilcox_p'}));
end
end
function y=tern(c,a,b), if c, y=a; else, y=b; end, end
