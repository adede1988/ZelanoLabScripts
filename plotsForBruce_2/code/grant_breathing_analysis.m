function grant_breathing_analysis(maxSess, nw)
if nargin<1||isempty(maxSess), maxSess=inf; end
if nargin<2||isempty(nw), nw=3; end   % worker count; 0 = serial
% GRANT_BREATHING_ANALYSIS  For the two breathing conditions (audiobook = "baseline",
% focusedBreathing = "ATB"), aligned to inhale onset + 200 ms:
%   (1) build group-average OB spectrograms for ALL / Dupi-S1 / Dupi-S2/3 / Control
%       (4 groups x 2 conditions = 8), display -100..+1000 ms re: the +200 ms point,
%       shared color scale (floor 4), 60 Hz line-noise band zeroed;
%   (2) extract per-breath mean within-frequency z-power across the gamma band (25-57 Hz)
%       in the 200-400 ms window (re: the +200 ms point) and test group x condition
%       differences at the session level.
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); repo=fileparts(proj);  % machine-agnostic
addpath(repo); addpath(codeDir);
T=fullfile(proj,'out','tables');
idx=readtable(fullfile(T,'session_index.csv'),'TextType','string');
mb =readtable(fullfile(T,'macbp_best.csv'),'TextType','string');
mbB=mb(mb.task=="breathingTask",:);

fs=500; F=(15:1:150).'; c1=3; ord=[3 30]; mult=1;
alignMs=200; alignSamp=round(alignMs/1000*fs);
compPre=1000; compPost=3000; padMs=600;
cpre=round(compPre/1000*fs); cpost=round(compPost/1000*fs); pad=round(padMs/1000*fs);
tMsF=(-cpre:cpost)/fs*1000;                 % relative to the +200 ms alignment point
tasks={'audiobook','focusedBreathing'};
gband=[25 57]; scwin=[200 400]; fLim=[25 60];

% ---------- build session list (unique breathing sessions) ----------
sids=strings(0);
for t=1:2
  f=dir(fullfile(proj,'out','gamma','perbreath',sprintf('*%s.csv',tasks{t})));
  for i=1:numel(f)
    nm=erase(f(i).name,sprintf('__%s.csv',tasks{t})); sids(end+1)=string(nm); %#ok<AGROW>
  end
end
sids=unique(sids);
S=struct('sessID',{},'finalPath',{},'bestLabel',{},'cohort',{},'sessNum',{},'participant',{},'abCSV',{},'fbCSV',{});
for i=1:numel(sids)
  s=sids(i);
  fp=idx.finalPath(idx.sessID==s & idx.task=="breathingTask"); if isempty(fp), continue; end
  bl=mbB.bestLabel(mbB.sessID==s); if isempty(bl), continue; end
  abc=fullfile(proj,'out','gamma','perbreath',char(s)+"__audiobook.csv");
  fbc=fullfile(proj,'out','gamma','perbreath',char(s)+"__focusedBreathing.csv");
  % read cohort/sessNum/participant from whichever CSV exists
  cc=''; sn=NaN; pp='';
  for cpath=[abc fbc]
    if isfile(cpath), tt=readtable(cpath,'TextType','string'); cc=char(tt.cohort(1)); sn=double(tt.sessNum(1)); pp=char(tt.participant(1)); break; end
  end
  k=numel(S)+1;
  S(k).sessID=char(s); S(k).finalPath=char(fp(1)); S(k).bestLabel=char(bl(1));
  S(k).cohort=cc; S(k).sessNum=sn; S(k).participant=pp;
  S(k).abCSV=char(abc); S(k).fbCSV=char(fbc);
end
nS=numel(S); if isfinite(maxSess), nS=min(nS,maxSess); S=S(1:nS); end
fprintf('sessions: %d\n', nS);

% ---------- extract (parfor over sessions) ----------
if nw>=2
  delete(gcp('nocreate')); parpool('Processes',nw);
  wait(parfevalOnAll(@() addpath(repo,codeDir), 0));
end
R=cell(nS,1);
parfor (i=1:nS, nw)
 try
  s=S(i);
  D=load(s.finalPath); vv=fieldnames(D); od=D.(vv{1}); D=[];
  labs=od.labels; ci=find(strcmpi(labs,s.bestLabel),1);
  isRsp=cellfun(@(x) ~isempty(regexpi(x,'rsp','once')), labs);
  ra=od.data(isRsp,:); ri=1; if isfield(od,'rspIDX'),ri=double(od.rspIDX);end
  rf=1; if isfield(od,'rspFlip'),rf=double(od.rspFlip);end
  sig=od.data(ci,:); rsp=ra(ri,:)*rf; N=numel(sig);
  od=[]; ra=[];   % free the big data matrix ASAP (some breathing finals are ~700 MB)
  out=struct();
  csvs={s.abCSV,s.fbCSV};
  for t=1:2
    tag=tasks{t}; cpath=csvs{t};
    sz=zeros(numel(F),numel(tMsF)); sr=zeros(1,numel(tMsF)); n=0; sc=[];
    if isfile(cpath)
      pb=readtable(cpath,'TextType','string'); pb=pb(pb.goodBreath==1 & isfinite(pb.onsetSample),:);
      for k=1:height(pb)
        eff=pb.onsetSample(k)+alignSamp; a=eff-cpre-pad; b=eff+cpost+pad;
        if a<1||b>N, continue; end
        sw=sig(a:b); rw=rsp(a:b);
        P=slt_power_cont(sw,fs,F,c1,ord,mult);
        Z=myChanZscore(P.').'; cr=(pad+1):(size(Z,2)-pad); Z=Z(:,cr); Z(~isfinite(Z))=0;
        sz=sz+Z; sr=sr+rw(cr); n=n+1;
        fmask=F>=gband(1)&F<=gband(2); tmask=tMsF>=scwin(1)&tMsF<=scwin(2);
        sc(end+1)=mean(Z(fmask,tmask),'all'); %#ok<AGROW>
      end
    end
    out.(tag)=struct('sumZ',sz,'sumRsp',sr,'n',n,'sc',sc);
  end
  % group assignment
  if strcmp(s.cohort,'OBE'), grp='control';
  elseif strcmp(s.cohort,'Dupi') && s.sessNum==1, grp='dupiS1';
  elseif strcmp(s.cohort,'Dupi') && ismember(s.sessNum,[2 3]), grp='dupiS23';
  else, grp='other'; end
  out.group=grp; out.participant=s.participant; out.sessID=s.sessID; out.sessNum=s.sessNum;
  R{i}=out;
  fprintf('  %-24s %-8s ab=%d fb=%d\n', s.sessID, grp, out.audiobook.n, out.focusedBreathing.n);
 catch ME
  R{i}=[]; fprintf('  FAIL %s: %s\n', S(i).sessID, ME.message);
 end
end
delete(gcp('nocreate'));
save(fullfile(T,'grant_breathing_R.mat'),'R','tMsF','F','tasks','S','-v7.3');  % safety: raw per-session sums

% ---------- combine into group maps + session-level scalar table ----------
groups={'all','control','dupiS1','dupiS23'};
M=struct(); sessTab=table();
for g=1:numel(groups), for t=1:2
  M.(groups{g}).(tasks{t})=struct('sumZ',zeros(numel(F),numel(tMsF)),'sumRsp',zeros(1,numel(tMsF)),'n',0,'parts',{{}});
end, end
for i=1:nS
  r=R{i}; if isempty(r)||strcmp(r.group,'other'), continue; end
  for t=1:2
    tag=tasks{t}; d=r.(tag); if d.n==0, continue; end
    for gname=[string(r.group) "all"]
      gg=char(gname);
      M.(gg).(tag).sumZ=M.(gg).(tag).sumZ+d.sumZ; M.(gg).(tag).sumRsp=M.(gg).(tag).sumRsp+d.sumRsp;
      M.(gg).(tag).n=M.(gg).(tag).n+d.n; M.(gg).(tag).parts{end+1}=r.participant;
    end
    sessTab=[sessTab; table(string(r.sessID),string(r.group),string(tag),string(r.participant),r.sessNum,d.n,mean(d.sc), ...
             'VariableNames',{'sessID','group','condition','participant','sessNum','nBreaths','gammaZ_200_400'})]; %#ok<AGROW>
  end
end
writetable(sessTab, fullfile(T,'grant_gamma200_400_session.csv'));

% save 8 group maps
for g=1:numel(groups), for t=1:2
  gg=groups{g}; tag=tasks{t}; d=M.(gg).(tag);
  meanZ=d.sumZ/max(d.n,1); meanRsp=d.sumRsp/max(d.n,1);
  nP=numel(unique(d.parts)); nTot=d.n;
  save(fullfile(T,sprintf('grant_group_meanZ_%s_%s.mat',gg,tag)),'meanZ','meanRsp','tMsF','F','nP','nTot');
end, end

% ---------- shared clim across all 8 (zero >57, floor 4) ----------
allv=[];
for g=1:numel(groups), for t=1:2
  d=M.(groups{g}).(tasks{t}); mz=d.sumZ/max(d.n,1); mz(F>57,:)=0;
  patch=mz(F>=fLim(1)&F<=fLim(2), tMsF>=-100&tMsF<=1000); allv=[allv; patch(:)]; %#ok<AGROW>
end, end
sharedUp=prctile(allv,99.5);
fprintf('\nshared clim upper=%.2f (floor 4)\n', sharedUp);

% ---------- render 8 ----------
xlab='Time from inhale (ms)';
for g=1:numel(groups), for t=1:2
  gg=groups{g}; tag=tasks{t};
  mp=fullfile(T,sprintf('grant_group_meanZ_%s_%s.mat',gg,tag));
  grant_group_render(mp, tag, [100 1000], [], fLim, 57, sharedUp, 4, gg, xlab);
end, end

% ---------- stats: group x condition on session-level 200-400 ms gamma z ----------
grant_breathing_stats(sessTab, T);
fprintf('\nDONE grant_breathing_analysis\n');
end
