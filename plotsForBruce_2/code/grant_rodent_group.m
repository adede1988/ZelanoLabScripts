function grant_rodent_group(taskRow, winMs, climPct, fLim, scope)
% GRANT_RODENT_GROUP  Group OB spectrogram in the single-trial rodent-analog style,
% averaged over ALL breaths of a task across ALL control (OBE) participants. Each breath
% is windowed on inhale onset, superlet power is within-frequency z-scored (same as the
% single-trial figure), and the z-maps are averaged. Reports breaths/participant and n.
%
%   taskRow  'cueTask' (default) | 'focusedBreathing' | 'audiobook'
%   winMs    display window [preMs postMs]  (default [100 1000])
%   climPct  color-scale upper percentile of the displayed patch (default 99.5) — set
%            per task so each condition's gamma is highlighted on its OWN scale.
%   fLim     display frequency limits Hz (default [25 60])
if nargin<1||isempty(taskRow), taskRow='cueTask'; end
if nargin<2||isempty(winMs), winMs=[100 1000]; end
if nargin<3||isempty(climPct), climPct=99.5; end
if nargin<4||isempty(fLim), fLim=[25 60]; end
if nargin<5||isempty(scope), scope='control'; end   % 'control' (OBE only) | 'all' (OBE+Dupi)
proj='C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\plotsForBruce_2';
repo='C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
addpath(repo); addpath(fullfile(proj,'code'));
T=fullfile(proj,'out','tables'); outdir=fullfile(proj,'out');
idx=readtable(fullfile(T,'session_index.csv'),'TextType','string');
mb =readtable(fullfile(T,'macbp_best.csv'),'TextType','string');
mbTask = taskRow; if any(strcmp(taskRow,{'focusedBreathing','audiobook'})), mbTask='breathingTask'; end

fs=500; F=(15:1:150).'; c1=3; ord=[3 30]; mult=1;
compPre=1000; compPost=3000; padMs=600;
cpre=round(compPre/1000*fs); cpost=round(compPost/1000*fs); pad=round(padMs/1000*fs);
tMsF=(-cpre:cpost)/fs*1000; nT=numel(tMsF);

glob = sprintf('*OBE*%s.csv',taskRow); if strcmp(scope,'all'), glob=sprintf('*%s.csv',taskRow); end
pbfiles=dir(fullfile(proj,'out','gamma','perbreath',glob));
sumZ=zeros(numel(F),nT); sumRsp=zeros(1,nT); nTot=0;
perP=containers.Map('KeyType','char','ValueType','double');
fprintf('=== extracting group %s (control) ===\n', taskRow);
for i=1:numel(pbfiles)
    pb=readtable(fullfile(pbfiles(i).folder,pbfiles(i).name),'TextType','string');
    pb=pb(pb.goodBreath==1 & isfinite(pb.onsetSample),:);
    if isempty(pb), continue; end
    s=char(pb.sessID(1)); part=char(pb.participant(1));
    fp=idx.finalPath(idx.sessID==string(s) & idx.task==string(mbTask)); if isempty(fp), continue; end
    D=load(char(fp(1))); vv=fieldnames(D); od=D.(vv{1}); clear D;
    labs=od.labels; blc=mb.bestLabel(mb.sessID==string(s) & mb.task==string(mbTask)); if isempty(blc), clear od; continue; end
    bl=char(blc(1)); ci=find(strcmpi(labs,bl),1); if isempty(ci), clear od; continue; end
    isRsp=cellfun(@(x) ~isempty(regexpi(x,'rsp','once')), labs);
    ra=od.data(isRsp,:); ri=1; if isfield(od,'rspIDX'),ri=double(od.rspIDX);end
    rf=1; if isfield(od,'rspFlip'),rf=double(od.rspFlip);end
    sig=od.data(ci,:); rsp=ra(ri,:)*rf; N=numel(sig); clear od;
    nsess=0;
    for k=1:height(pb)
        on=pb.onsetSample(k); a=on-cpre-pad; b=on+cpost+pad;
        if a<1 || b>N, continue; end
        sw=sig(a:b); rw=rsp(a:b);
        P=slt_power_cont(sw,fs,F,c1,ord,mult);
        Z=myChanZscore(P.').'; cr=(pad+1):(size(Z,2)-pad); Z=Z(:,cr); Z(~isfinite(Z))=0;
        sumZ=sumZ+Z; sumRsp=sumRsp+rw(cr); nTot=nTot+1; nsess=nsess+1;
    end
    if ~isKey(perP,part), perP(part)=0; end; perP(part)=perP(part)+nsess;
    fprintf('  %-22s (%s): %d breaths\n', s, part, nsess);
    clear sig rsp;
end
meanZ=sumZ/max(nTot,1); meanRsp=sumRsp/max(nTot,1);

ks=keys(perP); nP=numel(ks);
fprintf('\n==== GROUP: %s %s ====\n', scope, taskRow);
fprintf('participants: %d | total breaths: %d | mean breaths/participant: %.1f\n', nP, nTot, nTot/max(nP,1));
fprintf('breaths per participant:\n'); for j=1:nP, fprintf('  %-6s %d\n', ks{j}, perP(ks{j})); end
Rt=table(string(ks(:)), cell2mat(values(perP))', 'VariableNames',{'participant','nBreaths'});
writetable(Rt, fullfile(T,sprintf('grant_group_breathcounts_%s_%s.csv',scope,taskRow)));
matPath=fullfile(T,sprintf('grant_group_meanZ_%s_%s.mat',scope,taskRow));
save(matPath,'meanZ','meanRsp','tMsF','F','nP','nTot','taskRow','scope');

grant_group_render(matPath, taskRow, winMs, climPct, fLim, [], [], [], scope);
end
