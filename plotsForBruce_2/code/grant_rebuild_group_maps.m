function grant_rebuild_group_maps()
% Rebuild the 8 breathing group meanZ maps (ALL/control/dupiS1/dupiS23 x audiobook/focused)
% from grant_breathing_R.mat -- the per-session source of truth -- WITHOUT re-extraction,
% then re-render. Use this whenever the saved grant_group_meanZ_*.mat maps are stale or
% partial (e.g. a grant_breathing_analysis run that saved R but died before rewriting the
% maps, leaving an empty control map / undercounted groups). Aggregation matches
% grant_breathing_analysis exactly.
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); addpath(codeDir);
T=fullfile(proj,'out','tables');
L=load(fullfile(T,'grant_breathing_R.mat')); R=L.R; tMsF=L.tMsF; F=L.F;
groups={'all','control','dupiS1','dupiS23'}; tasks={'audiobook','focusedBreathing'};
M=struct();
for g=1:numel(groups), for t=1:2
  M.(groups{g}).(tasks{t})=struct('sumZ',zeros(numel(F),numel(tMsF)),'sumRsp',zeros(1,numel(tMsF)),'n',0,'parts',{{}});
end, end
for i=1:numel(R), r=R{i}; if isempty(r)||strcmp(r.group,'other'), continue; end
  for t=1:2, tag=tasks{t}; d=r.(tag); if d.n==0, continue; end
    for gname=[string(r.group) "all"], gg=char(gname);
      if ~isfield(M,gg), continue; end
      M.(gg).(tag).sumZ=M.(gg).(tag).sumZ+d.sumZ; M.(gg).(tag).sumRsp=M.(gg).(tag).sumRsp+d.sumRsp;
      M.(gg).(tag).n=M.(gg).(tag).n+d.n; M.(gg).(tag).parts{end+1}=r.participant;
    end
  end
end
for g=1:numel(groups), for t=1:2
  gg=groups{g}; tag=tasks{t}; d=M.(gg).(tag);
  meanZ=d.sumZ/max(d.n,1); meanRsp=d.sumRsp/max(d.n,1); nP=numel(unique(d.parts)); nTot=d.n;
  save(fullfile(T,sprintf('grant_group_meanZ_%s_%s.mat',gg,tag)),'meanZ','meanRsp','tMsF','F','nP','nTot');
  fprintf('  rebuilt %-8s %-16s nP=%d nTot=%d\n', gg, tag, nP, nTot);
end, end
fprintf('rebuilt 8 group maps from grant_breathing_R.mat; re-rendering...\n');
grant_rerender_breathing();
end
