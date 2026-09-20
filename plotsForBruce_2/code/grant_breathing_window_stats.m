function grant_breathing_window_stats(windows)
% Recompute the session-level group x condition test on the gamma-band (25-57 Hz) z-power
% for arbitrary post-alignment windows, directly from the saved per-session z-sum maps
% (grant_breathing_R.mat) -- no re-extraction. The session mean of per-breath window-means
% equals the window-mean of the session mean map (averaging is linear), so this is exact.
% windows: cell of [w1 w2] ms re: the +200 ms alignment point. Default the three windows.
if nargin<1||isempty(windows), windows={[200 400],[400 600],[600 800]}; end
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); T=fullfile(proj,'out','tables');
L=load(fullfile(T,'grant_breathing_R.mat')); R=L.R; tMsF=L.tMsF; F=L.F; tasks=L.tasks;
gband=[25 57]; fmask=F>=gband(1)&F<=gband(2);
allTab=table();
for wi=1:numel(windows)
  W=windows{wi}; tmask=tMsF>=W(1)&tMsF<=W(2);
  st=table();
  for i=1:numel(R)
    r=R{i}; if isempty(r)||strcmp(r.group,'other'), continue; end
    for t=1:2
      tag=tasks{t}; d=r.(tag); if d.n==0, continue; end
      val=mean(d.sumZ(fmask,tmask)/d.n,'all');
      st=[st; table(string(r.sessID),string(r.group),string(tag),val,W(1),W(2), ...
        'VariableNames',{'sessID','group','condition','gammaZ','win_lo','win_hi'})]; %#ok<AGROW>
    end
  end
  allTab=[allTab; st]; %#ok<AGROW>
  fprintf('\n############## WINDOW %d-%d ms (re +200), gamma 25-57 Hz ##############\n', W(1),W(2));
  G=unique(st.group); C=["audiobook";"focusedBreathing"];
  fprintf('%-9s %-16s  n   mean     SE\n','group','condition');
  for gi=1:numel(G), for ci=1:numel(C)
    v=st.gammaZ(st.group==G(gi)&st.condition==C(ci)); v=v(isfinite(v));
    fprintf('%-9s %-16s %2d  %6.3f  %6.3f\n', G(gi),C(ci),numel(v),mean(v),std(v)/sqrt(max(numel(v),1)));
  end, end
  w=unstack(st(:,{'sessID','group','condition','gammaZ'}),'gammaZ','condition'); grp=cellstr(w.group);
  both=w(~isnan(w.audiobook)&~isnan(w.focusedBreathing),:); d=both.focusedBreathing-both.audiobook;
  [~,pC,~,stC]=ttest(both.focusedBreathing,both.audiobook);
  fprintf('condition (focus-audiobook) paired n=%d: diff=%.3f, t(%d)=%.2f, p=%.4f\n',height(both),mean(d),stC.df,stC.tstat,pC);
  mn=mean([w.audiobook w.focusedBreathing],2,'omitnan'); [pG,tG]=anova1(mn,grp,'off');
  fprintf('group effect: F(%d,%d)=%.2f, p=%.4f\n',tG{2,3},tG{3,3},tG{2,5},pG);
  [pI,tI]=anova1(d,cellstr(both.group),'off');
  fprintf('interaction (focus-audiobook by group): F(%d,%d)=%.2f, p=%.4f\n',tI{2,3},tI{3,3},tI{2,5},pI);
  pA=anovan(st.gammaZ,{cellstr(st.group),cellstr(st.condition)},'model','interaction','varnames',{'group','condition'},'display','off');
  fprintf('2-way ANOVA (approx): group p=%.4f, condition p=%.4f, interaction p=%.4f\n',pA(1),pA(2),pA(3));
end
writetable(allTab, fullfile(T,'grant_gamma_windows_session.csv'));
fprintf('\nwrote grant_gamma_windows_session.csv\n');
end
