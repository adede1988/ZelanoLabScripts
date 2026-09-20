function grant_breathing_stats(sessTab, T)
% Session-level group x condition test on the 200-400 ms gamma-band (25-57 Hz) z-power.
% Groups: control / dupiS1 / dupiS23. Conditions: audiobook / focusedBreathing.
st=sessTab;
G=unique(st.group); C=["audiobook";"focusedBreathing"];
fprintf('\n==== 200-400 ms gamma z (25-57 Hz), session-level cell means ====\n');
fprintf('%-9s %-16s  n   mean     SE\n','group','condition');
for gi=1:numel(G), for ci=1:numel(C)
  v=st.gammaZ_200_400(st.group==G(gi)&st.condition==C(ci)); v=v(isfinite(v));
  fprintf('%-9s %-16s %2d  %6.3f  %6.3f\n', G(gi), C(ci), numel(v), mean(v), std(v)/sqrt(max(numel(v),1)));
end, end

% pivot to one row per session (audiobook, focusedBreathing)
w=unstack(st(:,{'sessID','group','condition','gammaZ_200_400'}),'gammaZ_200_400','condition');
hasAB=ismember('audiobook',w.Properties.VariableNames); hasFB=ismember('focusedBreathing',w.Properties.VariableNames);
if ~hasAB||~hasFB, fprintf('\n(missing a condition column; skipping paired tests)\n'); writetable(w,fullfile(T,'grant_gamma200_400_bySession.csv')); return; end
grp=cellstr(w.group);

% condition effect: paired focus - audiobook (sessions with both)
both=w(~isnan(w.audiobook)&~isnan(w.focusedBreathing),:); d=both.focusedBreathing-both.audiobook;
[~,pC,~,stC]=ttest(both.focusedBreathing, both.audiobook);
fprintf('\ncondition effect (focus - audiobook), paired n=%d: mean diff=%.3f, t(%d)=%.2f, p=%.4f\n', ...
  height(both), mean(d), stC.df, stC.tstat, pC);

% group effect: one-way ANOVA on per-session condition-mean
mn=mean([w.audiobook w.focusedBreathing],2,'omitnan');
if numel(unique(grp))>=2
  [pG,tG]=anova1(mn, grp, 'off');
  fprintf('group effect (ANOVA on per-session mean of conditions): F(%d,%d)=%.2f, p=%.4f\n', tG{2,3},tG{3,3},tG{2,5},pG);
  % interaction: (focus - audiobook) difference across groups
  if numel(unique(both.group))>=2
    [pI,tI]=anova1(d, cellstr(both.group), 'off');
    fprintf('interaction (focus-audiobook diff across groups): F(%d,%d)=%.2f, p=%.4f\n', tI{2,3},tI{3,3},tI{2,5},pI);
  end
  % 2-way ANOVA (approximate; condition treated as a fixed factor)
  pA=anovan(st.gammaZ_200_400, {cellstr(st.group), cellstr(st.condition)}, ...
    'model','interaction','varnames',{'group','condition'},'display','off');
  fprintf('2-way ANOVA (approx): group p=%.4f, condition p=%.4f, interaction p=%.4f\n', pA(1),pA(2),pA(3));
else
  fprintf('(only one group present; skipping group/interaction ANOVA)\n');
end

writetable(w, fullfile(T,'grant_gamma200_400_bySession.csv'));
fprintf('wrote grant_gamma200_400_session.csv + _bySession.csv\n');
end
