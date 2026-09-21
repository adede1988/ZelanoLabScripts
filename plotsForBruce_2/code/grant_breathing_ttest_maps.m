function grant_breathing_ttest_maps()
% Pixel-wise t-score maps (with BH-FDR) over the gamma-band time-frequency plane
% (25-57 Hz x -100..+1000 ms re: the +200 ms alignment point), from the saved per-session
% within-frequency z maps (grant_breathing_R.mat).
%   (1) Control: focused breathing vs audiobook  (paired, within session).
%   (2) Control (both conditions pooled) vs Dupi S1 (both conditions pooled)  (two-sample).
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); T=fullfile(proj,'out','tables'); Fg=fullfile(proj,'out');
L=load(fullfile(T,'grant_breathing_R.mat')); R=L.R; tMsF=L.tMsF; F=L.F;
fmask=F>=25 & F<=57; tmask=tMsF>=-100 & tMsF<=1000; Fb=F(fmask); tb=tMsF(tmask);

cA={}; cB={}; cC={}; d1={};   % control ab, control fb, control combined, dupiS1 combined
for i=1:numel(R)
  r=R{i}; if isempty(r), continue; end
  ab=r.audiobook; fb=r.focusedBreathing;
  comb=zeros(numel(F),numel(tMsF)); nt=0;
  if ab.n>0, comb=comb+ab.sumZ; nt=nt+ab.n; end
  if fb.n>0, comb=comb+fb.sumZ; nt=nt+fb.n; end
  if strcmp(r.group,'control')
    if ab.n>0 && fb.n>0, cA{end+1}=ab.sumZ(fmask,tmask)/ab.n; cB{end+1}=fb.sumZ(fmask,tmask)/fb.n; end %#ok<AGROW>
    if nt>0, cC{end+1}=comb(fmask,tmask)/nt; end %#ok<AGROW>
  elseif strcmp(r.group,'dupiS1')
    if nt>0, d1{end+1}=comb(fmask,tmask)/nt; end %#ok<AGROW>
  end
end
CA=cat(3,cA{:}); CB=cat(3,cB{:}); CC=cat(3,cC{:}); D1=cat(3,d1{:});
fprintf('control paired sessions=%d | control combined=%d | dupiS1 combined=%d\n', size(CA,3), size(CC,3), size(D1,3));

% ---- (1) paired: focus - audiobook (control) ----
[~,p1,~,s1]=ttest(CB,CA,'dim',3); t1=s1.tstat;
[sig1,crit1]=bh(p1(:),0.05); sig1=reshape(sig1,size(p1));
fprintf('\n(1) Control focus vs audiobook (paired, n=%d): min p=%.4g, BH crit p=%.4g, #sig(FDR .05)=%d of %d\n', ...
  size(CA,3), min(p1(:)), crit1, nnz(sig1), numel(p1));
render_tmap(tb,Fb,t1,sig1, sprintf('Control: focus - audiobook (paired, n=%d)',size(CA,3)), fullfile(Fg,'grant_tmap_ctrl_focus_vs_audiobook.png'));

% ---- (2) two-sample: control vs dupiS1 (both conditions pooled) ----
[~,p2,~,s2]=ttest2(CC,D1,'dim',3); t2=s2.tstat;
[sig2,crit2]=bh(p2(:),0.05); sig2=reshape(sig2,size(p2));
fprintf('\n(2) Control vs Dupi S1 (both cond. pooled, two-sample n=%d vs %d): min p=%.4g, BH crit p=%.4g, #sig(FDR .05)=%d of %d\n', ...
  size(CC,3), size(D1,3), min(p2(:)), crit2, nnz(sig2), numel(p2));
render_tmap(tb,Fb,t2,sig2, sprintf('Control - Dupi S1 (pooled, n=%d vs %d)',size(CC,3),size(D1,3)), fullfile(Fg,'grant_tmap_ctrl_vs_dupiS1.png'));

save(fullfile(T,'grant_ttest_maps.mat'),'t1','p1','sig1','t2','p2','sig2','Fb','tb');
fprintf('\nwrote grant_tmap_*.png + grant_ttest_maps.mat\n');
end

function [sig,critp]=bh(p,q)   % Benjamini-Hochberg FDR
p=p(:); [ps,ix]=sort(p); m=numel(ps); thr=(1:m)'/m*q; bel=ps<=thr;
if any(bel), k=find(bel,1,'last'); critp=ps(k); else, critp=0; end
sig=false(size(p)); sig(ix(ps<=critp))=true;
end

function render_tmap(tb,Fb,t,sig,ttl,outpng)
cmax=max(4, prctile(abs(t(:)),99));
f=figure('Position',[40 40 1260 520],'Color','w','Visible','off');
ax=axes(f,'Position',[0.15 0.26 0.68 0.60]);
imagesc(ax,tb,Fb,t,[-cmax cmax]); axis(ax,'xy'); colormap(ax,rb()); hold(ax,'on');
if any(sig(:)), contour(ax,tb,Fb,double(sig),[0.5 0.5],'k','LineWidth',3); end
xline(ax,0,'k--','LineWidth',2);
set(ax,'FontSize',30,'LineWidth',3.5); ylabel(ax,'Hz','FontWeight','bold','FontSize',42);
xlabel(ax,'Time from +200 ms (ms)','FontWeight','bold','FontSize',38);
cb=colorbar(ax); ax.Position=[0.15 0.26 0.68 0.60]; cb.Position=[0.845 0.26 0.024 0.60];
cb.LineWidth=3; cb.FontSize=26; cb.Label.String='t'; cb.Label.FontSize=32;
title(ax,ttl,'FontSize',24,'FontWeight','bold','Interpreter','none');
exportgraphics(f,outpng,'Resolution',200); close(f);
end

function c=rb()   % blue-white-red diverging colormap
n=256; x=linspace(0,1,n)'; r=min(1,2*x); b=min(1,2*(1-x)); g=1-abs(2*x-1); c=[r g b];
end
