function taskcmp_olf_grid()
% Parametric 2 bands x 4 windows re-slice of the cached olfactory-vs-non-olfactory
% z-maps (taskcmp_olf_blocks.mat; myChanZscore, superlet [3 30], control, CP dropped).
% NO re-extraction -- reads the stored per-block freq x time z-maps.
%   bands  : high = 40-58 Hz, low = 25-40 Hz   (same as before)
%   windows: 0-500, 500-1000, 1000-1500, 1500-2000 ms
% Writes out/tables/taskcmp_olf_roi_grid.csv (long: block x band x window) and
% redraws out/figs/taskcmp/olf_vs_nonolf_spectrograms.png with the 8-box overlay.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
L=load(fullfile(tdir,'taskcmp_olf_blocks.mat')); B=L.B; tShow=L.tShow; F=L.F;

BANDS={'high',[40 58]; 'low',[25 40]};
WINS =[0 500; 500 1000; 1000 1500; 1500 2000]; winLab={'0-500','500-1000','1000-1500','1500-2000'};

ROI=table('Size',[0 10],'VariableTypes', ...
   {'string','string','string','string','string','string','double','double','double','double'}, ...
   'VariableNames',{'sessID','participant','task','category','block','band','window','winMid','medZ','medZ_ninv'});
% attach nBreaths separately (double)
nBcol=[];
for i=1:numel(B)
  z=B(i).z; rn=sqrt(B(i).nBreaths);
  for bi=1:size(BANDS,1)
    fm = F>=BANDS{bi,2}(1) & F<=BANDS{bi,2}(2);
    for wi=1:size(WINS,1)
      tm = tShow>=WINS(wi,1) & tShow<=WINS(wi,2);
      m=median(z(fm,tm),'all');
      ROI=[ROI; table(string(B(i).sessID),string(B(i).participant),string(B(i).task), ...
           string(B(i).category),string(B(i).block),string(BANDS{bi,1}),wi,mean(WINS(wi,:)), ...
           m, m/rn,'VariableNames',ROI.Properties.VariableNames)]; %#ok<AGROW>
      nBcol(end+1,1)=B(i).nBreaths; %#ok<AGROW>
    end
  end
end
ROI.nBreaths=nBcol; ROI.window=winLab(ROI.window)';   % replace numeric index with label
writetable(ROI, fullfile(tdir,'taskcmp_olf_roi_grid.csv'));
fprintf('wrote taskcmp_olf_roi_grid.csv (%d rows = %d blocks x %d bands x %d windows)\n', ...
    height(ROI), numel(B), size(BANDS,1), size(WINS,1));

% ---- redraw spectrograms with 8-box overlay ----
olf=cat_mean(B,'olf'); non=cat_mean(B,'nonolf'); dif=olf-non;
fig=figure('Position',[20 20 1650 480],'Color','w','Visible','off');
climAbs=max([abs(prctile(olf(:),[2 98])) abs(prctile(non(:),[2 98]))]);
titles={sprintf('Olfactory (cue+thresh+O15), %d blocks',sum(strcmp({B.category},'olf'))), ...
        sprintf('Non-olfactory (audiobook+focus), %d blocks',sum(strcmp({B.category},'nonolf'))), ...
        'Difference (olfactory - non-olfactory)'};
maps={olf,non,dif};
for p=1:3
  subplot(1,3,p); imagesc(tShow, F, maps{p}); axis xy; hold on;
  ylim([25 58]); xlim([-500 2500]);
  if p<3, clim([-climAbs climAbs]); else, dc=max(abs(prctile(dif(:),[2 98]))); clim([-dc dc]); end
  colormap(gca, redblue()); colorbar;
  for wi=1:size(WINS,1)
    w=WINS(wi,:);
    rectangle('Position',[w(1) 40 diff(w) 18],'EdgeColor','k','LineWidth',1.3);  % high band
    rectangle('Position',[w(1) 25 diff(w) 15],'EdgeColor','k','LineWidth',1.3);  % low band
  end
  xline(0,'k-','LineWidth',1); yline(40,'k:');
  xlabel('time from onset (ms)','FontSize',12,'FontWeight','bold'); ylabel('Hz','FontSize',12,'FontWeight','bold');
  title(titles{p},'FontSize',12,'FontWeight','bold'); set(gca,'FontSize',11);
end
sgtitle('Control: olfactory vs non-olfactory gamma z (myChanZscore, [3 30]); boxes = 2 bands x 4 windows','FontSize',13,'FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'olf_vs_nonolf_spectrograms.png'),'Resolution',150); close(fig);
fprintf('redrew olf_vs_nonolf_spectrograms.png with 8-box overlay\n');
end
function M=cat_mean(B,catg)
idx=find(strcmp({B.category},catg)); Z=cat(3,B(idx).z); Z=min(max(Z,-10),10); M=mean(Z,3,'omitnan');
end
function c=redblue()
n=256; c=interp1([0;0.5;1], [0.23 0.30 0.75; 1 1 1; 0.75 0.20 0.20], linspace(0,1,n)');
end
