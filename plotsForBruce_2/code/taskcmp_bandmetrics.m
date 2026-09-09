function taskcmp_bandmetrics()
% Per breathing session x condition, from the spectro2 maps (no re-extraction):
%   theta band = mean/max over 4-8 Hz of the LOW-band z-map
%   gamma band = mean/max over 25-58 Hz of the MID-band z-map
% band-z = mean over freq (for time-to-peak + time series); band-max = max over freq
% (for the early/late window means the user asked for).
% Metrics: peakLatMs (argmax of 50-ms-smoothed band-z in [50,950]),
%          earlyBandMax (mean band-max 100-500 ms), lateBandMax (mean band-max 3000-4000 ms).
% Writes taskcmp_bandmetrics.csv (long: one row per session x condition x band) and the
% 3x2 group time-series figure (Control/DupiS1/DupiS2 x theta/gamma). CP excluded from Control.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
sdir=fullfile(proj,'out','spectro2'); tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
conds={'audiobook','focus'}; CLIP=10; EXCL_CTRL={'CP'}; fs=500;
bands={'theta','low',[4 8]; 'gamma','mid',[25 58]};

fa=dir(fullfile(sdir,'*__audiobook.mat')); tG=[]; R=table();
TC=struct();  % TC.(group).(band).(cond) = matrix subj x time ; plus *_subs
for i=1:numel(fa)
    base=strrep(fa(i).name,'__audiobook.mat','');
    for c=1:2
        f=fullfile(sdir,[base '__' tern(c==1,'audiobook','focusedBreathing') '.mat']);
        if ~exist(f,'file'), continue; end; o=load(f); o=o.out; if o.nBreaths==0, continue; end
        if isempty(tG), tG=o.tMs; end
        coh=char(string(o.cohort)); sn=o.sessNum; part=char(string(o.participant));
        if strcmpi(coh,'OBE'), grp='Control'; if any(strcmpi(part,EXCL_CTRL)), continue; end
        elseif strcmpi(coh,'Dupi'), grp=sprintf('DupiS%d',sn); else, continue; end
        for bi=1:2
            bn=bands{bi,1}; src=bands{bi,2}; rng=bands{bi,3};
            nb=o.nBreaths*o.nBaseSamp; mp=o.(src).sumPow/o.nBreaths; bmu=o.(src).baseSum/nb;
            bv=o.(src).baseSumSq/nb-bmu.^2; bv(bv<=0)=eps; z=(mp-bmu)./sqrt(bv/o.nBreaths); z=max(min(z,CLIP),-CLIP);
            fm=o.(src).freqs>=rng(1)&o.(src).freqs<=rng(2);
            bz=mean(z(fm,:),1); bmx=max(z(fm,:),[],1);
            % metrics
            sm=round(0.05*fs); bzs=movmean(bz,sm); win=tG>=50&tG<=950; yy=bzs; yy(~win)=-inf; [~,pk]=max(yy); pl=tG(pk);
            if pl<=52 || pl>=948, pl=NaN; end   % argmax hit the search-window edge => no genuine peak (invalid latency)
            eBM=mean(bmx(tG>=100&tG<=500)); lBM=mean(bmx(tG>=3000&tG<=4000));
            R=[R; table(string(o.sessID),string(part),string(coh),sn,string(grp),string(conds{c}),string(bn),pl,eBM,lBM, ...
              'VariableNames',{'sessID','participant','cohort','sessNum','group','condition','band','peakLatMs','earlyBandMax','lateBandMax'})]; %#ok<AGROW>
            % accumulate band-z tc for plot (participant-level later)
            key=matlab.lang.makeValidName([grp '_' bn '_' conds{c}]);
            if ~isfield(TC,key), TC.(key)=struct('subj',{{}},'tc',{{}}); end
            TC.(key).subj{end+1}=part; TC.(key).tc{end+1}=bz;
        end
    end
end
writetable(R, fullfile(tdir,'taskcmp_bandmetrics.csv'));
fprintf('wrote taskcmp_bandmetrics.csv (%d rows)\n', height(R));

% ---- 3x2 group time-series (participant-mean band-z) ----
groups={'Control','DupiS1','DupiS2'}; col=struct('audiobook',[.2 .4 .8],'focus',[.85 .3 .1]);
fig=figure('Position',[20 20 1350 1050],'Color','w','Visible','off'); p=0;
for gi=1:3, for bi=1:2
    p=p+1; subplot(3,2,p); hold on; bn=bands{bi,1}; g=groups{gi}; leg=[];
    for c=1:2
        key=matlab.lang.makeValidName([g '_' bn '_' conds{c}]);
        if ~isfield(TC,key), continue; end
        su=unique(TC.(key).subj,'stable'); M=nan(numel(su),numel(tG));
        for si=1:numel(su), m=strcmp(TC.(key).subj,su{si}); M(si,:)=mean(cell2mat(TC.(key).tc(m)'),1,'omitnan'); end
        mm=mean(M,1,'omitnan'); se=std(M,0,1,'omitnan')/sqrt(size(M,1));
        fill([tG fliplr(tG)],[mm+se fliplr(mm-se)],col.(conds{c}),'FaceAlpha',.15,'EdgeColor','none');
        h=plot(tG,mm,'-','Color',col.(conds{c}),'LineWidth',2); leg=[leg h]; %#ok<AGROW>
    end
    xline(0,'k-'); yline(0,'k:'); xlim([-1000 4000]);
    yl=ylim; patch([100 500 500 100],[yl(1) yl(1) yl(2) yl(2)],[.4 .8 .4],'FaceAlpha',.12,'EdgeColor','none');
    patch([3000 4000 4000 3000],[yl(1) yl(1) yl(2) yl(2)],[.9 .9 .4],'FaceAlpha',.15,'EdgeColor','none');
    title(sprintf('%s — %s band-z',g,bn),'FontWeight','bold'); ylabel([bn ' z']);
    if p==1 && numel(leg)==2, legend(leg,{'audiobook','focus'},'Location','northeast'); end
    if gi==3, xlabel('ms from inhale onset'); end; set(gca,'FontSize',10);
end, end
sgtitle('Theta & gamma band-z: audiobook vs focus (green=100-500 ms, yellow=3000-4000 ms)','FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'groups_theta_gamma_timecourse.png'),'Resolution',140); close(fig);
fprintf('wrote groups_theta_gamma_timecourse.png\n');
end
function y=tern(c,a,b), if c, y=a; else, y=b; end, end
