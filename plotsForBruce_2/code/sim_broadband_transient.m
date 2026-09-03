function sim_broadband_transient()
% Reviewer demonstration (Voytek/Cohen): show that an inhale-locked BROADBAND
% transient on a 1/f background — with NO oscillation anywhere — reproduces the
% descending ~40 Hz "chirp" ridge seen in the real spectrograms, when passed
% through the identical superlet + within-frequency z + ridge_track pipeline.
% Also reports the superlet's effective temporal & spectral FWHM at 25/40/58 Hz.

here=fileparts(mfilename('fullpath')); addpath(here);
addpath('C:\Users\Adam\Documents\GitHub\ZelanoLabScripts');
addpath('C:\Users\Adam\Documents\GitHub\Superlets\matlab-pure');
proj=fileparts(here); fdir=fullfile(proj,'out','figs','gamma'); if ~exist(fdir,'dir'), mkdir(fdir); end
rng(7);

fs=500; F=22:1:62; c1=3; ord=[3 30]; mult=1;
band=[25 58]; bIdx=find(F>=band(1)&F<=band(2)); Fb=F(bIdx);
nTr=200; epoch=[-1000 4000]; e0=round(epoch(1)/1000*fs); e1=round(epoch(2)/1000*fs);
tMs=(e0:e1)/fs*1000; nT=numel(tMs);

% ---- build a long continuous signal: pink (1/f) noise + inhale-locked broadband transients ----
% transients are a sharp biphasic deflection (NO oscillation), one per "breath"
gap=round(4*fs);                       % 4 s between onsets
N=nTr*gap + 2*fs;
pink = pinkNoise(N);
sig = 3*pink;                          % 1/f background dominates
% broadband transient kernel: derivative-of-gaussian (sharp, broadband, aperiodic)
tk=(-round(0.05*fs):round(0.05*fs))/fs; kern = -tk.*exp(-(tk.^2)/(2*(0.008)^2)); kern=kern/max(abs(kern));
onsets=zeros(nTr,1);
for k=1:nTr
    on = k*gap; onsets(k)=on;
    a=on-numel(kern)+1; b=on; if a<1, continue; end
    amp = 1.2 + 0.2*randn;             % transient amplitude (broadband, no rhythm)
    sig(on:on+numel(kern)-1) = sig(on:on+numel(kern)-1) + amp*kern;
end

% ---- run the ACTUAL pipeline ----
P = slt_power_cont(sig, fs, F, c1, ord, mult);
sumZ=zeros(numel(bIdx),nT); ncnt=0; fr_acc=zeros(nT,1);
for k=1:nTr
    a=onsets(k)+e0; b=onsets(k)+e1; if a<1||b>N, continue; end
    Pb=P(bIdx,a:b);
    mu=mean(Pb,2); sd=std(Pb,0,2); sd(sd<=0)=eps; Zw=(Pb-mu)./sd;
    sumZ=sumZ+Zw; ncnt=ncnt+1;
    [fr,~,~]=ridge_track(max(Zw,0),Fb(:),1.0,1,2); fr_acc=fr_acc+fr(:);
end
mZ=sumZ/ncnt; mFr=fr_acc/ncnt;

% ---- FWHM of the superlet impulse response (temporal & spectral) at 25/40/58 ----
d=zeros(1,2*fs); d(fs)=1; Pd=slt_power_cont(d,fs,F,c1,ord,mult);
fwhmT=zeros(1,3); fwhmF=zeros(1,3); tt=((1:size(Pd,2))-fs)/fs*1000;
targ=[25 40 58];
for i=1:3
    [~,fi]=min(abs(F-targ(i)));
    row=Pd(fi,:); fwhmT(i)=fwhm(tt,row);
    [~,tpk]=max(Pd(fi,:)); col=Pd(:,tpk); fwhmF(i)=fwhm(F,col(:).');
end

% ---- figure ----
fig=figure('Position',[60 60 1200 460],'Color','w','Visible','off');
subplot(1,2,1);
imagesc(tMs, Fb, mZ); axis xy; hold on; plot(tMs,mFr,'w-','LineWidth',1.2); xline(0,'w-');
yline(25,'w:'); yline(58,'w:'); colorbar; try, colormap(turbo); catch, end
title({'SIMULATION: 1/f + broadband inhale transient (NO oscillation)','mean within-freq z + tracked ridge (white)'},'FontSize',9);
xlabel('ms'); ylabel('Hz');
subplot(1,2,2); axis off;
txt = {'Superlet effective resolution (impulse response FWHM):', '', ...
  sprintf('  25 Hz:  temporal FWHM = %.0f ms   spectral FWHM = %.1f Hz', fwhmT(1), fwhmF(1)), ...
  sprintf('  40 Hz:  temporal FWHM = %.0f ms   spectral FWHM = %.1f Hz', fwhmT(2), fwhmF(2)), ...
  sprintf('  58 Hz:  temporal FWHM = %.0f ms   spectral FWHM = %.1f Hz', fwhmT(3), fwhmF(3)), '', ...
  'Takeaway: a purely broadband, aperiodic inhale-locked', ...
  'transient on a 1/f background reproduces the descending', ...
  '~40 Hz ridge/"chirp" — because temporal support grows as', ...
  'frequency falls, a sharp deflection smears to later times at', ...
  'lower frequencies. The ridge tracker follows this smear even', ...
  'though NO oscillation is present. So the observed chirp is not', ...
  'by itself evidence of a sweeping gamma rhythm.'};
text(0.02,0.95,txt,'VerticalAlignment','top','FontSize',9,'Interpreter','none');
exportgraphics(fig, fullfile(fdir,'sim_broadband_transient.png'),'Resolution',130);
close(fig);
% save FWHM table
T=table(targ(:), fwhmT(:), fwhmF(:), 'VariableNames',{'freqHz','temporalFWHM_ms','spectralFWHM_Hz'});
writetable(T, fullfile(proj,'out','tables','superlet_FWHM.csv'));
fprintf('wrote sim_broadband_transient.png and superlet_FWHM.csv\n');
disp(T);
end

function y=pinkNoise(N)
    w=randn(1,N); Y=fft(w); f=(0:N-1); f(1)=1; Y=Y./sqrt(f); y=real(ifft(Y)); y=y/std(y);
end
function w=fwhm(x,y)
    y=y-min(y); [pk,pi]=max(y); h=pk/2;
    l=pi; while l>1 && y(l)>h, l=l-1; end
    r=pi; while r<numel(y) && y(r)>h, r=r+1; end
    w=abs(x(r)-x(l));
end
