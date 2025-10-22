function [P, f] = trial_psd(data, Fs, varargin)
% TRIAL_PSD  Welch PSD per trial
%   [P,f] = TRIAL_PSD(data, Fs, 'WindowSec',1, 'Overlap',0.5, 'NFFT',[])
%   data: C x T x N   (Channels x Time x Trials)
%   Fs:   sampling rate (Hz)
%
%   Outputs:
%     P : C x F x N    (power spectral density [power/Hz], Welch-averaged within each trial)
%     f : F x 1        (frequencies in Hz)
%
%   Name-value options:
%     'WindowSec' (default = min(1, T/Fs))  length of Welch window in seconds
%     'Overlap'   (default = 0.5)           fraction overlap (0..<1)
%     'NFFT'      (default = nextpow2(windowLen))  FFT length
%     'Demean'    (default = true)          remove per-channel mean before PSD

% ----- parse inputs -----
p = inputParser;
p.addParameter('WindowSec', [], @(x) isempty(x) || (isscalar(x) && x>0));
p.addParameter('Overlap',   0.5, @(x) isscalar(x) && x>=0 && x<1);
p.addParameter('NFFT',      [], @(x) isempty(x) || (isscalar(x) && x>=1));
p.addParameter('Demean',    true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

[C,T,N] = size(data);
if isempty(opt.WindowSec)
    opt.WindowSec = min(1, T/Fs);  % default: 1 s if trial is long enough, else whole trial
end
winLen   = max(16, round(opt.WindowSec * Fs));       % at least 16 samples
noverlap = min(winLen-1, round(opt.Overlap * winLen));
if isempty(opt.NFFT)
    opt.NFFT = 2^nextpow2(winLen);
end
win = hamming(winLen, 'periodic');

% ----- compute PSD per trial (columns=channels) -----
P = []; f = [];
for tr = 1:N
    % Make it Time x Channels for pwelch (each column = one channel)
    X = squeeze(data(:,:,tr)).';    % T x C
    if opt.Demean
        X = X - mean(X,1,'omitnan');
    end
    % Welch PSD (one-sided for real signals), Pxx: F x C
    [Pxx, f] = pwelch(X, win, noverlap, opt.NFFT, Fs, 'psd'); 
    if tr == 1
        P = zeros(C, size(Pxx,1), N, 'like', Pxx);  % C x F x N
    end
    P(:,:,tr) = Pxx.';  % transpose to C x F
end
end
