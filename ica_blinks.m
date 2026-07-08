function out = ica_blinks(data, varargin)
% ICA_FLAG_SPIKES  Flag spike-like ICA components and reconstruct cleaned data.
%   out = ica_flag_spikes(data, 'Fs', Fs, 'Highpass',1, 'Notch',60, 'Method','auto', ...)
%
% INPUT
%   data:    C x T x N   (Channels x Time x Trials)
%   spikeDat: C X T X N (channels X time X trials) spike detections
%   spikeProm: C X T X N (channels X time X trials) spike prominence
% Name-Value (all optional)
%   'Fs'              : sampling rate (Hz). If provided, adds spectral features (gamma/hifreq ratios).
%   'Highpass'        : high-pass cutoff in Hz (default 1). Set [] to skip.
%   'Notch'           : line notch in Hz (60 or 50). Set [] to skip.
% % %   'GammaBand'       : [40 90] Hz (used if Fs provided)
% % %   'HiFreqBand'      : [100 300] Hz (used if Fs provided and Nyquist allows)
%   'Method'          : 'auto' (default), 'runica', or 'fastica'
%   --- Thresholds for flagging (tune to taste) ---
% % %   'KurtosisMin'     : 10         (components with very heavy tails)
% % %   'SparsityMin'     : 0.85       (Hoyer sparsity 0..1; spikes ~ high)
% % %   'OutlierSigma'    : 6          (#MAD for outlier definition |z|>k)
% % %   'OutlierFracMin'  : 1e-3       (fraction of samples beyond OutlierSigma*MAD)
% % %   'P2P_over_MADMin' : 20         (peak-to-peak relative to MAD)
% % %   'GammaFracMax'    : 0.5        (reject components with dominant gamma)
%
% OUTPUT (struct)
%   .A, .W, .S        : mixing, unmixing, activations (IC x samples)
%   .badICs           : indices of components flagged as spike-like
%   .features         : table of per-IC features (for inspection)
%   .Xclean           : cleaned 2-D data (C x (T*N)), bad ICs zeroed
%   .data_clean       : cleaned C x T x N array (reshaped from Xclean)
%   .X                : preprocessed input 2-D (C x (T*N)) used for ICA
%   .means            : channel means removed before ICA (added back after)
%
% NOTE: Requires EEGLAB's runica OR FastICA. If neither found, error.

% -------------------- Parse opts --------------------
p = inputParser;
p.addParameter('Fs', [], @(x) isempty(x) || (isscalar(x) && x>0));
p.addParameter('Highpass', 1, @(x) isempty(x) || (isscalar(x) && x>0));
p.addParameter('Notch', 60, @(x) isempty(x) || (isscalar(x) && x>0));
% p.addParameter('spikeDat', 21, @(x) isempty(x) || (isscalar(x) && x>0));
% p.addParameter('spikeThresh', 20, @(x) isempty(x) || (isscalar(x) && x>0));
% p.addParameter('GammaBand', [30 70], @(x)isnumeric(x)&&numel(x)==2);
% p.addParameter('HiFreqBand', [70 300], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Method', 'runica', @(s)ischar(s)||isstring(s));
p.addParameter('knownIC', [],  @(x) isempty(x) || (isscalar(x) && x>0))
p.addParameter('blinkChan', 1,  @(x) isempty(x) || (isscalar(x) && x>0))

% p.addParameter('KurtosisMin', 10, @isscalar);
% p.addParameter('SparsityMin', 0.85, @isscalar);
% p.addParameter('OutlierSigma', 6, @isscalar);
% p.addParameter('OutlierFracMin', 1e-3, @isscalar);
% p.addParameter('P2P_over_MADMin', 20, @isscalar);
% p.addParameter('GammaFracMax', 0.5, @isscalar);

p.parse(varargin{:});
S = p.Results;

% Accept 2D (C x T) or 3D (C x T x N)
origIs2D = ismatrix(data);
if origIs2D
    data    = reshape(data,    size(data,1), size(data,2), 1);
    
end


trialMeans = mean(data,2,'omitnan'); 
data = data - trialMeans;  %per trial mean center    

[C,T,N] = size(data);
X = reshape(data, C, T*N);              % C x samples
out = struct();

% -------------------- Preprocess --------------------
% remove channel DC
chMeans = mean(X,2,'omitnan');
X = X - chMeans;
Xtrain = X;
% High-pass
if ~isempty(S.Fs) && ~isempty(S.Highpass)
    [b,a] = butter(2, S.Highpass/(S.Fs/2), 'high');
    Xtrain = filtfilt(b,a, Xtrain.').';          % filter along time (columns), keep as C x samples
end

% Notch
if ~isempty(S.Fs) && ~isempty(S.Notch)
    w0 = S.Notch/(S.Fs/2);
    if w0 < 1
        bw = w0/35;                     % Q ~ 35 (adjust if needed)
        [bn,an] = iirnotch(w0, bw);
        Xtrain = filtfilt(bn,an, Xtrain.').';
    end
end

% -------------------- ICA --------------------
method = lower(string(S.Method));
hasRunica  = ~isempty(which('runica'));
hasFastica = ~isempty(which('fastica'));

if method=="auto"
    if hasRunica, method="runica";
    elseif hasFastica, method="fastica";
    else, error('Need EEGLAB runica or FastICA on path.'); end
end

switch method
case "runica"
    % EEGLAB expects (channels x frames)
    % Extended ICA helps for super-/sub-Gaussian comps
    [weights, sphere] = runica(Xtrain, 'pca', C, 'extended', 1, 'verbose', 'off');
    W = weights * sphere;               % unmixing (IC x C)
    Sact = W * X;                       % IC x samples
    A = pinv(W);                        % mixing (C x IC)
case "fastica"
    % FastICA returns S (signals) and A (mixing)
    % Note: FastICA usually whitens internally
    [Sact, A, W] = fastica(Xtrain, 'numOfIC', C, 'verbose', 'off', ...
                           'g','tanh','approach','defl');
otherwise
    error('Unknown method: %s', method);
end




K  = C;                       % number of ICs
Ssamples = T*N;



zd1 = data(S.blinkChan,:);
zd1 = (zd1 - mean(zd1)) / std(zd1); 
blinkGuess = arrayfun(@(x) corr(Sact(x,:)', data(S.blinkChan,:)'), 1:K);
plotK = sum(abs(blinkGuess)>.1);
plotidx = find(abs(blinkGuess)>.1); 
blinkGuess = abs(blinkGuess); 
zMax = (max(blinkGuess) - mean(blinkGuess)) ./ std(blinkGuess);
if zMax > 2
    [~, S.knownIC] = max(blinkGuess);
end

% Manual override for batch re-runs: force a specific blink IC via env var
% (used to reprocess a session whose auto-selection was ambiguous).
envIC = str2double(getenv('ZLP_BLINK_IC'));
if ~isnan(envIC) && envIC >= 1 && envIC <= K
    S.knownIC = round(envIC);
end
% -------------------- Plot all ICs & select one --------------------
%% ADDED: quick visualization and selection prompt
if isempty(S.knownIC)
    isBatch = batchStartupOptionUsed || ~usejava('desktop');
    if plotK > 0 && isBatch
        % Non-interactive run and no IC clearly matches the blink channel:
        % signal ambiguity so the caller can save a review plot and skip this
        % session (instead of hanging on input() or silently skipping blink
        % removal). Return the candidate activations for plotting.
        out.A = A; out.W = W; out.S = Sact;
        out.X = X; out.Xclean = A*Sact + chMeans;
        out.data_clean = reshape(out.Xclean, C, T, N) + trialMeans;
        out.means = chMeans;
        out.badICs = [];
        out.ambiguous = true;
        out.candIdx  = plotidx;
        out.candSact = Sact(plotidx,:);
        out.blinkSig = zd1;
        return;
    elseif plotK > 0
        figure('Name','ICA activations (z-scored)','Color','w');
        for ic = 1:plotK
            ax = subplot(plotK,1,ic);
            zic = (Sact(plotidx(ic),:) - mean(Sact(plotidx(ic),:))) / ...
                max(std(Sact(plotidx(ic),:)), eps);
            plot(zic); xlim([1 Ssamples]);
            hold on
            plot(zd1)

            ylabel(sprintf('IC %d',plotidx(ic)));
            if ic==1
                title('Choose blink-like IC (then close figure)');
            end
            if ic==plotK
                xlabel('Samples');
            end
            grid(ax,'on');
        end
        drawnow;
        linkaxes()
        badIC = input(sprintf(...
        'Enter the index of the blink-like IC (1..%d), or [] to skip: '...
                                                            ,K));
    else
        badIC = [];
    end
else
    badIC = S.knownIC;
end
if isempty(badIC) || badIC<1 || badIC>K
    % nothing to remove
    Sclean = Sact;
    Xclean = A * Sclean;
    Xclean = Xclean + chMeans;
    data_clean = reshape(Xclean, C, T, N);
    data_clean = data_clean + trialMeans;
    out.A = A; out.W = W; out.S = Sact;
    out.X = X; out.Xclean = Xclean;
    out.data_clean = data_clean;
    out.means = chMeans;
    out.badICs = [];
    return;
end


% -------------------- Cross-reference with chosen IC --------------------




% -------------------- Reconstruct cleaned data --------------------
Sclean = Sact;
Sclean(badIC,:) = 0; % removal of blink IC entirely 
Xclean = A * Sclean;                     % back to channel space
Xclean = Xclean + chMeans;               % restore means
data_clean = reshape(Xclean, C, T, N);   % C x T x N
data_clean = data_clean + trialMeans;    % restore trial means

% -------------------- Outputs --------------------
out.A = A; out.W = W; out.S = Sact;
out.badICs = badIC;
out.X = X; out.Xclean = Xclean;
out.data_clean = data_clean;
out.means = chMeans;

end




% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %code here to plot all 5 ICs and ask the user to select which (if any)
% %bad
% %put user index number into variable badIC: 
% badIC = 1; 
% 
% zbadIC = abs( (Sact(badIC,:) - mean(Sact(badIC,:))) / std(Sact(badIC,:))); 
% 
% 
% %threshold badIC to choose "real" spikes
% spikeIdx = spikeIdx(spikeIdx && zbadIC>3);
% 
% %create removal weight vector: 
% %start with a vector of 1s of length equal to spikeIdx
% mixVector = ones(size(spikeIdx)); 
% 
% %create a smooth deflection to zero around every spikeIdx point in the 
% %mixVector. Specifically, the points where spikeIdx is true should be 0
% %there should be a smooth gradient of points around zeros rising back 
% %to 1 across 100 time points
% mixVector = 
% 
% % -------------------- Reconstruct cleaned data --------------------
% Sclean = Sact;
% Sclean(badIC,:) = Sclean(badIC,:) .* mixVector; %surgical removal of spikes
% Xclean = A * Sclean;                     % back to channel space
% Xclean = Xclean + chMeans;               % restore means
% data_clean = reshape(Xclean, C, T, N);   % C x T x N
% data_clean = data_clean + trialMeans;            %restore trial means
% 
% 
% 
% end
