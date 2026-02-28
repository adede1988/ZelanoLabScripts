function out = ica_flag_spikes_targeted(data, spikeDat,spikeProm, varargin)
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
p.addParameter('Highpass', 2, @(x) isempty(x) || (isscalar(x) && x>0));
p.addParameter('Notch', 60, @(x) isempty(x) || (isscalar(x) && x>0));
% p.addParameter('spikeDat', 21, @(x) isempty(x) || (isscalar(x) && x>0));
% p.addParameter('spikeThresh', 20, @(x) isempty(x) || (isscalar(x) && x>0));
% p.addParameter('GammaBand', [30 70], @(x)isnumeric(x)&&numel(x)==2);
% p.addParameter('HiFreqBand', [70 300], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Method', 'runica', @(s)ischar(s)||isstring(s));
p.addParameter('knownIC', [],  @(x) isempty(x) || (isscalar(x) && x>0))

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
    spikeDat= reshape(spikeDat,size(spikeDat,1), size(spikeDat,2), 1);
    spikeProm = reshape(spikeProm,size(spikeProm,1), size(spikeProm,2), 1);
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



% ------------------Prepare spike targeted snippets for ICA-----------
spikeFlags = reshape(spikeDat, size(spikeDat,1), []);               % C x (T*N)
spikeMaskData = any(spikeFlags==1, 1);                     % 1 x (T*N), union over channels
spikeProm = reshape(spikeProm, size(spikeProm,1), []); 
spikeProm = max(spikeProm, [], 1); 
spikeMaskData(spikeProm<1) = false; 

spikeIDX = find(diff(spikeMaskData)>0);
keep = false(size(spikeIDX));
last = -inf;
for i = 1:numel(spikeIDX)
    if spikeIDX(i) - last >= 100
        keep(i) = true;
        last = spikeIDX(i);
    end
end
spikeIDX = spikeIDX(keep);
spikeIDX(spikeIDX<100 | spikeIDX>length(spikeMaskData)-100) = []; 

nonSpikeIDX = sample_nonspike_onsets(spikeIDX, length(spikeMaskData), 50, 50); 

allSnips = [spikeIDX, nonSpikeIDX']; 
allSnips = randsample(allSnips, length(allSnips), false);
allSnips(allSnips<100 | allSnips>length(spikeMaskData)-100)= []; 
allSnips = arrayfun(@(x) Xtrain(:,x-50:x+50), allSnips, ...
                                            'UniformOutput', false); 
allSnips = cat(3,allSnips{:}); 

allSnips = reshape(allSnips, C, []); 


% -------------------- ICA --------------------

minSamples = 30 * C;      % rule-of-thumb lower bound
haveSamples = size(allSnips, 2);

if haveSamples < minSamples
    warning('Not enough samples for ICA: have %d, need ≥ %d (30×rank).', ...
            haveSamples, minSamples);
    out = [];     % or out = struct('ok',false); etc.
    return
end


% EEGLAB expects (channels x frames)
% Extended ICA helps for super-/sub-Gaussian comps
[weights, sphere] = runica(allSnips, 'pca', C, 'extended', 1, 'verbose', 'off');
W = weights * sphere;               % unmixing (IC x C)
Sact = W * X;                       % IC x samples
A = pinv(W);                        % mixing (C x IC)





K  = C;                       % number of ICs
Ssamples = T*N;


% -------------------- Plot all ICs & select one --------------------
%% ADDED: quick visualization and selection prompt

zic = (Sact - mean(Sact, 2)) ./ std(Sact, [], 2);

spikeSnips = arrayfun(@(x) zic(:, x-50:x+50), spikeIDX, ...
    'UniformOutput', false);
spikeSnips = cat(3, spikeSnips{:}); 
spikeSnips = squeeze(max(spikeSnips, [], 2) - min(spikeSnips, [], 2));


nonSnips = arrayfun(@(x) zic(:, x-50:x+50), nonSpikeIDX, ...
    'UniformOutput', false);
nonSnips = cat(3, nonSnips{:}); 
nonSnips = squeeze(max(nonSnips, [], 2) - min(nonSnips, [], 2));


difs = mean(spikeSnips, 2) - mean(nonSnips, 2);

% figure 
% for ic = 1:K
%         ax = subplot(K,1,ic);
% 
%         histogram(spikeSnips(ic,:)); 
%         hold on 
%         histogram(nonSnips(ic,:)); 
% 
%         ylabel(sprintf('IC %d',ic));
%         if ic==1, title('Choose spike-like IC (then close figure)'); end
%         if ic==K, xlabel('max - min around spikes'); end
%         grid(ax,'on');
% end
%     drawnow;



if (max(difs) - min(difs)) / median(difs) > 3
    [~, S.knownIC] = max(difs);
else
    S.knownIC = -1; 
end

if isempty(S.knownIC)
    figure('Name','ICA activations (z-scored)','Color','w');
 
    for ic = 1:K
        ax = subplot(K,1,ic);
        zic = (Sact(ic,:) - mean(Sact(ic,:))) / max(std(Sact(ic,:)), eps);
     
        plot(zic); xlim([1 Ssamples]);
        hold on 
        % pos = find(spikeMaskData);
        plot(spikeMaskData*3)
        ylabel(sprintf('IC %d',ic));
        if ic==1, title('Choose spike-like IC (then close figure)'); end
        if ic==K, xlabel('Samples'); end
        grid(ax,'on');
    end
    drawnow;
    
    linkaxes
    badIC = input(sprintf('Enter the index of the spike-like IC (1..%d), or [] to skip: ',K));
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
    out.mixVector = ones(1, Ssamples);
    out.badICs = [];
    return;
end


% -------------------- Cross-reference with chosen IC --------------------
%% Split chosen IC into low/high frequency and only "surgically" edit high-freq part
x_full = double(Sact(badIC,:));   % 1 x (T*N)
Ssamples = numel(x_full);

if isempty(S.Fs)
    error('Fs must be provided to split IC into low/high frequency components.');
end

Fs = S.Fs;
splitFreq = 10;          % Hz cutoff between "low" and "high" components
hpOrder   = 4;           % 4th-order Butterworth for high-pass

% Design high-pass for the spike-y part (> splitFreq)
[b_hp, a_hp] = butter(hpOrder, splitFreq/(Fs/2), 'high');

% High-frequency component of the IC
x_high = filtfilt(b_hp, a_hp, x_full(:));   % column
% Low-frequency residual (everything not captured by high-pass)
x_low  = x_full(:) - x_high;                % column

% For convenience, keep row versions for later reconstruction
x_high_row = x_high(:).';
x_low_row  = x_low(:).';

% --- Spike detection based on high-frequency component only ---
medx = median(x_high);
madx = 1.4826 * median(abs(x_high - medx));
if madx == 0
    madx = std(x_high);
end
zbadIC = abs( (x_high - medx) / max(madx, eps) );   % robust |z|, column

% Keep times where BOTH the original data had spikes AND high-freq IC is large
spikeMask = spikeMaskData & (zbadIC(:).' > 2);      % logical 1 x (T*N)

% -------------------- Build smooth surgical mixVector for high-freq part only --------------------
halfWin = 50;                             % half-width (total ~100 points)
w = hann(2*halfWin + 1)';                 % 0 at edges, 1 at center
notch = 1 - w;                            % 1 at edges, 0 at center
mixVector = ones(1, Ssamples);            % start at 1 everywhere

spikeIdx = find(spikeMask);
for ii = 1:numel(spikeIdx)
    t0 = spikeIdx(ii);
    i1 = max(1, t0 - halfWin);
    i2 = min(Ssamples, t0 + halfWin);
    % align notch segment if we clipped at edges
    n1 = 1 + (i1 - (t0 - halfWin));
    n2 = (2*halfWin + 1) - ((t0 + halfWin) - i2);
    mixVector(i1:i2) = min(mixVector(i1:i2), notch(n1:n2));
end

% -------------------- Reconstruct cleaned data --------------------
Sclean = Sact;

% Only attenuate the high-frequency portion around spikes,
% keep the low-frequency part untouched.
Sclean(badIC,:) = x_low_row + x_high_row .* mixVector;

Xclean = A * Sclean;                     % back to channel space
Xclean = Xclean + chMeans;               % restore means
data_clean = reshape(Xclean, C, T, N);   % C x T x N
data_clean = data_clean + trialMeans;    % restore trial means

% -------------------- Outputs --------------------
out.A = A; out.W = W; out.S = Sact;
out.badICs = badIC;
out.data_clean = squeeze(data_clean);
out.means = chMeans;
out.mixVector = mixVector;

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
