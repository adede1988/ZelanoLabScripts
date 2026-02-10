function [] = singleChanPipeline(chanFiles, filei, subFiles, codePre)

%% set frequency parameters

    
frex = logspace(log10(.1),log10(200),300);
% bandWidth = logspace(log10(.1), log10(20), 300);
numfrex = length(frex); 
stds = logspace(log10(3),log10(15),numfrex);

standardPhaseLen = 50; 

% highfrex = linspace(70, 150, 81); 
% highnumfrex = length(highfrex); 
% highstds = logspace(log10(10),log10(20),highnumfrex)./(2*pi*highfrex);

%% load the data
folderName = chanFiles(filei).folder; 
stem = strsplit(folderName, 'CHANDAT'); 
stem = stem{1};
try %try loading the processed file
    
    chanDat = load([stem 'CHANDAT_processed/' chanFiles(filei).name]).chanDat; 
    QC = load([chanDat.QCFileDir '/' chanDat.QCFileName]).cleanDat;
catch
    chanDat = load([chanFiles(filei).folder '/' chanFiles(filei).name]).chanDat; % go raw if it's not working!
    QC = load([chanDat.QCFileDir '/' chanDat.QCFileName]).cleanDat;
end

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])


%% hardCode timing relative to sniff

%hard code that breath window is [-2:10) 

chanDat.tim = -2:1/chanDat.fs:...
    10 - 1/chanDat.fs;
relSamp = round(chanDat.tim * chanDat.fs);   % sample offsets relative to onset
nTime   = numel(chanDat.tim);
onsetSamp = chanDat.behDat.finalOnset;
onsetSamp = onsetSamp(:);
nTrials = numel(onsetSamp);
nSamp = numel(chanDat.data);

lowPassFilt = designfilt('lowpassiir','FilterOrder',4, ...
            'HalfPowerFrequency',2+1.7,'SampleRate',2000); 

lowRsp = filtfilt(lowPassFilt, chanDat.rsp);

trialDat   = nan(nTrials, nTime);
trialRsp   = nan(nTrials, nTime);
trialLowRsp = nan(nTrials, nTime);
trialSpike = nan(nTrials, nTime);
trialBlink = nan(nTrials, nTime);
trialBad   = nan(nTrials, nTime);

if strcmp(chanDat.task, 'breathingTask')
    trialRR   = nan(nTrials, nTime);
    trialTarg = nan(nTrials, nTime);
end

% --- Build trials
for tr = 1:nTrials
    idx = onsetSamp(tr) + relSamp;                 % absolute sample indices
    ok  = (idx >= 1) & (idx <= nSamp);             % guard (should mostly be all-true)

    if any(ok)
        trialDat(tr, ok)   = chanDat.data(idx(ok));
        trialRsp(tr, ok)   = chanDat.rsp(idx(ok));
        trialLowRsp(tr, ok)   = lowRsp(idx(ok));

        trialSpike(tr, ok) = QC.spikeCleanVec(idx(ok));
        trialBlink(tr, ok) = QC.blinkCleanVec(idx(ok));
        trialBad(tr, ok)   = QC.badTS(idx(ok));

        if strcmp(chanDat.task, 'breathingTask')
            trialRR(tr, ok)   = chanDat.RRint(idx(ok));
            trialTarg(tr, ok) = chanDat.targTrace(idx(ok));
        end
    end
end

% (optional) stash back into structs for convenience
chanDat.trial.tim  = chanDat.tim;
chanDat.trial.data = trialDat;
chanDat.trial.rsp  = trialRsp;
chanDat.trial.lowRsp = trialLowRsp; 

QC.trial.spikeCleanVec = trialSpike;
QC.trial.blinkCleanVec = trialBlink;
QC.trial.badTS         = trialBad;

if strcmp(chanDat.task, 'breathingTask')
    chanDat.trial.RRint     = trialRR;
    chanDat.trial.targTrace = trialTarg;
end



%% breath phase info per trial: 

if ~isfield(chanDat, 'targIDX')
   
    
    idx50 = breathPiecewiseTemplateIdx(chanDat);         % nBreaths x 50
    chanDat.targIDX = idx50; 
    
   

end

%% QC use/notuse breath by breath

use = ones(size(chanDat.trial.data,1),1); 
reasonEliminate = zeros(size(chanDat.trial.data,1), 4); 
%columns are reasons for eliminating: 
%col 1: bad breath
%col 2: couldn't template match
%col 3: bad EEG
%col 4: blink
for ii = 1:length(use)
    if chanDat.behDat.goodBreath(ii) == 0
        use(ii) = 0;
        reasonEliminate(ii,1) = 1; 
    end
    if sum(isnan(chanDat.targIDX(ii,:)))>10
        use(ii) = 0; 
        reasonEliminate(ii,2) = 1; 
    end
    startIdx = 500; 
    endIdx = round(chanDat.behDat.length(ii)+1000); 
    L = endIdx - startIdx; 
    if sum(QC.trial.badTS(ii, startIdx:endIdx)) / L > .25
        use(ii) = 0; 
        reasonEliminate(ii,3) = 1; 
    end

    if sum(QC.trial.blinkCleanVec(ii, startIdx:endIdx)) / L > .25
        use(ii) = 0; 
        reasonEliminate(ii,4) = 1; 
    end


    

end

chanDat.use = use; 
chanDat.reasonEliminate = reasonEliminate; 


%% long time series power/phase extraction
%this code needs to do everything for one frequency at a time to avoid
%problems with memory overload 
disp('working on tf extraction')
if ~isfield(chanDat, 'tf')
    %store into 50 timepoint samples per breath
    outPowZ = nan([size(chanDat.targIDX) numfrex]);
    outPow = nan([size(chanDat.targIDX) numfrex]);
    outPhase = nan([size(chanDat.targIDX) numfrex]);
    outBreathSeg = nan([size(chanDat.targIDX)]);
    %breaths X phases X frequencies 
    spectra = nan([size(chanDat.targIDX,1) 5 numfrex]);
    for fi = 1:length(frex)
        fi
      % extract power and phase information
        [phase,pow] = multiphasevec3(frex(fi),...
                        chanDat.data,chanDat.fs,stds(fi), true);
        pow = squeeze(pow); 
        phase = squeeze(phase); 
        % cut to trials: 
        curTrialPow = nan(size(chanDat.trial.data)); 
        curTrialPhase = curTrialPow; 
    
        for tr = 1:nTrials
            idx = onsetSamp(tr) + relSamp;      % absolute sample indices
            ok  = (idx >= 1) & (idx <= nSamp);  % guard (should mostly be all-true)
        
            if any(ok)
                curTrialPow(tr, ok)   = pow(idx(ok));
                curTrialPhase(tr, ok)   = phase(idx(ok));
               
            end
        end
        curTrialPow = curTrialPow'; 
        curTrialPhase = curTrialPhase';
        %z-score using only QC passed trials in the standardization
        powz = myChanZscore(curTrialPow, ...
            [find(chanDat.tim>=-.450,1), find(chanDat.tim>=-.050,1)],...
            use); %z-score
    
        for tt = 1:size(outPowZ, 1)
            idx = chanDat.targIDX(tt,:); 
            idx(isnan(idx)) = []; 
            outPow(tt, 1:length(idx),fi) = curTrialPow(idx,tt); 
            outPowZ(tt,1:length(idx),fi) = powz(idx,tt);
            outPhase(tt, 1:length(idx),fi) = curTrialPhase(idx,tt); 
            outBreathSeg(tt,1:length(idx)) = chanDat.trial.rsp(tt,idx); 
        end
    
    
        %get power spectra in each respiratory phase
         for tt = 1: size(spectra,1)
            if ~any(isnan(chanDat.targIDX(tt,1:10)))
                idx = chanDat.targIDX(tt,1:10); 
                spectra(tt, 1, fi) = ...
                    mean(curTrialPow(idx, tt), 1);
            end
            if ~any(isnan(chanDat.targIDX(tt,11:20)))
                idx = chanDat.targIDX(tt,11:20); 
                spectra(tt, 2, fi) = ...
                    mean(curTrialPow(idx, tt), 1);
            end
            if ~any(isnan(chanDat.targIDX(tt,21:30)))
                idx = chanDat.targIDX(tt,21:30); 
                spectra(tt, 3, fi) = ...
                    mean(curTrialPow(idx, tt, :), 1);
            end
            if ~any(isnan(chanDat.targIDX(tt,31:40)))
                idx = chanDat.targIDX(tt,31:40); 
                spectra(tt, 4, fi) = ...
                    mean(curTrialPow(idx, tt, :), 1);
            end
            if ~any(isnan(chanDat.targIDX(tt,41:50)))
                idx = chanDat.targIDX(tt,41:50); 
                spectra(tt, 5, fi) = ...
                    mean(curTrialPow(idx, tt, :), 1);
            end
        end
    
    % ---- save per-frequency trial-resolved pow/phase (drop-in) ----
    saveDir = fullfile(stem,'CHANDAT_processed','tf_files');
    if fi==1 && ~exist(saveDir,'dir'), mkdir(saveDir); end
    
    [~,baseName] = fileparts(chanFiles(filei).name); % strips .mat
    
    tfMeta = struct( ...
        'subID',chanDat.subID,'task',chanDat.task,'tim',chanDat.tim,'fs',chanDat.fs,'use',chanDat.use, ...
        'targIDX',chanDat.targIDX,'chi',chanDat.chi,'chanType',chanDat.chanType,'behDat',chanDat.behDat, ...
        'sessID',chanDat.sessID,'sessNum',chanDat.sessNum,'baseEmotion',chanDat.baseEmotion, ...
        'OGdataDir',chanDat.OGdataDir,'fi',fi,'frex',frex(fi));
    tfMeta.labels = chanDat.labels; 
    
    tfInfo = tfMeta; 
    tfInfo.pow   = curTrialPow;
    save(fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi)), 'tfInfo', '-v7.3');
    
    tfInfo = tfMeta; 
    tfInfo.phase = curTrialPhase;
    save(fullfile(saveDir, sprintf('%s_phase_fi%03d.mat', baseName, fi)), 'tfInfo', '-v7.3');
    
    
    
    
    end
    %store for saving out
    chanDat.tf.phase = outPhase; 
    chanDat.tf.pow = outPow; 
    chanDat.tf.powZ = outPowZ; 
    chanDat.tf.breathSeg = outBreathSeg; 
    chanDat.tf.frex = frex; 
    chanDat.tf.spectra = spectra; 
    
    clear trialBad trialBlink trialDat trialRR trialRsp trialSpike ...
        trialTarg ok onsetSamp idx tr ff cleanFields nTime nSamp ...
        nTrials relSamp lowRsp
    saveDir = fullfile(stem,'CHANDAT_processed'); 
    save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3'); 
end


%% remove long timeseries data, no longer needed
if isfield(chanDat, 'data')
    cleanFields = {'data', 'targTrace', 'RRint', 'rsp','spikeCleanVec', ...
        'blinkCleanVec', 'badTS'};
    for ff = 1:length(cleanFields)
        try
            chanDat = rmfield(chanDat, cleanFields{ff}); 
        catch
        end
        try
            QC = rmfield(QC, cleanFields{ff}); 
        catch
        end
    end
end

    

%% fooof
disp('working on fooof')
if ~isfield(chanDat, 'fooof')
    %--- Aperiodic subtraction per breath × epoch WITH KNEE ---
    % Requires: spectra (nBreaths x 5 x numfrex), frex (numfrex x 1 or 1 x numfrex)
    spectra = chanDat.tf.spectra; 
    frex = frex(:);
    
    % ----- USER SETTINGS -----
    fitRangeHz = [.5 150];
    excludeLineHz = [60 120];        % [] to disable
    excludeHalfWidthHz = 2;
    minFitPoints = 12;               
    
    fitMask = frex >= fitRangeHz(1) & frex <= fitRangeHz(2);
    if ~isempty(excludeLineHz)
        for f0 = excludeLineHz(:)'
            fitMask = fitMask & ~(frex > (f0 - excludeHalfWidthHz) & frex < (f0 + excludeHalfWidthHz));
        end
    end
    
    % Outputs
    spectra_flat_log10 = nan(size(spectra));              % log10(power) residual
    aperiodic_log10    = nan(size(spectra));              % fitted aperiodic in log10(power)
    aperiodic_params   = nan([size(spectra,1) size(spectra,2) 3]); % [offset exponent knee]
    
    % Knee aperiodic model in log10 space:
    % y = offset - log10(knee + f^exponent)
    % We parameterize exponent = exp(logExp) and knee = 10^(klog) to enforce positivity.
    aper_model = @(th,f) th(1) - log10( (10.^th(3)) + (f.^exp(th(2))) );
    
    useLSQ = exist('lsqcurvefit','file') == 2;
    
    if useLSQ
        opts = optimoptions('lsqcurvefit','Display','off');
        % Reasonable soft bounds in transformed space:
        lb = [-Inf, log(1e-3), -12];   % exponent >= 1e-3, knee >= 1e-12
        ub = [ Inf, log(50),   12];   % exponent <= 50,   knee <= 1e12
    end
    
    for tt = 1:size(spectra,1) %trials 
        for ee = 1:size(spectra,2) % epochs
    
            p = squeeze(spectra(tt,ee,:));  % numfrex x 1
            if all(isnan(p)), continue; end
    
            p(p <= 0) = eps;               % guard
            y = log10(p);
    
            m = fitMask & ~isnan(y);
            if nnz(m) < minFitPoints, continue; end
    
            % --- initial guesses from straight-line (no-knee) fit in log-log ---
            xf = log10(frex(m));
            cf = polyfit(xf, y(m), 1);         % y ≈ slope*log10(f) + intercept
            slope0    = cf(1);
            offset0   = cf(2);
            exponent0 = max(1e-3, -slope0);    % knee model exponent is positive
    
            fmin = min(frex(m));
            knee0 = max(1e-12, (fmin.^exponent0) / 10);  % small knee to start
    
            th0 = [offset0, log(exponent0), log10(knee0)];
    
            % --- fit ---
            if useLSQ
                th = lsqcurvefit(@(th,f) aper_model(th,f), th0, frex(m), y(m), lb, ub, opts);
            else
                obj = @(th) sum((aper_model(th, frex(m)) - y(m)).^2);
                th  = fminsearch(obj, th0, optimset('Display','off'));
            end
    
            yfit = aper_model(th, frex);
    
            offset   = th(1);
            exponent = exp(th(2));
            knee     = 10.^th(3);
    
            aperiodic_log10(tt,ee,:)    = yfit;
            spectra_flat_log10(tt,ee,:) = y - yfit;
            aperiodic_params(tt,ee,:)   = [offset exponent knee];
        end
    end
    
    % Optional linear-space versions:
    aperiodic_fit = 10.^aperiodic_log10;     % fitted 1/f(+knee) component in power units
    spectra_flat  = 10.^spectra_flat_log10;  % equals spectra ./ aperiodic_fit (ratio)
    
    f_knee = aperiodic_params(:,:,3).^(1 ./ aperiodic_params(:,:,2));  % Hz
    
    
    chanDat.fooof.aperiodic_fit = aperiodic_fit; 
    chanDat.fooof.spectra_flat = spectra_flat; 
    chanDat.fooof.aperiodic_params = aperiodic_params; 
    chanDat.fooof.f_knee = f_knee; 
    
    % --- Gamma peak frequency per breath × epoch (25–60 Hz) ---
    gammaMask = frex >= 25 & frex <= 60;
    
    tmp = spectra_flat(:,:,gammaMask);   % [nBreaths x 5 x nGamma]
    allNan = all(isnan(tmp),3);          % breaths/epochs missing
    
    tmp(isnan(tmp)) = -Inf;              % so max ignores NaNs
    [~, imax] = max(tmp, [], 3);         % argmax within gamma band
    gammaFreqs = frex(gammaMask);        % Hz values for gamma bins
    
    gamma_peak_freq = gammaFreqs(imax);  % [nBreaths x 5]
    gamma_peak_freq(allNan) = NaN;
    
    chanDat.fooof.gamma_peak_freq = gamma_peak_freq;

    % --- Low-frequency peak per breath × epoch (3–14 Hz; theta/alpha) ---
    lowMask = frex >= .5 & frex <= 14;
    
    tmp = spectra_flat(:,:,lowMask);     % [nBreaths x 5 x nLow]
    allNan = all(isnan(tmp),3);          % breaths/epochs missing
    
    tmp(isnan(tmp)) = -Inf;              % so max ignores NaNs
    [~, imax] = max(tmp, [], 3);         % argmax within low band
    lowFreqs = frex(lowMask);            % Hz values for low bins
    
    low_peak_freq = lowFreqs(imax);      % [nBreaths x 5]
    low_peak_freq(allNan) = NaN;
    
    chanDat.fooof.low_peak_freq = low_peak_freq;

    
    figure;
    plot(squeeze(mean(spectra_flat, 1,'omitnan'))')  % 5 x numfrex -> plotted as 5 lines
    xticks(20:20:300)
    xticklabels(round(frex(20:20:300)))

    legend({'inhale rise','inhale fall','exhale rise','exhale fall','pause'}, ...
           'Location','best','Box','off');
    
    
    % --- Peak frequency (Hz) per breath, per epoch, within 30–50 Hz ---
    % Inputs assumed in workspace:
    %   spectra_flat : [nBreaths x 5 x numfrex]
    %   frex         : [numfrex x 1] or [1 x numfrex]
    saveDir = fullfile(stem,'CHANDAT_processed'); 
    save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');
end



%% gamma bursts =======================
% STEP A: Gamma-burst event detection (1 peak per breath)
% Uses gamma_peak_freq(:,2) as the breath-wise gamma frequency (Hz)
% Loads the corresponding *_pow_fi###.mat file per breath, smooths with 100 ms FWHM Gaussian,
% then finds peak within [inhaleOnset-200ms, endOfBreath].
% Stores timing + QC metrics into chanDat.gammaBurst
%% =======================
disp('working on gammaBurst')
if ~isfield(chanDat, 'gammaBurst')
gamma_peak_freq = chanDat.fooof.gamma_peak_freq;
% Inputs assumed in workspace:
%   gamma_peak_freq : [nTrials x 5] (Hz)
%   frex            : [numfrex x 1] (Hz)
%   chanDat.behDat.length : [nTrials x 1] (sec)
%   chanFiles(filei).name, stem, filei
%   chanDat.fs (500), chanDat.tim (optional time vector aligned so tim(1001)=0)

fs = chanDat.fs;
nTrials = size(gamma_peak_freq,1);

% where the tf files live
[~, baseName] = fileparts(chanFiles(filei).name);
saveDir = fullfile(stem,'CHANDAT_processed','tf_files');

% constants
inhaleOnsetIdx = 1001;
preMs = 800;
startIdxFixed = inhaleOnsetIdx - round((preMs/1000)*fs);  % 200 ms pre-inhale
startIdxFixed = max(1, startIdxFixed);

% smoothing kernel: 100 ms FWHM gaussian
fwhm_sec  = 0.100;
fwhm_samp = fwhm_sec * fs;
sigma_samp = fwhm_samp / (2*sqrt(2*log(2))); % FWHM -> sigma
halfK = max(1, ceil(4*sigma_samp));
xk = (-halfK:halfK)';
gk = exp(-(xk.^2) ./ (2*sigma_samp^2));
gk = gk ./ sum(gk);

% allocate outputs
gammaBurst = struct();
gammaBurst.freqHz          = nan(nTrials,1);
gammaBurst.fi              = nan(nTrials,1);
gammaBurst.t0_idx          = nan(nTrials,1);
gammaBurst.t0_sec          = nan(nTrials,1);
gammaBurst.peak            = nan(nTrials,1);
gammaBurst.prominence      = nan(nTrials,1);
gammaBurst.width_samp      = nan(nTrials,1);
gammaBurst.width_sec       = nan(nTrials,1);
gammaBurst.edgeDist_samp   = nan(nTrials,1);
gammaBurst.edgeDist_sec    = nan(nTrials,1);
gammaBurst.snr             = nan(nTrials,1);

% (optional, but useful for debugging/QC later)
gammaBurst.searchStart_idx = nan(nTrials,1);
gammaBurst.searchEnd_idx   = nan(nTrials,1);
gammaBurst.baseline        = nan(nTrials,1);

% metadata about how detection was done
gammaBurst.smooth_fwhm_sec  = fwhm_sec;
gammaBurst.smooth_fwhm_samp = fwhm_samp;
gammaBurst.inhaleOnsetIdx   = inhaleOnsetIdx;
gammaBurst.search_preMs     = preMs;

% breath-wise gamma frequency (use ONLY epoch 2)
gammaFreqVec = gamma_peak_freq(:,2);

for tr = 1:nTrials
tr
    f = gammaFreqVec(tr);
    if isnan(f), continue; end
    gammaBurst.freqHz(tr) = f;

    % map frequency value -> index in frex (used in filename)
    idxExact = find(frex == f, 1, 'first');
    if ~isempty(idxExact)
        fi = idxExact;
    else
        % fallback: nearest frex bin (handles floating precision issues)
        [d, fi] = min(abs(frex - f));
        % (no thresholding requested; just proceed)
    end
    gammaBurst.fi(tr) = fi;

    % load the corresponding power file (raw power)
    powFile = fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi));
    S = load(powFile, 'tfInfo');
    powMat = S.tfInfo.pow;   % [6000 x nTrials]

    if tr > size(powMat,2)
        continue
    end

    powTS = powMat(:,tr);    % [6000 x 1]
    nT = numel(powTS);

    % smooth power time series with gaussian (100 ms FWHM)
    powSmooth = conv(powTS, gk, 'same');

    % define breath end index and search window
    endIdx = round(chanDat.behDat.length(tr)*fs + inhaleOnsetIdx);
    endIdx = min(nT, max(1, endIdx));

    startIdx = startIdxFixed;
    startIdx = min(startIdx, endIdx);  % guard

    gammaBurst.searchStart_idx(tr) = startIdx;
    gammaBurst.searchEnd_idx(tr)   = endIdx;

    seg = powSmooth(startIdx:endIdx);
    if all(isnan(seg)), continue; end

    % peak in smoothed data within search window
    [pk, relIdx] = max(seg);               % relIdx is within seg
    t0_idx = startIdx + relIdx - 1;

    gammaBurst.t0_idx(tr) = t0_idx;

    % time in seconds (prefer chanDat.tim if it exists and matches)
    if isfield(chanDat,'tim') && numel(chanDat.tim) >= t0_idx
        gammaBurst.t0_sec(tr) = chanDat.tim(t0_idx);
    else
        gammaBurst.t0_sec(tr) = (t0_idx - inhaleOnsetIdx) / fs;
    end

    gammaBurst.peak(tr) = pk;

    % baseline & prominence (baseline = median in search segment)
    base = median(seg, 'omitnan');
    gammaBurst.baseline(tr)   = base;
    gammaBurst.prominence(tr) = pk - base;

    % width (samples) at half-prominence above baseline, measured on smoothed seg
    if isfinite(base) && isfinite(pk) && pk > base
        halfLevel = base + 0.5*(pk - base);

        % expand left/right from the peak until dropping below halfLevel (NaNs treated as below)
        L = relIdx; R = relIdx;
        while L > 1
            if ~isfinite(seg(L-1)) || seg(L-1) < halfLevel, break; end
            L = L - 1;
        end
        while R < numel(seg)
            if ~isfinite(seg(R+1)) || seg(R+1) < halfLevel, break; end
            R = R + 1;
        end

        w_samp = R - L + 1;
    else
        w_samp = NaN;
    end

    gammaBurst.width_samp(tr) = w_samp;
    gammaBurst.width_sec(tr)  = w_samp / fs;

    % edge distance (how close peak is to either boundary of search window)
    distLeft  = relIdx - 1;
    distRight = numel(seg) - relIdx;
    ed_samp = min(distLeft, distRight);

    gammaBurst.edgeDist_samp(tr) = ed_samp;
    gammaBurst.edgeDist_sec(tr)  = ed_samp / fs;

    % SNR: (peak - baseline) / robust noise SD within search window
    resid = seg - base;
    resid = resid(isfinite(resid));
    if ~isempty(resid)
        noiseStd = 1.4826 * mad(resid, 1);         % robust SD estimate
        if noiseStd == 0
            noiseStd = std(resid, 0);              % fallback
        end
        if noiseStd > 0
            gammaBurst.snr(tr) = (pk - base) / noiseStd;
        else
            gammaBurst.snr(tr) = NaN;
        end
    end

end

chanDat.gammaBurst = gammaBurst;

  saveDir = fullfile(stem,'CHANDAT_processed'); 
    save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');

end




end














% % extract power and phase information
    % [phase,pow] = multiphasevec3(frex,chanDat.trial.data,chanDat.fs,6, 0);
    % pow = permute(pow, [3,1,2]);
    % powz = arrayfun(@(x) myChanZscore(pow(:,:,x), ...
    %     [find(chanDat.tim>=-.450,1), find(chanDat.tim>=-.050,1)],...
    %     use),...
    %     1:size(pow,3), 'UniformOutput',false ); %z-score
    % 
    % powz = cell2mat(powz); %organize
    % powz = reshape(powz, size(powz,1), size(powz,2)/numfrex, []); %organize
    % 
    % %store into 50 timepoint samples per breath
    % outPowZ = nan([size(chanDat.targIDX) numfrex]);
    % outPow = nan([size(chanDat.targIDX) numfrex]);
    % outPhase = nan([size(chanDat.targIDX) numfrex]);
    % outBreathSeg = nan([size(chanDat.targIDX)]);
    % for tt = 1:size(outPowZ, 1)
    %     idx = chanDat.targIDX(tt,:); 
    %     idx(isnan(idx)) = []; 
    %     outPow(tt, 1:length(idx),:) = pow(idx,tt,:); 
    %     outPowZ(tt,1:length(idx),:) = powz(idx,tt,:);
    %     outPhase(tt, 1:length(idx),:) = squeeze(phase(tt,:,idx))'; 
    %     outBreathSeg(tt,1:length(idx)) = chanDat.trial.rsp(tt,idx); 
    % end
    % 
    % %store for saving out
    % chanDat.phase = outPhase; 
    % chanDat.pow = outPow; 
    % chanDat.powZ = outPowZ; 
    % chanDat.breathSeg = outBreathSeg; 
    % chanDat.frex = frex; 
    

    
    % figure; imagesc(squeeze(mean(powz(500:end-500,use==1,:),2, 'omitnan'))')
    % set(gca, 'ydir', 'normal')
    % yticks([20:20:200])
    % yticklabels(round(frex(20:20:200)))
    % colorbar
    % caxis([-2,2])
    % yyaxis right
    % plot(mean(chanDat.trial.rsp(use==1,500:end-500),1))
    % 
    % figure; imagesc(squeeze(mean(outPowZ(use==1,:,:),1, 'omitnan'))')
    % set(gca, 'ydir', 'normal')
    % yticks([20:20:300])
    % yticklabels(round(frex(20:20:300)))
    % colorbar
    % caxis([-5,5])
    % yyaxis right
    % plot(mean(outBreathSeg(use==1,:),1, 'omitnan'))
    % 
    % figure; imagesc(squeeze(mean(outPow(use==1,:,:),1, 'omitnan'))')
    % set(gca, 'ydir', 'normal')
    % yticks([20:20:200])
    % yticklabels(round(frex(20:20:200)))
    % colorbar
    % caxis([0,200000])
    % yyaxis right
    % plot(mean(outBreathSeg(use==1,:),1, 'omitnan'))
    % 
    % figure; imagesc(squeeze(abs(mean(exp(1i.*outPhase(use==1,:,:)),1, 'omitnan')))')
    % set(gca, 'ydir', 'normal')
    % yticks([20:20:200])
    % yticklabels(round(frex(20:20:200)))
    % colorbar
    % caxis([0,.2])
    % yyaxis right
    % plot(mean(outBreathSeg(use==1,:),1, 'omitnan'))

    
    %create per breath per phase power spectra for peak analysis: 
    %breaths X phases X frequencies 
    % spectra = nan([size(chanDat.targIDX,1) 5 numfrex]);
    % for tt = 1: size(spectra,1)
    %     if ~any(isnan(chanDat.targIDX(tt,1:10)))
    %         idx = chanDat.targIDX(tt,1:10); 
    %         spectra(tt, 1, :) = mean(pow(idx, tt, :), 1);
    %     end
    %     if ~any(isnan(chanDat.targIDX(tt,11:20)))
    %         idx = chanDat.targIDX(tt,11:20); 
    %         spectra(tt, 2, :) = mean(pow(idx, tt, :), 1);
    %     end
    %     if ~any(isnan(chanDat.targIDX(tt,21:30)))
    %         idx = chanDat.targIDX(tt,21:30); 
    %         spectra(tt, 3, :) = mean(pow(idx, tt, :), 1);
    %     end
    %     if ~any(isnan(chanDat.targIDX(tt,31:40)))
    %         idx = chanDat.targIDX(tt,31:40); 
    %         spectra(tt, 4, :) = mean(pow(idx, tt, :), 1);
    %     end
    %     if ~any(isnan(chanDat.targIDX(tt,41:50)))
    %         idx = chanDat.targIDX(tt,41:50); 
    %         spectra(tt, 5, :) = mean(pow(idx, tt, :), 1);
    %     end
    % end

