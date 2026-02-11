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
if ~isfield(chanDat, 'trial')
    chanDat.tim = -2:1/chanDat.fs:...
        10 - 1/chanDat.fs;
    % sample offsets relative to onset
    relSamp = round(chanDat.tim * chanDat.fs);   
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
        idx = onsetSamp(tr) + relSamp;    % absolute sample indices
        % guard (should mostly be all-true)
        ok  = (idx >= 1) & (idx <= nSamp); 
    
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
end


%% breath phase info per trial: 

if ~isfield(chanDat, 'targIDX')
   
    
    idx50 = breathPiecewiseTemplateIdx(chanDat);         % nBreaths x 50
    chanDat.targIDX = idx50; 
    
   

end

%% QC use/notuse breath by breath

if ~isfield(chanDat, 'use')
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
end

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
    
    % ---- save per-frequency trial-resolved pow/phase (DOWNSAMPLED to 100 Hz) ----
    saveDir = fullfile(stem,'CHANDAT_processed','tf_files');
    if fi==1 && ~exist(saveDir,'dir'), mkdir(saveDir); end

    [~,baseName] = fileparts(chanFiles(filei).name); % strips .mat

    % === Downsample for disk output only ===
    fs_full = chanDat.fs;     % original sampling rate of curTrialPow/Phase
    fs_tf   = 100;            % desired saved sampling rate

    ds = fs_full / fs_tf;
    if abs(ds - round(ds)) > 1e-12
        error('Downsample requires integer factor: fs_full=%g, fs_tf=%g', fs_full, fs_tf);
    end
    ds = round(ds);

    % Downsample along time dimension (rows). This is toolbox-free and keeps alignment.
    powSave   = curTrialPow(1:ds:end, :);      % [nTime_ds x nTrials]
    phaseSave = curTrialPhase(1:ds:end, :);
    timSave   = chanDat.tim(1:ds:end);

    % Optional: cut file size ~2x more (highly recommended)
    powSave   = single(powSave);
    phaseSave = single(phaseSave);
    timSave   = single(timSave);

    % Save metadata (note fs/tim are for the SAVED files)
    tfMeta = struct( ...
        'subID',chanDat.subID,'task',chanDat.task,'tim',timSave,'fs',fs_tf,'use',chanDat.use, ...
        'targIDX',chanDat.targIDX,'chi',chanDat.chi,'chanType',chanDat.chanType,'behDat',chanDat.behDat, ...
        'sessID',chanDat.sessID,'sessNum',chanDat.sessNum,'baseEmotion',chanDat.baseEmotion, ...
        'OGdataDir',chanDat.OGdataDir,'fi',fi,'frex',frex(fi), ...
        'fs_full',fs_full,'dsFactor',ds);

    tfMeta.labels = chanDat.labels;

    % Save pow
    tfInfo = tfMeta;
    tfInfo.pow = powSave;
    save(fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi)), 'tfInfo', '-v7.3');

    % Save phase
    tfInfo = tfMeta;
    tfInfo.phase = phaseSave;
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
if ~isfield(chanDat, 'gammaBuasdrst')

gamma_peak_freq = chanDat.fooof.gamma_peak_freq;
nTrials = size(gamma_peak_freq,1);

[~, baseName] = fileparts(chanFiles(filei).name);
saveDir = fullfile(stem,'CHANDAT_processed','tf_files');

% breath-wise gamma frequency (use ONLY epoch 2)
gammaFreqVec = gamma_peak_freq(:,2);

% --- Preload one tf pow file to get tf sampling + timebase ---
tr0 = find(~isnan(gammaFreqVec), 1, 'first');
if isempty(tr0)
    warning('No valid gamma peak freqs; skipping gammaBurst.');
    chanDat.gammaBurst = struct();
    return
end

% map tr0 freq -> fi
f0 = gammaFreqVec(tr0);
idxExact = find(frex == f0, 1, 'first');
if isempty(idxExact)
    [~, fi0] = min(abs(frex - f0));
else
    fi0 = idxExact;
end

powFile0 = fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi0));
S0 = load(powFile0, 'tfInfo');

fs_tf  = double(S0.tfInfo.fs);     % SHOULD be 100 now
tfTim  = double(S0.tfInfo.tim(:)); % time vector for saved files
dsFact = double(S0.tfInfo.dsFactor);

% inhale onset index = first sample at/after t=0 in tfTim
inhaleOnsetIdx = find(tfTim >= 0, 1, 'first');
if isempty(inhaleOnsetIdx)
    % fallback if tfTim doesn't cross 0 for some reason
    inhaleOnsetIdx = round(abs(tfTim(1))*fs_tf) + 1;
end

% constants
preMs = 800;  % (your code uses 800; comment said 200 earlier)
startIdxFixed = inhaleOnsetIdx - round((preMs/1000)*fs_tf);
startIdxFixed = max(1, startIdxFixed);

% smoothing kernel: 100 ms FWHM gaussian (in tf samples)
fwhm_sec   = 0.100;
fwhm_samp  = fwhm_sec * fs_tf;
sigma_samp = fwhm_samp / (2*sqrt(2*log(2)));
halfK = max(1, ceil(4*sigma_samp));
xk = (-halfK:halfK)';
gk = exp(-(xk.^2) ./ (2*sigma_samp^2));
gk = gk ./ sum(gk);

% allocate outputs
gammaBurst = struct();
gammaBurst.freqHz          = nan(nTrials,1);
gammaBurst.fi              = nan(nTrials,1);

% Indices are now in TF (100 Hz) sample units
gammaBurst.t0_idx          = nan(nTrials,1);
gammaBurst.t0_sec          = nan(nTrials,1);

% Optional mapping back to original sampling (e.g., 500 Hz) indices
gammaBurst.t0_idx_full     = nan(nTrials,1);

gammaBurst.peak            = nan(nTrials,1);
gammaBurst.prominence      = nan(nTrials,1);
gammaBurst.width_samp      = nan(nTrials,1);
gammaBurst.width_sec       = nan(nTrials,1);
gammaBurst.edgeDist_samp   = nan(nTrials,1);
gammaBurst.edgeDist_sec    = nan(nTrials,1);
gammaBurst.snr             = nan(nTrials,1);
gammaBurst.searchStart_idx = nan(nTrials,1);
gammaBurst.searchEnd_idx   = nan(nTrials,1);
gammaBurst.baseline        = nan(nTrials,1);

% metadata
gammaBurst.fs_tf           = fs_tf;
gammaBurst.dsFactor        = dsFact;
gammaBurst.smooth_fwhm_sec = fwhm_sec;
gammaBurst.inhaleOnsetIdx  = inhaleOnsetIdx;
gammaBurst.search_preMs    = preMs;

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
        [~, fi] = min(abs(frex - f));
    end
    gammaBurst.fi(tr) = fi;

    % load the corresponding power file (now saved at fs_tf)
    powFile = fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi));
    S = load(powFile, 'tfInfo');
    powMat = S.tfInfo.pow;   % [nTime_tf x nTrials]

    if tr > size(powMat,2), continue; end

    powTS = double(powMat(:,tr));   % [nTime_tf x 1]
    nT = numel(powTS);

    % smooth
    powSmooth = conv(powTS, gk, 'same');
   
    % define breath end index and search window (convert sec -> tf samples)
    endIdx = round(chanDat.behDat.length(tr)*fs_tf + inhaleOnsetIdx);
    endIdx = min(nT, max(1, endIdx));

    startIdx = min(startIdxFixed, endIdx);

    gammaBurst.searchStart_idx(tr) = startIdx;
    gammaBurst.searchEnd_idx(tr)   = endIdx;

    seg = powSmooth(startIdx:endIdx);
    if all(isnan(seg)), continue; end

    [pk, relIdx] = max(seg);
    t0_idx = startIdx + relIdx - 1;

    gammaBurst.t0_idx(tr) = t0_idx;
    gammaBurst.t0_sec(tr) = tfTim(t0_idx);   % IMPORTANT: use tf file tim

    % map back to full-rate index if you want it
    % (works when dsFactor is integer; here it should be 5)
    gammaBurst.t0_idx_full(tr) = (t0_idx - 1)*dsFact + 1;

    gammaBurst.peak(tr) = pk;

    base = median(seg, 'omitnan');
    gammaBurst.baseline(tr)   = base;
    gammaBurst.prominence(tr) = pk - base;

    % width at half-prominence (in tf samples)
    if isfinite(base) && isfinite(pk) && pk > base
        halfLevel = base + 0.5*(pk - base);
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
    gammaBurst.width_sec(tr)  = w_samp / fs_tf;

    distLeft  = relIdx - 1;
    distRight = numel(seg) - relIdx;
    ed_samp = min(distLeft, distRight);

    gammaBurst.edgeDist_samp(tr) = ed_samp;
    gammaBurst.edgeDist_sec(tr)  = ed_samp / fs_tf;

    resid = seg - base;
    resid = resid(isfinite(resid));
    if ~isempty(resid)
        noiseStd = 1.4826 * mad(resid, 1);
        if noiseStd == 0, noiseStd = std(resid, 0); end
        if noiseStd > 0
            gammaBurst.snr(tr) = (pk - base) / noiseStd;
        end
    end
end

% --- Assign respiratory phase (1..5) by distance to nearest phase centroid ---
% targIDX is in 500 Hz sample indices; convert to TF (100 Hz) index space
targIDX_tf = floor((chanDat.targIDX - 1) ./ dsFact) + 1;  % dsFact should be 5
gammaBurst.phaseID = nan(nTrials,1);
gammaBurst.phaseDist_samp = nan(nTrials,1);               % distance (TF samples) to chosen centroid

for tr = 1:nTrials
    t0 = gammaBurst.t0_idx(tr);
    if isnan(t0), continue; end

    cent = nan(1,5);
    for ph = 1:5
        cols = (ph-1)*10 + (1:10);
        idxs = targIDX_tf(tr, cols);
        idxs = idxs(~isnan(idxs));
        if ~isempty(idxs)
            cent(ph) = mean(idxs);   % centroid of that phase's indices (in TF samples)
        end
    end

    if all(isnan(cent)), continue; end

    [dmin, phBest] = min(abs(cent - t0), [], 'omitnan');
    gammaBurst.phaseID(tr) = phBest;
    gammaBurst.phaseDist_samp(tr) = dmin;
end



chanDat.gammaBurst = gammaBurst;

saveDir = fullfile(stem,'CHANDAT_processed');
save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');

end



%% =========================================================
% Secondary gamma-burst surrogate:
% highest peak at least 500 ms away from the primary peak (no phaseID logic)
% Produces gammaBurstSecondary with same per-trial fields as gammaBurst
%% =========================================================

disp('working on gammaBurstSecondary')

[~, baseName] = fileparts(chanFiles(filei).name);
saveDir = fullfile(stem,'CHANDAT_processed','tf_files');

gammaBurst = chanDat.gammaBurst;  % primary
nTrials = numel(gammaBurst.freqHz);

% Rebuild smoothing kernel from stored metadata (TF sample units)
fs_tf   = gammaBurst.fs_tf;          % e.g., 100
dsFact  = gammaBurst.dsFactor;       % e.g., 5
inhaleOnsetIdx = gammaBurst.inhaleOnsetIdx;

preMs = gammaBurst.search_preMs;
startIdxFixed = inhaleOnsetIdx - round((preMs/1000)*fs_tf);
startIdxFixed = max(1, startIdxFixed);

fwhm_sec   = gammaBurst.smooth_fwhm_sec;
fwhm_samp  = fwhm_sec * fs_tf;
sigma_samp = fwhm_samp / (2*sqrt(2*log(2)));
halfK = max(1, ceil(4*sigma_samp));
xk = (-halfK:halfK)';
gk = exp(-(xk.^2) ./ (2*sigma_samp^2));
gk = gk ./ sum(gk);

% Exclusion radius: 500 ms away from primary peak
exclSec  = 0.500;
exclSamp = ceil(exclSec * fs_tf);   % inclusive exclusion zone

% Initialize secondary struct with same layout/metadata as primary
gammaBurstSecondary = gammaBurst;

% Reset per-trial fields that will be recomputed
resetFields = {'t0_idx','t0_sec','t0_idx_full','peak','prominence','width_samp','width_sec', ...
               'edgeDist_samp','edgeDist_sec','snr','searchStart_idx','searchEnd_idx','baseline'};

for ff = 1:numel(resetFields)
    if isfield(gammaBurstSecondary, resetFields{ff})
        gammaBurstSecondary.(resetFields{ff})(:) = NaN;
    end
end

% (Optional) If you want phaseID for secondary too and your primary struct has it,
% we will recompute it below using nearest-centroid. If you don't care, it will stay as-is.
if isfield(gammaBurstSecondary,'phaseID'), gammaBurstSecondary.phaseID(:) = NaN; end
if isfield(gammaBurstSecondary,'phaseDist_samp'), gammaBurstSecondary.phaseDist_samp(:) = NaN; end

% For optional phase assignment (nearest centroid), convert targIDX to TF indices
targIDX_tf = floor((chanDat.targIDX - 1) ./ dsFact) + 1;  % breaths x 50

for tr = 1:nTrials
tr
    % Need a valid primary peak and a valid frequency index for this breath
    t0_primary = gammaBurst.t0_idx(tr);
    fi = gammaBurst.fi(tr);
    if isnan(t0_primary) || isnan(fi), continue; end

    % Load the corresponding power file (inefficient but clean)
    powFile = fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi));
    S = load(powFile, 'tfInfo');
    powMat = S.tfInfo.pow;   % [nTime_tf x nTrials]
    if tr > size(powMat,2), continue; end

    powTS = double(powMat(:,tr));
    nT = numel(powTS);

    % Smooth
    powSmooth = conv(powTS, gk, 'same');

    % Define breath end index and search window (TF sample units)
    endIdx = round(chanDat.behDat.length(tr)*fs_tf + inhaleOnsetIdx);
    endIdx = min(nT, max(1, endIdx));
    startIdx = min(startIdxFixed, endIdx);

    gammaBurstSecondary.searchStart_idx(tr) = startIdx;
    gammaBurstSecondary.searchEnd_idx(tr)   = endIdx;

    seg = powSmooth(startIdx:endIdx);
    if all(isnan(seg)), continue; end

    % Mask out ±500 ms around the primary peak (in TF samples), within the search window
    seg2 = seg;
    relPrimary = round(t0_primary - startIdx + 1);  % primary location within seg (1..length(seg))
    if relPrimary >= 1 && relPrimary <= numel(seg)
        exLo = max(1, relPrimary - exclSamp);
        exHi = min(numel(seg), relPrimary + exclSamp);
        seg2(exLo:exHi) = NaN;
    end

    if all(isnan(seg2)), continue; end

    % Find highest remaining peak (secondary)
    [pk, relIdx] = max(seg2);
    if isnan(pk) || isempty(relIdx), continue; end

    t0_idx = startIdx + relIdx - 1;

    gammaBurstSecondary.t0_idx(tr)      = t0_idx;
    gammaBurstSecondary.t0_sec(tr)      = double(S.tfInfo.tim(t0_idx));
    gammaBurstSecondary.t0_idx_full(tr) = (t0_idx - 1)*dsFact + 1;

    gammaBurstSecondary.peak(tr) = pk;

    % Baseline/prominence computed on allowed samples (excluding masked region)
    base = median(seg2, 'omitnan');
    gammaBurstSecondary.baseline(tr)   = base;
    gammaBurstSecondary.prominence(tr) = pk - base;

    % Width at half-prominence (on allowed samples; NaNs stop width)
    if isfinite(base) && isfinite(pk) && pk > base
        halfLevel = base + 0.5*(pk - base);
        L = relIdx; R = relIdx;
        while L > 1
            if ~isfinite(seg2(L-1)) || seg2(L-1) < halfLevel, break; end
            L = L - 1;
        end
        while R < numel(seg2)
            if ~isfinite(seg2(R+1)) || seg2(R+1) < halfLevel, break; end
            R = R + 1;
        end
        w_samp = R - L + 1;
    else
        w_samp = NaN;
    end

    gammaBurstSecondary.width_samp(tr) = w_samp;
    gammaBurstSecondary.width_sec(tr)  = w_samp / fs_tf;

    % Edge distance relative to full search window (same definition as primary)
    distLeft  = relIdx - 1;
    distRight = numel(seg) - relIdx;
    ed_samp = min(distLeft, distRight);

    gammaBurstSecondary.edgeDist_samp(tr) = ed_samp;
    gammaBurstSecondary.edgeDist_sec(tr)  = ed_samp / fs_tf;

    % SNR on allowed samples
    resid = seg2 - base;
    resid = resid(isfinite(resid));
    if ~isempty(resid)
        noiseStd = 1.4826 * mad(resid, 1);
        if noiseStd == 0, noiseStd = std(resid, 0); end
        if noiseStd > 0
            gammaBurstSecondary.snr(tr) = (pk - base) / noiseStd;
        end
    end

    % Optional: assign respiratory phase for secondary peak by nearest centroid
    if isfield(gammaBurstSecondary,'phaseID')
        cent = nan(1,5);
        for ph = 1:5
            cols = (ph-1)*10 + (1:10);
            idxs = targIDX_tf(tr, cols);
            idxs = idxs(~isnan(idxs));
            if ~isempty(idxs), cent(ph) = mean(idxs); end
        end
        if ~all(isnan(cent))
            [dmin, phBest] = min(abs(cent - t0_idx), [], 'omitnan');
            gammaBurstSecondary.phaseID(tr) = phBest;
            if isfield(gammaBurstSecondary,'phaseDist_samp')
                gammaBurstSecondary.phaseDist_samp(tr) = dmin;
            end
        end
    end

end

chanDat.gammaBurstSecondary = gammaBurstSecondary;


saveDir = fullfile(stem,'CHANDAT_processed');
save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');



%% =========================================================
% STEP C (REWRITE): Frequency-first build of event-locked tensors + per-frequency null files
%   - Primary/Secondary tensors saved into chanDat.gammaLockTF
%   - Shuffled nulls saved as ONE FILE PER FREQUENCY:
%       fullfile(stem,'CHANDAT')processed'','burstNullFiles',[baseName sprintf('_BurstNull_fi%03d.mat',fi)])
%   - Each file contains:
%       nullOut.primary.powNull   [nShuf x nBreaths x 21] single
%       nullOut.primary.phaseNull [nShuf x nBreaths x 21] single
%       nullOut.secondary.powNull [nShuf x nBreaths x 21] single
%       nullOut.secondary.phaseNull [nShuf x nBreaths x 21] single
%       nullOut.meta  (small metadata)
%% =========================================================

disp('working on gammaLockTF (Step C) - per-frequency null files')

% ---- basic checks ----
if ~isfield(chanDat,'gammaBurst') || ~isfield(chanDat.gammaBurst,'t0_idx')
    error('chanDat.gammaBurst.t0_idx not found.');
end
if ~isfield(chanDat,'gammaBurstSecondary') || ~isfield(chanDat.gammaBurstSecondary,'t0_idx')
    error('chanDat.gammaBurstSecondary.t0_idx not found.');
end

t0_primary   = double(chanDat.gammaBurst.t0_idx(:));          % [nBreaths x 1] TF (100 Hz) indices
t0_secondary = double(chanDat.gammaBurstSecondary.t0_idx(:)); % [nBreaths x 1] TF (100 Hz) indices
nBreaths     = numel(t0_primary);

[~, baseName] = fileparts(chanFiles(filei).name);
saveDirTf   = fullfile(stem,'CHANDAT_processed','tf_files');
saveDirNull = fullfile(stem,'CHANDAT_processed','burstNullFiles');
if ~exist(saveDirNull,'dir'), mkdir(saveDirNull); end

% ---- frequency downsampling: take every 3rd frequency (300 -> 100) ----
frexSelIdx = 1:3:numel(frex);
frexSel    = frex(frexSelIdx);
nFsel      = numel(frexSelIdx);

% ---- time grid: 100 Hz TF -> 20 Hz window centered on t0 (decimation) ----
fs_tf  = double(chanDat.gammaBurst.fs_tf);   % should be 100
fs_out = 20;
ds_out = fs_tf / fs_out;
if abs(ds_out - round(ds_out)) > eps
    error('fs_tf (%g) not divisible by fs_out (%g).', fs_tf, fs_out);
end
ds_out = round(ds_out);  % should be 5

offsets_out = (-10:10);               % 21 points at 20 Hz
offsets_tf  = offsets_out * ds_out;   % e.g., -50..+50 in 100 Hz samples
centerIdx   = find(offsets_out==0,1,'first');

% ---- breath boundaries (pad outside breath with NaN) ----
inhaleOnsetIdx = double(chanDat.gammaBurst.inhaleOnsetIdx);  % TF index where t=0
breathStartIdx = inhaleOnsetIdx * ones(nBreaths,1);

breathEndIdx = round(double(chanDat.behDat.length(:)) * fs_tf + inhaleOnsetIdx);

% ---- allocate tensors in chanDat.gammaLockTF (store as single to reduce size) ----
gammaLockTF = struct();
gammaLockTF.fs_tf       = fs_tf;
gammaLockTF.fs_out      = fs_out;
gammaLockTF.ds_out      = ds_out;
gammaLockTF.win_sec     = 0.5;
gammaLockTF.offsets_out = offsets_out;
gammaLockTF.offsets_tf  = offsets_tf;
gammaLockTF.centerIdx   = centerIdx;

gammaLockTF.frexSelIdx  = frexSelIdx;
gammaLockTF.frexSel     = frexSel;

gammaLockTF.t0_primary   = t0_primary;
gammaLockTF.t0_secondary = t0_secondary;

gammaLockTF.pow_primary     = single(nan(nBreaths, numel(offsets_out), nFsel));
gammaLockTF.phase_primary   = single(nan(nBreaths, numel(offsets_out), nFsel));
gammaLockTF.pow_secondary   = single(nan(nBreaths, numel(offsets_out), nFsel));
gammaLockTF.phase_secondary = single(nan(nBreaths, numel(offsets_out), nFsel));

% ---- shuffle setup ----
nShuf = 1000;
useVec = double(chanDat.use(:)) == 1;

poolPrim = find(useVec & isfinite(t0_primary));
poolSec  = find(useVec & isfinite(t0_secondary));

if isempty(poolPrim)
    warning('No usable breaths for primary shuffle null (chanDat.use==1 & finite t0).');
end
if isempty(poolSec)
    warning('No usable breaths for secondary shuffle null (chanDat.use==1 & finite t0).');
end

% Precompute permutation maps (lightweight) so shuffles are consistent across frequencies
% permPrim(ss,:) gives indices into poolPrim (same length each time)
permPrim = [];
permSec  = [];
if ~isempty(poolPrim)
    permPrim = zeros(nShuf, numel(poolPrim), 'uint32');
    for ss = 1:nShuf
        permPrim(ss,:) = uint32(randperm(numel(poolPrim)));
    end
end
if ~isempty(poolSec)
    permSec = zeros(nShuf, numel(poolSec), 'uint32');
    for ss = 1:nShuf
        permSec(ss,:) = uint32(randperm(numel(poolSec)));
    end
end

% helper for fast sub2ind mapping
breathIdxMat = repmat((1:nBreaths)', 1, numel(offsets_tf)); % [nBreaths x 21]

% ======================================================
% MAIN LOOP: frequency (load each TF file once)
% ======================================================
for kk = 1:nFsel
    fi = frexSelIdx(kk);
    fprintf('gammaLockTF: %d/%d (orig fi=%d, f=%.3f Hz)\n', kk, nFsel, fi, frex(fi));

    % load pow + phase for this frequency
    powFile   = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat',   baseName, fi));
    phaseFile = fullfile(saveDirTf, sprintf('%s_phase_fi%03d.mat', baseName, fi));

    Sp = load(powFile,   'tfInfo');
    Sh = load(phaseFile, 'tfInfo');

    powMat   = double(Sp.tfInfo.pow);    % [nTime_tf x nBreaths]
    phaseMat = double(Sh.tfInfo.phase);  % [nTime_tf x nBreaths]
    nT       = size(powMat,1);

    breathEndClamped = min(max(1, breathEndIdx), nT);

    % ---------- PRIMARY extraction ----------
    idxMat = bsxfun(@plus, t0_primary, offsets_tf); % [nBreaths x 21]
    valid = isfinite(idxMat) & idxMat>=1 & idxMat<=nT & ...
            bsxfun(@ge, idxMat, breathStartIdx) & bsxfun(@le, idxMat, breathEndClamped);

    idxSafe = idxMat; idxSafe(~valid) = 1;
    lin = sub2ind([nT nBreaths], idxSafe(valid), breathIdxMat(valid));

    winPow = nan(nBreaths, numel(offsets_tf));
    winPh  = nan(nBreaths, numel(offsets_tf));
    tmpPow = powMat(:);
    tmpPh  = phaseMat(:);
    winPow(valid) = tmpPow(lin);
    winPh(valid)  = tmpPh(lin);

    gammaLockTF.pow_primary(:,:,kk)   = single(winPow);
    gammaLockTF.phase_primary(:,:,kk) = single(winPh);

    % ---------- SECONDARY extraction ----------
    idxMat = bsxfun(@plus, t0_secondary, offsets_tf);
    valid = isfinite(idxMat) & idxMat>=1 & idxMat<=nT & ...
            bsxfun(@ge, idxMat, breathStartIdx) & bsxfun(@le, idxMat, breathEndClamped);

    idxSafe = idxMat; idxSafe(~valid) = 1;
    lin = sub2ind([nT nBreaths], idxSafe(valid), breathIdxMat(valid));

    winPow = nan(nBreaths, numel(offsets_tf));
    winPh  = nan(nBreaths, numel(offsets_tf));
    winPow(valid) = tmpPow(lin);
    winPh(valid)  = tmpPh(lin);

    gammaLockTF.pow_secondary(:,:,kk)   = single(winPow);
    gammaLockTF.phase_secondary(:,:,kk) = single(winPh);

    % ======================================================
    % Build shuffled nulls IN MEMORY for this frequency, then save once
    % ======================================================
    powNull_primary   = single(nan(nShuf, nBreaths, numel(offsets_tf)));
    phaseNull_primary = single(nan(nShuf, nBreaths, numel(offsets_tf)));

    powNull_secondary   = single(nan(nShuf, nBreaths, numel(offsets_tf)));
    phaseNull_secondary = single(nan(nShuf, nBreaths, numel(offsets_tf)));

    % ---- PRIMARY shuffles ----
    if ~isempty(poolPrim)
        for ss = 1:nShuf
            t0s = nan(nBreaths,1);
            t0s(poolPrim) = t0_primary(poolPrim(double(permPrim(ss,:)))); % permute within pool

            idxMat = bsxfun(@plus, t0s, offsets_tf);
            valid = isfinite(idxMat) & idxMat>=1 & idxMat<=nT & ...
                    bsxfun(@ge, idxMat, breathStartIdx) & bsxfun(@le, idxMat, breathEndClamped);

            idxSafe = idxMat; idxSafe(~valid) = 1;
            lin = sub2ind([nT nBreaths], idxSafe(valid), breathIdxMat(valid));

            wPow = nan(nBreaths, numel(offsets_tf));
            wPh  = nan(nBreaths, numel(offsets_tf));
            wPow(valid) = tmpPow(lin);
            wPh(valid)  = tmpPh(lin);

            powNull_primary(ss,:,:)   = single(wPow);
            phaseNull_primary(ss,:,:) = single(wPh);
        end
    end

    % ---- SECONDARY shuffles ----
    if ~isempty(poolSec)
        for ss = 1:nShuf
            t0s = nan(nBreaths,1);
            t0s(poolSec) = t0_secondary(poolSec(double(permSec(ss,:)))); % permute within pool

            idxMat = bsxfun(@plus, t0s, offsets_tf);
            valid = isfinite(idxMat) & idxMat>=1 & idxMat<=nT & ...
                    bsxfun(@ge, idxMat, breathStartIdx) & bsxfun(@le, idxMat, breathEndClamped);

            idxSafe = idxMat; idxSafe(~valid) = 1;
            lin = sub2ind([nT nBreaths], idxSafe(valid), breathIdxMat(valid));

            wPow = nan(nBreaths, numel(offsets_tf));
            wPh  = nan(nBreaths, numel(offsets_tf));
            wPow(valid) = tmpPow(lin);
            wPh(valid)  = tmpPh(lin);

            powNull_secondary(ss,:,:)   = single(wPow);
            phaseNull_secondary(ss,:,:) = single(wPh);
        end
    end

    % ---- save ONE file for this frequency ----
    nullOut = struct();
    nullOut.primary = struct('powNull', powNull_primary, 'phaseNull', phaseNull_primary);
    nullOut.secondary = struct('powNull', powNull_secondary, 'phaseNull', phaseNull_secondary);

    nullOut.meta = struct();
    nullOut.meta.baseName     = baseName;
    nullOut.meta.fi           = fi;
    nullOut.meta.freqHz       = frex(fi);
    nullOut.meta.frexSelIdx   = frexSelIdx;
    nullOut.meta.fs_tf        = fs_tf;
    nullOut.meta.fs_out       = fs_out;
    nullOut.meta.ds_out       = ds_out;
    nullOut.meta.offsets_out  = offsets_out;
    nullOut.meta.offsets_tf   = offsets_tf;
    nullOut.meta.inhaleOnsetIdx = inhaleOnsetIdx;
    nullOut.meta.poolPrim     = poolPrim;
    nullOut.meta.poolSec      = poolSec;
    nullOut.meta.nShuf        = nShuf;

    nullFile = fullfile(saveDirNull, sprintf('%s_BurstNull_fi%03d.mat', baseName, fi));
    save(nullFile, 'nullOut', '-v7.3');

    clear powNull_primary phaseNull_primary powNull_secondary phaseNull_secondary nullOut
end

% ======================================================
% Compute ITPC (21 × 100) and phase preference at t=0 (freq-wise)
% ======================================================
Z = exp(1i * double(gammaLockTF.phase_primary));
gammaLockTF.itpc_primary = abs(squeeze(mean(Z, 1, 'omitnan')));
v0 = squeeze(mean(exp(1i * double(gammaLockTF.phase_primary(:,centerIdx,:))), 1, 'omitnan'));
gammaLockTF.phasePref_primary_angle = angle(v0);
gammaLockTF.phasePref_primary_r     = abs(v0);

Z = exp(1i * double(gammaLockTF.phase_secondary));
gammaLockTF.itpc_secondary = abs(squeeze(mean(Z, 1, 'omitnan')));
v0 = squeeze(mean(exp(1i * double(gammaLockTF.phase_secondary(:,centerIdx,:))), 1, 'omitnan'));
gammaLockTF.phasePref_secondary_angle = angle(v0);
gammaLockTF.phasePref_secondary_r     = abs(v0);

chanDat.gammaLockTF = gammaLockTF;

disp('gammaLockTF complete (primary/secondary tensors + per-frequency null files + ITPC/phasePref)')
saveDir = fullfile(stem,'CHANDAT_processed');
save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');

%% =========================================================
% STEP D: Lead/Lag analysis (template matching with Pearson correlation)
%
% Empirical:
%   For each breath, build a gamma template from that breath's gamma-peak
%   frequency power time series in a ±500 ms window (100 Hz).
%   Then for each TARGET frequency, slide the template across the low-freq
%   power time series and find the lag (sec) of maximum Pearson correlation.
%   Store lag (sec) + correlation strength into chanDat.leadLag.
%
% Null:
%   For chanDat.use==1 breaths only, compute 1000 circular-shift null lags
%   by circularly rotating ONLY the searchable portion of the low-freq signal.
%   Save null lag cube: [nShuf x nBreaths x nTargetFreq] into a file in:
%       fullfile(stem,'CHANDAT_processed','burstNullFiles')
%
% Rules for sliding window:
%   - Template length L = 2*halfL + 1 where halfL = round(0.5*fs_tf)
%   - Allowed start indices s = 1 ... sMax, where:
%         sMax = min( nT - L + 1, breathEndIdx_tf - halfL )
%     (i.e., stop when template center would align to end of breath OR
%      template end would align to end of trial, whichever comes first)
%   - Lag sign convention:
%         lag_sec = (bestCenterIdx - t0_idx) / fs_tf
%     Negative => low-frequency pattern leads gamma peak.
%
% Notes:
%   - No detrending / ramp removal (as requested)
%   - Uses raw power (as requested)
%   - Uses Pearson correlation via local window normalization
% =========================================================

disp('working on leadLag (Step D)')

% --- inputs ---
if ~isfield(chanDat,'gammaBurst') || ~isfield(chanDat.gammaBurst,'t0_idx')
    error('chanDat.gammaBurst.t0_idx not found.');
end
if ~isfield(chanDat,'gammaBurst') || ~isfield(chanDat.gammaBurst,'fi')
    error('chanDat.gammaBurst.fi not found. (Run gammaBurst block that stores fi.)');
end

t0_idx   = double(chanDat.gammaBurst.t0_idx(:));  % TF (100 Hz) index, breaths x 1
gammaFi  = double(chanDat.gammaBurst.fi(:));      % index into frex, breaths x 1
useVec   = double(chanDat.use(:)) == 1;

nBreaths = numel(t0_idx);

% paths
[~, baseName] = fileparts(chanFiles(filei).name);
saveDirTf   = fullfile(stem,'CHANDAT_processed','tf_files');
saveDirNull = fullfile(stem,'CHANDAT_processed','burstNullFiles');
if ~exist(saveDirNull,'dir'), mkdir(saveDirNull); end

% --- TF sampling rate and timebase (assume consistent across tf files) ---
% Use any valid gammaFi breath to load one file and get fs_tf and nT
tr0 = find(isfinite(gammaFi) & isfinite(t0_idx), 1, 'first');
if isempty(tr0)
    warning('No valid breaths for leadLag; skipping.');
    chanDat.leadLag = struct();
    return
end

powFile0 = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, gammaFi(tr0)));
S0 = load(powFile0, 'tfInfo');
fs_tf = double(S0.tfInfo.fs);         % should be 100
nT    = size(S0.tfInfo.pow,1);        % e.g., 1200
tfTim = double(S0.tfInfo.tim(:));     % optional, for metadata

% --- template length: ±500 ms around gamma peak at 100 Hz ---
halfWin_sec = 0.500;
halfL = round(halfWin_sec * fs_tf);   % 50 if fs_tf=100
L = 2*halfL + 1;                      % 101

% inhale onset and breath end indices (TF samples)
inhaleOnsetIdx = double(chanDat.gammaBurst.inhaleOnsetIdx);
breathEndIdx = round(double(chanDat.behDat.length(:))*fs_tf + inhaleOnsetIdx);
breathEndIdx = min(max(1, breathEndIdx), nT);

% --- target frequencies (Hz), map -> frex indices ---
targetHz = [.2 .6 1.0 1.8 2.6 3.4 4.2 5 6 7 8 9 10 11 12 13 15 17 ...
    19 21 24 27 30 35 40 45 50 55 65 75 85 95 105 115 130 145 160 175 190];
nF = numel(targetHz);

targetFi = nan(1,nF);
targetFrex = nan(1,nF);
for k = 1:nF
    [~, ii] = min(abs(frex - targetHz(k)));
    targetFi(k) = ii;
    targetFrex(k) = frex(ii);
end

% --- output arrays (empirical) ---
lag_sec  = nan(nBreaths, nF);
corr_max = nan(nBreaths, nF);

% --- null setup ---
nShuf = 1000;
lagNull = nan(nShuf, nBreaths, nF, 'single');  % only fill for use==1 breaths

% For reproducibility
rngState = rng;

% =========================================================
% Build gamma templates per breath (raw gamma power snippet)
% Store centered template gc and template std sg for fast correlation
% =========================================================
gammaTemplate_gc = nan(nBreaths, L);   % centered template (g - mean(g))
gammaTemplate_sg = nan(nBreaths, 1);   % std(g) with N-1 normalization

validTemplate = isfinite(t0_idx) & isfinite(gammaFi) & (t0_idx-halfL>=1) & (t0_idx+halfL<=nT);

uFi = unique(gammaFi(validTemplate));
uFi = uFi(isfinite(uFi));

for uf = uFi(:)'
    powFileG = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, uf));
    SG = load(powFileG, 'tfInfo');
    powG = double(SG.tfInfo.pow);  % [nT x nBreaths]

    breathsHere = find(validTemplate & gammaFi==uf);
    for tr = breathsHere(:)'
        g = powG( (t0_idx(tr)-halfL):(t0_idx(tr)+halfL), tr );
        if any(~isfinite(g)), continue; end
        mg = mean(g);
        gc = g - mg;
        sg = std(g,0);  % unbiased (N-1)
        if sg<=0 || ~isfinite(sg), continue; end
        gammaTemplate_gc(tr,:) = gc(:)';
        gammaTemplate_sg(tr)   = sg;
    end
end

haveTemplate = all(isfinite(gammaTemplate_gc),2) & isfinite(gammaTemplate_sg);

% =========================================================
% MAIN LOOP: target frequency (load each low-freq file once)
% =========================================================
onesL = ones(L,1);

for k = 1:nF
    fi = targetFi(k);
    fprintf('leadLag: %d/%d (fi=%d, f=%.3f Hz)\n', k, nF, fi, frex(fi));

    % load low-frequency power file
    powFile = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, fi));
    S = load(powFile, 'tfInfo');
    powLow = double(S.tfInfo.pow);   % [nT x nBreaths]
    if size(powLow,1) ~= nT
        % if length mismatch, clamp to min
        nTcur = min(nT, size(powLow,1));
        powLow = powLow(1:nTcur,:);
    end

    for tr = 1:nBreaths
        if ~haveTemplate(tr), continue; end

        t0 = t0_idx(tr);
        if ~isfinite(t0), continue; end

        % Sliding rule: s = 1 .. sMax
        sMax = min( (nT - L + 1), (breathEndIdx(tr) - halfL) );
        sMax = floor(sMax);
        if sMax < 1, continue; end

        % searchable portion length M so that linear starts 1..sMax use indices 1..M
        M = sMax + L - 1;
        if M > nT, M = nT; end

        ySeg = powLow(1:M, tr);
        if any(~isfinite(ySeg)), continue; end

        % Build circular-extended signal to compute Pearson corr for ALL circular starts 1..M
        yExt = [ySeg; ySeg(1:(L-1))];  % length M+L-1

        gc = gammaTemplate_gc(tr,:)';   % L x 1
        sg = gammaTemplate_sg(tr);      % scalar

        % Numerator for start s: sum(gc .* y_window)
        num = conv(yExt, flipud(gc), 'valid'); % length M

        % Local mean/std for each window (circular)
        sumY  = conv(yExt, onesL, 'valid');    % length M
        sumY2 = conv(yExt.^2, onesL, 'valid'); % length M

        % Unbiased local variance (N-1) form
        varY = (sumY2 - (sumY.^2)/L) / (L-1);
        varY(varY < 0) = 0;
        sdY = sqrt(varY);

        den = (L-1) * sg .* sdY;
        corr_circ = num ./ den;
        corr_circ(~isfinite(corr_circ)) = NaN;

        % Empirical: only starts 1..sMax (no wrap)
        cEmp = corr_circ(1:sMax);
        if all(isnan(cEmp)), continue; end

        [cBest, pos] = max(cEmp, [], 'omitnan'); % pos is best start index (1..sMax)
        if isempty(pos) || ~isfinite(cBest), continue; end

        bestCenter = pos + halfL; % trial index (since starts are in trial coords)
        lag_sec(tr,k)  = (bestCenter - t0) / fs_tf;
        corr_max(tr,k) = cBest;

        % Nulls: only for use==1 breaths
        if ~useVec(tr), continue; end

        % Precompute best-start for EVERY possible circular shift kShift = 0..M-1
        % For shift kShift, the correlation values for linear starts 1..sMax in the shifted signal are:
        %   corrShift(s) = corr_circ( mod(s - kShift - 1, M) + 1 )
        % We find argmax over s=1..sMax for each kShift.

        kVec = (0:(M-1))'; % M x 1
        idxMat = mod((1:sMax) - kVec - 1, M) + 1;     % [M x sMax]
        vals = corr_circ(idxMat);                      % [M x sMax]
        vals(isnan(vals)) = -Inf;

        [maxVal, posMat] = max(vals, [], 2);          % [M x 1]
        posMat(maxVal==-Inf) = NaN;                   % shifts with no valid corr

        % draw 1000 random circular shifts (0..M-1)
        kRand = randi([0 M-1], nShuf, 1);             % with replacement
        bestStartNull = posMat(kRand+1);              % [nShuf x 1] start indices (1..sMax)

        % Convert to lag (sec): center = start + halfL
        % (If bestStartNull is NaN, lag stays NaN)
        lagNull(:,tr,k) = single((bestStartNull + halfL - t0) / fs_tf);
    end
end

% =========================================================
% Store empirical outputs in chanDat.leadLag
% =========================================================
leadLag = struct();
leadLag.fs_tf = fs_tf;
leadLag.halfWin_sec = halfWin_sec;
leadLag.L = L;
leadLag.halfL = halfL;

leadLag.targetHz = targetHz;
leadLag.targetFi = targetFi;
leadLag.targetFrex = targetFrex;

leadLag.t0_idx = t0_idx;
leadLag.gammaFi = gammaFi;

leadLag.lag_sec = lag_sec;
leadLag.corr_max = corr_max;

chanDat.leadLag = leadLag;

% =========================================================
% Save null lags to disk (separate file)
% =========================================================
nullOut = struct();
nullOut.lagNull = lagNull;   % [nShuf x nBreaths x nF] single
nullOut.meta = struct();
nullOut.meta.baseName = baseName;
nullOut.meta.fs_tf = fs_tf;
nullOut.meta.halfWin_sec = halfWin_sec;
nullOut.meta.L = L;
nullOut.meta.halfL = halfL;
nullOut.meta.targetHz = targetHz;
nullOut.meta.targetFi = targetFi;
nullOut.meta.targetFrex = targetFrex;
nullOut.meta.useBreaths = find(useVec);
nullOut.meta.nShuf = nShuf;
nullOut.meta.rngState = rngState;

nullFile = fullfile(saveDirNull, sprintf('%s_LeadLagNull_circShift_n%d.mat', baseName, nShuf));
save(nullFile, 'nullOut', '-v7.3');

disp('leadLag complete: stored chanDat.leadLag and saved null lag cube')
saveDir = fullfile(stem,'CHANDAT_processed');
save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');

%%

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

