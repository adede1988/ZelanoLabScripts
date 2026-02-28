function [] = singleChanEEGPipeline(chanFiles, filei, subFiles, codePre)

%% set frequency parameters

    
frex = logspace(log10(.1),log10(55),150);
% bandWidth = logspace(log10(.1), log10(20), 300);
numfrex = length(frex); 
stds = logspace(log10(3),log10(6),numfrex);

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
    try
        QC = load([chanDat.QCFileDir '/' chanDat.QCFileName]).cleanDat;
    catch
        QCnameBits = strsplit(chanDat.QCFileName, 'cleaningVecs.mat'); 
        chanDat.QCFileName = [QCnameBits{1} char(chanDat.task) ...
                            '_cleaningVecs.mat'];
        QC = load([chanDat.QCFileDir '/' chanDat.QCFileName]).cleanDat;

    end
catch
    chanDat = load([chanFiles(filei).folder '/' chanFiles(filei).name]).chanDat; % go raw if it's not working!
    try
        QC = load([chanDat.QCFileDir '/' chanDat.QCFileName]).cleanDat;
    catch
        QCnameBits = strsplit(chanDat.QCFileName, 'cleaningVecs.mat'); 
        chanDat.QCFileName = [QCnameBits{1} char(chanDat.task) ...
                            '_cleaningVecs.mat'];
        QC = load([chanDat.QCFileDir '/' chanDat.QCFileName]).cleanDat;

    end
end

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])


%% hardCode timing relative to sniff

%hard code that breath window is [-2:10) 
if ~isfield(chanDat, 'trial')
    disp('making trials')
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
            if isempty(QC.blinkCleanVec)
                trialBlink(tr, ok) = zeros(sum(ok),1); 
            else
                trialBlink(tr, ok) = QC.blinkCleanVec(idx(ok));
            end
            if isempty(QC.badTS)
                trialBad(tr, ok) = zeros(sum(ok),1); 
            else
                trialBad(tr, ok)   = QC.badTS(idx(ok));
            end
    
            if strcmp(chanDat.task, 'breathingTask')
                trialRR(tr, ok)   = chanDat.RRint(idx(ok));
                if isempty(chanDat.targTrace)
                    trialTarg(tr, ok) = zeros(sum(ok),1); 
                else
                    trialTarg(tr, ok) = chanDat.targTrace(idx(ok));
                end
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
   disp('breath time indicies')
    
    [idx50, lm] = breathPiecewiseTemplateIdx(chanDat);         % nBreaths x 50
    chanDat.targIDX = idx50; 
    
    try
        chanDat.behDat.length
    catch
        lengthVals = (lm.winEnd - lm.onsetIdx) ./ chanDat.fs; 
        chanDat.behDat.length = lengthVals; 
    end
end

%% QC use/notuse breath by breath

if ~isfield(chanDat, 'use')
    disp('doing QC')
    use = ones(size(chanDat.trial.data,1),1); 
    reasonEliminate = zeros(size(chanDat.trial.data,1), 4); 
    %columns are reasons for eliminating: 
    %col 1: bad breath
    %col 2: couldn't template match
    %col 3: bad EEG
    %col 4: blink
    for ii = 1:length(use)
        try
        if chanDat.behDat.goodBreath(ii) == 0
            use(ii) = 0;
            reasonEliminate(ii,1) = 1; 
        end
        catch
        end
        if sum(isnan(chanDat.targIDX(ii,:)))>10
            use(ii) = 0; 
            reasonEliminate(ii,2) = 1; 
        end
        startIdx = 500; 
        endIdx = max(chanDat.targIDX(ii,:), [],'omitnan'); 
        L = endIdx - startIdx; 
        if ~isnan(L)
        if sum(QC.trial.badTS(ii, startIdx:endIdx)) / L > .50
            use(ii) = 0; 
            reasonEliminate(ii,3) = 1; 
        end
     
    
        if sum(QC.trial.blinkCleanVec(ii, startIdx:endIdx)) / L > .25
            use(ii) = 0; 
            reasonEliminate(ii,4) = 1; 
        end
        end
    
        
    
    end
    
    chanDat.use = use; 
    chanDat.reasonEliminate = reasonEliminate; 
end

%% long time series power/phase extraction
%this code needs to do everything for one frequency at a time to avoid
%problems with memory overload 
if ~isfield(chanDat, 'tf')
disp('working on tf extraction')
    %store into 50 timepoint samples per breath
    outPowZ = nan([size(chanDat.targIDX) numfrex]);
    outPow = nan([size(chanDat.targIDX) numfrex]);
    outPhase = nan([size(chanDat.targIDX) numfrex]);
    outBreathSeg = nan([size(chanDat.targIDX)]);
    outRawSeg = nan(size(chanDat.targIDX));   % breaths x 50 (or however many idx per breath)
    %breaths X phases X frequencies 
    spectra = nan([size(chanDat.targIDX,1) 5 numfrex]);
    for fi = 1:length(frex)
        
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
        rawTrial = chanDat.trial.data;
        rawTrial = movmean(rawTrial,100,2);
        for tt = 1:size(outPowZ, 1)
            idx = chanDat.targIDX(tt,:); 
            idx(isnan(idx)) = []; 
            outPow(tt, 1:length(idx),fi) = curTrialPow(idx,tt); 
            outPowZ(tt,1:length(idx),fi) = powz(idx,tt);
            outPhase(tt, 1:length(idx),fi) = curTrialPhase(idx,tt); 
            outBreathSeg(tt,1:length(idx)) = chanDat.trial.rsp(tt,idx); 
            outRawSeg(tt, 1:numel(idx)) = rawTrial(tt, idx);
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

    if isfield(chanDat, 'sessNum')
        curSess = chanDat.sessNum; 
    else
        curSess = 1; 
        chanDat.sessNum = 1; 
    end
    % Save metadata (note fs/tim are for the SAVED files)
    tfMeta = struct( ...
        'subID',chanDat.subID,'task',chanDat.task,'tim',timSave,'fs',fs_tf,'use',chanDat.use, ...
        'targIDX',chanDat.targIDX,'chi',chanDat.chi,'chanType',chanDat.chanType,'behDat',chanDat.behDat, ...
        'sessID',chanDat.sessID,'sessNum',curSess, ...
        'OGdataDir',chanDat.OGdataDir,'fi',fi,'frex',frex(fi), ...
        'fs_full',fs_full,'dsFactor',ds);

    tfMeta.labels = chanDat.labels;

    % Save pow
    % tfInfo = tfMeta;
    % tfInfo.pow = powSave;
    % save(fullfile(saveDir, sprintf('%s_pow_fi%03d.mat', baseName, fi)), 'tfInfo', '-v7.3');
    % 
    % % Save phase
    % tfInfo = tfMeta;
    % tfInfo.phase = phaseSave;
    % save(fullfile(saveDir, sprintf('%s_phase_fi%03d.mat', baseName, fi)), 'tfInfo', '-v7.3');

    
    
    
    
    end

    rawTrial = chanDat.trial.data;

    %store for saving out
    chanDat.tf.phase = outPhase; 
    chanDat.tf.pow = outPow; 
    chanDat.tf.powZ = outPowZ; 
    chanDat.tf.breathSeg = outBreathSeg; 
    chanDat.tf.frex = frex; 
    chanDat.tf.spectra = spectra; 
    chanDat.tf.rawSeg = outRawSeg;
    
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

    

%% fooof: 

% 
% if ~isfield(chanDat, 'fooof')
%     disp('working on fooof')
%     %--- Aperiodic subtraction per breath × epoch WITH KNEE ---
%     % Requires: spectra (nBreaths x 5 x numfrex), frex (numfrex x 1 or 1 x numfrex)
%     spectra = chanDat.tf.spectra; 
%     frex = frex(:);
% 
%     % ----- USER SETTINGS -----
%     fitRangeHz = [.5 150];
%     excludeLineHz = [60 120];        % [] to disable
%     excludeHalfWidthHz = 2;
%     minFitPoints = 12;               
% 
%     fitMask = frex >= fitRangeHz(1) & frex <= fitRangeHz(2);
%     if ~isempty(excludeLineHz)
%         for f0 = excludeLineHz(:)'
%             fitMask = fitMask & ~(frex > (f0 - excludeHalfWidthHz) & frex < (f0 + excludeHalfWidthHz));
%         end
%     end
% 
%     % Outputs
%     spectra_flat_log10 = nan(size(spectra));              % log10(power) residual
%     aperiodic_log10    = nan(size(spectra));              % fitted aperiodic in log10(power)
%     aperiodic_params   = nan([size(spectra,1) size(spectra,2) 3]); % [offset exponent knee]
% 
%     % Knee aperiodic model in log10 space:
%     % y = offset - log10(knee + f^exponent)
%     % We parameterize exponent = exp(logExp) and knee = 10^(klog) to enforce positivity.
%     aper_model = @(th,f) th(1) - log10( (10.^th(3)) + (f.^exp(th(2))) );
% 
%     useLSQ = exist('lsqcurvefit','file') == 2;
% 
%     if useLSQ
%         opts = optimoptions('lsqcurvefit','Display','off');
%         % Reasonable soft bounds in transformed space:
%         lb = [-Inf, log(1e-3), -12];   % exponent >= 1e-3, knee >= 1e-12
%         ub = [ Inf, log(50),   12];   % exponent <= 50,   knee <= 1e12
%     end
% 
%     for tt = 1:size(spectra,1) %trials 
%         for ee = 1:size(spectra,2) % epochs
% 
%             p = squeeze(spectra(tt,ee,:));  % numfrex x 1
%             if all(isnan(p)), continue; end
% 
%             p(p <= 0) = eps;               % guard
%             y = log10(p);
% 
%             m = fitMask & ~isnan(y);
%             if nnz(m) < minFitPoints, continue; end
% 
%             % --- initial guesses from straight-line (no-knee) fit in log-log ---
%             xf = log10(frex(m));
%             cf = polyfit(xf, y(m), 1);         % y ≈ slope*log10(f) + intercept
%             slope0    = cf(1);
%             offset0   = cf(2);
%             exponent0 = max(1e-3, -slope0);    % knee model exponent is positive
% 
%             fmin = min(frex(m));
%             knee0 = max(1e-12, (fmin.^exponent0) / 10);  % small knee to start
% 
%             th0 = [offset0, log(exponent0), log10(knee0)];
% 
%             % --- fit ---
%             if useLSQ
%                 th = lsqcurvefit(@(th,f) aper_model(th,f), th0, frex(m), y(m), lb, ub, opts);
%             else
%                 obj = @(th) sum((aper_model(th, frex(m)) - y(m)).^2);
%                 th  = fminsearch(obj, th0, optimset('Display','off'));
%             end
% 
%             yfit = aper_model(th, frex);
% 
%             offset   = th(1);
%             exponent = exp(th(2));
%             knee     = 10.^th(3);
% 
%             aperiodic_log10(tt,ee,:)    = yfit;
%             spectra_flat_log10(tt,ee,:) = y - yfit;
%             aperiodic_params(tt,ee,:)   = [offset exponent knee];
%         end
%     end
% 
%     % Optional linear-space versions:
%     aperiodic_fit = 10.^aperiodic_log10;     % fitted 1/f(+knee) component in power units
%     spectra_flat  = 10.^spectra_flat_log10;  % equals spectra ./ aperiodic_fit (ratio)
% 
%     f_knee = aperiodic_params(:,:,3).^(1 ./ aperiodic_params(:,:,2));  % Hz
% 
% 
%     chanDat.fooof.aperiodic_fit = aperiodic_fit; 
%     chanDat.fooof.spectra_flat = spectra_flat; 
%     chanDat.fooof.aperiodic_params = aperiodic_params; 
%     chanDat.fooof.f_knee = f_knee; 
% 
%     % --- Gamma peak frequency per breath × epoch (25–60 Hz) ---
%     gammaMask = frex >= 25 & frex <= 60;
% 
%     tmp = spectra_flat(:,:,gammaMask);   % [nBreaths x 5 x nGamma]
%     allNan = all(isnan(tmp),3);          % breaths/epochs missing
% 
%     tmp(isnan(tmp)) = -Inf;              % so max ignores NaNs
%     [~, imax] = max(tmp, [], 3);         % argmax within gamma band
%     gammaFreqs = frex(gammaMask);        % Hz values for gamma bins
% 
%     gamma_peak_freq = gammaFreqs(imax);  % [nBreaths x 5]
%     gamma_peak_freq(allNan) = NaN;
% 
%     chanDat.fooof.gamma_peak_freq = gamma_peak_freq;
% 
%     % --- Low-frequency peak per breath × epoch (3–14 Hz; theta/alpha) ---
%     lowMask = frex >= .5 & frex <= 14;
% 
%     tmp = spectra_flat(:,:,lowMask);     % [nBreaths x 5 x nLow]
%     allNan = all(isnan(tmp),3);          % breaths/epochs missing
% 
%     tmp(isnan(tmp)) = -Inf;              % so max ignores NaNs
%     [~, imax] = max(tmp, [], 3);         % argmax within low band
%     lowFreqs = frex(lowMask);            % Hz values for low bins
% 
%     low_peak_freq = lowFreqs(imax);      % [nBreaths x 5]
%     low_peak_freq(allNan) = NaN;
% 
%     chanDat.fooof.low_peak_freq = low_peak_freq;
% 
% 
%     figure;
%     plot(squeeze(mean(spectra_flat, 1,'omitnan'))')  % 5 x numfrex -> plotted as 5 lines
%     xticks(20:20:300)
%     xticklabels(round(frex(20:20:300)))
% 
%     legend({'inhale rise','inhale fall','exhale rise','exhale fall','pause'}, ...
%            'Location','best','Box','off');
% 
% 
%     % --- Peak frequency (Hz) per breath, per epoch, within 30–50 Hz ---
%     % Inputs assumed in workspace:
%     %   spectra_flat : [nBreaths x 5 x numfrex]
%     %   frex         : [numfrex x 1] or [1 x numfrex]
%     saveDir = fullfile(stem,'CHANDAT_processed'); 
%     save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');
% end

%% fooof (FOOOF_BASIC)

if ~isfield(chanDat, 'fooof')
    disp('working on fooof (fooof_basic)')
    % ---------------- USER SETTINGS ----------------
    fitRangeHz = [2 150];
    excludeLineHz = [60 120];        % [] to disable
    excludeHalfWidthHz = 2;

    % Inputs
    spectra = chanDat.tf.spectra;            % [nBreaths x 5 x numfrex] linear power
    frex    = frex(:);                       % [numfrex x 1]
    nBreaths = size(spectra,1);
    nEpochs  = size(spectra,2);
    nF       = numel(frex);

    % Use vector (fit only QC-passed breaths)
    if isfield(chanDat,'use') && ~isempty(chanDat.use)
        useVec = chanDat.use(:)==1;
    else
        useVec = true(nBreaths,1);
    end

    

    % Avoid interpreting extrapolation: set flattened outputs outside fit range to NaN
    nanOutsideFitRange = false;

    % fooof_basic options (tune as needed)
    fooofArgs = { ...
        'f_range', fitRangeHz, ...
        'aperiodic_mode', 'knee', ...
        'max_peaks', 3, ...
        'peak_thresh', 2.5, ...
        'peak_width_limits', [1 10], ...
        'smooth_w', 5, ...
        'max_clip_iters', 3, ...
        'verbose', false ...
    };

    % ---------------- Masks ----------------
    inFitRange = (frex >= fitRangeHz(1)) & (frex <= fitRangeHz(2)) & (frex > 0);

    lineMask = false(nF,1);
    if ~isempty(excludeLineHz)
        for f0 = excludeLineHz(:)'
            lineMask = lineMask | (frex > (f0 - excludeHalfWidthHz) & frex < (f0 + excludeHalfWidthHz));
        end
    end

    fitMask = inFitRange & ~lineMask;

    % ---------------- Preallocate outputs ----------------
    aperiodic_log10      = nan(nBreaths, nEpochs, nF, 'single');
    spectra_flat_log10   = nan(nBreaths, nEpochs, nF, 'single');
    aperiodic_fit        = nan(nBreaths, nEpochs, nF, 'single');
    spectra_flat         = nan(nBreaths, nEpochs, nF, 'single');
    gamma_peaks          = nan(nBreaths, nEpochs, 'single');
    beta_peaks           = nan(nBreaths, nEpochs, 'single');
    alpha_peaks          = nan(nBreaths, nEpochs, 'single');
    theta_peaks          = nan(nBreaths, nEpochs, 'single');
    delta_peaks          = nan(nBreaths, nEpochs, 'single');

    % [offset exponent kneeParam] where kneeParam is k in linear units (FOOOF knee)
    aperiodic_params     = nan(nBreaths, nEpochs, 3, 'single');
    f_knee               = nan(nBreaths, nEpochs, 'single');

    % Optional: store peak tables per breath×epoch
    peaks_cell           = cell(nBreaths, nEpochs);

    % ---------------- Main loop ----------------
    for tt = 1:nBreaths
        % if ~useVec(tt), continue; end

        for ee = 1:nEpochs
            p = squeeze(spectra(tt,ee,:));          % [nF x 1]
            if all(isnan(p)), continue; end

            % Guard power
            p = double(p);
            p(p <= 0) = eps;

            % Build PSD for fitting: exclude line bins by setting to NaN
            pFit = p;
            pFit(lineMask) = NaN;

            % Require enough points in fit mask
            if nnz(fitMask & isfinite(pFit)) < 16
                continue
            end

            % Run fooof_basic (knee), fallback to fixed if knee fails
            try
                out = fooof_basic(frex, pFit, fooofArgs{:});
            catch
                try
                    out = fooof_basic(frex, pFit, ...
                        'f_range', fitRangeHz, ...
                        'aperiodic_mode','fixed', ...
                        fooofArgs{3:end});
                catch
                    continue
                end
            end

            % Aperiodic model in log10 power (full length)
            ap_log10 = double(out.full.ap_log10(:));     % [nF x 1]

            % Flattened residual in log10: log10(p) - ap
            y = log10(p);
            flat_log10 = y - ap_log10;

            % Optional: blank outside the fit range (prevents end-range curvature dominating plots)
            if nanOutsideFitRange
                flat_log10(~inFitRange) = NaN;
                ap_log10(~inFitRange)   = NaN;
            end


          

            % Also blank line bins (so you don't interpret line-noise “structure”)
            flat_log10(lineMask) = NaN;
            ap_log10(lineMask)   = NaN;

            % Store
            aperiodic_log10(tt,ee,:)    = single(ap_log10);
            spectra_flat_log10(tt,ee,:) = single(flat_log10);

            ap_fit = 10.^ap_log10;                  % linear power
            flat_ratio = 10.^flat_log10;            % ratio p./ap (since flat_log10 = log10(p)-log10(ap))



            % ---- grab peak values from fooof output (store peak center freqs) ----
            % out.peaks: [nPeaks x 3] = [center_Hz, amp_log10, fwhm_Hz]
            % flat_ratio: [numfrex x 1] ratio spectrum aligned to frex (e.g., 10.^flattened_log10)
            
            if ~isempty(out.peaks) && size(out.peaks,1) > 0
            
                pkHz = out.peaks(:,1);   % centers in Hz
            
                % ----- GAMMA (25–60) -----
                inBand = (pkHz >= 25) & (pkHz < 60);
                if any(inBand)
                    idx = find(inBand);
                    if numel(idx) == 1
                        gamma_peaks(tt,ee) = single(pkHz(idx));
                    else
                        score = nan(numel(idx),1);
                        for ii = 1:numel(idx)
                            [~, fiPk] = min(abs(frex - pkHz(idx(ii))));
                            score(ii) = flat_ratio(fiPk);
                        end
                        [~, imax] = max(score);
                        gamma_peaks(tt,ee) = single(pkHz(idx(imax)));
                    end
                end
            
                % ----- BETA (14–25) -----
                inBand = (pkHz >= 14) & (pkHz < 25);
                if any(inBand)
                    idx = find(inBand);
                    if numel(idx) == 1
                        beta_peaks(tt,ee) = single(pkHz(idx));
                    else
                        score = nan(numel(idx),1);
                        for ii = 1:numel(idx)
                            [~, fiPk] = min(abs(frex - pkHz(idx(ii))));
                            score(ii) = flat_ratio(fiPk);
                        end
                        [~, imax] = max(score);
                        beta_peaks(tt,ee) = single(pkHz(idx(imax)));
                    end
                end
            
                % ----- ALPHA (8.5–14) -----
                inBand = (pkHz >= 8.5) & (pkHz < 14);
                if any(inBand)
                    idx = find(inBand);
                    if numel(idx) == 1
                        alpha_peaks(tt,ee) = single(pkHz(idx));
                    else
                        score = nan(numel(idx),1);
                        for ii = 1:numel(idx)
                            [~, fiPk] = min(abs(frex - pkHz(idx(ii))));
                            score(ii) = flat_ratio(fiPk);
                        end
                        [~, imax] = max(score);
                        alpha_peaks(tt,ee) = single(pkHz(idx(imax)));
                    end
                end
            
                % ----- THETA (4–8.5) -----
                inBand = (pkHz >= 4) & (pkHz < 8.5);
                if any(inBand)
                    idx = find(inBand);
                    if numel(idx) == 1
                        theta_peaks(tt,ee) = single(pkHz(idx));
                    else
                        score = nan(numel(idx),1);
                        for ii = 1:numel(idx)
                            [~, fiPk] = min(abs(frex - pkHz(idx(ii))));
                            score(ii) = flat_ratio(fiPk);
                        end
                        [~, imax] = max(score);
                        theta_peaks(tt,ee) = single(pkHz(idx(imax)));
                    end
                end
            
                % ----- DELTA (0.5–4) -----
                inBand = (pkHz >= 0.5) & (pkHz < 4);
                if any(inBand)
                    idx = find(inBand);
                    if numel(idx) == 1
                        delta_peaks(tt,ee) = single(pkHz(idx));
                    else
                        score = nan(numel(idx),1);
                        for ii = 1:numel(idx)
                            [~, fiPk] = min(abs(frex - pkHz(idx(ii))));
                            score(ii) = flat_ratio(fiPk);
                        end
                        [~, imax] = max(score);
                        delta_peaks(tt,ee) = single(pkHz(idx(imax)));
                    end
                end
            
            end


            aperiodic_fit(tt,ee,:) = single(ap_fit);
            spectra_flat(tt,ee,:)  = single(flat_ratio);

            % Params
            b   = out.ap.offset;
            chi = out.ap.exponent;
            k   = out.ap.knee;                      % 0 if fixed
            aperiodic_params(tt,ee,:) = single([b chi k]);
            f_knee(tt,ee) = single(out.ap.knee_freq);

            % Peaks (optional)
            peaks_cell{tt,ee} = out.peaks;          % [nPeaks x 3] = [centerHz amp_log10 fwhmHz]
        end
    end

    % ---------------- Save into chanDat ----------------
    chanDat.fooof = struct();
    chanDat.fooof.aperiodic_fit      = aperiodic_fit;
    chanDat.fooof.aperiodic_log10    = aperiodic_log10;

    chanDat.fooof.spectra_flat       = spectra_flat;
    chanDat.fooof.spectra_flat_log10 = spectra_flat_log10;

    chanDat.fooof.aperiodic_params   = aperiodic_params;
    chanDat.fooof.f_knee             = f_knee;

    % bookkeeping / reproducibility
    chanDat.fooof.fitRangeHz         = fitRangeHz;
    chanDat.fooof.excludeLineHz      = excludeLineHz;
    chanDat.fooof.excludeHalfWidthHz = excludeHalfWidthHz;
    chanDat.fooof.nanOutsideFitRange = nanOutsideFitRange;
    chanDat.fooof.fooofArgs          = fooofArgs;
    chanDat.fooof.delta_peaks        = delta_peaks; 
    chanDat.fooof.theta_peaks        = theta_peaks; 
    chanDat.fooof.alpha_peaks        = alpha_peaks; 
    chanDat.fooof.beta_peaks         = beta_peaks; 
    chanDat.fooof.gamma_peaks        = gamma_peaks; 

    % optional peaks
    chanDat.fooof.peaks              = peaks_cell;

    %macro specific stuff: 
    % % ---------------- Peak frequencies (same logic as you had) ----------------
    % % Gamma peak (25–60 Hz) using spectra_flat ratio (QC breaths were fit; others are NaN)
    % gammaMask = frex >= 25 & frex <= 60;
    % tmp = chanDat.fooof.spectra_flat(:,:,gammaMask);   % [nBreaths x 5 x nGamma]
    % allNan = all(isnan(tmp),3);
    % tmp(isnan(tmp)) = -Inf;
    % [~, imax] = max(tmp, [], 3);
    % gammaFreqs = frex(gammaMask);
    % gamma_peak_freq = gammaFreqs(imax);
    % gamma_peak_freq(allNan) = NaN;
    % chanDat.fooof.gamma_peak_freq = gamma_peak_freq;
    % 
    % % Low-frequency peak (0.5–14 Hz) using spectra_flat ratio
    % lowMask = frex >= .5 & frex <= 14;
    % tmp = chanDat.fooof.spectra_flat(:,:,lowMask);
    % allNan = all(isnan(tmp),3);
    % tmp(isnan(tmp)) = -Inf;
    % [~, imax] = max(tmp, [], 3);
    % lowFreqs = frex(lowMask);
    % low_peak_freq = lowFreqs(imax);
    % low_peak_freq(allNan) = NaN;
    % chanDat.fooof.low_peak_freq = low_peak_freq;

    % Save
    saveDir = fullfile(stem,'CHANDAT_processed');
    save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');
end

%% bring in the best macro channel for this participant: 

disp('bringing in the macro channel')

epochNames = ["inhale rise","inhale fall","exhale rise","exhale fall","pause"];
gammaBandHz     = [25 60];
baselineBandHz  = [10 200];
excludeHzAroundPeak = 5;

test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, '_macro'));
groupFiles = subFiles(test); 

chanSumm = struct([]);  % <-- NOT repmat(struct(),...)
cci = 1; 
for ci = 1:numel(groupFiles)
    fpath = fullfile(stem, 'CHANDAT_processed/', groupFiles(ci).name); 
    if exist(fpath, 'file')
    S = load(fpath, "chanDat");
    macChan = S.chanDat;

    s = summarize_channel_gamma(macChan, epochNames, ...
        gammaBandHz, baselineBandHz, excludeHzAroundPeak);

    s.filePath = string(fpath);
    s.fileName = string(groupFiles(ci).name);

    if cci == 1
        chanSumm = repmat(s, numel(groupFiles), 1);  % template has all fields
    end

    chanSumm(cci) = s;
    cci = cci + 1; 
    end
end

% Pick winner (robust z-scored composite within this group)
[winnerIdx, scoreTable, chanSumm] = pick_winner_channel(chanSumm);

winner = chanSumm(winnerIdx);
fprintf("Winner: %s  (score=%.3f)\n", winner.chanLabelSafe, winner.score);
sessNum = winner.sessNum;
taskName = chanDat.task;
prefix = sprintf("%s_sess%d_%s_%s", winner.subIDSafe, sessNum, winner.chanLabelSafe, taskName);


Sw = load(winner.filePath, "chanDat");
macChan = Sw.chanDat;




%% PAC between macro and EEG: 

rawDat = load([chanFiles(filei).folder '/' chanFiles(filei).name]).chanDat;

macRaw = load([chanFiles(filei).folder '/' macChan.sessID '_' ...
                    macChan.chanType '_' macChan.task '_' ...
                    num2str(macChan.chi) '.mat']).chanDat; 
deliberate error! 
%think here about adding in PAC relative to gamma Burst, but juice might
%not be worth the squeeze
if ~isfield(chanDat, 'pac')
    disp('doing pac')
    keyBreathIDX = chanDat.targIDX; 
    onsets = chanDat.behDat.finalOnset; 
    keyBreathIDX = keyBreathIDX + onsets; 
    keyBreathIDX(:,51) = macChan.gammaBurst.t0_idx_full + onsets;
    gamMed = median(macChan.fooof.gamma_peaks, 'all', 'omitnan');
    
    % X: [breaths x time]
    fs = chanDat.fs;               % Hz
    
    halfBW = 5;                    % +/- 5 Hz
    bpHz   = double([gamMed-halfBW, gamMed+halfBW]);
    
    PACfrex = logspace(log10(2), log10(25), 20); 
    
    
     [pacOut, meta] = pac_breathTemplate_timeResolvedPAC(rawDat,...
         macRaw, keyBreathIDX, gamMed, fs, bpHz, PACfrex);
    chanDat.pac = meta; 
    chanDat.pac.pac = pacOut; 
    
   
end

%% while the raw data are here, get better HRV measures! 

chanDat = addBreathHRV_fromRRint(chanDat, rawDat, fs);

%% Low frequency breathing lock in .11 to 2 Hz

out = helper_breathISPC(rawDat);
chanDat.breathLock = out; 
saveDir = fullfile(stem,'CHANDAT_processed');
save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');

%% go back for the respiration data and get phase locking with respiration: 

% gamMed = median(macChan.fooof.gamma_peaks, 'all', 'omitnan');
% 
% % X: [breaths x time]
% fs = chanDat.fs;               % Hz
% 
% halfBW = 5;                    % +/- 5 Hz
% bpHz   = double([gamMed-halfBW, gamMed+halfBW]);
% 
% 
% % --- Design + zero-phase filter (Butterworth IIR) ---
% filtOrder = 4; % filtfilt makes it effectively 8th-order magnitude response
% d = designfilt('bandpassiir', ...
%     'FilterOrder', 2*filtOrder, ...            % designfilt order is overall order
%     'HalfPowerFrequency1', bpHz(1), ...
%     'HalfPowerFrequency2', bpHz(2), ...
%     'SampleRate', double(fs), ...
%     'DesignMethod', 'butter');
% XMac = macRaw.data; 
% % filtfilt 
% X_bp = filtfilt(d, XMac);
% X_bp = log(abs(hilbert(X_bp)).^2); 



% % % 
% % % 
% % % % ============================================================
% % % % Breath-by-breath lead/lag correlations (raw time series)
% % % % Output: R is [nBreaths x nOffsets x nCenters]
% % % %   offsets = -250:25:250  (negative => macChan earlier than chanDat)
% % % %   window  = 501 samples, centered at centers = 250:50:5750
% % % % Notes:
% % % %   Any (offset,center) whose 501-sample window would go out of bounds
% % % %   in either signal is left as NaN (still preserves the requested size).
% % % % ============================================================
% % % 
% % % X = chanDat.trial.data;   % [breaths x time]
% % % Y = macChan.trial.data;   % [breaths x time]
% % % 
% % % offsets = -250:25:250;    % 21 offsets
% % % centers = 1:6000;    % 111 center points
% % % winLen  = 501;
% % % 
% % % R = leadLagCorrTrials(X, Y, offsets, centers, winLen);  % [breaths x 21 x 111]
% % % 
% % % % (optional) quick sanity:
% % % size(R)
% % % 
% % % Rnorm = nan([size(chanDat.targIDX) size(R,2)]);
% % % 
% % % for ii = 1:size(Rnorm,1)
% % %     if ~isnan(chanDat.targIDX(ii,1))
% % %         curIDX = chanDat.targIDX(ii,:); 
% % %         curIDX(isnan(curIDX)) = []; 
% % %         curMat = R(ii, :, curIDX); 
% % %         curMat = permute(curMat, [1 3 2]);
% % %         Rnorm(ii,1:length(curIDX), :) = curMat; 
% % %     end
% % % end
% % % 
% % % figure; imagesc([], offsets, ...
% % %                 squeeze(mean(Rnorm(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'audio'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan'))')
% % % 
% % % figure; imagesc([], offsets, ...
% % %                 squeeze(mean(Rnorm(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'focus'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan'))')
% % % 
% % % figure; imagesc([], offsets, ...
% % %                 squeeze(mean(Rnorm(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'shadow'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan'))')
% % % 
% % % 
% % % 
% % % X = chanDat.trial.data;   % [breaths x time]
% % % Y = macChan.trial.data;   % [breaths x time]
% % % 
% % % gamMed = median(macChan.fooof.gamma_peaks, 'all', 'omitnan');
% % % 
% % % % X: [breaths x time]
% % % fs = chanDat.fs;               % Hz
% % % 
% % % halfBW = 5;                    % +/- 5 Hz
% % % bpHz   = double([gamMed-halfBW, gamMed+halfBW]);
% % % 
% % % 
% % % % --- Design + zero-phase filter (Butterworth IIR) ---
% % % filtOrder = 4; % filtfilt makes it effectively 8th-order magnitude response
% % % d = designfilt('bandpassiir', ...
% % %     'FilterOrder', 2*filtOrder, ...            % designfilt order is overall order
% % %     'HalfPowerFrequency1', bpHz(1), ...
% % %     'HalfPowerFrequency2', bpHz(2), ...
% % %     'SampleRate', double(fs), ...
% % %     'DesignMethod', 'butter');
% % % 
% % % % filtfilt filters along columns, so transpose to [time x breaths]
% % % X_bp = filtfilt(d, double(X(1:end-1,:)).');
% % % 
% % % X_bp = (abs(hilbert(X_bp))).';
% % % X_bp(end+1,:) = nan(size(X_bp(end,:))); 
% % % 
% % % offsets = -250:25:250;    % 21 offsets
% % % centers = 1:6000;    % 111 center points
% % % winLen  = 501;
% % % 
% % % R = leadLagCorrTrials(X_bp, Y, offsets, centers, winLen);  % [breaths x 21 x 111]
% % % 
% % % % (optional) quick sanity:
% % % size(R)
% % % 
% % % Rnorm2 = nan([size(chanDat.targIDX) size(R,2)]);
% % % 
% % % for ii = 1:size(Rnorm2,1)
% % %     if ~isnan(chanDat.targIDX(ii,1))
% % %         curIDX = chanDat.targIDX(ii,:); 
% % %         curIDX(isnan(curIDX)) = []; 
% % %         curMat = R(ii, :, curIDX); 
% % %         curMat = permute(curMat, [1 3 2]);
% % %         Rnorm2(ii,1:length(curIDX), :) = curMat; 
% % %     end
% % % end
% % % 
% % % figure; imagesc([], offsets, ...
% % %                 squeeze(mean(Rnorm2(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'audio'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan'))')
% % % 
% % % figure; imagesc([], offsets, ...
% % %                 squeeze(mean(Rnorm2(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'focus'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan'))')
% % % 
% % % figure; imagesc([], offsets, ...
% % %                 squeeze(mean(Rnorm2(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'shadow'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan'))')
% % % 
% % % 
% % % 
% % % 
% % % figure; imagesc(centers, offsets, ...
% % %                 squeeze(mean(R(chanDat.use ==1 & ...
% % %                 cellfun(@(x) strcmp(x, 'audio'), chanDat.behDat.task),...
% % %                 :,:), 1, 'omitnan')))
% % % 
% % % 



















gamERP = nan(length(macChan.gammaEnv.evBreath), 2001); 
w = -1000:1000; 
us = chanDat.fs / macChan.gammaEnv.fs; 
for ii = 1:length(macChan.gammaEnv.evBreath)
    t0 = macChan.gammaEnv.evT0(ii); 
    bi = macChan.gammaEnv.evBreath(ii); 
    t0 = round(t0*us); 
    if t0>1001 && t0<4999
        gamERP(ii,:) = chanDat.trial.data(bi,t0+w);
    end



end

nBreath = length(macChan.gammaBurstSecondary.t0_idx); 

gamERP = nan(nBreath, 2001); 
w = -1000:1000; 
us = chanDat.fs / macChan.gammaBurstSecondary.fs_tf; 
for ii = 1:nBreath
    t0 = macChan.gammaBurstSecondary.t0_idx(ii); 
    t0 = round(t0*us); 
    if t0>1001 && t0<4999
        gamERP(ii,:) = chanDat.trial.data(ii,t0+w);
    end

end

test = mean(gamERP(:,1000:1500), 2, 'omitnan'); 

% -----------------------------
% 1000 random mean-ERP samples
% (<= 1 event per breath/sample)
% -----------------------------

nIter = 1000;
w     = -1000:1000;                          % window in chanDat samples
us    = chanDat.fs / macChan.gammaEnv.fs;    % upsample factor

data = chanDat.trial.data;
[nBreaths, nT] = size(data);

evBreath = macChan.gammaEnv.evBreath(:);     % breath index (row in data)
evT0     = macChan.gammaEnv.evT0(:);         % event time (in macChan.gammaEnv.fs)

% Convert event times to chanDat sample indices
t0s = round(evT0 * us);

%overwrite with values from gammaburst 
t0s = macChan.gammaBurst.t0_idx_full; 
evBreath = [1:length(t0s)]; 
t0s = t0s(:); 
evBreath = evBreath(:); 

% Keep only events whose snippet fits fully in bounds
valid = isfinite(evBreath) & isfinite(t0s) & ...
        evBreath >= 1 & evBreath <= nBreaths & ...
        (t0s + w(1)  >= 1) & (t0s + w(end) <= nT);

evBreath = evBreath(valid);
t0s      = t0s(valid);
useVec   = chanDat.use(valid) == 1; 

nEv = numel(evBreath);
if nEv == 0
    error('No valid events after bounds checking.');
end

% Extract ERP snippet for each valid event
gamERP = nan(nEv, numel(w));
for ii = 1:nEv
    gamERP(ii,:) = data(evBreath(ii), t0s(ii) + w);
end

% Group event indices by breath
eventsByBreath = accumarray(evBreath, (1:nEv)', [nBreaths 1], @(x){x}, {[]});
breathList     = find(~cellfun(@isempty, eventsByBreath));  % breaths with >=1 valid event

% Bootstrap: pick ONE event per breath, then average across breaths
meanERP_boot = nan(nIter, numel(w));
selEventIdx  = nan(nIter, numel(breathList)); % optional: which events were chosen each iteration

for it = 1:nIter
    sel = zeros(numel(breathList), 1);
    for jj = 1:numel(breathList)
        b = breathList(jj);
        idxs = eventsByBreath{b};
        sel(jj) = idxs(randi(numel(idxs)));   % choose 1 event from that breath
    end
    selEventIdx(it,:)  = sel;
    meanERP_boot(it,:) = mean(gamERP(sel,:), 1, 'omitnan');
end

% meanERP_boot is [1000 x 2001], each row is a random mean ERP

% ============================================================
% Condition-specific regular ERPs (no bootstrap)
% Conditions: audio / focus / shadow
% Uses only breaths where chanDat.use == 1
% Mean ± SEM across events (breaths) within each condition
% ============================================================

condNames = ["audio","focus","shadow"];

w  = -1000:1000;                           % window in chanDat samples
data = chanDat.trial.data;                % [breaths x time]
[nBreaths, nT] = size(data);

useBreath  = chanDat.use(:) == 1;         % [breaths x 1]
taskBreath = string(chanDat.behDat.task(:));

% Use gammaBurst indices (as in your original block)
t0s     = macChan.gammaBurst.t0_idx_full(:);   % per-breath event time (chanDat samples)
evBreath = (1:numel(t0s))';                    % breath index

% Valid breaths/events
valid = isfinite(evBreath) & isfinite(t0s) & ...
        evBreath>=1 & evBreath<=nBreaths & ...
        (t0s+w(1) >= 1) & (t0s+w(end) <= nT) & ...
        useBreath(evBreath);

tSec = w / chanDat.fs;
cols = lines(numel(condNames));

figure; hold on

for cc = 1:numel(condNames)
    cond = condNames(cc);

    inCond = valid & (taskBreath(evBreath) == cond);
    evB = evBreath(inCond);
    t0  = t0s(inCond);

    nEv = numel(evB);
    if nEv < 2
        warning('Too few events for condition "%s" (n=%d).', cond, nEv);
        continue
    end

    % Extract ERP snippets (events x time)
    gamERP = nan(nEv, numel(w));
    for ii = 1:nEv
        gamERP(ii,:) = data(evB(ii), t0(ii) + w);
    end

    % Same normalization you had (global over this condition's events)
    % gamERP = gamERP - prctile(gamERP(:), 2);
    % gamERP = gamERP ./ prctile(gamERP(:), 98);

    mu  = mean(gamERP, 1, 'omitnan');
    sem = std(gamERP, 0, 1, 'omitnan') ./ sqrt(nEv);

    % optional smoothing (keep consistent with your prior plot)
    muS  = movmean(mu, 20);
    semS = movmean(sem,20);

    x  = tSec(:);
    lo = (muS - semS).';
    hi = (muS + semS).';

    fill([x; flipud(x)], [lo; flipud(hi)], cols(cc,:), ...
        'FaceAlpha', 0.2, 'EdgeColor','none', 'HandleVisibility','off');
    plot(tSec, muS, 'LineWidth', 2, 'Color', cols(cc,:));
end

xlabel('Time (s)'); ylabel('ERP (uV)');
legend(condNames, 'Location','best'); box off
title('Condition-specific mean ERP (mean ± SEM; use==1)');
xlim([-.5 .5])

% --- Add null (t0-shuffle) lines: 100 shuffles -> 1 grey line per condition ---
nShuf = 100;

figure; hold on
cols = lines(numel(condNames));

for cc = 1:numel(condNames)
    cond = condNames(cc);

    inCond = valid & (taskBreath(evBreath) == cond);
    evB = evBreath(inCond);
    t0  = t0s(inCond);

    nEv = numel(evB);
    if nEv < 2
        warning('Too few events for condition "%s" (n=%d).', cond, nEv);
        continue
    end

    % -------- Observed ERP (mean ± SEM across events) --------
    gamERP = nan(nEv, numel(w));
    for ii = 1:nEv
        gamERP(ii,:) = data(evB(ii), t0(ii) + w);
    end
    % gamERP = gamERP - prctile(gamERP(:), 2);
    % gamERP = gamERP ./ prctile(gamERP(:), 98);

    mu  = mean(gamERP, 1, 'omitnan');
    sem = std(gamERP, 0, 1, 'omitnan') ./ sqrt(nEv);

    muS  = movmean(mu, 20);
    semS = movmean(sem,20);

    x  = tSec(:);
    lo = (muS - semS).';
    hi = (muS + semS).';

    fill([x; flipud(x)], [lo; flipud(hi)], cols(cc,:), ...
        'FaceAlpha', 0.2, 'EdgeColor','none', 'HandleVisibility','off');
    plot(tSec, muS, 'LineWidth', 2, 'Color', cols(cc,:));

    % -------- Null ERP: shuffle t0 across breaths (break breath-event linkage) --------
    nullMu = zeros(1, numel(w));
    for ss = 1:nShuf
        t0sh = t0(randperm(nEv));  % shuffle within condition

        tmp = nan(nEv, numel(w));
        for ii = 1:nEv
            tmp(ii,:) = data(evB(ii), t0sh(ii) + w);
        end
        tmp = tmp - prctile(tmp(:), 2);
        tmp = tmp ./ prctile(tmp(:), 98);

        nullMu = nullMu + mean(tmp, 1, 'omitnan');
    end
    nullMu = nullMu / nShuf;
    nullMuS = movmean(nullMu, 20);

    plot(tSec, nullMuS, 'LineWidth', 1.5, 'Color', [0.6 0.6 0.6], ...
        'HandleVisibility','off');
end

xlabel('Time (s)'); ylabel('ERP (uV)');
legend(condNames, 'Location','best'); box off
title('Condition-specific mean ERP (mean ± SEM) + t0-shuffle null (grey)');
xlim([-.5 .5])




%% 



%%
% --- after you load chanDat ---
chanDat = add_OB_gamma_features_to_behDat(chanDat);
chanDat.behDat.subID = repmat({chanDat.subID}, height(chanDat.behDat), 1); 
chanDat.behDat.type = repmat({chanDat.type}, height(chanDat.behDat), 1); 
chanDat.behDat.sessType = repmat({chanDat.task}, height(chanDat.behDat), 1); 
chanDat.behDat.sessNum = repmat({chanDat.sessNum}, height(chanDat.behDat), 1); 
% --- resave chanDat back to original file ---
saveDir = fullfile(stem,'CHANDAT_processed');
save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');

% --- save behDat out (mat + csv) ---
[~, baseName] = fileparts(chanFiles(filei).name);

saveDirBehDat = fullfile(stem,'CHANDAT_processed','BehFiles');
if ~exist(saveDirBehDat,'dir'), mkdir(saveDirBehDat); end

behFileCsv = fullfile(saveDirBehDat, sprintf('%s_behFile.csv', baseName));


writetable(chanDat.behDat, behFileCsv);



end


% ============================================================
% Helper function (place at end of your script, or in its own .m file)
% ============================================================
function R = leadLagCorrTrials(X, Y, offsets, centers, winLen)
    % X, Y: [nTrials x nTime]
    % offsets: vector of sample offsets applied to Y relative to X
    %          (idxY = idxX + offset; negative => Y earlier)
    % centers: vector of window centers (sample indices, 1-based)
    % winLen:  window length (must be odd; e.g. 501)
    %
    % Returns:
    %   R: [nTrials x nOffsets x nCenters] row-wise Pearson r

    if nargin < 5 || isempty(winLen), winLen = 501; end
    if mod(winLen,2) ~= 1
        error('winLen must be odd (e.g., 501).');
    end

    [nTrials, nT] = size(X);
    if ~isequal(size(Y), [nTrials, nT])
        error('X and Y must have the same size [nTrials x nTime].');
    end

    half = (winLen - 1)/2;
    offsets = offsets(:);
    centers = centers(:);

    nO = numel(offsets);
    nC = numel(centers);

    R = nan(nTrials, nO, nC);

    for ci = 1:nC
        c = centers(ci);
        idxX = (c-half):(c+half);

        % If X window itself is out of bounds, all offsets at this center are invalid
        if idxX(1) < 1 || idxX(end) > nT
            continue
        end

        segX = X(:, idxX);

        % Pre-center X once per center (saves a little time)
        mx = mean(segX, 2, 'omitnan');
        X0 = segX - mx;
        ssX = sum(X0.^2, 2, 'omitnan');  % [nTrials x 1]

        for oi = 1:nO
            d = offsets(oi);
            idxY = idxX + d;

            if idxY(1) < 1 || idxY(end) > nT
                continue
            end

            segY = Y(:, idxY);

            my  = mean(segY, 2, 'omitnan');
            Y0  = segY - my;
            ssY = sum(Y0.^2, 2, 'omitnan');

            num = sum(X0 .* Y0, 2, 'omitnan');
            den = sqrt(ssX .* ssY);

            r = num ./ den;
            r(den == 0) = NaN;

            R(:, oi, ci) = r;
        end
    end
end

function summ = summarize_channel_gamma(chanDat, epochNames, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Summarize gamma characteristics for channel selection

useVec = chanDat.use(:)==1;
nBreaths = numel(useVec);
nEpochs = numel(epochNames);

summ.subID   = string(chanDat.subID);
summ.subIDSafe = sanitize_for_filename(summ.subID);

if isfield(chanDat,'sessNum') && ~isempty(chanDat.sessNum) && isfinite(chanDat.sessNum)
    summ.sessNum = double(chanDat.sessNum);
else
    summ.sessNum = 1;
end

% Channel label for filenames
try
    cl = string(chanDat.labels{chanDat.chi});
catch
    cl = "chan" + string(chanDat.chi);
end
summ.chanLabel = cl;
summ.chanLabelSafe = sanitize_for_filename(cl);

% ----- gamma peak detection matrix -----
G = chanDat.fooof.gamma_peaks; % breaths x 5 (NaN = no peak)
G = double(G);

% Prevalence (QC breaths only)
prev = nan(1,nEpochs);
cnt  = zeros(1,nEpochs);
den  = sum(useVec);
for e = 1:nEpochs
    good = useVec & isfinite(G(:,e));
    cnt(e) = sum(good);
    prev(e) = cnt(e) / max(den,1);
end
summ.prev_by_epoch = prev;
summ.count_by_epoch = cnt;

% Frequency MAD (QC breaths only, pooled and per epoch)
madEpoch = nan(1,nEpochs);
for e=1:nEpochs
    x = G(useVec & isfinite(G(:,e)), e);
    madEpoch(e) = robust_mad(x);
end
summ.madFreq_by_epoch = madEpoch;

xAll = G(useVec & isfinite(G));
summ.madFreq_pooled = robust_mad(xAll);

% ----- prominence (QC breaths only; pooled + per epoch) -----
[promEpoch, promAll] = compute_gamma_prominence(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);
summ.prom_by_epoch = promEpoch;
summ.prom_pooled   = promAll;

% ----- burst quality (QC breaths only) -----
if isfield(chanDat,'gammaBurst') && isfield(chanDat.gammaBurst,'prominence')
    gbProm = double(chanDat.gammaBurst.prominence(:));
    gbSNR  = [];
    if isfield(chanDat.gammaBurst,'snr')
        gbSNR = double(chanDat.gammaBurst.snr(:));
    end
    have = useVec & isfinite(gbProm);
    summ.burstProm_median = median(gbProm(have), 'omitnan');
    summ.burstCoverage    = sum(have)/max(sum(useVec),1);
    if ~isempty(gbSNR)
        summ.burstSNR_median = median(gbSNR(useVec & isfinite(gbSNR)), 'omitnan');
    else
        summ.burstSNR_median = NaN;
    end
else
    summ.burstProm_median = NaN;
    summ.burstCoverage    = 0;
    summ.burstSNR_median  = NaN;
end

summ.score = NaN; % filled later
summ.filePath = "";
summ.fileName = "";

end


function [winnerIdx, T, chanSumm] = pick_winner_channel(chanSumm)
% Robust within-group scoring to select winner

n = numel(chanSumm);
chanLabel = strings(n,1);
subID = strings(n,1);
sessNum = nan(n,1);

prevPooled = nan(n,1);
promPooled = nan(n,1);
madFreq    = nan(n,1);
burstProm  = nan(n,1);
burstCov   = nan(n,1);

for i=1:n
    chanLabel(i) = chanSumm(i).chanLabelSafe;
    subID(i)     = chanSumm(i).subIDSafe;
    sessNum(i)   = chanSumm(i).sessNum;

    prevPooled(i) = mean(chanSumm(i).prev_by_epoch, 'omitnan');
    promPooled(i) = chanSumm(i).prom_pooled;
    madFreq(i)    = chanSumm(i).madFreq_pooled;
    burstProm(i)  = chanSumm(i).burstProm_median;
    burstCov(i)   = chanSumm(i).burstCoverage;
end

% Robust z (median/MAD)
zPrev  = robust_z(prevPooled);
zProm  = robust_z(promPooled);
zMadF  = robust_z(madFreq);
zBProm = robust_z(burstProm);
zBCov  = robust_z(burstCov);

score = zPrev + zProm;

% Package table sorted
T = table(chanLabel, subID, sessNum, prevPooled, promPooled, madFreq, burstProm, burstCov, ...
    zPrev, zProm, zMadF, zBProm, zBCov, score);

[~, ord] = sort(score, 'descend', 'MissingPlacement','last');
T = T(ord,:);

% Winner = first row
winnerChanLabel = T.chanLabel(1);

% map back to index in chanSumm
winnerIdx = find(strcmp(chanLabel, winnerChanLabel), 1, 'first');

% Save score into chanSumm (optional)
for i=1:n
    chanSumm(i).score = score(i);
end

end


function out = sanitize_for_filename(s)
s = string(s);
out = regexprep(s, '[^\w\-]+', '_');     % keep letters/numbers/_/-
out = regexprep(out, '_+', '_');
out = strip(out, "_");
if strlength(out)==0, out="NA"; end
end


function m = robust_mad(x)
x = x(:);
x = x(isfinite(x));
if isempty(x), m = NaN; return; end
med = median(x);
m = median(abs(x-med));
end

function [promEpoch, promPooled] = compute_gamma_prominence(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Returns median prominence per epoch and pooled median across all epochs
[promByEpoch, promAll] = compute_gamma_prominence_full(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);
nEpochs = numel(promByEpoch);

promEpoch = nan(1,nEpochs);
for e=1:nEpochs
    promEpoch(e) = median(promByEpoch{e}, 'omitnan');
end
promPooled = median(promAll, 'omitnan');
end


function [promByEpoch, promAll, promRaw] = compute_gamma_prominence_full(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Returns per-epoch prominence samples (cell), pooled sample vector, and raw per-breath×epoch matrix

Gdet = double(chanDat.fooof.gamma_peaks);        % NaN if no peak
Gmax = double(chanDat.fooof.gamma_peak_freq);    % fallback
flat = double(chanDat.fooof.spectra_flat_log10); % breaths x epochs x frex
frex = double(chanDat.tf.frex(:));

nBreaths = size(Gdet,1);
nEpochs  = size(Gdet,2);

% Choose peak frequency per breath/epoch: detected if possible else max
Guse = Gdet;
missing = ~isfinite(Guse);
Guse(missing) = Gmax(missing);

% Precompute masks for baseline band
maskBaseBand = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);

promRaw = nan(nBreaths, nEpochs);
for b=1:nBreaths
    if ~useVec(b), continue; end
    for e=1:nEpochs
        fpk = Guse(b,e);
        if ~isfinite(fpk), continue; end
        [~, fiPk] = min(abs(frex - fpk));

        % baseline mask excludes +/- excludeHzAroundPeak around peak
        maskEx = maskBaseBand & ~(frex >= (fpk-excludeHzAroundPeak) & frex <= (fpk+excludeHzAroundPeak));

        base = median(flat(b,e,maskEx), 3, 'omitnan');
        if ~isfinite(base), continue; end

        promRaw(b,e) = flat(b,e,fiPk) - base;
    end
end

promByEpoch = cell(1,nEpochs);
promAll = [];
for e=1:nEpochs
    x = promRaw(useVec, e);
    x = x(isfinite(x));
    promByEpoch{e} = x(:);
    promAll = [promAll; x(:)];
end

end


function z = robust_z(x)
x = x(:);
med = median(x,'omitnan');
m = robust_mad(x);
if ~isfinite(m) || m==0
    z = (x - med);
else
    z = (x - med) ./ (1.4826*m);
end
end

% % gammaEnv_build(chanDat, stem, chanFiles(filei).name, ...
% %     'W', 100, 'pctl', 98, 'nShuf', 1000, 'rngSeed', 0, ...
% %     'rawArtifactThresh', 100, 'verbose', false);
% % 
% % 
% % 
% % % paths
% % [~, baseName] = fileparts(chanFiles(filei).name);
% % saveDirTf   = fullfile(stem,'CHANDAT_processed','tf_files');
% % gamFoo   = chanDat.fooof.gamma_peaks;  
% % gamMax  = chanDat.fooof.gamma_peak_freq;   
% % 
% % useVec   = double(chanDat.use(:)) == 1;
% % 
% % nBreaths = size(gamFoo,1);
% % frex = chanDat.tf.frex; 
% % powFile = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, 1));
% % S = load(powFile, 'tfInfo');
% % gammaEnv = struct; 
% % 
% % ampFrex = logspace(log10(2), log10(15), 20); 
% % allPow = nan([flip(size(S.tfInfo.pow)), 20]); 
% % allPhase = allPow; 
% % allGam = nan(nBreaths, size(S.tfInfo.pow,1));
% % gamFreq = nan(nBreaths, 1); 
% % for ii = 1:nBreaths
% % 
% %     idx = find(isnan(gamFoo(ii,:)));
% %     gamFoo(ii,idx) = gamMax(ii,idx); 
% %     curGam = median(gamFoo(ii,:)); 
% %     gamFreq(ii) = curGam; 
% %     fi = find(min(abs(frex - curGam)) == abs(frex - curGam), 1);
% %     if ~isempty(fi)
% % 
% %     powFile = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, fi));
% %     S = load(powFile, 'tfInfo');
% %     powLow = double(S.tfInfo.pow);   % [nT x nBreaths]
% %     allGam(ii,:) = powLow(:,ii); 
% %     powLow = log(powLow); 
% % 
% %  % extract power and phase information
% %     powLow = powLow(:,ii) - mean(powLow(:,ii)); 
% %     [phase,pow] = multiphasevec3(ampFrex,...
% %                     powLow,S.tfInfo.fs,6, true);
% %     allPow(ii, :, :) = squeeze(pow)';
% %     allPhase(ii,:,:) = squeeze(phase)';
% % 
% %     end
% % 
% % end
% % 
% % gammaEnv.pow = nan([size(chanDat.targIDX), length(ampFrex)]); 
% % gammaEnv.phase = nan([size(chanDat.targIDX), length(ampFrex)]);
% % gammaEnv.gamEnv = nan([size(chanDat.targIDX)]); 
% % 
% % for ii = 1:nBreaths
% %     curIDX = round(chanDat.targIDX(ii,:) * (S.tfInfo.fs/ chanDat.fs));
% %     if ~isnan(curIDX(1))
% %         curPow   = squeeze(allPow(ii,:,:)); 
% %         curPhase = squeeze(allPhase(ii,:,:));
% %         curEnv   = squeeze(allGam(ii,:));
% %         gammaEnv.pow(ii,:,:)   = curPow(curIDX,:); 
% %         gammaEnv.phase(ii,:,:) = curPhase(curIDX,:);
% %         gammaEnv.gamEnv(ii,:,:) = curEnv(curIDX); 
% %     end
% % end
% % 
% % 
% % 
% % gammaEnv.gamFreq = gamFreq; 
% % gammaEnv.frex = ampFrex; 
% % gammaEnv.fs = S.tfInfo.fs; 
% % 
% % 
% % % ============================================================
% % % Gamma-power up-crossings @ 97th percentile -> ITPC heatmap
% % % (TF timebase is 100 Hz; chanDat.fs is 500 Hz)
% % % ============================================================
% % W = 100;                 % +/- samples around crossing (IN TF SAMPLES)
% % pctl = 98;              % percentile threshold
% % 
% % fs_tf = chanDat.fs / 5; % 100 Hz (timebase of S.tfInfo.pow, allGam, allPhase)
% % 
% % % ---- 1) threshold (global across all use==1 breaths & timepoints) ----
% % mUse = useVec(:)==1;
% % x = allGam(mUse, :);
% % x = x(isfinite(x));                 % drop NaNs/Infs
% % thr = prctile(x, pctl);
% % 
% % % ---- 2) find up-crossings per breath ----
% % evBreath = [];
% % evT0     = [];
% % 
% % nT = size(allGam,2);
% % for ii = 1:nBreaths
% %     if ~mUse(ii), continue; end
% %     g = allGam(ii,:);
% %     if all(~isfinite(g)), continue; end
% % 
% %     above = g > thr;
% % 
% %     % rising edge: (t-1)<=thr and (t)>thr
% %     t0 = find( above(2:end) & ~above(1:end-1) ) + 1;
% % 
% %     for tt = 1:length(t0)
% % 
% %         t1 = t0(tt) + find( above(t0(tt):end-1) & ~above(t0(tt)+1:end),1) - 1;
% %         if isempty(t1)
% %             t1 = length(g); 
% %         end
% %         [~, t1] = max(g(t0(tt):t1));  
% %         t0(tt) = t1 + t0(tt) - 1; 
% %     end
% %     % keep only events with full +/-W window in-bounds
% %     t0 = t0(t0 > W & t0 <= (nT - W));
% % 
% % 
% %     rawVals = chanDat.trial.data(ii, t0*5); 
% %     t0(abs(rawVals)>100) = [] ; 
% % 
% %     if ~isempty(t0)
% %         evBreath = [evBreath; repmat(ii, numel(t0), 1)];
% %         evT0     = [evT0;     t0(:)];
% %     end
% % end
% % 
% % fprintf('Found %d up-crossings above %g (97th pct).\n', numel(evT0), thr);
% % 
% % % ---- 3) extract phase snippets around each event ----
% % nE   = numel(evT0);
% % nF   = size(allPhase,3);            % should be 20 (ampFrex)
% % nOff = 2*W + 1;
% % offs = -W:W;
% % 
% % if nE == 0
% %     warning('No events found. Consider lowering percentile or checking allGam values.');
% % else
% %     phSnip = nan(nE, nOff, nF);
% %     powSnip= nan(nE, nOff, nF); 
% %     rawSnip= nan(nE, nOff); 
% %     for e = 1:nE
% %         ii  = evBreath(e);
% %         t0  = evT0(e);
% %         idx = (t0-W):(t0+W);
% %         phSnip(e,:,:) = allPhase(ii, idx, :);  % [1 x nOff x nF]
% %         powSnip(e,:,:)= allPow(  ii, idx, :); 
% %         curTrial = chanDat.trial.data(ii,1:5:end); 
% %         rawSnip(e,:)  = curTrial(idx); 
% %     end
% % 
% %     % ---- 4) ITPC across events ----
% %     % ITPC(t,f) = |mean(exp(1j*phase), events)|
% %     itpc = abs(squeeze(mean(exp(1j*phSnip), 1, 'omitnan')));  % [nOff x nF]
% %     powM = squeeze(mean(powSnip, 1, 'omitnan')); 
% %     % ---- 5) plot heatmap ----
% %     tOff_sec = offs ./ fs_tf;       % TF timebase (100 Hz)
% %     figure;
% %     subplot 121
% %     imagesc(tOff_sec, [], itpc'); axis xy;
% %     yticks(2:2:20)
% %     yticklabels(ampFrex(2:2:20))
% %     colorbar;
% %     xlabel('Time offset from up-crossing (s)');
% %     ylabel('Modulation frequency (Hz)');
% %     title(sprintf('ITPC aligned to gamma-power %dth%% up-crossings (n=%d)', pctl, nE));
% %     subplot 122
% %     imagesc(tOff_sec, [], powM'); axis xy;
% %     yticks(2:2:20)
% %     yticklabels(ampFrex(2:2:20))
% %     colorbar;
% %     xlabel('Time offset from up-crossing (s)');
% %     ylabel('Modulation frequency (Hz)');
% %     title(sprintf('pow of modulation aligned to gamma-power %dth%% up-crossings (n=%d)', pctl, nE));
% % end
% % 
% % % ---- optional: save the event indices for later ----
% % gammaEnv.phSnip   = phSnip; 
% % gammaEnv.powSnip  = powSnip;
% % gammaEnv.thr      = thr;
% % gammaEnv.evBreath = evBreath;
% % gammaEnv.evT0     = evT0;
% % gammaEnv.W        = W;
% % gammaEnv.offs     = offs;
% % if exist('itpc','var'), gammaEnv.itpc = itpc; end
% % 
% % % ============================================================
% % % Surrogate ITPC distribution by shuffling detected timepoints
% % % across breaths, then z-score observed ITPC vs shuffles
% % % Assumes you already computed:
% % %   evBreath, evT0, W, offs, allPhase, ampFrex, useVec
% % % And that TF timebase is fs_tf = chanDat.fs/5 (100 Hz)
% % % ============================================================
% % 
% % fs_tf  = chanDat.fs / 5;     % 100 Hz TF sampling rate
% % nShuf  = 1000;
% % rng(0);                      % reproducible
% % 
% % % ---- quick guards ----
% % nBreaths = size(allPhase,1);
% % nT       = size(allPhase,2);
% % nF       = size(allPhase,3);
% % nE       = numel(evT0);
% % nOff     = 2*W + 1;
% % offs     = -W:W;
% % 
% % if nE == 0
% %     error('No events (evT0) available. Lower threshold or verify detection.');
% % end
% % 
% % % ---- precompute time-index matrix for the window around each event ----
% % idxMat = evT0(:) + offs;  % [nE x nOff] via implicit expansion (R2016b+)
% % % If older MATLAB, use: idxMat = bsxfun(@plus, evT0(:), offs);
% % 
% % % ---- helper: compute ITPC given a breath vector bVec (length nE) ----
% % % Uses only a loop over frequencies (20), no per-event loop.
% % compute_itpc = @(bVec) local_compute_itpc_from_allPhase(allPhase, bVec, idxMat, nBreaths, nT, nF);
% % 
% % % ---- 1) observed ITPC (using original event breath assignments) ----
% % itpc_obs = compute_itpc(evBreath(:));     % [nOff x nF]
% % 
% % % ---- 2) shuffles ----
% % itpc_shuf = nan(nOff, nF, nShuf, 'single');
% % 
% % for ss = 1:nShuf
% %     if mod(ss,100)==0, fprintf('Shuffle %d / %d\n', ss, nShuf); end
% % 
% %     % shuffle breath assignments across events, keep the detected timepoints fixed
% %     bShuf = evBreath(randperm(nE));
% % 
% %     itpc_shuf(:,:,ss) = single(compute_itpc(bShuf));
% % end
% % 
% % mu  = mean(itpc_shuf, 3, 'omitnan');
% % sig = std(itpc_shuf, 0, 3, 'omitnan');
% % 
% % % avoid divide-by-zero
% % sig(sig < 1e-6) = 1e-6;
% % 
% % z_itpc = (itpc_obs - double(mu)) ./ double(sig);   % [nOff x nF]
% % 
% % % ---- 3) plots ----
% % tOff_sec = offs ./ fs_tf;
% % 
% % figure;
% % 
% % % Observed ITPC
% % subplot(1,2,1);
% % imagesc(tOff_sec, [], itpc_obs'); axis xy;
% % yticks(2:2:20)
% % yticklabels(ampFrex(2:2:20))
% % colorbar;
% % xlabel('Time offset from up-crossing (s)');
% % ylabel('Modulation frequency (Hz)');
% % title(sprintf('Observed ITPC (n=%d events)', nE));
% % 
% % % Z-scored ITPC vs shuffle distribution
% % subplot(1,2,2);
% % imagesc(tOff_sec, [], z_itpc'); axis xy;
% % yticks(2:2:20)
% % yticklabels(ampFrex(2:2:20))
% % colorbar;
% % xlabel('Time offset from up-crossing (s)');
% % ylabel('Modulation frequency (Hz)');
% % title(sprintf('ITPC z-score vs %d shuffles', nShuf));
% % 
% % % ---- optional: stash results ----
% % gammaEnv.itpc_obs   = itpc_obs;
% % % gammaUpCross.itpc_shuf  = itpc_shuf;   % can be large; comment out if you prefer
% % gammaEnv.itpc_mu    = mu;
% % gammaEnv.itpc_sig   = sig;
% % gammaEnv.itpc_z     = z_itpc;
% % gammaEnv.nShuf      = nShuf;
% % gammaEnv.tim        = S.tfInfo.tim; 

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

