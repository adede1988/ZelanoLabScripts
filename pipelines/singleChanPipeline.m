function [] = singleChanPipeline(chanFiles, filei, subFiles, datPre)

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
    try
        QC = load([datPre  'cleanFiles/' chanDat.QCFileName]).cleanDat;
    catch
        QCnameBits = strsplit(chanDat.QCFileName, 'cleaningVecs.mat'); 
        chanDat.QCFileName = [QCnameBits{1} char(chanDat.task) ...
                            '_cleaningVecs.mat'];
        QC = load([datPre 'cleanFiles/' chanDat.QCFileName]).cleanDat;

    end
catch
    chanDat = load([chanFiles(filei).folder '/' chanFiles(filei).name]).chanDat; % go raw if it's not working!
    try
        QC = load([datPre 'cleanFiles/' chanDat.QCFileName]).cleanDat;
    catch
        QCnameBits = strsplit(chanDat.QCFileName, 'cleaningVecs.mat'); 
        chanDat.QCFileName = [QCnameBits{1} char(chanDat.task) ...
                            '_cleaningVecs.mat'];
        QC = load([datPre 'cleanFiles/' chanDat.QCFileName]).cleanDat;

    end
end

 try
        chanDat.behDat.length;
    disp('length checked')
 catch
        lengthVals = (lm.winEnd - lm.onsetIdx) ./ chanDat.fs; 
        chanDat.behDat.length = lengthVals; 
        disp('length added')
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

    % ---------------- USER SETTINGS ----------------
    fitRangeHz = [2 150];
    excludeLineHz = [60 120];        % [] to disable
    excludeHalfWidthHz = 2;

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

    % ---------------- Peak frequencies (same logic as you had) ----------------
    % Gamma peak (25–60 Hz) using spectra_flat ratio (QC breaths were fit; others are NaN)
    gammaMask = frex >= 25 & frex <= 60;
    tmp = chanDat.fooof.spectra_flat(:,:,gammaMask);   % [nBreaths x 5 x nGamma]
    allNan = all(isnan(tmp),3);
    tmp(isnan(tmp)) = -Inf;
    [~, imax] = max(tmp, [], 3);
    gammaFreqs = frex(gammaMask);
    gamma_peak_freq = gammaFreqs(imax);
    gamma_peak_freq(allNan) = NaN;
    chanDat.fooof.gamma_peak_freq = gamma_peak_freq;

    % Low-frequency peak (0.5–14 Hz) using spectra_flat ratio
    lowMask = frex >= .5 & frex <= 14;
    tmp = chanDat.fooof.spectra_flat(:,:,lowMask);
    allNan = all(isnan(tmp),3);
    tmp(isnan(tmp)) = -Inf;
    [~, imax] = max(tmp, [], 3);
    lowFreqs = frex(lowMask);
    low_peak_freq = lowFreqs(imax);
    low_peak_freq(allNan) = NaN;
    chanDat.fooof.low_peak_freq = low_peak_freq;

    % Save
    saveDir = fullfile(stem,'CHANDAT_processed');
    save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');
end


%% gamma bursts =======================
% STEP A: Gamma-burst event detection (1 peak per breath)
% Uses gamma_peak_freq(:,2) as the breath-wise gamma frequency (Hz)
% Loads the corresponding *_pow_fi###.mat file per breath, smooths with 100 ms FWHM Gaussian,
% then finds peak within [inhaleOnset-200ms, endOfBreath].
% Stores timing + QC metrics into chanDat.gammaBurst
% =======================

if ~isfield(chanDat, 'gammaBurst')
disp('working on gammaBurst')
baselineBandHz = [25 58]; gammaBandHz = [25 60]; 
[nBreaths, nEpochs, nF] = size(chanDat.tf.spectra);

prom = nan(nBreaths, nEpochs);     % spectral prominence (log10 units)
promF= nan(nBreaths, nEpochs); 
fiPk = nan(nBreaths, nEpochs);     % nearest frex index for the peak
pkHz = nan(nBreaths, nEpochs);     % snapped peak Hz

baseBandMask = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);
gammaMaskAll = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);
excludeHzAroundPeak = 2; 
useVec = chanDat.use == 1; 
G = chanDat.fooof.gamma_peak_freq;
flat_log10 = chanDat.fooof.spectra_flat_log10;
for e = 1:nEpochs
    for b = 1:nBreaths
        if ~useVec(b), continue; end
        f0 = G(b,e);
        if ~isfinite(f0), continue; end

        [~,ii] = min(abs(frex - f0));
        fiPk(b,e) = ii;
        pkHz(b,e) = frex(ii);

        exclMask = frex >= (pkHz(b,e)-excludeHzAroundPeak) & frex <= (pkHz(b,e)+excludeHzAroundPeak);
        bmask = baseBandMask & ~exclMask;

        y = squeeze(flat_log10(b,e,:));
        if ~any(isfinite(y(gammaMaskAll))), continue; end

        base = median(y(bmask), 'omitnan');
        if ~isfinite(base), continue; end

        prom(b,e) = y(ii) - base;
    end
end

GF = chanDat.fooof.gamma_peaks;
for e = 1:nEpochs
    for b = 1:nBreaths
        if ~useVec(b), continue; end
        f0 = GF(b,e);
        if ~isfinite(f0), continue; end

        [~,ii] = min(abs(frex - f0));
        fiPk(b,e) = ii;
        pkHz(b,e) = frex(ii);

        exclMask = frex >= (pkHz(b,e)-excludeHzAroundPeak) & frex <= (pkHz(b,e)+excludeHzAroundPeak);
        bmask = baseBandMask & ~exclMask;

        y = squeeze(flat_log10(b,e,:));
        if ~any(isfinite(y(gammaMaskAll))), continue; end

        base = median(y(bmask), 'omitnan');
        if ~isfinite(base), continue; end

        promF(b,e) = y(ii) - base;
    end
end

%choose the peak gamma frequency 
%choose the peak gamma frequency
gamma_peak_freq = nan(nBreaths,1);
inputPhaseID    = nan(nBreaths,1);
fooofBased      = false(nBreaths,1);   % true if taken from GF, false if from G (or no peak)

for b = 1:nBreaths
    % if ~useVec(b)
    %     continue
    % end

    % ---------- 1) Try FOOOF-based peaks (GF / promF) ----------
    gf = GF(b,:);                       % 1 x nEpochs
    ok = isfinite(gf);

    if any(ok)
        if nnz(ok) == 1
            eSel = find(ok, 1, 'first');
        else
            cand = find(ok);
            p    = promF(b,cand);
            p(~isfinite(p)) = -Inf;     % so NaN prominences never win
            [~,k] = max(p);
            eSel = cand(k);
        end

        gamma_peak_freq(b) = gf(eSel);
        inputPhaseID(b)    = eSel;
        fooofBased(b)      = true;
        continue
    end

    % ---------- 2) Fall back to non-FOOOF peaks (G / prom) ----------
    g  = G(b,:);
    ok = isfinite(g);

    if any(ok)
        if nnz(ok) == 1
            eSel = find(ok, 1, 'first');
        else
            cand = find(ok);
            p    = prom(b,cand);
            p(~isfinite(p)) = -Inf;
            [~,k] = max(p);
            eSel = cand(k);
        end

        gamma_peak_freq(b) = g(eSel);
        inputPhaseID(b)    = eSel;
        fooofBased(b)      = false;
    end
end


nTrials = size(gamma_peak_freq,1);

[~, baseName] = fileparts(chanFiles(filei).name);
saveDir = fullfile(stem,'CHANDAT_processed','tf_files');


gammaFreqVec = gamma_peak_freq;

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
gammaBurst.fooofBased      = fooofBased;
gammaBurst.inputPhase      = inputPhaseID; 

% metadata
gammaBurst.fs_tf           = fs_tf;
gammaBurst.dsFactor        = dsFact;
gammaBurst.smooth_fwhm_sec = fwhm_sec;
gammaBurst.inhaleOnsetIdx  = inhaleOnsetIdx;
gammaBurst.search_preMs    = preMs;

% breath-wise gamma frequency (use ONLY epoch 2)
gammaFreqVec = gamma_peak_freq;

for tr = 1:nTrials
   
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
%=========================================================

if ~isfield(chanDat, 'gammaBurstSecondary')
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
    base = median(seg, 'omitnan');
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
end


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
% =========================================================

if ~isfield(chanDat, 'gammaLockTF')
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
gammaLockTF.powZ_primary    = single(nan(nBreaths, numel(offsets_out), nFsel));
gammaLockTF.phase_primary   = single(nan(nBreaths, numel(offsets_out), nFsel));
gammaLockTF.pow_secondary   = single(nan(nBreaths, numel(offsets_out), nFsel));
gammaLockTF.powZ_secondary  = single(nan(nBreaths, numel(offsets_out), nFsel));
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
    powz = myChanZscore(powMat, ...
            [find(Sp.tfInfo.tim>=-.450,1), find(Sp.tfInfo.tim>=-.050,1)],...
            chanDat.use); %z-score
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
    winPowZ= nan(nBreaths, numel(offsets_tf));
    winPh  = nan(nBreaths, numel(offsets_tf));
    tmpPow = powMat(:);
    tmpPowZ= powz(:); 
    tmpPh  = phaseMat(:);
    winPow(valid) = tmpPow(lin);
    winPowZ(valid)= tmpPowZ(lin); 
    winPh(valid)  = tmpPh(lin);

    gammaLockTF.pow_primary(:,:,kk)   = single(winPow);
    gammaLockTF.powZ_primary(:,:,kk)  = single(winPowZ); 
    gammaLockTF.phase_primary(:,:,kk) = single(winPh);

    % ---------- SECONDARY extraction ----------
    idxMat = bsxfun(@plus, t0_secondary, offsets_tf);
    valid = isfinite(idxMat) & idxMat>=1 & idxMat<=nT & ...
            bsxfun(@ge, idxMat, breathStartIdx) & bsxfun(@le, idxMat, breathEndClamped);

    idxSafe = idxMat; idxSafe(~valid) = 1;
    lin = sub2ind([nT nBreaths], idxSafe(valid), breathIdxMat(valid));

    winPow = nan(nBreaths, numel(offsets_tf));
    winPowZ= nan(nBreaths, numel(offsets_tf));
    winPh  = nan(nBreaths, numel(offsets_tf));
    winPow(valid) = tmpPow(lin);
    winPowZ(valid)= tmpPowZ(lin); 
    winPh(valid)  = tmpPh(lin);

    gammaLockTF.pow_secondary(:,:,kk)   = single(winPow);
    gammaLockTF.powZ_secondary(:,:,kk)  = single(winPowZ); 
    gammaLockTF.phase_secondary(:,:,kk) = single(winPh);

    % ======================================================
    % Build shuffled nulls IN MEMORY for this frequency, then save once
    % ======================================================
    powNull_primary   = single(nan(nShuf, nBreaths, numel(offsets_tf)));
    powZNull_primary   = single(nan(nShuf, nBreaths, numel(offsets_tf)));
    phaseNull_primary = single(nan(nShuf, nBreaths, numel(offsets_tf)));

    powNull_secondary   = single(nan(nShuf, nBreaths, numel(offsets_tf)));
    powZNull_secondary   = single(nan(nShuf, nBreaths, numel(offsets_tf)));
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
            wPowZ= nan(nBreaths, numel(offsets_tf));
            wPh  = nan(nBreaths, numel(offsets_tf));
            wPow(valid) = tmpPow(lin);
            wPowZ(valid)= tmpPowZ(lin);
            wPh(valid)  = tmpPh(lin);

            powNull_primary(ss,:,:)   = single(wPow);
            powZNull_primary(ss,:,:)  = single(wPowZ);
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
    
            wPow  = nan(nBreaths, numel(offsets_tf));
            wPowZ = nan(nBreaths, numel(offsets_tf));
            wPh   = nan(nBreaths, numel(offsets_tf));
    
            wPow(valid)  = tmpPow(lin);
            wPowZ(valid) = tmpPowZ(lin);
            wPh(valid)   = tmpPh(lin);
    
            powNull_secondary(ss,:,:)   = single(wPow);
            powZNull_secondary(ss,:,:)  = single(wPowZ);
            phaseNull_secondary(ss,:,:) = single(wPh);
        end
    end


    % ---- save ONE file for this frequency ----
    nullOut = struct();
    nullOut.primary = struct('powNull', powNull_primary,...
                             'powNullZ', powZNull_primary, ...
                             'phaseNull', phaseNull_primary);
    nullOut.secondary = struct('powNull', powNull_secondary, ...
                               'powNullZ', powZNull_secondary, ...
                               'phaseNull', phaseNull_secondary);

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

    % nullFile = fullfile(saveDirNull, sprintf('%s_BurstNull_fi%03d.mat', baseName, fi));
    % save(nullFile, 'nullOut', '-v7.3');

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
end
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
% 
% if ~isfield(chanDat, 'leadLag')
% disp('working on leadLag (Step D)')
% % --- inputs ---
% if ~isfield(chanDat,'gammaBurst') || ~isfield(chanDat.gammaBurst,'t0_idx')
%     error('chanDat.gammaBurst.t0_idx not found.');
% end
% if ~isfield(chanDat,'gammaBurst') || ~isfield(chanDat.gammaBurst,'fi')
%     error('chanDat.gammaBurst.fi not found. (Run gammaBurst block that stores fi.)');
% end
% 
% t0_idx   = double(chanDat.gammaBurst.t0_idx(:));  % TF (100 Hz) index, breaths x 1
% gammaFi  = double(chanDat.gammaBurst.fi(:));      % index into frex, breaths x 1
% useVec   = double(chanDat.use(:)) == 1;
% 
% nBreaths = numel(t0_idx);
% 
% % paths
% [~, baseName] = fileparts(chanFiles(filei).name);
% saveDirTf   = fullfile(stem,'CHANDAT_processed','tf_files');
% saveDirNull = fullfile(stem,'CHANDAT_processed','burstNullFiles');
% if ~exist(saveDirNull,'dir'), mkdir(saveDirNull); end
% 
% % --- TF sampling rate and timebase (assume consistent across tf files) ---
% % Use any valid gammaFi breath to load one file and get fs_tf and nT
% tr0 = find(isfinite(gammaFi) & isfinite(t0_idx), 1, 'first');
% if isempty(tr0)
%     warning('No valid breaths for leadLag; skipping.');
%     chanDat.leadLag = struct();
%     return
% end
% 
% powFile0 = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, gammaFi(tr0)));
% S0 = load(powFile0, 'tfInfo');
% fs_tf = double(S0.tfInfo.fs);         % should be 100
% nT    = size(S0.tfInfo.pow,1);        % e.g., 1200
% tfTim = double(S0.tfInfo.tim(:));     % optional, for metadata
% 
% % --- template length: ±500 ms around gamma peak at 100 Hz ---
% halfWin_sec = 0.500;
% halfL = round(halfWin_sec * fs_tf);   % 50 if fs_tf=100
% L = 2*halfL + 1;                      % 101
% 
% % inhale onset and breath end indices (TF samples)
% inhaleOnsetIdx = double(chanDat.gammaBurst.inhaleOnsetIdx);
% breathEndIdx = round(double(chanDat.behDat.length(:))*fs_tf + inhaleOnsetIdx);
% breathEndIdx = min(max(1, breathEndIdx), nT);
% 
% % --- target frequencies (Hz), map -> frex indices ---
% targetHz = [.2 .6 1.0 1.8 2.6 3.4 4.2 5 6 7 8 9 10 11 12 13 15 17 ...
%     19 21 24 27 30 35 40 45 50 55 65 75 85 95 105 115 130 145 160 175 190];
% nF = numel(targetHz);
% 
% targetFi = nan(1,nF);
% targetFrex = nan(1,nF);
% for k = 1:nF
%     [~, ii] = min(abs(frex - targetHz(k)));
%     targetFi(k) = ii;
%     targetFrex(k) = frex(ii);
% end
% 
% % --- output arrays (empirical) ---
% lag_sec  = nan(nBreaths, nF);
% corr_max = nan(nBreaths, nF);
% 
% % --- null setup ---
% nShuf = 1000;
% lagNull = nan(nShuf, nBreaths, nF, 'single');  % only fill for use==1 breaths
% 
% % For reproducibility
% rngState = rng;
% 
% % =========================================================
% % Build gamma templates per breath (raw gamma power snippet)
% % Store centered template gc and template std sg for fast correlation
% % =========================================================
% gammaTemplate_gc = nan(nBreaths, L);   % centered template (g - mean(g))
% gammaTemplate_sg = nan(nBreaths, 1);   % std(g) with N-1 normalization
% 
% validTemplate = isfinite(t0_idx) & isfinite(gammaFi) & (t0_idx-halfL>=1) & (t0_idx+halfL<=nT);
% 
% uFi = unique(gammaFi(validTemplate));
% uFi = uFi(isfinite(uFi));
% 
% for uf = uFi(:)'
%     powFileG = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, uf));
%     SG = load(powFileG, 'tfInfo');
%     powG = double(SG.tfInfo.pow);  % [nT x nBreaths]
% 
%     breathsHere = find(validTemplate & gammaFi==uf);
%     for tr = breathsHere(:)'
%         g = powG( (t0_idx(tr)-halfL):(t0_idx(tr)+halfL), tr );
%         if any(~isfinite(g)), continue; end
%         mg = mean(g);
%         gc = g - mg;
%         sg = std(g,0);  % unbiased (N-1)
%         if sg<=0 || ~isfinite(sg), continue; end
%         gammaTemplate_gc(tr,:) = gc(:)';
%         gammaTemplate_sg(tr)   = sg;
%     end
% end
% 
% haveTemplate = all(isfinite(gammaTemplate_gc),2) & isfinite(gammaTemplate_sg);
% 
% % =========================================================
% % MAIN LOOP: target frequency (load each low-freq file once)
% % =========================================================
% onesL = ones(L,1);
% 
% for k = 1:nF
%     fi = targetFi(k);
%     fprintf('leadLag: %d/%d (fi=%d, f=%.3f Hz)\n', k, nF, fi, frex(fi));
% 
%     % load low-frequency power file
%     powFile = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, fi));
%     S = load(powFile, 'tfInfo');
%     powLow = double(S.tfInfo.pow);   % [nT x nBreaths]
%     if size(powLow,1) ~= nT
%         % if length mismatch, clamp to min
%         nTcur = min(nT, size(powLow,1));
%         powLow = powLow(1:nTcur,:);
%     end
% 
%     for tr = 1:nBreaths
%         if ~haveTemplate(tr), continue; end
% 
%         t0 = t0_idx(tr);
%         if ~isfinite(t0), continue; end
% 
%         % Sliding rule: s = 1 .. sMax
%         sMax = min( (nT - L + 1), (breathEndIdx(tr) - halfL) );
%         sMax = floor(sMax);
%         if sMax < 1, continue; end
% 
%         % searchable portion length M so that linear starts 1..sMax use indices 1..M
%         M = sMax + L - 1;
%         if M > nT, M = nT; end
% 
%         ySeg = powLow(1:M, tr);
%         if any(~isfinite(ySeg)), continue; end
% 
%         % Build circular-extended signal to compute Pearson corr for ALL circular starts 1..M
%         yExt = [ySeg; ySeg(1:(L-1))];  % length M+L-1
% 
%         gc = gammaTemplate_gc(tr,:)';   % L x 1
%         sg = gammaTemplate_sg(tr);      % scalar
% 
%         % Numerator for start s: sum(gc .* y_window)
%         num = conv(yExt, flipud(gc), 'valid'); % length M
% 
%         % Local mean/std for each window (circular)
%         sumY  = conv(yExt, onesL, 'valid');    % length M
%         sumY2 = conv(yExt.^2, onesL, 'valid'); % length M
% 
%         % Unbiased local variance (N-1) form
%         varY = (sumY2 - (sumY.^2)/L) / (L-1);
%         varY(varY < 0) = 0;
%         sdY = sqrt(varY);
% 
%         den = (L-1) * sg .* sdY;
%         corr_circ = num ./ den;
%         corr_circ(~isfinite(corr_circ)) = NaN;
% 
%         % Empirical: only starts 1..sMax (no wrap)
%         cEmp = corr_circ(1:sMax);
%         if all(isnan(cEmp)), continue; end
% 
%         [cBest, pos] = max(cEmp, [], 'omitnan'); % pos is best start index (1..sMax)
%         if isempty(pos) || ~isfinite(cBest), continue; end
% 
%         bestCenter = pos + halfL; % trial index (since starts are in trial coords)
%         lag_sec(tr,k)  = (bestCenter - t0) / fs_tf;
%         corr_max(tr,k) = cBest;
% 
%         % Nulls: only for use==1 breaths
%         if ~useVec(tr), continue; end
% 
%         % Precompute best-start for EVERY possible circular shift kShift = 0..M-1
%         % For shift kShift, the correlation values for linear starts 1..sMax in the shifted signal are:
%         %   corrShift(s) = corr_circ( mod(s - kShift - 1, M) + 1 )
%         % We find argmax over s=1..sMax for each kShift.
% 
%         kVec = (0:(M-1))'; % M x 1
%         idxMat = mod((1:sMax) - kVec - 1, M) + 1;     % [M x sMax]
%         vals = corr_circ(idxMat);                      % [M x sMax]
%         vals(isnan(vals)) = -Inf;
% 
%         [maxVal, posMat] = max(vals, [], 2);          % [M x 1]
%         posMat(maxVal==-Inf) = NaN;                   % shifts with no valid corr
% 
%         % draw 1000 random circular shifts (0..M-1)
%         kRand = randi([0 M-1], nShuf, 1);             % with replacement
%         bestStartNull = posMat(kRand+1);              % [nShuf x 1] start indices (1..sMax)
% 
%         % Convert to lag (sec): center = start + halfL
%         % (If bestStartNull is NaN, lag stays NaN)
%         lagNull(:,tr,k) = single((bestStartNull + halfL - t0) / fs_tf);
%     end
% end
% 
% % =========================================================
% % Store empirical outputs in chanDat.leadLag
% % =========================================================
% leadLag = struct();
% leadLag.fs_tf = fs_tf;
% leadLag.halfWin_sec = halfWin_sec;
% leadLag.L = L;
% leadLag.halfL = halfL;
% 
% leadLag.targetHz = targetHz;
% leadLag.targetFi = targetFi;
% leadLag.targetFrex = targetFrex;
% 
% leadLag.t0_idx = t0_idx;
% leadLag.gammaFi = gammaFi;
% 
% leadLag.lag_sec = lag_sec;
% leadLag.corr_max = corr_max;
% 
% chanDat.leadLag = leadLag;
% 
% % =========================================================
% % Save null lags to disk (separate file)
% % =========================================================
% nullOut = struct();
% nullOut.lagNull = lagNull;   % [nShuf x nBreaths x nF] single
% nullOut.meta = struct();
% nullOut.meta.baseName = baseName;
% nullOut.meta.fs_tf = fs_tf;
% nullOut.meta.halfWin_sec = halfWin_sec;
% nullOut.meta.L = L;
% nullOut.meta.halfL = halfL;
% nullOut.meta.targetHz = targetHz;
% nullOut.meta.targetFi = targetFi;
% nullOut.meta.targetFrex = targetFrex;
% nullOut.meta.useBreaths = find(useVec);
% nullOut.meta.nShuf = nShuf;
% nullOut.meta.rngState = rngState;
% 
% nullFile = fullfile(saveDirNull, sprintf('%s_LeadLagNull_circShift_n%d.mat', baseName, nShuf));
% save(nullFile, 'nullOut', '-v7.3');
% 
% disp('leadLag complete: stored chanDat.leadLag and saved null lag cube')
% saveDir = fullfile(stem,'CHANDAT_processed');
% save(fullfile(saveDir, chanFiles(filei).name), 'chanDat', '-v7.3');
% end
%% get info for gamma envelope. This turned out not to be very helpful

if ~isfield(chanDat, 'gammaEnv')

    %get information about the gamma envelope
chanDat.gammaEnv = gammaEnv_build(chanDat, chanFiles, filei, stem, 100, 98, 1000);


end


%% respiration PAC
if ~isfield(chanDat, 'pac')
    disp('working on pac')
  rawDat = load([chanFiles(filei).folder '/' chanFiles(filei).name]).chanDat;
  keyBreathIDX = chanDat.targIDX; 
    onsets = chanDat.behDat.finalOnset; 
    keyBreathIDX = keyBreathIDX + onsets; 
    keyBreathIDX = keyBreathIDX(:,[ 1 25 50]);
  gamMed = median(chanDat.fooof.gamma_peaks, 'all', 'omitnan');
    fs = chanDat.fs;               % Hz
    
    halfBW = 5;                    % +/- 5 Hz
    bpHz   = double([gamMed-halfBW, gamMed+halfBW]);
    
    PACfrex = logspace(log10(.05), log10(2), 50); 
    
  [pacOut, meta] = pac_breathTemplate_timeResolvedPAC(rawDat, [], keyBreathIDX, gamMed, fs, bpHz, PACfrex, ...
    'targetFs', 20);
   meta = pac_addBreathDiagnostics(meta, pacOut, chanDat.behDat);
   meta.pac = pacOut; 
   chanDat.pac = meta; 
   out = helper_gammaExtremaNearPrefPhase(rawDat, meta, PACfrex);
   chanDat.pac_peaks = out; 

end

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

