function [] = singleChanPipeline_general(chanFiles, idx, subFiles, type)

%% set frequency parameters
frex = logspace(log10(2),log10(150),200);
numfrex = length(frex); 
stds = logspace(log10(3),log10(15),numfrex)./(2*pi*frex);

%HFB
highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = logspace(log10(10),log10(20),highnumfrex)./(2*pi*highfrex);

%% load the data
folderName = split(chanFiles(idx).folder, 'CHANDAT'); 
folderName = folderName{1}; 
try %try loading the processed file
    
    chanDat = load([folderName '/CHANDAT_processed/' ...
        chanFiles(idx).name]).chanDat; 
catch
    chanDat = load([chanFiles(idx).folder '/' ...
        chanFiles(idx).name]).out; % go raw if it's not working!
end

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])

%hard code that breath window is [-2:end) 
if ~isfield(chanDat, 'tim')
    chanDat.tim = -2:1/chanDat.fs:...
        size(chanDat.trialDat,1)/chanDat.fs - 2 - 1/chanDat.fs;
end

if strcmp(chanDat.task, 'breathing')
    chanDat.rsp(:,:,chanDat.behDat(:,13)==0) = []; 
    chanDat.RRints(:,chanDat.behDat(:,13)==0) = []; 
    chanDat.targTrace(:,chanDat.behDat(:,13)==0) = []; 
    chanDat.trialDat(:,chanDat.behDat(:,13)==0) = []; 
    chanDat.behDat(chanDat.behDat(:,13)==0,:) = []; 
end


%% breath phase info per trial: 

if ~isfield(chanDat, 'phaseindices')
    rsp = chanDat.rsp;
    if chanDat.rspIDX == 999
        rsp = squeeze(rsp(1,:,:));
    else
        rsp = squeeze(rsp(chanDat.rspIDX,:,:)); 
    end
    if strcmp(chanDat.task, 'breathing')
        rsp(:,chanDat.behDat(:,13)==0) = []; 
    end
    rspSmooth = smoothdata(rsp, 1, 'gaussian', round(chanDat.fs/2)); 
    rspPhase = angle(hilbert(rspSmooth));
    idx0        = find(chanDat.tim >= 0, 1, 'first');
    ph0         = rspPhase(idx0, :);     % phases in [-pi, pi]
    startPhase  = atan2(mean(sin(ph0)), mean(cos(ph0)));  % circular mean
    if isnan(startPhase)
        error('nan values in respiration')
    end
    n = 150; %length of breath cycle in standardized phases
    theta_open      = startPhase + (0:n-1)*(2*pi/n);
    targPhases = angle(exp(1j*theta_open));
    [~, phaseIndices] = match_phase_indices(rspPhase, idx0, ...
                                startPhase, targPhases); 

    chanDat.targPhases = targPhases; 
    chanDat.phaseindices = phaseIndices; 

end


%% fooof
if ~isfield(chanDat, 'bands')
    % Welch per trial, then average
    win     = round(2.0*chanDat.fs);         % 1 s window (tune as needed)
    nover   = round(0.5*win);
    nfft    = max(2^nextpow2(win), 2^nextpow2(chanDat.fs));  % decent freq grid
    
    for tr = 1:size(chanDat.trialDat,2)
        x = detrend(chanDat.trialDat(:,tr),'constant');  % remove DC
        [pxx(:,tr), f] = pwelch(x, hamming(win), nover, ...
                nfft, chanDat.fs, 'psd');
    end
    psd_mean = mean(pxx, 2);             % linear PSD
    
    % fits for low and high frequency to avoid non linearity in aperiodic psd
    res_low  = fooof_basic(f, psd_mean, ...
        'peak_thresh', 1, ...
        'f_range', [1 40], ...
        'aperiodic_mode','fixed', ...
        'max_peaks', 10, ...
        'peak_width_limits', [0.5 8], ...
        'verbose', false);
    
    res_high = fooof_basic(f, psd_mean, ...
        'peak_thresh', 1, ...
        'f_range', [20 200], ...
        'aperiodic_mode','fixed', ...
        'max_peaks', 6, ...
        'peak_width_limits', [2 20], ...
        'verbose', false);
    
    % --- band definitions (Hz) ---
    bands.theta = [4 8];
    bands.alpha = [8 13];
    bands.beta  = [13 30];
    bands.gamma1 = [30 58];
    bands.gamma2 = [62 80];
    bands.high  = [80 min(200, max(f))];   
    
    % --- extract one peak per band (largest by amplitude) ---
    bp = extract_band_peaks(res_low, res_high, bands);
    
    chanDat.bands = bp; 

end


%% single trial time frequency

if ~isfield(chanDat, 'TF')
    disp('working on TF')
    %to keep size down, don't put large variables into the chanDat struct!
    TFout = struct;
    pow = getChanTrialTF_highPrec(chanDat.trialDat, frex, ...
                 stds, chanDat.fs, 'UseParfor', true, 'WidthMode','std');

    %%% should there be some kind of noise trial detection and rejection
    %%% here? 

    phase = angle(pow); 
    pow = abs(pow).^2; 
    powz = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), ...
                    'UniformOutput',false ); %z-score
    powz = cell2mat(powz); %organize
    powz = reshape(powz, size(powz,1), size(powz,2)/numfrex, []); %organize
    
    %grab 50 Hz sampled output: 
    sampRatio = (chanDat.fs / 50); %convert to 50 Hz for output
    TFout.powz = powz(1:round(sampRatio):size(powz,1), :, :); 
    TFout.phase = phase(1:round(sampRatio):size(powz,1), :, :); 
    TFout.TF_tim = chanDat.tim(1:round(sampRatio):size(powz,1));

    TFout.pow = pow(1:round(sampRatio):size(pow,1), :, :); 

    %grab breath phase aligned data: 
    
    f = @(dat, idx1, idx2) dat(idx1(~isnan(idx1)), idx2, :); 
    TFout.breathPowz = arrayfun(@(x) f(powz, ...
                                      chanDat.phaseindices(:,x), ...
                                      x),...
            1:size(chanDat.phaseindices,2), 'uniformoutput', false);
    f = nan([length(chanDat.targPhases), 1, numfrex]); 
   
    TFout.breathPowz = cellfun(@(x) set_slice(f, x), TFout.breathPowz, ...
        'UniformOutput', false); 
    TFout.breathPowz = cat(2, TFout.breathPowz{:}); 

    f = @(dat, idx1, idx2) dat(idx1(~isnan(idx1)), idx2, :); 
    TFout.breathPow = arrayfun(@(x) f(pow, ...
                                      chanDat.phaseindices(:,x), ...
                                      x),...
            1:size(chanDat.phaseindices,2), 'uniformoutput', false);
    f = nan([length(chanDat.targPhases), 1, numfrex]); 
   
    TFout.breathPow = cellfun(@(x) set_slice(f, x), TFout.breathPow, ...
        'UniformOutput', false); 
    TFout.breathPow = cat(2, TFout.breathPow{:}); 


    chanDat.TF = TFout; 

    
    
    
    clear TFout 
     disp('attempting saving')
    save([folderName 'CHANDAT_processed/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name])
else
    disp('TF already done, skipping')
end





%% High frequency Broadband 


if ~isfield(chanDat, 'reactiveRes') && type == 'macro'
    HFB = getHFB(chanDat, highfrex); 

    chanDat.HFB = HFB; 
    chanDat.reactiveRes = reactiveTest_100(chanDat.HFB);
    % downsample for size
    HFB = chanDat.HFB;
    HFB_names = fieldnames(HFB); 
    
    datidx = {[1,2], [5,6], [9,10,11,12], [15,16,17,18]}; 
    timidx = [3,7,13,19]; 

    for dati = 1:4
        curtim = HFB.(HFB_names{timidx(dati)});
        HFB.(HFB_names{timidx(dati)}) = curtim(1:5:end);
        curDi = datidx{dati}; 
        for fi = 1:length(curDi)
            cur = HFB.(HFB_names{curDi(fi)});

            HFB.(HFB_names{curDi(fi)}) = cur(1:5:end,:);
        end
    end
    chanDat.HFB = HFB; 
    clear HFB 
    % done downsample 
    disp('attempting saving')
    save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name])

else
    disp('HFB already done')
end

%% get the HFB latencies and save them in order to reference other variables to them

if ~isfield(chanDat, 'HFB_lat')
    HFB_lat = struct; 
   

    %encoding
    HFB_lat.subHit = gausLat(chanDat.HFB.subHit, ...
        chanDat.HFB.encMulTim, ...
        chanDat.encInfo(chanDat.use & chanDat.hits, 4), ...
        1);

    HFB_lat.subMiss = gausLat(chanDat.HFB.subMiss, ...
        chanDat.HFB.encMulTim, ...
        chanDat.encInfo(chanDat.use & chanDat.misses, 4), ...
        1);

    %retrieval
    triali = chanDat.retInfo(:,1)==1;
    HFB_lat.retHit = gausLat(chanDat.HFB.hit_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);

    triali = chanDat.retInfo(:,1)==2;
    HFB_lat.retMiss = gausLat(chanDat.HFB.miss_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);
    
    triali = chanDat.retInfo(:,1)==3;
    HFB_lat.retCR = gausLat(chanDat.HFB.cr_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);

    triali = chanDat.retInfo(:,1)==4;
    HFB_lat.retFA = gausLat(chanDat.HFB.fa_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);

    chanDat.HFB_lat = HFB_lat; 



end











%% get ISPC and PPC values 

if ~isfield(chanDat, 'ISPC')
    disp('working on ISPC')
    ISPCout = struct; 
    %store the downsample index (di) 
    multim = chanDat.HFB.encMulTim; 
    ISPCout.encdi = arrayfun(@(x) find(x<=chanDat.enctim,1), multim);
    
    multim = chanDat.HFB.encRT_tim; 
    ISPCout.encRdi = arrayfun(@(x) find(x<=chanDat.enctimRT,1), multim);

    multim = chanDat.HFB.onMulTim; 
    ISPCout.ondi = arrayfun(@(x) find(x<=chanDat.retOtim,1), multim);
    
    multim = chanDat.HFB.rtMulTim; 
    ISPCout.rtdi = arrayfun(@(x) find(x<=chanDat.retRtim,1), multim);
    
    %preallocate: 
    %channels X time X frequencies X ISPC/PPC
    frex = logspace(log10(2), log10(25), 20); 
    numfrex = length(frex); 
    stds = logspace(log10(3),log10(5),numfrex)./(2*pi*frex);

    ISPCout.subMiss = zeros(length(subFiles), length(ISPCout.encdi), length(frex), 4); 
    ISPCout.subHit = zeros(length(subFiles), length(ISPCout.encdi), length(frex), 4); 
    ISPCout.subMissRT = zeros(length(subFiles), length(ISPCout.encRdi), length(frex), 4);
    ISPCout.subHitRT = zeros(length(subFiles), length(ISPCout.encRdi), length(frex), 4);
    
    
    ISPCout.hit_on = zeros(length(subFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.miss_on = zeros(length(subFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.cr_on = zeros(length(subFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.miss_onRT = zeros(length(subFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.hit_onRT = zeros(length(subFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.cr_onRT = zeros(length(subFiles), length(ISPCout.rtdi), length(frex), 4);




    %will need to loop channels
    %NOTE: all trial types must have at least two trials! 
    for chan = 1:length(subFiles)
        tic
       
        chanDat2 = load([subFiles(chan).folder '/' subFiles(chan).name]).chanDat; 

        %ENCODING DATA: ***********************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMiss(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, ...
            chanDat.use & chanDat.misses, chanDat.HFB_lat.subMiss, ...
            chanDat.HFB.encMulTim, true);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHit(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, ...
            chanDat.use & chanDat.hits, chanDat.HFB_lat.subHit, ...
            chanDat.HFB.encMulTim, true);
        end

     
        %ENCODING DATA BEHAVIOR RESPONSE: *************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMissRT(chan,:,:,:) = getChanISPC(chanDat.encRT, chanDat2.encRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encRdi, ...
            chanDat.use & chanDat.misses, chanDat.HFB_lat.subMiss, ...
            chanDat.HFB.encRT_tim, false);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHitRT(chan,:,:,:) = getChanISPC(chanDat.encRT, chanDat2.encRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encRdi, ...
            chanDat.use & chanDat.hits, chanDat.HFB_lat.subHit, ...
            chanDat.HFB.encRT_tim, false);
        end



        %RETRIEVAL STIM ONSET: ****************************************************
        if sum(chanDat.retInfo(:,1)==1) > 1
        ISPCout.hit_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, ...
            chanDat.retInfo(:,1)==1, chanDat.HFB_lat.retHit, ...
            chanDat.HFB.onMulTim, true);
        end
        if sum(chanDat.retInfo(:,1)==2) > 1
        ISPCout.miss_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, ...
            chanDat.retInfo(:,1)==2, chanDat.HFB_lat.retMiss, ...
            chanDat.HFB.onMulTim, true);
        end
        if sum(chanDat.retInfo(:,1)==3) > 1
        ISPCout.cr_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, ...
            chanDat.retInfo(:,1)==3, chanDat.HFB_lat.retCR, ...
            chanDat.HFB.onMulTim, true);
        end
     
        
        %RETRIEVAL BEHAVIOR RESPONSE: ****************************************************
        if sum(chanDat.retInfo(:,1)==1) > 1
        ISPCout.hit_onRT(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, ...
            chanDat.retInfo(:,1)==1, chanDat.HFB_lat.retHit, ...
            chanDat.HFB.rtMulTim, false);
        end
        if sum(chanDat.retInfo(:,1)==2) > 1
        ISPCout.miss_onRT(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, ...
            chanDat.retInfo(:,1)==2, chanDat.HFB_lat.retMiss, ...
            chanDat.HFB.rtMulTim, false);
        end
        if sum(chanDat.retInfo(:,1)==3) > 1
        ISPCout.cr_onRT(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, ...
            chanDat.retInfo(:,1)==3, chanDat.HFB_lat.retCR, ...
            chanDat.HFB.rtMulTim, false);
        end
     
        

        
        disp(['channel: ' num2str(chan) ' took ' num2str(round(toc/60,1)) ' minutes'])
    end
    chanDat.ISPC = ISPCout; 
    disp('attempting saving')
    save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name])



else
    disp('connectivity already done, skipping')
end




%% final save out

save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name]) 
% delete([chanFiles(idx).folder '/' chanFiles(idx).name]) 







end















