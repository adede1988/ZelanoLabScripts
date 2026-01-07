function outTTLs = detect_ttls_EmotionalMovie(raw, P)

    idx = cellfun(@(x) contains(x, 'event'), raw.labels);
    photoDiode = raw.data(idx, :); 
   
    wo = 60/(raw.fs_raw/2);                 % normalized center freq
    bw = wo/35;                     % Q=35 ~ 1.7 Hz 3-dB bandwidth at 60 Hz
    [b,a] = iirnotch(wo, bw);
    sig = filtfilt(b, a, double(photoDiode));   

    wo = 120/(raw.fs_raw/2);                 % normalized center freq
    bw = wo/35;                     % Q=35 ~ 1.7 Hz 3-dB bandwidth at 120 Hz
    [b,a] = iirnotch(wo, bw);
    sig = filtfilt(b, a, double(sig));  

    wo = 180/(raw.fs_raw/2);                 % normalized center freq
    bw = wo/35;                     % Q=35 ~ 1.7 Hz 3-dB bandwidth at 180 Hz
    [b,a] = iirnotch(wo, bw);
    sig = filtfilt(b, a, double(sig));  

    
    photoDiode = (sig - mean(sig)) / std(sig);
    

    TTLs = find(photoDiode(1:end-1)>P.pd.zthresh &...
                photoDiode(2:end)<P.pd.zthresh);
    diffVals = diff(TTLs); 
    TTLs(diffVals<100) = []; 
    

    outTTLs = table;
    outTTLs.timStamp = zeros(180,1); 
    outTTLs.type = zeros(180,1); 
    ii = 1;
    ti = 1; 
    while true
        
        curTTL = TTLs(ii); 
        outTTLs.timStamp(ti) = curTTL;
        if TTLs(ii+2) - curTTL < P.pd.searchWin
            outTTLs.type(ti) = P.pd.numNeg; %negative
            ii = ii + P.pd.numNeg; 
        elseif TTLs(ii+1) - curTTL < P.pd.searchWin
            outTTLs.type(ti) = P.pd.numPos; %positive
            ii = ii + P.pd.numPos; 
        else
            outTTLs.type(ti) = P.pd.numNeu; %neutral 
            ii = ii + P.pd.numNeu; 
        end
        ti = ti + 1; 
        if ii >= length(TTLs)-2
            break
        end
    
    end

    outTTLs(outTTLs.timStamp==0,:) = []; 
   

    figure('visible', false, 'position', [0,0,1000,500])
    plot(photoDiode)
    xline(outTTLs.timStamp(outTTLs.type==1), 'color', 'green')
    xline(outTTLs.timStamp(outTTLs.type==2), 'color', 'blue')
    xline(outTTLs.timStamp(outTTLs.type==3), 'color', 'red')
    title([raw.sessID ' TTLs'], 'interpreter', 'none')


    outTTLs.timStamp = round(outTTLs.timStamp ./ ...
                        (raw.fs_raw / P.fs_target));

end