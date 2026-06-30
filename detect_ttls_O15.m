function [outTTLs, raw] = detect_ttls_O15(raw, P)

    idx = cellfun(@(x) contains(x, 'event'), raw.labels);
    photoDiode = raw.data(idx, :); 
    

 % Find missing samples in the photodiode channel
        nanidx = isnan(photoDiode);
        nNan   = sum(nanidx);
        
        if nNan > 4000
            error('too many missing values!');
        end
        
        % Get raw data for this trial (rows = channels, cols = time)
        rawData = raw.data;
        
        if nNan > 0
            % --- 1) Interpolate photodiode (1-D vector) ---
            % Fill internal NaNs by linear interpolation
            photoDiode = fillmissing(photoDiode, 'linear');
            % Any leading/trailing NaNs get replaced by nearest neighbor
            photoDiode = fillmissing(photoDiode, 'nearest');
        
            % --- 2) Interpolate rawData along time (dimension 2) ---
            % Linear interpolation across time for each channel
            rawData = fillmissing(rawData, 'linear', 2);
            % Nearest neighbor for any remaining edge NaNs
            rawData = fillmissing(rawData, 'nearest', 2);
        end

        raw.data = rawData;





    photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

    downs = find(photoDiode(1:end-1) > P.pd.zthresh & ...
                photoDiode(2:end)    < P.pd.zthresh);
    ups = find(photoDiode(1:end-1) < P.pd.zthresh & ...
               photoDiode(2:end)   > P.pd.zthresh);
    difVals = ups - downs; %difVals is length of TTL pulses
    downs(difVals>P.pd.maxPulseSamp) = []; 
    ups(difVals>P.pd.maxPulseSamp) = []; 
    difVals(difVals>P.pd.maxPulseSamp) = []; 

    starti = 1; 
    for di = 5:length(downs)
        if downs(di) - downs(di-4) < 3500
            starti = di; 
            if strcmp(raw.sessID, '250929_Dupi_NMH_GH_2') || ...
                strcmp(raw.sessID, '260105_OBE_NWU_ZF_1') || ... 
                strcmp(raw.sessID, '251120_Dupi_NMH_JL_2')
                break
            end
        end
    end
    downs(1:starti) = []; 
    difVals(1:starti) = []; 
    downs(difVals<P.pd.minPulseSamp) = []; 
    difVals(difVals<P.pd.minPulseSamp) = []; 

    trialMarks = downs(difVals < P.pd.trialSplitSamp);
    sniffMarks = downs(difVals > P.pd.trialSplitSamp); 

    trialMarks(diff(trialMarks)<raw.fs_raw) = []; 

    if ~isempty(P.ttl.removeTrialMarksIdx)
        trialMarks(P.ttl.removeTrialMarksIdx) = []; %aberant extra TTL  
    end
    
    
 

    if length(trialMarks) ~= P.ttl.expectedTrialCount
        error('wrong trial count!')
        
    end

    confirmMarks = []; 
    for tt = 1:length(trialMarks)
        if tt < length(trialMarks)
            idx = find(sniffMarks>trialMarks(tt) & ...
                sniffMarks<trialMarks(tt+1));
        else
            idx = find(sniffMarks>trialMarks(tt));
        end
        if length(idx) == 2
            if sniffMarks(idx(2))- sniffMarks(idx(1)) < raw.fs_raw
                confirmMarks = [confirmMarks sniffMarks(idx(1))];
                sniffMarks(idx) = []; 
            end
        end
    end
    
    % confirmMarks = sniffMarks(diff(sniffMarks)<raw.fs_raw);
    % idx = find(diff(sniffMarks)<raw.fs_raw);
    % sniffMarks([idx, idx+1]) = []; 




    figure  %('visible', false, 'position', [0,0,1000,500])
    plot(photoDiode)
    xline(trialMarks)
    xline(confirmMarks, 'color', 'red')
    xline(sniffMarks, 'color', 'green')
    title([raw.sessID ' TTLs'])



    trialStarts = trialMarks(1:2:30); 
    xline(trialStarts, 'color', 'magenta')
    buttonPresses = trialMarks(2:2:30); 
    saveas(gcf,fullfile(raw.paths.fig, 'TTLs.jpg'));


    
    %store all TTLs into one matrix: 
    %col 1: trialStarts
    %col 2: buttonPress 
    %col 3: confirmatory sniff (trial end)
    %col 4: free sniff 1
    %col 5: free sniff 2 
    %      ....
    %col 20: free sniff 17
    TTLs = nan(15, 20); 
    TTLs(:,1) = trialStarts; 
    TTLs(:,2) = buttonPresses; 
    TTLs(:,3) = confirmMarks;  
    for triali = 1:15
        idx = sniffMarks;
        idx = idx(idx>trialStarts(triali) & ...
                            idx < buttonPresses(triali));
        for sniffi = 1:length(idx)
            TTLs(triali,3+sniffi) = idx(sniffi);
        end
    end
    TTLs =round(TTLs ./ (raw.fs_raw / P.fs_target));

    outTTLs = table; 
    outTTLs.trialStart = TTLs(:,1); 
    outTTLs.buttonPress = TTLs(:,2); 
    outTTLs.confirmSniff = TTLs(:,3); 
    outTTLs.free1 = TTLs(:,4); 
    outTTLs.free2 = TTLs(:,5); 
    outTTLs.free3 = TTLs(:,6); 
    outTTLs.free4 = TTLs(:,7); 
    outTTLs.free5 = TTLs(:,8); 
    outTTLs.free6 = TTLs(:,9); 
    outTTLs.free7 = TTLs(:,10); 
    outTTLs.free8 = TTLs(:,11); 
    outTTLs.free9 = TTLs(:,12); 
    outTTLs.free10= TTLs(:,13); 
    outTTLs.free11= TTLs(:,14); 
    outTTLs.free12= TTLs(:,15); 
    outTTLs.free13= TTLs(:,16); 
    outTTLs.free14= TTLs(:,17); 
    outTTLs.free15= TTLs(:,18); 
    outTTLs.free16= TTLs(:,19); 
    outTTLs.free17= TTLs(:,20); 


end