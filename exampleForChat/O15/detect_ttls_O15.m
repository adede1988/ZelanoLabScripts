function outTTLs = detect_ttls_O15(raw, P)

    idx = cellfun(@(x) contains(x, 'event'), raw.labels);
    photoDiode = raw.data(idx, :); 
    
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
    
    confirmMarks = sniffMarks(diff(sniffMarks)<raw.fs_raw);
    idx = find(diff(sniffMarks)<raw.fs_raw);
    sniffMarks([idx, idx+1]) = []; 

    %THIS WON'T BE CORRECT FOR THS SHORT EXAMPLE DATA SO REMOVE: 
    % if length(trialMarks) ~= P.ttl.expectedTrialCount
    %     error('wrong trial count!')
    % 
    % end

    figure('visible', false, 'position', [0,0,1000,500])
    plot(photoDiode)
    xline(trialMarks)
    xline(confirmMarks, 'color', 'red')
    xline(sniffMarks, 'color', 'green')
    title([raw.sessID ' TTLs'])



    trialStarts = trialMarks(1:2:end); 
    xline(trialStarts, 'color', 'magenta')
    buttonPresses = trialMarks(2:2:end); 
    saveas(gcf,fullfile(raw.paths.fig, 'TTLs.jpg'));


    
    %store all TTLs into one matrix: 
    %col 1: trialStarts
    %col 2: buttonPress 
    %col 3: confirmatory sniff (trial end)
    %col 4: free sniff 1
    %col 5: free sniff 2 
    %      ....
    %col 20: free sniff 17
    TTLs = nan(length(trialStarts), 20); 
    TTLs(:,1) = trialStarts; 
    TTLs(:,2) = buttonPresses; 
    TTLs(:,3) = confirmMarks;  
    for triali = 1:length(trialStarts)
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