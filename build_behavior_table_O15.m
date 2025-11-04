function behDat = build_behavior_table_O15(sniffs, rawBeh)

    %integrate information from the behavioral data with the sniff onset
    %information
   
    
    sniffTypes = {"start", "free", "confirm"};     
    
    behDat = table; 
    behDat.sniffOnset = sniffs(:,1); %sniff onset
    behDat.n = sniffs(:,2); %trial num
    behDat.wiTriali = sniffs(:,3); %sniff within the trial 
    behDat.TTLoffSet = sniffs(:,4); %how far off from the TTL
    behDat.sniffType = sniffs(:,6); %what kind of sniff is it? 
    behDat.sniffLabel = [sniffTypes{sniffs(:,6)}]'; %readable sniff label
    for ii = 1:length(rawBeh.target)
        idx = find(behDat.n == ii);
        for jj = 1:length(idx)
            behDat.target(idx(jj)) = string(rawBeh.target{ii});
            behDat.response(idx(jj)) = string(rawBeh.response{ii});
            behDat.expScore(idx(jj)) = rawBeh.expScore(ii);
        end
    
    end














end