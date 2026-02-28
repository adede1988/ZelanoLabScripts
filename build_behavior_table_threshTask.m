function behDat = build_behavior_table_threshTask(sniffs, rawBeh)

    %integrate information from the behavioral data with the sniff onset
    %information
   
    
    sniffTypes = {"cued"};     
    trialTypes = {'air', 'low', 'med'}; 
    behDat = table; 
    behDat.sniffOnset = sniffs(:,1); %sniff onset
    behDat.n = sniffs(:,2); %trial num
    behDat.wiTriali = sniffs(:,3); %sniff within the trial 
    behDat.TTLoffSet = sniffs(:,4); %how far off from the TTL
    behDat.sniffType = sniffs(:,6); %what kind of sniff is it? 
    behDat.sniffLabel = [sniffTypes{sniffs(:,6)}]'; %readable sniff label
    for ii = 1:length(rawBeh.trialNum)
        
        idx = find(behDat.n == ii);
        for jj = 1:length(idx)
            
            behDat.odor(idx(jj)) = rawBeh.Odor(ii);
            behDat.pleasantness(idx(jj)) = rawBeh.pleasantness(ii);
            behDat.intensity(idx(jj)) = rawBeh.intensity(ii);
           
            behDat.type(idx(jj)) = string(trialTypes{rawBeh.Odor(ii)});
        end
    
    end














end