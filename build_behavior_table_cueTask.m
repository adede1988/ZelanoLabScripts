function behDat = build_behavior_table_cueTask(sniffs, rawBeh)

    %integrate information from the behavioral data with the sniff onset
    %information
   
    
    sniffTypes = {"cued"};     
    
    behDat = table; 
    behDat.sniffOnset = sniffs(:,1); %sniff onset
    behDat.n = sniffs(:,2); %trial num
    behDat.wiTriali = sniffs(:,3); %sniff within the trial 
    behDat.TTLoffSet = sniffs(:,4); %how far off from the TTL
    behDat.sniffType = sniffs(:,6); %what kind of sniff is it? 
    behDat.sniffLabel = [sniffTypes{sniffs(:,6)}]'; %readable sniff label
    for ii = 1:length(rawBeh.n)
        idx = find(behDat.n == ii);
        for jj = 1:length(idx)
            behDat.cue(idx(jj)) = rawBeh.cue(ii);
            behDat.odor(idx(jj)) = rawBeh.odor(ii);
            behDat.response(idx(jj)) = rawBeh.response(ii);
            if ~isempty(rawBeh.response_str{ii})
                behDat.respString(idx(jj)) = string(rawBeh.response_str(ii));
            else
                behDat.respString(idx(jj)) = "SKIP"; 
            end
            behDat.type(idx(jj)) = string(rawBeh.type(ii));
        end
    
    end














end