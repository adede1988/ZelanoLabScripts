function [outSniffs] = detect_sniffs_from_TTLs(R, P, outDat)

    %outDat.TTL      table with col1: trialStart
    %                           col2: buttonpress
    %                           col3-20: sniffs

    outSniffs = zeros(30, 7); 
    TTL = outDat.TTL; 
    %col 1: sniff onset index into data
    %col 2: trial number
    %col 3: sniff within trial number
    %col 4: off from TTL by
    %col 6: sniff type 1 = start, 2 = free, 3 = confirm
    %col 7: adjustment for phase align
    oi = 1; %index variable for out sniffs

    for triali = 1:size(TTL,1)
        targets = TTL(triali, 3:end);
        targets = table2array(targets); 
        targets(isnan(targets)) = [];
        targets = sort(targets); 
        sniffi = 1; 
        for tt = 1:length(targets)
            if tt == 1 || tt == length(targets) 
                %first and last sniff! should be after TTL
                startSearch = targets(tt) - 30 - P.cuedBackBuff; 
                endSearch = targets(tt)+2000 - 30; %search window of 4s
                val = find(R.testSig(startSearch:endSearch-1)<...
                                                    P.respThresh  & ... 
                    R.testSig(startSearch+1:endSearch) > P.respThresh, 1);
                if tt == 1
                    type = 1; 
                else
                    type = 3; 
                end
             
            else
                %free sniffs! should be before TTL
                startSearch = targets(tt)- 30 - 1000; %search window of 2s
                endSearch = targets(tt) -30 + 300; 

                val = find(R.testSig(startSearch:endSearch-1)<...
                                                    P.respThresh & ... 
                    R.testSig(startSearch+1:endSearch) > P.respThresh);
                % val(val<16) = 16; 
                % val(val>1285) = 1285; 
                if length(val)>1
                    testSpace = R.testSig(startSearch:endSearch-1);
                    L = length(testSpace); 
                    stopAt = arrayfun(@(x) min([x+15, L]),... 
                                        val); 
                    startAt= arrayfun(@(x) max([x-15, 1]), val); 
                    sizes = testSpace(stopAt) - testSpace(startAt); 
                    [~, maxSniff] = max(sizes); 
                    val = val(maxSniff); 
                end
                type = 2; 
            end
            
                
                %was a sniff found? if yes then find overall idx
                %val is how far from the start search the sniff happened
            if ~isempty(val)
                targidx = val + startSearch; 
              
                %col 1: sniff onset index
                outSniffs(oi,1) = targidx; 
                %col 2: trial number
                outSniffs(oi,2) = triali; 
                %col 3: sniff within trial number
                outSniffs(oi,3) = sniffi; 
                sniffi = sniffi + 1; 
                %col 4: off from TTL by
                if type == 2
                    outSniffs(oi,4) = - 30 - 1000 +val; 
                else
                    outSniffs(oi,4) = val - 30 - P.cuedBackBuff;
                end
                %col 6: sniff type 1 = start, 2 = free, 3 = confirm
                outSniffs(oi,6) = type;
                oi = oi + 1; 
            end
        end
    end


    TTL = table2array(TTL); 


    for triali  = 1:size(TTL,1)
        figure('visible', false)
        idx = TTL(triali, 1); 
        endidx = TTL(triali,3)+P.fs_target*4; 
        rspTrial = R.smoothR(idx:endidx); 
        testHere = R.testSig(idx+30:endidx+30);
     
        plot(rspTrial )
       
        idx = TTL(triali,3:end) - TTL(triali, 1);
        idx(isnan(idx)) = []; 
        xline(idx, 'linewidth', 2)
        title([outDat.sessID ...
            ' backBuff: ' num2str(P.cuedBackBuff) ...
            ' rspThresh:' num2str(P.respThresh) ...
            ' trial:' num2str(triali)], 'interpreter', 'none')
        curSniffs = outSniffs(outSniffs(:,2) == triali, 1); 
        if ~isempty(curSniffs)
            xline(curSniffs - TTL(triali, 1), ...
                'color', 'green', 'linewidth', 2); 
        end
        yyaxis right
        plot(testHere)
        saveas(gcf,fullfile(outDat.figs, ...
                    ['sniffsTrial' num2str(triali) '.jpg']));
    end




    





end