function outDat = flagBadBreaths(outDat)
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    
    %choose which one looks right
    % rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    
    idx = outDat.rspIDX; 
    rspDat = squeeze(rspDat(idx, :));  
    
    %flip signal
    rspDat = rspDat .* outDat.rspFlip;

    %get the heartbeat interval time series
    idx = cellfun(@(x) contains(x, 'RRint'), outDat.labels);
    rrDat = outDat.data(idx, :); 


    %new values to add: 
    %Good Breaths indicator 1 = good 0 = bad
    %RRmin
    %RRmax
    %RRmax - RRmin
    for bb = 1:size(outDat.behDat, 1)
       
        idx = outDat.behDat.finalOnset(bb);
        bL = outDat.behDat.length(bb) ; 
        %check breath isn't too near start or end of recording
        %also not more than 18 seconds long
        endPnt = idx + outDat.fs*20  - outDat.fs*2;
        startPnt = endPnt - outDat.fs*20 + 1;  
        if startPnt > 0 && ...
           endPnt < length(outDat.data(1,:)) && ...
           bL < 18

          
           curRsp = rspDat(startPnt:endPnt); 
           curRR = rrDat(startPnt:endPnt); 
           %quality check: 
           %end point should be higher than minimum
           %start point should be higher than minimum
           %one max
           %one min
           smoothRsp = smoothdata(curRsp, 'gaussian', outDat.fs/2);
           bonset = round(outDat.fs*2); 
           boffset = round(outDat.fs*2+bL*outDat.fs); 
           startVal = smoothRsp(bonset); 
           endVal = smoothRsp(boffset);
            
           lowThresh = prctile(smoothRsp(bonset:boffset), 5); 
           highThresh = prctile(smoothRsp(bonset:boffset), 95); 
           
           upPeakIdx = find(smoothRsp(bonset:boffset)<highThresh & ...
                 smoothRsp(bonset+1:boffset+1)>highThresh) ;
            downPeakIdx = find(smoothRsp(bonset:boffset)>highThresh & ...
                               smoothRsp(bonset+1:boffset+1)<highThresh) ;
            
            upTroughIdx = find(smoothRsp(bonset:boffset)<lowThresh & ...
                               smoothRsp(bonset+1:boffset+1)>lowThresh) ;
            downTroughIdx = find(smoothRsp(bonset:boffset)>lowThresh & ...
                                 smoothRsp(bonset+1:boffset+1)<lowThresh) ;
            
            % --- Helper function: true if any pair is 1s apart ---
            tooFar = @(idx) numel(idx) > 1 && any(diff(idx) > outDat.fs);
            
            % --- Reject condition ---
            if startVal < highThresh && ...
               endVal   < highThresh && ...  
               startVal > lowThresh  && ...
               endVal   > lowThresh  && ...
               ~tooFar(upPeakIdx)    && ...
               ~tooFar(downPeakIdx)  && ...
               ~tooFar(upTroughIdx)  && ...
               ~tooFar(downTroughIdx) || ...
                strcmp(outDat.behDat.task(bb), 'cyclicSigh') %just take all cyclic sighs, might have to do similar for slow breathing

               %good breath! 
               outDat.behDat.goodBreath(bb) = 1; 

               %get RR variability: 
                %col15: RRmin
                %col16: RRmax
                %col17: RRmax - RRmin
               outDat.behDat.maxRR(bb) = max(curRR( ...
                                            bonset:boffset));
               outDat.behDat.minRR(bb) = min(curRR( ...
                                            bonset:boffset));
               outDat.behDat.RR_max_min(bb) = outDat.behDat.maxRR(bb) -...
                                            outDat.behDat.minRR(bb);

           else
                outDat.behDat.goodBreath(bb) = 0; 

               %get RR variability: 
                %col15: RRmin
                %col16: RRmax
                %col17: RRmax - RRmin
               outDat.behDat.maxRR(bb) = max(rrDat( ...
                                            bonset:boffset));
               outDat.behDat.minRR(bb) = min(rrDat( ...
                                            bonset:boffset));
               outDat.behDat.RR_max_min(bb) = outDat.behDat.maxRR(bb) -...
                                            outDat.behDat.minRR(bb);
                
           end
        end

    end









end