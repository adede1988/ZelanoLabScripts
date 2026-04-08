function outDat = process_respiration_breathing(outDat, P)

    % respiration  
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx,:); 
    rspDat = rspDat(P.rspIDX,:);
    rspDat = rspDat .* P.rspFlip;
    bmObj = breathTemplates4(rspDat, outDat.fs);
    %col 1: onset Y value
    %col 2: onset tim
    %col 3: peak Y value
    %col 4: peak tim
    %col 5: end Y value
    %col 6: end tim
    %col 7: length (end tim - onset tim)
    %col 8: amp (peak Y - avg of two ends)
    %col 9: idx of peak in rspSig2
    %col10: exhale peak Y value
    %col11: exhale peak tim
    %col12: condition
    %col13: empty
    %col14: index
    fs = outDat.fs; 
    if ismember('cndName', outDat.behDat.Properties.VariableNames)
        orderIdx = arrayfun(@(x) find(outDat.behDat.order == x, 1),...
                            unique(outDat.behDat.order));
        orderIdx(1) = []; %get rid of the pre ratings
        
        if length(orderIdx) == length(outDat.TTL) %sanity check! 
            for ii = 1:length(orderIdx)
                if ismember('cyclicSigh', ...
                        outDat.behDat.cndName(orderIdx(ii)))
                    if ii == length(orderIdx)
                        endIdx = length(rspDat); 
                    else
                        endIdx = outDat.TTL(outDat.behDat.order(orderIdx(ii+1))); 
                    end
                    startIdx = outDat.TTL(outDat.behDat.order(orderIdx(ii))); 

                    cndSeg = rspDat(startIdx:endIdx-1); 
                    timSeg = 1/fs:1/fs:(endIdx-startIdx)/fs;
                    figure; 
                    hold on 
                    plot(timSeg, cndSeg)
                    starts = bmObj(:,2); 
                    bmIdx = find(starts>startIdx/fs & starts < endIdx/fs); 
                    starts = starts(starts>startIdx/fs & starts < endIdx/fs); 
                    starts = starts - startIdx/fs; %adjust to 0 start
                    deleteIdx = []; 
                    for bb = 1:length(starts)-1
                    %correct cyclic sighs are when participant has two
                    %consecutive inhales near each other and no more for
                    %several seconds
                        if starts(bb+1) - starts(bb) > 3
                            deleteIdx = [deleteIdx bmIdx(bb)]; 
                            % scatter(starts(bb), 0, 30, 'red', 'filled')
                        else
                            % scatter(starts(bb), 0, 30, 'green', 'filled')
                            
                            %col 1: onset Y value
                            %col 2: onset tim
                            %col 3: peak Y value
                            [bmObj(bmIdx(bb),3), ti] = max([bmObj(bmIdx(bb),3)...
                                                      bmObj(bmIdx(bb+1),3)]);
                            %col 4: peak tim
                            tmp = [bmObj(bmIdx(bb),4) bmObj(bmIdx(bb+1),4)];
                            bmObj(bmIdx(bb),4) = tmp(ti);
                            %col 5: end Y value
                            bmObj(bmIdx(bb),5) = bmObj(bmIdx(bb+1),5); 
                            %col 6: end tim
                            bmObj(bmIdx(bb),6) = bmObj(bmIdx(bb+1),6); 
                            %col 7: length (end tim - onset tim)
                            bmObj(bmIdx(bb),7) = bmObj(bmIdx(bb),6) - ...
                                                 bmObj(bmIdx(bb),2); 
                            %col 8: amp (peak Y - avg of two ends)
                            bmObj(bmIdx(bb),8) = bmObj(bmIdx(bb), 3) - ...
                                          mean([bmObj(bmIdx(bb),1) ...
                                                bmObj(bmIdx(bb),5)]); 
                                                 
                            %col 9: idx of peak in rspSig2
                            tmp = [bmObj(bmIdx(bb),9) bmObj(bmIdx(bb+1),9)];
                            bmObj(bmIdx(bb),9) = tmp(ti);
                            %col10: exhale peak Y value
                            tmp = [bmObj(bmIdx(bb),10) bmObj(bmIdx(bb+1),10)];
                            [bmObj(bmIdx(bb),10), ti] = min(tmp);
                            %col11: exhale peak tim
                             tmp = [bmObj(bmIdx(bb),11) bmObj(bmIdx(bb+1),11)];
                             bmObj(bmIdx(bb),11) = tmp(ti); 
                            %col12: condition
                            %col13: empty
                            %col14: index
                        end
                    end
                    bmObj(deleteIdx,:) = []; 
                else
                     if ii == length(orderIdx)
                        endIdx = length(rspDat); 
                    else
                        endIdx = outDat.TTL(outDat.behDat.order(orderIdx(ii+1))); 
                    end
                    startIdx = outDat.TTL(outDat.behDat.order(orderIdx(ii)));
                    if startIdx == 0 
                        startIdx = 1; 
                    end
                    % figure; 
                    % hold on 
                    % plot(1/500:1/500:(endIdx-startIdx)/500, rspDat(startIdx:endIdx-1))
                    % starts = bmObj(:,2); 
                    % bmIdx = find(starts>startIdx/500 & starts < endIdx/500); 
                    % starts = starts(starts>startIdx/500 & starts < endIdx/500); 
                    % starts = starts - startIdx/500; 
                    % scatter(starts, zeros(size(starts)))
                end
            end
        end
    end
        



    outDat.bmObj = bmObj; 










end