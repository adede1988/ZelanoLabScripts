function [outDat, P] = paramCheck(outDat, P)

    set(0, 'defaultfigurewindowstyle', 'docked')

% is the respiration index correct? 
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 

    
    figure; 
    plot(rspDat', 'color', 'k')
    hold on 
    plot(rspDat(P.rspIDX,:).*P.rspFlip, 'color', 'red')
    xlim([10000 100000])

    userData = input('Press 1 to accept; 0 to reject: ');
    ii = 1; 
    while userData == 0
        figure; 
        plot(rspDat(ii,:), 'color', 'k')
        xlim([10000 100000])
        userData = input('Press 1 to accept; 0 to reject: ');
        if userData == 1
            P.rspIDX = ii; 
            userData = input('should it be flipped? 0 = no; 1=yes');
            if userData == 0
                P.rspFlip = 1; 
            else 
                P.rspFlip = -1; 
            end
            userData = 1; 
        end
        ii = ii+1;
        if ii == 4
            ii = 1; 
        end
    end


% should we do spike cleaning or reject any macros? 
    figure;
    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    if sum(idx)>0
        macDat = outDat.data(idx,:); 
        
        hold on 
        for ii = 1:size(macDat,1)
            plot(macDat(ii,:)+ii*100)
        end
        legend()
        xlim([10000 30000])
        userData = input('remove any macros? enter as [# #] to indicate channels or [] to indicate none');
        
        P.macroRemove = userData;
        
        userData = input('spike removal, should be avoided if macro channels are being removed? 0=no; 1=yes');

        P.spikeClean = userData; 
    end









    set(0, 'defaultfigurewindowstyle', 'normal')















end