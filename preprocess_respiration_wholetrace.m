function [R] = preprocess_respiration_wholetrace(outDat)


    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    rspDat = squeeze(rspDat(outDat.rspIDX, :)); 
    %flip signal
    rspDat = rspDat .* outDat.rspFlip;
    
    smoothRsp = smoothdata(rspDat, 'gaussian', round(outDat.fs*.6)); 
    win = round(outDat.fs * .06); 
    deflectionDetect = (smoothRsp(win:end) - smoothRsp(1:end-win+1)).^2 .* ...
                                    smoothRsp(win:end); 
    deflectionDetect(deflectionDetect<0) = 0; 
    deflectionDetect(deflectionDetect>10000) = 10000; 
    deflectionDetect = smoothdata(deflectionDetect, 'gaussian', 500);
   

    analytic = hilbert(rspDat); 
    analytic = smoothdata(analytic, 1, 'gaussian', 500);
    analytic = lowpass(analytic, 1, outDat.fs); 
    rspPhase = angle(analytic);

    R = struct; 
    R.rspDat = rspDat; 
    R.smoothR = smoothRsp; 
    R.testSig = deflectionDetect; 
    R.analytic = analytic; 
    R.rspPhase = rspPhase; 
    









end