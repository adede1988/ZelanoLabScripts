function [fridge, ridgeIdx, ridgePow] = ridge_track(E, freqs, penalty, numRidges, BW)
% RIDGE_TRACK  Forward-backward penalized frequency-ridge tracker.
% Faithful MATLAB port of Frequency_ridge_tracking/tfridge_tracking/
% ridge_tracking.py (extract_fridges, D. Bondesson) — the repo the project
% specifies. Operates on a real, nonnegative time-frequency ENERGY/POWER
% matrix (we pass the superlet power, optionally z-scored+clipped).
%
%   [fridge, ridgeIdx, ridgePow] = ridge_track(E, freqs, penalty, numRidges, BW)
%     E        : nFreq x nTime nonnegative energy/power
%     freqs    : nFreq x 1 frequency scale (Hz), used for the jump penalty
%     penalty  : scalar jump penalty (freq^2 distance is multiplied by this)
%     numRidges: number of ridges to extract (default 1)
%     BW       : # freq bins zeroed around a found ridge before next (default 2)
%   fridge   : nTime x numRidges ridge frequency (Hz)
%   ridgeIdx : nTime x numRidges ridge freq-bin index
%   ridgePow : nTime x numRidges energy/power ON the ridge (from ORIGINAL E)

    if nargin<4 || isempty(numRidges), numRidges=1; end
    if nargin<5 || isempty(BW), BW=2; end
    [nF,nT] = size(E);
    fs = freqs(:);
    penMat = (fs - fs').^2 * penalty;      % nF x nF, penMat(a,b)=(f_a-f_b)^2*penalty
    ep = eps;

    fridge=zeros(nT,numRidges); ridgeIdx=zeros(nT,numRidges); ridgePow=zeros(nT,numRidges);
    Ework = E;                              % gets peeled between ridges
    for r = 1:numRidges
        emax = max(Ework,[],1); emax(emax<=0)=ep;
        En = -log(Ework./emax + ep);        % neg-log normalized cost (low = high energy)

        % ---- forward accumulation ----
        pe = En;
        for t = 2:nT
            prev = pe(:,t-1).';             % 1 x nF
            % cost(a) = min_b( pe(b,t-1) + penMat(a,b) )
            [mincost] = min(penMat + prev, [], 2);   % nF x 1
            pe(:,t) = pe(:,t) + mincost;
        end
        [~, ridge] = min(pe, [], 1);        % 1 x nT forward ridge

        % ---- backward refinement ----
        for t = nT-1:-1:1
            nx  = ridge(t+1);
            val = pe(nx,t+1) - En(nx,t+1);
            cand = pe(:,t) + penMat(nx,:).';           % pe(b,t)+penalty(nx,b)
            [~,bi] = min(abs(val - cand));
            ridge(t) = bi;
        end

        ind = sub2ind([nF,nT], ridge, 1:nT);
        fridge(:,r)   = fs(ridge);
        ridgeIdx(:,r) = ridge(:);
        ridgePow(:,r) = E(ind).';           % power from ORIGINAL (unpeeled) E

        % ---- peel BW bins around ridge for next ridge ----
        for t = 1:nT
            lo = max(1, ridge(t)-BW); hi = min(nF, ridge(t)+BW);
            Ework(lo:hi,t) = 0;
        end
    end
end
