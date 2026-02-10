function [idx_start_near, idx_seq] = match_phase_indices(phases, idx0, ...
    startPhase, targPhases, rsp)
% phases:   [T x N] matrix of phases in radians (typically in [-pi, pi])
% idx0:     scalar starting time index (1-based)
% startPhase: scalar phase (rad)
% targPhases: [M x 1] (or 1xM) vector of target phases
%
% Returns:
%   idx_start_near : [1 x N] index per column for the closest phase to startPhase
%   idx_seq        : [M x N] forward-only indices per column for each targPhase
%                    Row 1 repeats idx_start_near (the match for startPhase).
%                    Unmatched (ran out of time) are NaN.

if size(phases, 1) < size(phases, 2)
    phases = phases';
    rsp = rsp';
end

[T, N] = size(phases); %Time X Trials 
targPhases = targPhases(:);          % ensure column
M = numel(targPhases);

idx_start_near = nan(1, N);
idx_seq        = nan(M, N); %num target phases X trials

for c = 1:N
    ph = phases(:, c);

    figure; plot(rsp(:,c))
    hold on 
    plot(ph)
    % --- Step 1: nearest (in time) best match to startPhase by searching
    %             both backward [1:idx0] and forward [idx0:T] ---
    dAll = abs(angle(exp(1i*(ph - startPhase))));  % circular absolute difference

    %forward search
    curDif = dAll(idx0); 
    nexDif = dAll(idx0+1); 
    iFwd = idx0; 
    while curDif > nexDif
        iFwd = iFwd + 1; 
        curDif = dAll(iFwd); 
        nexDif = dAll(iFwd+1); 
    end
    minFwd = curDif;
    %back search
    curDif = dAll(idx0); 
    nexDif = dAll(idx0-1); 
    iBack = idx0; 
    while curDif > nexDif
        iBack = iBack - 1; 
        curDif = dAll(iBack); 
        nexDif = dAll(iBack-1); 
    end
    minBack = curDif; 

    if minBack < minFwd
        i0 = iBack;
    elseif minFwd < minBack
        i0 = iFwd;
    else
        % tie: pick the one closer (in time) to idx0
        if abs(iBack - idx0) <= abs(iFwd - idx0)
            i0 = iBack;
        else
            i0 = iFwd;
        end
    end

    idx_start_near(1, c) = i0;
    idx_seq(1, c)        = i0;

    % --- Step 2: forward-only matches for subsequent target phases ---
    prevIdx = i0;
    for k = 2:M
        if prevIdx >= T
            idx_seq(k, c) = NaN;     % no room left in time
            continue
        end
        seg = ph(prevIdx+1:T);
        idx_k = find(seg>=targPhases(k), 1);
        if isempty(idx_k)
            idx_seq(k, c) = NaN;     % no room left in time
            continue
        else
            idx_k       = prevIdx + idx_k;
          
            idx_seq(k,c)= idx_k;
            prevIdx     = idx_k;
        end
    end
end
end
