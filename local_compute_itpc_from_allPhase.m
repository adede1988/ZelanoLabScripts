% ============================================================
% Local function (keep at end of script, or convert to separate .m)
% ============================================================
function itpc = local_compute_itpc_from_allPhase(allPhase, bVec, idxMat, nBreaths, nT, nF)
    % bVec: [nE x 1] breath index for each event
    % idxMat: [nE x nOff] time indices for each event window
    nE   = numel(bVec);
    nOff = size(idxMat,2);

    itpc = nan(nOff, nF);

    % linear indices into [nBreaths x nT] for each event x offset
    % sub2ind supports matrices; returns [nE x nOff]
    lin = sub2ind([nBreaths, nT], repmat(bVec(:),1,nOff), idxMat);

    for f = 1:nF
        ph2 = allPhase(:,:,f);          % [nBreaths x nT]
        ph  = ph2(lin);                 % [nE x nOff]
        itpc(:,f) = abs(mean(exp(1j*ph), 1, 'omitnan')).';  % [nOff x 1]
    end
end
