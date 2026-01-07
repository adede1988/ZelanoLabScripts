function ax = plotICATopo(ica, icIdx, theta, phi)
% plotICATopo  Quick scalp topo for one ICA component.
%
% ax = plotICATopo(ica, icIdx, theta, phi)
%   ica    : struct with fields A (nChan x nIC) or W (nIC x nChan)
%   icIdx  : component index (1-based)
%   theta  : per-channel polar angle from +Z (colatitude). rad or deg.
%   phi    : per-channel azimuth around Z. rad or deg.
%
% Notes:
% - Uses A(:,icIdx) if available; else uses pinv(W)(:,icIdx).
% - theta/phi may be in degrees or radians (auto-detected).
% - Assumes channels lie on a unit sphere; projects with:
%       X = cos(phi).*sin(theta);  Y = sin(phi).*sin(theta)

    % --- get spatial weights for the component
    if isfield(ica,'A') && ~isempty(ica.A)
        w = ica.A(:, icIdx);
    elseif isfield(ica,'W') && ~isempty(ica.W)
        A = pinv(ica.W);
        w = A(:, icIdx);
    else
        error('plotICATopo:NoAorW', 'Need ica.A or ica.W in the struct.');
    end
    w = w(:);

    % --- normalize angles to radians
    if max(abs(theta)) > 2*pi || max(abs(phi)) > 2*pi
        theta = deg2rad(theta);
        phi   = deg2rad(phi);
    end
    theta = theta(:); phi = phi(:);

    % --- 3D->2D projection on unit disk
    X = cos(phi).*sin(theta);
    Y = sin(phi).*sin(theta);

    % --- grid + interpolation on the disk
    n = 200;
    [xg, yg] = meshgrid(linspace(-1,1,n), linspace(-1,1,n));
    mask = (xg.^2 + yg.^2) <= 1.0;

    F = scatteredInterpolant(X, Y, w, 'natural', 'nearest');
    zg = nan(size(xg));
    zg(mask) = F(xg(mask), yg(mask));

    % --- plot
    figure('Color','w');
    ax = axes; hold(ax,'on'); axis(ax,'equal','tight'); axis(ax,[-1 1 -1 1]);
    himg = imagesc(ax, linspace(-1,1,n), linspace(-1,1,n), zg);
    set(himg,'AlphaData', ~isnan(zg)); set(ax,'YDir','normal');

    % scalp outline, nose, ears
    th = linspace(0,2*pi,360);
    plot(ax, cos(th), sin(th), 'k', 'LineWidth', 1.2);
    plot(ax, [0 .06 -0.06 0], [1 1.10 1.10 1], 'k', 'LineWidth', 1);    % nose
    plot(ax, [1 1.10 1 1], [.1 0 -.1 .1], 'k');                         % right ear
    plot(ax, [-1 -1.10 -1 -1], [.1 0 -.1 .1], 'k');                     % left ear

    % electrodes
    plot(ax, X, Y, 'k.', 'MarkerSize', 12);

    title(ax, sprintf('ICA %d scalp map (A(:,%d))', icIdx, icIdx));
    colorbar(ax); colormap(ax, parula);
    axis(ax,'off');
end
