function [Zi, Xi, Yi, hSurf] = miniTopo(m, x, y, plotit)
% miniTopo
% Minimal EEGLAB-style topomap:
%   - auto-squeezes input coords so the outermost sensor sits at r=0.5
%   - interpolates on a cartesian grid
%   - masks outside head circle
%
% Inputs:
%   m      [nChan x 1] intensity
%   x,y    [nChan x 1] coordinates (any common scale; will be auto-squeezed)
%   plotit (0/1) optional (default 1)
%
% Outputs:
%   Zi, Xi, Yi  interpolated grid + mesh
%   hSurf       surface handle (if plotted)

if nargin < 4 || isempty(plotit), plotit = 1; end

% ---- constants (EEGLAB-ish) ----
rmax      = 0.5;          % anatomically correct head radius
GRID_SCALE= 67;
AXHEADFAC = 1.05;
CIRCGRID  = 201;
HEADCOLOR = [0 0 0];
HLINEWIDTH= 1.7;
RINGW     = 0.03;         % blanking ring width (cosmetic)

% ---- columnize + drop bad chans ----
m = double(m(:));
x = double(x(:));
y = double(y(:));
good = isfinite(m) & isfinite(x) & isfinite(y);
m = m(good); x = x(good); y = y(good);

% ---- squeeze coords so data fills head nicely ----
rd = hypot(x, y);
squeezefac = rmax / (max(rd)*1.02);
x = x * squeezefac;
y = y * squeezefac;

% ---- grid ----
xi = linspace(-rmax, rmax, GRID_SCALE);
yi = linspace(-rmax, rmax, GRID_SCALE);
[Xi, Yi] = meshgrid(xi, yi);

% ---- interpolate (smooth fill) ----
Zi = griddata(x, y, m, Xi, Yi, 'v4');

% clamp to data range to prevent v4 overshoot halos
mn = min(m); mx = max(m);
Zi = min(max(Zi, mn), mx);

% ---- mask outside head ----
mask = sqrt(Xi.^2 + Yi.^2) <= rmax;
Zi(~mask) = NaN;

% ---- plot ----
hSurf = [];
if plotit > 0
    ax = gca;
    hold(ax,'on');

    hSurf = surface(ax, Xi, Yi, zeros(size(Zi)), Zi, ...
        'EdgeColor','none', 'FaceColor','interp');

    axis(ax,'equal'); axis(ax,'off');
    set(ax,'XLim',[-rmax rmax]*AXHEADFAC, 'YLim',[-rmax rmax]*AXHEADFAC);

    % stable color scaling
    % caxis(ax, [mn mx]);
  

    % blanking ring to hide grid edge jaggies
    circ = linspace(0, 2*pi, CIRCGRID);
    rx = cos(circ); ry = sin(circ);
    bg = get(gcf,'color');
    rout = rmax -.05 + RINGW;
    ringx = [rx*rout, fliplr(rx*rmax)];
    ringy = [ry*rout, fliplr(ry*rmax)];
    patch(ax, ringx, ringy, bg, 'EdgeColor','none', 'HitTest','off');
    plot(ax, rx*rout, ry*rout, 'color', 'k', 'linewidth', 3);

   % ---- add nose + ears (drop into miniTopo after head outline is drawn) ----
    % assumes you have: ax, rmax (head radius), HEADCOLOR, HLINEWIDTH already defined
    
    ELECTRODE_HEIGHT = 2;  % just to lift lines above the surface
    
    % Nose (EEGLAB-style triangle)
    base  = rmax - 0.0046;
    basex = 0.18*rmax;
    tip   = 1.15*rmax;
    tiphw = 0.04*rmax;
    tipr  = 0.01*rmax;
    
    plot3(ax, [basex; tiphw; 0; -tiphw; -basex], ...
             [base;  tip-tipr; tip; tip-tipr; base] - .02, ...
             ELECTRODE_HEIGHT*ones(5,1), ...
             'Color', HEADCOLOR, 'LineWidth', HLINEWIDTH, 'HitTest','off');
    ylim([-rmax tip])
    % Ears (EEGLAB ear template, scaled to rmax)
    q = 0.04;
    EarX = [.497-.005 .510 .518 .5299 .5419 .54 .547 .532 .510 .489-.005]; % template for rmax=0.5
    EarY = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    
    sf = rmax/0.5;   % ear template assumes rmax=0.5
    earOff = rmax;   % EEGLAB ear horizontal offset (in "rmax=0.5 units")
    plot3(ax,  (EarX*sf - rmax + min(EarX*sf)),  (EarY*sf), ELECTRODE_HEIGHT*ones(size(EarY)),...
          'Color', HEADCOLOR, 'LineWidth', HLINEWIDTH);
    plot3(ax, -(flip(EarX*sf)) + rmax +max(-(flip(EarX*sf))),  flip(EarY*sf), ELECTRODE_HEIGHT*ones(size(EarY)), ...
          'Color', HEADCOLOR, 'LineWidth', HLINEWIDTH, 'HitTest','off');
    xlim([min(-(flip(EarX*sf))) max((EarX*sf - rmax + min(EarX*sf)))])
        hold(ax,'off');
end
end