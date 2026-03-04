function [Zi, Xi, Yi] = quickTopo2(m, x, y, plotit)

    headrad = 0.70;
    plotrad = .8;
    intrad  = 1;

    GRID_SCALE = 67;
    AXHEADFAC  = 1.05;
    BLANKINGRINGWIDTH = .035;
    HEADRINGWIDTH     = .007;
    SHADING = 'interp';
    CIRCGRID = 201;
    HEADCOLOR = [0 0 0];
    HLINEWIDTH = 1.7;

    % --- ensure column vectors ---
    m = double(m(:));
    x = double(x(:));
    y = double(y(:));

    % --- channels to interpolate (radial) ---
    intchans = find(isfinite(m) & isfinite(x) & isfinite(y) & (sqrt(x.^2 + y.^2) <= intrad));

    intValues = m(intchans);
    intx = x(intchans);
    inty = y(intchans);

    squeezefac = headrad/plotrad;
    intx = intx * squeezefac;
    inty = inty * squeezefac;

    xmin = min(-headrad, min(intx)); xmax = max(headrad, max(intx));
    ymin = min(-headrad, min(inty)); ymax = max(headrad, max(inty));

    xi = linspace(xmin, xmax, GRID_SCALE);
    yi = linspace(ymin, ymax, GRID_SCALE);
    [Xi, Yi] = meshgrid(xi, yi);

    % --- interpolation: natural neighbor (less ringing than 'v4') ---
    if numel(intValues) >= 3 && range(intValues(isfinite(intValues))) > 0
        F  = scatteredInterpolant(intx, inty, intValues, 'natural', 'none');
        Zi = F(Xi, Yi);
    else
        Zi = nan(size(Xi));
    end

    % --- mask outside the head circle ---
    mask = (sqrt(Xi.^2 + Yi.^2) <= headrad);
    Zi(~mask) = NaN;

    if plotit>0
        set(gca,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC)

        surface(Xi, Yi, zeros(size(Zi)), Zi, 'EdgeColor','none', 'FaceColor','interp');
        axis off
        axis equal

        % lock color scaling to data range (stabilizes "gradient")
        if any(isfinite(intValues))
            caxis([min(intValues(isfinite(intValues))) max(intValues(isfinite(intValues)))]);
        end

        %% Plot filled ring to mask jagged grid boundary
        hwidth = HEADRINGWIDTH;
        hin    = squeezefac*headrad*(1- hwidth/2);

        if strcmp(SHADING,'interp')
            rwidth = BLANKINGRINGWIDTH*1.3;
        else
            rwidth = BLANKINGRINGWIDTH;
        end
        rin = headrad*(1-rwidth/2);
        if hin>rin
            rin = hin;
        end

        circ = linspace(0,2*pi,CIRCGRID);
        rx = sin(circ);
        ry = cos(circ);

        ringx = [[rx(:)' rx(1)]*(rin+rwidth) [rx(:)' rx(1)]*rin];
        ringy = [[ry(:)' ry(1)]*(rin+rwidth) [ry(:)' ry(1)]*rin];
        patch(ringx, ringy, 0.01*ones(size(ringx)), get(gcf,'color'), ...
            'edgecolor','none','hittest','off'); hold on

        %% Plot cartoon head ring
        headx = [[rx(:)' rx(1)]*(hin+hwidth) [rx(:)' rx(1)]*hin];
        heady = [[ry(:)' ry(1)]*(hin+hwidth) [ry(:)' ry(1)]*hin];
        patch(headx, heady, ones(size(headx)), HEADCOLOR, ...
            'edgecolor',HEADCOLOR,'hittest','off'); hold on

        % Plot ears and nose
        base  = headrad-.0046;
        basex = 0.18*headrad;
        tip   = 1.15*headrad;
        tiphw = .04*headrad;
        tipr  = .01*headrad;

        q = .04;
        EarX  = [.497-.005 .510 .518 .5299 .5419 .54 .547 .532 .510 .489-.005];
        EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];

        sf = headrad/plotrad;

        plot3([basex;tiphw;0;-tiphw;-basex]*sf, [base;tip-tipr;tip;tip-tipr;base]*sf, ...
            2*ones(5,1), 'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');

        plot3(EarX*sf+.18,  EarY*sf, 2*ones(size(EarX)), 'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')
        plot3(-EarX*sf-.18, EarY*sf, 2*ones(size(EarY)), 'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')
    end
end