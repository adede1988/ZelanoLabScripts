function [m2, isBad] = topo_fix_isolated_electrodes(m, x, y, varargin)
% topo_fix_isolated_electrodes
% Detect electrodes whose value is far from their spatial neighbors and replace
% them with the (simple) mean of neighboring values.
%
% NEW: optional input 'isBad' lets you reuse a detection mask from one map
% and apply the same replacements to another map.
%
% Inputs:
%   m [n x 1] values
%   x [n x 1] x-coords
%   y [n x 1] y-coords
%
% Options (name-value):
%   'K'        (default 6)    number of nearest neighbors
%   'Z'        (default 4)    robust z-threshold vs neighbor median/MAD (used if isBad not provided)
%   'AbsThresh'(default [])   optional absolute threshold fallback (same units as m)
%   'isBad'    (default [])   logical mask (n x 1) or index vector of electrodes to replace
%
% Outputs:
%   m2    corrected values
%   isBad logical mask of corrected electrodes (detected or passed in)

p = inputParser;
p.addParameter('K', 6, @(v)isscalar(v)&&v>=2);
p.addParameter('Z', 4, @(v)isscalar(v)&&v>0);
p.addParameter('AbsThresh', [], @(v)isempty(v)||(isscalar(v)&&v>0));
p.addParameter('isBad', [], @(v) isempty(v) || islogical(v) || isnumeric(v));
p.parse(varargin{:});

K        = p.Results.K;
Z        = p.Results.Z;
AbsThresh= p.Results.AbsThresh;
isBadIn  = p.Results.isBad;

m = double(m(:)); x = double(x(:)); y = double(y(:));
n = numel(m);
m2 = m;

% Pairwise distances (no toolboxes)
D = hypot(x - x.', y - y.');
D(1:n+1:end) = inf; % ignore self

% Convert isBad input to logical mask if provided
if ~isempty(isBadIn)
    if islogical(isBadIn)
        isBad = isBadIn(:);
    else
        isBad = false(n,1);
        ii = isBadIn(:);
        ii = ii(isfinite(ii) & ii>=1 & ii<=n);
        isBad(ii) = true;
    end
else
    isBad = false(n,1);

    % Global robust scale fallback (in case local MAD ~ 0)
    gmed = median(m, 'omitnan');
    gmad = median(abs(m - gmed), 'omitnan');
    gmad = max(gmad, eps);

    % --- detect outliers vs neighbors ---
    for i = 1:n
        if ~isfinite(m(i)), continue; end

        [~, ord] = sort(D(i,:), 'ascend');
        nn = ord(1:min(K, n-1));
        v  = m(nn);
        v  = v(isfinite(v));

        if numel(v) < 3, continue; end

        medv = median(v, 'omitnan');
        madv = median(abs(v - medv), 'omitnan');
        scale = madv;
        if scale < 1e-12
            scale = gmad;
        end

        rz = abs(m(i) - medv) / (1.4826*scale); % MAD->sigma approx
        bad = rz > Z;

        if ~isempty(AbsThresh)
            bad = bad | (abs(m(i) - medv) > AbsThresh);
        end

        if bad
            isBad(i) = true;
        end
    end
end

% --- replace bad electrodes using mean of nearby GOOD electrodes ---
goodMask = ~isBad & isfinite(m);

for i = find(isBad(:))'
    % neighbor candidates: prefer good electrodes
    cand = find(goodMask);
    if isempty(cand)
        continue;
    end

    % distances from i to candidates
    di = D(i, cand);
    [~, ord] = sort(di, 'ascend');

    nn = cand(ord(1:min(K, numel(cand))));
    v  = m(nn);
    v  = v(isfinite(v));

    if isempty(v)
        continue;
    end

    m2(i) = mean(v, 'omitnan');  % simple average of surrounding values
end

end