function y = tile_or_truncate(x, L)
    if isempty(x), y = zeros(L,1); return; end
    n = numel(x);
    if n >= L
        y = x(1:L);
    else
        r = ceil(L/n);
        y = repmat(x(:), r, 1);
        y = y(1:L);
    end
end

