function y = local_minmax01(v)
    v = v(:);
    mn = min(v(isfinite(v)));
    mx = max(v(isfinite(v)));
    if isempty(mn) || isempty(mx) || mx <= mn
        y = nan(size(v));
        if any(isfinite(v)), y(isfinite(v)) = 0; end
        return
    end
    y = (v - mn) ./ (mx - mn);
end