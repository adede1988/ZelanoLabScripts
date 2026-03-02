
function m = robust_mad(x)
x = x(:);
x = x(isfinite(x));
if isempty(x), m = NaN; return; end
med = median(x);
m = median(abs(x-med));
end
