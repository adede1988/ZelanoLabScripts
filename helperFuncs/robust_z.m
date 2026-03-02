

function z = robust_z(x)
x = x(:);
med = median(x,'omitnan');
m = robust_mad(x);
if ~isfinite(m) || m==0
    z = (x - med);
else
    z = (x - med) ./ (1.4826*m);
end
end
