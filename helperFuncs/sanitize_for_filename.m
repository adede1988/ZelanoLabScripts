
function out = sanitize_for_filename(s)
s = string(s);
out = regexprep(s, '[^\w\-]+', '_');     % keep letters/numbers/_/-
out = regexprep(out, '_+', '_');
out = strip(out, "_");
if strlength(out)==0, out="NA"; end
end

