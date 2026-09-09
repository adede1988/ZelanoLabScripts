function g = v2_grouplabel(type, sessID)
% V2_GROUPLABEL  Analysis group from cohort type + session number.
%   Dupi is split by visit/session number (trailing _N in the sessID): DupiS1/DupiS2/DupiS3
%   (others -> DupiSother). OBE (controls) stay a single group. Accepts either the raw
%   outDat.type ('Dupi'/'OBE') OR an already-derived 'Dupi'/'DupiSx' label as `type`.
    t = char(string(type));
    if ~isempty(regexpi(t,'^dupi','once'))
        tok = regexp(char(string(sessID)), '_(\d+)$', 'tokens', 'once');
        n = 0; if ~isempty(tok), n = str2double(tok{1}); end
        if ismember(n,[1 2 3]), g = sprintf('DupiS%d', n); else, g = 'DupiSother'; end
    else
        g = 'OBE';
    end
end
