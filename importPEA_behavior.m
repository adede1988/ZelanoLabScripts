function T = importPEA_behavior(fp)
%IMPORT_TXT_TABLE  Import a tab-delimited .txt as a table with headers preserved.
%   T = IMPORT_TXT_TABLE(fp)
%   - Preserves original column headers
%   - Handles consecutive tabs / ragged rows
%   - Keeps text columns as string
%   - Strips non-printable control chars from text columns

    arguments
        fp (1,1) string
    end

    % Build import options
    opts = detectImportOptions(fp, ...
        'FileType','text', ...
        'Delimiter','\t', ...
        'Whitespace','');                % don't auto-trim

    opts.PreserveVariableNames = true;
    opts.ExtraColumnsRule = 'ignore';
    opts.ConsecutiveDelimitersRule = 'join';
    % Preserve whitespace / handle empties for all vars
    opts = setvaropts(opts, opts.VariableNames, 'WhitespaceRule','preserve');
    opts = setvaropts(opts, opts.VariableNames, 'EmptyFieldRule','auto');

    % Make text-like vars string type (not char/cellstr)
    vt = opts.VariableTypes;
    textMask = ismember(vt, {'char','string'});
    if any(textMask)
        opts = setvartype(opts, opts.VariableNames(textMask), 'string');
    end

    % Read table
    T = readtable(fp, opts);

    % Clean non-printable characters from string variables
    for vn = T.Properties.VariableNames
        col = T.(vn{1});
        if iscellstr(col), col = string(col); end
        if isstring(col)
            col = regexprep(col, '[^\x20-\x7E]', ''); % keep printable ASCII only
            T.(vn{1}) = col;
        end
    end
end
