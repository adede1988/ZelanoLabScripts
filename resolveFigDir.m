function figDir = resolveFigDir(S)
%RESOLVEFIGDIR  SHARED: the session's figure directory (created if needed).
%   figDir = S.fig if set, else fullfile(S.figPath, S.id).

    if isfield(S, 'fig') && ~isempty(S.fig)
        figDir = S.fig;
    elseif isfield(S, 'figPath') && ~isempty(S.figPath)
        figDir = fullfile(S.figPath, S.id);
    else
        error('resolveFigDir:noFigDir', 'S needs .fig or .figPath.');
    end
    if ~isfolder(figDir), mkdir(figDir); end
end
