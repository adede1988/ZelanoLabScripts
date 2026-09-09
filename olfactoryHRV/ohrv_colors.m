function [cmap, subs] = ohrv_colors()
%OHRV_COLORS  Fixed patient colour assignment, shared across grant figures.
%   The seven patients with a matched session 1 -> 2 interval. Colours are
%   fixed here so a given patient is the same colour in every figure.
subs = ["JH","AB","GH","DB","KS","JL","PC"];
cmap = [ 0.05 0.40 0.36      % JH  deep teal
         0.85 0.42 0.06      % AB  orange
         0.25 0.47 0.78      % GH  blue
         0.62 0.15 0.42      % DB  magenta
         0.55 0.62 0.10      % KS  olive
         0.35 0.30 0.62      % JL  violet
         0.78 0.18 0.20 ];   % PC  red
end
