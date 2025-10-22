


% Load standard electrode locations
stdLocs = readtable("G:\My Drive\GitHub\AB_pic_analysis/standardEEGlocs.csv");

% Your electrode list
myLabels = {'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7',...
    'TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8',...
    'TP10','CP6','CP2','Cz','C4','T8','FT10','FC6','FC2','F4','F8','Fp2'};



theta = nan(length(myLabels),1);
phi   = nan(length(myLabels),1);

for i = 1:length(myLabels)
    idx = strcmpi(stdLocs.Electrode, myLabels{i});
    if any(idx)
        theta(i) = stdLocs.theta(idx);
        phi(i)   = stdLocs.phi(idx);
    end
end

% Interpolate missing electrodes directly in theta/phi space
% FT9 ~ midpoint of F7 & T7, small phi shift forward
F7 = [theta(strcmpi(myLabels,'F7')), phi(strcmpi(myLabels,'F7'))];
T7 = [theta(strcmpi(myLabels,'T7')), phi(strcmpi(myLabels,'T7'))];
FT9 = mean([F7; T7],1); FT9(2) = FT9(2) - 5; % adjust phi forward
theta(strcmpi(myLabels,'FT9')) = FT9(1);
phi(strcmpi(myLabels,'FT9'))   = FT9(2);

% FT10 ~ midpoint of F8 & T8, small phi shift forward
F8 = [theta(strcmpi(myLabels,'F8')), phi(strcmpi(myLabels,'F8'))];
T8 = [theta(strcmpi(myLabels,'T8')), phi(strcmpi(myLabels,'T8'))];
FT10 = mean([F8; T8],1); FT10(2) = FT10(2) - 5;
theta(strcmpi(myLabels,'FT10')) = FT10(1);
phi(strcmpi(myLabels,'FT10'))   = FT10(2);

% TP9 ~ midpoint of P7 & T7, small phi shift backward
P7 = [theta(strcmpi(myLabels,'P7')), phi(strcmpi(myLabels,'P7'))];
TP9 = mean([P7; T7],1); TP9(2) = TP9(2) + 5; % adjust phi posterior
theta(strcmpi(myLabels,'TP9')) = TP9(1);
phi(strcmpi(myLabels,'TP9'))   = TP9(2);

% TP10 ~ midpoint of P8 & T8, small phi shift backward
P8 = [theta(strcmpi(myLabels,'P8')), phi(strcmpi(myLabels,'P8'))];
TP10 = mean([P8; T8],1); TP10(2) = TP10(2) + 5;
theta(strcmpi(myLabels,'TP10')) = TP10(1);
phi(strcmpi(myLabels,'TP10'))   = TP10(2);

% Combine into final table
outTable = table(myLabels(:), theta, phi, ...
                 'VariableNames', {'Label','Theta','Phi'});

% Save as CSV
writetable(outTable,'G:\My Drive\GitHub\ZelanoLabScripts\myEEGcoords_thetaPhi.csv');