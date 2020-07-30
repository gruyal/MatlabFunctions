function plotModel9CondProfile(cellDat, normFlag, axh)

% function plotModelCondProfile(cellDat)
%
% This function plots the E and I conductance profile between positions -6
% to 6
% it assumes a certain structure to the parameters vector.
% if its length is 12 - it will plot 2 gaussians
% if the length is 21 - it will plot a gaussian for E and position values
% for I
%
% INPUT 
% axh -             axes handle to plot in 
% cellDat -         generated with organizingClusterData
% normFlag -        (optional) if true normalizes the E and I cond
% (default)


if nargin < 2
    normFlag = 1;
end

if nargin < 3
    fh  = figure;
    axh  = gca; 
else
    fh = axh.Parent; 
end

parLen = length(cellDat.params);
assert(parLen ==14, 'model has wrong number of parameters')
relXX = -6:6; 

relCol = cbrewer('qual', 'Set1', 3); 

ampE = cellDat.params(6);
muE = cellDat.params(7);
sigE = cellDat.params(8);

tt = cellDat.params(9); vv = cellDat.params(10); tau = cellDat.params(11); lam = cellDat.params(12); ampI = cellDat.params(13);

yyE = ampE * normpdf(relXX, muE, sigE); 
yyI = EI_profile(relXX, tt, vv, tau, lam, ampI); 


if normFlag
    yyE = yyE ./ max(yyE);
    yyI = yyI ./ max(yyI); 
end

fh.CurrentAxes = axh; 

hold on 
plot(relXX, yyE, 'linewidth', 2, 'color', relCol(1, :))

plot(relXX, yyI, 'linewidth', 2, 'color', relCol(2, :))
% plot(relXX, yyE-yyI, 'linewidth', 2, 'color', relCol(3, :))
hold off

legend('E profile', 'I profile')

end




