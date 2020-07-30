function plotModelCondProfile(cellDat, normFlag, axh)

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
relXX = -6:6; 

relCol = cbrewer('qual', 'Set1', 3); 

if parLen == 12

    ampE = cellDat.params(4);
    muE = cellDat.params(5);
    sigE = cellDat.params(6);
    
    tt = cellDat.params(7); vv = cellDat.params(8); tau = cellDat.params(9); lam = cellDat.params(10); ampI = cellDat.params(11);
    
%     ampI = cellDat.params(7);
%     muI = cellDat.params(9);
%     sigI = cellDat.params(11);

    yyE = ampE * normpdf(relXX, muE, sigE); 
%     yyI = ampI * normpdf(relXX, muI, sigI); 
    yyI = I_profile(relXX, tt, vv, tau, lam, ampI); 

elseif parLen == 21

    ampE = cellDat.params(5);
    muE = cellDat.params(6);
    sigE = cellDat.params(7);
    
    profI = cellDat.params(8:20);

    yyE = ampE * normpdf(relXX, muE, sigE); 
    yyI = profI; 
    
elseif parLen == 13

    ampE = cellDat.params(5);
    muE = cellDat.params(6);
    sigE = cellDat.params(7);
    
    tt = cellDat.params(8); vv = cellDat.params(9); tau = cellDat.params(10); lam = cellDat.params(11); ampI = cellDat.params(12);

    yyE = ampE * normpdf(relXX, muE, sigE); 
    yyI = I_profile(relXX, tt, vv, tau, lam, ampI); 
    
elseif parLen == 14

    
    tE = cellDat.params(4); vE = cellDat.params(5); tauE = cellDat.params(6); lamE = cellDat.params(7); ampE = cellDat.params(8);
    tI = cellDat.params(9); vI = cellDat.params(10); tauI = cellDat.params(11); lamI = cellDat.params(12); ampI = cellDat.params(13);

    yyE = EI_profile(relXX, tE, vE, tauE, lamE, ampE);
    yyI = EI_profile(relXX, tI, vI, tauI, lamI, ampI);
    

else 
    
    error('wrong number of parameters')
end

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




