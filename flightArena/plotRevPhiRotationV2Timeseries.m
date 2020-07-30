function plotRevPhiRotationV2Timeseries(processedDataFileName)

% function plotRevPhiRotationV2Timeseries(direcName)
%
% This function plots the timeseries of the LmR channel and the freq
% channel from the protocol revPhiRotationV202-26-20_10-24-23 
% (contains 2bar widths and phi and revPhi motion)
%
% processedDataFileName is the file name containing experimental results. 
% This includesmatrix timeseries, timestamps, channelNames etc. (generated
% at the end of the experiment by the processing function) 


load(processedDataFileName, 'timeseries', 'timestamps')

relCh = [2,3,6]; % LmR and freq
yyLim = [0, 7.5; 0, 7.5; -6.1, 6.1];
xxLim = [0, 5]; 

winSize = 21; % for smoothing 

matSiz = size(timeseries);

plotTit = {'phi', 'revPhi'};

numReps = matSiz(3);
numStim = matSiz(2); 
numCh = matSiz(1); 

relCol = cbrewer('seq', 'YlOrRd', numReps+2);

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.01, 0.02, [2, numStim/2]);
axh = gobjects([length(relCh), size(posCell)]); 

for ff = 1:length(relCh)
    
    fh = figure('position', [1400, 200, 700, 1100]);
    chInd = relCh(ff); 
    
    for ss=1:numStim
        
        xInd = 2-mod(ss, 2);
        yInd = ceil(ss/2);
        
        axh(ff, xInd, yInd) = axes('position', posCell{xInd, yInd}); 
        hold on
        chData = squeeze(timeseries(chInd, ss, :, :)); 
        
        filtChData = smoothdata(chData, 2, 'movmean', winSize); 
        
        plot(timestamps, filtChData)
        
        axh(ff,xInd,yInd).Box = 'off'; 
        
        lineH = flipud(findobj(axh(ff,xInd,yInd).Children, 'LineStyle', '-'));
        
        for ll=1:length(lineH)
            lineH(ll).Color = relCol(ll+2, :);
        end
        
        plot(timestamps, mean(filtChData), 'linewidth', 2, 'color', 'k')
        
        line(xxLim, [0,0], 'linewidth', 1, 'color', 'k')
        hold off
        
        axh(ff,xInd,yInd).XLim = xxLim; 
        
        axh(ff,xInd,yInd).YLim = yyLim(ff, :); 
        
        if xInd>1 && ff ==1
            axh(ff,xInd,yInd).YColor = 'none'; 
        end
        
        if yInd < size(axh,3)
            axh(ff,xInd,yInd).XColor = 'none'; 
        end
        
        if yInd == 1
            axh(ff, xInd,yInd).Title.String = plotTit{xInd}; 
        end
            
        axh(ff, xInd,yInd).Children = flipud(axh(ff, xInd,yInd).Children);
        
    end
    
end
        




end