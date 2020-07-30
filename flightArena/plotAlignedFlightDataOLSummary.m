function varargout = plotAlignedFlightDataOLSummary(alignExpStruct)

% function varargout = plotAlignedFlightDataOL(alignExpStruct)
% 
% This function uses the output from alignFlightArenaExpStruct to plot the
% data from stimulus stimNum. Plot include recalculated L-R, L+R and Freq
% data for individual trials and mean. 
%
% 

numOL = size(alignExpStruct.dataOL,1);

relTab = [];
for ii=1:numOL
    tempTab = alignExpStruct.dataOL(ii,1).table;
    relTab = vertcat(relTab, tempTab);
end

relTags = {'onLev', 'freq'};

numSpeeds = length(unique(relTab{:, relTags{2}}));
numONLev = length(unique(relTab{:, relTags{1}}));
xxLim = [-1, 11];
yyLim = [-5, 4];
yyLim2 = [4,10];
% yyLim3 = [-0.1, 2.5];
posSwitch = 0:2.5:10; 

close all

timeCh = 1; lCh = 3; rCh = 4; freqCh = 5;  % ch 2 is L-R from the flight controller (is used for closed loop but less accurate)

rCols = cbrewer('qual', 'Paired', 8);

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.02, 0.02, [numSpeeds, numONLev]); 
axh = gobjects(3, numSpeeds, numONLev);



fh = figure('name', 'L-R', 'position', [1100, 400, 1000, 800]); 
fh2 = figure('name', 'L+R', 'position', [1100, 400, 1000, 800]); 
fh3 = figure('name', 'Freq', 'position', [1100, 400, 1000, 800]); 

for ii=1:numSpeeds
    
    for jj=1:numONLev
        
        stimNum = numONLev*(ii-1) + jj;
        
        plotTag = [];
        for kk=1:length(relTags)
            plotTag = [plotTag, alignExpStruct.dataOL(stimNum,1).table{:,relTags{kk}}];
        end
        
        relData = alignExpStruct.dataOL(stimNum, end);
        

        excludedT = ~relData.usefulFlag;

        datTime = relData.alignData(:,:,timeCh);
        lMinusR = relData.alignData(:,:, lCh) - relData.alignData(:,:, rCh);
        lPlusR = relData.alignData(:,:, lCh) + relData.alignData(:,:, rCh);
        freqDat = relData.alignData(:,:, freqCh); 


        meanTime = relData.meanData(:,timeCh);
        meanLminR = relData.meanData(:,lCh) - relData.meanData(:, rCh);
        meanLplusR = relData.meanData(:,lCh) + relData.meanData(:, rCh);
        meanFreq = relData.meanData(:,freqCh);

        exFlag = 0;
        if sum(excludedT > 0)
            exTime = datTime(:, excludedT);
            exLMR = lMinusR(:,excludedT);
            exLPR = lPlusR(:,excludedT);
            exFreq = freqDat(:, excludedT);
            exFlag = 1;
        end

        posData = relData.expPos;
        set(0, 'currentFigure', fh)

        axh(1, ii, jj) = axes('position', posCell{ii,jj});
        hold on 
        line(xxLim, [0, 0], 'color', [1,1,1]*0.8)
        plot(datTime, lMinusR, 'linewidth', 2, 'color', rCols(1, :))
        
        title(num2str(plotTag))
        
        for kk=1:length(posSwitch)-1
            if mod(kk,2)
                stCol = 'k';
            else
                stCol = [1,1,1]*0.75;
            end
            
            posTempInds = posData(:,1) > posSwitch(kk) & posData(:,1) < posSwitch(kk+1);
            
            plot(posData(posTempInds,1), posData(posTempInds,2) ./10 + yyLim(1) + 0.5, ...
                 'color', stCol, 'linewidth', 2)
        end
            
        
        plot(meanTime, meanLminR, 'linewidth', 2, 'color', rCols(2, :))

        if exFlag
            plot(exTime, exLMR, 'linewidth', 1, 'color', [1,1,1]*0.6);
        end

        hold off

        axh(1,ii,jj).YLim = yyLim;
        axh(1,ii,jj).XLim = xxLim;
        
        set(0, 'currentFigure', fh2)
        
        axh(2, ii, jj) = axes('position', posCell{ii,jj});
        hold on 
        line(xxLim, [0, 0], 'color', [1,1,1]*0.8)
        plot(datTime, lPlusR, 'linewidth', 2, 'color', rCols(3, :))

        for kk=1:length(posSwitch)-1
            if mod(kk,2)
                stCol = 'k';
            else
                stCol = [1,1,1]*0.75;
            end
            line([posSwitch(kk), posSwitch(kk+1)], [yyLim2(1)+0.5, yyLim2(1)+0.5], ...
                 'color', stCol, 'linewidth', 4)
        end
            
        
        plot(meanTime, meanLplusR, 'linewidth', 2, 'color', rCols(4, :))

        if exFlag
            plot(exTime, exLPR, 'linewidth', 1, 'color', [1,1,1]*0.6);
        end

        hold off

        axh(2,ii,jj).YLim = yyLim2;
        axh(2,ii,jj).XLim = xxLim;
        
        
        set(0, 'currentFigure', fh3)
        
        axh(3, ii, jj) = axes('position', posCell{ii,jj});
        hold on 
        plot(datTime, freqDat, 'linewidth', 2, 'color', rCols(5, :))

        for kk=1:length(posSwitch)-1
            if mod(kk,2)
                stCol = 'k';
            else
                stCol = [1,1,1]*0.75;
            end
%             line([posSwitch(kk), posSwitch(kk+1)], [yyLim3(1)+0.5, yyLim3(1)+0.5], ...
              line([posSwitch(kk), posSwitch(kk+1)], [2, 2], ...
                 'color', stCol, 'linewidth', 4)
        end
            
        plot(meanTime, meanFreq, 'linewidth', 2, 'color', rCols(6, :))

        if exFlag
            plot(exTime, exFreq, 'linewidth', 1, 'color', [1,1,1]*0.6);
        end

        hold off

%         axh(3,ii,jj).YLim = yyLim3;
        axh(3,ii,jj).XLim = xxLim;
        
        
        
    end
    
end


fh4 = figure('name', 'Fixation', 'position', [1100, 400, 1000, 800]); 

plotAlignFlightCLNew(alignExpStruct)


if nargout > 0
    
    varargout{1} = axh;
    
end



end



