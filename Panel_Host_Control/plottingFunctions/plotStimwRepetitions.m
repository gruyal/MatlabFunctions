function plotStimwRepetitions(pStruct, getStimIndsInput)

% function plotStimwRepatitions(pStruct, getStimIndsInput)
%
% This function plots a overlaid plot for each stimulus (all repetitions on
% one plot). It uses getStimInds for the input, and generates a figure for
% each separate stim. 
%
% INPUT
% pStruct -             protocolStruct w/ stim field and data in it
% getStimIndsInput -    input for getStimInds. 1X4 vector describing the
%                       desired relInds (for more see getStimInds)

indsSt = getStimInds(pStruct, getStimIndsInput);

relCh = 2; 

for ii=1:length(indsSt)
    
    figure
   
    set(gcf, 'NumberTitle', 'off','Name', num2str(indsSt(ii).val))
    hold on
    
    numRep = length(indsSt(ii).inds);
    linCol = cbrewer('seq', 'Greys', numRep+2); 
    
    for jj=1:numRep
        plot(pStruct.stim(indsSt(ii).inds(jj)).data{1}(:,relCh), 'color', linCol(jj+1, :));
    end
    
    legend(arrayfun(@num2str, indsSt(ii).inds, 'uniformoutput', 0))
    
end


end
