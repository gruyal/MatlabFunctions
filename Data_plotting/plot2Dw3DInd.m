function varargout = plot2Dw3DInd(plotDat, relInds)

% function varrgout = plot2Dw3DInd(plotDat)
%
% This function plots 2d data with a plot3 and rotates it to show only xy
% the Z axis is used to index into the order of the original 2d Data
%
% For now only plots with 'o' option
%
% INPUT
% plotDat -     2d data will be plotted as x y plot. can be either 2XN or NX2
% relInds -     (optional) relevant indices to present using cursor. This
%               option is to be used when the data is not index
%               sequentially (e.g. when extracted from a bigger table)
%               If not given the sequence is 1:length(plotDat)
%
%
% OUTPUT
% axh -         optional. axes handles for the plot


dSiz = size(plotDat);

if dSiz(1) == 2 && dSiz(2) > 2
    dat1 = plotDat(1, :);
    dat2 = plotDat(2, :); 
    relS = dSiz(2);
elseif dSiz(2) == 2 && dSiz(1) > 2
    dat1 = plotDat(:,1);
    dat2 = plotDat(:,2); 
    relS = dSiz(1);
else
    error('data not compatible with function')
end

plotInd = 1:relS; 

if nargin == 2
    assert(length(plotInd) == length(relInds), 'relInds are not compatible with plotDat - wrong length')
    plotInd = relInds; 
end
 

lh = plot3(dat1, dat2, plotInd, 'o', 'markeredgecolor', [1,1,1]*0.8, 'markerfacecolor', [1,1,1]*0.8);
axh = lh.Parent; 
view(0, 90)
datacursormode on
dcH = datacursormode(axh.Parent);

set(dcH,'UpdateFcn',@myUpdateFcn)

if nargout == 1
    varargout{1} = axh;
end


end
    

% function txt = myUpdateFcn(~, event_obj)
% % Customizes text of data tips
% 
% pos = get(event_obj,'Position');
% txt = {['x: ',num2str(exp(pos(1)), '%4.3f')], ...
% 	   ['y: ',num2str(exp(pos(2)), '%4.3f')], ...
%        ['Ind:', num2str(pos(3), '%d')]};
%       
% end

function txt = myUpdateFcn(~, event_obj)
% Customizes text of data tips

pos = get(event_obj,'Position');
txt = {['x: ',num2str(pos(1), '%4.3f')], ...
	   ['y: ',num2str(pos(2), '%4.3f')], ...
       ['Ind:', num2str(pos(3), '%d')]};
      
end


