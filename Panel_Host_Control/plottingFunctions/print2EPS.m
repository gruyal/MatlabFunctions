function print2EPS(fileName, fh)

% function print2EPS(fileName)
%
% This function prints the current figure to a nicely sized EPS figure
% fh is optional figure handle (if not given uses gcf) 
% modified from print2PDF

if nargin < 2
    fh = gcf;
end

set(fh,'Units','inches');
screenposition = get(fh,'Position');
set(fh,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print(fh, fileName, '-depsc2', '-painters', '-r0') 


end