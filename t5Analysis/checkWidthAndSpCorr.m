
function [newMmProt, protW, protSp] = checkWidthAndSpCorr(mmProt)

% this is a quick function to facilitate combining minMot protocols
% it simply checks that the width of the first and second bars are the same
% and extracts width and spCorr values (since they are not in gratingTable)
%
% Also removes NaNs from fAppear and fDisappear in the table since they do
% not allow proper table comparisons (setdiff treats them as different 
% replaces them with -1

protW = mmProt.inputParams.fBarWid;
assert(protW == mmProt.inputParams.sBarWid, 'bar widths are not the same')
protSp = mmProt.inputParams.speedCor;

nanInds = isnan(mmProt.gratingTable.fAppear);

if any(nanInds) % since they are the same for fAppear and fDisappear
    
    mmProt.gratingTable.fAppear(nanInds) = -1;
    mmProt.gratingTable.fDisappear(nanInds) = -1;
    
    fprintf('repleaced %d NaNs with -1 \n', sum(nanInds))
else
    fprintf('no NaNs were replaced \n')
    
end

newMmProt = mmProt; 

end