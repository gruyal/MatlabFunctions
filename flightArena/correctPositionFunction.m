function corrFunc = correctPositionFunction(func, pfnParam)

% function corrFunc = correctPositionFunction(func, pfnParam)
%
% this function is designed to correct position function so that motion is
% continuous between different segments
%
% INPUT
% func      - vector. position function output from Function_Maker_G4
% pfnParam  - struct. parameters structure that was used as input to Function_Maker_G4
%
%
% OUTPUT 
% corrFunc  - vector. same length as func only with no jumps between
%             transitions
%
% checked in sawtooth functions

corrFunc = func; 

if pfnParam.gs_val == 4
    sampFreq = 500; 
else
    sampFreq = 1000; 
end

secFlips = pfnParam.flip;
secHVals = pfnParam.high;
secLVals = pfnParam.low;
secTimes = cumsum(pfnParam.dur * sampFreq); 
secTypeInds = find(~ismember(pfnParam.section, {'static'}));

assert(sum(secTypeInds) > 1,  'only one non-static segment - functio aborted')

for ii = secTypeInds(1:end-1)
    
    relFlip = secFlips(ii);
    endV = corrFunc(secTimes(ii));
    startVNext = corrFunc(secTimes(ii)+1);
    
    if abs(endV - startVNext) ~= 1 % checks that there is a jump
        
        % which position the next section should begin
        if ~relFlip 
            relInc = -1;
        else
            relInc = 1;
        end
        
        relVecSec = corrFunc(secTimes(ii)+1:secTimes(ii+1));
        [spRVSInds, spRVSVals] = SplitVec(relVecSec, 'equal', 'bracket', 'firstval');
        newVecSec = zeros(size(relVecSec));
        currVal = endV;
        
        for jj=1:length(spRVSVals)
           
            newVal = currVal + relInc; 
            if newVal < secLVals(ii)
                newVal = secHVals(ii);
            elseif newVal > secHVals(ii)
                newVal = secLVals(ii);
            end
            
            newVecSec(spRVSInds(jj, 1):spRVSInds(jj, 2)) = newVal; 
            
            currVal = newVal; 

        end
       
        corrFunc(secTimes(ii)+1:secTimes(ii+1)) = newVecSec; 
        
    end
    
end
    



end