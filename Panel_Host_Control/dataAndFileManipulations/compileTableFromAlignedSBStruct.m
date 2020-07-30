function sbTab = compileTableFromAlignedSBStruct(alignSBSt)

% This is a convenience function to go into the structure and generate a
% complied table from all the individual tables in the structure. 
% it adds the fields counter (which is equivalent to index - already
% taken), and inds - which are the indices into the structure)

dSiz = size(alignSBSt); 

if length(dSiz) == 3
    dSiz = [dSiz , 1]; 
end
    
sbTab = [];
counter = 0;

for ii=1:dSiz(1)
    
    for jj=1:dSiz(2)
        
        for kk=1:dSiz(3)
            
            for mm=1:dSiz(4)
                
                if ~isempty(alignSBSt(ii,jj,kk,mm).data)
                    
                    counter=counter+1;
                    tempTab = alignSBSt(ii,jj,kk,mm).data.table; 
                    
                    inds = [ii,jj,kk,mm];
                    
                    sbTab = [sbTab; [table(counter), tempTab, table(inds)]];
                    
                end
                
            end
            
        end
        
    end
    
end

maxEPos = alignSBSt(end, end, end, end).maxExtPosVal;
maxIPos = alignSBSt(end, end, end, end).minInhPosVal;

relPD = sign(maxIPos - maxEPos);

normPos = nan(height(sbTab), 1);

for pp=1:length(normPos)

    if relPD == 1
        normPos(pp) = sbTab.position(pp) - maxEPos;
    else
        normPos(pp) = (sbTab.position(pp) - maxEPos - sbTab.width(pp) + 1) * relPD; 
    end 
end

sbTab = [sbTab, table(normPos)];



end
                    
                    
                    
                