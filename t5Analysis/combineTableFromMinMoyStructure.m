function [fullTab, emptyTab] = combineTableFromMinMoyStructure(minMotLinSt)

% function [fullTab, emptyTab] = combineTableFromMinMoyStructure(minMotLinSt)
%
% This function is an accessory function designed to work with the output
% of calcMinMotExtLinCompDiffWandV. it goes over all the structure and
% combines it back into 2 tables. One where linear data has been
% reconstructed and one where it hasn't 
%
% INPUT 
%
% minMotLinSt -         output from calcMinMotExtLinCompDiffWandV. 6D
%                       struture
%
% OUTPUT
%
% fullTab -             table from all the components of the table where a
%                       linear comparison has been made
% emptyTab -            complementary table

dSiz = size(minMotLinSt); 


fullTab = table;
fullTab2 = table;
emptyTab = table;
emptyTab2 = table;

for ii=1:dSiz(1)
    
    for jj=1:dSiz(2)
        
        for kk=1:dSiz(3)
            
            for mm=1:dSiz(4)
                
                for nn=1:dSiz(5)
                    
                    for tt=1:dSiz(6)
                        
                        if minMotLinSt(ii,jj,kk,mm,nn,tt).empty && isempty(minMotLinSt(ii,jj,kk,mm,nn,tt).data)
                            
                            emptyTab2 = vertcat(emptyTab2, table(ii,jj,kk,mm,nn,tt));
                            
                        elseif minMotLinSt(ii,jj,kk,mm,nn,tt).empty && ~isempty(minMotLinSt(ii,jj,kk,mm,nn,tt).data)
                            
                            emptyTab = vertcat(emptyTab, minMotLinSt(ii,jj,kk,mm,nn,tt).data.table);
                            
                        elseif minMotLinSt(ii,jj,kk,mm,nn,tt).empty == 0
                            
                            tempTab = minMotLinSt(ii,jj,kk,mm,nn,tt).data.table;
                            
                            if width(tempTab) == 14
                            
                                fullTab = vertcat(fullTab, tempTab);
                                
                            else
                                
                                fullTab2 = vertcat(fullTab2, tempTab);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end






end
                        
                        
                        
                        
                        