function newTable = addColumnToGratingTable(oldTable, varName, varValue)

% function newTable = addColumnToGratingTable(oldTable, varName, varValue)
%
% This function adds a variable to a table, either as a single value for
% the whole table.
%
%
% INPUT
% oldTable -        table with no variable with the same name as varName
% varName -         string. name for the column.
% varValue -        single value the will be applied for the whole table
%
% OUTPUT
% newTable -        same as old with new var added 

allNames = oldTable.Properties.VariableNames;
assert(~ismember(varName, allNames), 'variable name already exists in the table')

tabH = height(oldTable); 

newCol = repmat(varValue, tabH, 1);

tempTab = table(newCol, 'VariableNames', {varName});

newTable = [oldTable, tempTab];







end


