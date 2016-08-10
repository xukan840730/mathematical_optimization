function hashModel = createHashTable(myStruct,FieldsName,inputArgName)
%CREATEHASHTABLE Create a Java hash table for 'myStruct'.
%   myStruct is a problem or options structure used in Optimization area GUIs.
%   FieldsName are the fields of myStruct for which a Java hash table will be created.
%   The optional argument 'FieldsName' must be a cell array of strings.
%   The optional argument 'inputArgName' is the 'variable name' used to denote 'myStruct'
%   when passed as the input to OPTIMTOOL.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2008/05/23 15:36:36 $

if nargin < 3,
    inputArgName = '';
    if nargin < 2
        FieldsName = fieldnames(myStruct);
    end
end
% Create an empty hashtable
hashModel = java.util.Hashtable;
% Insert key and value in the hashtable
for i = 1:length(FieldsName)
    if isfield(myStruct,FieldsName{i}) && ~isempty(myStruct.(FieldsName{i}))
        rhsValue = getStringValue(myStruct.(FieldsName{i}),FieldsName{i},inputArgName);
        % remove string quotes
        rhsValue(find(rhsValue == '''')) = [];
        hashModel.put(FieldsName{i},rhsValue);
    end
end

% Special code to handle states for random number generators (problem
% structures only). Need to remove 'rngstate' and insert 
% appropriate 'randchoice' key to the hashModel
if isfield(myStruct,'solver') && any(strcmpi(myStruct.solver,{'ga', ...
        'patternsearch','simulannealbnd','threshacceptbnd','gamultiobj'})) && ...
        isfield(myStruct,'rngstate')
    hashModel = setoptimrandstates(myStruct,hashModel,true);
end

%----------------------------------------
function rhsValue = getStringValue(value,FieldName,inputArgName)
% Function to wrap around value2RHS so that large matrices are not
% displayed
MAX_NUM_ELEMENT_SHOW = max(2,10);
if(numel(value) > MAX_NUM_ELEMENT_SHOW) && isnumeric(value)
    rhsValue = sprintf('%s.%s',inputArgName,FieldName);
else
    rhsValue = value2RHS(value);
end
