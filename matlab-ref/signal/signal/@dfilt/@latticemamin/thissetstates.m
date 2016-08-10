function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2004/04/12 23:58:39 $

if isempty(S),
    S = [nullstate1(Hd.filterquantizer);nullstate1(Hd.filterquantizer)];
else
    % Check data type, quantize if needed
    S = validatestates(Hd.filterquantizer, S);
    S = [S;zeros(2,size(S,2))];
end

Hd.HiddenStates = S;
S = [];

