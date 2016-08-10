function NL = nodelist(N)
%NODE Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/14 15:13:35 $

error(nargchk(0,1,nargin,'struct'));

NL = filtgraph.nodelist;

if nargin > 0
    NL.nodeCount = N;
end

if  NL.nodeCount> 0
    for I = 1:NL.nodeCount
        X(I) = filtgraph.node('DUMMY');
        X(I).setindex(I);
    end
else
    X = [];
end

NL.nodes = X;
end
