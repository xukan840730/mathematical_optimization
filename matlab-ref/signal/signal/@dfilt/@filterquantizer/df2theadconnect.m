function [NL, NextIPorts, NextOPorts, mainparams]=df2theadconnect(q,NL,H,mainparams)
%DF2THEADCONNECT specifies connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/12/26 22:06:09 $


% specify the qparam
% gain
set(NL.nodes(2),'qparam','double');
set(NL.nodes(4),'qparam','double');
% sum
set(NL.nodes(3),'qparam','double');


% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(3,1,4,1);
NL.connect(4,1,5,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[filtgraph.nodeport(3,2)];
NextOPorts=[filtgraph.nodeport(1,1) filtgraph.nodeport(4,1)];

