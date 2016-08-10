function [NL, NextIPorts, NextOPorts, mainparams]=df1headconnect(q,NL,H,mainparams)
%DF1HEADCONNECT specifies the blocks, connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/12/26 22:05:39 $


% specify the qparam

set(NL.nodes(1),'qparam','double');
set(NL.nodes(2),'qparam','double');


% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(3,1,1,1);
NL.connect(2,1,4,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[filtgraph.nodeport(5,1) filtgraph.nodeport(2,1)];
NextOPorts=[filtgraph.nodeport(3,1) filtgraph.nodeport(1,1) filtgraph.nodeport(5,1) filtgraph.nodeport(2,1)];

