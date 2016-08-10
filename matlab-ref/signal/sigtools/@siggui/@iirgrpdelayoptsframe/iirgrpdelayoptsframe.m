function h = iirgrpdelayoptsframe(varargin)
%IIRGRPDELAYOPTSFRAME  Constructor for the options frame
%
%   DENSITYFACTOR   -   Value for the density factor
%   MAXPOLERADIUS   -   Value for Max Pole Radius
%   INITDEN         -   Initial guess at denominator
%   NAME            -   Name

%   Author(s): Z. Mecklai
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2.4.1 $  $Date: 2003/03/02 10:27:40 $

%  Builtin-in constructor
h = siggui.iirgrpdelayoptsframe;

% Set the version and tag
set(h, 'version', 1.0);
settag(h);

% [EOF]
