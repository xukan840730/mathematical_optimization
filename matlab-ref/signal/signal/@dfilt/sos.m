%SOS  Convert to second-order-sections (for IIR filters only).
%   Hsos = SOS(Hd) converts IIR discrete-time filter Hd to second-order section
%   form. Hd can be of class dfilt.df1, dfilt.df1t, dfilt.df2 or
%   dfilt.df2t.
%
%   SOS(Hd,DIR_FLAG) specifies the ordering of the 2nd order sections. If
%   DIR_FLAG is equal to 'UP', the first row will contain the poles closest
%   to the origin, and the last row will contain the poles closest to the
%   unit circle. If DIR_FLAG is equal to 'DOWN', the sections are ordered
%   in the opposite direction. The zeros are always paired with the poles
%   closest to them. DIR_FLAG defaults to 'UP'.
%
%     % Example:
%     [b,a] = butter(8,.5);
%     Hd = dfilt.df2(b,a);
%     Hsos = sos(Hd,'up',inf)
%   
%     % Example 2: Note that if the Filter Design Toolbox is available, it
%     % is preferable todesign directly in SOS form rather than to construct
%     % the transfer function and then convert which can introduce roundoff
%     % error
%     f = fdesign.lowpass('N,F3db',8,0.5);
%     Hsos = design(f,'butter','FilterStructure','df1sos');
%
%   See also DFILT/SCALE in the Filter Design Toolbox.

%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.4.2 $  $Date: 2006/06/27 23:33:05 $ 

% Help for the p-coded SOS method of DFILT.DF1, DFILT.DF1T, DFILT.DF2 and
% DFILT.DF2T classes.
