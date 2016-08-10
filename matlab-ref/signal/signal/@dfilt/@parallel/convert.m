function Hd2 = convert(Hd,newstruct)
%CONVERT Convert structure of DFILT object.
%  Hd2 = CONVERT(Hd,NEWSTRUCT) converts all stages of multisection DFILT
%  object Hd to the structure defined by string NEWSTRUCT.
%
%  EXAMPLE:
%          Hd1 = dfilt.df1 + dfilt.df2;
%          Hd2 = convert(Hd1,'df1');
%          % returns Hd2 with all stages as direct-form 1 discrete-time filters.
%  
%   See also DFILT.

%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.2.4.4 $  $Date: 2008/09/13 07:14:39 $ 

error(nargchk(2,2,nargin,'struct'))

Hd2 = convert(Hd.Stage(1),newstruct);
for k=2:nstages(Hd)
  Hd2 = parallel(Hd2, convert(Hd.Stage(k),newstruct));
end
