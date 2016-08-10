function msg = fdatoolmessage(id, varargin)
%FDATOOLMESSAGE Translates the string from an ID.

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/05/23 08:16:23 $

[~, catalogInfo] = DAStudio.getMessageCatalogInfo('signal', 'fdatool');

msg = feval(catalogInfo.fcnHndl, id);

if nargin > 1
    msg = sprintf(msg, varargin{:});
end

% [EOF]
