function delete(obj)

%   Author: T. Krauss
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.7.4.1 $

obj = struct(obj);
for i = 1:length(obj)
    fig = get(obj(i).h,'parent');
    ud = get(fig,'userdata');
    
    ind = find(obj(i).h == [ud.Objects.fdmeas.h]);
    if isempty(ind)
        error(generatemsgid('InvalidParam'),'Not a valid object')
    end
        
    objud = get(obj(i).h,'userdata');
    
    delete(objud.hlabel)
    delete(obj(i).h)
    
    ud.Objects.fdmeas(ind) = [];
    
    set(fig,'userdata',ud)
end

