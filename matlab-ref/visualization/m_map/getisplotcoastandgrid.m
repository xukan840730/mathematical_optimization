%% Get m_map plotting information, based on which ice sheet to plot
m_coast('color','k','linewidth',1);
is=strcmp(isname{isn},'ais');
if is==1;  
   m_grid('XaxisLocation', 'top', 'xtick',12,'tickdir','out','ytick',[])
else
    is=strcmp(isname{isn},'gis');
    if is==1;
       m_grid('XaxisLocation', 'bottom', 'xtick',12,'tickdir','out','YaxisLocation', 'middle','ytick',[])
    else
        display('ERROR NO ICE SHEET DEFINED')
    end
end
