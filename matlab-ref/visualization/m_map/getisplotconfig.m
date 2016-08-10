%% Get m_map plotting information, based on which ice sheet to plot
lat_array=linspace(-90.,90);lon_array=linspace(0,360); 
[Plg,Plt]=meshgrid(lon_array,lat_array);
is=strcmp(isname{isn},'ais');
if is==1;
    m_proj('stereographic','longitude', 0,'latitude',-90,'radius', 30);
else
    is=strcmp(isname{isn},'gis');
    if is==1;
       m_proj('stereographic','longitude', 310,'latitude',73,'radius', 15);  
    else
        display('ERROR NO ICE SHEET DEFINED')
    end
end
