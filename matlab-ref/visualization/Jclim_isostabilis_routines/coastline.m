%%
clear all

%%
ncload('is.08000.01.01.nc');
tempsize = size(output_field1);
cum(1:tempsize(1),:,:) = output_field1;
tot(1:tempsize(1),:,:) = output_field2;

clear output_field1 output_field2

%%

ncload('is.01950.01.01.nc');
tempsize = size(output_field1);
cum((end+1):(end+tempsize(1)),:,:) = output_field1;
tot((end+1):(end+tempsize(1)),:,:) = output_field2;

clear output_field1 output_field2 tempsize

%%
ncload('kmt_orig.nc');

tempsize = size(kmt);
kmt = kmt(2:(tempsize(1)-1),2:(tempsize(2)-1));

clear xt yt tempsize


%%
array_dim = size(cum);
mask=zeros(array_dim(1), array_dim(2));

for t=1:array_dim(1)
for i=1:array_dim(2);
for j=1:array_dim(3);

    %more than melt period over land
  if kmt(i,j) < 0.5;
    if cum(t,i,j) > 0.5;
       mask(t,i,j) = 1;
    end
  end
  
    %more than melt period over ocean
  if kmt(i,j) > 0.5;
    if cum(t,i,j) > 0.5;
       mask(t,i,j) = 2;
    end
  end
  
    %less than melt period over land
  if kmt(i,j) < 0.5;
    if cum(t,i,j) < 0.5;
       mask(t,i,j) = 3;
    end
  end
  
    %less than melt period over ocean
  if kmt(i,j) > 0.5;
    if cum(t,i,j) < 0.5;
       mask(t,i,j) = 4;
    end
  end 
  
    
end
end
end

%%

figure

tempsize = size(mask)
for n=1:tempsize(1);
  date = 1850+n
  date_text = num2str(date);
  plot_var(:,:)=mask(n,:,:);
  surface(plot_var), shading flat;
  text(5,5,date_text)
  M(n) = getframe;
  clf
end

movie2avi(M,'cumulative.avi','FPS',5)





