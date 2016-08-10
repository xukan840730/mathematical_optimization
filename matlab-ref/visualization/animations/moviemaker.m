function moviemaker(var,fname)

arr_sz = size(var);
for n=1:arr_sz(1);
  plot_var(:,:)=var(n,:,:);
  figure
  surface(plot_var), shading flat;
  M(n) = getframe;
  close all
end
hello = 'hello'
movie(M)

movie2avi(M,fname)