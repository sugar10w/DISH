function result_show(result)
[x,y,z,H] = reducevolume(result,[1,1,1]);
f1 = isosurface(x,y,z,H,0.8,'verbose');
p1=patch(f1,'FaceColor','red','EdgeColor','none');
isonormals(x,y,z,H,p1);
f2=isocaps(x,y,z,H,0.8);
p2=patch(f2,'FaceColor','interp','EdgeColor','none');                                                                           
view(45,45);
axis tight;
daspect([300,300,300]);
colormap(gray(100));
camlight;

lighting phong;
lightangle(90,30);
end

