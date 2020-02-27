function FEA_field_plot_3D(node,element,Cdata)
%% this function plot the 3D quanticity for hexahedral element

p=patch('Vertices',node,'Faces',element(:,2:5));
set(p,'FaceColor','interp','FaceVertexCData',Cdata,'CDataMapping','scaled')

s1=element(:,[2 3 7 6]);
p=patch('Vertices',node,'Faces',s1);
set(p,'FaceColor','interp','FaceVertexCData',Cdata,'CDataMapping','scaled')

s1=element(:,[3 4 8 7]);
p=patch('Vertices',node,'Faces',s1);
set(p,'FaceColor','interp','FaceVertexCData',Cdata,'CDataMapping','scaled')

s1=element(:,[4 5 9 8]);
p=patch('Vertices',node,'Faces',s1);
set(p,'FaceColor','interp','FaceVertexCData',Cdata,'CDataMapping','scaled')

s1=element(:,[5 2 6 9]);
p=patch('Vertices',node,'Faces',s1);
set(p,'FaceColor','interp','FaceVertexCData',Cdata,'CDataMapping','scaled')

p=patch('Vertices',node,'Faces',element(:,6:9));
set(p,'FaceColor','interp','FaceVertexCData',Cdata,'CDataMapping','scaled')

%% 
daspect([1 1 1]);
view(30,50); 
grid on;
%camlight; lighting gouraud;
alpha(.75)
axis([-60 60 -10 70 -5 15]);
xlabel('X');ylabel('Y');zlabel('Z');colorbar;



end