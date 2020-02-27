function [displacement strain stress von_stress c dc] = Mechanical_analysis(element,node,boundary_displace,point_force,F0,surface_traction,xPhys,loadedele,penal,E0,Emin,mu)

num_element=size(element,1);
num_node=node(end,1);

C=zeros(6,6);
C(1:3,1:3)=[1-mu  mu   mu
             mu 1-mu   mu
             mu   mu 1-mu];
C(4:6,4:6)=[(1-2*mu)/2          0          0
                     0 (1-2*mu)/2          0 
                     0          0 (1-2*mu)/2 ];
C=C*E0/((1+mu)*(1-2*mu));  
                 
b = zeros(3*num_node,1);    
K = sparse(3*num_node,3*num_node);

[Gaussian_point weight]=Gaussian_quadrature_3D('second');              

%% start to calculate the stiffness matrix for all the element 
disp('calculate stiffness matrix')
for i=1:num_element
      node_num   = element(i,2:9);
     XYZ_element = node(node_num,2:4);
     Ke=zeros(24,24);
 
  for j=1:length(weight)
      [B N N_diff J H J_inv x y z] = master_element_linear_elasticity_3D('linear',Gaussian_point(j,1),Gaussian_point(j,2),Gaussian_point(j,3),XYZ_element);     
      Ke = Ke+weight(j)*(N_diff'*J_inv'*H'*C*H*J_inv*N_diff)*det(J(1:3,1:3));
  end
% assembling 
     assembling_index=zeros(24,1);
     assembling_index(1:3:24,1)=3*(node_num-1)+1;
     assembling_index(2:3:24,1)=3*(node_num-1)+2;
     assembling_index(3:3:24,1)=3*node_num;
     K(assembling_index, assembling_index)=K(assembling_index,assembling_index)+Ke*(Emin+xPhys(i).^penal*(E0-Emin))/E0;
end

% add point force
for i=1:size(point_force,1)
  index=3*(point_force(i,1)-1)+point_force(i,2);
  b(index)=b(index)+point_force(i,3);
end

%%  begin to set essential boundary condition
disp('set boundary conditions and solve the linear system')
M=1e20;
for i=1:size(boundary_displace,1)
 index=3*(boundary_displace(i,1)-1)+boundary_displace(i,2);
 K(index,index)=max(abs(K(index,:)))*M;
 b(index)=M*boundary_displace(i,3)*K(index,index);
end

% slove the linear system 
displacement = K\b;
u_displacement=displacement(1:3:end);
v_displacement=displacement(2:3:end);
w_displacement=displacement(3:3:end);

disp('calculate sensitivity')
c = 0; dc = zeros(num_element, 1);
for i=1:num_element
    node_num   = element(i,2:9);
    Ue = [u_displacement(node_num)'; v_displacement(node_num)'; w_displacement(node_num)'];
    Fe = zeros(1,24);
    if F0(3*node_num(1)-2) ~= 0
        Fe(1) = 0.5;  Fe(10) = 0.5; Fe(13) = 0.5; Fe(22) = 0.5;
        Fe(4) = -0.5;  Fe(7) = -0.5; Fe(16) = -0.5; Fe(19) = -0.5;
    elseif F0(3*node_num(7)-2) ~= 0
        Fe(1) = 0.5;  Fe(10) = 0.5; Fe(13) = 0.5; Fe(22) = 0.5;
        Fe(4) = -0.5;  Fe(7) = -0.5; Fe(16) = -0.5; Fe(19) = -0.5;
    elseif F0(3*node_num(1)-1) ~= 0
        Fe(2) = 0.5;  Fe(5) = 0.5; Fe(14) = 0.5; Fe(17) = 0.5;
        Fe(8) = -0.5;  Fe(11) = -0.5; Fe(20) = -0.5; Fe(23) = -0.5;
    elseif F0(3*node_num(7)-1) ~= 0
        Fe(2) = 0.5;  Fe(5) = 0.5; Fe(14) = 0.5; Fe(17) = 0.5;
        Fe(8) = -0.5;  Fe(11) = -0.5; Fe(20) = -0.5; Fe(23) = -0.5;   
    elseif F0(3*node_num(7)) ~= 0 
        Fe(15) = 0.5;  Fe(18) = 0.5; Fe(21) = 0.5; Fe(24) = 0.5;
        Fe(3) = -0.5;  Fe(6) = -0.5; Fe(9) = -0.5; Fe(12) = -0.5;  
    elseif F0(3*node_num(1)) ~= 0 
        Fe(15) = 0.5;  Fe(18) = 0.5; Fe(21) = 0.5; Fe(24) = 0.5;
        Fe(3) = -0.5;  Fe(6) = -0.5; Fe(9) = -0.5; Fe(12) = -0.5;  
    end
    Fe = Fe * 2 * max(abs(F0));
    c = c + 0.5*xPhys(i)^penal*Ue(:)'*Ke*Ue(:);    
    if ~isempty(intersect(element(i,1),loadedele(:,1)))
        dc(i) = 0.5*xPhys(i)^(penal-1)*(Ue(:)'*Ke*Ue(:) + Fe*Ue(:));
    else
        dc(i) = 0.5*xPhys(i)^(penal-1)*Ue(:)'*Ke*Ue(:);
    end
end   
            
% plot amplified displacement 
new_node=node;
amplify=100;
new_node(:,2)=new_node(:,2)+u_displacement*amplify;
new_node(:,3)=new_node(:,3)+v_displacement*amplify;
new_node(:,4)=new_node(:,4)+w_displacement*amplify;
elements=element;
nodes=new_node;
xmin = min(node(:,2)); xmax = max(node(:,2));
ymin = min(node(:,3)); ymax = max(node(:,3));
zmin = min(node(:,4)); zmax = max(node(:,4));
figure;plot_mesh;title('deformed mesh with displacements be amplified by 100 times '); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
% saveas(gcf,'deformed mesh','png')
%% post processing(calculate strain and stress field)
disp('post processing')
[strain stress] = post_processing_3D(node,element,displacement,C,xPhys,penal,E0,Emin);
von_stress = sqrt(0.5*((stress(:,1)-stress(:,2)).^2 + (stress(:,2)-stress(:,3)).^2 + ...
                  (stress(:,3)-stress(:,1)).^2 + 6*(stress(:,4).^2 + stress(:,5).^2 + stress(:,6).^2)));
%% plot field 
disp('strain and stress field')
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,u_displacement);title('u'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,v_displacement);title('v'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,w_displacement);title('w'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');

figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,strain(:,1));title('\epsilon_x_x'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,strain(:,2));title('\epsilon_y_y'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,strain(:,3));title('\epsilon_z_z'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,strain(:,4));title('\epsilon_x_y'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,strain(:,5));title('\epsilon_y_z'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse'); 
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,strain(:,6));title('\epsilon_x_z'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse'); 

figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,stress(:,1));title('\sigma_x_x'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,stress(:,2));title('\sigma_y_y'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,stress(:,3));title('\sigma_z_z'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse'); 
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,stress(:,4));title('\sigma_x_y'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,stress(:,5));title('\sigma_y_z'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,stress(:,6));title('\sigma_x_z'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,von_stress);title({'$\overline{\sigma}$'},'interpreter','latex'); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');
% 
% color = zeros(num_node,1);
% count = zeros(num_node,1);
% for i = 1:num_element
%     node_num = element(i,2:9);
%     color(node_num) = color(node_num) + thickness(element(i,1));
%     count(node_num) = count(node_num) + 1;
% end
% color = color./count;
% figure('windowstyle','docked');FEA_field_plot_3D(new_node(:,2:4),element,color*40); xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]); set(gca,'YDir','reverse');