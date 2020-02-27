function [strain stress]=post_processing_3D(node,element,displacement,C,xPhys,penal,E0,Emin)
%% this function is used to calcualte strain and stress field for 3D object at each node position 
 master_xyz=zeros(8,3);
 master_xyz(1,:)=[1 1 -1];master_xyz(2,:)=[-1 1 -1];master_xyz(3,:)=[-1 -1 -1];master_xyz(4,:)=[1 -1 -1];
 master_xyz(5,:)=[1 1 1];master_xyz(6,:)=[-1 1 1];master_xyz(7,:)=[-1 -1 1];master_xyz(8,:)=[1 -1 1];
 
 length_node=size(node,1);
 strain       = zeros(length_node,6);
  stress=zeros(length_node,6);
 strain_index = zeros(length_node,1);

 for i=1:size(element,1) 
 node_num   = element(i,2:9);
 XYZ_element = node(node_num,2:4);
 dis_index  = zeros(24,1);
 dis_index(1:3:24,1)  = 3*node_num-2;
 dis_index(2:3:24,1)  = 3*node_num-1; 
 dis_index(3:3:24,1)  = 3*node_num; 
    for j=1:8
     [B N N_diff J H J_inv x y z] = master_element_linear_elasticity_3D('linear',master_xyz(j,1),master_xyz(j,2),master_xyz(j,3),XYZ_element);    
     strain_element=(H*J_inv*N_diff*displacement(dis_index));
     stress_element=(C*(Emin+xPhys(i).^penal*(E0-Emin))/E0*H*J_inv*N_diff*displacement(dis_index));
     strain(node_num(j),:)=strain(node_num(j),:) + strain_element';
     stress(node_num(j),:)=stress(node_num(j),:) + stress_element'; 
     strain_index(node_num(j),1)=strain_index(node_num(j),1)+1;
    end
 end
 strain=strain./repmat(strain_index,1,6);
 stress=stress./repmat(strain_index,1,6);
end