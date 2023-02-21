function DGPart_2 = new_Err_interface(node,vertice,index1, swap1, index2,swap2,  coef1,coef2,n ,Po ,Lege_ind,a , sigma)



dim_elem =size(Lege_ind,1); % number of basis for each element.


DGPart_2=0;

% generating quadrature points and weights

%[weights,ref_Qpoints] = Golub_Welsch(ceil(Po/2));

[weights,ref_Qpoints] = quad_GL(ceil((Po+3)*0.5));

%index1

[quad_x1,quad_y1] = quadpreference(ref_Qpoints,index1,swap1);
  

[quad_x2,quad_y2] = quadpreference(ref_Qpoints,index2,swap2);

% Change to reference quads

 [quad_x1,quad_y1] = rstoab(quad_x1,quad_y1);

 % the direction of the nodes is the opposite
 
 [quad_x2,quad_y2] = rstoab(quad_x2,quad_y2);
 

%change the quadrature nodes from reference domain to physical domain.

mid= sum(node)./2;   tanvec = 0.5* (node(2,:)-node(1,:));

C = kron(mid,ones(size(ref_Qpoints,1),1));

P_Qpoints = kron(ref_Qpoints,tanvec) + C;  De = norm((node(2,:)-node(1,:))).*0.5;

 weights  =  weights.*De;
 
% penalty term

a_center = a(mid);

%lamda = eigs(reshape(a_center,2,2),1);


temp = [n(1).^2,n(1).*n(2),n(1).*n(2),n(2).^2 ];

lamda = dot(a_center,temp);

%correct measure

h_per = max(abs(sum((vertice-kron(node(1,:),ones(size(vertice,1),1))).*kron(n,ones(size(vertice,1),1)),2)));

sigma = sigma.*lamda/h_per; %


 
    % construct the matrix for all the local basis function

    
    P1 = zeros(size(P_Qpoints,1) ,dim_elem);
    

    
    for i =1:dim_elem
        
        %P1(:,i)= tensor_leg(P_Qpoints,m1,h1,Lege_ind(i,:));
        
        P1(:,i) = Simplex2DP(quad_x1,quad_y1,Lege_ind(i,1),Lege_ind(i,2));
        
 
        
    end
  
    u_DG_val1 = P1*coef1;   %DG solution on Kappa1;
    
   
    
     P2 = zeros(size(P_Qpoints,1) ,dim_elem);

     
    
    for i =1:dim_elem
        
        %P2(:,i)= tensor_leg(P_Qpoints,m2,h2,Lege_ind(i,:));

        P2(:,i) = Simplex2DP(quad_x2,quad_y2,Lege_ind(i,1),Lege_ind(i,2));
        
    end
  
    u_DG_val2 = P2*coef2;   %DG solution on Kappa2;
    
    
    
     
    
     
    
     
     % Part 2 DG norm error   
     % sigma*jump{u_DG}^2+2a(grad_u - aver{grad_U_DG})\cdot jump{u_DG} 
     
     jump_u_DG = (u_DG_val1-u_DG_val2);
     
     
   
     
     
     t = sigma.*jump_u_DG.^2;
     

     DGPart_2= DGPart_2+ dot((t),weights);

end