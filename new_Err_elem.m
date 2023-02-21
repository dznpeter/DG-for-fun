function [DGPart_1, L2part ,H1part] = new_Err_elem(node, coef ,Po ,Lege_ind,  u_exact, grad_u,a)



dim_elem =size(Lege_ind,1); % number of basis for each element.

L2part=0; H1part=0;  DGPart_1=0;

[w_x,x] = quad_GL(ceil((Po+3)*0.5));  [w_y,y] = quad_GJ1(ceil((Po+3)*0.5));


quad_x = kron(x,ones(size(w_y,1),1)); quad_y = kron(ones(size(w_x,1),1),y);
 
 
weights = kron(w_x,w_y) ;

shiftpoints = [(1+quad_x).*(1-quad_y).*0.5-1, quad_y ];
 
 ref_points = 0.5.*shiftpoints+0.5;  

    
     % For each triangle
    
    %B is the affine matrix
    
    B = 0.5*[node(2,:)-node(1,:); node(3,:)-node(1,:)];
    
     De_tri = abs(det(B));    inv_B = inv(B); weights = weights.*0.5.*De_tri;
    
        
    phy = reference_to_physical_t3 (node([1,2,3],:)', size(ref_points,1), ref_points' );
    
    P_Qpoints  = phy';

     
    % data for quadrature
    
    %%% numerical intergration by Duffy Transfomatio
    
     u_val = u_exact(P_Qpoints);    a_val = a(P_Qpoints);
    
     grad_u_val = grad_u(P_Qpoints);
    
    % construct the matrix for all the local basis function

    
    P = zeros(size(P_Qpoints,1) ,dim_elem);
    
    Px = zeros(size(P_Qpoints,1) ,dim_elem);
    
    Py = zeros(size(P_Qpoints,1) ,dim_elem);
    
    for i =1:dim_elem
        
       % P(:,i)= tensor_leg(P_Qpoints,m,h,Lege_ind(i,:));
        
        
        P(:,i) = Simplex2DP(quad_x,quad_y,Lege_ind(i,1),Lege_ind(i,2));
                        
        t = GradSimplex2DP(quad_x,quad_y,Lege_ind(i,1),Lege_ind(i,2),inv_B);
        
        Px(:,i) = t(:,1); Py(:,i) = t(:,2); 
    end
  
    u_DG_val = P*coef;   %DG solution;
    
    grad_u_DG = [Px*coef , Py*coef];   %gradient of DG

        
      % L_2 norm error  int_\kappa (u - u_DG)^2 dx
    
     t1 = (u_val - u_DG_val).^2;
        
     L2part = L2part+ dot((t1),weights);
     
     % H1 semi norm error  int_\kappa (grad_u - grad_u_DG)^2 dx
     
     t2 = sum((grad_u_val - grad_u_DG).^2,2);
     
     H1part = H1part+ dot((t2),weights);
     
     
     % Part 1 DG L_2 norm error  int_\kappa a(u - u_DG)^2 dx
    
     
     grad1 = grad_u_val - grad_u_DG;    grad2 = grad_u_val - grad_u_DG;
        
     grad = [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2),...
            grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2)];
     
     
     t3 = sum((grad).*a_val,2);
     
     
     DGPart_1 = DGPart_1+ dot((t3),weights);
    
end




