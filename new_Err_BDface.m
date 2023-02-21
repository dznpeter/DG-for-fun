function DGPart_3 = new_Err_BDface(node, vertice,index1,swap1,coef, n ,Po ,Lege_ind, u_exact,a , sigma)



dim_elem =size(Lege_ind,1); % number of basis for each element.

DGPart_3=0;

% generating quadrature points and weights

%[weights,ref_Qpoints] = Golub_Welsch(ceil(Po/2));



[weights,ref_Qpoints] = quad_GL(ceil((Po+3)*0.5));

[quad_x1,quad_y1] = quadpreference(ref_Qpoints,index1,swap1);

% Change to reference quads

 [quad_x1,quad_y1] = rstoab(quad_x1,quad_y1);

%change the quadrature nodes from reference domain to physical domain.

mid= sum(node)./2;   tanvec = 0.5* (node(2,:)-node(1,:));

C = kron(mid,ones(size(ref_Qpoints,1),1));

P_Qpoints = kron(ref_Qpoints,tanvec) + C;  De = norm((node(2,:)-node(1,:))).*0.5;

 weights  =  weights.*De;
 
% penalty term

a_center = a(mid);



temp = [n(1).^2,n(1).*n(2),n(1).*n(2),n(2).^2 ];

lamda = dot(a_center,temp);

%correct measure

h_per = max(abs(sum((vertice-kron(node(1,:),ones(size(vertice,1),1))).*kron(n,ones(size(vertice,1),1)),2)));

sigma = sigma.*lamda/h_per; 

 % data for quadrature, function value b and normal vector n_vec
    

      
u_val = u_exact(P_Qpoints); 
 
    % construct the matrix for all the local basis function

    
   P = zeros(size(P_Qpoints,1) ,dim_elem);
   
  
    for i =1:dim_elem
        
        
        P(:,i) = Simplex2DP(quad_x1,quad_y1,Lege_ind(i,1),Lege_ind(i,2));
        
       
        
    end
  
    u_DG_val = P*coef;   %DG solution;
    
    
    
    
    
     
     % Part 3 DG norm error  
     % sigma*(u-u_DG)^2-2a(grad_u - grad_U_DG)\cdot n* (u-u_DG) 
     
     t = sigma.*(u_val-u_DG_val).^2;
     

     DGPart_3= DGPart_3+ dot(t,weights);

end