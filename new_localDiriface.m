function z = new_localDiriface(node, vertice,T1,index1,swap1, n ,Po ,Lege_ind, sigma,a)

dim_elem =size(Lege_ind,1); % number of basis for each element.



%[weights,ref_Qpoints] = Golub_Welsch(ceil(Po/2));

[weights,ref_Qpoints] = quad_GL(ceil((Po+1)*0.5));

[quad_x1,quad_y1] = quadpreference(ref_Qpoints,index1,swap1);

% Change to reference quads

 [quad_x1,quad_y1] = rstoab(quad_x1,quad_y1);

 
%change the quadrature nodes from reference domain to physical domain.

mid= sum(node)./2;   tanvec = 0.5* (node(2,:)-node(1,:));

C = kron(mid,ones(size(ref_Qpoints,1),1));

P_Qpoints = kron(ref_Qpoints,tanvec) + C;   De = norm((node(2,:)-node(1,:))).*0.5;



a_val = a(P_Qpoints);   weights  =  weights.*De;

% penalty term

a_center = a(mid);

%lamda = eigs(reshape(a_center,2,2),1);


temp = [n(1).^2,n(1).*n(2),n(1).*n(2),n(2).^2 ];

lamda = dot(a_center,temp);


%correct measure

h_per = max(abs(sum((vertice-kron(node(1,:),ones(size(vertice,1),1))).*kron(n,ones(size(vertice,1),1)),2)));

sigma = sigma.*lamda/h_per; 

 % data for quadrature
    
 n_vec =kron(n,ones(size(ref_Qpoints,1),1));

 
 
 % inverse of B

B1 = 0.5*[T1(2,:)-T1(1,:); T1(3,:)-T1(1,:)];


inv_B1 = inv(B1); 


 
first = zeros(dim_elem,dim_elem); 


% i is row which is basis of u and j is v. Only half of the matrix is
% calculated



for i = 1:dim_elem
   
   
    
    %%symetric term
    
    for j=i:dim_elem
    
      % % aver{grad{u}}.jump(v)+aver{grad{v}}.jump(u) - sigma.*jump(u)\cdot jump(v);
      
      
      coe = [(a_val(:,1).*n_vec(:,1)+a_val(:,2).*n_vec(:,2)) , (a_val(:,3).*n_vec(:,1)+a_val(:,4).*n_vec(:,2))];
      
       
    
       %gradu1 =  gradtensor_leg(P_Qpoints,m,h,Lege_ind(i,:));  
       
       gradu1 = GradSimplex2DP(quad_x1,quad_y1,Lege_ind(i,1),Lege_ind(i,2),inv_B1);
       
       %gradv1 =  gradtensor_leg(P_Qpoints,m,h,Lege_ind(j,:));
       
       gradv1 = GradSimplex2DP(quad_x1,quad_y1,Lege_ind(j,1),Lege_ind(j,2),inv_B1);
        
       
       U1x = gradu1(:,1);    U1y = gradu1(:,2);
       
       V1x = gradv1(:,1);    V1y = gradv1(:,2);
      
       %U1 = tensor_leg(P_Qpoints,m,h,Lege_ind(i,:));
       
       U1 = Simplex2DP(quad_x1,quad_y1,Lege_ind(i,1),Lege_ind(i,2));
       
       %V1 = tensor_leg(P_Qpoints,m,h,Lege_ind(j,:));
       
       V1 = Simplex2DP(quad_x1,quad_y1,Lege_ind(j,1),Lege_ind(j,2));
       
       
       
       t = (coe(:,1).*(U1x.*V1+V1x.*U1) ...
           +coe(:,2).*(U1y.*V1+V1y.*U1))...
               - sigma.*U1.*V1;
       
       
        
        first(j,i) = dot((t),weights);
        
    end
end


%z =  first;

%%symetric term

z = first + (tril(first,-1))';


end