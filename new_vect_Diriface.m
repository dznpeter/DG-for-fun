function z = new_vect_Diriface(node, vertice,T1,index1,swap1,n ,Po ,Lege_ind, sigma, g,a)


dim_elem =size(Lege_ind,1); % number of basis for each element.

z = zeros(dim_elem,1);  % initialize the local vector

% generating quadrature points and weights

%[weights,ref_Qpoints] = Golub_Welsch(ceil(Po/2));

[weights,ref_Qpoints] = quad_GL(ceil((Po+3)*0.5));

[quad_x1,quad_y1] = quadpreference(ref_Qpoints,index1,swap1);

% Change to reference quads

 [quad_x1,quad_y1] = rstoab(quad_x1,quad_y1);
 

%change the quadrature nodes from reference domain to physical domain.

mid= sum(node)./2;   tanvec = 0.5* (node(2,:)-node(1,:));

C = kron(mid,ones(size(ref_Qpoints,1),1));

P_Qpoints = kron(ref_Qpoints,tanvec) + C;   De = norm((node(2,:)-node(1,:))).*0.5;

% penalty term

a_center = a(mid); weights  =  weights.*De;

%lamda = eigs(reshape(a_center,2,2),1);


temp = [n(1).^2,n(1).*n(2),n(1).*n(2),n(2).^2 ];

lamda = dot(a_center,temp);


%correct measure

h_per = max(abs(sum((vertice-kron(node(1,:),ones(size(vertice,1),1))).*kron(n,ones(size(vertice,1),1)),2)));

sigma = sigma.*lamda/h_per; 

 % data for quadrature
    
 
 
 % inverse of B

B1 = 0.5*[T1(2,:)-T1(1,:); T1(3,:)-T1(1,:)];


inv_B1 = inv(B1);
 
 g_val = g(P_Qpoints);   n_vec =kron(n,ones(size(ref_Qpoints,1),1));
 
 
 a_val = a(P_Qpoints);
 
first = zeros(dim_elem,1);


% i is row which is basis of u and j is v. 

   
    for j = 1:dim_elem
    
      % first term g_D(grad{v} cdot n - sigma v)
      
       coe = [(a_val(:,1).*n_vec(:,1)+a_val(:,2).*n_vec(:,2)) , (a_val(:,3).*n_vec(:,1)+a_val(:,4).*n_vec(:,2))]; 
      
        
       
       
       gradv1 = GradSimplex2DP(quad_x1,quad_y1,Lege_ind(j,1),Lege_ind(j,2),inv_B1);
        
       
       V1x = gradv1(:,1);    V1y = gradv1(:,2);
      
      
       
      
        V1 = Simplex2DP(quad_x1,quad_y1,Lege_ind(j,1),Lege_ind(j,2)); 
       
       t = g_val.*(coe(:,1).*V1x +coe(:,2).*V1y -sigma*V1);
         
       
        
        first(j) = dot((t),weights);
         
        
    end



z = z + first;


end
    