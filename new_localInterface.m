function  z = new_localInterface(node, vertice,T1,T2,index1,swap1, index2,swap2, n ,Po ,Lege_ind, sigma,a)


dim_elem =size(Lege_ind,1); % number of basis for each element.

% generating quadrature points and weights

%[weights,ref_Qpoints] = Golub_Welsch(ceil(Po/2));

[weights,ref_Qpoints] = quad_GL(ceil((Po+1)*0.5));

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

a_val = a(P_Qpoints);  weights  =  weights.*De;

% penalty term

a_center = a(mid);

%lamda = eigs(reshape(a_center,2,2),1);


temp = [n(1).^2,n(1).*n(2),n(1).*n(2),n(2).^2 ];

lamda = dot(a_center,temp);

%correct measure

h_per = max(abs(sum((vertice-kron(node(1,:),ones(size(vertice,1),1))).*kron(n,ones(size(vertice,1),1)),2)));

sigma = sigma.*lamda/h_per; %

 % data for quadrature, function value b and normal vector n_vec
    
n_vec =kron(n,ones(size(ref_Qpoints,1),1));


% inverse of B

B1 = 0.5*[T1(2,:)-T1(1,:); T1(3,:)-T1(1,:)];  
    
B2 = 0.5*[T2(2,:)-T2(1,:); T2(3,:)-T2(1,:)];    


inv_B1 = inv(B1); inv_B2 = inv(B2);

    
 % for splitting the caculation from 2 parts

auxilary_grad_1 = zeros(dim_elem,dim_elem); auxilary_grad_2 = zeros(dim_elem,dim_elem);

auxilary_grad_3 = zeros(dim_elem,dim_elem); auxilary_grad_4 = zeros(dim_elem,dim_elem);

auxilary_sigma_1 = zeros(dim_elem,dim_elem); auxilary_sigma_2 = zeros(dim_elem,dim_elem);

auxilary_sigma_3 = zeros(dim_elem,dim_elem); auxilary_sigma_4 = zeros(dim_elem,dim_elem);

% i is row which is basis of u and j is v. 1 is kappa and 2 is neighbour




for i = 1:dim_elem
   

    
    %%symetric term
    
    for j=i:dim_elem
    
    
        
% Define some values      
        
        coe = [(a_val(:,1).*n_vec(:,1)+a_val(:,2).*n_vec(:,2)) , (a_val(:,3).*n_vec(:,1)+a_val(:,4).*n_vec(:,2))];
        
      
        
        %gradu1 =  gradtensor_leg(P_Qpoints,m1,h1,Lege_ind(i,:));
        
        gradu1 = GradSimplex2DP(quad_x1,quad_y1,Lege_ind(i,1),Lege_ind(i,2),inv_B1);
        
        %gradu2 =  gradtensor_leg(P_Qpoints,m2,h2,Lege_ind(i,:));
        
        gradu2 = GradSimplex2DP(quad_x2,quad_y2,Lege_ind(i,1),Lege_ind(i,2),inv_B2);
        
        %gradv1 =  gradtensor_leg(P_Qpoints,m1,h1,Lege_ind(j,:));
        
        gradv1 = GradSimplex2DP(quad_x1,quad_y1,Lege_ind(j,1),Lege_ind(j,2),inv_B1);
        
        %gradv2 =  gradtensor_leg(P_Qpoints,m2,h2,Lege_ind(j,:));
        
        gradv2 = GradSimplex2DP(quad_x2,quad_y2,Lege_ind(j,1),Lege_ind(j,2),inv_B2);
        
        U1x = gradu1(:,1);    U1y = gradu1(:,2);
        
        U2x = gradu2(:,1);    U2y = gradu2(:,2);
           
        V1x = gradv1(:,1);    V1y = gradv1(:,2);
        
        V2x = gradv2(:,1);    V2y = gradv2(:,2);
                 
        U1 = Simplex2DP(quad_x1,quad_y1,Lege_ind(i,1),Lege_ind(i,2));
      
        U2 = Simplex2DP(quad_x2,quad_y2,Lege_ind(i,1),Lege_ind(i,2));
         
        V1 = Simplex2DP(quad_x1,quad_y1,Lege_ind(j,1),Lege_ind(j,2));
        
        V2 = Simplex2DP(quad_x2,quad_y2,Lege_ind(j,1),Lege_ind(j,2));
        
  % first term a/2*[(grad(u^+)\cdot n)*(v^+)+(grad(v^+)\cdot n)*(u^+)]
  %            - sigma*(u^+)*(v^+)      
        
        
        
 t1 = 0.5.*(coe(:,1).*(U1x.*V1+V1x.*U1) ...
                 +coe(:,2).*(U1y.*V1+V1y.*U1));
        
       
auxilary_grad_1(j,i) = dot((t1),weights);


s1 = - (sigma.*U1.*V1) ;

auxilary_sigma_1(j,i) = dot((s1),weights);             
       
        
        % second term 1/2*[(grad(u^-)\cdot n)*(v^+)-(grad(v^+)\cdot n)*(u^-)]
        %            + sigma*(u^-)*(v^+)
        
               
        
t2 = 0.5.*(coe(:,1).*(U2x.*V1-V1x.*U2) ...
               +coe(:,2).*(U2y.*V1-V1y.*U2));
        
       
auxilary_grad_2(j,i) = dot((t2),weights);


s2 = sigma.*U2.*V1 ;

auxilary_sigma_2(j,i) = dot((s2),weights);             
             
        
        % third term 1/2*[-(grad(u^+)\cdot n)*(v^-)+(grad(v^-)\cdot n)*(u^+)]
        %            + sigma*(u^+)*(v^-)
        

        
t3 = 0.5.*(coe(:,1).*(-U1x.*V2+V2x.*U1) ...
          +coe(:,2).*(-U1y.*V2+V2y.*U1));
        
       
auxilary_grad_3(j,i) = dot((t3),weights);


s3 = sigma.*U1.*V2;

auxilary_sigma_3(j,i) = dot((s3),weights);          
         
         
         % fourth term -1/2*[(grad(u^-)\cdot n)*(v^-)+(grad(v^-)\cdot n)*(u^-)]
         %           - sigma*(u^-)*(v^-)
        

        
  t4 = -0.5.*(coe(:,1).*(U2x.*V2+V2x.*U2) ...
                +coe(:,2).*(U2y.*V2+V2y.*U2));
        
       
auxilary_grad_4(j,i) = dot((t4),weights);


s4 = - (sigma.*U2.*V2) ;

auxilary_sigma_4(j,i) = dot((s4),weights);         
        
    end
end

% sigma is always symmetric but grad is sometimes symmetric but sometimes
% skewsymmetric


local_1 = auxilary_grad_1+ (tril(auxilary_grad_1,-1))'...
         +auxilary_sigma_1+(tril(auxilary_sigma_1,-1))';

     
% skewsymmetric

local_2 = auxilary_grad_2+ (tril(auxilary_grad_3,-1))'...
         +auxilary_sigma_2+(tril(auxilary_sigma_3,-1))';

% skewsymmetric
     
local_3 = auxilary_grad_3+ (tril(auxilary_grad_2,-1))'...
         +auxilary_sigma_3+(tril(auxilary_sigma_2,-1))';

local_4 = auxilary_grad_4+ (tril(auxilary_grad_4,-1))'...
         +auxilary_sigma_4+(tril(auxilary_sigma_4,-1))';


z = [local_1, local_2; local_3,local_4];

end