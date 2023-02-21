%%%%%% The second order elliptic PDE %%%%%%%%%%%

%% for different NO of elements

Number_Elements = [32,128,512,2048,8192];

Polynomial_order =[1:4];

for NO_elem = 1: size(Number_Elements,2)

 %% initializing all the data for DG from VEM
    
load([num2str(Number_Elements(NO_elem)) ' triangle Elements.mat'])

%% for different order of Basis!

for  O_B = 1:size(Polynomial_order,2)
    
disp([num2str(Number_Elements(NO_elem)) ' element with polynomial basis order ' num2str(Polynomial_order(O_B))]);

%%Innitialize the information about the PDE%%%%%%%%%%%

a = @(x)diffusive(x); % Diffusive function

g = @(x)u_true(x); % the Dirichelet boundary condition 

f = @(x)forcing(x); % the forcing founction

u = @(x)u_true(x); % the analytic solution 

grad_u = @(x)grad_u_true(x); % the gradient of analytic solution

Dimension = 2; % 2D space

Polydegree = Polynomial_order(O_B);   % polynomial degree of DGFEM    

Po = 2*Polydegree+1 ;         % quadrature order 

penalty = 10; % Pernalty parameter    

sigma = penalty*Polydegree.^2;   % \sigma = penalty * p^2/h, h is inside the face...

%% Set up different boundary condition

DiribdEdge=bdEdge;   Diribd_edge2elem=bd_edge2elem;  Diribdnorvec = bdnorvec;

%% Genrating Sparse Matrix A, we are in two D geometry, so the number od basis for each lement is nchoose(p+d,d)

dim_elem = nchoosek(Polydegree+Dimension,Dimension);    dim_FEM = nchoosek(Polydegree+Dimension,Dimension) .*NT;

%legendre polynomial order

Lege_ind = Basis_index_generator(Polydegree,Dimension);


%% Contruct the bilinear form with Global matrix..
%% Part 1

% contribution from element


i =zeros(dim_elem.^2,NT ); j =zeros(dim_elem.^2,NT); s =zeros(dim_elem.^2,NT ); 

parfor t =1:NT
          
    [elem, BDbox] = Elem{t,:}; 
    
    local = new_localstiff(Node(elem,:), Po ,Lege_ind,@(x)a(x));
    
   local = local';   
   
   ind = (t-1)*dim_elem+1 : t*dim_elem;   ind = ind'; 
    
   i(:,t) = kron(ind , ones(dim_elem,1)) ;
   
   j(:,t) = kron(ones(dim_elem,1),ind )  ;
   
   s(:,t) = local(:);
    
end

B1 = sparse(i(:),j(:) ,s(:) ,dim_FEM,dim_FEM );

disp('Bilinear form Part 1 Element is over');


 %% part 2 is the comtribution from interface
 


i = zeros(4.*dim_elem.^2,size(intEdge,1) ); j =zeros(4.*dim_elem.^2,size(intEdge,1) ); s =zeros(4.*dim_elem.^2,size(intEdge,1) );
 


parfor t =1:size(intEdge,1)
    
    
    elem_oneface =int_edge2elem(t,:);
    
    tri1 = Elem{elem_oneface(1),1};   tri2 = Elem{elem_oneface(2),1};        
    
    [index1, swap1] = wrap(intEdge(t,:),tri1);
    
    [index2, swap2] = wrap(intEdge(t,:),tri2);
    
     T1 = Node(Elem{elem_oneface(1),1},:);  T2 = Node(Elem{elem_oneface(2),1},:); 
         
     
     vertice = [T1;T2]; 
    
    
    interface = new_localInterface(Node(intEdge(t,:),:),vertice ,T1,T2 ,index1,swap1,index2,swap2,intnorvec(t,:) ,Po ,Lege_ind,sigma,@(x)a(x));
    
    interface = interface';
     
     ind1 = (elem_oneface(1)-1)*dim_elem+1 : elem_oneface(1)*dim_elem;
     
     ind2 = (elem_oneface(2)-1)*dim_elem+1 : elem_oneface(2)*dim_elem;
     
     ind1 = ind1';   ind2 = ind2';
     
     ind = [ind1 ;ind2];
     
     
     i(:,t) = kron(ind , ones(2*dim_elem,1)) ;
   
     j(:,t) = kron(ones(2*dim_elem,1),ind )  ;
   
     s(:,t) = interface(:);
    
end

B2 = sparse(i(:),j(:) ,s(:) ,dim_FEM,dim_FEM );

disp('Bilinear form Part 2 Interface is over');




 %% part 3 is the comtribution from Dirichelte boundary
 

i = zeros(dim_elem.^2,size(DiribdEdge,1) ); j = zeros(dim_elem.^2,size(DiribdEdge,1) );
s = zeros(dim_elem.^2,size(DiribdEdge,1) ); 

parfor t =1:size(DiribdEdge,1)
    
    
    Dirielem_bdface =Diribd_edge2elem(t,:);
    
    
    tri1 = Elem{Dirielem_bdface,1};   %tri2 = Elem{elem_oneface(2),1};
    
    
    
    [index1, swap1] = wrap(bdEdge(t,:),tri1);
    
     T1 = Node(Elem{Dirielem_bdface,1},:);  
    
     vertice = [T1]; 
          
     
    DiriBDface = new_localDiriface(Node(DiribdEdge(t,:),:), vertice,T1,index1,swap1,Diribdnorvec(t,:),Po ,Lege_ind,sigma,@(x)a(x));
    
   
   ind = (Dirielem_bdface-1)*dim_elem+1 : Dirielem_bdface*dim_elem;
   
   ind = ind'; 
    
   i(:,t) = kron(ind , ones(dim_elem,1)) ;
   
   j(:,t) = kron(ones(dim_elem,1),ind )  ;
   
   s(:,t) = DiriBDface(:);
    
end

B3 = sparse(i(:),j(:) ,s(:) ,dim_FEM,dim_FEM );

disp('Bilinear form Part 3 Boundary Face is over');






 
%% The global stiffness matrix is the summation of Part1-3

B = B1-B2-B3 ;

disp('The condition NO of Sparse matrix B is ');

Condition_NO = condest(B)





%% The right handside of the matrix, the linear functional;   
    
%Part 1

% contribution from element

i =zeros(dim_elem,NT ); j =ones(dim_elem,NT ); s =zeros(dim_elem,NT ); 


parfor t =1:NT
   
    [elem, BDbox] = Elem{t,:}; 
    
   
    
    localvec = new_vect_forcing(Node(elem,:),  Po, Lege_ind, @(x)f(x));
    
   
     ind = (t-1)*dim_elem+1 : t*dim_elem;   ind = ind'; 
    
   i(:,t) = ind ;
   
   s(:,t) = localvec;
    
end

L1 = sparse(i(:),j(:) ,s(:)  ,dim_FEM,1 );

disp('Linear form Part 1 Element is over');






%Part 2

% contribution from Dirichelet boundary  


i = zeros(dim_elem,size(DiribdEdge,1) );  j =ones(dim_elem,size(DiribdEdge,1)); 

s = zeros(dim_elem,size(DiribdEdge,1) ); 


parfor t =1:size(DiribdEdge,1)
    
 
    
    Dirielem_Diric =Diribd_edge2elem(t,:);
    
    
    tri1 = Elem{Dirielem_Diric,1};   
    
    
    
    [index1, swap1] = wrap(bdEdge(t,:),tri1);
    
     T1 = Node(Elem{Dirielem_Diric,1},:); 
         
    
     
     vertice = [T1]; 
     
     
    Diriface = new_vect_Diriface(Node(DiribdEdge(t,:),:),vertice,T1,index1,swap1,Diribdnorvec(t,:) ,Po ,Lege_ind,sigma ,@(x)g(x),@(x)a(x));
    
    
   ind = (Dirielem_Diric-1)*dim_elem+1 : Dirielem_Diric*dim_elem;
   
   ind = ind'; 
    
   i(:,t) = ind;
   
   s(:,t) = Diriface;
    
end

L2 = sparse(i(:),j(:) ,s(:) ,dim_FEM,1 );

disp('Linear form Part 2 Dirichelet boundary is over');


%% The global vector of linear functional is the summation of Part1-2

L=L1-L2;

%% The coefficient of global basis is U

U = B\L;


%% Test error in L_2 norm and DG norm 


%Part 1

%Error from element  (L_2 norm)

L2_err = NaN(NT,1);  H1_err = NaN(NT,1);  DG_1=NaN(NT,1);

parfor t =1:NT
   
    [elem] = Elem{t,:};   coef = U((t-1)*dim_elem+1 : 1 :t*dim_elem ,1); % c is coefficeint
   
    
    [Part_1, z,y] = new_Err_elem(Node(elem,:), coef, Po ,Lege_ind , @(x)u(x) ,@(x)grad_u(x) ,@(x)a(x));
         
    
    L2_err(t) = z;
    
    H1_err(t) = y;
    
    
    DG_1(t) = Part_1;
    
end



L2_err = sqrt(sum(L2_err)); 

H1_err = sqrt(sum(H1_err));

DG_1 = sum(DG_1);  



%Part 2 from the interface  on Face 1/2*|b \cdot n|jump(e)^2

DG_2=NaN(size(intEdge,1),1);

parfor t =1:size(intEdge,1)
    
    elem_oneface =int_edge2elem(t,:);
    
    
    tri1 = Elem{elem_oneface(1),1};   tri2 = Elem{elem_oneface(2),1};
    
    [index1, swap1] = wrap(intEdge(t,:),tri1);
    
    
    [index2, swap2] = wrap(intEdge(t,:),tri2);
    
    
    coef1 = U((elem_oneface(1)-1)*dim_elem+1 : 1 :elem_oneface(1)*dim_elem ,1);
    
    coef2 = U((elem_oneface(2)-1)*dim_elem+1 : 1 :elem_oneface(2)*dim_elem ,1);        
    
    vertice = [Node(Elem{elem_oneface(1),1},:);Node(Elem{elem_oneface(2),1},:)]; 
    
    Part_2 = new_Err_interface(Node(intEdge(t,:),:), vertice,index1,swap1, index2,swap2,coef1,coef2,intnorvec(t,:) ,Po ,Lege_ind  ,@(x)a(x),sigma);    
    
    
    DG_2(t) =Part_2;
end



DG_2 = sum(DG_2);


%Part 3 from the interface  on Face 1/2*|b \cdot n|jump(e)^2
DG_3=NaN(size(DiribdEdge,1),1); 

parfor t =1:size(DiribdEdge,1)
    
     Dirielem_bdface =Diribd_edge2elem(t,:); 
     
     tri1 = Elem{Dirielem_bdface(1),1};
     
     [index1, swap1] = wrap(bdEdge(t,:),tri1);
     
     coef = U((Dirielem_bdface-1)*dim_elem+1 : 1 :Dirielem_bdface*dim_elem ,1); % c is coefficeint
             
     vertice = Node(Elem{Dirielem_bdface(1),1},:); 
     
     Part_3 = new_Err_BDface(Node(DiribdEdge(t,:),:), vertice, index1,swap1,coef,Diribdnorvec(t,:) ,Po ,Lege_ind, @(x)u(x),@(x)a(x),sigma);
    
    
    
     DG_3(t) = Part_3;
end



DG_3 = sum(DG_3);


disp('L2 norm error is');

L2_err


disp('H1 semi-norm error is');

H1_err



disp('DG norm error is');

DG_err = sqrt(DG_1+DG_2+DG_3); 

DG_err

savefile = ['Error ' num2str(NT) ' triangle Elements penalty ' num2str(penalty) ' P' num2str(Polydegree) ' basis.mat'];


save(savefile, 'U','L2_err', 'H1_err','DG_err','dim_FEM','Polydegree','Condition_NO','-v7.3');


end

end

