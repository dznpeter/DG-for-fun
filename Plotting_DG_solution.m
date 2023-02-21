%%plot DG solution 


Dimension =2;

Polydegree =1; 

Basis_type = 'P';%Q

NT = 2048;

penalty = 10;

load([num2str(NT) ' triangle Elements.mat']);

% true solution 
U_true = u_true(Node);



%%DG Solution 
load(['Error ' num2str(NT) ' triangle Elements penalty ' num2str(penalty) ' P' num2str(Polydegree) ' basis.mat'])


if Polydegree == 0
    
    Lege_ind = [ 0 0];
else
    
    Lege_ind = [0 0];
    
for i = 1 : Polydegree 

Lege_ind = [Lege_ind ;multiindex( 2 , i +2) - 1] ;

end

end


dim_elem = size(Lege_ind,1); %number of basis for each element.


figure 

hold on

for t =1:NT
   
    [elem, BDbox] = Elem{t,:};   node = Node(elem,:);
    
    coef = U((t-1)*dim_elem+1 : 1 :t*dim_elem ,1); % c is coefficeint
       
    % information about the bounding box

    h = (BDbox(2,:)-BDbox(1,:))./2;  

    m = 0.5.*sum(BDbox);

    P = zeros(size(node,1) ,dim_elem);  
       
quad_x = [-1,1,-1]'; quad_y = [-1,-1,1]';
    
    for i =1:dim_elem
        
        P(:,i)= Simplex2DP(quad_x,quad_y,Lege_ind(i,1),Lege_ind(i,2));             
                
        
    end
  
    u_DG_val = P*coef;   %DG solution;

    % plot 
    
    c = (u_DG_val-min(U_true))./max(U_true);
    
    fill3(node(:,1),node(:,2),u_DG_val,c) ; 
end



hold off;view(3);

%title([num2str(NT) ' polygonal meshes with ' num2str(Basis_type) num2str(Polydegree) ' basis'],'FontSize',20);

set(gca,'FontSize',20);set(gca,'FontSize',20);