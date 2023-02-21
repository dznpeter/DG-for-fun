%%Initialize the Polymesh generator to construct the information about the
%%underlying FEM space and the geometry of the computational domain.



clear all; 

%% Generate the triangles for L-shaped domain.

%% triangle domain

[node,elem] = squaremesh([-1,1,-1,1],1);

for k=1:5

[node,elem]=uniformrefine(node,elem);

end


figure(1); showmesh(node,elem);


Node=node;  Element=elem;  

clear  node elem;

%% Generate bounding box for each polygonal element

% move to the correct place.

n = 8;   Node = round(10^n.*Node)/10^n;




%  and  finding all the edges. Label element to edges 



N = size(Node,1); NT = size(Element,1); 

Elem = cell(NT,2);  P = NaN(NT,2);

totalEdge = NaN(3*NT,2);  elemperedge = 3*ones(NT,1);




for i =1: NT
    
    % barycenter
    
    P(i,:) = sum(Node(Element(i,:)',:))./size(Node(Element(i,:)',:),1);
    
    % calculate all the bounding box
    
    Elem{i,1}= Element(i,:)';
   
    
    bdbox_x = sort(Node(Elem{i,1},1));
    
    bdbox_y = sort(Node(Elem{i,1},2));
    
    Elem{i,2} = [bdbox_x([1,end]) ,  bdbox_y([1,end])]; %bounding box is [x_min y_min
                                                        %                 x_max y_max]
                                                                                                        
    % take out all the edges  
    
    totalEdge((i-1)*3+1:3*i,:) = [[Elem{i,1}, [Elem{i,1}(2:end); Elem{i,1}(1) ] ] ] ;                                                       
                    

end

elemtotaledge = cumsum(elemperedge);


%{

%ploting the bounding box to cover  element!

if NT <= 300

figure;
hold on;
for i = 1: NT
 %% plot element
 
b = Node(Elem{i,1},:);
plot([b(:,1) ;b(1,1)],[b(:,2); b(1,2)],'k-','LineWidth',1)    
 %% bounding box 
 
a=[Elem{i,2}(1,:);[Elem{i,2}(1,1) Elem{i,2}(2,2)]  ;Elem{i,2}(2,:);  [Elem{i,2}(2,1) Elem{i,2}(1,2)] ] ;
plot([a(:,1) ;a(1,1)],[a(:,2); a(1,2)],'o-','LineWidth',1,'MarkerSize',5)

text(P(i,1),P(i,2), num2str(i));

end; 


end;

 %}

%% Classify all the edges

% totalEdge is all the edges from summation of each element
% edge is the all the edges from triagulation 
% bdEdge is the boundary edge and intEdge is the interior edge

totalEdge = sort(totalEdge, 2);

[i , j ,s ] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));

edge = [j,i]; bdEdge = [j(s==1), i(s==1)]; intEdge = [j(s==2), i(s==2)];

%% The relation between edge to element and

% internal edge is shared by 2 elements and boundary edge is only used by 1
% edge

int_edge2elem = NaN(size(intEdge,1),2); bd_edge2elem = NaN(size(bdEdge,1),1);

for i = 1: size(intEdge,1)
   
    edge2elem = elem_share_edge(intEdge(i,:),totalEdge,elemtotaledge);
    
    int_edge2elem(i,:) =   edge2elem';
    
end

for i = 1: size(bdEdge,1)
   
    edge2elem = elem_share_edge(bdEdge(i,:),totalEdge,elemtotaledge);
    
    bd_edge2elem(i) =   edge2elem;
    
end


%% Plot for bdedge


% figure; hold on;
% 
% for i = 1 : size(bdEdge,1) 
%     
%    
%    
%   t=Node(bdEdge(i,:),:);
%    
%   plot(t(:,1), t(:,2),'o-','LineWidth',1,'MarkerSize',5)
%        
% end



%% Outward Normal vectors for all boundary edges

% For bdedge, only one normal vector
% For internal edge there two normnal vector, but we only use one.

% tangert vector of the edge%%%%%%
 bdtan_vec = Node(bdEdge(:,1),:) - Node(bdEdge(:,2),:);
 
inttan_vec = Node(intEdge(:,1),:) - Node(intEdge(:,2),:);
 
  % normal vector of the edge%%%%%%
 bdnorvec = [bdtan_vec(:,2)./sqrt(bdtan_vec(:,1).^2 +bdtan_vec(:,2).^2) , -bdtan_vec(:,1)./sqrt(bdtan_vec(:,1).^2 +bdtan_vec(:,2).^2)];

intnorvec = [inttan_vec(:,2)./sqrt(inttan_vec(:,1).^2 +inttan_vec(:,2).^2) , -inttan_vec(:,1)./sqrt(inttan_vec(:,1).^2 +inttan_vec(:,2).^2)];
 
% outward normal vector of the edge%%%%%%
  bdoutward = Node(bdEdge(:,1),:)-P(bd_edge2elem(:),:);
  
  intoutward = Node(intEdge(:,1),:)-P(int_edge2elem(:,1),:); % the first element
 
  
 
 bdindex =  max(sum(bdnorvec.*bdoutward,2),0);
 
 intindex =  max(sum(intnorvec.*intoutward,2),0);
 
 [i ,j ,s] = find(bdindex==0);   [m ,n ,k] = find(intindex==0);
 
 bdnorvec(i,:) = - bdnorvec(i,:);
  
 intnorvec(m,:) = - intnorvec(m,:);

 

 %% Data for Oliver!

% elements = Element; vertices = Node; boundaryVertices=unique(bdEdge(:));
% 
% 
% savefile = [ num2str(NT) ' rectangle Elements for Oliver.mat'];
% 
% save(savefile, 'elements','vertices','boundaryVertices','-v7.3');

 


%save all the useful data

%savefile = [ num2str(NT) ' radial triangle Elements.mat'];

savefile = [ num2str(NT) ' triangle Elements.mat'];
save(savefile, 'Elem','Node','N','NT','bdEdge','bd_edge2elem','intEdge','int_edge2elem','bdnorvec','intnorvec','-v7.3');

%%%%%%%%%%%%%%%%%% The intialization of Mesh is finished %%%%%%%%%%%%%%%%%%%%%

