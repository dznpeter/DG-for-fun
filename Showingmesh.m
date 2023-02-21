%% Show the mesh

%clear all;


%load('16 polygonal Elements.mat')


%figure;
hold on;
for i = 1: NT
 %% plot element
 
b = Node(Elem{i,1},:);
plot([b(:,1) ;b(1,1)],[b(:,2); b(1,2)],'k-','LineWidth',1)    

%barycenter = sum(Node(Elem{i,1},:))./size(Node(Elem{i,1},:),1);

%text(barycenter(1),barycenter(2),num2str(i)); 
end; 



for j = 1 : size(Node,1)
    
   
    text(Node(j,1),Node(j,2),num2str(j)); 
    
end





%{

edge = DiribdEdge;

%edge = outflowEdge;

for t = 1:size(edge,1) 

 t1 = [Node(edge(t,1),:); Node(edge(t,2),:)];

plot(t1(:,1), t1(:,2),'*r-','LineWidth',1,'MarkerSize',5)
 
xlim([0 ,1]); ylim([0 ,1]);

barycenter = (Node(edge(t,1),:) + Node(edge(t,2),:)).*0.5;

text(barycenter(1),barycenter(2),num2str(t)); 

end
%}
