% plotting edge


figure; hold on;


edge = bdEdge;


%edge = DiribdEdge;

%edge = outflowEdge;

for t = 1:size(edge,1) 

 t1 = [node(edge(t,1),:); node(edge(t,2),:)];

plot(t1(:,1), t1(:,2),'*r-','LineWidth',1,'MarkerSize',5)
 
xlim([0 ,1]); ylim([0 ,1]);

end

%{


NT =16

for i = 1: NT
 %% plot element
 
b = Node(Elem{i,1},:);
plot([b(:,1) ;b(1,1)],[b(:,2); b(1,2)],'k-','LineWidth',1)    
 %% bounding box 
 
a=[Elem{i,2}(1,:);[Elem{i,2}(1,1) Elem{i,2}(2,2)]  ;Elem{i,2}(2,:);  [Elem{i,2}(2,1) Elem{i,2}(1,2)] ] ;
plot([a(:,1) ;a(1,1)],[a(:,2); a(1,2)],'o-','LineWidth',1,'MarkerSize',5)



end

%}