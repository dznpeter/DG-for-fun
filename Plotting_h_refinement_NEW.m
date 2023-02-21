
NO_Elements = [75,300,1200,4800 ]; %, 19200,76800];


color = {'b','r','g','k'};

L2=figure;


%% plotting errors for the L2 norm for p refinement

%%P basis


Polynomial_degree = 3;    penalty=10;

NO_elem = NO_Elements;   

err_P = NaN(length(NO_elem),1); dof_P = NaN(length(NO_elem),1);


for i=1 :length(NO_elem)

%load(['Error ' num2str(NO_elem(i)) ' rectangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
    
    
load(['Error ' num2str(NO_elem(i)) ' radial rectangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
        
    
%load(['Error ' num2str(NO_elem(i)) ' polygonal Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
    
    
%load(['Error ' num2str(NO_elem(i)) ' radial polygonal Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])

err_P(i,1)=L2_err;         dof_P(i) = dim_FEM;

       

end



%slope for rect

logerr_P1 = abs(log(err_P(:,1))); 

slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(log(dof_P(2:end))-log(dof_P(1:end-1))));

%slope_P1 = sort(slope_P1); 

slope_Rect =  mean(slope_P1(end));



loglog(dof_P,err_P(:,1),['-s' num2str(color{Polynomial_degree})],'LineWidth',1.5,'MarkerSize',10);

hold on;






legend(['rect P' num2str(Polynomial_degree) ' slope ' num2str(slope_Rect(1))]);

xlabel('Dof','FontSize',20);

ylabel('||u-u_h ||_{L^2(\Omega)}','FontSize',20);

set(gca,'FontSize',15)

title(['L2 norm error under h-refinement penalty ' num2str(penalty) ],'FontSize',20)

%saveas(L2,['radial_rect_Poisson_L_shpae_L2_norm_error_h_refine_penalty_' num2str(penalty) ],'fig');

%print(L2,'-depsc',['radial_rect_Poisson_L_shpae_L2_norm_error_h_refine_penalty_' num2str(penalty) '.eps']);





H1=figure;


%% plotting errors for the L2 norm for p refinement

%%P basis


Polynomial_degree = 2;    penalty=10;

NO_elem = NO_Elements;   

err_P = NaN(length(NO_elem),1); dof_P = NaN(length(NO_elem),1);


for i=1 :length(NO_elem)

%load(['Error ' num2str(NO_elem(i)) ' rectangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
    
    
load(['Error ' num2str(NO_elem(i)) ' radial rectangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
        
    
%load(['Error ' num2str(NO_elem(i)) ' polygonal Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
    
    
%load(['Error ' num2str(NO_elem(i)) ' radial polygonal Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])

err_P(i,1)=H1_err;         dof_P(i) = dim_FEM;

       

end



%slope for rect

logerr_P1 = abs(log(err_P(:,1))); 

slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(log(dof_P(2:end))-log(dof_P(1:end-1))));

%slope_P1 = sort(slope_P1); 

slope_Rect =  mean(slope_P1(end));



loglog(dof_P,err_P(:,1),['-s' num2str(color{Polynomial_degree})],'LineWidth',1.5,'MarkerSize',10);

hold on;






legend(['rect P' num2str(Polynomial_degree) ' slope ' num2str(slope_Rect(1))]);

xlabel('Dof','FontSize',20);

ylabel('|u-u_h |_{H^1(\Omega)}','FontSize',20);

set(gca,'FontSize',20)

title(['H1 seminorm error under h-refinement penalty ' num2str(penalty) ],'FontSize',20)

%saveas(H1,['radial_rect_Poisson_L_shpae_H1_semi_norm_error_h_refine_penalty_' num2str(penalty) ],'fig');

%print(H1,'-depsc',['radial_rect_Poisson_L_shpae_H1_semi_norm_error_h_refine_penalty_' num2str(penalty) '.eps']);


