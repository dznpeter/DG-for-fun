%% plotting errors for the DG for h refinement

%{
clear all;



Polynomial_degree = [3];      penalty=10;

NO_elem = [12,48,192,768,3072,12288, 49152 ,196608]';

%errL_2 =NaN(length(NO_elem),1); dof = NaN(length(NO_elem),1);

errL_2 =NaN(length(NO_elem),1); dof = NaN(length(NO_elem),1);


for i=1 :length(NO_elem)


load(['Error ' num2str(NO_elem(i)) ' rectangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])


errL_2(i)=L2_err;         dof(i) = dim_FEM;

end

logerrL_2 = abs(log(errL_2)); 

slope = abs((logerrL_2(2:end)-logerrL_2(1:end-1))./(log(dof(2:end).^0.5)-log(dof(1:end-1).^0.5)));

slope

%slope = sort(slope);

slope = mean(slope(end-2:end));
%slope = max(slope);


loglog(dof.^0.5,errL_2,'-h','LineWidth',1,'MarkerSize',8);

%plot(log((dof(end:-1:1)).^0.5), log(errL_2(end:-1:1).^0.5),'-h','LineWidth',1,'MarkerSize',8);

xlabel('Dof^{1/2}','FontSize',18);

ylabel('||u-u_h ||_{L^2{(\Omega)}}','FontSize',18);

%ylabel('||| u-u_{h}|||_{DG}','FontSize',18);

legend(['P' num2str(Polydegree) ' basis with slope ' num2str(slope)]);

title(['P' num2str(Polynomial_degree) ' basis under h refinement'],'FontSize',18)


%}



figure


%% plotting errors for the DG for p refinement

Polynomial_degree = 1:8;    penalty=10;

NO_elem = [128];   Basis_type = 'P'; % P, Q

errL_2 =NaN(length(Polynomial_degree),1); dof = NaN(length(Polynomial_degree),1);


for i=1 :length(Polynomial_degree)
    
    

load(['Error ' num2str(NO_elem) ' triangle Elements penalty ' num2str(penalty) ' ' num2str(Basis_type) num2str(Polynomial_degree(i)) ' basis.mat'])


errL_2(i)=H1_err;         dof(i) = dim_FEM;

   

end

%{

%semilogy(Polynomial_degree,errL_2,'-h','LineWidth',1,'MarkerSize',8);

semilogy(dof.^0.5,errL_2,'-h','LineWidth',1,'MarkerSize',8);

legend([num2str(NO_elem) ' rect penalty ' num2str(penalty)])

%xlabel('Polynomial order of basis','FontSize',18);

xlabel('Dof^{1/2}','FontSize',18);

ylabel('||u-u_h ||_{L^2{(\Omega)}}','FontSize',18);

%ylabel('||| u-u_{h}|||_{DG}','FontSize',18);

title('P refinement','FontSize',18)
%}


logerrL_2 = abs(log(errL_2)); 

Polynomial_degree =Polynomial_degree';

slope = ((logerrL_2(2:end)-logerrL_2(1:end-1))./(log(Polynomial_degree(2:end))-log(Polynomial_degree(1:end-1))));

%slope = ((logerrL_2(2:end)-logerrL_2(1:end-1))./(log(dof(2:end).^0.5)-log(dof(1:end-1).^0.5)));


slope

%slope = sort(slope);

slope = mean(slope(end-3:end));

%slope = median(slope);

%slope = max(slope);


%semilogy(Polynomial_degree,errL_2,'-s','LineWidth',1,'MarkerSize',8);

%loglog(dof.^0.5,errL_2,'-s','LineWidth',1,'MarkerSize',8);

switch Basis_type 
    
    case 'P' 

semilogy(Polynomial_degree,errL_2,'-sr','LineWidth',1,'MarkerSize',8);

%loglog(dof.^0.5,errL_2,'-sr','LineWidth',1,'MarkerSize',8);


    case 'Q' 

loglog(Polynomial_degree,errL_2,'-sb','LineWidth',1,'MarkerSize',8);

%loglog(dof.^0.5,errL_2,'-sb','LineWidth',1,'MarkerSize',8);


end


legend([num2str(NO_elem) ' Rect ' num2str(Basis_type) ' basis slope ' num2str(slope)])

xlabel('Polynomial order of basis','FontSize',18); 

%xlabel('Dof^{1/2}','FontSize',18);

%ylabel('||u-u_h ||_{L^2{(\Omega)}}','FontSize',18);

ylabel('||| u-u_{h}|||_{DG}','FontSize',18); 

title('DG error under p-refinement','FontSize',18)

set(gca,'FontSize',18);
