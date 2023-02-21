function val = gradtensor_leg(x,m,h,order)

val = NaN(size(x,1),2);


%val(:,1)=shift_legP(x(:,1),m(1),h(1),order(1)).*shift_leg(x(:,2),m(2),h(2),order(2));

%val(:,2)=shift_leg(x(:,1),m(1),h(1),order(1)).*shift_legP(x(:,2),m(2),h(2),order(2));

val(:,1)=shift_leg_derivative(x(:,1),m(1),h(1),order(1),1).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),0);

val(:,2)=shift_leg_derivative(x(:,1),m(1),h(1),order(1),0).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),1);


end