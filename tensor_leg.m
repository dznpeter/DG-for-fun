function val = tensor_leg(x,m,h,order)

val=shift_leg_derivative(x(:,1),m(1),h(1),order(1),0).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),0);

end