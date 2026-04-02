function f = rhs_isoth_euler(t,y,para,dx)
num_f_aux = euler_flux(y,para);
f = -[ num_f_aux(1:end/2) -  circshift(num_f_aux(1:end/2),1); num_f_aux(end/2 + 1 : end) -  circshift(num_f_aux(1 + end/2 : end),1) ] / dx;
end