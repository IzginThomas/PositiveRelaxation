function [gamma, y_gamma, eta_ygamma] = positive_relax_cons_2nd_order_step(y_gamma,yn, gamma, eta,  eta_ygamma,etaold,eta_prime,y_func,M_func)
aux_vec = (y_gamma - yn) / gamma;
y_gammaPrime = M_func(gamma) \ aux_vec;
rPrime = (eta_prime(y_gamma))' * y_gammaPrime;
gamma = gamma - (eta_ygamma - etaold) / (rPrime + realmin);
if gamma >0
    y_gamma = y_func(gamma);
    eta_ygamma = eta(y_gamma);
end