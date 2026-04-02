function dygamma = get_derivative(yn, ygamma, gamma,dsigmagamma,sigmagamma,Mgamma)
vgamma = ygamma .* dsigmagamma ./ (sigmagamma + realmin);
% aux = 0;
% for j=1:length(b)
%     aux = aux + b(j)*(P{j} - diag((Rminus{j} + sum(D{j},2)))) * vgamma;
% end
rhs = ygamma - yn - gamma * vgamma;
dygamma = (1 / gamma) * (Mgamma \ rhs) + vgamma;
end