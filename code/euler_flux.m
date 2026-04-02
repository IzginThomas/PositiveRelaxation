function fl = euler_flux(u,c)
Nx = length(u)/2;
fl = zeros(size(u));
u1 = u(1:Nx);
u2 = u(Nx+1:2*Nx);

fl(1 : Nx) = log_mean(u1,@(u) u) .* arith_mean(u2./u1);
fl(Nx + 1 : end) = fl(1:Nx) .* (arith_mean(u2./u1)) + arith_mean(c^2*u1);
end