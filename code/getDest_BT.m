function DD = getDest_BT(y)
test_BT(); % getting d, A, Nx and dx here
vecminus = -min(0,num_flux(y)) / dx;
vecplus = max(0,num_flux(y)) / dx;
d = length(A);
DD = zeros(d*Nx,d*Nx);
for k=1:d
    ids = (k-1)*Nx+1:k*Nx;
    DD(ids,ids) = circshift(spdiags(vecplus(ids),0,Nx,Nx), [0 1]) + circshift(spdiags(vecminus(ids),0,Nx,Nx), [1 0]);
    % ids_inner = (k-1)*Nx + 2: k*Nx-1;
  %  DD(ids_inner,ids_inner) = D(2:end-1,2:end-1); % no-flux Boundary conditions 
end