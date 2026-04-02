%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linear advection with speed 1                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem,Rminus] = PD_isothEuler(t,y,varargin)
getpara(); % getting phys_flux, num_flux, Nx, dx, bc


vecminus = -min(0,num_flux(y,para)) / dx;
vecplus = max(0,num_flux(y,para)) / dx;

switch bc
    case 'periodic'
        D =circshift(spdiags(vecplus(1:Nx),0,Nx,Nx), [0 1]) + circshift(spdiags(vecminus(1:Nx),0,Nx,Nx), [1 0]);
end
% P = spdiags(y_vec,[-1,Nx-1], Nx, Nx);
%D = DD(y);
P = D';



if nargout > 2
    switch bc
        case 'periodic'
             f = @(t,y) rhs_isoth_euler(t,y,para,dx);
            % f = @(t,y) sum((DD(y))' - DD(y), 2);
    end
end
% if max(abs(sum(P-D,2)-f(t,y))) > 1e-14
 %    disp(max(abs(sum(P-D,2)-f(t,y))));
 %end
if nargout > 3 %% Need to be done
 S=zeros(length(y));
  
  
   r = @(y) zeros(length(y),1);

end

if nargout > 5 
Rrem = zeros(length(y)/2,1); % remainder
Rminus = Rrem;
end