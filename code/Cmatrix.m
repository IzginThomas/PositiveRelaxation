%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 04/08/2019                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = Cmatrix(S,r,y1)

N = length(y1);
C = 0;
aux = r(y1);
for j = 1:size(S,2)
  rj = aux(j);
  
  s_p = max(S(:,j),zeros(N,1));
  s_m = min(S(:,j),zeros(N,1));
  
  w = 1/norm(s_p,1)*s_p;
  
  for i = 1:N
    if s_m(i) ~= 0
      e = zeros(N,1);
      e(i) = 1;
      C = C + rj/(y1(i)+realmin)*s_m(i)*(e - w)*e';
    end
  end
end
