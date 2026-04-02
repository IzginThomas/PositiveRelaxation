function mean = arith_mean(u,varargin)
mean = (circshift(u,-1) + u)/2;
end