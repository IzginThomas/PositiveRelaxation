function f = log_mean(u,phys_flux)
pfluxu = phys_flux(u);
uplus = circshift(u,-1);
fuplus = circshift(pfluxu,-1);
f = (fuplus - pfluxu) ./ (log(uplus+realmin) - log(u+realmin)); % =NaN iff uplus = u
 f(isnan(f)) = pfluxu(isnan(f)); % for consistency
ind = find(abs(uplus - u) < 1e-14);
f(ind) = (pfluxu(ind) + fuplus(ind))/2; 
end