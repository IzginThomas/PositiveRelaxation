function f = harm_mean(u,phys_flux)
uplus = circshift(u,-1);
f = 2.*u.*uplus./(u + uplus);
end