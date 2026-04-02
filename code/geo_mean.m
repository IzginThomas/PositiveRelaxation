function f = geo_mean(u,phys_flux)
uplus = circshift(u,-1);
f = sqrt(u.*uplus);
end