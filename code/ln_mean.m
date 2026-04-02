function out = ln_mean(x,y)
epsilon = 1e-4;
f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y + realmin);
if f2 < epsilon
    out = (x + y) / (2 + 2/3 * f2 + 2/5 * f2^2 + 2/7 * f2^3 + realmin);
else
    out = (y - x) / (log(y / x) + realmin);
end