function [c,flag_relax] = bisectionMethod(f,a,b,fa,error,max_it) %(ensure change of sign between a and b)
c=(a+b)/2;
i = 0;
fc = f(c);
while abs(fc) > error  && i <= max_it && abs(b-a) > error
    if fc<0 && fa<0
        a=c;
        fa = fc;
    else
        b=c;
    end
    c=(a+b)/2;
    fc = f(c);
    i = i + 1;
end
if abs(fc) <= error && c > 0
    flag_relax = true;
else
    flag_relax = false;
end