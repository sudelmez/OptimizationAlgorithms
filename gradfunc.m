function ypr = gradfunc(x)
    syms x1 x2
    F=1 + sin(x1)^2 + sin(x2)^2 - 0.1 * exp(-x1^2 - x2^2);
    g=gradient(F);
    ypr=double(subs(g,[x1 x2], [x(2) x(2)]));
end