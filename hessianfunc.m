function yprpr = hessianfunc(x)
    syms x1 x2
    F= 1 + sin(x1)^2 + sin(x2)^2 - 0.1 * exp(-x1^2 - x2^2);
    H=hessian(F);
    yprpr=double(subs(H,[x1 x2], [x(2) x(2)]));
end