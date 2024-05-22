function y = func(x)
  x1=x(1);
  x2=x(2);
  y= 1 + sin(x1)^2 + sin(x2)^2 - 0.1 * exp(-x1^2 - x2^2);
end