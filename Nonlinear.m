clear all
clc

X1=-10:0.01:10;
X2=-10:0.01:10;
[x1,x2]=meshgrid(X1,X2);
F = 1 + sin(x1).^2 + sin(x2).^2 - 0.1 * exp(-x1.^2 - x2.^2);
realFMin=min(min(F))
mesh(x1,x2,F)

figure
contourf(x1,x2,F)
hold on

algorithms = {'Newton-Raphson', 'Hestenes-Stiefel', 'Polak-Ribiere', 'Fletcher-Reeves'};
executionTimes = zeros(1, 4);

%% Newton-Raphson
fprintf('Newton-Raphson Algorithm\n');
x0 = -10 + 20 * rand(2, 1); % -> [-10, 10]
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
plot(x(1),x(2),'r.')

x_next= x - inv(hessianfunc(x))*gradfunc(x);

fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')
k=3;
%while(norm(func(x_next)-func(x))>epsilon)
while(abs(gradfunc(x_next))>epsilon)
    %-
    x=x_next;
    x_next= x - inv(hessianfunc(x))*gradfunc(x);
    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k=k+1;
    %-
end
elapsed1=toc
executionTimes(1) = elapsed1;
title('Newton-Raphson Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Hestenes-Stiefel Algorithm
figure
contourf(x1,x2,F)
hold on

fprintf('Hestenes-Stiefel Algorithm\n');
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))

plot(x(1),x(2),'r.')
g=gradfunc(x);
d=-g;

%alpha argmin procedure
alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);
%end of alpha procedure

x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta=(g_next'*(g_next-g))/(d'*(g_next-g));
d_next=-g_next+beta*d;

fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')
k=3;
while(norm(gradfunc(x_next))>epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    %alpha argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    %end of alpha procedure

    x_next=x+alpha*d;
    g_next=gradfunc(x_next);
    beta=(g_next'*(g_next-g))/(d'*(g_next-g));
    d_next=-g_next+beta*d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k=k+1;

end
elapsed2=toc
executionTimes(2) = elapsed2;
title('Hestenes-Stiefel Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);


%% Polak-Ribi`ere Algorithm
figure
contourf(x1,x2,F)
hold on
fprintf('Polak-Ribi Algorithm\n');

x = x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))

plot(x(1),x(2),'r.')
g = gradfunc(x);
d = -g;

%alpha argmin procedure
alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);
%end of alpha procedure

x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta= (g_next'*(g_next-g))/(g'*g);
d_next=-g_next+beta*d;

fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')

k = 3;
while (norm(gradfunc(x_next)) > epsilon) && (abs(func(x_next) - func(x)) > epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    %alpha argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    %end of alpha procedure

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);
    beta = (g_next'*(g_next-g))/(g'*g);
    d_next = -g_next + beta * d;
    
    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n', k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)))
    plot(x_next(1), x_next(2), 'r*')
    k = k + 1;
end

elapsed3=toc
executionTimes(3) = elapsed3;
title('Polak-Ribi Algorithm')
set(gca, 'fontsize', 35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);


%% Fletcher-Reeves Algorithm
figure
contourf(x1,x2,F)
hold on
fprintf('Fletcher-Reeves Algorithm\n');

x = x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))

plot(x(1),x(2),'r.')
g = gradfunc(x);
d = -g;

%alpha argmin procedure
alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);
%end of alpha procedure

x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta= ((g_next)'*(g_next))/(g'*g);
d_next=-g_next+beta*d;

fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')

k = 3;

while (norm(gradfunc(x_next)) > epsilon) && (abs(func(x_next) - func(x)) > epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    %alpha argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    %end of alpha procedure

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);
    beta= ((g_next)'*(g_next))/(g'*g);
    d_next = -g_next + beta * d;
    
    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n', k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)))
    plot(x_next(1), x_next(2), 'r*')
    k = k + 1;
end

elapsed4=toc
executionTimes(4) = elapsed4;
title('Fletcher-Reeves Algorithm')
set(gca, 'fontsize', 35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

figure;
bar(executionTimes);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('Execution Time (seconds)');
title('Execution Time of Different Optimization Algorithms');
set(gca, 'fontsize', 20);