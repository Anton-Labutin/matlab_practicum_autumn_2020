% 
% Лабораторная работа 2
% Работа с графикой
% Вариант 2
%
%%
%
% Двумерная графика, часть 2
%
%%
% Task 9
%
% Используя функцию fmincon, написать функцию supportLebesgue(f, opts), 
% которая выдает (приближенную) опорную функцию множества X, которую 
% можно использовать в предыдущем задании. Функция f предполагается 
% выпуклой. Параметры opts являются параметрами функции fmincon.

f = @(x) (x(1)^2 + x(2)^2 - 9);

opts = struct();
opts.x0 = [0, -0.5];
opts.A = [[0, 1]; [0, 0]];
opts.b = zeros(2, 1);
opts.Aeq = zeros(2);
opts.beq = zeros(2, 1);
opts.lb = [];
opts.ub = [];

[support_func, support_x] = supportLebesgue(f, opts);
l = [0; -2];
support_x(l)
support_func(l)


function [support_func, support_x] = supportLebesgue(f, opts)
    scalar_product = @(x, l) (x(1) * l(1) + x(2) * l(2));
    
    function [c, ceq] = mycon(x)
        c = f(x);
        ceq = [];
    end

    sup_x = @(l) fmincon(@(x) -scalar_product(x, l), ...
        opts.x0, opts.A, opts.b, opts.Aeq, opts.beq, opts.lb, opts.ub, ...
        @mycon);
        
    support_x = @(l) sup_x(l);
    support_func = @(l) scalar_product(support_x(l), l);
end