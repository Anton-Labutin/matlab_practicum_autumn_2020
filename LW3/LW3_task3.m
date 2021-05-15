% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Для функции нарисовать график на отрезке [−a, a]: по оси абсцисс — 
% начальное приближение, по оси ординат — ближайший к начальному 
% приближению корень функции, найденный с помощью fzero.

a = 1; 
step = 0.001;
x = -a : step : a; 


figure('Position', [100 100 1000 1000]);
plot(x, f(x), 'r');
xlabel('x');
ylabel('y');
title('y = f(x)');
grid on;


figure('Position', [100 100 1000 1000]);
hold on;
for i = 1 : numel(x)
    [x0, y0] = fzero(@f, x(i));
    plot(x(i), x0, 'ro', x(i), y0, 'bo', 'MarkerSize', 3);
end
hold off;
legend("root","approximation value");
xlabel('x');
ylabel('y');
title('Roots and approximation values');
grid on;


function [y] = f(x)
    eps = 0.0001;
    y = x;
    y(abs(x) < eps) = 1;
    y = y .* cos(log(abs(y)));
    y(abs(x) < eps) = 0;
end