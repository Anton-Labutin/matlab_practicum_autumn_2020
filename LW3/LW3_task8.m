% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 8
% 
% С помощью функции bvp4c решить численно краевую задачу. Сравнить решение 
% с аналитическим в L2-норме и C-норме.

y = @ (x, c) pi * cos(x) + c * sin(x) + 2 * x - pi;

x_left = 0;
x_right = pi;

xmesh = linspace(x_left, x_right, 5);
solinit = bvpinit(xmesh, @guess);

sol = bvp4c(@bvpfcn, @bcfcn, solinit);
% sol.x
% sol.y

% let y(pi/2) = 1  =>  c = 1
figure('Position', [100 100 1000 1000]);
plot(sol.x, sol.y(1, :), '-o', sol.x, y(sol.x, 1), '-*');
legend('numerical solution', 'exact solution: pi * cos(x) + c * sin(x) + 2 * x - pi');

L2_norm = sqrt(trapz(sol.x, (y(sol.x, 1) - sol.y(1, :)).^2))
c_norm = max(abs(y(sol.x, 1) - sol.y(1, :)))


function dydx = bvpfcn(x, y)
    dydx = zeros(2, 1);
    dydx = [y(2); -y(1) + 2 * x - pi];
end


function res = bcfcn(ya, yb)
    res = [ya(1); yb(1)];
end


function g = guess(x)
    g = [sin(x); cos(x)];
end