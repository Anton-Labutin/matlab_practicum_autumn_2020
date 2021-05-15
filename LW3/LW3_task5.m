% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 5
% Рассмотреть систему двух тел на плоскости. Решить систему численно. 
% Нарисовать на плоскости анимацию движения траекторий x1(t), x2(t). 
% Подобрать параметры так, что бы продемонстрировать движение двух типов: 
% по пересекающимся орбитам («восьмёрка») и вокруг общего центра. 

clear
clc

% движение по пересекающимся орбитам
% x10 = [-1; 0];
% x20 = [1; 0];
% dx1_dt0 = [0; 2];
% dx2_dt0 = [0; -2];
% 
% G = 100;
% m1 = 1;
% m2 = 1;


% движение по окружности
x10 = [-2; 0];
x20 = [2; 0];
dx1_dt0 = [0; 5];
dx2_dt0 = [0; -5];

G = 200;
m1 = 1;
m2 = 1;

% движение по винтовой линии
% x10 = [0; 0];
% x20 = [0; 0.4];
% dx1_dt0 = [5; 5];
% dx2_dt0 = [-5; -5];
% 
% G = 0.367;
% m1 = 120;
% m2 = 100;

save lw3_task5.mat G m1 m2;

y0 = [x10; x20; dx1_dt0; dx2_dt0];

t0 = 0;
t1 = 10;

[t, y] = ode45(@f, [t0, t1], y0);

figure('Position', [100 100 1000 1000]);
axis_x_min = min(min(y(:, 1)), min(y(:, 3)));
axis_x_max = max(max(y(:, 1)), max(y(:, 3)));
axis_y_min = min(min(y(:, 2)), min(y(:, 4)));
axis_y_max = max(max(y(:, 2)), max(y(:, 4)));
axis([axis_x_min axis_x_max axis_y_min axis_y_max]);

hold on
% max_time_point_cnt = 1000;
% n = min(size(t, 1), max_time_point_cnt);
n = size(t, 1);
for i = 1 : n
    plot(y(1 : i, 1), y(1 : i, 2), 'r', y(1 : i, 3), y(1 : i, 4), 'b');
    pause(0.01);
end
legend({'x1(t)', 'x2(t)'});
xlabel('x\_1');
ylabel('x\_2');
hold off


function [dxdt] = f(t, x)
    load lw3_task5.mat G m1 m2;
    
    % x1 = [x11, x12]
    x11 = x(1);
    x12 = x(2);
    % x2 = [x21, x22]
    x21 = x(3);
    x22 = x(4);
    % dx1/dt = [dx11/dt, dx12/dt]
    dx11_dt = x(5);
    dx12_dt = x(6);
    % dx2/dt = [dx21/dt, dx22/dt]
    dx21_dt = x(7);
    dx22_dt = x(8);
    
    % norm = @(x) sqrt(x(1)^2 + x(2)^2);
    norm12_3 = (norm([x11 - x21; x12 - x22]))^3;
    
    % 1st diff eq
    d2x11_dt2 = G * m2 * (x21 - x11) / norm12_3;
    d2x12_dt2 = G * m2 * (x22 - x12) / norm12_3;
    % 2nd diff eq
    d2x21_dt2 = G * m1 * (x11 - x21) / norm12_3;
    d2x22_dt2 = G * m1 * (x12 - x22) / norm12_3;
    
    dxdt = [ dx11_dt; dx12_dt; ... 
             dx21_dt; dx22_dt; ...
             d2x11_dt2;  d2x12_dt2; ...
             d2x21_dt2;   d2x22_dt2 ];
end