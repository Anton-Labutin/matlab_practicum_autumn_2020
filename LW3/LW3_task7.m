% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 7
% 
% Для систем
%   x' = x^3 - y, x из R,
%   y' = x + y^3, y из R,
% и
%   x' = 2*y^3 - x^5, x из R,
%   y' = -x - y^3 + y^5, y из R,
% исследовать на устойчивость нулевое положение равновесия, построив 
% функцию Ляпунова и применив теоремы Ляпунова или Четаева. Нарисовать 
% фазовый портрет системы и линии уровня функции Ляпунова. Траектории 
% нарисовать меняющими цвет в соответствии со значением функции Ляпунова 
% вдоль траектории (например, чем больше — тем краснее).

clc
clear

%   x' = x^3 - y, x из R,
%   y' = x + y^3, y из R,

% функция Ляпунова
V = @(x, y) x.^2 + y.^2;

% V - непрерывно-дифференцируемая
% V > 0 в области U = { x из R^2: 0 < |x| < eps }
% V = 0 на границе области (точка (0,0))
% dV/dt = 2 * x * (x^3 - y) + 2 * y * (x + y^3) = 2 * (x^4 + y^4)
% dV/dt > 0 в области U
% в области  U(a) = {x из R^2: sqrt(a) < |x| < eps} dV/dt >= (c * a^2) > 0, 
% где с = const
% => по теореме Четаева нулевое решение неустойчиво
    
dxdt = @(x, y) x.^3 - y;
dydt = @(x, y) x + y.^3;

% начальные данные
t0 = 0;
t1 = 1;
time = [t0, t1];

n = 30;
angle = linspace(0, 2 * pi, n);
start = [cos(angle); sin(angle)];

border = 8;
grid_step = 0.02;
level_line_cnt = 20;

save lw3_task7.mat border;


% линии уровня функции Ляпунова
figure('Position', [100 100 1000 1000]);
axis([-border border -border border]);

hold on;
[X, Y] = meshgrid(-border : grid_step : border, -border : grid_step : border);
Z = V(X, Y);
contour(X, Y, Z, level_line_cnt);


% фазовый портрет
phase_portrait(dxdt, dydt, @ode_ls_1, time, start, -border : 15 * grid_step : border, 'System 1');
hold off
    
 
% x' = 2*y^3 - x^5, x из R,
% y' = -x - y^3 + y^5, y из R,

% функция Ляпунова
V = @(x, y) x.^2 + y.^4;
% V > 0 при х != 0 и V = 0 на границе области (точка (0,0)) 
% dV/dt = -2 * x^6 - 4 * y^6 * (1 - y^2) <= 2 * x^6 - 4 * y^6 < 0
% dV/dt < 0 в области U = {x из R2, t из R: |x| < eps, t > t0} 
% => по теореме Ляпунова нулевое решение устойчиво
   
dxdt = @(x, y) 2 * y.^3 - x.^5;
dydt = @(x, y) -x - y.^3 + y.^5;

% начальные данные
t0 = -10;
t1 = 10;
time = [t0, t1];

n = 50;
angle = linspace(0, 2 * pi, n);
r = 6;
start = r * [cos(angle); sin(angle)];

border = 20;
grid_step = 0.02;
level_line_cnt = 20;

save lw3_task7.mat border;

% линии уровня
figure('Position', [100 100 1000 1000]);
axis([-border border -border border]);

hold on;
[X, Y] = meshgrid(-border : grid_step : border, -border : grid_step : border);
Z = V(X, Y);
contour(X, Y, Z, level_line_cnt);
    
% фазовый портрет
phase_portrait(dxdt, dydt, @ode_ls_2, time, start, -border : 25 * grid_step : border, 'System 2');
hold off
 

function dydt = ode_ls_1(t, y)
    % x' = x^3 - y, x из R,
    % y' = x + y^3, y из R,
    
    dydt = zeros(2, 1);
    dydt(1) = y(1)^3 - y(2); ...
    dydt(2) = y(1) + y(2)^3;
end


function dydt = ode_ls_2(t, y)
    dydt = zeros(2, 1);
    dydt(1) = 2 * y(2)^3 - y(1)^5;
    dydt(2) = -y(1) - y(2)^3 + y(2)^5;
end


function phase_portrait(dxdt, dydt, odefun, time, start, grid, title_name)
    options = odeset('Events', @events);
    for start_idx = 1 : size(start, 2)
        x0 = start(1, start_idx);
        y0 = start(2, start_idx);
        
        [t, y] = ode45(odefun, time, [x0; y0], options);
        
        for i = 1 : size(y, 1) - 1
            line = plot(y(i : i+1, 1), y(i : i+1, 2), 'LineWidth', 1); 
            line.Color = [1, 1 - i / (size(y, 1) - 1), 0];
        end 
    end
    
    [X, Y] = meshgrid(grid, grid);
    quiver(X, Y, dxdt(X, Y), dydt(X, Y));
    daspect([1 1 1]);
    
    title(title_name);
    legend('линии уровня', 'траектория системы');
    xlabel('x');
    ylabel('y');
end


function [value, isterminal, direction] = events(t, x)
    % value(i), isterminal(i), direction(i) corresponds to i-th event
    % value(i) - the value of the ith event function
    % isterminal(i) = 1 if the integration is to terminate at a zero of 
    %   this event function. Otherwise, it is 0.
    % direction(i) = 0 if all zeros are to be located (the default). 
    %   A value of +1 locates only zeros where the event function is 
    %   increasing, and -1 locates only zeros where the event function 
    %   is decreasing.
    
    load lw3_task7.mat border;
    
    value = [ (border - x(1)) > 1, (x(1) + border) > 1, ... 
        (border - x(2)) > 1, (x(2) + border) > 1 ];
    isterminal = [1, 1, 1, 1];
    direction = [0, 0, 0, 0];
end