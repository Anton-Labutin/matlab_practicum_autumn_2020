% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 6
% 
% Для линейной системы дифференциальных уравнений второго порядка 
% построить фазовый портрет. Подобрать системы таким образом, чтобы 
% проиллюстрировать различные виды особых точек (узел, дикритический 
% узел, седло, фокус, центр). На траекториях должно быть видно 
% направление движения (стрелки).


angle = linspace(0, 2 * pi, 100);
start = [cos(angle); sin(angle)];

time = linspace(-0.1, 0.1, 10);


% неустойчивый узел
% dx/dt = 2x + 3y,
% dy/dt = x + 4y
% lyambda_1 = 1, lyambda2 = 5 - собственные значения
% h1 = [-3; 1], h2 = [1; 1] - собственные векторы
% решение:
% x(t) = -3 * c1 * exp(t) + c2 * exp(5t)
% y(t) = c1 * exp(t) + c2 * exp(5t)

dxdt = @(x, y) 2 * x + 3 * y;
dydt = @(x, y) x + 4 * y;

phase_portrait(dxdt, dydt, @odefun1, time, start, 'Неустойчивый узел');


% устойчивый узел
% dx/dt = -6x + 3y,
% dy/dt = -2x - y
% lyambda_1 = -3, lyambda2 = -4 - собственные значения
% h1 = [1; 1], h2 = [3; 2] - собственные векторы
% решение:
% x(t) = c1 * exp(-3t) + 3 * c2 * exp(-4t)
% y(t) = c1 * exp(-3t) + 2 * c2 * exp(-4t)

dxdt = @(x, y) -6 * x + 3 * y;
dydt = @(x, y) -2 * x - y;

phase_portrait(dxdt, dydt, @odefun2, time, start, 'Устойчивый узел');


% неустойчивый дикритический узел
% dx/dt = 3x,
% dy/dt = 3y
% lyambda_1 = 3 - собственное значение кратности 2
% h1 = [1; 0], h2 = [0; 1] - собственные векторы
% решение:
% x(t) = c1 * exp(3t)
% y(t) = c2 * exp(3t)

dxdt = @(x, y) 3 * x;
dydt = @(x, y) 3 * y;

phase_portrait(dxdt, dydt, @odefun3, time, start, 'Неусточивый дикритический узел');


% устойчивый дикритический узел
% dx/dt = -2x,
% dy/dt = -2y
% lyambda_1 = -2 - собственное значение кратности 2
% h1 = [1; 0], h2 = [0; 1] - собственные векторы
% решение:
% x(t) = c1 * exp(-2t)
% y(t) = c2 * exp(-2t)

dxdt = @(x, y) -2 * x;
dydt = @(x, y) -2 * y;

phase_portrait(dxdt, dydt, @odefun4, time, start, 'Усточивый дикритический узел');


% седло
% dx/dt = 3x + 4y,
% dy/dt = 2x + y;
% lyambda_1 = -1, lyambda = 5 - собственные значения
% h1 = [-1; 1], h2 = [2; 1] - собственные векторы
% решение:
% x(t) = -c1 * exp(-t) + 2 * c2 * exp(5t)
% y(t) = c1 * exp(-t) + c2 * exp(5t)

dxdt = @(x, y) 3 * x + 4 * y;
dydt = @(x, y) 2 * x + y;

phase_portrait(dxdt, dydt, @odefun5, time, start, 'Седло');


% Неустойчивый фокус
% dx/dt = y,
% dy/dt = -2x + y;
% lyambda_1 = 0.5 + sqrt(7)/2 * i, lyambda = 0.5 - sqrt(7)/2 * i - собственные значения
% собственные векторы:
% h1 = sqrt(e) * [cos(sqrt(7)/2 * t) + sqrt(7) * sin(sqrt(7)/2 * t); 
%                   4 * cos(sqrt(7)/2 * t)], 
% h2 = sqrt(e) * [sqrt(7) * cos(sqrt(7)/2 * t) - sin(sqrt(7)/2 * t); 
%                   -4 * sin(sqrt(7)/2 * t)] 
% решение:
% x(t) = c1 * h11 + c2 * h12,
% y(t) = c1 * h21 + c2 * h22

dxdt = @(x, y) y;
dydt = @(x, y) -2 * x + y;

phase_portrait(dxdt, dydt, @odefun6, time, start, 'Неустойчивый фокус');


% Устойчивый фокус
% dx/dt = -y,
% dy/dt = -2x + y;
% lyambda_1 = -0.5 + sqrt(7)/2 * i, lyambda = -0.5 - sqrt(7)/2 * i - собственные значения
% собственные векторы:
% h1 = e(-1/2) * [cos(sqrt(7)/2 * t) - sqrt(7) * sin(sqrt(7)/2 * t); 
%                   -4 * cos(sqrt(7)/2 * t)], 
% h2 = e^(-1/2) * [sqrt(7) * cos(sqrt(7)/2 * t) + sin(sqrt(7)/2 * t); 
%                   -4 * sin(sqrt(7)/2 * t)] 
% решение:
% x(t) = c1 * h11 + c2 * h12,
% y(t) = c1 * h21 + c2 * h22

dxdt = @(x, y) y;
dydt = @(x, y) -2 * x - y;

phase_portrait(dxdt, dydt, @odefun7, time, start, 'Устойчивый фокус');


% Центр
% dx/dt = -2x - 5y,
% dy/dt = 2x + 2y;
% lyambda_1 = sqrt(6) * i, lyambda = -sqrt(6) * i - собственные значения
% собственные векторы:
% h1 = [2 * cos(sqrt(6) * t) + sqrt(6) * sin(sqrt(6) * t); 
%                   -2 * cos(sqrt(6) * t)], 
% h2 = [sqrt(6) * cos(sqrt(6) * t) - 2 * sin(sqrt(6) * t); 
%                   4 * sin(sqrt(6) * t)] 
% решение:
% x(t) = c1 * h11 + c2 * h12,
% y(t) = c1 * h21 + c2 * h22

dxdt = @(x, y) -2 * x - 5 * y;
dydt = @(x, y) 2 * x + 2 * y;

phase_portrait(dxdt, dydt, @odefun8, time, start, 'Центр');


function dydt = odefun1(t, y)
    % dx/dt = 2x + 3y,
    % dy/dt = x + 4y

    dydt = zeros(2, 1);
    dydt(1) = 2 * y(1) + 3 * y(2);
    dydt(2) = y(1) + 4 * y(2);
end


function dydt = odefun2(t, y)
    % dx/dt = -6x + 3y,
    % dy/dt = -2x - y

    dydt = zeros(2, 1);
    dydt(1) = -6 * y(1) + 3 * y(2);
    dydt(2) = -2 * y(1) - y(2);
end


function dydt = odefun3(t, y)
    % dx/dt = 3x,
    % dy/dt = 3y

    dydt = zeros(2, 1);
    dydt(1) = 3 * y(1);
    dydt(2) = 3 * y(2);
end


function dydt = odefun4(t, y)
    % dx/dt = -2x,
    % dy/dt = -2y

    dydt = zeros(2, 1);
    dydt(1) = -2 * y(1);
    dydt(2) = -2 * y(2);
end


function dydt = odefun5(t, y)
    % dx/dt = 3x + 4y,
    % dy/dt = 2x + y;

    dydt = zeros(2, 1);
    dydt(1) = 3 * y(1) + 4 * y(2);
    dydt(2) = 2 * y(1) + y(2);
end


function dydt = odefun6(t, y)
    % dx/dt = y,
    % dy/dt = -2x + y;

    dydt = zeros(2, 1);
    dydt(1) = y(2);
    dydt(2) = -2 * y(1) + y(2);
end


function dydt = odefun7(t, y)
    % dx/dt = -y,
    % dy/dt = -2x + y;

    dydt = zeros(2, 1);
    dydt(1) = y(2);
    dydt(2) = -2 * y(1) - y(2);
end


function dydt = odefun8(t, y)
    % dx/dt = -2x - 5y,
    % dy/dt = 2x + 2y;

    dydt = zeros(2, 1);
    dydt(1) = -2 * y(1) - 5 * y(2);
    dydt(2) = 2 * y(1) + 2 * y(2);
end


function phase_portrait(dxdt, dydt, odefun, time, start, title_name)
    figure('Position', [100 100 1000 1000]);
    hold on;
    for start_idx = 1 : size(start, 2)
        x0 = start(1, start_idx);
        y0 = start(2, start_idx);
    
        [t, y] = ode45(odefun, time, [x0; y0]);
        plot(y(:, 1), y(:, 2));
        quiver(y(:, 1), y(:, 2), dxdt(y(:, 1), y(:, 2)), dydt(y(:, 1), y(:, 2)));
    end
    title(title_name);
    hold off;
end