% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 4
% Движение шарика на плоскости описывается уравнением x'' = −αx, x ∈ R2. 
% Реализовать моделирование (см. ode45) движения шарика внутри участка, 
% окруженного четырьмя перегородками, параллельными осям координат. 
% При попадании на перегородку шарик от нее упруго отскакивает 
% (так, при ударе о вертикальную стенку в момент t′ вертикальная 
% компонента скорости меняет знак: x1(t′+0) = x1(t′−0), x2(t′+0) = 
% −x2(t′−0), и так далее). Нарисовать анимацию, изображающую движение 
% шарика с ненулевой начальной скоростью.

clear
clc

% начальные данные
t0 = 0;
tf = 20;

x0 = [5; 6];
dx_dt0 = [2; 3];
y0 = [x0; dx_dt0];

alpha = 0.1;
left_lower_bound = [0, 0];
right_upper_bound = [10, 10];
save('lw3_task4_temp.mat', 'alpha', 'left_lower_bound', 'right_upper_bound');

eps = 1e-6;

options = odeset('Events', @events);

% результаты
tout = t0;
yout = y0';
    
while t0 < tf
    [t, y, te, ye, ie] = ode45(@f, [t0, tf], y0, options);
    % t - time
    % y - solution
    % te - time of the events
    % ye - solution at 'te'
    % ie - the index of the triggered event
    
    tout = [tout; t];
    yout = [yout; y];
    
    if((t(end) >= tf) ... time is out
        || (isempty(ie) == 1)) % stop inside the field
        break;
    end
    t0 = te;
    % может выйти так, что за время вычисления траектории значения
    % несколько раз вышли за ограничения, поэтому делаем проверки
    
    % ye = [x_1, x_2, dx_1/dt, dx_2/dt]
    % если вышли за левый край или дошли до левого края, то
    if (ye(end, 1) - left_lower_bound(1)) <= eps && ye(end, 3) < 0
        y0(3) = -ye(end, 3);
        y0(1) = left_lower_bound(1) + eps;
    else
        % если вышли за правый край или дошли до правого края, то
        if (right_upper_bound(1) - ye(end, 1)) <= eps && ye(end, 3) > 0
            y0(3) = -ye(end, 3);
            y0(1) = right_upper_bound(1) - eps;
        else
            y0(3) = ye(end, 3);
            y0(1) = ye(end, 1);
        end
    end
    % если вышли за нижний край или дошли до нижнего края, то
    if (ye(end, 2) - left_lower_bound(2)) <= eps && ye(end, 4) < 0
        y0(4) = -ye(end, 4);
        y0(2) = left_lower_bound(2) + eps;
    else
        % если вышли за верхний край или дошли до верхнего края, то
        if (right_upper_bound(2) - ye(end, 2)) <= eps && ye(end, 4) > 0
            y0(4) = -ye(end, 4);
            y0(2) = right_upper_bound(2) - eps;
        else
            y0(4) = ye(end, 4);
            y0(2) = ye(end, 2);
        end
    end
end

% траектория мяча в промежутке времени [t0, tf]
figure('Position', [100 100 1000 1000]);
axis([left_lower_bound(1), right_upper_bound(1), ... 
    left_lower_bound(2), right_upper_bound(2)]);

hold on
for i = 1 : size(yout, 1)
    plot(yout(1 : i, 1), yout(1 : i, 2), 'b');
    pause(0.1);
end
xlabel('x');
ylabel('y');
hold off

% график зависимости координат х и у от времени t
% figure('Position', [100 100 1000 1000]);
% plot(tout, yout(:, 1), 'r', tout, yout(:, 2), 'g');
% xlabel('time');
% ylabel('переменные х и у');
% legend('x = x(t)', 'y = y(t)');
% title('Зависимость переменных х и у от времени');


function dxdt = f(t, x)
    % x'' = -alpha * x 
    % <=>
    % x' = y,
    % y' = -alpha * x
    
    load lw3_task4_temp.mat alpha; 
    dxdt = [x(3); x(4); -alpha * x(1); -alpha * x(2)];
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
    
    load lw3_task4_temp.mat left_lower_bound right_upper_bound;
    
    value = [ (x(1) - left_lower_bound(1)) * (right_upper_bound(1) - x(1)), ...
              (x(2) - left_lower_bound(2)) * (right_upper_bound(2) - x(2)) ];
    isterminal = [1, 1];
    direction = [0, 0];
end