% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 9
% 
% Реализовать функцию, ищущую минимум функции многих переменных методом 
% покоординатного спуска. (Функцию, её частные производные и начальное 
% приближение задаёт пользователь.) Для функции двух переменных построить 
% набор линий уровня, на которых отметить шаги алгоритма. Сравнить 
% результат работы с функцией fminbnd.

eps = 10^(-5);

f = @(x) 4 * (x(1) - 5).^2 + (x(2) - 6).^2;
init_approx = [15; 15];
borders = [-20 20];

[X, Y] = meshgrid(linspace(borders(1), borders(2), 100), ...
     linspace(borders(1), borders(2), 100));
Z = 4 * (X - 5).^2 + (Y - 6).^2;

x_min = zeros(size(init_approx, 1), 1);
start_flag = 1;
temp = x_min;

figure('Position', [100 100 1000 1000]);
hold on;
iter_idx = 1;
while norm(f(x_min) - f(temp)) >= eps || start_flag == 1
    for i = 1 : size(init_approx, 1)
        if norm(f(x_min) - f(temp)) < eps && start_flag == 0
            break
        end
        disp(x_min);
        contour(X, Y, Z, [f(x_min) f(x_min)], 'ShowText', 'on');
        text(x_min(1), x_min(2), num2str(iter_idx), 'FontSize', 10, 'Color', 'red');
        iter_idx = iter_idx + 1;
        temp = x_min;
        x_min(i) = golden_section(f, i, x_min, borders, eps);
    end
    start_flag = 0;
end
hold off;

x_min = x_min
f_x_min = f(x_min)



f = @(x) sin(x);
init_approx = [pi/2];
borders = [-pi, pi];
% borders = [pi/4, 5*pi/4]

x_min = coordinate_descent(f, init_approx, borders, eps)
f_x_min = f(x_min)

fminbnd_x_min = fminbnd(f, borders(1), borders(2))
f_fminbnd_x_min = f(fminbnd_x_min)




function x_min = coordinate_descent(f, init_approx, borders, eps)
    x_min = zeros(size(init_approx, 1), 1);
    start_flag = 1;
    temp = x_min;
    
    while norm(f(x_min) - f(temp)) >= eps || start_flag == 1
        for i = 1 : size(init_approx, 1)
            temp = x_min;
            x_min(i) = golden_section(f, i, x_min, borders, eps);
        end
        start_flag = 0;
    end
end


function x_min_i = golden_section(f, f_arg_idx, x_min, segment, eps)
    phi = (1 + sqrt(5)) / 2;
    a = segment(1);
    b = segment(2);
    x1 = b - (b - a) / phi;
    x2 = a + (b - a) / phi;
    temp1 = x_min;
    temp2 = x_min;
    
    while (b - a) / 2 >= eps
        temp1(f_arg_idx) = x1;
        temp2(f_arg_idx) = x2;
        
        if f(temp1) > f(temp2)
            a = x1;
            x1 = x2;
            x2 = b - (x1 - a);
        elseif f(temp1) < f(temp2)
            b = x2;
            x2 = x1;
            x1 = a + (b - x2);
        end
    end
    
    x_min_i = (a + b) / 2;
end