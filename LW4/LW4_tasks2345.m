% 
% Лабораторная работа 4
% Численные методы (продолжение). Интеграция с языками С/С++.
% Вариант 2
%
%% Задание 2
% Реализовать mex-фунцию [A, B, C, D] = createspline_c(x, f), 
% рассчитывающую коэффициенты кубического сплайна по вектору значений 
% функции f, заданных на узлах сетки x. Реализовать аналогичную функцию 
% [A, B, C, D] = createspline_m(x, f) простейшими средствами Matlab 
% (циклы; оператором двоеточия пользоваться нельзя).

f = @(x) sin(x);
x_left = -pi;
x_right = pi;
x = linspace(x_left, x_right, 20)';

[A1, B1, C1, D1] = createspline_cpp(x, f(x));
% A1
% B1
% C1
% D1
[A2, B2, C2, D2] = createspline_m(x, f(x));
% A2
% B2
% C2
% D2
spl = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;

figure('Position', [100 100 1000 1000]);
hold on;
plot(x, f(x));
for i = 2 : size(x, 1)
    grid = linspace(x(i - 1), x(i), 10);
    plot(grid, spl(A1(i - 1), B1(i - 1), C1(i - 1), D1(i - 1), grid, x(i)), 'r');
    plot(grid, spl(A2(i - 1), B2(i - 1), C2(i - 1), D2(i - 1), grid, x(i)), 'g');
    % (A[i] + B[i] * (x[i] - x[i + 1]) + 0.5 * C[i] * pow(x[i] - x[i + 1], 2) + (1/6) * D[i] * pow(x[i] - x[i + 1], 3)) - f[i]
    % (A(i - 1) + B(i - 1) * (x(i - 1) - x(i)) + 0.5 * C(i - 1) * (x(i - 1) - x(i))^2 + (1/6) * D(i - 1) * (x(i - 1) - x(i))^3) - y(x(i - 1))
end
legend({'sin', 'spline\_cpp', 'spline\_m'});
xlabel('x');
ylabel('y(x)');
hold off;


%% Задание 3
% Сравнить точность функций interp1 (с ключом spline), spline 
% (стандартные матлабовские функции), createspline_cpp, createspline_m 
% для сеток различной длины, построив соответствующие графики.

f = @(x) cos(x);

x_left = -pi;
x_right = pi;
start_point_cnt = 10;
x = linspace(x_left, x_right, start_point_cnt);

grid_dim = start_point_cnt * [2, 4, 6, 8, 10, 12, 14, 16, 18, 20];
abs_error = zeros(4, size(grid_dim, 2));

spl = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;

for i = 1 : size(grid_dim, 2)
    grid = linspace(x_left, x_right, grid_dim(i));
    
    f_interp1 = interp1(x, f(x), grid, 'spline');
    abs_error(1, i) = max(abs(f(grid) - f_interp1));
    
    f_spline = spline(x, f(x), grid);
    abs_error(2, i) = max(abs(f(grid) - f_spline));
    
    [A1, B1, C1, D1] = createspline_cpp(x', (f(x))');
    [A2, B2, C2, D2] = createspline_m(x', (f(x))');
    
    for j = 2 : size(x, 2)
        current_grid = grid(grid > x(j - 1) & grid < x(j));
        % createspline_cpp
        current_abs_error_1 = max(abs(f(current_grid) - ...
            spl(A1(j - 1), B1(j - 1), C1(j - 1), D1(j - 1), current_grid, x(j))));
        % createspline_m
        current_abs_error_2 = max(abs(f(current_grid) - ...
            spl(A2(j - 1), B2(j - 1), C2(j - 1), D2(j - 1), current_grid, x(j))));
        
        if (j > 2) 
            if current_abs_error_1 > abs_error(3, i)
                abs_error(3, i) = current_abs_error_1;
            end
            if current_abs_error_2 > abs_error(4, i)
                abs_error(4, i) = current_abs_error_2;
            end
        else
            abs_error(3, i) = current_abs_error_1; % createspline_cpp
            abs_error(4, i) = current_abs_error_2; % createspline_m
        end
    end
end

figure('Position', [100 100 1000 1000]);
hold on;
plot(grid_dim, abs_error(1, :), 'o-', 'Color', 'r'); % interp1
plot(grid_dim, abs_error(2, :), 'o-', 'Color', 'g'); % spline
plot(grid_dim, abs_error(3, :), 'o-', 'Color','b'); % createspline_cpp
plot(grid_dim, abs_error(4, :), 'o-', 'Color', 'k'); % createspline_m
legend({'interp1', 'spline', 'spline\_cpp', 'spline\_m'});
xlabel('dots count');
ylabel('abs\_error');
hold off;

%% Задание 4
% Сравнить быстродействие функций interp1, spline, createspline_cpp, 
% createspline_m для сеток различной размерности, построив соответствующие 
% графики.

f = @(x) sin(x) + cos(x);

x_left = -pi;
x_right = pi;
start_point_cnt = 10;
x = linspace(x_left, x_right, start_point_cnt);
grid_dim = start_point_cnt * [2, 4, 6, 8, 10, 12, 14, 16, 18, 20];

exec_time = zeros(4, size(grid_dim, 2));

calc_time_repeat_cnt = 10;

spl = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;

for i = 1 : size(grid_dim, 2)
    grid = linspace(x_left, x_right, grid_dim(i));
    
    tic
    for k = 1 : calc_time_repeat_cnt
        f_interp1 = interp1(x, f(x), grid, 'spline');
    end
    exec_time(1, i) = toc;
    exec_time(1, i) = exec_time(1, i) / calc_time_repeat_cnt;
    
    tic
    for k = 1 : calc_time_repeat_cnt 
        f_spline = spline(x, f(x), grid);
    end
    exec_time(2, i) = toc;
    exec_time(2, i) = exec_time(2, i) / calc_time_repeat_cnt;
    
    tic
    for k = 1 : calc_time_repeat_cnt
        [A1, B1, C1, D1] = createspline_cpp(x', (f(x))');
    end
    exec_time(3, i) = toc;
    exec_time(3, i) = exec_time(3, i) / calc_time_repeat_cnt;
    
    tic
    for k = 1 : calc_time_repeat_cnt
        [A2, B2, C2, D2] = createspline_m(x', (f(x))');
    end
    exec_time(4, i) = toc;
    exec_time(4, i) = exec_time(4, i) / calc_time_repeat_cnt;
end

figure('Position', [100 100 1000 1000]);
hold on;
plot(grid_dim, exec_time(1, :), 'o-', 'Color', 'r'); % interp1
plot(grid_dim, exec_time(2, :), 'o-', 'Color', 'g'); % spline
plot(grid_dim, exec_time(3, :), 'o-', 'Color','b'); % createspline_cpp
plot(grid_dim, exec_time(4, :), 'o-', 'Color', 'k'); % createspline_m
legend({'interp1', 'spline', 'spline\_cpp', 'spline\_m'});
xlabel('dots cnt');
ylabel('execution time');
hold off;

%% Задание 5
% Обозначим T_s(n) время работы методов из предыдущего пункта на матрицах 
% порядка n (s = spline, createspline_c, ...). Написать функцию, которая, 
% используя линейную регрессию, аппроксимирует эти функции с помощью 
% многочленов степени не выше заданной.
pol_deg = 4;
pol_coef = zeros(4, pol_deg + 1);
pol_coef(1, :) = polyfit(grid_dim, exec_time(1, :), pol_deg); % interp1
pol_coef(2, :) = polyfit(grid_dim, exec_time(2, :), pol_deg); % spline
pol_coef(3, :) = polyfit(grid_dim, exec_time(3, :), pol_deg); % createspline_cpp
pol_coef(4, :) = polyfit(grid_dim, exec_time(4, :), pol_deg); % createspline_m

methods = {'interp1', 'spline', 'createspline\_cpp', 'createspline\_m'};

for i = 1 : 4
    figure('Position', [100 100 1000 1000]);
    hold on;
    plot(grid_dim, exec_time(i, :), 'o-', 'Color', 'b');
    plot(grid_dim, polyval(pol_coef(i, :), grid_dim), 'o-', 'Color', 'r');
    title([methods(i), 'approximation by a polynomial of the', pol_deg, 'degree']);
    legend({'tic-toc', 'approx'});
    xlabel('dots cnt');
    ylabel('execution time');
    hold off
end

%% createspline_m
function [A, B, C, D] = createspline_m(x, f)
    h = circshift(x, -1) - x;
    h(size(x, 1)) = [];
    
    
    A1 = h;
    A1(size(h, 1)) = [];
    
    B1 = h;
    B1(1) = [];
   
    C1 = 2.0 * (A1 + B1);
    
    
    f1 = f;
    f1(1) = [];
    f1(1) = [];
    
    f2 = f;
    f2(size(f2, 1)) = [];
    f2(1) = [];
    
    f3 = f;
    f3(size(f3, 1)) = [];
    f3(size(f3, 1)) = [];
    
    h12 = h;
    h12(1) = [];
    
    h23 = h;
    h23(size(h23, 1)) = [];
    
    F = 6.0 * ((f1 - f2) ./ h12 - (f2 - f3) ./ h23);
    
    
    tmpT = [circshift(A1, -1), C1, circshift(B1, 1)];
    T = spdiags(tmpT, [-1, 0, 1], size(A1, 1), size(A1, 1));
    
    
    % solve T * y = F;
    y = T \ F;
    
    
    tmpA = circshift(f, -1);
    tmpA(size(tmpA, 1)) = [];
    A = tmpA;
    
    C = y;
    C(size(C, 1) + 1) = 0.0;
    
    
    D = (C - circshift(C, 1)) ./ h;
    
    tmp1B = f;
    tmp1B(1) = [];
    tmp2B = f;
    tmp2B(size(tmp2B, 1)) = [];
    
    B = (0.5 * h .* C) - ((1/6) * (h.^2) .* D) + (tmp1B - tmp2B) ./ h;
end