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
% Task 7
%
% Написать функцию getEqual(f, g, t0, t1, N), которая принимает на вход 
% две функции, описывающие параметрическую кривую на плоскости, и 
% возвращает N точек таких, что ... Продемонстрировать работу на 
% фигурах Лиссажу.

A = 2;
B = 1;
a = 3;
b = 1;
d = 2;

x = @(t) A * sin(a * t + d);
y = @(t) B * sin(b * t);

getEqual(x, y, 1, 10, 50);
axis equal;

function result = getEqual(f, g, t0, t1, N)
    t_point_cnt = 3000;
    t = linspace(t0, t1, t_point_cnt);
    x = f(t);
    y = g(t);
    
    N_points_grid = 1 : N;
    N_points_grid(N) = t_point_cnt;
    
    flag = 1;
    max_norm_idx = 0;
    while flag
        max_norm = 0;
        for i = 1 : N-1
            dist = norm([x(N_points_grid(i + 1)) - x(N_points_grid(i)), y(N_points_grid(i + 1)) - y(N_points_grid(i))]);
            
            if dist >= max_norm
               max_norm = dist;
               max_norm_idx = i;
            end    
            
        end
        
        if max_norm_idx == 1
            flag = 0;
        else
            N_points_grid(max_norm_idx) = N_points_grid(max_norm_idx) + 1;
        end  
    end
    
    figure('Position', [100 100 1000 1000]);
    hold on
    
    % равномерная параметризация по t
    dt = (t1 - t0) / (N - 1);
    T = zeros(1, N);
    for i = 1 : N
        T(i) = t0 + (i - 1) * dt;
    end
    plot(f(T), g(T), 'b*', 'MarkerSize', 10);
    
    
    plot(x, y);
    plot(f(t(N_points_grid)), g(t(N_points_grid)), 'r.','MarkerSize', 10);
    plot(f(t(1)), g(t(1)), 'r.','MarkerSize', 20);
    plot(f(t(t_point_cnt)) ,g(t(t_point_cnt)), 'r.', 'MarkerSize', 20);
    
    % for i = 1 : N
    %   text(f(t(N_points_grid(i))), g(t(N_points_grid(i))), num2str(i), 'FontSize', 10);
    %   text(f(T(i)), g(T(i)), num2str(i), 'Color', 'blue', 'FontSize', 12); 
    %end
  
    hold off
    
    % среднее расстояние между точками
    norm_N = 0;
    norm_T_uniform = 0;
    
    for i = 1 : N - 1
        norm_N = norm_N + ...
            norm([f(t(N_points_grid(i + 1))) - f(t(N_points_grid(i))), ...
            g(t(N_points_grid(i + 1))) - g(t(N_points_grid(i)))]);
        
        norm_T_uniform = norm_T_uniform + ...
            norm([f(T(i + 1)) - f(T(i)), g(T(i + 1)) - g(T(i))]);
    end
    
    disp('result: ');
    disp(norm_N / N);
    disp(norm_T_uniform / N);
end