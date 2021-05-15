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
% Task 10
% 
% Написать функцию drawPolar(rho,N), которая получает на вход опорную функцию
% [support_value, support_point] = rho(x), возвращающую значение опорной функции некоторого 
% плоского множества X в направлении x и соответствующий опорный вектор. 
% Функция drawPolar рисует поляру множества X и само X . Подобрать 
% 3–4 примера, когда вид поляры известен заранее, в том числе, когда 
% 0 не ∈ X .

drawPolar(@drawPolar_withZero_ellipse, 100);

%@drawPolar_withZero_square  с 0
%@drawPolar_withoutZero_square  без 0
%@drawPolar_withZero_ellipse без 0
%@drawPolar_withZero_rhomb с 0
%@drawPolar_withoutZero_rhomb без 0 

function result = drawPolar(rho, N)
    t = linspace(0, 2 * pi * (N - 1) / N, N);
    y1 = cos(t);
    y2 = sin(t);
    
    dir_vec = [y1(1), y2(1)]; % [1, 0]
    [support_value, support_point] = rho(dir_vec);
    point_border = dir_vec ./ abs(support_value);
    
    point1 = point_border;
    
    polar_x = zeros(1, N);
    polar_y = zeros(1, N);
    polar_x(1) = point_border(1); % 1
    polar_y(1) = point_border(2); % 0
    
    hold on
    
    for i = 2 : N
        dir_vec = [y1(i), y2(i)];
        [support_value, support_point] = rho(dir_vec);
        point2 = dir_vec ./ abs(support_value);
        polar_x(i) = point2(1);
        polar_y(i) = point2(2);
    end
   
    K = convhull(polar_x, polar_y); % convex hull of polar
    patch(polar_x(K),polar_y(K),'green');
    drawSet(rho, 50);
    axis equal;
    hold off
end


function [] = drawSet(rho, N)
    for i = 0 : N - 1
        angle = (2 * pi / N) * i;
        next_angle = (2 * pi / N) * (i + 1);
        vector_i(1) = cos(angle);
        vector_i(2) = sin(angle);
        vector_i_plus_1(1) = cos(next_angle);
        vector_i_plus_1(2) = sin(next_angle);
        
        [supp_func_i_val, supp_func_i_point] = rho(vector_i);
        
        [supp_func_i_plus_1_val, supp_func_i_plus_1_point] = rho(vector_i_plus_1);
     
        hold on;
        % external approximation
        
        % solve the linear system (vec - intersection of hyperplanes)
        % <vec, vector_i> = supp_func_i_val and <vec, vector_i_plus_1> =
        % supp_func_i_plus_1_val
        
        det = vector_i(1) * vector_i_plus_1(2) - vector_i_plus_1(1) * vector_i(2);
        if det ~= 0
            intersection_x = (supp_func_i_val * vector_i_plus_1(2) - supp_func_i_plus_1_val * vector_i(2)) / det;
            intersection_y = (vector_i(1) * supp_func_i_plus_1_val - vector_i_plus_1(1) * supp_func_i_val) / det;
            
            line([supp_func_i_point(1), intersection_x], [supp_func_i_point(2), intersection_y], 'Color', 'blue');
            line([intersection_x, supp_func_i_plus_1_point(1)], [intersection_y, supp_func_i_plus_1_point(2)], 'Color', 'blue');         
        else
            line([ supp_func_i_point(1), supp_func_i_plus_1_point(1) ], [ supp_func_i_point(2), supp_func_i_plus_1_point(2) ], 'Color', 'blue');
        end
        
        % inner approximation
        line([supp_func_i_point(1), supp_func_i_plus_1_point(1)], [supp_func_i_point(2), supp_func_i_plus_1_point(2)], 'Color', 'red');
    end
    hold off;
end


function [support_value, support_point] = supportDrawPolar(f, opts)
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    
    scalar_product = @(x, l) (x(1) * l(1) + x(2) * l(2));
    
    function [c, ceq] = mycon(x)
        ceq = f(x);
        c = [];
    end

    nonlcon = @mycon;
    
    sup_x = @(l) fmincon(@(x) -scalar_product(x, l), ...
        opts.x0, A, b, Aeq, beq, lb, ub, ...
        @mycon);
     
    support_point = sup_x([opts.l1, opts.l2]);
    support_value = scalar_product(support_point, [opts.l1, opts.l2]);
end


function [support_value, support_point] = drawPolar_withZero_ellipse(y)
    f = @(x) (6 * x(1))^2 + (3 * x(2))^2 - 1;
    opts = struct('l1', y(1), 'l2', y(2), 'x0', [3, 0.25]);
    [support_value, support_point] = supportDrawPolar(f, opts);
end

function [support_value, support_point] = drawPolar_withZero_square(y)
    f = @(x) 2 - max(abs(x(1)), abs(x(2)));
    opts = struct('l1', y(1), 'l2', y(2), 'x0', [0, 0]);
    [support_value, support_point] = supportDrawPolar(f, opts);
end

function [support_value, support_point] = drawPolar_withoutZero_square(y)
    f = @(x) 1 - max(abs(x(1) + 2), abs(x(2) - 2));
    opts = struct('l1', y(1), 'l2', y(2), 'x0', [-2, 2]);
    [support_value, support_point] = supportDrawPolar(f, opts);
end

function [support_value, support_point] = drawPolar_withZero_rhomb(y)
    f = @(x) 1 - abs(x(1)) - abs(x(2));
    opts = struct('l1',y(1),'l2',y(2),'x0',[0, 0]);
    [support_value, support_point] = supportDrawPolar(f,opts);
end

function [support_value, support_point] = drawPolar_withoutZero_rhomb(y)
    f = @(x) 1 - abs(x(1) - 2) - abs(x(2));
    opts = struct('l1', y(1), 'l2', y(2), 'x0', [2, 0]);
    [support_value, support_point] = supportDrawPolar(f,opts);
end