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
% Task 8 
% 
% Написать функцию drawSet(rho, N) которая получает на вход 
% опорную функцию некоторого плоского множества [val, supp_func_i_point] = rho(x),
% возвращающую значение опорной функции в направлении x и соответствующий 
% опорный вектор. Функция drawSet рисует внутреннюю и внешнюю 
% кусочно-линейные аппроксимации границы множества с N точками. 
% Функцией convhull пользоваться нельзя. Подготовить 3 примера с 
% аналитически рассчитанными опорными функциями: эллипс, квадрат, ромб 
% (в каждом случае центр не обязательно нулевой; центр и полуоси являются 
% параметрами).
% При помощи системы Latex в отдельном pdf-файле выписать вывод опорных 
% функций для трех указанных выше множеств.
%

% Центр ромба
rhomb_c = [2, 1];
% первая диагональ
rhomb_diag1 = 10;
% вторая диагональ
rhomb_diag2 = 5;
% Ромб:
rhomb_func = @(x) ( 2 * abs(x(1) - rhomb_c(1)) / rhomb_diag1) + (2 * abs(x(2) - rhomb_c(2)) / rhomb_diag2) - 1;
% Функция rho для ромба:
rhomb_supp_func = @(l) [ ...
    l(1) * rhomb_c(1) + l(2) * rhomb_c(2) + max( (rhomb_diag1 / 2) * abs(l(1)), (rhomb_diag2 / 2) * abs(l(2)) ), ... % опорная функция
    rhomb_c(1) + ((rhomb_diag1 / 2) * abs(l(1)) >= (rhomb_diag2 / 2) * abs(l(2))) * (rhomb_diag1 / 2) * sign(l(1)), ... % 1я координата опорного вектора
    rhomb_c(2) + ((rhomb_diag1 / 2) * abs(l(1)) < (rhomb_diag2 / 2) * abs(l(2))) * (rhomb_diag2 / 2) * sign(l(2)) ... % 2я координата опорного вектора
    ];

% Центр квадрата
square_c = [1, 1];
% Сторона квадрата
square_side = 8;
% Квадрат:
square_func = @(x) max(abs(x(1) - square_c(1)), abs(x(2) - square_c(2))) - square_side / 2;
% Функция rho для квадрата:
square_supp_func = @(l) [ ...
    l(1) * square_c(1) + l(2) * square_c(2) + (square_side / 2) * (abs(l(1)) + abs(l(2))), ... 
    square_c(1) + (l(1) > 0) * square_side / 2 + (l(1) <= 0) * (-square_side / 2), ... 
    square_c(2) + (l(2) > 0) * square_side / 2 + (l(2) <= 0) * (-square_side / 2) ...
    ];

% Центр эллипса
ellipse_c = [2, 3];
% Полуоси эллипса
a = 7;
b = 4;
% Эллипс:
ellipse_func = @(x) (x(1) - ellipse_c(1)) ^ 2 / a ^ 2 + (x(2) - ellipse_c(2)) ^ 2 / b ^ 2 - 1;
% Функция rho для эллипса
ellipse_supp_func = @(l) [ ...
    l(1) * ellipse_c(1) + l(2) * ellipse_c(2) + sqrt( (l(1) * a) ^ 2 + (l(2) * b) ^ 2  ), ...
    ellipse_c(1) + l(1) * (a ^ 2) / sqrt((l(1) ^ 2) * (a ^ 2) + (l(2) ^ 2) * (b ^ 2)), ...
    ellipse_c(2) + l(2) * (b ^ 2) / sqrt((l(1) ^ 2) * (a ^ 2) + (l(2) ^ 2) * (b ^ 2)) ... 
    ];


drawSet(ellipse_supp_func, 30)


function [] = drawSet(rho, N)
    for i = 0 : N - 1
        angle = (2 * pi / N) * i;
        next_angle = (2 * pi / N) * (i + 1);
        vector_i(1) = cos(angle);
        vector_i(2) = sin(angle);
        vector_i_plus_1(1) = cos(next_angle);
        vector_i_plus_1(2) = sin(next_angle);
        
        supp_func_i = rho(vector_i);
        supp_func_i_val = supp_func_i(1);
        supp_func_i_point = supp_func_i(2:3);
        
        supp_func_i_plus_1 = rho(vector_i_plus_1);
        supp_func_i_plus_1_val = supp_func_i_plus_1(1);
        supp_func_i_plus_1_point = supp_func_i_plus_1(2:3); 
        
        segment_left = -10;
        segment_right = 10;
        step = 0.01;
        grid = segment_left : step : segment_right;
     
        hold on;
        axis([segment_left segment_right segment_left segment_right]);
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