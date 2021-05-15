%%
% Лабораторная работа 2 
% Работа с графикой
% Вариант 2
%%
% Двумерная графика, часть вторая
%%
% Task 6
% Создать скрипт, в первом блоке которого задаётся переменная-функция 
% и некоторая сетка. Второй блок рисует график функции на этой сетке, 
% отмечает на ней все точки локального минимума, отмечает один глобальный 
% максимум, и запускает от него до ближайшего минимума комету (команда 
% comet). Подобрать примеры функций со многими экстремумами. В разных 
% примерах комета должна иметь возможность двигаться как вправо, так 
% и влево.
%%
grid = linspace(0, 6 * pi, 100);
f = @(x) sin(x);
%f = @(x) exp(-x.^2) + sin(x);
%grid = linspace(-20, 20, 1000);
y = f(grid);
%%
f_local_min = y(islocalmin(y));
x_local_min = grid(islocalmin(y));

[f_global_max, global_max_idx] = max(y);

figure('Position', [100 100 1000 1000]);
plot(grid, y, x_local_min, f_local_min, '*', grid(global_max_idx), f_global_max, '*', 'MarkerSize', 30);
hold on

local_min_idx = find(islocalmin(y) == 1);
[tmp, nearest_local_min_idx] = min(abs(grid(local_min_idx) - grid(global_max_idx)));

if local_min_idx(nearest_local_min_idx) < global_max_idx
    comet_x = grid(local_min_idx(nearest_local_min_idx) : global_max_idx);
    comet_y = f(comet_x);
    comet(fliplr(comet_x), fliplr(comet_y));
else
    comet_x = grid(global_max_idx : local_min_idx(nearest_local_min_idx));
    comet_y = f(comet_x);
    comet(comet_x, comet_y);
end