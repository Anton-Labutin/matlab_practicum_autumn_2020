%%
% Лабораторная работа 2 
% Работа с графикой
% Вариант 2
%%
% Трёхмерная графика
% 
% Task 13
% 
% Колонисты установили на поверхности Марса N независимых антенных станций 
% (поля от них складываются), каждая из которых генерирует вокруг себя 
% беспроводную сеть с уровнем сигнала V / (1 + d(pk, p)), где pk — точка, 
% где находится k-я станция, p — точка, где проводится замер, d(·, ·) — 
% евклидово расстоние, V — исходный уровень сигнала. Для уверенной работы 
% марсохода требуется сигнал с уровнем, не меньшим L. Написать функцию 
% viewPossible(points, P, L), принимающую массив из координат N точек
% на плоскости, уровни сигналов L и P, и выводящую на экран область, в 
% которой можно уверенно управлять марсоходом. Определить, будет ли 
% полученная область односвязной.
%

points = [-1, -1; 5, 8; 13, 8; -4, 9; -2.5, 5; 0, -10; -5, -5; 10, -10];
V = 10;
L = 4;

viewPossible(points, V, L);

function viewPossible(points, V, L)
    x = linspace(min(points(:, 1)) - 2 * max(points(:, 1)), 3 * max(points(:, 1)), 1000);
    y = linspace(min(points(:, 2)) - 2 * max(points(:, 2)), 3 * max(points(:, 2)), 1000);
    [X, Y] = meshgrid(x, y);
    
    sygLevel = @(pk_x, pk_y, p_x, p_y) (ones(size(X, 1), size(X, 2)) .* V) ./ ...
        (1 + sqrt((pk_x - p_x) .^ 2 + (pk_y - p_y) .^ 2));
    sygLevelArea = zeros(size(X, 1), size(X, 2));
    
    for i = 1 : size(points, 1)
        sygLevelArea = sygLevelArea + ...
            sygLevel(X, Y, points(i, 1), points(i, 2));
    end
    
    M = contourf(X, Y, sygLevelArea, [L, L], 'blue');
    hold on;
    plot(points(:, 1), points(:, 2), '.', 'Color', 'red', 'MarkerSize', 15);
    hold off;
   
    
    if M(2, 1) + 1 == size(M, 2)
        title('Область односвязная');
    else
        title('Область не односвязная');
    end
end