% 
% Лабораторная работа 2
% Работа с графикой
% Вариант 2
%
%%
%
% Четырёхмерная графика
%
%%
% Task 14
%
% Написать функцию drawBall(alpha, params), 
% которая создаёт трехмерную сетку и рисует на ней линию уровня функции 
% f(x, y, z). В params (это может быть список параметров или структура) 
% передать параметры отрисовки (как минимум, цвет, диапазоны изменения 
% переменных и число точек в сетке). Если исследуемое множество пустое, 
% то функция должна выводить соответствующее сообщение об ошибке.
% 

params = struct();
params.x_min = -10;
params.x_max = -8;
params.y_min = -10;
params.y_max = 10;
params.z_min = -10;
params.z_max = 10;

params.face_color = 'red';
params.edge_color = 'blue';
params.face_alpha = 0.3;
params.line_style = '-';
params.marker = 'o';
params.marker_size = 3;
params.marker_face_color = 'green';

params.point_per_grid = 100;

drawBall(2, params);

%%
% Task 15
%
% Написать функцию drawManyBalls(alphas, colors, edges), которая 
% рисует единичные шары в метриках, указанных в векторе alphas, 
% с цветами, задаваемыми в векторе colors. Параметр edges отвечает 
% за цвет граней и может принимать значение ’None’

alphas = 1:7;

% colors = zeros(length(alphas), 3);
% edges = zeros(length(alphas), 3);
% for i = 1 : length(alphas)
%    colors(i, :) = [ i / length(alphas), (i - 0.25) / length(alphas), (i - 0.5) / length(alphas) ];
%    edges(i, :) = [ (i - 0.5) / length(alphas), (i - 0.75) / length(alphas), (i - 0.25) / length(alphas) ];
%end

colors = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];
edges = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];
drawManyBalls(alphas, colors, edges);

%%
function drawBall(alpha, params)
    if isinf(alpha)
        f = @(x, y, z) max(abs(x), max(abs(y), abs(z)));
    else
        f = @(x, y, z) (abs(x)) .^ alpha + (abs(y)) .^ alpha + (abs(z)) .^ alpha;
    end
    
    x = linspace(params.x_min, params.x_max, params.point_per_grid);
    y = linspace(params.y_min, params.y_max, params.point_per_grid);
    z = linspace(params.z_min, params.z_max, params.point_per_grid);
    [X, Y, Z] = meshgrid(x, y, z);
    
    V = f(X, Y, Z);
    fv = isosurface(X, Y, Z, V, 1);
    
    if isempty(fv.vertices) == true || isempty(fv.faces) == true
        error("fv is empty");
    end
    
    figure('Position', [100 100 1000 1000]);
    p = patch(fv);
    
    p.FaceColor = params.face_color;
    p.EdgeColor = params.edge_color; 
    p.FaceAlpha = params.face_alpha; 
    p.LineStyle = params.line_style;
    p.Marker = params.marker;
    p.MarkerFaceColor = params.marker_face_color; 
    p.MarkerSize = params.marker_size; 
    daspect([1 1 1]); 
    view(3); 
    axis tight; 
    camlight('right');
    lighting gouraud;
end


function drawManyBalls(alphas, colors, edges)
    x = linspace(-1, 1, 100);
    y = linspace(-1, 1, 100);
    z = linspace(-1, 1, 100);
    
    f = @(x, y, z, i) (abs(x).^alphas(i) + abs(y).^(alphas(i)) + abs(z).^(alphas(i))).^(1 / alphas(i));
    
    [X, Y, Z] = meshgrid(x, y, z);
    
    for i = 1 : length(alphas)
        V = f(X, Y, Z, i);
        fv = isosurface(X, Y, Z, V, 1);
        p = patch(fv, 'FaceColor', colors(i), 'FaceAlpha', (length(alphas) - i) / length(alphas));
        
        if ~strcmp(edges, 'None')
            p.EdgeColor = edges(i);
        end 
    end
    
    axis tight; 
    camlight('right');
    lighting gouraud;
    view(3); 
end