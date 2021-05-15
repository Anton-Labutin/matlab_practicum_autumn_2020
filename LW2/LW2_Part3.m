%%
% Лабораторная работа 2 
% Работа с графикой
% Вариант 2
%%
% Трёхмерная графика
% 
% Создать блочный скрипт. В первом блоке задается функция, зависяшая от 
% двух переменных и скалярного параметра, двухмерная сетка (см. команду 
% meshgrid) и границы изменения параметра.
%
%%
%
% Task 11
% Создать блок, рисующий анимацию с эволюцией поверхности по параметру 
% (см. surf) и сохраняющий анимацию в переменную. На каждом кадре 
% неоходимо отметить локальные максимумы и минимумы, подобрать примеры, 
% где их несколько и где они с течением времени меняют своё положение и 
% число. В следующем блоке воспроизвести эту анимацию командой movie. 
% Написать блок, где фиксируется некоторое значение параметра и при 
% помощи команды contour строится проекция сечения функции на некотором 
% фиксированном уровне на плоскость Oxy.
%
%%
%
[X, Y] = meshgrid(0:pi/20:(2 * pi), 0:pi/20:(2 * pi));

nFrames = 100;
mov(1:nFrames) = struct('cdata', [], 'colormap', []);



for i = 1:nFrames
    f(X, Y, (nFrames - i + 1) / nFrames);
    mov(i) = getframe();
end
%%
movie(mov);
%%
param = 0.5;
height = 0.5;
Z = param * cos(X) + (1 - param) * sin(Y);
contour(X, Y, Z, [height, height]);
%%
% Task 12
% Создать блоки, сохраняющие анимацию в файл на диске в форматах .mat и
% .avi
v = VideoWriter('animation.avi'); % 'Uncompressed AVI'
v.FrameRate = 10;
open(v);
writeVideo(v, mov);
close(v);
save('animation.mat', 'mov', '-nocompression');
load('animation.mat', '-mat');
%%
function f(X, Y, param)
    Z = param * cos(X) + (1 - param) * sin(Y);
    
    zLocalMin = Z .* islocalmin(Z .* islocalmin(Z));
    xLocalMin = X .* islocalmin(Z .* islocalmin(Z)); 
    yLocalMin = Y .* islocalmin(Z .* islocalmin(Z));
    
    zLocalMax = Z .* islocalmax(Z .* islocalmax(Z));
    xLocalMax = X .* islocalmax(Z .* islocalmax(Z)); 
    yLocalMax = Y .* islocalmax(Z .* islocalmax(Z));
    
    surf(X, Y, Z);
    hold on
        plot3(xLocalMin, yLocalMin, zLocalMin, 'g.', 'MarkerSize', 30);
        plot3(xLocalMax, yLocalMax, zLocalMax, 'r.', 'MarkerSize', 30);
    hold off
end