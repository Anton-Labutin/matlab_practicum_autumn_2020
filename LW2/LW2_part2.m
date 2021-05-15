%%
% Лабораторная работа 2 
% Работа с графикой
% Вариант 2
%%
% Двумерная графика, часть вторая
%
% Task 4
% Написать функцию convergenceFunc(fn, f, a, b, n, convType)), принимающую 
% на вход аргументы: функцию fn(n, x), такую, что fn(n, x) = f_n(x), 
% и функцию f, считающуюся пределом последовательности fn(x) на [a,b] 
% в смысле, задаваемом аргументом-строкой convType: это может быть 
% поточечная сходимость, равномерная сходимость, среднеквадратичная 
% сходимость. Функция рисует анимацию из n кадров, на каждом i-м из 
% которых нарисованы fi и f на отрезке [a, b]. В заглавии графика стоит 
% вывести значении метрики разности для всех сходимостей, кроме поточечной. 
% Подготовить несколько примеров, когда есть один вид сходимости, 
% но нет другого.

% Поточечная и среднеквадратичная сходимость, но не равномерная 
f1 = @(x) 1 .* (x == 1) + 0 * (x < 1);
f1n = @(x, n) x .^ n;
a1 = 0;
b1 = 1;

% Равномерная сходимость (и поточечная, и среднеквадратичная)
f2 = @(x) 0;
f2n = @(x, n) sin(n * x) ./ n;
a2 = 0;
b2 = 2 * pi;

% Сходится поточечно, но не в среднем и не равномерно
f3 = @(x) 0;
f3n = @(x, n) n .* (x > 0 & x < 1 / n)
a3 = 0;
b3 = 1;

n = 100;

convergenceFunc(f2n, f2, a2, b2, n, 'uniform'); 
% pointwise - поточечная, uniform - равномерная, MS - среднеквадратичная
%%
% Task 5
% Написать функции fourierApprox(f,a,b, n, meth)), принимающаю на вход 
% аргументы: функцию f, и рисующую анимацию из n кадров, на каждом из 
% которых рисуется i-я частичная сумма ряда Фурье для этой функции на 
% [a, b] по системе функций, задаваемой параметром meth. Должна быть 
% реализована стандартная тригонометрчиеская система функций, система 
% многолченов Чебышёва и любая другая полная система функций, не 
% встречающаяся ни в одном из вариантов этого задания. Каждая система 
% должна порождаться функцией вида getFunc(n), возвращающую анонимную 
% функцию номер n в той или иной системе.

f = @(x) exp(x) + (x - 1) .^ 2

a_tfs = -pi;
b_tfs = pi;
n_tfs = 100;

a_cfs = -1;
b_cfs = 1;
n_cfs = 50;

a_lfs = -1;
b_lfs = 1;
n_lfs = 10;

fourierApprox (f, a_lfs, b_lfs, n_lfs, 'lfs');
% tfs = system of trigonometric functions
% cfs = function system by Chebyshev
% lfs = function system by Lezhandr

function [] = fourierApprox (f, a, b, n, meth)
    x = linspace(a, b, 1000);
    sum = zeros(1, size(x, 2)); % частичная сумма ряда Фурье
    
    y = zeros(n + 1, size(x, 2));
    
    for i = 1 : n
        if meth == 'tfs'
            left = -pi;
            right = pi;
            funcIdx = getFuncTFS(i);
            helpFunc = @(x) ones(1, size(x, 2));
        elseif meth == 'cfs'
            left = -1;
            right = 1;
            funcIdx = getFuncCFS(i);
            helpFunc = @(x) ones(1, size(x, 2)) ./ sqrt(1 - x .^ 2);
        elseif meth == 'lfs'
            left = -1;
            right = 1;
            func = getFuncLFS(i);
            funcIdx = @(x) func(x) ./ sqrt(2 / (2 * i - 1));
            helpFunc = @(x) ones(1, size(x, 2));
        end
        
        func_find_lyambda = @(x) (funcIdx(x) .^ 2) .* helpFunc(x);
        func_find_c = @(x) f(x) .* funcIdx(x) .* helpFunc(x);
        lyambdaIdx = integral(func_find_lyambda, left, right);
        cIdx = integral(func_find_c, left, right) / lyambdaIdx;
        sum = sum + funcIdx(x) .* cIdx;
        y(i, :) = sum;
    end
    
    y(n + 1, :) = f(x);
    min_y = min(min(y, [], 2));
    max_y = max(max(y, [], 2));
    
    mov(1 : n) = struct('cdata', [], 'colormap', []);
    
    for i = 1 : n
        plot(x, f(x), 'r');
        hold on;
        
        axis( [a b  min_y max_y] );
        plot(x, y(i, :), 'b');
        hold off;
        
        title(['Частичная сумма ряда Фурье по', newline, 'системе ', meth, ' с номером ', num2str(i)]);
        
        mov(i) = getframe();
        pause(0.1);
    end
end



function [f] =  getFuncTFS(i)
% 1 / sqrt(2 * pi), sin(x) / sqrt(2 * pi), cos(x) / (sqrt(2 * pi), ...
    if i == 1
        f = @(x) ones(1, size(x, 2)) / sqrt(2 * pi);
    elseif mod(i, 2) == 1
        f = @(x) sin(floor(i / 2) * x) / sqrt(pi);
    else
        f = @(x) cos((i / 2) * x) / sqrt(pi);
    end
end



function [f] =  getFuncCFS(i)
% система многочленов Чебышёва первого рода T_n(z) = cos(n * arccos(z))
    if i == 1
        f = @(x) cos((i - 1) * acos(x)) ./ sqrt(pi);
    else
        f = @(x) cos((i - 1) * acos(x)) ./ sqrt(pi / 2);
    end
end



function [f] =  getFuncLFS(i)
% система многочленов Лежандра
    if i == 1
        f = @(x) ones(1, size(x,2));
    elseif i == 2
        f = @(x) x;
    else
        lezh1 = getFuncLFS(i - 1);
        lezh2 = getFuncLFS(i - 2);
        f = @(x)  x .* (2 * i - 3) ./ (i - 1) .* lezh1(x) - lezh2(x) .* (i - 2) ./ (i - 1); 
    end
end

%%
function [] = convergenceFunc(fn, f, a, b, n, convType)
    x = linspace(a, b, 1000);
    y = ones(n, 1000);
    diff = ones(1, n);
    
    for seqIdx = 1:n 
        if strcmp(convType, 'pointwise')
            y(seqIdx, :) = fn(x, seqIdx);
        end

        if strcmp(convType, 'uniform')
            y(seqIdx, :) = fn(x, seqIdx);
            diff(seqIdx) = max(abs(f(x) - y(seqIdx, :)));
        end

        if strcmp(convType, 'MS')
            y(seqIdx, :) = fn(x, seqIdx);
            d = f(x) - y(seqIdx, :);
            d = d .^ 2;
            diff(seqIdx) = sqrt(trapz(x, d));
        end
    end
    
    frame(1:n) = struct('cdata', [], 'colormap', []);
    
    max_y = max(max(y, [], 2));
    min_y = min(min(y, [], 2));
    
    for seqIdx = 1 : n
        plot(x, f(x), 'r');
        hold on;
        
        if strcmp(convType, 'uniform')
            title( ['Значение метрики разности для равномерной сходимости: ', num2str(diff(seqIdx), '%.3f')] );
        end
        
        if strcmp(convType, 'MS')
            title( ['Значение метрики разности для среднеквадратичной сходимости: ', num2str(diff(seqIdx), '%.3f')] );
        end
        
        axis([a b min_y max_y]);
        plot(x, y(seqIdx, :), 'b');
        hold off;

        frame(seqIdx) = getframe();
        pause(0.01);
    end
end