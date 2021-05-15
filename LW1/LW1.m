%% Task 1 
%
% Задать два вещественных числа (a и b), натуральное число n и равномерную сетку на [a, b] с n точками.
% Задать функцию f(x) = sin(ln(1+|x|)+x^2)
% Нарисовать график её значений на сетке, отметить отдельно максимальные и минимальные значения.
%
% Input data
a = input('Input a: ');
b = input('Input b: ');
n = input('Input n: ');

% function 
grid = linspace(a, b, n);
val = sin(log(1 + abs(grid)) + grid .^ 2);

% find local maximums and minimums
min_idx = islocalmin(val);
max_idx = islocalmax(val);

% plot
plot(grid, val, grid(min_idx), val(min_idx), '*', grid(max_idx), val(max_idx), '*');
%
%
%% Task 2
%
% Запросить у пользователя ввод числа n. Проверить, что введенное число —
% натуральное. 
% 1. Создать вектор из всех нечетных чисел, делящихся на 9, из промежутка от 1 до n. 
% 2. Построить матрицу размера n × n, все элементы i–й строки которой
% равны i.
% 3. Создать матрицу B n × (n + 1) вида: B = [[1, 2, 3]; [4, 5, 6]]
% Вытянуть матрицу B в вектор c. Присвоить переменной D последние 2 столбца матрицы B.
%
% Input data
n = input('Enter number n: ');

% Check the number
if n > 0 & floor(n) == ceil(n)
    'n is integer'
else
    'n is not integer'
end

% create vector

% i = 1;
% for x = 1:n
%    if rem(x, 2) == 1 & rem(x, 9) == 0
%        vec(i) = x;
%        i = i + 1;
%    end
% end 

vec = 9:18:n

% create matrix A

% A = ones(n, n);
%for i = 2:n
%    A(i,:) = i;
%end

A = (1:n)';
A = repmat(A, 1, n)
%%
% create matrix B

b = 1 : n * (n + 1);
B = reshape(b, n + 1, n)'

% transform natrix to vector

C = B';
c = C(:)

% create D
D = B(:, n:(n + 1))
%
%
%% Task 3
% 
% Создать матрицу размера 7 × 7, состоящую из случайных элементов 
% с равномерным распределением среди натуральных чисел от 1 до 315, 
% найти максимальный элемент на диагонали этой матрицы, 
% найти максимальное и минимальное отношение произведения к сумме 
% для строк этой матрицы, отсортировать строки матрицы 
% в лексико-графическом порядке (то есть строка [a1, a2, a3, . . . , an] 
% стоит в матрице выше строки [b1, b2, b3, . . . , bn], если ai = bi при
% i = 1,...,k−1 и ak < bk для некоторого k).
%
A = floor(1 + (316 - 1) .* rand(7))
diagMax = max(diag(A))
maxRel = max(prod(A, 2) ./ sum(A, 2))
minRel = min(prod(A, 2) ./ sum(A, 2))
sortrows(A)
%
%% Task 4
%
% Построить таблицу умножения всевозможных пар элементов таких, 
% что первый — элемент вектора X, а второй — вектора Y.
%
% input
n = input('Input n: ');
m = input('Input m: ');
% generate vectors
x = 1 + (10 - 1) * rand(n, 1)
y = 1 + (10 - 1) * rand(m, 1)

kronProd = kron(y', x) 
% 2й способ
% y' .* x
%
%% Task 5
%
% Запросить у пользователя ввод числа n. Проверить, что введенное число — 
% простое. Создать случайную матрицу A ∈ R(n×n) и вектор b ∈ R(n×1), 
% в случае, если A не вырождена, решить уравнение Ax = b (решить задачу 
% не менее чем двумя способами и вставить проверку возможности решения и 
% правильности решения). 
%
% input
n = input('Input n: ');
% check primary
if isprime(n)
    'n is prime'
else 
    'n is not prime'
end
% create matrices and vector
A = rand(n)
b = rand(n, 1)
if abs(det(A)) > 10 ^ (-6)
    % 1st method
    x1 = inv(A) * b
    if abs(A * x1 - b) < 10 ^ (-6)
        'x1 is the solution of the linear system A * x = b'
    end
    % 2nd method: Gauss' method
    x2 = A \ b
    if abs(A * x2 - b) < 10 ^ (-6)
        'x2 is the solution of the linear system A * x = b'
    end
    % 3d method: LU
    [L, U] = lu(A);
    y = L \ b;
    x3 = U \ y
    if abs(A * x3 - b) < 10 ^ (-6)
        'x3 is the solution of the linear system A * x = b'
    end
    % 4th method: QR
    [Q, R] = qr(A);
    x4 = R \ (Q' * b)
    if abs(A * x4 - b) < 10 ^ (-6)
        'x4 is the solution of the linear system A * x = b'
    end
    % 5th method
    x5 = linsolve(A, b)
    if abs(A * x5 - b) < 10 ^ (-6)
        'x5 is the solution of the linear system A * x = b'
    end
else 
    'det(A) = 0'
end
%
%
%% Task 6
% 
% Даны вектора a размерности n и b размерности m. Найти, используя только 
% арифметические операции и команды max и min, максимум функции |ai − bj|, 
% где ai — элемент вектора a, bj — элемент вектора b. Функцию abs и 
% дополнительную память не использовать.
%
% input 
n = input('Input n: ');
m = input('Input m: ');
% create vectors
a = 1 + (10 - 1) * rand(n, 1)
b = 1 + (10 - 1) * rand(m, 1)
% calc maximum of |a_i - b_j|
maxVal = max( max(a) - min(b), max(b) - min(a) );
%
%% Task 7
%
% Пусть у нас задано n точек в пространстве Rk в виде матрицы double[n,k]. 
% Требуется построить матрицу double[n,n] расстояний между каждой парой 
% точек. Пользоваться командами pdist и squareform нельзя.
%
% input
n = input('Input number of points: ');
k = input('Input dimension of the space: ');
% create points
points = 1 + (10 - 1) * rand(n, k)
% calculate pairwise distances

A_i_times_j = points * points'; % point_i * point_j, n * n
temp_i_squared = diag(A_i_times_j); % [...point_i ^ 2...]
A_i_squared_in_row = repmat(temp_i_squared, 1, n); % [[...], ..., [point_i ^ 2], ..., [...]], n * n
A_i_squared_in_col = A_i_squared_in_row'; % [[...]', ..., [point_i ^ 2]', ..., [...]']

pwDist = sqrt(A_i_squared_in_row - 2 * A_i_times_j + A_i_squared_in_col)

% pwDist = zeros(n, n);
% for i = 1:n
%     for j = (i + 1):n
%         for l = 1:k
%             pwDist(i, j) = pwDist(i, j) + (points(l, i) - points(l, j)) ^ 2;
%         end
%         pwDist(i, j) = sqrt(pwDist(i, j));
%    end
% end
% pwDist = pwDist + pwDist'
%


%% Task 8
% 
% Построить матрицу, в которой по строкам записаны все n-мерные бинарные 
% векторы. Натуральное число n задается пользователем.
% 
% input 
n = input('Input n: ');
% output all n-digit binary vectors

% binNumbers = zeros(2 ^ n, n);
% for num = 1:(2 ^ n - 1)
%     binNumbers(num + 1, :) = de2bi(num, n, 'left-msb');
% end

numbers = (0:n)';
binNumbers = de2bi(numbers, 'left-msb')
%
%% Task 9
% 
% Реализовать функцию C = my_multiply(A, B), которая выполняет расчет 
% значения C = AB по определению («строка на столбец»). Сравнить 
% быстродействие этой функции и стандартного умножения матриц для матриц 
% различной размерности. Построить график времени работы в зависимости от 
% размера матриц (в случае квадратных матриц).
%
n = input('Input maximum matrix size: ');
repeatCnt = input('Input count of repeatitions for each matrix size: ');
% calc times
execTime1 = zeros(n);
execTime2 = zeros(n);
totalTime1 = 0;
totalTime2 = 0;
for matrSz = 1:n
    A = 1 + (100 - 1) * randn(matrSz);
    B = 1 + (100 - 1) * randn(matrSz);
    
    tic
    for i = 1:repeatCnt
        C1 = A * B;
    end
    totalTime1 = toc;
    
    tic
    for i = 1:repeatCnt
        C2 = my_multiply(A, B); % function is defined at the end of the script!
    end
    totalTime2 = toc;
    
    execTime1(matrSz) = totalTime1 / repeatCnt;
    execTime2(matrSz) = totalTime2 / repeatCnt;
    totalTime1 = 0;
    totalTime2 = 0;
end
% plot
plot(1:n, execTime1, 'g');
legend({'A * B', 'my multiply'});
plot(1:n, execTime2, 'r');
legend({'A * B', 'my multiply'});
xlabel('matrix size');
ylabel('execution time');
%
%
%% Task 10
%
% Напишите функцию, которая находит средние значения (по одному 
% направлению) с учётом NaN элементов матрицы. Команду nanmean использовать
% нельзя.
%
A = [[NaN, 1, 2]; [NaN, 0, 6]; [1, 5, NaN]]
vec = calc_mean(A) % function is defined at the end of the script!
%
%
%% Task 11
%
% Сгенерировать вектор из n случайных величин с нормальным распределением 
% N(a, σ2). Проверить «правило трёх сигм»: вывести долю элементов вектора, 
% находящихся в интервале [a − 3σ, a + 3σ].
%
clear
n = input('Input the size of the vector: ');
a = input('Input mathematical experience: ');
d = input('Input dispersion: ');
vec = normrnd(a * ones(1, n), sqrt(d));

n_part = sum(vec >= a - 3 * sqrt(d) & vec <= a + 3 * sqrt(d)) * 100 / n
%
%% Task 12
%
% По аналогии с функцией trapz реализовать аналогичные функции rectangles 
% (интегрирование методом прямоугольников) и simpson (методом Симпсона). 
% С помощью них построить график первообразной функции f(x) = = sin(x) / x. 
% Сравнить внутреннюю скорость сходимости при использовании всех трёх 
% методов (внутренняя скорость сходимости определяется с помощью сравнения 
% разностей решений при шаге h и h/2, нарисовать график этой ошибки
% в зависимости от h). Сравнить время вычисления.
%
% input
a = input('Input a: ');
b = input('Input b: ');
pointCnt = input('Input count of points: ');
% func
step = (b - a) / pointCnt;
X = a:step:b;
f = @(x) sin(x) ./ x;
% antideriv
y_rect = zeros(1, pointCnt + 1);
y_trapz = zeros(1, pointCnt + 1);
y_simp = zeros(1, pointCnt + 1);

for ind = 1:pointCnt
    x = a:step:(a + ind * step);
    y = f(x);
    
    y_rect(1, ind) = rectangles(x, y);
    y_trapz(1, ind) = trapz(x, y);
    y_simp(1, ind) = simpson(x, y);
end

figure(1);
plot(X, y_rect, 'r', X, y_trapz, 'g', X, y_simp, 'b')
title("Antiderivative");
legend({'rectangles', 'trapz', 'simpson'});
xlabel('X');
ylabel('Y');

h = zeros(1, pointCnt);
antideriv_trap_h = zeros(1, pointCnt);
antideriv_rect_h = zeros(1, pointCnt);
antideriv_simp_h = zeros(1, pointCnt);
antideriv_trap_h_half = zeros(1, pointCnt);
antideriv_rect_h_half = zeros(1, pointCnt);
antideriv_simp_h_half = zeros(1, pointCnt);

% convergence rate
for ind = 1:pointCnt
    h(ind) = (b - a) / (2 * ind);
    x1 = a:h(ind):b;
    y1 = f(x1);
    
    antideriv_trap_h(ind) = trapz(x1, y1);
    antideriv_rect_h(ind) = rectangles(x1, y1);
    antideriv_simp_h(ind) = simpson(x1, y1);
    
    h_half = h(ind) / 2;
    x2 = a:h_half:b;
    y2 = f(x2);
    
    antideriv_trap_h_half(ind) = trapz(x2, y2);
    antideriv_rect_h_half(ind) = rectangles(x2, y2);
    antideriv_simp_h_half(ind) = simpson(x2, y2);
end

antideriv_trap = zeros(1, pointCnt);
antideriv_rect = zeros(1, pointCnt);
antideriv_simp = zeros(1, pointCnt);

antideriv_trap(:) = abs(antideriv_trap_h(:) - antideriv_trap_h_half(:));
antideriv_rect(:) = abs(antideriv_rect_h(:) - antideriv_rect_h_half(:));
antideriv_simp(:) = abs(antideriv_simp_h(:) - antideriv_simp_h_half(:));

figure(2);
plot(h, antideriv_rect, 'r', h, antideriv_trap, 'g', h, antideriv_simp, 'b')
title('Convergence rate');
legend({'rectangles', 'trapz', 'simpson'});
xlabel('h');
ylabel('conv');
%
%% Task 13
%
% Задать формулу для некоторой функции и её производной. На одном графике 
% в логарифмическом масштабе (loglog) вывести модули разностей между 
% точным значением производной в некоторой точке и правой и центральной 
% разностной производной в зависимости от шага численного дифференцирования.
%
x = pi / 4;
y = func(x);
% derivations
f_deriv = func_deriv(x);
h = logspace(-10,-1, 100);
f_right_deriv = (func(x + h) - func(x)) ./ h;
f_central_deriv = (func(x + h) - func(x - h)) ./ (2 * h);
% plots
loglog(h, abs(f_deriv - f_right_deriv), 'r', h, abs(f_deriv - f_central_deriv), 'g')
legend({'right', 'central'});
clear;
%
%% Functions
% Task 9
function C = my_multiply(A, B)
    sizesA = size(A);
    sizesB = size(B);
    C = zeros(sizesA(1), sizesB(2));
    for i = 1:sizesA(1)
        for j = 1:sizesB(2)
             C(i, j) = A(i, :) * B(:, j);
             % for k = 1:sizesA(2)
             %     C(i, j) = C(i, j) + A(i, k) * B(k, j);
             % end
        end
    end
end

% Task 10
function colMeanVal = calc_mean(A)
    sizesA = size(A);
    colMeanVal = zeros(1, sizesA(2));
    colMeanVal = sum(A, 'omitnan') ./ (sizesA(2) - sum(isnan(A)));
end
%

% Task 12
function I = rectangles(x, y)
    I = sum(y) * (x(2) - x(1));
end

function I = simpson(x, y)
    segCnt = size(x, 2) - 1; 
    segStart = x(1); 
    segEnd = x(segCnt); 
    step = (segEnd - segStart) / segCnt; 

    sum_1 = y(1 : 2 : segCnt - 2); 
    sum_2 = y(2 : 2 : segCnt - 1); 
    sum_3 = y(3 : 2 : segCnt); 

    I = (step / 3) * sum( sum_1(:) + 4 * sum_2(:) + sum_3(:) ); 
end
%

% Task 13
function y = func(x, y)
    y = sin(x);
end

function y = func_deriv(x)
    y = cos(x);
end
%