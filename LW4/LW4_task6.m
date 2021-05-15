%
% Лабораторная работа 4
% Численные методы (продолжение). Интеграция с языками С/С++.
% Вариант 2
%
% Задание 6
% Дана краевая задача. Для неё рассматривается разностная схема. 
% Реализовать численный метод и подобрать примеры. Проверить корректность 
% работы численного алгоритма.

%% Основной пример
M = 40;
N = 30;

fGiven = @(x,y) (5 + x.^2) .* exp(3 * x) - 2 * y .* sin(5 * y);

u1Zero = 1;
u2Zero = 1;
mu1 = 1;
val = uNumerical(u1Zero, u2Zero, mu1, M, N);

x = linspace(0, 1 - 1 / M, M);
y = linspace(0, 1 - 1 / N, N);
[X, Y] =  meshgrid(y,x);

val_analyt = uAnalytical(X, Y, u1Zero, u2Zero, mu1);

figure('Position', [100 100 1000 1000]);
hold on;
surf(X, Y, real(val), 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
surf(X, Y, val_analyt, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none');
legend("Численное решение", "Аналитическое решение");
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
view(3);

figure('Position', [100 100 1000 1000]);
hold on;
surf(X, Y, abs(val_analyt - real(val)), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
zlabel('z');
legend("Модуль разности аналитического и численного решений");
hold off;
view(3);

max(abs(val_analyt - real(val)), [], 'all')

%% Пример 1: решение x * (x - 1) + y * (y^3 - 1)
M = 30;
N = 40;

mu3 = 3;
fHandle3 = @(x,y) 2 + 12 * y.^2 - mu3 * (x .* (x - 1) + y .* (y.^3 - 1));
xiHandle3 = @(x) x .* (x - 1);
etaHandle3 = @(y) y .* (y.^3 - 1);

% решение
val3 = solveDirichlet(fHandle3,xiHandle3,etaHandle3,mu3,M,N);

x = linspace(0, 1 - 1 / M, M);
y = linspace(0, 1 - 1 / N, N);
[X, Y] =  meshgrid(y, x);

f_1 = @(x,y) x .* (x - 1) + y .* (y.^3 - 1);
val_f1 = f_1(X, Y);

figure('Position', [100 100 1000 1000]);
hold on;
surf(X, Y, real(val3), 'FaceAlpha', 0.5, 'FaceColor', 'b');
surf(X, Y, val_f1, 'FaceAlpha', 0.5, 'FaceColor', 'r');
legend('My solve','Real solve');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
view(3);

figure('Position', [100 100 1000 1000]);
surf(X,Y,abs(real(val_f1) - val3));
legend('abs (My solve-Real solve)');
xlabel('x');
ylabel('y');
zlabel('z');
hold off
view(3);

disp('Максимальная погрешность');
disp(max(max(abs(real(val3)-val_f1))))
%% Пример 2: решение есть 3 * sin(pi * x) + sin(pi * y)
M = 30;
N = 35;

mu1 = 1;
fHandle2 = @(x,y) -(pi^2 + mu1) * (3 * sin(pi * x) + sin(pi * y));
xiHandle2 = @(x) 3 * sin(pi * x);
etaHandle2 = @(y) sin(pi * y);

% решение
val = solveDirichlet(fHandle2, xiHandle2, etaHandle2, mu1, M, N);

x = linspace(0, 1 - 1 / M, M);
y = linspace(0, 1 - 1 / N, N);
[X, Y] =  meshgrid(y, x);

f_1 = @(x, y) 3 * sin(pi * x) + sin(y * pi);
val_f1 = f_1(X, Y);

figure('Position', [100 100 1000 1000]);
hold on;
surf(X, Y, real(val), 'FaceAlpha', 0.5, 'FaceColor', 'b');
surf(X, Y, val_f1, 'FaceAlpha', 0.5, 'FaceColor', 'r');
legend('My solve','Real solve');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
view(3);

figure('Position', [100 100 1000 1000]);
hold on;
surf(X, Y, abs(real(val_f1) - (val)));
legend('abs (My solve-Real solve)');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
view(3);

disp('Максимальная погрешность');
disp(max(max(abs(real(val) - val_f1))));


% решение КЗ с помощью разностной схемы
function res =  solveDirichlet(fHandle, ... 
    xiHandle, ... // краевое условие по x
    etaHandle, ... // краевое условие по y
    mu, ... 
    M, ... 
    N)

    res = zeros(M, N);
    
    h_x = 1 / M;
    h_y = 1 / N;
    x = linspace(0, 1 - h_x, M);
    y = linspace(0, 1 - h_y, N);
    [X, Y] =  meshgrid(y, x);
    
    f = zeros(M, N); % значения функции на двумерной сетке
    f = fHandle(X, Y);

    f(:, 1) = 0;
    f(1, :) = 0;

    xi = zeros(1, N); % значения краевого условия по x на сетке по y
    xi = xiHandle(y);
    
    eta = zeros(1, M); % значения краевого условия по y на сетке по x
    eta = etaHandle(x);
    
    alpha = ifft(xi);
    beta = ifft(eta);

    C = ones(M, N);
    C = -mu * C;
    for p = 1 : M
        C(p, 1 : N) = C(p, 1 : N) - (4 / (h_x^2)) * (sin(pi * (p - 1) / M)).^2;
        C(p, 1 : N) = C(p, 1 : N) - (4 / (h_y^2)) * (sin(pi * ((1 : N) - 1) / N)).^2;
    end
    
    D_mat = ifft2(f);
    
    Q = zeros(M + N - 1, M + N - 1);
    t = zeros(1 , M + N - 1);
    
    % заполняем первыми M числами
    for p = 1 : M
        % Вычислем значение правой части
        D = sum((1 ./ C(p, 1 : N)) .* D_mat(p, 1 : N), 2);
        
        t(p) = beta(p) - D; 
        
        inv_c_vec = zeros(1, N);
        inv_c_vec(1 : N) = 1 ./ C(p, 1 : N);
        
        A = ifft(inv_c_vec) / M;
        % первый элемент заполняем отдельно
        Q(p, 2 : N) = A(2 : N);
        
        % Заполняем первые N элементов p-ого столбца
        Q(p, 1) = Q(p, 1) + sum(inv_c_vec) ./ (N * M);
        
        % Вычисляем (N + 1)( ..., (M + N - 1) коэффициенты системы
        Q(p, N + 1 : M + N - 1) = (sum(inv_c_vec) ./ (N * M)) .* ... 
            exp((2 * pi * 1j * (p - 1) * (1 : M - 1)) / M);
    end
    
    % Заполняем вторую часть матрицы
    for q = 2 : N
        D = sum((1 ./ C(1 : M, q)) .* D_mat(1 : M, q), 1);
        
        t(q + M - 1) = alpha(q) - D; 
        
        % вычисляем правую часть
        inv_c_vec = 1 ./ C(1 : M, q);
        
        %вычисляем вспомогательные значения
        Q(q + M - 1, 1) = Q(q + M - 1, 1) + sum(inv_c_vec) ./ (N * M);
        Q(q + M - 1, 2 : N) = (sum(inv_c_vec) ./ (N * M)) .* ...
                 exp((2 * pi * 1j * (q - 1) * (1 : N - 1)) / N);
        
        % Заполняем первые N элементов
        B = ifft(inv_c_vec) / N;
       
        Q(q + M - 1, N + 1 : N + M - 1) = B(2 : M);
    end
    
    % Решаем систему
    eps = 1e-7;
    max_it_cnt = 1000;
    f_NM = bicg(Q, t', eps, max_it_cnt);
    
    % Заполняем матрицу f полученнымим значениями
    f(1, 1 : N) = f_NM(1 : N);
    f(2 : M, 1) = f_NM(N + 1 : N + M - 1);
    
    b_pq = ifft2(f);

    a_pq = zeros (M, N);
    a_pq = b_pq ./ C;
    
    res = real(fft2(a_pq));
end

% аналитическое решение КЗ с правой частью f_1
function res = uAnalytical(xMat, ... % матрица x
    yMat, ... % матрица y
    u1Zero, ... % u_10
    u2Zero, ... % u2_0
    mu)

    v1 = @(z) exp(-(mu.^(1/2)) * z) - exp((mu.^(1/2)) * z);
    v2 = @(z) exp((2 - z) * mu.^(1/2)) - exp((mu.^(1/2)) * z);
    W = 2 * (mu.^(1/2)) * (exp(2 * mu.^(1/2)) - 1);
    
    I1 = v1fx;
    I11 = I1(mu, u1Zero, xMat);
    I11 = v2(xMat) .* I11 ./ W;
    
    I2 = v2fx;
    I21 = I2(mu, u1Zero, xMat);
    I21 = v1(xMat) .* I21 ./ W;
    
    I3 = v1fy;
    I31 = I3(mu, u2Zero, yMat);
    I31 = v2(yMat) .* I31 ./ W;
    
    I4 = v2fy;
    I41 = I4(mu, u2Zero, yMat);
    I41 = v1(yMat) .* I41 ./ W;
    
    res = I11 + I21 + I31 + I41 + u1Zero + u2Zero;
    
%     fx = @(s) (5 + s.^3) .* exp(3 * s) + mu * u1Zero;
%     fy = @(s) -2 * s .* sin(5 * s) + mu * u2Zero;
%     M = size(xMat, 1);
%     N = size(xMat, 2);
%     set_dim = 1000;
%     val = zeros(size(xMat));
%     for i = 1 : M
%         for j = 1 : N
%             x = xMat(i, j);
%             y = yMat(i, j);
%             
%             set_x_1 = linspace(0, x, set_dim);
%             I_x_1 = trapz(set_x_1, v_1(set_x_1) .* f_x(set_x_1));
%             I_x_1 = (v_2(x) / W) .* I_x_1;
%             
%             set_x_2 = linspace(x, 1, set_dim);
%             I_x_2 = trapz(set_x_2, v_2(set_x_2) .* f_x(set_x_2));
%             I_x_2 = (v_1(x) / W) .* I_x_2;
%             
%             set_y_1 = linspace(0, y, set_dim);
%             I_y_1 = trapz(set_y_1, v_1(set_y_1) .* f_y(set_y_1));
%             I_y_1 = (v_2(y) / W) .* I_y_1;
%             
%             set_y_2 = linspace(y, 1, set_dim);
%             I_y_2 = trapz(set_y_2, v_2(set_y_2) .* f_y(set_y_2));
%             I_y_2 = (v_1(y) / W) .* I_y_2;
%             
%             val(i, j) = I_x_1 + I_x_2 + I_y_1 + I_y_2 + u1Zero + u2Zero;
%         end
%     end
%     
%     res = val;
end


function res = v1fx(mu, u1zero, x) 
    res =  @(mu, u1zero, x) ((exp(x .* 3.0 + sqrt(mu) .* x) - 1.0) .* - 5.0) ./ ...
        (sqrt(mu) + 3.0) - ((exp(x .* 3.0 - sqrt(mu) .* x) - 1.0) .* 5.0) ./ ...
        (sqrt(mu) - 3.0) + 1.0 ./ (sqrt(mu) - 3.0).^4 .* 6.0 - 1.0 ./ ...
        (sqrt(mu) + 3.0).^4 .* 6.0 - sqrt(mu) .* u1zero .* (exp(sqrt(mu) .* x) ...
        -1.0) - sqrt(mu) .* u1zero .* (exp(-sqrt(mu) .* x) - 1.0) - exp(x .* 3.0 ...
        + sqrt(mu) .* x) .* 1.0 ./ (sqrt(mu) + 3.0).^4 .* (x .* (sqrt(mu) + 3.0) ...
        .* 6.0 - x .^2 .* (sqrt(mu) + 3.0) .^ 2 .* 3.0 + x .^3 .* (sqrt(mu) + 3.0) ...
        .^3 - 6.0) - exp(x .* 3.0 - sqrt(mu) .* x) .* 1.0 ./ (sqrt(mu) - 3.0) ...
        .^ 4 .* (x .* (sqrt(mu) - 3.0) .* 6.0 + x .^ 2 .* (sqrt(mu) - 3.0) .^2 ...
        .* 3.0 + x .^ 3 .* (sqrt(mu) - 3.0) .^ 3 + 6.0);
end


function res = v2fx(mu, u1zero, x)
    res =  @(mu, u1zero, x) (exp(x .* 3.0 - sqrt(mu) .* x + sqrt(mu) .* 2.0) ...
        .* 5.0) ./ (sqrt(mu) - 3.0) + ((exp(x .* 3.0 + sqrt(mu) .* x) - ...
        exp(sqrt(mu) + 3.0)) .* 5.0) ./ (sqrt(mu) + 3.0) - (exp(sqrt(mu) + 3.0) ...
        .* 5.0) ./ (sqrt(mu) - 3.0) - exp(sqrt(mu) + 3.0) .* 1.0 ./ (sqrt(mu) ...
        - 3.0) .^ 3 .* 1.5e+1 + exp(x .* 3.0 + sqrt(mu) .* x) .* 1.0 ./ ...
        (sqrt(mu) + 3.0) .^ 4 .* (x .* (sqrt(mu) + 3.0) .* 6.0 - x .^ 2 .* ...
        (sqrt(mu) + 3.0) .^ 2 .* 3.0 + x .^ 3 .* (sqrt(mu) + 3.0) .^ 3 - 6.0) ...
        + exp(x .* 3.0 - sqrt(mu) .* x + sqrt(mu) .* 2.0) .* 1.0 ./ (sqrt(mu) ...
        - 3.0) .^ 4 .* (mu .^ (3.0 ./ 2.0) .* x .^ 3 + mu .* x .^ 2 .* 3.0 ...
        - mu .* x .^3 .* 9.0 - x .^ 2 .* 2.7e+1 + x .^ 3 .* 5.4e+1 + 6.0) ...
        - exp(sqrt(mu) + 3.0) .* 1.0 ./ (sqrt(mu) - 3.0) .^ 4 .* (mu .* - ...
        6.0 + mu .^ (3.0 ./ 2.0) + 3.3e+1) - exp(sqrt(mu) + 3.0) .* 1.0 ./ ...
        (sqrt(mu) + 3.0) .^ 4 .* ((sqrt(mu) + 3.0) .^ 2 .* - 3.0 + ...
        (sqrt(mu) + 3.0) .^ 3 + sqrt(mu) .* 6.0 + 1.2e+1) + sqrt(mu) .* ...
        u1zero .* (exp(sqrt(mu) .* x) - exp(sqrt(mu))) - sqrt(mu) .* u1zero ...
        .* (exp(sqrt(mu)) - exp(-sqrt(mu) .* (x - 2.0))) + x .* exp(x .* 3.0 ...
        - sqrt(mu) .* x +sqrt(mu) .* 2.0) .* 1.0 ./ (sqrt(mu) - 3.0) .^ 3 .* ...
        (x .* - 6.0 + x .^ 2 .* 9.0 + 2.0) .* 3.0;
end


function res = v1fy(mu, u2zero, y)
    res = @(mu, u2zero, y) sqrt(mu) .* 1.0 ./ (mu + 2.5e+1) .^ 2 .* -4.0e+1 ...
        - sqrt(mu) .* u2zero .* sinh((sqrt(mu) .* y) ./ 2.0) .^ 2 .* 4.0 + ...
        exp(sqrt(mu) .* y) .* 1.0 ./ (mu + 2.5e+1) .^ 2 .* (sin(y .* 5.0) .* ...
        2.5e+1 - y .* cos(y .* 5.0) .* 1.25e+2 - mu .* sin(y .* 5.0) ...
        + sqrt(mu) .* cos(y .* 5.0) .* 1.0e+1 - mu .* y .* cos(y .* 5.0) ...
        .* 5.0 + sqrt(mu) .* y .* sin(y .* 5.0) .* 2.5e+1 + mu .^ (3.0 ./ 2.0) ...
        .* y .* sin(y .* 5.0)) .* 2.0 + exp(-sqrt(mu) .* y) .* 1.0 ./ (mu ...
        + 2.5e+1) .^ 2 .* (sin(y .* 5.0) .* -2.5e+1 + y .* cos(y .* 5.0) ...
        .* 1.25e+2 + mu .* sin(y .* 5.0) + sqrt(mu) .* cos(y .* 5.0) .* ...
        1.0e+1 + mu .* y .* cos(y .* 5.0) .* 5.0 + sqrt(mu) .* y .* ... 
        sin(y .* 5.0) .* 2.5e+1 + mu .^ (3.0 ./ 2.0) .* y .* sin(y .* 5.0)) .* 2.0;
end


function res = v2fy(mu, u2zero, y)
    res =  @(mu, u2zero, y) exp(-sqrt(mu) .* y + sqrt(mu) .* 2.0) .* 1.0 ...
        ./ (mu + 2.5e+1) .^ 2 .* (sin(y .* 5.0) .* -2.5e+1 + y .* cos(y .* 5.0) ...
        .* 1.25e+2 + mu .* sin(y .* 5.0) + sqrt(mu) .* cos(y .* 5.0) .* ...
        1.0e+1 + mu .* y .* cos(y .* 5.0) .* 5.0 + sqrt(mu) .* y .* ...
        sin(y .* 5.0) .* 2.5e+1 + mu .^ (3.0 ./ 2.0) .* y .* sin(y .* 5.0)) ...
        .* -2.0 - sqrt(mu) .* u2zero .* (exp(sqrt(mu)) - exp(-sqrt(mu) .* ...
        y + sqrt(mu) .* 2.0)) - exp(sqrt(mu) .* y) .* 1.0 ./ (mu + 2.5e+1) ...
        .^ 2 .* (sin(y .* 5.0) .* 2.5e+1 - y .* cos(y .* 5.0) .* 1.25e+2 ...
        - mu .* sin(y .* 5.0) + sqrt(mu) .* cos(y .* 5.0) .* 1.0e+1 - ...
        mu .* y .* cos(y .* 5.0) .* 5.0 + sqrt(mu) .* y .* sin(y .* 5.0) ...
        .* 2.5e+1 + mu .^ (3.0 ./ 2.0) .* y .* sin(y .* 5.0)) .* 2.0 + ...
        sqrt(mu) .* u2zero .* (exp(sqrt(mu) .* y) - exp(sqrt(mu))) + ...
        exp(sqrt(mu)) .* 1.0 ./ (mu + 2.5e+1) .^ 2 .* (cos(5.0) .* 1.25e+2 ...
        - sin(5.0) .* 2.5e+1 + mu .* cos(5.0) .* 5.0 + mu .* sin(5.0) ...
        + sqrt(mu) .* cos(5.0) .* 1.0e+1 + sqrt(mu) .* sin(5.0) .* 2.5e+1 ...
        + mu .^ (3.0 ./ 2.0) .* sin(5.0)) .* 2.0 + exp(sqrt(mu)) .* 1.0 ./ ...
        (mu + 2.5e+1) .^ 2 .* (cos(5.0) .* -1.25e+2 + sin(5.0) .* 2.5e+1 ...
        - mu .* cos(5.0) .* 5.0 - mu .* sin(5.0) + sqrt(mu) .* cos(5.0) .* ...
        1.0e+1 + sqrt(mu) .* sin(5.0) .* 2.5e+1 + mu .^ (3.0 ./ 2.0) .* sin(5.0)) .* 2.0;
end
   

% вызов solveDirichlet
function res = uNumerical(u1Zero, u2Zero, mu, M, N)
    fHandle = @(x,y) (5 + x.^2) .* exp(3 * x) - 2 * y .* sin(5 * y);
    %fHandle = @(x,y) 2*x.^2.*exp(2*x)+0*y;
    
    xiHandle = @(x) uAnalytical(x, zeros(size(x)), u1Zero, u2Zero, mu);
    etaHandle = @(y) uAnalytical(zeros(size(y)), y, u1Zero, u2Zero, mu);
    
    res = solveDirichlet(fHandle, xiHandle, etaHandle, mu, M, N);
end