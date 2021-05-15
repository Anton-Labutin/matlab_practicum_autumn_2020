% 
% Лабораторная работа 3
% Численные методы
% Вариант 2
%
% Задание 10
% Получить аппроксимацию преобразования Фурье F(λ) для каждой функции f(t) 
% из набора при помощи быстрого преобразования Фурье (БПФ), выбирая 
% различные шаги дискретизации исходной функции и различные окна, 
% ограничивающие область определения f(t). Построить графики F(λ). 
% Для первых двух функций f(t) вычислить F(λ) в явном виде и сравнить 
% графики F(λ), полученные из аналитического представления F (λ) и из 
% аппроксимации через БПФ.

figure_handle = gcf; % получаем handle фигуры для вывода графиков
step = 0.001;
inpLimVec = [-100, 100];
outLimVec = [-10, 10];

% @funci <-> @ftfunci, i = 1, 2
% @funci <-> [], i = 3, 4
plotFT(figure_handle, @func4, [], step, inpLimVec, outLimVec);


function [structure] = plotFT (hFigure, ... % handle фигуры, в которую осуществялется вывод графиков
    fHandle, ... % handle функции
    fFTHandle, ... % handle аналитическго преобразования Фурье для func1 и func2; [] - для func3 и func4
    step, ... % шаг дискретизации delta_t > 0
    inpLimVec, ... % [a, b] - окно для func_i (a < b)
    outLimVec) % [c, d] - окно для вывода преобразования Фурье, c < d
    
    % обработка свойств фигуры -------------------------------------------
    
    % UserData - свойство фигуры, хранящее метаинформацию
    SPlotInfo = get(hFigure, 'UserData'); % SPlotInfo - структура, хранящая метаинформацию
    
    if ~isempty(SPlotInfo)
        % очистка старых осей 
        delete(SPlotInfo.first_axises);
        delete(SPlotInfo.second_axises);
    else
        % инициализация метаинформации и осей
        SPlotInfo = struct('first_axises', [], 'second_axises', []);
    end
    
    % вычисление параметров ----------------------------------------------
    
    a = inpLimVec(1);
    b = inpLimVec(2);
    T = b - a; % длина окна вывода ПФ
    N = T / step; % число узов сеточной функции
    N = round(N);
    step = T / N; % перевычислим шаг дискретизации
    grid = a : step : b;
    grid_func_val = fHandle(grid);
    
    % продлим fHandle по периоду T ---------------------------------------
    
    T_cnt = 0;
    if a > 0 && b > 0
        while T_cnt * T < a
            T_cnt = T_cnt + 1;
        end
    elseif a < 0 && b < 0
        while T_cnt * T > b 
            T_cnt = T_cnt - 1;
        end
    end
    new_limit = T_cnt * T;
    
    grid_sep_idx = 1;
    size(grid);
    while grid_sep_idx < size(grid, 2) && grid(grid_sep_idx + 1) <= new_limit
        grid_sep_idx = grid_sep_idx + 1;
    end
    
    % сдвиг функции с [a, b] на [0, T]
    func_on_grid_0T = zeros(1, N + 1); % N - число узлов сеточной функции
    func_on_grid_0T(1 : (N + 1) - grid_sep_idx) = grid_func_val(grid_sep_idx + 1 : N + 1);
    func_on_grid_0T((N + 1) - grid_sep_idx + 1 : N + 1) = grid_func_val(1 : grid_sep_idx);
    
    % бысрое преобразование Фурье ----------------------------------------
    fourier_on_grid_0T = step * fft(func_on_grid_0T);

    delta_lyambda = (2 * pi) / T;
    
    new_T = delta_lyambda * N;
    left = -new_T;
    right = new_T;
    
    c = outLimVec(1);
    d = outLimVec(2);
    
    newT_in_cd_cnt = 2;
    while c < left || d > right
        left = left - new_T;
        right = right + new_T;
        newT_in_cd_cnt = newT_in_cd_cnt + 2;
    end
    
    new_grid = left : delta_lyambda : right;
    fourier_on_new_grid = repmat(fourier_on_grid_0T(2 : end), 1, newT_in_cd_cnt);
    
    % строим графики -----------------------------------------------------
    % вещественная часть преобразования Фурье
    subplot(2, 1, 1);
    hold on;
    plot(new_grid(2 : end), real(fourier_on_new_grid), 'b');
    
    if ~isempty(fFTHandle)
        grid_cd = outLimVec(1) : step : outLimVec(2);
        plot(grid_cd, real(fFTHandle(grid_cd)) , 'r');
        legend({'Вещественная часть аппроксимации F(\lambda) через БПФ', 'Вещественная часть аналитически посчитанного F(\lambda)'});  
    else
        legend('Вещественная часть аппроксимации F(\lambda) через БПФ');    
    end
    
    axis([outLimVec(1) outLimVec(2) min(real(fourier_on_new_grid)) max(real(fourier_on_new_grid))]);
    xlabel('\lambda');
    ylabel('Re F(\lambda)');
    SPlotInfo.first_axises = gca; % сохраняем текущие оси
    hold off;
    
    % мнимая часть преобразования Фурье
    subplot(2, 1, 2);
    hold on;
    plot(new_grid(2 : end), imag(fourier_on_new_grid), 'b');
    
    if ~isempty(fFTHandle)
        grid_cd = outLimVec(1) : step : outLimVec(2); 
        plot(grid_cd, imag(fFTHandle(grid_cd)) , 'r');
        
        % plot(grid_cd + (2 * pi / step), imag(fFTHandle(grid_cd)), 'r');
        % plot(grid_cd - (2 * pi / step), imag(fFTHandle(grid_cd)), 'r');
        % plot(grid_cd + (4 * pi / step), imag(fFTHandle(grid_cd)), 'r');
        % plot(grid_cd - (4 * pi / step), imag(fFTHandle(grid_cd)), 'r');
        
        legend({'Мнимая часть аппроксимации F(\lambda) через БПФ', 'Мнимая часть аналитически посчитанного F(\lambda)'});
    else 
        legend('Мнимая часть аппроксимации F(\lambda) через БПФ');
    end
    
    axis([outLimVec(1) outLimVec(2) min(imag(fourier_on_new_grid)) max(imag(fourier_on_new_grid))]);
    xlabel('\lambda');
    ylabel('Im F(\lambda)');
    SPlotInfo.second_axises = gca; % сохраняем текущие оси
    hold off;
    
    % сохранить структуру structure с метаинформацией в свойство фигуры
    % UserData
    structure = struct('nPoints', N, 'step', step, 'inpLimVec', inpLimVec, 'outLimVec', outLimVec);
    set(hFigure, 'UserData', SPlotInfo);
   
    %figure(2);
    %hold on;
    %true_fourier_on_big_gr = repmat(fFTHandle(func_on_grid_0T(2:end)), 1, newT_in_cd_cnt);
    %plot(new_grid(2:end), imag(true_fourier_on_big_gr)); 
    %hold off;
end


function func_value = func1(t)
    func_value = t .* exp(-2 * t.^2);
end


function func_value = func2(t)
    func_value = atan(3 * t) - atan(2 * t);
end


function func_value = func3(t)
   for i = 1 : length(t)
        if t(i) ~= 0
            func_value(i) = (-cos((t(i))^3) + 1) ./ (t(i))^5;
        else
            func_value(i) = 0;
        end
    end
end


function func_value = func4(t)
    func_value = exp(-2 * abs(t)) ./ (1 + (sin(t)).^2);
end


function func_value = ftfunc1(lyambda)
    Re = zeros(1, length(lyambda));
    Im = (-sqrt(2 * pi) / 8) * lyambda .* exp((-lyambda.^2) / 8);
    func_value = Re + 1i * Im;
end


function func_value = ftfunc2(lyambda)
    Re = zeros(1, length(lyambda));
    Im = (pi ./ lyambda) .* (exp(-abs(lyambda) / 2) - exp(-abs(lyambda) / 3));
    func_value = Re + 1i * Im;
end