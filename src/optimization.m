% =========================================================================
% MOM_Richmond_Antenna_FINAL.m
% Реализация метода моментов Ричмонда для тонкопроволочной антенны
% MATLAB R2019b, без использования Integral, sinint, cosint и других тулбоксов
% =========================================================================

clc; clear; close all;

%% Параметры расчёта
freq = 300e6;               % 300 МГц
c = 3e8;
lambda = c / freq;          % 1 м
k = 2*pi/lambda;
Z0 = 120*pi;                % волновое сопротивление свободного пространства
radius = lambda / 200;      % радиус проводника (тонкий)
max_seg_len = lambda / 15;  % максимальная длина сегмента

% Точка питания (центр нижнего вибратора)
feed_point = [0, 0, 0];

%% Функция построения узлов антенны (8 узлов)
% Вход: L - длина вибратора, d - расстояние между вибраторами, off - смещение короткозамыкателей
build_nodes = @(L, d, off) [
    -L/2, 0, 0;      % 1
     L/2, 0, 0;      % 2
    -L/2, d, 0;      % 3
     L/2, d, 0;      % 4
    -off, 0, 0;      % 5
    -off, d, 0;      % 6
     off, 0, 0;      % 7
     off, d, 0       % 8
];

% Топология соединений [индекс1, индекс2, радиус]
segments_ref = [
    1,2, radius;   % вибратор 1 (нижний)
    3,4, radius;   % вибратор 2 (верхний)
    5,6, radius;   % левый короткозамыкатель
    7,8, radius;   % правый короткозамыкатель
    1,5, radius;   % левая часть нижнего
    5,2, radius;   % правая часть нижнего
    3,6, radius;   % левая часть верхнего
    6,4, radius    % правая часть верхнего
];

%% ===================== 1. ОПТИМИЗАЦИЯ ГЕОМЕТРИИ (ПЕРЕБОР) =====================
L_vals = linspace(0.3, 0.7, 4) * lambda;
spacing_vals = linspace(0.03, 0.3, 4) * lambda;
off_vals = linspace(0.1, 0.4, 3) * lambda;

fprintf('Оптимизация (перебор %d вариантов)...\n', length(L_vals)*length(spacing_vals)*length(off_vals));
best_gain = -inf;
best_params = [0.5*lambda, 0.1*lambda, 0.25*lambda];
cnt = 0;

for L = L_vals
    for d = spacing_vals
        for off = off_vals
            cnt = cnt + 1;
            nodes = build_nodes(L, d, off);
            [~, ~, basis] = discretize_and_basis(nodes, segments_ref, max_seg_len);
            if length(basis) < 5; continue; end
            N = length(basis);
            Zmat = complex(zeros(N,N));
            for m = 1:N
                for n = 1:N
                    Zmat(m,n) = mutual_impedance_richmond(basis{m}, basis{n}, k, Z0);
                end
            end
            % Поиск сегмента, ближайшего к точке питания
            feed_idx = 1;
            min_dist = inf;
            for m = 1:N
                d_seg = norm(basis{m}.center - feed_point);
                if d_seg < min_dist
                    min_dist = d_seg;
                    feed_idx = m;
                end
            end
            V = zeros(N,1);
            V(feed_idx) = 1;
            I = Zmat \ V;
            gain_dB = estimate_gain_from_currents(basis, I, k, Z0);
            if gain_dB > best_gain
                best_gain = gain_dB;
                best_params = [L, d, off];
                fprintf('  Вариант %d: L=%.3f d=%.3f off=%.3f -> Gain=%.2f дБ (лучший)\n', cnt, L, d, off, gain_dB);
            elseif mod(cnt,5)==0
                fprintf('  Вариант %d: Gain=%.2f дБ\n', cnt, gain_dB);
            end
        end
    end
end

L_opt = best_params(1);
spacing_opt = best_params(2);
short_opt = best_params(3);
fprintf('\n=== Оптимальные параметры ===\n');
fprintf('L = %.3f м (%.2f ?)\n', L_opt, L_opt/lambda);
fprintf('spacing = %.3f м (%.2f ?)\n', spacing_opt, spacing_opt/lambda);
fprintf('short_offset = %.3f м (%.2f ?)\n', short_opt, short_opt/lambda);
fprintf('Gain = %.2f дБ\n', best_gain);

%% ===================== 2. ФИНАЛЬНЫЙ РАСЧЁТ ДЛЯ ОПТИМАЛЬНОЙ АНТЕННЫ =====================
nodes_opt = build_nodes(L_opt, spacing_opt, short_opt);
[wire_nodes, wire_segments, basis_opt] = discretize_and_basis(nodes_opt, segments_ref, max_seg_len);
N = length(basis_opt);
fprintf('\nКоличество сегментов: %d\n', N);

Zmat = complex(zeros(N,N));
fprintf('Расчёт матрицы импедансов (MoM Ричмонда)...\n');
for m = 1:N
    if mod(m,10)==0; fprintf('  строка %d из %d\n', m, N); end
    for n = 1:N
        Zmat(m,n) = mutual_impedance_richmond(basis_opt{m}, basis_opt{n}, k, Z0);
    end
end

% Находим сегмент питания
feed_idx = 1;
min_dist = inf;
for m = 1:N
    d_seg = norm(basis_opt{m}.center - feed_point);
    if d_seg < min_dist
        min_dist = d_seg;
        feed_idx = m;
    end
end
V = zeros(N,1);
V(feed_idx) = 1;
I = Zmat \ V;
Zin = 1 / I(feed_idx);
fprintf('Входное сопротивление: %.1f + j%.1f Ом\n', real(Zin), imag(Zin));

% Диаграмма направленности и усиление
[Gain_dB, Theta_deg, Phi_deg] = far_field_richmond(basis_opt, I, k, Z0);
fprintf('Максимальный коэффициент усиления: %.2f дБ\n', max(Gain_dB(:)));

%% ===================== 3. ВИЗУАЛИЗАЦИЯ =====================
% 3.1 Геометрия антенны
figure('Name', 'Геометрия антенны', 'NumberTitle', 'off');
hold on; axis equal; grid on;
xlabel('X, м'); ylabel('Y, м'); zlabel('Z, м');
title('Проволочная модель (красный – сегмент питания)');
for i = 1:size(wire_segments,1)
    i1 = wire_segments(i,1); i2 = wire_segments(i,2);
    P1 = wire_nodes(i1,:); P2 = wire_nodes(i2,:);
    if i == feed_idx
        plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], 'r-', 'LineWidth', 3);
    else
        plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], 'b-', 'LineWidth', 1.5);
    end
end
plot3(feed_point(1), feed_point(2), feed_point(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
view(3); hold off;

% 3.2 3D диаграмма направленности (сферическая проекция)
G_lin = 10.^(Gain_dB/10);
[Theta, Phi] = meshgrid(Theta_deg, Phi_deg);
theta_rad = Theta * pi/180;
phi_rad = Phi * pi/180;
X = G_lin .* sin(theta_rad) .* cos(phi_rad);
Y = G_lin .* sin(theta_rad) .* sin(phi_rad);
Z = G_lin .* cos(theta_rad);
figure('Name', '3D диаграмма направленности', 'NumberTitle', 'off');
surf(X, Y, Z, Gain_dB, 'EdgeColor', 'none');
colormap(jet); colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Коэффициент усиления, дБ');
axis equal; grid on; view(3);

% 3.3 Сечения E-plane (phi = 0° и 180°)
idx0 = find(Phi_deg == 0, 1);
idx180 = find(Phi_deg == 180, 1);
if ~isempty(idx0) && ~isempty(idx180)
    figure('Name', 'E-plane', 'NumberTitle', 'off');
    pat0 = interp1(Theta_deg, Gain_dB(idx0,:), Theta_deg);
    pat180 = interp1(Theta_deg, Gain_dB(idx180,:), Theta_deg);
    polarplot(Theta_deg*pi/180, 10.^(pat0/10), 'b-', 'LineWidth', 2);
    hold on;
    polarplot(Theta_deg*pi/180, 10.^(pat180/10), 'r--', 'LineWidth', 2);
    legend('phi=0°','phi=180°', 'Location','best');
    title('E-plane (нормированная ДН)');
end

% 3.4 Сечение H-plane (phi = 90°)
idx90 = find(Phi_deg == 90, 1);
if ~isempty(idx90)
    figure('Name', 'H-plane', 'NumberTitle', 'off');
    pat90 = 10.^(Gain_dB(idx90,:)/10);
    polarplot(Theta_deg*pi/180, pat90, 'g-', 'LineWidth', 2);
    title('H-plane (нормированная ДН)');
end

%% ==================== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ====================

% -------------------------------------------------------------------------
% Дискретизация проводников и создание базисных функций (с id)
function [wire_nodes, wire_segments, basis] = discretize_and_basis(nodes, segments, max_len)
    wire_nodes = [];
    wire_segments = [];
    for i = 1:size(segments,1)
        P1 = nodes(segments(i,1),:);
        P2 = nodes(segments(i,2),:);
        L_seg = norm(P2-P1);
        n_div = max(1, ceil(L_seg / max_len));
        t = linspace(0, 1, n_div+1);
        pts = (1-t')*P1 + t'*P2;
        start_idx = size(wire_nodes,1) + 1;
        wire_nodes = [wire_nodes; pts];
        for j = 1:n_div
            wire_segments = [wire_segments; start_idx+j-1, start_idx+j, segments(i,3)];
        end
    end
    M = size(wire_segments,1);
    basis = cell(M,1);
    for m = 1:M
        i1 = wire_segments(m,1);
        i2 = wire_segments(m,2);
        P1 = wire_nodes(i1,:);
        P2 = wire_nodes(i2,:);
        L = norm(P2-P1);
        basis{m}.id = m;
        basis{m}.P1 = P1;
        basis{m}.P2 = P2;
        basis{m}.L = L;
        basis{m}.center = (P1+P2)/2;
        basis{m}.dir = (P2-P1)/L;
        basis{m}.radius = wire_segments(m,3);
    end
end

% -------------------------------------------------------------------------
% Вычисление взаимного импеданса между двумя диполями (метод Ричмонда)
function Zmn = mutual_impedance_richmond(basis_m, basis_n, k, Z0)
    if basis_m.id == basis_n.id
        % Собственный импеданс: численное интегрирование с регуляризацией
        % (избегаем сингулярности, заменяя R на sqrt(R^2 + a^2))
        a = basis_m.radius;
        n_pts = 8;
        [xq, wq] = lgwt(n_pts, -0.5, 0.5);
        integral = 0;
        for i = 1:n_pts
            s_m = xq(i);
            t_m = s_m * basis_m.L;
            r_m = basis_m.center + t_m * basis_m.dir;
            J_m = sin(k*(basis_m.L/2 - abs(t_m)));
            if J_m < 0; J_m = 0; end   % защита от численных погрешностей
            for j = 1:n_pts
                s_n = xq(j);
                t_n = s_n * basis_m.L;
                r_n = basis_m.center + t_n * basis_m.dir;
                J_n = sin(k*(basis_m.L/2 - abs(t_n)));
                if J_n < 0; J_n = 0; end
                R = norm(r_m - r_n);
                R_reg = sqrt(R^2 + a^2);   % регуляризация
                G = exp(-1j*k*R_reg) / (4*pi*R_reg);
                integral = integral + wq(i)*wq(j) * J_m * J_n * G;
            end
        end
        integral = integral * (basis_m.L * basis_m.L);
        Zmn = -1j * k * Z0 * integral;
        return;
    end
    
    % Разные сегменты: обычное интегрирование
    n_pts = 6;
    [xq, wq] = lgwt(n_pts, -0.5, 0.5);
    integral = 0;
    for i = 1:n_pts
        s_m = xq(i);
        t_m = s_m * basis_m.L;
        r_m = basis_m.center + t_m * basis_m.dir;
        J_m = sin(k*(basis_m.L/2 - abs(t_m)));
        for j = 1:n_pts
            s_n = xq(j);
            t_n = s_n * basis_n.L;
            r_n = basis_n.center + t_n * basis_n.dir;
            J_n = sin(k*(basis_n.L/2 - abs(t_n)));
            R = norm(r_m - r_n);
            if R < 1e-9
                G = 1/(4*pi*basis_m.radius);
            else
                G = exp(-1j*k*R) / (4*pi*R);
            end
            integral = integral + wq(i)*wq(j) * J_m * J_n * G;
        end
    end
    integral = integral * (basis_m.L * basis_n.L);
    Zmn = -1j * k * Z0 * integral;
end

% -------------------------------------------------------------------------
% Квадратура Гаусса-Лежандра (узлы и веса на интервале [a,b])
function [x, w] = lgwt(N, a, b)
    N = N - 1;
    N1 = N+1; N2 = N+2;
    xu = linspace(-1,1,N1)';
    y = cos((2*(0:N)'+1)*pi/(2*N+2));
    iter = 1;
    while max(abs(y-xu)) > eps && iter < 20
        p_prev = 1; p_curr = y;
        for k = 1:N
            p_next = ((2*k+1)*y.*p_curr - k*p_prev)/(k+1);
            p_prev = p_curr; p_curr = p_next;
        end
        pp = (N+1)*(p_prev - y.*p_curr)./(1-y.^2);
        y = y - p_curr./pp;
        iter = iter + 1;
    end
    x = y;
    w = 2./((1-y.^2).*pp.^2);
    x = 0.5*(b-a)*x + 0.5*(b+a);
    w = 0.5*(b-a)*w;
end

% -------------------------------------------------------------------------
% Оценка усиления по токам (используется во время оптимизации)
function gain_dB = estimate_gain_from_currents(basis, I, k, Z0)
    Theta_deg = 0:10:180;
    Phi_deg = 0:20:360;
    [Theta, Phi] = meshgrid(Theta_deg, Phi_deg);
    theta_rad = Theta * pi/180;
    phi_rad = Phi * pi/180;
    ux = sin(theta_rad) .* cos(phi_rad);
    uy = sin(theta_rad) .* sin(phi_rad);
    uz = cos(theta_rad);
    E_max = 0;
    for m = 1:length(basis)
        center = basis{m}.center;
        dir = basis{m}.dir;
        Lm = basis{m}.L;
        Im = I(m);
        factor = sin(0.5*k*Lm) / (0.5*k*Lm);
        I_eff = Im * Lm * factor;
        phase = exp(-1j*k*(center(1)*ux + center(2)*uy + center(3)*uz));
        costheta_dir = dir(1)*ux + dir(2)*uy + dir(3)*uz;
        sin_theta_dir = sqrt(max(0, 1 - costheta_dir.^2));
        E_m = abs(1j * Z0 * k * I_eff .* sin_theta_dir .* phase / (4*pi));
        E_max = max(E_max, max(E_m(:)));
    end
    if E_max > 0 && sum(abs(I).^2) > 0
        gain_lin = (4*pi) * (E_max^2) / (sum(abs(I).^2) * Z0);
        gain_dB = 10*log10(gain_lin);
    else
        gain_dB = -30;
    end
end

% -------------------------------------------------------------------------
% Полный расчёт диаграммы направленности и коэффициента усиления
function [Gain_dB, Theta_deg, Phi_deg] = far_field_richmond(basis, I, k, Z0)
    Theta_deg = 0:10:360;
    Phi_deg = 0:10:360;
    [Theta, Phi] = meshgrid(Theta_deg, Phi_deg);
    theta_rad = Theta * pi/180;
    phi_rad = Phi * pi/180;
    ux = sin(theta_rad) .* cos(phi_rad);
    uy = sin(theta_rad) .* sin(phi_rad);
    uz = cos(theta_rad);
    E_total = zeros(size(theta_rad));
    for m = 1:length(basis)
        center = basis{m}.center;
        dir = basis{m}.dir;
        Lm = basis{m}.L;
        Im = I(m);
        factor = sin(0.5*k*Lm) / (0.5*k*Lm);
        I_eff = Im * Lm * factor;
        phase = exp(-1j*k*(center(1)*ux + center(2)*uy + center(3)*uz));
        costheta_dir = dir(1)*ux + dir(2)*uy + dir(3)*uz;
        sin_theta_dir = sqrt(max(0, 1 - costheta_dir.^2));
        E_m = 1j * Z0 * k * I_eff .* sin_theta_dir .* phase / (4*pi);
        E_total = E_total + E_m;
    end
    E_abs = abs(E_total);
    E_max = max(E_abs(:));
    if E_max > 0
        Gain_lin = (E_abs / E_max).^2 * 4*pi;
    else
        Gain_lin = zeros(size(E_abs));
    end
    Gain_dB = 10*log10(Gain_lin + eps);
end


