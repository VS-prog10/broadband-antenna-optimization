% ========================================================================
% ПРОЕКТ ДВА: МЕТОД РИЧМОНДА ДЛЯ ДВУХВИБРАТОРНОЙ АНТЕННЫ (ПО РИС. 2.1)
% ========================================================================
clear; clc; close all;

%% 1. ПАРАМЕТРЫ АНТЕННЫ (ИЗ ДОКУМЕНТА, РИС. 2.1)
f0 = 300e6;                     % частота, Гц
c = 299792458;                  % скорость света
lambda = c / f0;                % длина волны, м

L = 0.5 * lambda;               % длина каждого вибратора
d = 0.1 * lambda;               % расстояние между вибраторами
short_pos = 0.25 * L;           % положение короткозамыкателей (от центра)
a_rad = 0.005;                  % радиус провода, м

fprintf('=== Параметры антенны (проект два) ===\n');
fprintf('f0 = %.1f МГц, ? = %.3f м\n', f0/1e6, lambda);
fprintf('L = %.3f м, d = %.3f м, short_pos = %.3f м\n', L, d, short_pos);
fprintf('Радиус провода = %.4f м\n', a_rad);

%% 2. ПОСТРОЕНИЕ ГЕОМЕТРИИ (СПИСОК ОТРЕЗКОВ)
segments = [];
% Левый вибратор (сплошной, питание в центре)
segments = [segments; 0, 0, -L/2,   0, 0,  L/2];
% Правый вибратор
segments = [segments; d, 0, -L/2,   d, 0,  L/2];
% Верхний короткозамыкатель
segments = [segments; 0, 0,  short_pos,   d, 0,  short_pos];
% Нижний короткозамыкатель
segments = [segments; 0, 0, -short_pos,   d, 0, -short_pos];
fprintf('Построено %d отрезков\n', size(segments,1));

%% 3. ДИСКРЕТИЗАЦИЯ НА СЕГМЕНТЫ (кусочно-синусоидальные токи)
seg_len_max = lambda / 15;
dipoles = [];   % [x_c, y_c, z_c, half_len, ux, uy, uz]
for s = 1:size(segments,1)
    P1 = segments(s,1:3);
    P2 = segments(s,4:6);
    len = norm(P2-P1);
    Nsub = max(1, round(len / seg_len_max));
    for i = 1:Nsub
        t1 = (i-1)/Nsub;
        t2 = i/Nsub;
        start = P1 + t1*(P2-P1);
        endp  = P1 + t2*(P2-P1);
        center = (start + endp)/2;
        half_len = norm(endp-start)/2;
        dir = (endp-start) / (2*half_len);
        dipoles = [dipoles; center, half_len, dir];
    end
end
N = size(dipoles,1);
fprintf('Количество сегментов (диполей): %d\n', N);

%% 4. ОПРЕДЕЛЕНИЕ ТОЧКИ ПИТАНИЯ (ближайший сегмент к (0,0,0) на левом вибраторе)
min_dist = inf;
feed_idx = 1;
for i = 1:N
    if abs(dipoles(i,1)) < 1e-6   % x ? 0
        dist = norm(dipoles(i,1:3));
        if dist < min_dist
            min_dist = dist;
            feed_idx = i;
        end
    end
end
fprintf('Точка питания: сегмент %d, центр (%.4f, %.4f, %.4f)\n', ...
        feed_idx, dipoles(feed_idx,1), dipoles(feed_idx,2), dipoles(feed_idx,3));

%% 5. ВЫЧИСЛЕНИЕ МАТРИЦЫ ИМПЕДАНСОВ
k = 2*pi / lambda;
omega = 2*pi*f0;
mu = 4*pi*1e-7;
eta0 = 120*pi;

Z = zeros(N, N);
fprintf('Формирование матрицы импедансов (может занять время)...\n');
for m = 1:N
    for n = 1:N
        if m == n
            Z(m,n) = selfImpedance(dipoles(m,4), a_rad, k, omega, mu);
        else
            Z(m,n) = mutualImpedance(dipoles(m,:), dipoles(n,:), a_rad, k, omega, mu);
        end
    end
end

%% 6. РЕШЕНИЕ СИСТЕМЫ
V = zeros(N,1);
V(feed_idx) = 1.0;
I = Z \ V;

%% 7. ВХОДНОЕ СОПРОТИВЛЕНИЕ
Zin = V(feed_idx) / I(feed_idx);
fprintf('\n=== РЕЗУЛЬТАТЫ МЕТОДА РИЧМОНДА ===\n');
fprintf('Входное сопротивление: Zin = %.1f + j%.1f Ом\n', real(Zin), imag(Zin));

%% 8. ДИАГРАММА НАПРАВЛЕННОСТИ И КНД
theta = linspace(0, pi, 181);
E_abs = zeros(size(theta));
for n = 1:N
    x = dipoles(n,1);   z = dipoles(n,3);
    half_len = dipoles(n,4);
    uz = dipoles(n,7);
    In = I(n);
    if abs(uz) > 0.9   % учитываем только вертикальные токи
        for ti = 1:length(theta)
            th = theta(ti);
            if abs(sin(th)) > 1e-6
                F = (cos(k*half_len*cos(th)) - cos(k*half_len)) / sin(th);
            else
                F = 0;
            end
            phase = exp(-1j * k * (x*sin(th) + z*cos(th)));
            E_theta = 1j * eta0 * k * In * half_len * F * phase / (4*pi);
            E_abs(ti) = E_abs(ti) + abs(E_theta);
        end
    end
end
E_norm = E_abs / max(E_abs);
U = E_abs.^2 / (2*eta0);
P_rad = 2*pi * trapz(theta, U .* sin(theta));
D = 4*pi * max(U) / P_rad;
D_dBi = 10*log10(D);
fprintf('КНД: %.2f дБи\n', D_dBi);

%% 9. ВИЗУАЛИЗАЦИЯ
% ДН
figure('Name', 'Диаграмма направленности');
polarplot(theta, E_norm, 'LineWidth', 2);
title(sprintf('ДН (метод Ричмонда), КНД = %.1f дБи', D_dBi));
grid on;
% Геометрия
figure('Name', 'Геометрия антенны');
hold on; axis equal; grid on; view(3);
xlabel('X, м'); ylabel('Y, м'); zlabel('Z, м');
title('Двухвибраторная антенна с короткозамыкателями');
for s = 1:size(segments,1)
    X = [segments(s,1), segments(s,4)];
    Y = [segments(s,2), segments(s,5)];
    Z = [segments(s,3), segments(s,6)];
    plot3(X, Y, Z, 'b-', 'LineWidth', 2);
end
plot3(dipoles(feed_idx,1), dipoles(feed_idx,2), dipoles(feed_idx,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
legend('Проводники', 'Питание');
hold off;

%% ========== ЛОКАЛЬНЫЕ ФУНКЦИИ (каждая закрыта end) ==========

function Zii = selfImpedance(half_len, a, k, omega, mu)
    Nint = 32;
    [z_int, w_int] = lgwt(Nint, -half_len, half_len);
    sum_val = 0;
    for i = 1:Nint
        z1 = z_int(i);
        I1 = sin(k*(half_len - abs(z1))) / sin(k*half_len);
        for j = 1:Nint
            z2 = z_int(j);
            I2 = sin(k*(half_len - abs(z2))) / sin(k*half_len);
            R = sqrt(a^2 + (z1 - z2)^2);
            if R < 1e-12, R = a; end
            G = exp(-1j*k*R) / (4*pi*R);
            sum_val = sum_val + w_int(i)*w_int(j) * I1 * I2 * G;
        end
    end
    Zii = 1j * omega * mu * (2*half_len)^2 * sum_val;
end

function Zmn = mutualImpedance(dip_m, dip_n, a, k, omega, mu)
    dm = dip_m(4); dn = dip_n(4);
    cm = dip_m(1:3); cn = dip_n(1:3);
    um = dip_m(5:7); un = dip_n(5:7);
    [xq, wq] = lgwt(8, -1, 1);
    sum_val = 0;
    for i = 1:length(xq)
        t_m = xq(i) * dm;
        pm = cm + um * t_m;
        Im = sin(k*(dm - abs(t_m))) / sin(k*dm);
        for j = 1:length(xq)
            t_n = xq(j) * dn;
            pn = cn + un * t_n;
            In = sin(k*(dn - abs(t_n))) / sin(k*dn);
            R = norm(pm - pn);
            if R < a, R = a; end
            G = exp(-1j*k*R) / (4*pi*R);
            sum_val = sum_val + wq(i)*wq(j) * Im * In * G;
        end
    end
    Zmn = 1j * omega * mu * (dm*dn) * sum_val;
end

function [x, w] = lgwt(N, a, b)
    N = N-1; N1 = N+1; N2 = N+2;
    xu = linspace(-1,1,N1)';
    y = cos((2*(0:N)'+1)*pi/(2*N+2)) + (0.27/N1)*sin(pi*xu*N/N2);
    for iter = 1:10
        L = zeros(N1, N2);
        Lp = zeros(N1, N2);
        L(:,1) = 1; L(:,2)=y;
        Lp(:,1)=0; Lp(:,2)=1;
        for k=2:N1
            L(:,k+1) = ((2*k-1)*y.*L(:,k) - (k-1)*L(:,k-1))/k;
            Lp(:,k+1) = ((2*k-1)*(L(:,k) + y.*Lp(:,k)) - (k-1)*Lp(:,k-1))/k;
        end
        y_new = y - L(:,N2)./Lp(:,N2);
        if max(abs(y_new-y)) < 1e-12, break; end
        y = y_new;
    end
    x = (a*(1-y) + b*(1+y))/2;
    w = (b-a) ./ ((1-y.^2).*Lp(:,N2).^2) * (N2/2)^2;
end
