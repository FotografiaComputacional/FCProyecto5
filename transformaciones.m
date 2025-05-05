
% ------------------------------------------------------------
% Funciones
function T = get_afin(xy, uv)
    
    x = xy(1, :)'; 
    y = xy(2, :)';
    u = uv(1, :)';
    v = uv(2, :)';
    
    H_u = [x, y, ones(3, 1)];
    H_v = [x, y, ones(3, 1)];
    
    coefs_u = H_u \ u;
    coefs_v = H_v \ v;
    
    T = [coefs_u'; coefs_v'; 0 0 1];
end

function [uv] = aplica_T(xy, T)
    % Paso 1
    xy_homogeneas = [xy; ones(1, size(xy, 2))];
    % Paso 2
    resultado = T * xy_homogeneas;
    % Paso 3
    uv = resultado(1:2, :) ./ repmat(resultado(3, :), 2, 1);
end

function T = get_proy(xy, uv)
    x = xy(1, :)';
    y = xy(2, :)';
    u = uv(1, :)';
    v = uv(2, :)';
    
    A = zeros(8, 8);
    b = zeros(8, 1);
    
    for i = 1:4
        % Ecuaciones para u
        A(2*i-1, 1) = x(i);
        A(2*i-1, 2) = y(i);
        A(2*i-1, 3) = 1;
        A(2*i-1, 7) = -u(i)*x(i);
        A(2*i-1, 8) = -u(i)*y(i);
        b(2*i-1) = u(i);
        
        % Ecuaciones para v
        A(2*i, 4) = x(i);
        A(2*i, 5) = y(i);
        A(2*i, 6) = 1;
        A(2*i, 7) = -v(i)*x(i);
        A(2*i, 8) = -v(i)*y(i);
        b(2*i) = v(i);
    end
    
    params = A \ b;
    
    % Parametros
    a = params(1);
    b = params(2);
    c = params(3);
    d = params(4);
    e = params(5);
    f = params(6);
    g = params(7);
    h = params(8);
    
    T = [a b c; d e f; g h 1];
end

function [im2, RU, RV] = warp_img(im, T)
    % 1) Convertir la imagen a double (valores entre 0 y 1) y hallar su alto N y ancho M
    im = double(im) / 255;
    [N, M, ~] = size(im);
    
    % 2) Construir los vectores RU y RV con el rango de coordenadas de la imagen en el espacio destino
    x = [1, M, M, 1];
    y = [1, 1, N, N];
    xy = [x; y];
    
    uv = aplica_T(xy, T);
    
    u_min = floor(min(uv(1,:)));
    u_max = ceil(max(uv(1,:)));
    v_min = floor(min(uv(2,:)));
    v_max = ceil(max(uv(2,:)));
    
    RU = u_min:u_max;
    RV = v_min:v_max;
    
    % 3) Reservar memoria para la imagen de salida im2
    N_new = length(RV);
    M_new = length(RU);
    im2 = zeros(N_new, M_new, 3);

    X = zeros(N_new, M_new);
    Y = zeros(N_new, M_new);
    
    % 4) Calcular T^-1 la matriz inversa de T (para el warping "inverso")
    T_inv = inv(T);
    
    % 5) Hacer un bucle para barrer todas las filas de la imagen destino
    for k = 1:N_new
        % a) Crear matriz uv con las coordenadas u y v de esa fila
        u = RU;
        v = RV(k) * ones(size(RU));
        uv = [u; v];
        % b) Obtener las correspondientes coordenadas xy en la imagen de partida
        xy = aplica_T(uv, T_inv);
        % c) Guardar las coordenadas
        X(k,:) = xy(1,:);
        Y(k,:) = xy(2,:);
    end
    
    % 6) Usar interp2 para interpolar los valores de la imagen en las coordenadas (X,Y)
    for c = 1:3
        im2(:,:,c) = interp2(1:M, 1:N, im(:,:,c), X, Y, 'linear');
    end
    end


% ------------------------------------------------------------
% Parte 1
disp('------------------------ Parte 1 ------------------------')
% --- Obtención de la matriz T en una transformación afín ---
disp('Obtención de la matriz T en una transformación afín')


xy = [200 400 600; 375 125 375];
uv = [87 319 600; 445 69 375];

T = get_afin(xy, uv);

% Mostrar la matriz T resultante
disp('La matriz de transformacion afín T es:');
disp(T);

% Comprobacion
uv2 = aplica_T(xy, T);
dif_uv = uv - uv2;
disp('Diferencias entre coordenadas originales y transformadas:');
disp(dif_uv);

desv_std = std2(dif_uv(:));
fprintf('Desviacion estándar de las diferencias: %e\n', desv_std);

TT = get_afin(uv, xy);
disp('La matriz de transformacion inversa TT es:');
disp(TT);

T_inv = inv(T);
disp('La inversa de la matriz T es:');
disp(T_inv);


% Comprobar que TT = inv(T)
dif_T = abs(TT - T_inv);
max_dif_T = max(dif_T(:));
fprintf('Diferencia máxima entre TT e inv(T): %e\n', max_dif_T);

if max_dif_T < 1e-4
    disp(['TT = inv(T), son iguales']);
else
    disp('Hay diferencias significativas entre TT e inv(T)');
end


% --- Obtención de la matriz T en una transformación proyectiva ---
disp('Obtención de la matriz T en una transformación proyectiva')

xy_proy = [200 400 600 300; 375 125 375 200];
uv_proy = [87 319 600 200; 445 69 375 180];

T_proy = get_proy(xy_proy, uv_proy);

% Mostrar la matriz T
disp('La matriz de transformacion proyectiva T es:');
disp(T_proy);

TT_proy = get_proy(uv_proy, xy_proy);
disp('La matriz de transformacion proyectiva inversa iT es:');
disp(TT_proy);

T_proy_inv = inv(T_proy);
disp('La inversa inv(T) de la matriz T es:');
disp(T_proy_inv);

% Comparar iT con inv(T)
dif_T_proy = abs(TT_proy - T_proy_inv);
max_dif_T_proy = max(dif_T_proy(:));
fprintf('Diferencia máxima entre iT e inv(T): %e\n', max_dif_T_proy);
if max_dif_T_proy < 1e-4
    disp('iT = inv(T), son iguales');
else
    disp('Hay diferencias significativas entre iT e inv(T)');
end

% Comprobación
uv2_proy = aplica_T(xy_proy, T_proy);
dif_uv_proy = uv_proy - uv2_proy;
disp('Diferencias entre coordenadas deseadas y obtenidas:');
disp(dif_uv_proy);
desv_std_proy = std2(dif_uv_proy(:));
fprintf('Desviación estándar de las diferencias: %e\n', desv_std_proy);


% --- Aplicación de la matriz T de una transformación para deformar una imagen ---
disp('Aplicación de la matriz T de una transformación para deformar una imagen')
im = imread('foto.jpg');

T_ejemplo = [3.5 1.2 -1000; 1.0 4.5 -500; 0.0005 0.003 1];

[im_deformada, RU, RV] = warp_img(im, T_ejemplo);

figure;
imshow(im_deformada);
title('Imagen deformada');

[N_def, M_def, ~] = size(im_deformada);
fprintf('Tamaño de la imagen deformada: %d x %d\n', M_def, N_def);
fprintf('Rango de coordenadas u: [%d, %d]\n', min(RU), max(RU));
fprintf('Rango de coordenadas v: [%d, %d]\n', min(RV), max(RV));
