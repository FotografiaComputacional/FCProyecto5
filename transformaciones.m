
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

%PARTE NO NECESARIA
function [xy_n, T] = prepara(xy)
    x = xy(1,:); 
    y = xy(2,:);
    mx = mean(x);
    my = mean(y);
    sx = std(x);
    sy = std(y);
    
    %restar a x yy sus medias
    x0 = x - mx;    
    y0 = y - my;

    % dividir por sus desviaciones
    x1 = x0 / sx;     
    y1 = y0 / sy;

    xy_n = [x1; y1];
    
    T = [ 1/sx,    0,   -mx/sx;
        0,    1/sy,  -my/sy;
        0,      0,      1   ];
end

[xy_n, Txy] = prepara(xy);
disp('Matriz de normalización Txy:');
disp(Txy);

[uv_n, Tuv] = prepara(uv);
disp('Matriz de normalización Tuv:');
disp(Tuv);

%------------------------------------------

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


% ------------------------------------------------------------
% Parte 2
disp('------------------------ Parte 2 ------------------------')
% --- Ejemplo de transformaciones no lineales ---
disp('Ejemplo de transformaciones no lineales')

p = 1.15;

% Paso 1
im = imread('pano.jpg');
im = double(im) / 255;
[N, M, numCanales] = size(im);
%R = N;
R=round(N/1);

% Paso 2
destino = zeros(2*R, 2*R, numCanales);
X = zeros(2*R, 2*R);
Y = zeros(2*R, 2*R);

% Paso 3
for v = 1:2*R
    for u = 1:2*R
        u_centrado = u - R;
        v_centrado = v - R;

        r = sqrt(u_centrado^2 + v_centrado^2);
        theta = atan2(v_centrado, u_centrado);
        
        r_norm = r / (2*R);
        r_norm = r_norm^p;
        theta_norm = mod(theta, 2*pi) / (2*pi);
        
        x = 1 + theta_norm * (M - 1);
        y = N - r_norm * (N - 1);
        
        X(v, u) = x;
        Y(v, u) = y;
    end
end

% Paso 4
for c = 1:numCanales
    destino(:,:,c) = interp2(1:M, 1:N, im(:,:,c), X, Y, 'linear', 0);
end

% Paso 5
figure;
imshow(destino);
title('Transformacion de Panorama');



figure;
subplot(1,2,1);
imshow(im);
title('Panorama de Partida');

subplot(1,2,2);
imshow(destino);
title('Resultado Final');

% ------------------------------------------------------------
% Parte 3
disp('------------------------ Parte 3 ------------------------')
disp('Seam Carving')

im = im2double(imread('img.jpg'));
[Ny,Nx,~] = size(im);

% Energia
function E = energia(im)
    filtrog    = fspecial('gaussian',[7 7],1.51);
    w_gauss  = imfilter(im, filtrog,'symmetric');
    dif = abs(im - w_gauss);
    E    = 0.30*dif(:,:,1) + 0.55*dif(:,:,2) + 0.15*dif(:,:,3);
end

E0      = energia(im);                          
E_mean  = mean(E0(:));
E_max   = max(E0(:));
fprintf('Media(E)  = %.4f\n',E_mean);
fprintf('Máx(E)    = %.4f\n',E_max);

figure;
imagesc(E0), axis image
colormap('hot'), colorbar('vert')
title('Mapa de energía E')
truesize

Nx_convertida = round(Ny*3/2);              
n_fuera  = Nx - Nx_convertida;               
fprintf('Quitamos %d columnas (%dx%d  →  %dx%d)\n',...
        n_fuera,Nx,Ny,Nx_convertida,Ny);

% Calcula_M
function M = calcula_M(E)
    [Ny,Nx] = size(E);
    M = Inf*E;              
    M(1,2:Nx-1) = E(1,2:Nx-1);        

    for i = 2:Ny
        for j = 2:Nx-1
            M(i,j) = E(i,j) + min([M(i-1,j-1), M(i-1,j), M(i-1,j+1)]);
        end
    end
end

% los 3x4 valores de M en la esquina inferiorr derecha
M0        = calcula_M(E0);                    
esq      = M0(end-2:end , end-3:end);        % 3 últimas filas × 4 últimas columnas
fprintf('\nM 3×4 última esquina:\n');
for r = 1:3
    fprintf('%10.4f  %10.4f  %10.4f  %10.4f\n', esq(r,:));
end

figure, imagesc(M0), axis image
colormap('jet'), colorbar('vert')
title('Matriz acumulada M (paleta jet)')
truesize

figure
plot(1:size(M0,2), M0(end,:),'LineWidth', 1.3)
title('Energía acumulada en la última fila de M')
xlabel('columna'); ylabel('M(ALT0, :)')

%find_seam
function J = find_seam(M)
    [Ny,Nx] = size(M);
    J = zeros(Ny,1);                 
    [~,J(Ny)] = min(M(Ny,:)); % J(ALT0)=col_final      
    
    %retrcedemos fila a fila
    for i = Ny-1:-1:1
        prev   = J(i+1);
        cols   = prev-1 : prev+1; %indice solo +-1 columnas
        cols   = cols(cols>=1 & cols<=Nx);
        [~,k]  = min(M(i,cols));
        J(i)   = cols(k);
    end
end

J = find_seam(M0);
figure;
imshow(im), hold on
plot(J, 1:Ny, 'g', 'LineWidth',1.3)
title('Costura mínima superpuesta (verde)')

rows = (1:Ny).';                      
cols = J;                             
idx  = sub2ind(size(E0), rows, cols); 
acum = sum(E0(idx));       %energía acumulada
fprintf('\nCostura – 1er punto: (fila 1 , col %d)\n', J(1));
fprintf('Costura – último punto: (fila %d , col %d)\n', Ny, J(end));
fprintf('Energía acumulada = %.4f\n', acum);

% borrar 1 costura
function im2 = remove_seam(im, J)

    [Ny, Nx, Ch] = size(im);
    im2 = zeros(Ny, Nx-1, Ch);

    for i = 1:Ny                    % recorremos filas
        cols = [1:J(i)-1 , J(i)+1:Nx];   % columnas que se mantienen
        for c = 1:Ch              
            im2(i,:,c) = im(i, cols, c);
        end
    end
end

Nx_ = 972;
n_remove = Nx - Nx_;

im_1 = im;                             
for k = 1:n_remove
    E = energia(im_1);                 
    M = calcula_M(E);                   
    J = find_seam(M);                   
    im_1 = remove_seam(im_1,J);       
end

im_2 = imresize(im,[Ny Nx_]);

escala   = Nx_ / size(im,2);   
im_3  = imresize(im, escala); 
faltan   = Ny - size(im_3,1);      % filas que faltan
pad_arr  = floor(faltan/2);                % mitad arriba
pad_aba  = faltan - pad_arr;               % mitad abajo

im_bands = padarray(im_3,[pad_arr 0],0,'pre');   % negro arriba
im_bands = padarray(im_bands,[pad_aba 0],0,'post'); % negro abajo

montaje = [im_1 ; im_2 ; im_3]; 
figure, imshow(montaje)
title('1) Seam–carving   |   2) imresize   |   3) bandas negras')


%horizontal
Ny__ = round(Nx/2);
n_remove = Ny - Ny__;
fprintf('Eliminar %d filas (%dx%d → %dx%d) para 2:1\n\n',...
        n_remove, Ny, Nx, Ny__, Nx);

% primera costura horizontal verde
imT   = permute(im,[2 1 3]);            % transponer 
E0T   = energia(imT);                 
M0T   = calcula_M(E0T);            
Jt    = find_seam(M0T);              

figure, imshow(im), hold on
plot(1:Nx, Jt, 'g', 'LineWidth',1.3)
title('1ª costura horizontal (verde)')

% coordenadas y energía 
fprintf('Primer punto: (fila %d, col 1)\n', Jt(1));
fprintf('Último punto : (fila %d, col %d)\n', Jt(end), Nx);
acum = sum( E0T(sub2ind(size(E0T),(1:Nx).', Jt)) );
fprintf('Energía acumulada = %.4f\n\n', acum);

% eliminar costuras horizontales
im_h = im;
for k = 1:n_remove
    imT = permute(im_h,[2 1 3]);
    E   = energia(imT);
    M   = calcula_M(E);
    Jt  = find_seam(M);
    imT = remove_seam(imT, Jt);
    im_h = permute(imT,[2 1 3]);
end

im_resize_h = imresize(im,[Ny__ Nx]);

figure
subplot(1,2,1), imshow(im_h)
title(sprintf('Seam-horizontal (%dx%d)',size(im_h,1),size(im_h,2)))
subplot(1,2,2), imshow(im_resize_h)
title(sprintf('imresize       (%dx%d)',size(im_resize_h,1),size(im_resize_h,2)))
