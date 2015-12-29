% Representa las funciones de 
% coste del modelo mejorado

% Límites con las especificaciones 
% de producción de energía
x1 = [1:0.01:700]';
x2 = [0:0.01:500]';
x3 = [0:0.01:500]';

% Constantes de la función objetivo 
E1 = 0.0578082;
E2 = 0.0016517;
E3 = 0.0412916;

% Representamos cada función individualmente,
% una encima de la otra
scf(3);
clf(3); 
plot(x1, E1*log(x1.^100), 'r');
plot(x2, E2*x2.^2, 'b');
plot(x3, E3*x3, 'black');

% Representamos la leyenda para identificar cada función
h1 = legend(['0.0578082*log(x^100)';'0.0016517*x^2';'0.0412916*x']);
