// Se define la función objetivo, su gradiente y su Hessiano

function [f,g,H] = central(x);
  
E1=0.0578082;
E2=0.0016517;
E3=0.0412916;

f = log(x(1)^100)*E1+x(2)^2*E2+x(3)*E3;

g = zeros(3,1);
g(1) = 100/x(1)*E1;
g(2) = 2*x(2)*E2;
g(3) = E3;

H = zeros(3,3);
H(1,1) = -100/x(1)^2*E1;
H(2,2) = 2*E2;
H(3,3)=  0;

endfunction //Fin de la funcion central

// Se define la restriccion de igualdad y su Jacobiano

function [c,A] = demanda(x);

c = x(1,1)+x(2,1)+x(3,1)-1000;
A = [1 1 1];

endfunction //Fin de la funcion demanda

// Se implementa una función con el método de puntos interiores

function [x,f] = interior(x,lb,ub);

  // METODO DE PUNTOS INTERIORES

  // Variables duales de las cotas inferiores

  n=length(x);

  zlb=zeros(n,1);
  for i=1:n
   zlb(i)=1/(x(i) - lb(i));
  end;

  // Variables duales de las cotas superiores

  zub=zeros(n,1);
  for i=1:n
   zub(i)=1/(ub(i) - x(i));
  end;

  // Calculamos el valor inicial de los multiplicadores

  [f,g,H] = central(x);
  [c,A] = demanda(x);
  m=length(c);
  lam=A'\(g-zlb+zub);

  // Inicializamos el contador de iteraciones

  k=0;

  // Condiciones de complementariedad

  colb=zeros(n,1);
  for i=1:n
    colb(i)=(x(i) - lb(i))'*zlb(i);
  end;

  coub=zeros(n,1);
  for i=1:n
    coub(i)=(ub(i) - x(i))'*zub(i);
  end;

  // Inicializamos el parámetro de barrera mu

  mu=1;

  // Comienza el bucle principal de la funcion interior

  while norm([(g-A'*lam-zlb+zub)' c' colb' coub'])>1.e-5
  
  // Construimos la matriz del sistema
  
    Xlb=zeros(n,n);
    for i=1:n
      Xlb(i,i)=1/(x(i)-lb(i));
    end;
    Xub=zeros(n,n);
    for i=1:n
      Xub(i,i)=1/(ub(i)-x(i));
    end;
    Zlb=diag(zlb);
    Zub=diag(zub);
    Hd=H+Xlb*Zlb+Xub*Zub;
  
    // Hacemos la matriz H definida positiva
  
    [U,D]=spec(Hd);
    for i=1:n
      D(i,i)=max(abs(D(i,i)),1.e-3);
    end
    Hd=U*D*U';
    K=[Hd A'
       A 0];
   
    // Calculamos el lado derecho del sistema de ecuaciones
   
    der=[-g+mu*diag(Xlb)-mu*diag(Xub)+A'*lam
         -c];
   
    // Resolvemos el sistema de ecuaciones
   
    d=K\der;
   
    // Calculamos las direcciones de movimiento de las variables primales
   
    dx=d(1:n,1);
    dl=-d(n+1:n+m,1);
   
    // Direccion de movimiento de las variables duales

    dzlb = mu*diag(Xlb)-zlb-Xlb*Zlb*dx;
    dzub = mu*diag(Xub)-zub+Xub*Zub*dx;
   
    // Calculamos la longitud de paso de las variables duales

    tau=0.99995;                        
    azlb=1;
    for i=1:n'
     if dzlb(i)<0;
      azlb=min(tau*(-zlb(i))/dzlb(i),azlb); 
     end;
    end;
                       
    azub=1;
    for i=1:n'
     if dzub(i)<0;
      azub=min(tau*(-zub(i))/dzub(i),azub); 
     end;
    end;
    alfaz=min([azlb azub]);
  
    // Actualizamos las variables duales  
  
    zlb=zlb+alfaz*dzlb;
    zub=zub+alfaz*dzub;

    // Calculamos la longitud de paso de las variables primales   
   
    alb=1;
    aub=1;
            
    for i=1:n
     if dx(i)<0;
       alb=min(tau*(lb(i)-x(i))/dx(i),alb); 
     end;
    end;

    for i=1:n
     if dx(i)>0;
      aub=min(tau*(ub(i)-x(i))/dx(i),aub); 
     end;
    end;
    
    alfax=min([alb aub]);

    // Actualizamos las variables primales y los multiplicadores

    x=x+alfax*dx;
    lam=lam+alfax*dl;

    // Actualizamos los valores del problema

    [f,g,H]= central(x);
    [c,A] = demanda(x);

    // Actualizamos el parámetro barrera

    mu=0.9*mu;

    // Actualizamos el contador de iteraciones

    k=k+1;

    // Actualizamos las condiciones de complementariedad

    for i=1:n
      colb(i)=(x(i) - lb(i))'*zlb(i);
    end;
    for i=1:n
      coub(i)=(ub(i) - x(i))'*zub(i);
    end;
 
    // Limitamos el numero de iteraciones a 1000
  
    if k>1000 then
      disp('Demasiadas iteraciones')  
      break
    end
 
  end

endfunction //Fin de la funcion interior

// PROGRAMA PRINCIPAL

// Se define el vector de cotas inferiores de cada calentador

lb=[1 0 0]';

// Se define el vector de cotas superiores de cada calentador

ub=[700 500 500]';

// Se toma un punto inicial contenido entre las cotas

x=[300 200 300]';

// Resolvemos el problema

[x,f]=interior(x,lb,ub);

// Sacamos por pantalla la solución del problema

disp(f,'El coste óptimo del GW/hora (en miles de euros) es:')
disp(x(1),'Los MW/hora a producir por el primer calentador son:')
disp(x(2),'Los MW/hora a producir por el segundo calentador son:')
disp(x(3),'Los MW/hora a producir por el tercer calentador son:')

// FIN DEL PROGRAMA PRINCIPAL
