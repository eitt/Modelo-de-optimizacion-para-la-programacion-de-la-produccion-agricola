%Este es el c�digo para solucionar de manera exacta un modelo de
%optimizaci�n multiobjetivo (como aproximaci�n se trabajan las dos
%funciones de manera independiente), el cual hace parte del proyecto de
%investigaci�n titulado: Modelo de optimizaci�n multiobjetivo para la
% programaci�n de la producci�n agr�cola a peque�a escala en Santander,
% Colombia, desarrollado por Leonardo Talero, bajo la direcci�n de Henry
% Lamos & Edwin Garavito.

% El Script se divide en 3 partes que se describen a continuaci�n:

% 1)______________________________________________________________________
% Se ingresan los par�metros del modelo, los cuales dependen de
% instancias generadas, por tanto, existir�n m�ltiples archivos .mat, los
% cuales tienen almacenadas las inctancias y ser�n cargadas al principio
% del documento. Adem�s, en esta secci�n se construyen los vectores que se
% utilizar�n para generar las restricciones. (vectores b)
% 2)______________________________________________________________________
% Se construyen las matrices de restricci�n, el problema ha sido abordado
% como un modelo de optimizaci�n entera mixta, en esta secci�n se construye
% tambi�n el vector soluci�n.
% 3)______________________________________________________________________
% Se eval�a el modelo, mediante la aplicaci�n de una funci�n ya existente
% en MATLAB, se generan las salidas del modelo.

%Limpieza de variables
clearvars
clc
%Contador del tiempo de c�mputo
tic
% ========================================================================
%=============================SECCI�N # 1=================================
% ========================================================================

%Cargar par�metros de los productos: Duraci�n promedio de maduraci�n,
%periodos cuando se puede recoger el producto, etc. Todos los productos
%est�n organizados alfab�ticamente en orden descendentes (filas) con el fin
%de evitar recurrir a par�metros como el nombre o un identificador.

% load('productos_parametros.mat')
% Conjunto_s: Vector con los instantes donde puede sembrarse cada producto.
% q:          Familia bot�nica a la que pertenece cada producto.
% Ni:         N�mero de periodos que se demora el producto k en "Madurar".
% Pkt:        Precio de venta del producto K en el instante t
% Rkl:        Rendimiento en Kg/m2 de cada producto (col) por lote (fil)
% Al:         �rea en metros de cada lote
% Gv          Grupo de productos sustitutos para satisfacer una demanda

%====Variables para la instancia de prueba

% productos = datasample(1:21,3,'Replace',false);%Selecciono al azar 3 productos
% lotes=datasample(1:length(Rkl(:,1)),2,'Replace',false);%Selecciono al azar 2 lotes
% q=q(productos);
% Conjunto_s=Conjunto_s(productos,:);
% Ni=Ni(productos);
% Pkt=Pkt(productos,:);
% Rkl=Rkl(lotes,productos);
% Al=Al(lotes');
load('instancia_prueba.mat')
B=100000;
%Cantidad de periodos
T=96;
% Cantidad de productos
K=length(Ni);
% Cantidad de lotes
L=length(Al);
% Estimar la matriz de Covarianza de los precios;
Covkkp=cov((Pkt)');
% Se crea el conjunto de instantes de recogida.
Conjunto_r=Conjunto_s*0;
Restriccion_re=zeros(1,T*K*L);
% Se calcula el conjunto de recogida para cada producto:

% Bucle que recorre cada producto
for k=1:K
    % Ubicar los posibles periodos de siembra, posterior a ello, calcular
    % el periodo de madurez Ni, de tal manera que re=s+Ni para cada
    % producto. Como el conjunto de siembra y pro tanto el de recogida
    % tiene dimensiones k*T, estas deben ser transformadas a un vector, con
    % el fin de utilizarlas como restricciones.
    Conjunto_r(k,find(Conjunto_s(k,:)==1)+Ni(k))=1;
    % El conjunto de recogida es almacenado para un �nico lote y debe ser
    % completado en un bucle de lotes (para evitar la creaci�n de bucles
    % extras)
    Restriccion_re(((T*(k-1))+1:(T*(k-1))+T))=Conjunto_r(k,1:T);
end
% Se acota el conjunto de recogida para futuros c�lculos
Conjunto_r=Conjunto_r(:,1:T);
% Se crea el conjunto para determinar el �rea m�xima de siembra para cada
% periodo, este conjunto funciona como restricci�n.
Restriccion_Al=zeros(1,T*K*L);
for l=1:L
    % Se almacena el valor del �rea Al(l) para cada conjunto de tama�o k*T
    Restriccion_Al(((T*K)*(l-1)+1:(T*K*l)))=Al(l);
    % Se aprovecha el bucle de lotes para completar el vector de recogidas
    Restriccion_re(((T*K)*(l-1)+1:(T*K*l)))=Restriccion_re(1:(T*K));
end
%

% ========================================================================
%=============================SECCI�N # 2=================================
% ========================================================================

% ========================================================================
% Familia de restricciones #1: No negatividad
% ========================================================================
% Como son restricciones de no negatividad y se trabaja con la funci�n
% intlinprog, la no negatividad se modelo mediante el vector LB (Lower
% Bound, L�mite Inferior).
% ========================================================================

% ========================================================================
% Familia de restricciones #2: Periodos de recogida
% ========================================================================

% Se genera una matriz identidad en las variables Y (dicot�mica de recogida
% del producto):
A1=eye(T*K*L,T*K*L*4);

% El vector b1 es la transpuesta del vector Restriccion_re
b1=Restriccion_re';
% ========================================================================

% La matriz relaciona las variables de decisi�n Y (dicot�mica) con Z
% (cont�nua) a partir de la premultiplicaci�n por un n�mero muy grande, por
% tanto, la matriz de restricciones se compone de dos diagonales y un
% espacio vac�o:
A2=cat(2,eye(T*K*L,T*K*L)*(-B),zeros(T*K*L,T*K*L),eye(T*K*L,T*K*L),zeros(T*K*L,T*K*L));

% El vector b2 es igual a un vector 0
b2=b1*0;
% ========================================================================


% ========================================================================
% Familia de restricciones #3: �rea destinada a siembra
% ========================================================================

% La matriz indica el �rea m�xima destinada a cda cultivo (variable Z), por
% tanto, se modela como una diagonal igualada a la restricci�n:
A3=cat(2,zeros(T*K*L,T*K*L),zeros(T*K*L,T*K*L),eye(T*K*L,T*K*L),zeros(T*K*L,T*K*L));

% El vector b3 es la transpuesta del vector Restriccion_Al
b3=Restriccion_Al';
% ========================================================================

% ========================================================================
% Familia de restricciones #4: Satisfacci�n de la demanda
% ========================================================================
% Para el caso particular, se considera que existen productos sustitutos y,
% por tanto, estos pueden satisacer una misma demanda. Para agrupar los
% productos se cuenta con el par�metro Gv, el cual indica a cu�l familia de
% productos sustitutos pertenecen los productos k, como es posible que no
% existan todos las familias a la vez (dependiendo de la instancia), es
% necesario determinar cu�les familias existen y crear una familia de
% restricciones (matriz) por cada familia:

% Determinar la cantidad de familias
F_Gv = unique(Gv);
% Se actualiza la demanda para la cantidad de familias
Demanda_Gv=Demanda_Gv(F_Gv);

% Como la demanda cambia periodo a periodo, el tama�o de la matriz A4 es
% igual a la cantidad de familias a satisfacer en cada tiempo (en ls filas)
% (length(F_Gv)*T) por la cantidad de variables de decisi�n:

A4=zeros(length(F_Gv)*T,T*K*L*4);
% por su parte, el vectyor b4 tendr� la misma cantidad de filas y se
% subdividir� en tantas demandas como familias Gv existan:
b4=zeros(length(F_Gv)*T,1);
% Como se construye tantas submatrices como familias, se genera un bucle:
% Desde la primera familia hasta la �ltima
for fgv=1:length(F_Gv)
    % Como se debe recorrer todos los productos para determinar si
    % pertenecen o no a la misma familia Gv, se genera un bucle y se asigna
    % una submatriz identidad a aquellos productos que se van a relacionar
    % y una submatriz de ceros a aquellos que no pertenecen, se definen a
    % continuaci�n:
    
    M_identidad=eye(T,T);
    M_ceros=zeros(T,T);
    
    % Adem�s, se asignan los valores de demanda para el vector b4:
    b4(((fgv-1)*T+1:(fgv-1)*T+T))=Demanda_Gv(fgv);
    
    % Para cada uno de los lotes
    for l=1:L
        % Para cada uno de los productos
        for k=1:K
            if Gv(k)==F_Gv(fgv)
                A4(((fgv-1)*T+1:(fgv-1)*T+T),((T*K*L*2)+(T*K)*(l-1)+T*(k-1)+1:(T*K*L*2)+(T*K)*(l-1)+T*(k-1)+T))=M_identidad;
            else
                A4(((fgv-1)*T+1:(fgv-1)*T+T),((T*K*L*2)+(T*K)*(l-1)+T*(k-1)+1:(T*K*L*2)+(T*K)*(l-1)+T*(k-1)+T))=M_ceros;
            end
        end
    end
end
% ========================================================================

% ========================================================================
% Familia de restricciones #5: Evitar solapamiento
% ========================================================================
% Para evitar el solapamiento de los productos, se tiene en cuenta un
% producto k, el cual puede ser recogido en un instante de tiempo t, todos
% los productos (incluyendo al producto k) no pueden ser recogidos mientras
% dure el periodo de madurez del producto k, por tanto, se debe generar una
% submatriz de cargas que tenga como l�mite superior el instante de
% recogida -1 (t-1), y abarque el periodo de madurez del producto K. Por
% comodidad, se genera una matriz general (a�n cuando el producto k no
% pueda ser recogido en el instante t=1). La cantidad de restricciones en
% este caso es igual al n�mero de variables de decisi�n, por tanto, la
% matriz A5 inicial es igual a:

A5=zeros(K*T*L,K*T*L*4);
% Se hace un recorrido para cada lote
for l=1:L
    %     Para evitar excesos de bucles se calcula todo en el primer lote y se
    %     reproduce para los lotes posteriores
    if l==1
        %     Se realiza un recorrido para cada producto
        for k=1:K
            %         Se ubican las filas correspondientes al intervalo de
            %         madurez del producto:
            [row col] = (find(eye(T,T).*Conjunto_r(k,:)==1));
            %         Se define el intervalo de madurez previo a la
            %         recogida de cada producto
            intervalo_madurez=[find(Conjunto_r(k,1:T)==1)-Ni(k);find(Conjunto_r(k,1:T)==1)-1];
            %         Se contruye una matriz para registar los periodos de
            %         madurez
            M_madurez=zeros(T,T);
            %             Se realiza un bucle para cada uno de los posibles
            %             periodos de madurez
            for i=1:length(intervalo_madurez)
                M_madurez(row(i),intervalo_madurez(1,i):intervalo_madurez(2,i))=1;
            end
            
            %            fila1= (T)*(k-1)+(T*K)*(l-1)+1
            %             fila2=(T)*(k-1)+(T*K)*(l-1)+T
            %             col1=1
            %             col2=(T*K)
            %         Se repite el patr�n para todos los productos
            A5((T)*(k-1)+(T*K)*(l-1)+1:(T)*(k-1)+(T*K)*(l-1)+T,1:(T*K))=repmat(M_madurez,1,K);
            %         Se Agrega la diagonal en el producto correspondiente:
            A5((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)=A5((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)+eye(T,T).*Conjunto_r(k,:);
        end
    else
        % Se replica la distribuci�n de variables para todos los lotes
        A5(K*T*(l-1)+1:K*T*(l),K*T*(l-1)+1:K*T*(l))=A5(1:K*T,1:K*T);
    end
end
% Se crea el vector de restricci�n B5, el cual posee valores de 1:
b5=ones(K*T*L,1);
% ========================================================================

% ========================================================================
% Familia de restricciones #6: Rotaci�n de cultivos (k' diferente k)
% ========================================================================
% para la rotaci�n de productos es necesario que el terreno descanse en
% f^kk'o o^kk' periodos, si pertenecen o no a la misma familia bot�nica
% (F). Teniendo en cuenta que, durante el presente caso se parte del
% supuesto que todos los productos deben descansar la misma cantidad de
% periodos, se construye una �nica restricci�n, la cual tiene en
% consideraci�n dos familias de variables de decisi�n: Y (variable binaria
% utilziada para relacionar la recogida del producto) y la variable V
% (binaria que indica la activaci�n y desactivaci�n de restricciones). El
% tama�o total de la matriz A6 es equivalente a 2*K*T*L filas por la K*T*L
% columnas.

A6=zeros(2*T*K*L,4*T*K*L);
% Posteriormente, se construye unas submatrices para cada periodo, estas se
% componen de dos estilos de matrices, la primera es la de tiempos de
% recogida, la cual corresponde a la secci�n de la restriccio�n que incluye
% a la variable Y, premultiplicada por el tiempo y por la familia bot�nica
% del producto. La segunda submatriz tiene en consideraci�n el tiempo menos
% el tiempo Ni de maduraci�n de cada producto. Con el fin de evitar excesos
% de bucles, se realiza un �nico bucle para los lotes:

% para cada lote
for l=1:L
    if l==1
        %         Para cada producto (filas)
        for k1=1:K
            %         se realiza otro bucle en productos para completar la matriz
            %         Para cada producto (columnas)
            for k2=1: K
                
                %             se diligencia la diagonal cuando k1=k2
                if k1==k2
                    %                   Se crea el conjunto o submatriz de recolecci�n
                    M_recoleccion=ones(T,T).*(Conjunto_r(k1,:).*(1:T)*q(k1)).*Conjunto_r(k1,:)';
%                     filas=(T)*(k1-1)+1:(T)*(k1-1)+T;
%                     columnas=(T)*(k1-1)+1:(T)*(k1-1)+T;
                    A6((T)*(k1-1)+1:(T)*(k1-1)+T,(T)*(k1-1)+1:(T)*(k1-1)+T)=M_recoleccion;
                    %                   Avanzando por columnas
                elseif k2>k1
                    %                   Se crea el conjunto o submatriz de siembra
                    M_siembra=-ones(T,T).*(Conjunto_r(k2,:).*((1:T)-Ni(k2))*q(k2)).*Conjunto_r(k2,:)';
%                     filas=(T)*(k1-1)+1:(T)*(k1-1)+T;
%                     columnas=(T)*(k2-1)+1:(T)*(k2-1)+T;
                    A6((T)*(k1-1)+1:(T)*(k1-1)+T,(T)*(k2-1)+1:(T)*(k2-1)+T)=M_siembra;
                else
                    %                   Se crea el conjunto o submatriz de siembra
                    M_siembra=-ones(T,T).*(Conjunto_r(k2,:).*((1:T)-Ni(k2))*q(k2)).*Conjunto_r(k2,:)';
%                     columnas=(T)*(k2-1)+1:(T)*(k2-1)+T
%                     filas=(T)*(k1-1)+1:(T)*(k1-1)+T
                    A6((T)*(k1-1)+1:(T)*(k1-1)+T,(T)*(k2-1)+1:(T)*(k2-1)+T)=M_siembra;
                end
            end
            
        end
    else
        % Se replica la distribuci�n de variables para todos los lotes
        A6(K*T*(l-1)+1:K*T*(l),K*T*(l-1)+1:K*T*(l))=A6(1:K*T,1:K*T);
    end
end
% Posterior a la creaci�n de las restricciones Y, se construye la
% restricci�n de la variable V, la cual va desde 1 hasta K*L*T. Su valor es
% el negativo de un n�mero muy grande

A6(1:K*T*L,K*T*L+1:2*K*T*L)=-B*eye(K*T*L,K*T*L);
% Finalmente, se replica la matriz de carga, pero con valores negativos,
% con el fin de generar cada par de variables excluyentes:

[row,col] = size(A6);
A6(row/2+1:row,:)=-A6(1:row/2,:);

% para construir el vector b6, es neceario determinar en cuales instantes
% ning�n producto puede ser cosechado, con el fin de excluir las
% restricciones, para ello se parte de identificar en el vector Conjunto_r,
% en cu�les posiciones (para los tres productos) no es posible recoger,
% igualar ese valor a cero y escribir -1 para los dem�s casos:

% se crea un vector par almacenar la restricci�n
vector=((1:T)*0)';
% Se iguala a 1 las posiciones donde se pueda recoger alg�n producto:
vector(find(sum(Conjunto_r(:,1:length(Conjunto_r))==1)~=0))=-1;
% Se crea el vector b6
b6=zeros(2*K*T*L,1);
% Se asigna y replica el vector de restricci�n en la primera mitad del
% vector b6
b6(1:K*T*L)=repmat(vector,2*K,1);
% Se agrega el complemento de la retricci�n (sumar valor B)
b6(K*T*L+1:2*K*T*L)=(-b6(1:K*T*L).*B)+b6(1:K*T*L);
% ========================================================================


% ========================================================================
% Familia de restricciones #7 Espacio entre siembra del mismo cultivo
% ========================================================================
% Como caso especial, cuando no se realiza la rotaci�n de un producto, la
% recogida (y respectiva siembra) del siguiente producto no puede
% realizarse en un periodo inferior de descanso de al menos fkk semanas.
% Por tanto, se construye una matriz de restricci�n A7 con tantas filas
% como varialbes Y y que abarca una cantidad de columnas igual a la
% cantidad de variables de decisi�n:

% Se construye el la matriz:
A7=zeros(T*K*L,4*T*K*L);

% Se hace un recorrido para cada lote
for l=1:L
    %     Para evitar excesos de bucles se calcula todo en el primer lote y se
    %     reproduce para los lotes posteriores
    if l==1
        %     Se realiza un recorrido para cada producto
        for k=1:K
            %         Se ubican las filas correspondientes al intervalo de
            %         madurez del producto:
            [row col] = (find(eye(T,T).*Conjunto_r(k,:)==1));
            %         Se define el intervalo de madurez posterior a la
            %         recogida de cada producto m�s la holgura o descanso
%              [val,loc] = min((find(Conjunto_r(k,1:T)==1)+Ni(k)+1), ones(1,(length(find(Conjunto_r(k,1:T)==1)+Ni(k)+1)))*T)
            intervalo_madurez=[find(Conjunto_r(k,1:T)==1)+1;min((find(Conjunto_r(k,1:T)==1)+Ni(k)+1), ones(1,(length(find(Conjunto_r(k,1:T)==1)+Ni(k)+1)))*T)];
            %         Se contruye una matriz para registar los periodos de
            %         madurez
            M_madurez=zeros(T,T);
            %             Se realiza un bucle para cada uno de los posibles
            %             periodos de madurez
            for i=1:length(intervalo_madurez)
                M_madurez(row(i),intervalo_madurez(1,i):intervalo_madurez(2,i))=1;
            end
            
            %            fila1= (T)*(k-1)+(T*K)*(l-1)+1
            %             fila2=(T)*(k-1)+(T*K)*(l-1)+T
            %             col1=1
            %             col2=(T*K)
            %         Se repite el patr�n para todos los productos
            A7((T)*(k-1)+(T*K)*(l-1)+1:(T)*(k-1)+(T*K)*(l-1)+T,1:(T*K))=repmat(M_madurez,1,K);
            %         Se Agrega la diagonal en el producto correspondiente:
            A7((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)=A7((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)+eye(T,T).*Conjunto_r(k,:);
        end
    else
        % Se replica la distribuci�n de variables para todos los lotes
        A7(K*T*(l-1)+1:K*T*(l),K*T*(l-1)+1:K*T*(l))=A7(1:K*T,1:K*T);
    end
end
% Se crea el vector de restricci�n B5, el cual posee valores de 1:
b5=ones(K*T*L,1);
% ========================================================================


% ========================================================================
% Familia de restricciones #6
% ========================================================================


% ========================================================================


% ========================================================================
% Familia de restricciones #6
% ========================================================================

% ========================================================================


% ========================================================================
%=============================SECCI�N # 3=================================
% ========================================================================


toc