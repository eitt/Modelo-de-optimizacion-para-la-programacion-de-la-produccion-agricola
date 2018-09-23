% Este es el código para solucionar de manera exacta un modelo de
% optimización multiobjetivo (como aproximación se trabajan las dos
% funciones de manera independiente), el cual hace parte del proyecto de
% investigación titulado: Modelo de optimización multiobjetivo para la
% programación de la producción agrícola a pequeña escala en Santander,
% Colombia, desarrollado por Leonardo Talero, bajo la dirección de Henry
% Lamos & Edwin Garavito.

% El Script se divide en 3 partes que se describen a continuación:

% 1)______________________________________________________________________
% Se ingresan los parámetros del modelo, los cuales dependen de
% instancias generadas, por tanto, existirán múltiples archivos .mat, los
% cuales tienen almacenadas las inctancias y serán cargadas al principio
% del documento. Además, en esta sección se construyen los vectores que se
% utilizarán para generar las restricciones. (vectores b)
% 2)______________________________________________________________________
% Se construyen las matrices de restricción, el problema ha sido abordado
% como un modelo de optimización entera mixta, en esta sección se construye
% también el vector solución.
% 3)______________________________________________________________________
% Se evalúa el modelo, mediante la aplicación de una función ya existente
% en MATLAB, se generan las salidas del modelo.

%Limpieza de variables
clearvars
clc
%Contador del tiempo de cómputo
tic
% ========================================================================
%=============================SECCIÓN # 1=================================
% ========================================================================

%Cargar parámetros de los productos: Duración promedio de maduración,
%periodos cuando se puede recoger el producto, etc. Todos los productos
%están organizados alfabéticamente en orden descendentes (filas) con el fin
%de evitar recurrir a parámetros como el nombre o un identificador.
% 
% load('productos_parametros.mat')
% Conjunto_s: Vector con los instantes donde puede sembrarse cada producto.
% q:          Familia botánica a la que pertenece cada producto.
% Ni:         Número de periodos que se demora el producto k en "Madurar".
% Pkt:        Precio de venta del producto K en el instante t
% Rkl:        Rendimiento en Kg/m2 de cada producto (col) por lote (fil)
% Al:         Área en metros de cada lote
% Gv          Grupo de productos sustitutos para satisfacer una demanda

% ====Variables para la instancia de prueba
% 
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
T=length(Conjunto_s);
% Cantidad de productos
K=length(Ni);
% Cantidad de lotes
L=length(Al);
% Estimar la matriz de Covarianza de los precios;
% Covkkp=cov((Pkt)');
Covkkp=cov(diff(log(Pkt)',1));
% Determinar la cantidad de variables de decisión para cada una de las
% familias de variables:

% {Y_t }^{l,k}
Cant_Y=K*T*L;
% {V_t }^{l,k}
Cant_V=K*T*L;
% {Z_t }^{l,k}
Cant_Z=K*T*L;
% {_y }^{k,k}
Cant_U=K*K*T*L;
Cant_var=Cant_Y+Cant_V+Cant_Z+Cant_U;

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
    % El conjunto de recogida es almacenado para un único lote y debe ser
    % completado en un bucle de lotes (para evitar la creación de bucles
    % extras)
    Restriccion_re(((T*(k-1))+1:(T*(k-1))+T))=Conjunto_r(k,1:T);
end
% Se acota el conjunto de recogida para futuros cálculos
Conjunto_r=Conjunto_r(:,1:T);
% Se crea el conjunto para determinar el área máxima de siembra para cada
% periodo, este conjunto funciona como restricción.
Restriccion_Al=zeros(1,T*K*L);
for l=1:L
    % Se almacena el valor del área Al(l) para cada conjunto de tamaño k*T
    Restriccion_Al(((T*K)*(l-1)+1:(T*K*l)))=Al(l);
    % Se aprovecha el bucle de lotes para completar el vector de recogidas
    Restriccion_re(((T*K)*(l-1)+1:(T*K*l)))=Restriccion_re(1:(T*K));
end
%

% ========================================================================
%=============================SECCIÓN # 2=================================
% ========================================================================

% ========================================================================
% Familia de restricciones #1: No negatividad
% ========================================================================
% Como son restricciones de no negatividad y se trabaja con la función
% intlinprog, la no negatividad se modelo mediante el vector LB (Lower
% Bound, Límite Inferior).
% ========================================================================

% ========================================================================
% Familia de restricciones #2: Periodos de recogida
% ========================================================================

% Se genera una matriz identidad en las variables Y (dicotómica de recogida
% del producto):
A1=eye(T*K*L,Cant_var);

% El vector b1 es la transpuesta del vector Restriccion_re
b1=Restriccion_re';
% ========================================================================

% La matriz relaciona las variables de decisión Y (dicotómica) con Z
% (contínua) a partir de la premultiplicación por un número muy grande, por
% tanto, la matriz de restricciones se compone de dos diagonales y un
% espacio vacío:
A2=cat(2,eye(Cant_Y,Cant_Y)*(-B),zeros(Cant_V,Cant_V),eye(Cant_Z,Cant_Z),zeros(Cant_Y,Cant_U));

% El vector b2 es igual a un vector 0
b2=b1*0;
% ========================================================================


% ========================================================================
% Familia de restricciones #3: Área destinada a siembra
% ========================================================================

% La matriz indica el área máxima destinada a cda cultivo (variable Z), por
% tanto, se modela como una diagonal igualada a la restricción:
A3=cat(2,zeros(Cant_Y,Cant_Y),zeros(Cant_V,Cant_V),eye(Cant_Z,Cant_Z),zeros(Cant_Z,Cant_U));

% El vector b3 es la transpuesta del vector Restriccion_Al
b3=Restriccion_Al';
% ========================================================================

% ========================================================================
% Familia de restricciones #4: Satisfacción de la demanda
% ========================================================================
% Para el caso particular, se considera que existen productos sustitutos y,
% por tanto, estos pueden satisacer una misma demanda. Para agrupar los
% productos se cuenta con el parámetro Gv, el cual indica a cuál familia de
% productos sustitutos pertenecen los productos k, como es posible que no
% existan todos las familias a la vez (dependiendo de la instancia), es
% necesario determinar cuáles familias existen y crear una familia de
% restricciones (matriz) por cada familia:

% Determinar la cantidad de familias
F_Gv = unique(Gv);
% Se actualiza la demanda para la cantidad de familias
Demanda_Gv=Demanda_Gv(F_Gv);

% Como la demanda cambia periodo a periodo, el tamaño de la matriz A4 es
% igual a la cantidad de familias a satisfacer en cada tiempo (en ls filas)
% (length(F_Gv)*T) por la cantidad de variables de decisión:

A4=zeros(length(F_Gv)*T,Cant_var);
% por su parte, el vectyor b4 tendrá la misma cantidad de filas y se
% subdividirá en tantas demandas como familias Gv existan:
b4=zeros(length(F_Gv)*T,1);
% Como se construye tantas submatrices como familias, se genera un bucle:
% Desde la primera familia hasta la última
for fgv=1:length(F_Gv)
    % Como se debe recorrer todos los productos para determinar si
    % pertenecen o no a la misma familia Gv, se genera un bucle y se asigna
    % una submatriz identidad a aquellos productos que se van a relacionar
    % y una submatriz de ceros a aquellos que no pertenecen, se definen a
    % continuación:
    
    M_identidad=eye(T,T);
    M_ceros=zeros(T,T);
    
    % Además, se asignan los valores de demanda para el vector b4:
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
% submatriz de cargas que tenga como límite superior el instante de
% recogida -1 (t-1), y abarque el periodo de madurez del producto K. Por
% comodidad, se genera una matriz general (aún cuando el producto k no
% pueda ser recogido en el instante t=1). La cantidad de restricciones en
% este caso es igual al número de variables de decisión, por tanto, la
% matriz A5 inicial es igual a:

A5=zeros(Cant_Y,Cant_var);
% Se hace un recorrido para cada lote
for l=1:L
    %     Para evitar excesos de bucles se calcula todo en el primer lote y se
    %     reproduce para los lotes posteriores
    if l==1
        %     Se realiza un recorrido para cada producto
        for k=1:K
            %         Se ubican las filas correspondientes al intervalo de
            %         madurez del producto:
            [row, col] = (find(eye(T,T).*Conjunto_r(k,:)==1));
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
            %         Se repite el patrón para todos los productos
            A5((T)*(k-1)+(T*K)*(l-1)+1:(T)*(k-1)+(T*K)*(l-1)+T,1:(T*K))=repmat(M_madurez,1,K);
            %         Se Agrega la diagonal en el producto correspondiente:
            A5((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)=A5((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)+eye(T,T).*Conjunto_r(k,:);
        end
    else
        % Se replica la distribución de variables para todos los lotes
        A5(K*T*(l-1)+1:K*T*(l),K*T*(l-1)+1:K*T*(l))=A5(1:K*T,1:K*T);
    end
end
% Se crea el vector de restricción B5, el cual posee valores de 1:
b5=ones(Cant_Y,1);
% ========================================================================

% ========================================================================
% Familia de restricciones #6: Rotación de cultivos (k' diferente k)
% ========================================================================
% para la rotación de productos es necesario que el terreno descanse en
% f^kk'o o^kk' periodos, si pertenecen o no a la misma familia botánica
% (F). Teniendo en cuenta que, durante el presente caso se parte del
% supuesto que todos los productos deben descansar la misma cantidad de
% periodos, se construye una única restricción, la cual tiene en
% consideración dos familias de variables de decisión: Y (variable binaria
% utilziada para relacionar la recogida del producto) y la variable V
% (binaria que indica la activación y desactivación de restricciones). El
% tamaño total de la matriz A6 es equivalente a 2*K*T*L filas por la K*T*L
% columnas.

A6=zeros(Cant_Y*2,Cant_var);
% Posteriormente, se construye unas submatrices para cada periodo, estas se
% componen de dos estilos de matrices, la primera es la de tiempos de
% recogida, la cual corresponde a la sección de la restriccioón que incluye
% a la variable Y, premultiplicada por el tiempo y por la familia botánica
% del producto. La segunda submatriz tiene en consideración el tiempo menos
% el tiempo Ni de maduración de cada producto. Con el fin de evitar excesos
% de bucles, se realiza un único bucle para los lotes:


% para cada Lote
for l=1:L
%     Se hace la solución para el primer lote y luego se replica
    if l==1
%         para cada producto k
        for k1=1:K
%             Para cada producto k'
            for k2=1:K
%                 Se calcula la diagonal de la matriz
                if k1==k2                    
                    % Se crea el conjunto o submatriz de recolección
                    M_recoleccion=zeros(T,T);
%                     Para cada periodo
                    for t=1:T
%                         Se diligencia el producto k' sembrado en el
%                         instante t' y recogido en t
% el rango aplica sólo si t>1
                        if t>1
%                             el rango va desde mínimo t=1, hasta el
%                             instante previo a sembrar el product k en el
%                             instante t
                            RANGO_ANTES=(max(1,t-Ni(k2)):t-1);
%                             Se almacena el conjunto en la submatriz
                            M_recoleccion(t,RANGO_ANTES)=-(RANGO_ANTES-Ni(k2))*q(k2);
                        end
%                         Se diligencia el producto k, recogido en el
%                         instante t+Ni(k)
                        RANGO_DIAGONAL=(t:min(T,t+Ni(k1)-1));
%                         Se almacena el conjunto en la submatriz
                        M_recoleccion(t,RANGO_DIAGONAL)=(RANGO_DIAGONAL)*q(k1);
%                         Se calcula el conjunto para el producto k'
%                         recogido en el instante t', posterior a t t'>t
% Se verifica que no exceda el tiempo máximo del subconjunto
                        if t+Ni(k1)<T
%                             se crea el periodo de maduración.
                            RANGO_DESPUES=(t+Ni(k1):min(T,t+Ni(k1)+Ni(k2)-1));
%                             Se almacena el conjunto en la submatriz
                            M_recoleccion(t,RANGO_DESPUES)=-(RANGO_DESPUES-Ni(k2))*q(k2);
                        end
                    end
%                     se almacena el subconjunto en la matriz de cargas
                    A6((T)*(k1-1)+1:(T)*(k1-1)+T,(T)*(k2-1)+1:(T)*(k2-1)+T)=M_recoleccion;
%                 se replica el proceso cuando k'>k (posición en el conjunto K)    
                elseif k2>k1

                    M_recoleccion=zeros(T,T);
                    for t=1:T
                        %                         Diagonal antes
                        if t>1
                            RANGO_ANTES=(max(1,t-Ni(k2)):t-1);
                            M_recoleccion(t,RANGO_ANTES)=-(RANGO_ANTES-Ni(k2))*q(k2);
                        end
                        RANGO_DIAGONAL=(t:min(T,t+Ni(k1)-1));
                        M_recoleccion(t,RANGO_DIAGONAL)=(RANGO_DIAGONAL)*q(k1);
                        if t+Ni(k1)<T
                            RANGO_DESPUES=(t+Ni(k1):min(T,t+Ni(k1)+Ni(k2)-1));
                            M_recoleccion(t,RANGO_DESPUES)=-(RANGO_DESPUES-Ni(k2))*q(k2);
                        end
                    end
                    A6((T)*(k1-1)+1:(T)*(k1-1)+T,(T)*(k2-1)+1:(T)*(k2-1)+T)=M_recoleccion;
                    %                 se replica el proceso cuando k'<k (posición en el conjunto K)
                else
                    % Se crea el conjunto o submatriz de recolección
                    M_recoleccion=zeros(T,T);
                    for t=1:T
                        %                         Diagonal antes
                        if t>1
                            RANGO_ANTES=(max(1,t-Ni(k2)):t-1);
                            M_recoleccion(t,RANGO_ANTES)=-(RANGO_ANTES-Ni(k2))*q(k2);
                        end
                        RANGO_DIAGONAL=(t:min(T,t+Ni(k1)-1));
                        M_recoleccion(t,RANGO_DIAGONAL)=(RANGO_DIAGONAL)*q(k1);
                        if t+Ni(k1)<T
                            RANGO_DESPUES=(t+Ni(k1):min(T,t+Ni(k1)+Ni(k2)-1));
                            M_recoleccion(t,RANGO_DESPUES)=-(RANGO_DESPUES-Ni(k2))*q(k2);
                        end
                    end
                    A6((T)*(k1-1)+1:(T)*(k1-1)+T,(T)*(k2-1)+1:(T)*(k2-1)+T)=M_recoleccion;
                    
                end
            end
        end
    else
       % Se replica la distribución de variables para todos los lotes
        A6(K*T*(l-1)+1:K*T*(l),K*T*(l-1)+1:K*T*(l))=A6(1:K*T,1:K*T);  
    end
end

% Se construye la matriz de cargas de la variable V
A6(1:K*T*L,K*T*L+1:2*K*T*L)=-B*eye(K*T*L,K*T*L);
% Finalmente, se replica la matriz de carga, pero con valores negativos,
% con el fin de generar cada par de variables excluyentes:

[row,~] = size(A6);
A6(row/2+1:row,:)=-A6(1:row/2,:);

% Posterior a la creación de las restricciones Y, se construye la
% restricción de la variable V, la cual va desde 1 hasta K*L*T. Su valor es
% el negativo de un número muy grande

A6(1:K*T*L,K*T*L+1:2*K*T*L)=-B*eye(K*T*L,K*T*L);

% para la construcción del vector restricción, se tiene en cuenta la
% diferencia entre grupos F^k,k' o o^k,k'; en este caso, ambas igual a 1,
% por tanto:
b6=cat(1,-ones(K*L*T,1),(-ones(K*L*T,1)+B));
% ========================================================================

% ========================================================================
% Familia de restricciones #7 Espacio entre siembra del mismo cultivo
% ========================================================================
% Como caso especial, cuando no se realiza la rotación de un producto, la
% recogida (y respectiva siembra) del siguiente producto no puede
% realizarse en un periodo inferior de descanso de al menos fkk semanas.
% Por tanto, se construye una matriz de restricción A7 con tantas filas
% como varialbes Y y que abarca una cantidad de columnas igual a la
% cantidad de variables de decisión:

% Se construye el la matriz:
A7=zeros(Cant_Y,Cant_var);

% Se hace un recorrido para cada lote
for l=1:L
    %     Para evitar excesos de bucles se calcula todo en el primer lote y se
    %     reproduce para los lotes posteriores
    if l==1
        %     Se realiza un recorrido para cada producto
        for k=1:K
            %         Se ubican las filas correspondientes al intervalo de
            %         madurez del producto:
            [row, col] = (find(eye(T,T).*Conjunto_r(k,:)==1));
            %         Se define el intervalo de madurez posterior a la
            %         recogida de cada producto más la holgura o descanso
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
            %         Se repite el patrón para todos los productos
            A7((T)*(k-1)+(T*K)*(l-1)+1:(T)*(k-1)+(T*K)*(l-1)+T,1:(T*K))=repmat(M_madurez,1,K);
            %         Se Agrega la diagonal en el producto correspondiente:
            A7((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)=A7((T)*(k-1)+1:(T)*(k-1)+T,(T)*(k-1)+1:(T)*(k-1)+T)+eye(T,T).*Conjunto_r(k,:);
        end
    else
        % Se replica la distribución de variables para todos los lotes
        A7(K*T*(l-1)+1:K*T*(l),K*T*(l-1)+1:K*T*(l))=A7(1:K*T,1:K*T);
    end
end
% Se crea el vector de restricción b7, el cual posee valores de 1 en los
% periodos donde es posible realizar la recolección:

b7=b1;
% ========================================================================


% ========================================================================
% Familia de restricciones #8: Proporción de áreas recogidas (para cada
% par)
% ========================================================================

% La octava familia de restricciones relacionan dos variables de decisión:
% Uk,k' (la cual premultiplica l amatriz de varianza-covarianza) con la
% variable de decisión Z l,k,t y el rendimiento (calculando así el volumen 
% de producción) Como la variable U no posee índice t, la matriz de cargas
% para Z tiene (k^2 )*T filas (para estimar cada interacción) y L*K*T
% columnas.

A8=zeros(Cant_U,Cant_var);

% se contruye el vector b8 de restricciones, el cual relaciona las áreas
% máximas donde está cultivado cada producto:
b8=zeros(Cant_U,1);
% para reducir la cantidad de bucles, la asignación de las cargas comienza
% por los productos (conjuntos que más se repiten) y posteriormente se
% recalcula el valor para el respectivo lote:

for k1=1:K
    for k2=1:K
        % se genera el bucle para cada lote:
        for l=1:L
            if k1==k2
            %         Se crea una matriz diagonal a la cual se premultiplican las
            %         áreas:

            M_areas=eye(T,T)*(Rkl(l,k2)*Al(l)*Rkl(l,k1)*2);
            %         Se asigna el valor de la submatriz en la matriz de cargas A:
            A8((K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+1:(K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+T,Cant_Y+Cant_V+(K*T)*(l-1)+T*(k1-1)+1:Cant_Y+Cant_Z+(K*T)*(l-1)+T*(k1-1)+T)=M_areas;
            else
                M_areask1=eye(T,T)*(Rkl(l,k1)*Al(l)*Rkl(l,k1));
                M_areask2=eye(T,T)*(Rkl(l,k1)*Al(l)*Rkl(l,k2));
                A8((K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+1:(K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+T,Cant_Y+Cant_V+(K*T)*(l-1)+T*(k1-1)+1:Cant_Y+Cant_Z+(K*T)*(l-1)+T*(k1-1)+T)=M_areask1;
                A8((K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+1:(K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+T,Cant_Y+Cant_V+(K*T)*(l-1)+T*(k2-1)+1:Cant_Y+Cant_Z+(K*T)*(l-1)+T*(k2-1)+T)=M_areask2;
            end
        end
    end
end

% Respecto a la variable U, se crea una diagonal positiva
A8(:,Cant_Y+Cant_V+Cant_Z+1:Cant_var)=-eye(Cant_U,Cant_U);

% Se diligencia el vector b8

for l=1:L
    b8((Cant_U/L)*(l-1)+1:(Cant_U/L)*(l-1)+(Cant_U/L))=Al(l);
end
b8=b8.^2;
col1=b8*0;
col2=col1;

for l=1:L
    for k1=1:K
        recorrido1=(K*K*T)*(l-1)+(K*T)*(k1-1);
%         b8(recorrido1+1:recorrido1+K*T)=b8(recorrido1+1:recorrido1+K*T)*Rkl(l,k1);
        col1(recorrido1+1:recorrido1+K*T)=Rkl(l,k1);
        for k2=1:K
            recorrido2=recorrido1+T*(k2-1);
%             b8(recorrido2+1:recorrido2+T)=b8(recorrido2+1:recorrido2+T)*Rkl(l,k2);
            col2(recorrido2+1:recorrido2+T)=Rkl(l,k2);
        end
    end    
end


% ========================================================================


% ========================================================================
% Familia de restricciones #9: Proporción de Áreas U^{k,k'} Cantidad de
% producto k, recogido en cada instante t.
% ========================================================================
% Se contruye la matriz de cargas para la familia de restricciones
A9=zeros(Cant_U,Cant_var);
% Se construye el vector de restricciones b9:
b9=zeros(Cant_U,1);

% para reducir la cantidad de bucles, la asignación de las cargas comienza
% por los productos (conjuntos que más se repiten) y posteriormente se
% recalcula el valor para el respectivo lote:



for l=1:L
  A9((Cant_U/L)*(l-1)+1:(Cant_U/L)*(l-1)+(Cant_U/L),Cant_Y+Cant_V+(K*T)*(l-1)+1:Cant_Y+Cant_V+(K*T)*(l-1)+K*T)=repmat(eye(K*T),K,1);
  b9((Cant_U/L)*(l-1)+1:(Cant_U/L)*(l-1)+(Cant_U/L))=Al(l);
end

A9(:,Cant_Y+Cant_V+1:Cant_U)=A9(:,Cant_Y+Cant_V+1:Cant_U).*b9.*col1.*col2*-1;
A9(:,Cant_var-Cant_U+1:Cant_var)=eye(Cant_U,Cant_U);
b9=b9*0;
% ========================================================================

% ========================================================================
% Familia de restricciones #10: Proporción de Áreas U^{k,k'}?Cantidad de
% producto k, recogido en cada instante t.
% ========================================================================
% Se contruye la matriz de cargas para la familia de restricciones
A10=zeros(Cant_U,Cant_var);
% Se construye el vector de restricciones b9:
b10=zeros(Cant_U,1);

% para reducir la cantidad de bucles, la asignación de las cargas comienza
% por los productos (conjuntos que más se repiten) y posteriormente se
% recalcula el valor para el respectivo lote:


for k1=1:K
    for k2=1:K
        % se genera el bucle para cada lote:
        for l=1:L
            if k1==k2
            %         Se crea una matriz diagonal a la cual se premultiplican las
            %         áreas:
            M_areas=-eye(T,T)*(Rkl(l,k1)*Al(l)*Rkl(l,k2));
            %         Se asigna el valor de la submatriz en la matriz de cargas A:
            A10((K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+1:(K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+T,Cant_Y+Cant_V+(K*T)*(l-1)+T*(k1-1)+1:Cant_Y+Cant_Z+(K*T)*(l-1)+T*(k1-1)+T)=M_areas;
            else
%                 M_areask1=eye(T,T)*(-Rkl(l,k2)*Al(l));
                M_areask1=-eye(T,T)*(Rkl(l,k1)*Al(l)*Rkl(l,k2));
                A10((K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+1:(K*K*T)*(l-1)+(K*T)*(k1-1)+T*(k2-1)+T,Cant_Y+Cant_V+(K*T)*(l-1)+T*(k1-1)+1:Cant_Y+Cant_Z+(K*T)*(l-1)+T*(k1-1)+T)=M_areask1;
                
            end
        end
    end
end
A10(:,Cant_Y+Cant_V+Cant_Z+1:Cant_var)=eye(Cant_U,Cant_U);
% ========================================================================



% ========================================================================
%=============================SECCIÓN # 3=================================
% ========================================================================
% 
% % Se integran todas las restricciones, comenzando por la matriz de cargas:
% A=cat(1,A1,A2,A3,A4,A5,A6,A7,A8,A9);
% % Seguida del vector de desigualdades
% b=cat(1,b1,b2,b3,b4,b5,b6,b7,b8,b9);
%
% Se integran todas las restricciones, comenzando por la matriz de cargas:
A=cat(1,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10);
% Seguida del vector de desigualdades
b=cat(1,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10);



cova=zeros(1,Cant_U/L);
for k1= 1:K
    for k2=1:K
        recorrido=((K*T)*(k1-1))+T*(k2-1);
        cova(recorrido+1:recorrido+T)=Covkkp(k1,k2);
    end
end
cova=repmat(cova,1,L);

f2=cat(1,zeros(Cant_var-Cant_U,1),cova');

f = -repmat(Pkt(:),L,1);
f=cat(1,zeros(Cant_Y+Cant_V,1),f,zeros(Cant_U,1));
intcon = 1:Cant_Y+Cant_V;

Aeq = [];
beq = [];

lb = zeros(length(f),1);
% for l=1:L
%     b88((Cant_U/L)*(l-1)+1:(Cant_U/L)*(l-1)+(Cant_U/L))=Al(l);
% end
% ub =  cat(1,ones(Cant_Y+Cant_V,1),repmat(Inf,Cant_Z,1),b88');
ub =  cat(1,ones(Cant_Y+Cant_V,1),Inf(Cant_Z+Cant_U,1));
x=intlinprog(f,intcon',A,b,[],[],lb,ub);
x2=intlinprog(f2,intcon',A,b,[],[],lb,ub);
ponderacion=0.95;
f3=f*(1-ponderacion)+f2*ponderacion;
x3=intlinprog(f3,intcon',A,b,[],[],lb,ub);
% intcon2=[]
% x=intlinprog(f,intcon2,A,b,[],[],lb,ub)
toc