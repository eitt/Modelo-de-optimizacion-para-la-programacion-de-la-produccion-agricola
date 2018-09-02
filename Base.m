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
% del documento.
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

% ========================================================================
%=============================SECCI�N # 2=================================
% ========================================================================


% ========================================================================
%=============================SECCI�N # 3=================================
% ========================================================================

Covkkp=cov((Pkt)')%Covarianza de los precios;

toc