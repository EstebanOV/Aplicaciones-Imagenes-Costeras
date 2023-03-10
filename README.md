# Aplicaciones Imágenes Costeras 

![](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/main/Portada-Banner.jpg)

Este repositorio contiene los códigos utilizados para la obtención del campo de velocidades medias en la zona de la rompiente a partir de imágenes grabadas desde un dron. 

Estos códigos fueron generados como parte de mi trabajo de memoria para optar el título de Ingeniero Civil, y se enmarca dentro del proyecto *"FONDECYT 1170415-Quantification of Two Dimensional Wave Breaking Dissipation in the Surf Zone"* liderado por el profesor Patricio Catalán del Departamento de Ingeniería Civil de la Universidad Técnica Federico Santa María.


El repositorio se encuentra organizado en 3 grupos: 

- **Rectificación:** Códigos utilizados para la rectificación de imágenes. 

- **Batimetría:**  Códigos de cBathy adaptados al proyecto para la obtención de batimetría en la zona.

- **Corrientes:** Códigos generados para la obtención de corrientes medias longitudinales y transversales. 

Cada uno de estos grupos de códigos se puede utilizar de forma independiente y adaptar a las condiciones requeridas.

<!-- TABLE OF CONTENTS -->
## Indice
- [1. Rectificación](#1-rectificacion)
- [2. Batimetría cBathy](#2-batimetria_cbathy)
- [3. Corrientes](#3-corrientes)
- [4. Documentación](#4-documentacion)
- [5. Contacto](#5-contacto)


<!-- ABOUT THE PROJECT -->
## 1. Rectificacion
Estos scripts de matlab corresponden a una modificación y adaptación del paquete de códigos del CIRN (Coastal Imaging Research Network). Este paquete se encuentra almacenado en el siguiente repositorio: [UAV-Processing-Toolbox](https://github.com/Coastal-Imaging-Research-Network/UAV-Processing-Toolbox "UAV-Processing-Toolbox").

Previo al proceso de la ortorectificación, se deben tener en cuenta algunas recomendaciones en la grabación de imágenes en terreno. En el siguiente archivo se resumen algunos puntos a tener en cuenta: 

[Guía de Pre-Rectificación](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/61438c70ad05e72d21f3ade8688130c404e66538/Guia%20Pre-Rectificacion.pdf "Guía de Pre-Rectificación")

En el respositorio del CIRN se puede encontrar una descripción más extensa de los criterios a tener en cuenta en la grabación de imágenes.

### Inputs
Información de entrada:
- Archivo de video  (Tener clara la resolución, fps y duración)
- Perfil de calibración del lente (LCP)
- Puntos GCP

### A-VideoToImagen.m
Este script [A-VideoToImagen.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/61438c70ad05e72d21f3ade8688130c404e66538/1-Rectificacion/A-VideoToImagen.m "A-VideoToImagen.m") extrae los frames de un video a una razón de fps especificada. 

Escribir el nombre del archivo `nVid` y el fps de extracción `fpsR` 

```matlab
nVid = 'DJI_0001.mov'
fpsR = 2
```

Las imágenes se irán guardando en la carpeta `Outputs\1-Frames\`.

Generalmente en estos estudios se considera una frecuencia de muestreo de 2 Hz.

El video debe estar ubicado en la carpeta actual de matlab

### B-UTMtoLocal.m
Convierte las coordenadas de los puntos GCP de un sistema global (UTM) a un sistema local definido por el usuario.

Ingrese nombre y coordenadas de los puntos con coordenadas conocidas en UTM. Ejemplo con 3 puntos (pueden ser más):

```matlab
GCP_UTM.Name(1,:) = "P1";
GCP_UTM.Name(2,:) = "P2";
GCP_UTM.Name(3,:) = "P3";
GCP_UTM.CoordUTM(1,:) = [256888.121 6289629.746]';
GCP_UTM.CoordUTM(2,:) = [256867.711 6289684.774]';
GCP_UTM.CoordUTM(3,:) = [256904.663 6289669.813]';
```
Especifique origen del sistema de referencia local en UTM (punto de referencia).
```matlab
Re = [256884.448 6289741.483]';
```
Indique el ángulo de rotación (en grados) del sistema local con respecto al sistema de coordenadas UTM.

```matlab
theta =150;
```
Los puntos quedarán guardados en la carpeta `Outputs` con el nombre `refPOINT`

### C-ConfguracionUAV.m
El script [C-ConfiguracionUAV](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/07a42b5b102765b87e33f8f1ca37075fc3c22ecf/1-Rectificacion/C-ConfiguracionUAV.m "C-ConfiguracionUAV")  permite editar los valores de configuración del algoritmo de rectificación.

[makeLCPP3.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/7f3cc3ac0d49cdd1d7c2e5dd888b9942fe669e1c/1-Rectificacion/UAV-Processing-Toolbox-master/makeLCPP3.m "makeLCPP3.m"): Crea un perfil de calibración del lente

[demoInstsFile.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/07a42b5b102765b87e33f8f1ca37075fc3c22ecf/1-Rectificacion/UAV-Processing-Toolbox-master/demoInstsFile.m "demoInstsFile.m"):  Crea instrumentos de pixeles (En este caso no se utiliza)

[demoInputFile_12V1P1](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/7f3cc3ac0d49cdd1d7c2e5dd888b9942fe669e1c/1-Rectificacion/Inputs/demoInputFile_12V1P1.m "demoInputFile_12V1P1"): Configuración valores de entrada del análsis principal.

- **makeLCPP3.m**

Modificar de acuerdo con los valores de la calibración del lente (parámetros intrínsecos). Puede modificar sólo un caso asegurandose que en el archivo de configuración [demoInputFile_12V1P1](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/7f3cc3ac0d49cdd1d7c2e5dd888b9942fe669e1c/1-Rectificacion/Inputs/demoInputFile_12V1P1.m "demoInputFile_12V1P1") se indique el caso correcto (para este caso se usa 'Aerielle').

```matlab
lcp.NU = NU;
lcp.NV = NV;
lcp.c0U = 1957.13;       
lcp.c0V = 1088.21;
lcp.fx = 2298.59;        
lcp.fy = 2310.87;
lcp.d1 = -0.14185;  % radial distortion coefficients
lcp.d2 =  0.11168;
lcp.d3 = 0.0;
lcp.t1 = 0.00369;   % tangential distortion coefficients
lcp.t2 = 0.002314;
lcp.r = 0:0.001:1.5;
lcp = makeRadDist(lcp);
lcp = makeTangDist(lcp);    % add tangential dist template
```

- **demoInputFile_12V1P1.m**

Los principales valores que se deben modificar son los siguientes:

**Paths, names and time stamp info**

| Variable | Descripción  | 
| :------------: |:------------|
| stationStr      | ''Aerielle” (Relacionado con los parámetros intrínsecos y calibración de la cámara, revisar función makeLCPP3.m )   |
| dateVect     | Fecha del primer frame  [aaaa mm dd hh mm ss]    | 
| dt | Espaciamiento temporal. Si su frecuencia de muestreo es de 2 [Hz], no cambiar        |  
| dateVect     | Fecha del primer frame  [aaaa mm dd hh mm ss]    | 

Ejemplo:

```matlab
inputs.stationStr = 'Aerielle';  
inputs.dateVect = [2018 11 12 13 08 43];       % date/time of first frame
inputs.dt = 0.5/(24*3600);           % delta t (s) converted to datenums
inputs.frameFn = 'demoClip';            % root of frame names
inputs.gcpFn = ['.\Outputs\refPOINT.mat'];   % File that contains the names and locations of all the possible GCPs 
inputs.instsFn = ['./UAV-Processing-Toolbox-master/demoInstsFile'];            % instrument m-file location
```

**Geometry solution Inputs**

Parametros extrínsecos de la cámara  $\rightarrow$  [xCam yCam zCam Azimuth Tilt Roll]

| Variable | Descripción  | 
| :------------------: |:------------|
| knownFlags  |  1  $\rightarrow$ si se conoce la variable |
|       |  0 $\rightarrow$ si no se conoce la variable |
| xyCam, zCamz azTilt, roll  |  si su valor es conocido se debe anotar, si no es conocido, este valor sirve como primer punto de iteración |

Ejemplo:

```matlab
inputs.knownFlags = [0 0 0 0 0 0];
inputs.xyCam = [ -107.4604 265.6772];
inputs.zCam = 80;             % based on last data run                
inputs.azTilt = [95 60] / 180*pi;          % first guess
inputs.roll = 0 / 180*pi; 
```
**GCP Info**

| Variable | Descripción  | 
| :------------------: |:------------|
| gcpList |  Anotar el número de los puntos GCP que utilizará en el análisis, son los que se ven en el video |
| nRefs |  Número de puntos virtuales que utilizará, son los que se distinguen fácilmente en la imagen. Pueden o no ser GCP. |

Ejemplo:

```matlab
inputs.gcpList = [1 2 3];      % use these gcps for init beta soln
inputs.nRefs = 3;                    % number of ref points for stabilization
inputs.zRefs = 2;                    % assumed z level of ref points
```

**Processing Parameters**

| Variable | Descripción  | 
| :------------------: |:------------|
| rectxy |  Grilla de rectificación en X e Y, [xmin dx xmax ymin dy ymax] |
| rectz |  Nivel vertical para la rectificación (Generalmente es el nivel medio del mar).|

Ejemplo:

```matlab
inputs.doImageProducts = 1;                    % usually 1.
inputs.showFoundRefPoints = 0;                 % to display ref points as check
inputs.showInputImages = 1;                    % display each input image as it is processed
inputs.rectxy = [-50 0.5 400 -300 0.5 300];     % rectification specs
inputs.rectz = 0;                              % rectification z-level
```
### D-RectImagenes-Part1.m

Este script [D-RectImagenes-Part1.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/c0ffd452c3347aa4881b4d49cdbdf69b210a4a16/1-Rectificacion/D-RectImagenes-Part1.m "D-RectImagenes-Part1.m") permite determinar los parámetros extrínsecos de la cámara con respecto al sistema de referencia local, para lo cual se utilizan los puntos GCP y su ubicación en la imagen.

- Es recomendable ir ejecutando la lineas de código por sección, e ir comprobando que se van obteniendo los resultados esperados.

Modifique el nombre del archivo de configuración de acuerdo con el nombre asignado:

https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/c0ffd452c3347aa4881b4d49cdbdf69b210a4a16/1-Rectificacion/D-RectImagenes-Part1.m#L14-L17

Seleccione los GCP en la imagen. Es recomendable tener una foto de referencia (como apoyo en una seguna pantalla por ejemplo) que indique claramente donde está ubicado cada punto, tal como se muestra en la siguiente imagen:

![](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/656d2c0b76168bf10eda1ac20f6b9256f4e42eed/1-Rectificacion/13v2.jpg)

Elija los gcp virtuales. Estos tienen que ser puntos idealmente blancos o muy claros, con la idea de que se diferencien bastante de su entorno, suelen servir los techos de las casas.

![](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/45b5be5e54fe93521b7cc06ec8edce67d724f02a/1-Rectificacion/Ingres_%20GCP_Virtuales.png)

**Obs:** *No es necesario que se conozcan las coordenadas reales de estos puntos. Tampoco hay problema si existen puntos que cumplen tanto como GCP como GCP virtuales.*

Una vez elegido el recuadro que encierra al GCP virtual, se debe fijar un umbral de intensidad. Debe asegurarse que este umbral sea tal, que la forma (limite) que define el GCP virtual no se altere durante el análisis de todos los frames.

**Obs:** *El valor de este umbral está entre los 180 y 220 normalamente.*

### E-RectImagenes-Part2.m

Este script [E-RectImagenes-Part2.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/26e0344d6fc7bdd50e92d18032631cd94c9b7262/1-Rectificacion/E-RectImagenes-Part2.m "E-RectImagenes-Part2.m")  rectifica todas las imágenes ubicadas en `Outputs/1-Frames` y las guarda en `Outputs/1-FramesRect` en dos formatos, como matriz de matlab .mat y comprimido en JPG .

## 2. Batimetria_cBathy

Estos scripts corresponden a una aplicación o adaptación del algoritmo [cBathy-Toolbox](https://github.com/Coastal-Imaging-Research-Network/cBathy-Toolbox "cBathy-Toolbox")  del CIRN para obtener la batimetría del lugar estudiado.

### A_CrearInstrumentoPixel.m
Este script genera un instrumento pixel que consiste en una región rectangular discretizada  (con una resolución especfificada por el usuario) sobre la cual se van almacenando los pixeles rectificados y georeferenciados.

Especifique resolución espacial de la grilla (se considera la misma resolución en las dos direcciones x e y) y el nombre con el que se guardará el instrumento:

```matlab
res = 5          %Resolucion [m]
oname = ['15V3_' num2str(res) 'm'];
```
Luego escriba la direccion (carpeta) donde se encuentren las imágenes rectificadas. Deben ser en formato de matriz de matlab (.mat), en el caso de ser un formato distinto, realizar las modificaciones correspondientes.

Además, ingresar la ubicación del archivo `DataFrames` generado en el proceso de rectificación de las imágenes.

```matlab
%Imagenes Rectificadas
imageDirectory{1} =  'C:\Memoria2020\Resultados\Rectificados\15V3\MAT';

%Info Imagenes Rectificadas

DireccionDF = './Inputs/DataFrames_15V3';
```
Ingrese los límites de la grilla `xlim` y `ylim`

```matlab
pixInst(1).type = 'Grid';
pixInst(1). dx  = res;
pixInst(1). dy =res;
pixInst(1). xlim = [80 400];
pixInst(1). ylim = [-300 300];
pixInst(1).z = {};
```
El instrumento quedará guardado en la carpeta `/Outputs/1-PixelInstruments/`

### B_CreacionArchivoCbathy.m
Este script compatibiliza el formato y nomenclaturas de archivos para utilizar cBathy-Toolbox.

Sólo debe ingresar la dirección de los archivos `DataFrames` y el `Instrumento Pixel` generado anteriormente.

```matlab
load('./Outputs/1-PixelInstruments/15V3_5m_pixInst.mat')
load('./Inputs/DataFrames_15V3.mat')
```
**Obs:** Se considera por defecto un $\Delta$ t de 0.5 segundos.

En la sección 3, especificar nombre del archivo generado: 
```matlab
% Output Name
oname=['15V3_Archivo_cBathy-5m'];

% OutPut Directory
odir=['./Outputs/2-Archivos-cBathy/'];
```
### C_EstimarBathy.m
Este script ejecuta la caja de códigos cBathy-Toolbox para estimar batimetría.

Primero seleccione el archivo de entrada:

```matlab
load('15V3_Archivo_cBathy-5m.mat')
```
Luego edite el archivo de configuración `argus02a`:

```matlab
edit argus02a
```

En la siguiente tabla se muestra una breve descripción de los principales parámetros. Para mayor detalles consultar documentación en el repositorio de cBathy-Toolbox.

| Parámetro | Descripción  | 
| :------------------: |:------------|
| dxm y dym |  Espaciamiento deseado para los puntos de análisis, espaciamiento de salida |
| xyMinMax|  Mínimo y máximo de x e y|
| MINDEPTH y MAXDEPTH|  Límites para la profundidad|
| Lx y Ly|  Parámetros que definen la vecindad que se analiza por cada punto|
| fB|  Rango de frecuencias de análisis.|
| nKeep | N° de frecuencias que se consideran en el análisis (las que poseen mayor coherencia dentro de fB)|

Ejemplo:

```matlab
%%% Site-specific Inputs

params.stationStr = 'argus02a';
                
params.dxm = 10;  %4 %5                % analysis domain spacing in x
params.dym = 10;  %8 %10               % analysis domain spacing in y
params.xyMinMax = [80 400 -300 300];   % min, max of x, then y
                                      % default to [] for cBathy to choose
params.tideFunction = 'cBathyTide';   % tide level function for evel

%%%%%%%   Power user settings from here down   %%%%%%%
params.MINDEPTH = 0.25;             % for initialization and final QC
params.minValsForBathyEst = 4;      % min num f-k pairs for bathy est.

params.QTOL = 0.5;                  % reject skill below this in csm
params.minLam = 10;                 % min normalized eigenvalue to proceed
params.Lx = 3*params.dxm;           % tomographic domain smoothing
params.Ly = 3*params.dym;           % 
params.kappa0 = 2;                  % increase in smoothing at outer xm
params.DECIMATE = 1;                % decimate pixels to reduce work load.
params.maxNPix = 80;                % max num pixels per tile (decimate excess)

% f-domain etc.
params.fB = [1/18: 1/50: 1/4];		% frequencies for analysis (~40 dof)
params.nKeep = 4;                   % number of frequencies to keep
```

Ejecute todas las secciones y al finalizar especificar nombre del archivo de salida.

### D_KalmanFiltering.m
Este script permite suavizar y rellenar las estimaciones de batimetría a partir de estimaciones anteriores, utilizando un filtro de Kalman. 

Se aconseja revisar documentación de cBathy-Toolbox para una correcta ejecución.

## 3. Corrientes
Estos scripts permiten obtener las corrientes medias longitudinales y transversales en la zona de la rompiente.

Para la obtención de las corrientes longitudinales se utiliza el código creado por Chris Chickadel denominado [Video-Currents-Toolbox](https://github.com/Coastal-Imaging-Research-Network/Video-Currents-Toolbox "Video-Currents-Toolbox") disponible en el repositorio del CIRN.

Para la obtención de las velocidades en dirección transversal a la línea de costa, se aplica una continuidad de flujos, lo que implica necesariamente disponer de la batimetría (es un valor de entrada necesario).

Los parámetros principales de los algoritmos implementados se deben modificar en el script [Configuracion.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/c44be5d18da9d52d4940f0f3b167f289088472f4/3-Corrientes/Configuracion.m "Configuracion.m").

Ejemplo:

```matlab
NameV = '15V3'; %Vuelo a analizar

%Parámetros
dy = 10 ;  %[m] Longitud de los instrumentos vBar
dx = 10 ;  %[m] Separación instrumentos en dirección x
Tw = 32 ;  %[s] Duración ventana de análisis
Ts = 16 ;  %[s] Tiempo de muestreo
km = 0.25; %[1/m] Número de onda mínimo
P  = 0.9 ; %Límite para el nivel de confianza de ajuste al modelo
Lci= 0.2 ; %[m/s] Umbral para el intervalo de confianza de la velocidad media
Ir = 40  ; %Umbral para el rango de intensidad

%Límites grilla
Xlim=[70,370];
Ylim=[-270,260];
```

Además, para la creación de instrumentos a partir de imágenes rectificadas, debe especificar la dirección donde se encuentran las imágenes en formato .mat, y la ubicación del archivo `DataFrames` con la información de las imágenes en el script  [A_CreacionInstrumento.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/5fbdf9764821787608e8264bf845b4c93a0b4688/3-Corrientes/Toolbox-Corrientes/A_CreacionInstrumento.m "A_CreacionInstrumento.m"). Ejemplo:

https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/5fbdf9764821787608e8264bf845b4c93a0b4688/3-Corrientes/Toolbox-Corrientes/A_CreacionInstrumento.m#L11-L16

Posteriormente debe ejecutar el script [A_DemoCorrientes.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/c44be5d18da9d52d4940f0f3b167f289088472f4/3-Corrientes/A_DemoCorrientes.m "A_DemoCorrientes.m")

Todos los resultados quedarán guardados en la carpeta `Outputs`.

## 4. Documentacion
La base teórica de los algoritmos, así como su aplicación a una playa ubicada en la zona central de Chile (Playa Las Cruces) se encuentra explicada en mi memoria:

- [Determinación de campos de velocidades en la zona de la rompiente desde drones - Esteban Opazo Verdugo](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/852a275e61784a37778b0f2de36781110a8ebe9a/Determinacion%20de%20Velocidades%20-%20Memoria%20Esteban%20Opazo.pdf "Determinación de campos de velocidades en la zona de la rompiente desde drones - Esteban Opazo Verdugo")

Además comparto unos apuntes que escribí con una explicación más detallada de los algoritmos de rectificación, incluyendo cbathy y la obtención de velocidades:

- [Apuntes Rectificación de Imágenes](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/852a275e61784a37778b0f2de36781110a8ebe9a/Apuntes%20Rectificacion%20de%20Imagenes%20-%20EOV.pdf "Apuntes Rectificación de Imágenes")

Otros:

- [Guía rápida - Rectificación de imágenes](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/534b992e53cd3672b9da5c059dffd823398dc7c1/Guia%20Rapida-Rectificacion%20de%20Imagenes.pdf "Guía rápida - Rectificación de imágenes") (Versión desactualizada)

## 5. Contacto

Mail: esteban.opazo.13@sansano.usm.cl
