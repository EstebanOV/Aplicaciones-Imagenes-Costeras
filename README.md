# Aplicaciones Imágenes Costeras 
Este repositorio contiene los códigos utilizados para la obtención del campo de velocidades medias en la zona de la rompiente a partir de imágenes grabadas desde un dron. 

Estos códigos fueron generados como parte de mi trabajo de memoria para optar el título de Ingeniero Civil, y se enmarca dentro del proyecto *"FONDECYT 1170415-Quantification of Two Dimensional Wave Breaking Dissipation in the Surf Zone"* liderado por el profesor Patricio Catalán del Departamento de Ingeniería Civil de la Universidad Técnica Federico Santa María.


El repositorio se encuentra organizado en 3 grupos: 

- **Rectificación:** Códigos utilizados para la rectificación de imágenes. 

- **Batimetría:**  Códigos de cBathy adaptados al proyecto para la obtención de batimetría en la zona.

- **Corrientes:** Códigos generados para la obtención de corrientes medias longitudinales y transversales. 

Cada uno de estos grupos de códigos se puede utilizar de forma independiente y adaptar a las condiciones requeridas.

## 1- Rectificación
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

Paths, names and time stamp info

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

Geometry solution Inputs
