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
