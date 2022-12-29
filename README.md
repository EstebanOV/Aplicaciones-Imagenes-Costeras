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

#### Inputs
Información de entrada:
- Archivo de video  (Tener clara la resolución, fps y duración)
- Perfil de calibración del lente (LCP)
- Puntos GCP

#### A-VideoToImagen.m
Este script [A-VideoToImagen.m](https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/61438c70ad05e72d21f3ade8688130c404e66538/1-Rectificacion/A-VideoToImagen.m "A-VideoToImagen.m") extrae los frames de un video a una razón de fps especificada. 

Escribir el nombre del archivo (nVid) y el fps de extracción (fpsR) 

https://github.com/EstebanOV/Aplicaciones-Imagenes-Costeras/blob/61438c70ad05e72d21f3ade8688130c404e66538/1-Rectificacion/A-VideoToImagen.m#L11-L21

Las imágenes se irán guardando en la carpeta ** *Outputs\1-Frames\* **.

Obs: generalmente en estos estudios se considera una frecuencia de muestreo de 2 Hz.
