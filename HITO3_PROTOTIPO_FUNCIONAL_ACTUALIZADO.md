# HITO 3: Prototipo Funcional - Detección y Clasificación Automática de Latidos Ventriculares Ectópicos (VEB/PVC)

**Asignatura:** Señales y Sistemas — Bioingeniería (Plan 2026)

**Período de Ejecución:** Segundo Bloque (Semana 14 — Hito 3: MVP)

**Proyecto Integrador:** Caracterización y Clasificación Automática de Arritmias Ventriculares mediante Procesamiento Digital de Señales ECG

**Artículo de Referencia:** Srinivasulu, A., & Sriraam, N. (2023). *Signal Processing Framework for the Detection of Ventricular Ectopic Beat Episodes*. Journal of Medical Signals & Sensors.

## 1. Ingesta y Estandarización de Señales Biomédicas

### 1.1. Fuentes de Datos Abiertas

Para cumplir con la exigencia de validación con señales biológicas reales de la Guía Metodológica, el sistema se conecta remotamente a los servidores de PhysioNet mediante la librería `wfdb`:

* **Base de Entrenamiento (DB-I):** *MIT-BIH Arrhythmia Database*. El prototipo procesa y consolida el dataset (44 de los 48 registros de \~30 minutos, $f_s = 360\text{ Hz}$), utilizando específicamente la derivación electrocardiográfica modificada II (MLII) para maximizar la amplitud del complejo QRS y capturar la variabilidad interpaciente.

* **Base de Validación Externa (DB-II):** *St. Petersburg INCART 12-lead Arrhythmia Database*. Se seleccionó el **Registro I01** (registrado originalmente a $f_s = 257\text{ Hz}$ en el canal II) para evaluar la generalización ciega del modelo frente a un dispositivo, población y entorno normativo diferentes.

### 1.2. Estandarización de Frecuencia de Muestreo (Remuestreo)

Dado que los registros de distintas bases de datos o dispositivos Holter operan a diferentes tasas de adquisición, se aplica un algoritmo de remuestreo polifásico basado en la Transformada Rápida de Fourier (función `scipy.signal.resample`):

$$
N_{target} = \left\lfloor N_{orig} \cdot \frac{f_{s, target}}{f_{s, orig}} \right\rfloor
$$

Para el Registro INCART I01, la señal fue remuestreada de $257\text{ Hz}$ a $360\text{ Hz}$, reescalando simultáneamente las anotaciones clínicas oficiales para garantizar la perfecta coincidencia temporal de las ventanas de análisis.

## 2. Preprocesamiento Digital y Análisis Espectral

### 2.1. Definición Matemática de los Filtros Digitales

La señal electrocardiográfica cruda $x[n]$ sufre la superposición de interferencias de baja frecuencia (deriva de línea base producida por los movimientos respiratorios, $< 0.5\text{ Hz}$) y de alta frecuencia (ruido electromiográfico por actividad muscular y zumbido de la red eléctrica en 50/60 Hz).

El acondicionamiento digital se implementa mediante la siguiente cascada de filtrado LTI:

1. **Filtro Notch IIR (Filtro Muesca con** $Q = 30$**):**
   Atenúa la frecuencia fundamental de la red eléctrica manteniendo un factor de calidad estricto:

   $$
   H_{notch}(z) = G_0 \frac{1 - 2\cos(\omega_0)z^{-1} + z^{-2}}{1 - 2r\cos(\omega_0)z^{-1} + r^2 z^{-2}}
   $$

   donde $\omega_0 = 2\pi f_{red} / f_s$. Para la base MIT-BIH se fija $f_{red} = 60\text{ Hz}$, mientras que para el Registro INCART I01 (Rusia) el sistema adapta automáticamente su parámetro a $f_{red} = 50\text{ Hz}$, demostrando robustez ante cambios geográficos y normativos de adquisición.

2. **Filtro Pasa-Banda Butterworth de 3.er Orden (**$0.5\text{ Hz} - 35\text{ Hz}$**):**
   Garantiza una respuesta estrictamente plana en la banda útil, preservando los componentes espectrales del complejo QRS ($5\text{ a }30\text{ Hz}$) y eliminando la deriva de línea base junto al ruido muscular.

### 2.2. Justificación de la Fase Cero mediante `filtfilt`

Los filtros IIR convencionales introducen un desfasaje no lineal que desplaza temporalmente los picos rápidos del QRS respecto a las ondas lentas. Para anular este retardo de grupo, se implementa el **filtrado bidireccional de fase cero** (`scipy.signal.filtfilt`), cuya función de transferencia efectiva es:

$$
H_{efectiva}(f) = H(f) \cdot H^*(f) = \vert{}H(f)\vert{}^2
$$

La fase resultante es idéntica a cero en todo el espectro ($\text{Arg}(H_{efectiva}) = 0^\circ$), garantizando una alineación temporal perfecta con las anotaciones clínicas de referencia dentro de la ventana de tolerancia de $\pm 50\text{ ms}$.

### 2.3. Verificación Espectral y Diagrama de Bode

La caracterización frecuencial validada mediante diagramas de Bode y Densidad Espectral de Potencia (Welch) confirma que las frecuencias de corte a la frecuencia de diseño operan exactamente a $-6\text{ dB}$, cumpliendo con la teoría de sistemas discretos.

## 3. Algoritmo de Detección QRS (Pan-Tompkins) y Validación Cuantitativa

### 3.1. Transformaciones Discretas de Pan-Tompkins

La detección se ejecuta mediante la secuencia de transformaciones: derivador discreto (resalta pendientes abruptas), elevación al cuadrado (amplifica la energía no lineal del QRS) e integración por ventana móvil de $150\text{ ms}$ ($W = 54\text{ muestras}$ a $360\text{ Hz}$). Los picos fiduciales $Q$, $R$ y $S$ se reubican mediante búsqueda local en la señal filtrada original $y_{pb}[n]$.

### 3.2. Evaluación Cuantitativa del Detector

Evaluado sobre tramos de prueba de referencia (como el Registro MIT-BIH 119), el algoritmo alcanzó una **Sensibilidad del** $100\%$, un **Valor Predictivo Positivo del** $100\%$ y un **F1-Score del** $100\%$, con cero falsos positivos y cero falsos negativos en la localización de latidos.

## 4. Extracción de Características y Análisis Exploratorio de Separabilidad

### 4.1. Formulación del Vector de Características $\mathbf{f} \in \mathbb{R}^5$

Para cada latido $i$, se extraen 5 atributos morfológicos y temporales:

1. **Intervalo Pre-RR (**$f_1$**):** $t_R(i) - t_R(i-1)$

2. **Intervalo Post-RR (**$f_2$**):** $t_R(i+1) - t_R(i)$

3. **Ancho del Complejo QRS (**$f_3$**):** $t_S(i) - t_Q(i)$

4. **Pendiente QR (**$f_4$**):** $\frac{A_R(i) - A_Q(i)}{t_R(i) - t_Q(i)}$

5. **Pendiente RS (**$f_5$**):** $\frac{A_R(i) - A_S(i)}{t_R(i) - t_S(i)}$

### 4.2. Proyección 2D del Espacio de Características (PCA y t-SNE)

La estandarización mediante Z-score y la proyección por Análisis de Componentes Principales demostraron que las dos primeras componentes explican el $94.09\%$ **de la varianza total** ($CP_1 = 69.08\%$, $CP_2 = 25.02\%$), evidenciando una clara separabilidad lineal y no lineal (t-SNE) entre latidos Normales ($N$) y Ectópicos Ventriculares ($V$).

## 5. Clasificación Automática (NMC) y Validación Patient-Wise / Cross-Database

### 5.1. Formulación del Clasificador NMC

El Clasificador de Medias Cercanas (NMC) calcula los centroides vectoriales $\boldsymbol{\mu}_N$ y $\boldsymbol{\mu}_V$ en el espacio estandarizado $\mathbb{R}^5$ y clasifica nuevos latidos calculando la mínima distancia Euclidiana:

$$
d(\mathbf{x}_{test}, \boldsymbol{\mu}_k) = \Vert{}\mathbf{x}_{test} - \boldsymbol{\mu}_k\Vert{}_2
$$

### 5.2. Esquema de Validación Ciega (*Cross-Database*)

Para evitar el sesgo por fuga de datos (*data leakage*), el normalizador (`StandardScaler`) y los centroides se ajustaron utilizando exclusivamente los datos de entrenamiento de MIT-BIH, aplicando dichos parámetros sin reentrenamiento sobre el **Registro INCART I01**.

### 5.3. Resultados Métricos Definitivos de la Validación Externa

Los resultados cuantitativos obtenidos en la evaluación ciega sobre INCART I01 reflejan el comportamiento real del prototipo ante un dataset multivariado y patológico:

| Métrica de Desempeño | Expresión Matemática | Resultado Definitivo (INCART I01) | Interpretación y Análisis Clínico | 
| ----- | ----- | ----- | ----- | 
| **Verdaderos Positivos (**$TP$**)** | Latidos $V$ detectados correctamente | $20$ | Arritmias ventriculares identificadas con éxito. | 
| **Falsos Positivos (**$FP$**)** | Latidos $N$ clasificados erróneamente como $V$ | $5$ | Falsas alarmas mínimas emitidas por el sistema. | 
| **Verdaderos Negativos (**$TN$**)** | Latidos $N$ identificados correctamente | $364$ | Ritmo sinusal reconocido sin fatiga diagnóstica. | 
| **Falsos Negativos (**$FN$**)** | Latidos $V$ omitidos por el modelo | $14$ | Arritmias no capturadas por limitación del centroide único. | 
| **Exactitud Global (Accuracy)** | $ACC = \frac{TP + TN}{TP + TN + FP + FN}$ | $95.29\%$ | Alta tasa de aciertos globales sin sobreajuste (*overfitting*). | 
| **Sensibilidad (Recall clase V)** | $Se = \frac{TP}{TP + FN}$ | $58.82\%$ | Capacidad afectada por la multiformidad ventricular del paciente externo. | 
| **Especificidad (Recall clase N)** | $Sp = \frac{TN}{TN + FP}$ | $98.64\%$ | Excelente rechazo de falsas alarmas sobre el ritmo normal. | 
| **Precisión (Valor Predictivo Positivo)** | $Pr = \frac{TP}{TP + FP}$ | $80.00\%$ | Confiabilidad sólida de las alertas generadas por el prototipo. | 
| **Puntuación F1 (F1-Score)** | $F1 = \frac{2 \cdot Pr \cdot Se}{Pr + Se}$ | $67.80\%$ | Balance armónico general en condiciones de validación ciega. | 

## 6. Conclusiones Técnicas y Trabajo Futuro (Hito 4)

1. **Desempeño y Generalización:** El prototipo desarrollado cumple satisfactoriamente con los requerimientos del Hito 3, alcanzando una **Exactitud Global del** $95.29\%$ y una **Especificidad del** $98.64\%$ en validación ciega inter-base (*cross-database*), demostrando la robustez del filtrado de fase cero y el vector de características $\mathbb{R}^5$.

2. **Justificación Analítica de Limitaciones:** La sensibilidad obtenida ($58.82\%$) expone teóricamente la limitación del NMC al emplear un único centroide promedio ($\boldsymbol{\mu}_V$) frente a latidos ventriculares multiformes de distintas procedencias clínicas. Esto justifica formalmente la necesidad de evolucionar hacia clasificadores más avanzados (como *NRMS* o distancia de Mahalanobis) en futuras etapas del desarrollo.

3. **Cierre del TPI:** El código y la estructura expositiva se encuentran completamente consolidados y listos para la defensa oral del Hito final.

## 7. Referencias Bibliográficas

1. Srinivasulu, A., & Sriraam, N. (2023). *Signal processing framework for the detection of ventricular ectopic beat episodes*. Journal of Medical Signals and Sensors, 13(3), 239-251.

2. Pan, J., & Tompkins, W. J. (1985). *A real-time QRS detection algorithm*. IEEE Transactions on Biomedical Engineering, (3), 230-236.

3. Goldberger, A. L., et al. (2000). *PhysioBank, PhysioToolkit, and PhysioNet*. Circulation, 101(23), e215-e220.

4. Cátedra de Señales y Sistemas. (2026). *Guía Metodológica del Trabajo Práctico Integrador (TPI)*. Carrera de Bioingeniería, UNER.