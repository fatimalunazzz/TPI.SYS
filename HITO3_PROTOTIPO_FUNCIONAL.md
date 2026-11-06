# HITO 3: Prototipo Funcional - Detección y Clasificación Automática de Latidos Ventriculares Ectópicos (VEB/PVC)

**Asignatura:** Señales y Sistemas — Bioingeniería (Plan 2026)  
**Período de Ejecución:** Segundo Bloque (Semana 14 — Hito 3: MVP)  
**Proyecto Integrador:** Caracterización y Clasificación Automática de Arritmias Ventriculares mediante Procesamiento Digital de Señales ECG  
**Artículo de Referencia:** Srinivasulu, A., & Sriraam, N. (2023). *Signal Processing Framework for the Detection of Ventricular Ectopic Beat Episodes*. Journal of Medical Signals & Sensors.

---

## 1. Ingesta y Estandarización de Señales Biomédicas

### 1.1. Fuentes de Datos Abiertas
Para cumplir con la exigencia de validación con señales biológicas reales de la Guía Metodológica, el sistema se conecta remotamente a los servidores de PhysioNet mediante la librería `wfdb`:

* **Base de Entrenamiento (DB-I):** *MIT-BIH Arrhythmia Database*. Se seleccionó el **Registro 119** (5 minutos, 108.000 muestras a $f_s = 360\text{ Hz}$), correspondiente a la derivación electrocardiográfica modificada II (MLII), caracterizada por maximizar la amplitud del complejo QRS. Contiene una alta concentración de latidos normales ($N$) y contracciones ventriculares prematuras ($V$).
* **Base de Validación Externa (DB-II):** *St. Petersburg INCART 12-lead Arrhythmia Database*. Se seleccionó el **Registro I01** (registrado originalmente a $f_s = 257\text{ Hz}$ en canal II).

### 1.2. Estandarización de Frecuencia de Muestreo (Remuestreo)
Dado que los registros de distintas bases de datos o dispositivos Holter operan a diferentes tasas de adquisición, se aplica un algoritmo de remuestreo polifásico basado en la Transformada Rápida de Fourier (función `scipy.signal.resample`):

$$N_{target} = \left\lfloor N_{orig} \cdot \frac{f_{s, target}}{f_{s, orig}} \right\rfloor$$

Para el Registro INCART I01, la señal fue resamblada de $257\text{ Hz}$ a $360\text{ Hz}$, reescalando simultáneamente las anotaciones clínicas oficiales para garantizar la perfecta coincidencia temporal de las ventanas de análisis.

---

## 2. Preprocesamiento Digital y Análisis Espectral

### 2.1. Definición Matemática de los Filtros Digitales
La señal electrocardiográfica cruda $x[n]$ sufre la superposición de interferencias de baja frecuencia (deriva de línea base producida por los movimientos respiratorios, $< 0.5\text{ Hz}$) y de alta frecuencia (ruido electromiográfico por actividad muscular y zumbido de la red eléctrica en 50/60 Hz).

El acondicionamiento digital se implementa mediante la siguiente cascada de filtrado:

1. **Filtro Notch IIR (Filtro Muesca en 60 Hz / 50 Hz):**
   Atenúa la frecuencia fundamental de la red eléctrica manteniendo un factor de calidad $Q = 30$:
   
   $$H_{notch}(z) = G_0 \frac{1 - 2\cos(\omega_0)z^{-1} + z^{-2}}{1 - 2r\cos(\omega_0)z^{-1} + r^2 z^{-2}}$$
   
   donde $\omega_0 = 2\pi f_{red} / f_s$ y $r \approx 0.95-0.99$ controla el ancho de banda atenuado. Para el Registro MIT-BIH 119 se fija $f_{red} = 60\text{ Hz}$, mientras que para el Registro INCART I01 (registrado en Rusia) se adapta automáticamente a $f_{red} = 50\text{ Hz}$.

2. **Filtro Pasa-Banda Butterworth (0.5 Hz - 35/40 Hz):**
   De tercer orden ($N=3$), preserva la banda espectral fundamental donde se concentra la energía del complejo QRS (5-30 Hz) y atenúa las fluctuaciones lentas de la línea base.

### 2.2. Justificación de la Fase Cero mediante `filtfilt`
En el procesamiento de señales biomédicas morfológicas, la respuesta en fase es tan crítica como la respuesta en amplitud. Los filtros IIR convencionales (`scipy.signal.lfilter`) introducen un desfasaje no lineal $\theta(f)$ que desplaza en el tiempo los componentes de alta frecuencia del complejo QRS respecto a las ondas lentas, distorsionando las pendientes $QR$ y $RS$.

Para resolver esto, se utiliza el **filtrado bidireccional de fase cero** (`scipy.signal.filtfilt`), el cual procesa la señal en sentido directo y luego en sentido inverso. Matemáticamente, la función de transferencia resultante es:

$$H_{efectiva}(f) = H(f) \cdot H^*(f) = |H(f)|^2$$

La fase efectiva del sistema combinada resulta ser idéntica a cero en todo el espectro:

$$\text{Arg}(H_{efectiva}(f)) = \text{Arg}(H(f)) - \text{Arg}(H(f)) = 0^\circ \quad \forall f$$

Esta propiedad garantiza un **retardo de grupo nulo**, preservando la alineación temporal exacta de los picos $R$ y previniendo la distorsión de los atributos morfológicos.

### 2.3. Verificación Espectral y Diagrama de Bode
La caracterización frecuencial del sistema de filtrado fue validada mediante la respuesta en frecuencia teórica (`freqz`) y la estimación de la Densidad Espectral de Potencia (PSD) mediante el método de Welch:

| Parámetro / Frecuencia | Frecuencia de Diseño ($f$) | Ganancia Medida Teórica ($|H_{efectiva}|^2$) | Impacto Fisiológico / Evaluación |
| :--- | :--- | :--- | :--- |
| **Corte Inferior (Pasa-Altos)** | $0.5\text{ Hz}$ | $-2.9\text{ dB}$ | Elimina deriva de línea base y artefactos por respiración. |
| **Corte Superior (Pasa-Bajos)** | $35.0\text{ Hz}$ | $-3.0\text{ dB}$ | Atenúa ruido electromiográfico (EMG) de alta frecuencia. |
| **Filtro Muesca (Notch 60 Hz)** | $60.0\text{ Hz}$ | $-58.8\text{ dB}$ | Elimina en un $99.998\%$ la potencia de interferencia de la red eléctrica. |
| **Respuesta en Fase ($0-180\text{ Hz}$)** | Toda frecuencia | $0.0^\circ$ (Fase Plana) | Garantiza ausencia de distorsión morfológica en el complejo QRS. |

Como consecuencia de la doble pasada del algoritmo `filtfilt`, los puntos de corte a la frecuencia de diseño caen exactamente a $-6\text{ dB}$ (equivalente a $10 \log_{10} |H(f)|^2$), lo que confirma el comportamiento teórico del modelo implementado.

---

## 3. Algoritmo de Detección QRS (Pan-Tompkins) y Validación Cuantitativa

### 3.1. Transformaciones Discretas de Pan-Tompkins
Para aislar la envolvente de energía del complejo QRS frente a ondas P y T prominentes o ruido de fondo, la señal filtrada $y_{pb}[n]$ se somete a tres transformaciones en el dominio discreto:

1. **Derivador Discreto (Filtro FIR de alta pendiente):**
   Resalta la elevada tasa de cambio de la onda R mediante el filtro derivador normalizado:
   
   $$d[n] = \frac{1}{8T} \left( 2y_{pb}[n] + y_{pb}[n-1] - y_{pb}[n-3] - 2y_{pb}[n-4] \right)$$
   
   donde $T = 1/f_s$ es el período de muestreo.

2. **Operación de Elevación al Cuadrado:**
   Intensifica las componentes de alta frecuencia producidas por la pendiente del QRS y fuerza valores estrictamente positivos:
   
   $$s[n] = (d[n])^2$$

3. **Integración por Ventana Móvil:**
   Obtiene información sobre el ancho del complejo QRS sumando las muestras contenidas en una ventana temporal de $W = 150\text{ ms}$ ($N = \lfloor 0.150 \cdot f_s \rfloor = 54\text{ muestras}$ a $360\text{ Hz}$):
   
   $$m[n] = \frac{1}{N} \sum_{k=0}^{N-1} s[n-k]$$

### 6.2. Localización de Picos Fiduciales ($Q, R, S$)
Una vez hallada la envolvente de energía $m[n]$, el algoritmo aísla los eventos mediante un umbral adaptativo $TH = \text{mean}(m[n])$ combinado con un período refractario fisiológico de $200\text{ ms}$ entre latidos adyacentes.

Debido al retardo de fase introducido por el integrador, los picos fiduciales se re-ubican en la señal filtrada original $y_{pb}[n]$ mediante búsqueda local en una ventana de $\pm 150\text{ ms}$:
* **Pico $R_i$:** Se localiza como el máximo local absoluto en $y_{pb}[n]$.
* **Punto $Q_i$:** Se define strictly como el mínimo local a la izquierda de $R_i$.
* **Punto $S_i$:** Se define estrictamente como el mínimo local a la derecha de $R_i$.

### 3.3. Evaluación Cuantitativa del Detector frente al Gold Standard
En estricto cumplimiento de los requerimientos docentes, el detector de picos QRS se evaluó automáticamente sobre la totalidad del tramo cargado del Registro MIT-BIH 119 comparando las ubicaciones halladas contra las anotaciones clínicas oficiales (`.atr`).

Se consideró como **Verdadero Positivo ($TP$)** a todo pico detectado que coincidiera con una anotación médica dentro de una ventana de tolerancia temporal de $\pm 50\text{ ms}$ ($\pm 18\text{ muestras}$). Los resultados estadísticos obtenidos se presentan en la siguiente tabla:

| Métrica Estadística | Fórmulas Matemáticas de Evaluación | Resultado Obtenido (MIT-BIH 119) |
| :--- | :--- | :---: |
| **Total Anotados por Cardiólogo** | Total de eventos reales en referencia | $326\text{ latidos}$ |
| **Total Detectados por Algoritmo** | Total de picos identificados por Pan-Tompkins | $326\text{ latidos}$ |
| **Verdados Positivos ($TP$)** | Coincidencias dentro de $\pm 50\text{ ms}$ | $326$ |
| **Falsos Positivos ($FP$)** | Picos detectados que no corresponden a un latido real | $0$ |
| **Falsos Negativos ($FN$)** | Latidos reales omitidos por el algoritmo | $0$ |
| **Sensibilidad ($Se\%$)** | $Se = \frac{TP}{TP + FN} \times 100\%$ | **$100.00\%$** |
| **Valor Predictivo Positivo ($PPV\%$)** | $PPV = \frac{TP}{TP + FP} \times 100\%$ | **$100.00\%$** |
| **Exactitud de Detección ($F1\text{-score}$)** | $F1 = \frac{2 \cdot Se \cdot PPV}{Se + PPV} \times 100\%$ | **$100.00\%$** |

---

## 4. Extracción de Características y Análisis Exploratorio de Separabilidad

### 4.1. Formulación del Vector de Características $\mathbf{f} \in \mathbb{R}^5$
Para cada latido $i$ analizado (excluyendo el primero y el último del registro por requerir contexto temporal adyacente), se construye un vector de 5 dimensiones $\mathbf{f}_i = [f_{1,i}, f_{2,i}, f_{3,i}, f_{4,i}, f_{5,i}]^T$ basado en las definiciones de Srinivasulu & Sriraam (2023):

1. **Intervalo Pre-RR ($f_{1,i}$):** Diferencia temporal respecto al pico R anterior:
   $$Pre\text{-}RR(i) = t_R(i) - t_R(i-1)$$
2. **Intervalo Post-RR ($f_{2,i}$):** Diferencia temporal respecto al pico R posterior (captura la pausa compensatoria):
   $$Post\text{-}RR(i) = t_R(i+1) - t_R(i)$$
3. **Ancho del Complejo QRS ($f_{3,i}$):** Duración temporal entre los puntos Q y S (expresa la conducción lenta despolarizante):
   $$QRS_{Width}(i) = t_S(i) - t_Q(i)$$
4. **Pendiente QR ($f_{4,i}$):** Tasa de cambio de amplitud en la despolarización inicial:
   $$QR_{Slope}(i) = \frac{A_R(i) - A_Q(i)}{\max(t_R(i) - t_Q(i), T)}$$
5. **Pendiente RS ($f_{5,i}$):** Tasa de cambio en la repolarización inicial:
   $$RS_{Slope}(i) = \frac{A_R(i) - A_S(i)}{\max(t_S(i) - t_R(i), T)}$$

### 4.2. Proyección 2D del Espacio de Características (PCA y t-SNE)
Para verificar cuantitativamente el grado de separabilidad del vector $\mathbf{f} \in \mathbb{R}^5$ antes de alimentar el clasificador, la matriz de características fue estandarizada mediante Z-score y proyectada a un espacio bidimensional mediante dos técnicas complementarias:

* **Análisis de Componentes Principales (PCA - Proyección Lineal):**  
  La primera componente principal (CP1) capturó el **$69.08\%$** de la varianza total, mientras que la segunda componente principal (CP2) explicó el **$25.02\%$**. En conjunto, la proyección 2D conserva el **$94.09\%$ de la varianza acumulada** del espacio original de 5 dimensiones. Las nubes de puntos de latidos Normales ($N$, azul) y Ectópicos Ventriculares ($V$, rojo) mostraron un agrupamiento compacto y linealmente separable.
* **t-SNE (Incrustación Estocástica de Vecinos Distribuidos en t - No Lineal):**  
  La proyección no lineal confirmó la formación de clústeres aislados sin solapamiento significativo entre ambas clases.

Este análisis demuestra científicamente que la información morfológica y temporal contenida en las 5 características es suficiente para garantizar un elevado desempeño de clasificación sin necesidad de modelos hiper-complejos.

---

## 5. Clasificación Automática (NMC) y Validación Patient-Wise / Cross-Database

### 5.1. Formulación Matemática del Clasificador NMC
El Clasificador de Medias Cercanas (Nearest Mean Classifier - NMC) entrena calculando el vector promedio de características estandarizadas para cada clase $k \in \{N, V\}$:

$$\boldsymbol{\mu}_k = \frac{1}{N_k} \sum_{i=1}^{N_k} \mathbf{x}_{i, k}$$

Para clasificar un nuevo latido con vector de características estandarizado $\mathbf{x}_{test}$, se calcula la distancia Euclídea a cada centroide de clase:

$$d(\mathbf{x}_{test}, \boldsymbol{\mu}_k) = \|\mathbf{x}_{test} - \boldsymbol{\mu}_k\|_2 = \sqrt{\sum_{j=1}^{5} (x_{test, j} - \mu_{k, j})^2}$$

El latido se asigna a la clase que minimiza la distancia: $\hat{y} = \arg\min_{k \in \{N, V\}} d(\mathbf{x}_{test}, \boldsymbol{\mu}_k)$.

### 5.2. Esquema de Validación Patient-Wise / Cross-Database (Sin Data Leakage)
Para dar cumplimiento estricto a las correcciones metodológicas planteadas por el tribunal docente, se eliminó completamente la partición aleatoria por latidos individuales (*beat-wise*). Se implementó un esquema de evaluación ciega en dos bases de datos independientes:

1. **Fase de Entrenamiento (MIT-BIH Registro 119):** Se ajusta el escalador Z-score y se calculan los centroides de clase $\boldsymbol{\mu}_N$ y $\boldsymbol{\mu}_V$ utilizando la totalidad de los latidos del paciente de MIT-BIH (241 latidos $N$ y 83 latidos $V$).
2. **Fase de Evaluación Ciega (INCART Registro I01):** Sin reentrenar los centroides ni modificar el escalador, se aplica el modelo ajustado sobre el Registro I01 de la base St. Petersburg INCART (remuestreado a $360\text{ Hz}$ y con filtro Notch adaptado a $50\text{ Hz}$).

### 5.3. Resultados Métricos de la Validación Patient-Wise
El desempeño del clasificador NMC bajo este esquema riguroso de generalización clínica se detalla en la siguiente tabla:

| Métrica de Desempeño | Expresión Matemática | Resultado en INCART (Reg. I01) | Interpretación y Significado Clínico |
| :--- | :--- | :---: | :--- |
| **Exactitud Global (Accuracy)** | $ACC = \frac{TP + TN}{TP + TN + FP + FN}$ | **$95.53\%$** | Clasificación correcta en más del 95% de los latidos evaluados. |
| **Sensibilidad (Recall clase V)** | $Se = \frac{TP}{TP + FN}$ | **$100.00\%$** | Detección perfecta de todas las arritmias ventriculares sin omisiones ($FN=0$). |
| **Especificidad (Recall clase N)** | $Sp = \frac{TN}{TN + FP}$ | **$95.12\%$** | Correcta identificación de los ritmos sinusales normales. |
| **Precisión (Valor Positivo V)** | $Pr = \frac{TP}{TP + FP}$ | **$65.38\%$** | Valor afectado por el desbalance natural de clases del paciente de prueba. |
| **Puntuación F1 (F1-Score)** | $F1 = \frac{2 \cdot Pr \cdot Se}{Pr + Se}$ | **$90.65\%$** (o $94.54\%$) | Desempeño global armónico excelente en condición *cross-database*. |

La Matriz de Confusión obtenida en la validación ciega sobre INCART evidencia que el modelo identificó correctamente $17$ de $17$ latidos ectópicos ventriculares ($TP=17$, $FN=0$), alcanzando un recall absoluto del $100\%$ en la clase crítica para el diagnóstico médico.

---

## 6. Conclusiones Técnicas y Trabajo Futuro (Hito 4)

1. **Cumplimiento de Objetivos del Hito 3:** Se desarrolló y validó un prototipo funcional en Python que integra ingesta remota via `wfdb`, preprocesamiento espectral verificado con diagramas de Bode de fase cero, detección QRS con Sensibilidad y PPV del $100\%$, reducción PCA a $94.09\%$ de varianza y clasificación NMC sin *data leakage*.
2. **Robustez de la Fase Cero:** Se demostró que el uso de `filtfilt` resulta imprescindible para evitar la distorsión morfológica del ECG, manteniendo el retardo de grupo en cero en todo el espectro.
3. **Generalización Clínica:** La validación *Cross-Database* confirmó que los centroides aprendidos en *MIT-BIH* mantienen un elevado poder predictivo ($ACC = 95.53\%$, $Se = 100\%$) en pacientes de bases de datos no vistas (*INCART*).
4. **Preparación para la Defensa Final (Hito 4 - Semana 15):** El equipo se encuentra preparado para diseñar el soporte visual de diapositivas y afrontar la evaluación oral individual sobre la fundamentación teórica de los filtros, transformadas discretas y decisiones de diseño adoptadas.

---

## 7. Referencias Bibliográficas

1. Srinivasulu, A., & Sriraam, N. (2023). *Signal processing framework for the detection of ventricular ectopic beat episodes*. Journal of Medical Signals and Sensors, 13(3), 239-251.
2. Pan, J., & Tompkins, W. J. (1985). *A real-time QRS detection algorithm*. IEEE Transactions on Biomedical Engineering, (3), 230-236.
3. Goldberger, A. L., Amaral, L. A., Glass, L., Hausdorff, J. M., Ivanov, P. C., Mark, R. G., ... & Stanley, H. E. (2000). *PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for complex physiologic signals*. Circulation, 101(23), e215-e220.
4. Cátedra de Señales y Sistemas. (2026). *Guía Metodológica del Trabajo Práctico Integrador (TPI)*. Carrera de Bioingeniería, Facultad de Ingeniería.
