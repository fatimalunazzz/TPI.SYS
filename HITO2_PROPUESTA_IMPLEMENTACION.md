# Hito 2: Propuesta de Implementación del Sistema

**Asignatura:** Señales y Sistemas - Bioingeniería (Plan 2026)  
**Trabajo Práctico Integrador (TPI)** | **Semana 12**  
**Título del Proyecto:** Detección y Clasificación Automática de Latidos Ventriculares Ectópicos (VEB/PVC) mediante Procesamiento Digital de Señales de ECG  
**Base de Datos Seleccionada:** *MIT-BIH Arrhythmia Database* (PhysioNet)  
**Artículo de Referencia:** Srinivasulu, A., & Sriraam, N. (2023). *Signal processing framework for the detection of ventricular ectopic beat episodes*. Journal of Medical Signals & Sensors, 13(3), 239-251.


---

## 1. Justificación de la Elección del Artículo de Referencia

Para el desarrollo del TPI se seleccionó el trabajo de **Srinivasulu & Sriraam (2023)** por las siguientes razones fundamentales:

1. **Accesibilidad y Claridad Metodológica:** El artículo presenta un flujo de procesamiento directo y transparente, detallando en forma explícita cada ecuación matemática y etapa de procesamiento [JMSS-13-239].
2. **Pertinencia Teórica:** Centra el análisis en conceptos del programa analítico de *Señales y Sistemas* (filtrado digital pasa-banda/notch, derivación discreta, integración en ventana móvil y análisis temporal/morfométrico), evitando dependencias excesivas de arquitecturas opacas de aprendizaje profundo.
3. **Eficiencia en la Extracción de Características:** Utiliza únicamente **5 características fisiológicas y temporales** derivadas del complejo QRS, lo que simplifica la computación y permite interpretar clínicamente el comportamiento de la señal electrocardiográfica [JMSS-13-239].
4. **Validación sobre Datos Abiertos:** Emplea la base estándar pública *MIT-BIH Arrhythmia Database* de PhysioNet, garantizando la reproducibilidad del prototipo [JMSS-13-239].

---

## 2. Selección y Descripción de la Base de Datos

* **Nombre:** MIT-BIH Arrhythmia Database (PhysioNet)
* **Acceso:** PhysioNet Open Access Repository (MITDB 1.0.0)
* **Características Técnicas de los Registros:**
  * **Registros:** 48 registros de ECG de 30 minutos de duración cada uno.
  * **Frecuencia de Muestreo ($f_s$):** $360\text{ Hz}$ por canal.
  * **Resolución:** 11 bits en un rango de $\pm 5\text{ mV}$.
  * **Derivación Principal:** Derivación modificada de extremidades II (MLII).
  * **Anotaciones:** Marcas de latidos anotadas por cardiólogos clasificando latidos normales ($N$) y latidos ectópicos ventriculares ($V$ / PVC).

* **Criterios de Inclusión/Exclusión de Registros:**  
  Siguiendo las pautas del artículo de referencia, se excluirán los registros con marcapasos cardíacos (102, 104, 107 y 217) y los segmentos con severo ruido de movimiento no recuperable [JMSS-13-239].

---

## 3. Diagrama en Bloques del Sistema Propuesto

El flujo de procesamiento desde la señal bruta hasta la clasificación del latido se compone de 5 bloques funcionales:

```
+-------------------------------------------------------------------+
| 1. Carga e Ingesta de Señal ECG Cruda (PhysioNet WFDB / MLII)     |
+-------------------------------------------------------------------+
                                  |
                                  v
+-------------------------------------------------------------------+
| 2. Preprocesamiento y Acondicionamiento Digital de la Señal       |
|    - Filtro Pasa-Banda FIR/IIR (0.5 Hz - 40 Hz) [Deriva/EMG]      |
|    - Filtro Notch IIR (50/60 Hz) [Interferencia de Red]          |
+-------------------------------------------------------------------+
                                  |
                                  v
+-------------------------------------------------------------------+
| 3. Detección de Eventos y Puntos Fiduciales (Algoritmo Pan-Tompkins) |
|    - Derivación discreta + Elevación al cuadrado                   |
|    - Integración por ventana móvil + Umbral adaptativo             |
|    - Localización precisa de picos R, Q y S                       |
+-------------------------------------------------------------------+
                                  |
                                  v
+-------------------------------------------------------------------+
| 4. Extracción del Vector de 5 Características Morfométricas       |
|    - [Pre-RR, Post-RR, Ancho QRS, Pendiente QR, Pendiente RS]     |
+-------------------------------------------------------------------+
                                  |
                                  v
+-------------------------------------------------------------------+
| 5. Clasificación Automática de Latidos (NMC / NRMS)              |
|    - Comparación de distancia Euclídea respecto a los centroides   |
|    - Salida: Latido Normal (N) vs. Ectópico Ventricular (VEB/V)   |
+-------------------------------------------------------------------+
```

---

## 4. Fundamentación Teórica y Definición Matemática

### 4.1 Preprocesamiento y Filtrado Digital
La señal electrocardiográfica $x[n]$ capturada sufre contaminación por ruido de baja frecuencia (deriva de línea base por respiración, $< 0.5\text{ Hz}$) e interferencias de alta frecuencia (ruido electromiográfico EMG y red eléctrica de $50/60\text{ Hz}$).

1. **Filtro Pasa-Banda (0.5 Hz - 40 Hz):**  
   Atenúa la deriva y el ruido EMG manteniendo la energía fundamental del complejo QRS. Se implementa mediante un filtro digital cuya respuesta en frecuencia atenúa las componentes fuera del rango deseado:
   $$y_{pb}[n] = h_{pb}[n] * x[n] = \sum_{k=0}^{M} h_{pb}[k] x[n-k]$$

2. **Filtro Notch (Filtro Muesca en 50/60 Hz):**  
   Elimina el zumbido de la red eléctrica mediante pares de ceros complejos conjugados sobre la circunferencia de radio unidad en $e^{\pm j \omega_0}$ donde $\omega_0 = 2\pi f_{red} / f_s$:
   $$H_{notch}(z) = G \frac{1 - 2\cos(\omega_0)z^{-1} + z^{-2}}{1 - 2r\cos(\omega_0)z^{-1} + r^2 z^{-2}}$$
   con $r \approx 0.95 - 0.99$ para controlar el ancho de banda de la muesca.

---

### 4.2 Algoritmo de Detección del Complejo QRS (Pan-Tompkins)
Para localizar los picos $R$ y determinar los puntos $Q$ y $S$, se aplican las siguientes transformaciones discretas sobre la señal filtrada $y[n]$:

1. **Derivador Discreto (Diferenciación):**  
   Resalta la alta pendiente característica de la onda R:
   $$d[n] = \frac{1}{8T} \Big(2 y[n] + y[n-1] - y[n-3] - 2 y[n-4]\Big)$$

2. **Operación de Elevación al Cuadrado:**  
   Intensifica las componentes de alta frecuencia producidas por la pendiente del QRS y fuerza valores positivos:
   $$s[n] = \big(d[n]\big)^2$$

3. **Integración por Ventana Móvil:**  
   Obtiene información sobre el ancho del complejo QRS sumando puntos en una ventana de duración $W \approx 150\text{ ms}$ ($N = \lfloor 0.150 \cdot f_s \rfloor$ muestras):
   $$m[n] = \frac{1}{N} \sum_{k=0}^{N-1} s[n-k]$$

4. **Detección de Picos Fiduciales ($Q, R, S$):**  
   Tras aplicar un umbral adaptativo sobre $m[n]$ para aislar el complejo, el pico $R_i$ se re-ubica buscando el máximo local en la señal filtrada original. El punto $Q_i$ se define como el mínimo local a la izquierda de $R_i$, y el punto $S_i$ como el mínimo local a la derecha de $R_i$.

---

### 4.3 Extracción de Características Temporales y Morfométricas
Para cada latido $i$, se calcula un vector de características $\mathbf{f}_i \in \mathbb{R}^5$ basado en las propiedades anatómicas descritas por Srinivasulu & Sriraam (2023) [JMSS-13-239]:

1. **Intervalo Pre-RR:** Tiempo transcurrido desde el latido previo.
   $$\text{Pre-RR}(i) = t_{R}(i) - t_{R}(i-1)$$

2. **Intervalo Post-RR:** Tiempo transcurrido hasta el latido siguiente.
   $$\text{Post-RR}(i) = t_{R}(i+1) - t_{R}(i)$$

3. **Ancho del Complejo QRS:** Duración del complejo $QRS$.
   $$\text{Ancho}_{\text{QRS}}(i) = t_{S}(i) - t_{Q}(i)$$

4. **Pendientes Morfológicas ($QR_{slope}$ y $RS_{slope}$):**  
   Inclinación de los flancos ascendente y descendente del pico R:
   $$\text{Pendiente}_{\text{QR}}(i) = \frac{A_R(i) - A_Q(i)}{t_R(i) - t_Q(i)}$$
   $$\text{Pendiente}_{\text{RS}}(i) = \frac{A_R(i) - A_S(i)}{t_R(i) - t_S(i)}$$

---

### 4.4 Clasificador Replicado: Nearest Mean Classifier (NMC)
De acuerdo con las recomendaciones recibidas para acotar el alcance al marco de *Señales y Sistemas*, se implementará el **Nearest Mean Classifier (NMC)** por su transparencia lineal y bajo costo computacional ($CT \approx 0.0031\text{ s}$) [JMSS-13-239].

1. **Etapa de Entrenamiento (Cálculo de Centroides):**  
   A partir del conjunto de entrenamiento con etiquetas conocidas $k \in \{\text{Normal (N)}, \text{Ectópico (V)}\}$, se calcula el vector promedio de características $\mathbf{\bar{f}}_k \in \mathbb{R}^5$:
   $$\mathbf{\bar{f}}_k = \frac{1}{N_k} \sum_{j \in \text{Clase } k} \mathbf{f}_j$$

2. **Etapa de Regla de Decisión (Clasificación de un Latido Nuevo):**  
   Dado un latido de prueba con vector de características $\mathbf{f}_{\text{test}}$, se calcula la distancia Euclídea hacia cada centroide:
   $$d(\mathbf{f}_{\text{test}}, \mathbf{\bar{f}}_k) = \sqrt{\sum_{m=1}^{5} \big(f_{\text{test}}[m] - \bar{f}_k[m]\big)^2}$$

3. **Asignación de Clase:**  
   $$\text{Clase Predicha} = \arg\min_{k \in \{N, V\}} d(\mathbf{f}_{\text{test}}, \mathbf{\bar{f}}_k)$$

---

## 5. Referencias Bibliográficas

1. **Srinivasulu, A., & Sriraam, N. (2023).** *Signal processing framework for the detection of ventricular ectopic beat episodes*. Journal of Medical Signals & Sensors, 13(3), 239–251.
2. **Pan, J., & Tompkins, W. J. (1985).** *A real-time QRS detection algorithm*. IEEE Transactions on Biomedical Engineering, (3), 230–236.
3. **Goldberger, A. L., et al. (2000).** *PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals*. Circulation, 101(23), e215–e220.
