# Hito 1: Definición y Selección del Problema

**Asignatura:** Señales y Sistemas - Bioingeniería  
**Trabajo Práctico Integrador (TPI)** | **Semana 11**  
**Tema General:** Procesamiento de Señales Electrocardiográficas (ECG)  

---

## 1. Datos del Grupo y Asignación de Roles

* **Grupo N°:** [Insertar Número de Grupo]

* **Integrantes y Roles:**
  * **[Fátima Luna Latorre]** - *Moderador/Coordinador* (Gestión de reuniones y cronograma)
  * **[Nombre Apellido]** - *Responsable Técnico / Lead Coder* (Estructura de código en Python)
  * **[Kiara Catalina Furlán]** - *Analista / Documentador (Scribe)* (Registro de decisiones e informe)
  * **[Nombre Apellido]** - *Analista / Documentador (Scribe)* (Registro de decisiones e informe)
  * **[Nicolás Nista]** - *Desarrollador / Procesamiento de Señales*
  * **[Francisco Nicolás Butus]** - *Desarrollador / Análisis de Datos*
  * **[Nombre Apellido]** - *Desarrollador / Análisis de Datos*

---

# Propuestas de problemas potenciales (ECG)

## Propuesta 1: Filtrado digital y supresión de artefactos fisiológicos y de instrumentación

Las señales de electrocardiograma (ECG) pueden verse afectadas por interferencias de la red eléctrica, la respiración y el movimiento muscular, lo que dificulta su interpretación.

Para abordar este problema, se propone diseñar e implementar filtros digitales en Python, como filtros *Notch* y pasa-banda. Se analizarán sus respuestas en frecuencia y fase, con el objetivo de reducir el ruido sin distorsionar la morfología del complejo QRS.

## Propuesta 2: Detección automática del complejo QRS y análisis de la variabilidad de la frecuencia cardíaca (HRV)

Esta propuesta busca identificar los complejos QRS y los intervalos RR para estudiar el ritmo cardíaco y analizar la variabilidad de la frecuencia cardíaca (*Heart Rate Variability*, HRV), un indicador relacionado con la actividad del sistema nervioso autónomo.

Para ello, se propone implementar el algoritmo de Pan-Tompkins para la detección de picos R y analizar la serie de intervalos RR mediante parámetros temporales y herramientas de análisis espectral, como la Transformada Rápida de Fourier (FFT) y la densidad espectral de potencia (PSD).

## Propuesta 3: Detección y clasificación de latidos ectópicos ventriculares mediante el análisis de señales ECG

Los latidos ectópicos ventriculares son alteraciones del ritmo cardíaco que presentan características diferentes de las de los latidos normales. Su detección automática puede facilitar el análisis de registros electrocardiográficos y contribuir a la identificación de posibles irregularidades cardíacas.

Se propone analizar señales ECG mediante técnicas de procesamiento digital para identificar y diferenciar latidos normales de latidos ectópicos ventriculares. Para ello, se podrán explorar métodos de filtrado, extracción de características temporales y morfológicas, y algoritmos de detección y clasificación.

## Resultado de la selección

Tras debatir y analizar conjuntamente las propuestas con el profesor, se determinó que la **Propuesta 3: Detección y clasificación de latidos ectópicos ventriculares mediante el análisis de señales ECG** es la más pertinente y adecuada para el alcance del curso.
