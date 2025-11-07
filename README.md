## Descripción del Proyecto

Todo el código de las funciones utilizadas para ajustar los modelos se encuentra en el archivo **`FQLGR.R`**.

Cada documento **RMarkdown** genera los gráficos correspondientes y los guarda automáticamente en las carpetas que comienzan con `plots_...`.

> **Nota:**  
> El proyecto `fake-audio-discriminator` utiliza los datos de la carpeta comprimida.  
> Debido al peso de los datos, no fue posible ejecutar todos los gráficos; por lo tanto, la carpeta **`plots_fakeaudio`** está desactualizada y no incluye los resultados del modelo  
> \( \Upsilon^{(3)} \).

---

## Modelos

Los modelos estimados están basados en la siguiente formulación general:

$$
\eta_i = \alpha + \boldsymbol{\xi}_i^{(i)} \Upsilon_\beta^{(i)} + Q(\boldsymbol{\xi}_i^{(i)}) \Upsilon_\gamma^{(i)}
$$

donde:

- $ \beta = H_\beta^{(i)}(\Upsilon^{(i)}) $
- $ \gamma = H_\gamma^{(i)}(\Upsilon^{(i)}) $
- $ \xi^{(i)} = S^{(i)}(X) $ son las proyecciones de las variables originales bajo distintas transformaciones según el modelo.

---

## Operador Cuadrático

El operador $Q(X)$ toma una matriz $X$ y genera una nueva matriz cuyas columnas corresponden a las multiplicaciones entre todos los pares únicos de columnas de $X$:

$$
Q(X) = [\, X_{\cdot,i} \cdot X_{\cdot,j} \,]_{1 \leq i \leq j}
$$

---

## Optimización

Todos los modelos se resuelven mediante un **algoritmo de descenso por gradiente penalizado por suavidad**,  
donde las penalizaciones dependen de las transformaciones $$ H^{(i)}(X) $$ y $$ S^{(i)}(X) $$ empleadas en cada caso.
