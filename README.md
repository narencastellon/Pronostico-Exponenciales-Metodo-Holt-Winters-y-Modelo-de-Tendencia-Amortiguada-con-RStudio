# Pronostico-Exponenciales-Metodo-Holt-Winters-y-Modelo-de-Tendencia-Amortiguada-con-RStudio

## **5. Método estacional de Holt-Winters : técnica para datos con tendencia y estacionalidad.**

Para hacer predicciones utilizando datos con tendencia y estacionalidad, recurrimos al método estacional de Holt-Winters. Este método se puede implementar con una estructura "Aditiva" o una estructura "Multiplicativa", donde la elección del método depende del conjunto de datos. El modelo aditivo se utiliza mejor cuando la tendencia estacional es de la misma magnitud en todo el conjunto de datos, mientras que se prefiere el modelo multiplicativo cuando la magnitud de la estacionalidad cambia a medida que aumenta el tiempo.

Dado que los datos de Google no tienen estacionalidad, usaremos los datos `qcement` que configuramos en la **sección Replicación** para demostrarlo. Estos datos tienen estacionalidad y tendencia; sin embargo, no está claro si la estacionalidad es aditiva o multiplicativa. Usaremos el método de Holt-Winters para identificar el modelo que mejor se ajusta.
```{r}
autoplot(decompose(qcement), main="Descomposición aditiva de Series de tiempo")
```

#### **5.1 Aditivo**

Para el modelo aditivo, la forma de la ecuación regular es:
$$\hat{y}_{t+1}=L_t+kT_t+S_{t+k-m} $$
El nivel $L-t$, tendencia $T_t$ y temporada $S_t$ se actualizan a través de un par de ecuaciones de actualización, que es donde se ve la presencia de los tres parámetros de suavizado:

$$L_t=\alpha(y_t-S_{t-m})+(1->alpha)(L_{t-1}+T_{t-1} ) $$
$$T_t=\beta(L_{t-1}-T_{t-1} )+(1-\beta)T_{t-1} $$
$$S_t=\gamma(y_t-L_t)+(1-\gamma)S_{t-m} $$

dónde $\alpha$, $\beta$ y $\gamma$ son los tres parámetros de suavizado para tratar el patrón de nivel, la tendencia y la estacionalidad, respectivamente. De manera similar a SES y al método de Holt, los tres parámetros están restringidos a 0-1. Las ecuaciones de los componentes son las siguientes:
$$\hat{y}_{t+h}=l_t+hb_t+s_{t-m-hm} $$
$$l_t=\alpha(y_t-s_{t-m})+(1-\alpha)(l_{t-1}+b_{t-1}) $$
$$b_t=\beta(l_t-l_{t-1})+(1-\beta)b_{t-1} $$
$$s_t=\gamma(y_t-l_{t-1}-b_{t-1})+(1-\gamma)s_{t-m} $$

Para aplicar el método Holt-Winters, presentaremos una nueva función, `ets` que significa error, tendencia y estacionalidad. Lo importante que hay que entender sobre el modelo `ets` es cómo seleccionar el `model =`parámetro. En total, tiene 36 opciones de modelos para elegir. La configuración de los parámetros en el siguiente código ( `model = "AAA"`) representa un modelo con error aditivo, tendencia aditiva y estacionalidad aditiva.
```{r}
qcement.hw <- ets(qcement.train, model = "AAA")
autoplot(forecast(qcement.hw), main="Pronóstico del Modelo ETS(A,A,A)")
```

Por lo tanto, al especificar el tipo de modelo, siempre especifica el error, la tendencia y luego la estacionalidad (de ahí "ets"). Las opciones que puede especificar para cada componente son las siguientes:

* **error** : aditivo ("A"), multiplicativo ("M"), desconocido ("Z")
* **tendencia** : ninguna ("N"), aditiva ("A"), multiplicativa ("M"), desconocida ("Z")
* **estacionalidad** : ninguna ("N"), aditivo ("A"), multiplicativo ("M"), desconocido ("Z")

En consecuencia, si desea aplicar un modelo de Holt en el que el error y la tendencia son aditivos y no existe estacionalidad, debe seleccionar `model = "AAN"`. Si desea aplicar un modelo de Holt-Winters donde hay un error aditivo, una tendencia exponencial (multiplicativa) y una estacionalidad aditiva, debe seleccionar `model = "AMA"`. Si no está seguro del tipo de componente, utilice "Z". Entonces, si no estaba seguro de los componentes o si desea que el modelo seleccione la mejor opción, puede usar `model = "ZZZ"` y se seleccionará el modelo "óptimo".

Si evaluamos nuestro modelo aditivo podemos ver que  $\alpha=0.6208$, $\beta=0.0001$ y $\gamma=0.1913$
```{r}
summary(qcement.hw)
```

Si revisamos nuestros residuos, vemos que los residuos crecen con el tiempo. Esto puede sugerir que una tasa de error multiplicativa puede ser más apropiada.
```{r}
checkresiduals(qcement.hw)
```

Si comprobamos la precisión predictiva, vemos que nuestra precisión de predicción es de aproximadamente un 2,9% (según el MAPE).
```{r}
# pronosticar los próximos 5 trimestres
qcement.f1 <- forecast(qcement.hw, h = 5)
```

```{r}
# Revisando la precisión
accuracy(qcement.f1, qcement.test)
```

#### **5.2 Multiplcativo**

Como se indicó anteriormente, es posible que tengamos características multiplicativas para nuestro **método Holt-Winters**. Si tenemos estacionalidad multiplicativa, entonces nuestra forma de ecuación cambia a:
$$\hat{y}_{T+1}=(L_t+kT_t)S_{t+-m} $$

El nivel $L_t$, tendencia $T_t$ y temporada $S_t$ se actualizan a través de un par de ecuaciones de actualización, que es donde se ve la presencia de los tres parámetros de suavizado:

$$L_t=\alpha y_t/S_{t-m}+(1-\alpha)(L_{t-1}+T_{t-1}) $$
$$T_t=\beta(L_t-L_{t-1})+(1-\beta)T_{t-1} $$
$$S_t=\gamma(y_t/L_t)+(1-\gamma)S_{t-m} $$

Si aplicamos un **modelo de estacionalidad multiplicativa**, nuestro parámetro `model` se convierte en `model = "MAM"`(aquí, en realidad estamos aplicando un modelo de error multiplicativo y estacionalidad). Vemos que los residuos ilustran menos cambios de magnitud a lo largo del tiempo. Todavía tenemos un problema con la autocorrelación con errores, pero lo abordaremos en tutoriales posteriores.
```{r}
qcement.hw2 <- ets(qcement.train, model = "MAM")
checkresiduals(qcement.hw2)
```

Para comparar la precisión predictiva de nuestros modelos, comparemos cuatro modelos diferentes. Vemos que el primer modelo (error aditivo, tendencia y estacionalidad) da como resultado el RMSE y MAPE más bajos en nuestro conjunto de datos de prueba.
```{r}
# error aditivo, tendencia y estacionalidad
qcement.hw1 <- ets(qcement.train, model = "AAA")
qcement.f1 <- forecast(qcement.hw1, h = 5)
accuracy(qcement.f1, qcement.test)
```

```{r}
# error multiplicativo, tendencia aditiva y estacionalidad
qcement.hw2 <- ets(qcement.train, model = "MAA")
qcement.f2 <- forecast(qcement.hw2, h = 5)
accuracy(qcement.f2, qcement.test)
```

```{r}
# error aditivo y tendencia y estacionalidad multiplicativa
qcement.hw3 <- ets(qcement.train, model = "AAM", restrict = FALSE)
qcement.f3 <- forecast(qcement.hw3, h = 5)
accuracy(qcement.f3, qcement.test)
```

```{r}
# error multiplicativo, tendencia aditiva y estacionalidad multiplicativa
qcement.hw4 <- ets(qcement.train, model = "MAM")
qcement.f4 <- forecast(qcement.hw4, h = 5)
accuracy(qcement.f4, qcement.test)
```

Si tuviéramos que comparar esto con un modelo no especificado donde dejamos seleccionar `ets` el modelo óptimo, vemos que ets selecciona una especificación de modelo de error **multiplicativo, tendencia aditiva y estacionalidad multiplicativa (“MAM”)**. Esto es equivalente a nuestro cuarto modelo anterior. Se supone que este modelo es "óptimo" porque minimiza RMSE, AIC y BIC en el conjunto de datos de **entrenamiento** , pero no necesariamente minimiza los errores de predicción en el conjunto de prueba.

```{r}
qcement.hw5 <- ets(qcement.train, model = "ZZZ")
summary(qcement.hw5)
```

Como hicimos en la sección del **método SES y Holt** , podemos optimizar la $\gamma$ parámetro en nuestro **modelo de Holt-Winters**. Aquí, utilizamos el modelo de error aditivo, tendencia y estacionalidad que minimizó nuestros errores de predicción anteriores e identificamos el parámetro $\gamma$ que minimiza los errores de pronóstico. En este caso vemos que $\gamma=0.21$ minimiza la tasa de error.
 
```{r}
gamma <- seq(0.01, 0.85, 0.01)
RMSE <- NA

for(i in seq_along(gamma)) {
  hw.expo <- ets(qcement.train, "AAA", gamma = gamma[i])
  future <- forecast(hw.expo, h = 5)
  RMSE[i] = accuracy(future, qcement.test)[2,2]
}

error <- data.frame(gamma, RMSE)
minimum <- filter(error, RMSE==min(RMSE))
ggplot(error, aes(gamma, RMSE)) + geom_line() + geom_point(data=minimum, color = "blue", size = 2) + ggtitle("Pronósticos de los errores gamma", subtitle = "gamma = 0.21 minimizes RMSE")
```

Si actualizamos nuestro modelo con este $\gamma$ "óptimo" vemos que reducimos nuestra tasa de error de pronóstico del 2,88% al 2,76%. Se trata de una pequeña mejora, pero a menudo las pequeñas mejoras pueden tener grandes implicaciones comerciales.
```{r}
# modelo anterior con error aditivo, tendencia y estacionalidad
accuracy(qcement.f1, qcement.test)
```

```{r}
# Nuevo modelo con óptimo parámetro gamma
qcement.hw6 <- ets(qcement.train, model = "AAA", gamma = 0.21)
qcement.f6 <- forecast(qcement.hw6, h = 5)
accuracy(qcement.f6, qcement.test)
```

Con este nuevo modelo óptimo podemos obtener nuestros valores predichos:
```{r}
qcement.f6
```

y también visualice estos valores predichos:
```{r}
autoplot(qcement.f6)
```


## **6.Métodos de tendencia amortiguada : técnica para las tendencias que se cree que se vuelven más conservadoras o “planas” con el tiempo.**

Un último elemento para discutir es la idea de "amortiguar" su pronóstico. Los pronósticos amortiguados utilizan un coeficiente de amortiguamiento denotado $\phi$ para estimar de forma más conservadora la tendencia prevista. Básicamente, si cree que su tendencia aditiva o multiplicativa se está desacelerando o disminuirá ("línea plana") en el futuro cercano, entonces está asumiendo que disminuirá.

La forma de la ecuación para un modelo aditivo con un coeficiente de amortiguamiento es
$$\hat{y}_{t+h}=L_t+(\phi+\phi^2+...+\phi^h)\beta_t $$
$$L_t=\alpha y_t+\alpha(1-\alpha)(L_{t-1}+\phi\beta_{t-1}) $$
$$\beta_t=\beta(L_t-L_{t-1})+(1-\beta)\phi\beta_{t-1}. $$

donde $0< \phi <1$. Cuando $\phi=1$ el método es el mismo que el modelo aditivo de Holt. Como $\phi$ se acerca a 0, la tendencia se vuelve más conservadora y se convierte en una constante en un futuro más cercano. El resultado final de este método es que los pronósticos a corto plazo siguen teniendo tendencia, mientras que los pronósticos a largo plazo son constantes.

Para ilustrar el efecto de un pronóstico amortiguado usaremos el `fpp2::ausair` conjunto de datos. Aquí, creamos varios modelos (aditivo, aditivo + amortiguado, multiplicativo, multiplicativo + amortiguado). En el gráfico, puede ver que los modelos amortiguados (líneas discontinuas) tienen líneas de tendencia más conservadoras y si las pronosticamos lo suficientemente en el futuro, veríamos esta tendencia plana.

```{r}
# Modelo lineal (aditivo) de Holt 
fit1 <- ets(ausair, model = "ZAN", alpha = 0.8, beta = 0.2)
pred1 <- forecast(fit1, h = 5)

# Modelo lineal (aditivo) de Holt 
fit2 <- ets(ausair, model = "ZAN", damped = TRUE, alpha = 0.8, beta = 0.2, phi = 0.85)
pred2 <- forecast(fit2, h = 5)

# Modelo exponencial(Multiplicativo) de Holt 
fit3 <- ets(ausair, model = "ZMN", alpha = 0.8, beta = 0.2)
pred3 <- forecast(fit3, h = 5)

# modelo exponencial (multiplicativo) de holt amortiguado
fit4 <- ets(ausair, model = "ZMN", damped = TRUE, alpha = 0.8, beta = 0.2, phi = 0.85)
pred4 <- forecast(fit4, h = 5)

autoplot(ausair) +
  autolayer(pred1$mean, color = "blue") +
  autolayer(pred2$mean, color = "blue", linetype = "dashed") +
  autolayer(pred3$mean, color = "red") +
  autolayer(pred4$mean, color = "red", linetype = "dashed")
```

Los modelos anteriores fueron solo para fines ilustrativos. Aplicaría el mismo proceso que vio en secciones anteriores para identificar si un modelo amortiguado predice con mayor precisión que un modelo no amortiguado. Incluso puede aplicar los enfoques que vio anteriormente para ajustar este parámetro para identificar el $\phi$ óptimo coeficiente.

