# Proyecto Diseño en Acero (CIV-336), Grupo 5 (L = 10 m, H = 3,0 m)

| Archivo | Contenido |
|---|---|
| `G5_Edificio3D.s2k` | Modelo 3D completo para SAP2000 (formato de tablas v19.2.1) |
| `G5_Eje2_Pushover.s2k` | Modelo 2D del eje 2 para el análisis plástico (pushover) |
| `Memoria_Tarea1_G5.tex` | Memoria de la Tarea 1 en LaTeX. Es autocontenida: se pega en Overleaf y compila con pdfLaTeX |
| `Memoria_Tarea1_G5.pdf` | La memoria ya compilada, para revisarla |
| `generar_modelo.py` + `memoria_latex.py` | Script que genera todo lo anterior (`python3 generar_modelo.py`, requiere numpy) |
| `resultados_G5.json` | Resultados del análisis propio, usados para verificar el modelo |

## Cómo obtener el `.sdb`
1. Abrir SAP2000 → **File → Import → SAP2000 .s2k Text File…** → elegir `G5_Edificio3D.s2k`.
2. Revisar el log de importación, ejecutar el análisis (F5) y guardar como `G5_Edificio3D.sdb`.
3. Repetir lo mismo con `G5_Eje2_Pushover.s2k`. Luego asignar las rótulas M3 y crear los casos no lineales
   (CGNL y PUSHOVER) como se indica en el Anexo B de la memoria.

## Supuestos del modelo
- Las columnas están empotradas en la base. Las uniones viga-columna son rígidas. Las riostras están rotuladas (M2 y M3 liberados).
- Las columnas van giradas en 90°, de modo que el eje fuerte trabaja en los marcos rígidos 2, 3 y 4 (dirección Y).
- La losa se modela como diafragma rígido por nivel. Las fuerzas sísmicas se aplican en un nudo de CM (199, 299, …, 599).
- Las cargas de losa se reparten a las vigas en dos direcciones, con líneas de rotura a 45°.
- Patrones de carga: `PP` (peso propio × 1,10 por conexiones), `SCP` (losa de 450 + tabiques de 100 kgf/m²), `L` (5 kPa), `LR` (techo) y `EX`/`EY`.
- Se usan las combinaciones NCh3171 (C1 a C11) y una `ENVOLVENTE`.
- En las elevaciones de los ejes 1 y 5 se supuso que el lado izquierdo del dibujo corresponde al eje A.

## Pendientes / por verificar
- **Sobrecarga de techo transitable.** Se usó 2,0 kPa. Hay que confirmarla con la Tabla 4 de la NCh1537 vista en clases.
  Si es otro valor, se cambia `Q_TECHO_KPA` en `generar_modelo.py` y se vuelve a correr el script. Así se regeneran el .s2k y la memoria.
- Los esfuerzos de la memoria salen del análisis matricial propio del script, que reproduce el modelo SAP
  (el equilibrio global está verificado). Después de correr SAP2000 conviene comparar y, si hay diferencias, usar los valores de SAP.
- Hay que completar la portada (integrantes y profesor) y agregar las capturas de diagramas en el Anexo C.
