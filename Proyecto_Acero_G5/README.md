# Proyecto Diseño en Acero (CIV-336), Grupo 5 (L = 10 m, H = 3,0 m)

Todos los archivos `.s2k` están en formato de tablas de **SAP2000 v27**.
Usan la unidad `Tonf, m, C` y el código de diseño de acero AISC 360-22.

Hay dos configuraciones de perfiles:

| Configuración | Columnas | Vigas X / Vigas Y | Riostras X / Riostras Y | Acero |
|---|---|---|---|---|
| **AISC** (sin sufijo) | W14X90 | W21X50 / W21X62 | HSS8X8X1/2 / HSS10X10X1/2 | A572 Gr.50, A500 Gr.C |
| **Chilena** (`_CL`) | HN40x154 | IN50x74 / IN50x80 | CAJ200x200x12 / CAJ250x250x12 | NCh203 A270ES |

La configuración chilena usa **solo perfiles soldados armados con planchas**, ninguno conformado en frío:
- HN e IN, con la nomenclatura tipo ICHA (serie, altura en cm × peso en kgf/m).
- Cajones soldados de 4 planchas para las riostras.

| Archivo | Contenido |
|---|---|
| `G5_Edificio3D[_CL].s2k` | Modelo 3D completo |
| `G5_Eje2_Pushover[_CL].s2k` | Modelo 2D del eje 2 para el análisis plástico (pushover) |
| `Memoria_Tarea1_G5[_CL].tex` | Memoria de la Tarea 1 en LaTeX. Se pega en Overleaf y compila con pdfLaTeX |
| `Memoria_Tarea1_G5[_CL].pdf` | La memoria ya compilada, para revisarla |
| `generar_modelo.py` + `memoria_latex.py` | Generan todo: `python3 generar_modelo.py` (ambas configuraciones) o `python3 generar_modelo.py CL` (solo una). Requiere numpy |
| `resultados_G5[_CL].json` | Resultados del análisis propio, usados para verificar el modelo |

## Cómo obtener el `.sdb` en SAP2000 v27
1. **File → Import → SAP2000 .s2k Text File…** y elegir `G5_Edificio3D_CL.s2k` (o el de la versión AISC).
2. Revisar el log de importación, ejecutar el análisis (F5) y guardar como `.sdb`.
3. Repetir lo mismo con `G5_Eje2_Pushover_CL.s2k`. Luego asignar las rótulas M3 y crear los casos no lineales
   (CGNL y PUSHOVER) como se indica en el Anexo B de la memoria.

Nota: el enunciado exige que el `.s2k` se pueda abrir en SAP2000 19.2.1. Un `.s2k` exportado desde v27 puede traer
tablas o campos que v19 no reconoce. Antes de entregar, conviene probar que el archivo se abre en el PC10.

## Supuestos del modelo
- Las columnas están empotradas en la base. Las uniones viga-columna son rígidas. Las riostras están rotuladas (M2 y M3 liberados).
- Las columnas van giradas en 90°, de modo que el eje fuerte trabaja en los marcos rígidos 2, 3 y 4 (dirección Y).
- La losa se modela como diafragma rígido por nivel. Las fuerzas sísmicas se aplican en un nudo de CM (199, …, 599).
- Las cargas de losa se reparten a las vigas en dos direcciones, con líneas de rotura a 45°.
- Patrones de carga: `PP` (peso propio × 1,10 por conexiones), `SCP` (550 kgf/m²), `L` (5 kPa), `LR` (techo) y `EX`/`EY`.
- Se usan las combinaciones NCh3171 (C1 a C11) y una `ENVOLVENTE`.
- En las elevaciones de los ejes 1 y 5 se supuso que el lado izquierdo del dibujo corresponde al eje A.

## Pendientes / por verificar
- **Sobrecarga de techo transitable.** Se usó 2,0 kPa. Hay que confirmarla con la Tabla 4 de la NCh1537.
  Si es otro valor, se cambia `Q_TECHO_KPA` en `generar_modelo.py` y se vuelve a correr el script.
- **Perfiles soldados.** Las dimensiones se eligieron con espesores comerciales de plancha. Hay que verificar el perfil
  equivalente en el catálogo ICHA usado en clases o especificarlo como perfil a pedido.
  Para cambiarlos, se editan `CL_I_DIMS` y `CL_BOX_DIMS` en `generar_modelo.py` y se vuelve a correr el script.
- **Origen de los esfuerzos.** Los esfuerzos de la memoria salen del análisis matricial propio del script, que reproduce
  el modelo SAP (el equilibrio global está verificado). Después de correr SAP2000 conviene comparar.
- Hay que completar la portada (integrantes y profesor) y agregar las capturas de diagramas en el Anexo C.
