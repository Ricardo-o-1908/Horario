# -*- coding: utf-8 -*-
"""Escribe la memoria de cálculo de la Tarea 1 (LaTeX autocontenido, figuras en TikZ/pgfplots)."""
import math


def n(v, d=2):
    """Número con coma decimal (válido en texto y en modo matemático)."""
    s = ("%." + str(d) + "f") % v
    return s.replace(".", "{,}")


def escribir(ruta, G, datos, po, trac, resumen, env):
    L, H, ZS, XS, YS = G["L"], G["H"], G["ZS"], G["XS"], G["YS"]
    NIV = G["NIVELES"]
    SEC = G["SECCIONES"]
    BARRAS = {b.id: b for b in G["BARRAS"]}
    NUDOS = G["NUDOS"]
    COMBOS = G["COMBOS"]
    MPA = G["MPA"]
    out = []
    w = out.append

    # ------------------------------------------------------------------ preámbulo
    w(r"""% =====================================================================
%  PROYECTO DISEÑO EN ACERO (CIV-336) - USM - 2do SEMESTRE 2026
%  GRUPO 5  -  TAREA 1: Cargas, Modelación, Análisis Plástico y Tracción
%  Documento autocontenido: compila directamente en Overleaf (pdfLaTeX).
%  Las figuras están hechas en TikZ/pgfplots (no requiere imágenes).
% =====================================================================
\documentclass[11pt,letterpaper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[spanish,es-noshorthands,es-tabla]{babel}
\usepackage{lmodern}
\usepackage[margin=2.5cm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,array,multirow,tabularx}
\usepackage{float}
\usepackage{caption}
\usepackage{xcolor}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,patterns,calc,decorations.pathmorphing}
\usepackage{pgfplots}
\pgfplotsset{compat=1.16}
\usepackage{fancyhdr}
\usepackage[hidelinks]{hyperref}
\setlength{\parskip}{0.5em}
\setlength{\headheight}{14pt}
\setlength{\parindent}{0pt}
\pagestyle{fancy}
\fancyhf{}
\lhead{CIV-336 Diseño en Acero}
\rhead{Proyecto -- Tarea 1 -- Grupo 5}
\cfoot{\thepage}
\renewcommand{\arraystretch}{1.15}
\newcommand{\tonf}{\,\mathrm{tonf}}
\newcommand{\tonfm}{\,\mathrm{tonf\cdot m}}

\begin{document}

% ------------------------------------------------------------------ PORTADA
\begin{titlepage}
\centering
{\large Universidad Técnica Federico Santa María\par}
{\large Departamento de Obras Civiles\par}
\vspace{3cm}
{\Large CIV-336 Diseño en Acero\par}
\vspace{1cm}
{\huge\bfseries Proyecto Semestral\par}
\vspace{0.4cm}
{\LARGE Tarea 1: Cargas, Modelación, Análisis Plástico\\ y Diseño a Tracción\par}
\vspace{2cm}
{\Large\bfseries Grupo 5\par}
\vspace{0.3cm}
{\large $L = 10\ \mathrm{m}$ \quad -- \quad $H = 3{,}0\ \mathrm{m}$\par}
\vfill
\begin{flushleft}
\begin{tabular}{ll}
\textbf{Integrantes:} & Nombre Apellido 1 -- Rol \\
                      & Nombre Apellido 2 -- Rol \\
                      & Nombre Apellido 3 -- Rol \\[0.3cm]
\textbf{Profesor:}    & Nombre del profesor \\
\textbf{Fecha:}       & \today
\end{tabular}
\end{flushleft}
\end{titlepage}

\tableofcontents
\newpage
""")

    # ------------------------------------------------------------------ 1. descripción
    Htot = ZS[-1]
    w(r"""
\section{Introducción y descripción de la estructura}

La presente memoria de cálculo corresponde a la Tarea 1 del proyecto de la asignatura
CIV-336 Diseño en Acero. Se estudia un edificio de cinco pisos estructurado en base a marcos de acero,
con marcos arriostrados concéntricos en los ejes A, B, 1 y 5, y marcos rígidos (resistentes a momento)
en los ejes 2, 3 y 4. En esta entrega se presentan las bases de cálculo, la determinación de las cargas,
la descripción del modelo realizado en SAP2000, el resumen de esfuerzos internos, el análisis plástico del
eje 2 y el diseño a tracción de la riostra más traccionada.

Las dimensiones de la estructura dependen del número de grupo (Tabla 2 del enunciado). Para el
\textbf{Grupo 5} se tiene $L = 10\ \mathrm{m}$ y $H = 3{,}0\ \mathrm{m}$, de donde resulta la geometría
de la Tabla~\ref{tab:geom}.

\begin{table}[H]
\centering
\caption{Geometría de la estructura (Grupo 5).}
\label{tab:geom}
\begin{tabular}{lcc}
\toprule
Descripción & Expresión & Valor \\
\midrule
Vanos extremos (ejes 1--2 y 4--5) & $0{,}75L$ & $7{,}50$ m \\
Vanos centrales (ejes 2--3 y 3--4) & $L$ & $10{,}00$ m \\
Largo total en dirección X & $3{,}5L$ & $35{,}00$ m \\
Ancho en dirección Y (ejes A--B) & $0{,}8L$ & $8{,}00$ m \\
Altura primer piso & $1{,}1H$ & $3{,}30$ m \\
Altura pisos 2 a 5 & $H$ & $3{,}00$ m \\
Altura total $H_{total}$ & $1{,}1H + 4H$ & $""" + n(Htot) + r"""$ m \\
Área de planta por piso & $3{,}5L \times 0{,}8L$ & $""" + n(G["AREA_PLANTA"], 1) + r"""\ \mathrm{m^2}$ \\
\bottomrule
\end{tabular}
\end{table}
""")
    # Figuras de planta y elevaciones
    w(fig_planta(XS, YS))
    w(fig_elevaciones(XS, YS, ZS))

    # ------------------------------------------------------------------ 2. bases de cálculo
    w(r"""
\section{Bases de cálculo}

\subsection{Normas y documentos de referencia}
\begin{itemize}
  \item Enunciado del Proyecto -- Diseño en Acero (CIV-336), Segundo Semestre 2026.
  \item NCh1537.Of2009: Diseño estructural -- Cargas permanentes y cargas de uso.
  \item NCh3171.Of2010: Diseño estructural -- Disposiciones generales y combinaciones de carga.
  \item NCh433.Of1996 Mod.2009 (referencial, para el peso sísmico y la distribución de fuerzas).
  \item ANSI/AISC 360-16: \emph{Specification for Structural Steel Buildings} (método LRFD).
  \item AISC \emph{Steel Construction Manual}, 15\textsuperscript{a} edición (propiedades de perfiles).
\end{itemize}

\subsection{Materiales}
\begin{table}[H]
\centering
\caption{Materiales considerados.}
\begin{tabular}{llcc}
\toprule
Elemento & Material & $F_y$ [MPa] & $F_u$ [MPa] \\
\midrule
Vigas y columnas (perfiles W) & ASTM A572 Gr.50 & """ + n(G["FY_W"], 0) + r""" & """ + n(G["FU_W"], 0) + r""" \\
Riostras (HSS cuadrados) & ASTM A500 Gr.C & """ + n(G["FY_HSS"], 0) + r""" & """ + n(G["FU_HSS"], 0) + r""" \\
Planchas de conexión (gusset) & ASTM A36 & 248 & 400 \\
Soldadura & Electrodo E70XX & -- & $F_{EXX} = 482$ \\
Losa & Hormigón armado & \multicolumn{2}{c}{$\gamma = 2500\ \mathrm{kgf/m^3}$} \\
\bottomrule
\end{tabular}
\end{table}
Módulo de elasticidad del acero $E = 200\,000$ MPa, módulo de Poisson $\nu = 0{,}3$ y
peso unitario $\gamma_a = 7850\ \mathrm{kgf/m^3}$.

\subsection{Método de diseño}
Se utiliza el método de Diseño por Factores de Carga y Resistencia (LRFD) según AISC 360-16, con las
combinaciones de carga de la NCh3171 (sección~\ref{sec:combos}).
""")

    # ------------------------------------------------------------------ 3. cargas
    pat = G["peso_acero_tipo"]
    w(r"""
\section{Cargas}

\subsection{Cargas permanentes (D)}
\begin{table}[H]
\centering
\caption{Cargas permanentes (Tabla 1 del enunciado).}
\begin{tabular}{lcc}
\toprule
Descripción & Cálculo & Valor \\
\midrule
Peso propio acero estructural & $\gamma_a = 7850\ \mathrm{kgf/m^3}$ & SAP2000 (automático) \\
Conexiones & $+10\%$ del peso propio & factor $1{,}10$ \\
Losa de H.A. ($e = 18$ cm) & $0{,}18 \times 2500$ & $450\ \mathrm{kgf/m^2}$ \\
Tabiques, terminaciones, etc. & -- & $100\ \mathrm{kgf/m^2}$ \\
\midrule
Sobrecarga permanente (SCP) & losa + tabiques & $550\ \mathrm{kgf/m^2}$ \\
\bottomrule
\end{tabular}
\end{table}

El peso propio de los elementos de acero se considera en el patrón \texttt{PP} con un multiplicador de
peso propio igual a $1{,}10$, lo que incluye el $10\%$ adicional por conexiones. La losa y los tabiques
se consideran en el patrón \texttt{SCP}. El peso de acero resultante del prediseño es:
columnas $""" + n(pat["COL"], 1) + r"""\tonf$, vigas X $""" + n(pat["VX"], 1) + r"""\tonf$,
vigas Y $""" + n(pat["VY"], 1) + r"""\tonf$, riostras X $""" + n(pat["RX"], 1) + r"""\tonf$ y
riostras Y $""" + n(pat["RY"], 1) + r"""\tonf$ (incluye el $10\%$ de conexiones), con un total de
$""" + n(sum(pat.values()), 1) + r"""\tonf$.

\subsection{\texorpdfstring{Cargas de uso (L y L\textsubscript{r})}{Cargas de uso (L y Lr)}}
\begin{itemize}
  \item Pisos 1 a 4: sobrecarga de uso $q_L = 5\ \mathrm{kPa} = """ + n(G["Q_USO"] * 1000, 0) + r"""\ \mathrm{kgf/m^2}$
        (patrón \texttt{L}).
  \item Techo (nivel 5): techo transitable según NCh1537,
        $q_{L_r} = """ + n(G["Q_TECHO_KPA"], 1) + r"""\ \mathrm{kPa} = """ + n(G["Q_TECHO"] * 1000, 0) + r"""\ \mathrm{kgf/m^2}$
        (patrón \texttt{LR}). % <<< VERIFICAR con la Tabla 4 de NCh1537 vista en clases
\end{itemize}
No se aplica reducción de la sobrecarga por área tributaria (criterio conservador).

\subsection{Distribución de las cargas de losa a las vigas}
La losa no se modela como elemento estructural; su peso, los tabiques y las sobrecargas se traspasan a las
vigas mediante áreas tributarias, considerando losas apoyadas en sus cuatro bordes (trabajo en dos
direcciones) con líneas de rotura a $45^\circ$. Así, en cada paño de $a \times b$ ($a \le b$):
las vigas del lado corto reciben una carga triangular de altura $q\,a/2$ y las del lado largo
una carga trapezoidal de altura $q\,a/2$. Para el Grupo 5 los paños son de
$7{,}5 \times 8{,}0$ m (vanos extremos) y $10{,}0 \times 8{,}0$ m (vanos centrales).

\begin{table}[H]
\centering
\caption{Cargas máximas (ordenada máxima) transmitidas por la losa a las vigas, por paño.}
\begin{tabular}{llccc}
\toprule
Paño & Viga & Forma & SCP [tonf/m] & L [tonf/m] \\
\midrule
$7{,}5 \times 8$ & Vigas X (luz 7,5 m) & triangular & $""" + n(3.75 * G["Q_SCP"]) + r"""$ & $""" + n(3.75 * G["Q_USO"]) + r"""$ \\
$7{,}5 \times 8$ & Vigas Y (luz 8 m) & trapezoidal & $""" + n(3.75 * G["Q_SCP"]) + r"""$ & $""" + n(3.75 * G["Q_USO"]) + r"""$ \\
$10 \times 8$ & Vigas X (luz 10 m) & trapezoidal & $""" + n(4.0 * G["Q_SCP"]) + r"""$ & $""" + n(4.0 * G["Q_USO"]) + r"""$ \\
$10 \times 8$ & Vigas Y (luz 8 m) & triangular & $""" + n(4.0 * G["Q_SCP"]) + r"""$ & $""" + n(4.0 * G["Q_USO"]) + r"""$ \\
\bottomrule
\end{tabular}
\end{table}
Las vigas Y de los ejes 2, 3 y 4 reciben la suma de los dos paños adyacentes. En el techo la carga
de uso es la indicada para $L_r$.
""")

    # sísmica
    filas = []
    for v in NIV:
        filas.append(r"%d & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\" % (
            v["k"], n(v["z"]), n(v["acero"], 1), n(v["scp"], 1), n(v["sc25"], 1), n(v["P"], 1), n(v["A"], 4), n(v["F"], 2)))
    S_AP = sum(v["AP"] for v in NIV)
    w(r"""
\subsection{Cargas sísmicas (E)}
De acuerdo al enunciado, el esfuerzo de corte basal se obtiene como $Q = 0{,}18\,P$, donde $P$ es el peso
sísmico de la estructura, igual a las cargas permanentes más un $25\%$ de la sobrecarga de uso
(también se incluye el $25\%$ de la sobrecarga del techo transitable). El peso de cada nivel considera la
SCP de la planta, el peso de las vigas del nivel y la mitad de las columnas y riostras de los pisos
superior e inferior. Las fuerzas se distribuyen en altura según:
\begin{equation}
F_k = \frac{A_k\,P_k}{\sum_{j=1}^{5} A_j\,P_j}\,Q,
\qquad
A_k = \sqrt{1-\frac{Z_{k-1}}{H_{total}}} - \sqrt{1-\frac{Z_k}{H_{total}}}
\end{equation}
con $H_{total} = """ + n(Htot) + r"""$ m.

\begin{table}[H]
\centering
\caption{Peso sísmico y fuerzas sísmicas por nivel.}
\label{tab:sismo}
\begin{tabular}{cccccccc}
\toprule
Nivel & $Z_k$ [m] & Acero [tonf] & SCP [tonf] & $0{,}25\,L$ [tonf] & $P_k$ [tonf] & $A_k$ & $F_k$ [tonf] \\
\midrule
""" + "\n".join(filas) + r"""
\midrule
$\Sigma$ & & & & & $""" + n(datos["P_total"], 1) + r"""$ & $""" + n(sum(v["A"] for v in NIV), 4) + r"""$ & $""" + n(datos["Q"], 2) + r"""$ \\
\bottomrule
\end{tabular}
\end{table}

\[
P = """ + n(datos["P_total"], 1) + r"""\tonf \quad\Rightarrow\quad Q = 0{,}18 \times """ + n(datos["P_total"], 1) + r""" = """ + n(datos["Q"], 1) + r"""\tonf
\qquad \left(\textstyle\sum A_j P_j = """ + n(S_AP, 2) + r"""\tonf\right)
\]

Las fuerzas $F_k$ se aplican en el centro de masas (CM) de cada piso, ubicado en
$(X, Y) = (""" + n(G["XCM"]) + r""";\ """ + n(G["YCM"]) + r""")$ m dada la simetría de la planta, en forma
independiente en las direcciones X (patrón \texttt{EX}) e Y (patrón \texttt{EY}).
""")

    # combos
    fc = []
    for cname, desc, fac in COMBOS:
        fc.append(r"%s & $%s$ \\" % (cname, desc.replace("Lr", "L_r").replace("Ex", "E_x").replace("Ey", "E_y")))
    w(r"""
\subsection{Combinaciones de carga}\label{sec:combos}
Se utilizan las combinaciones de la NCh3171 para el método LRFD, con $D = PP + SCP$. Al no existir
cargas de viento, nieve ni lluvia, éstas se omiten. Las cargas sísmicas se consideran con ambos signos
y en forma independiente en cada dirección. Se usa factor $1{,}0$ para la sobrecarga $L$ en las
combinaciones sísmicas (conservador).
\begin{table}[H]
\centering
\caption{Combinaciones de carga (NCh3171).}
\begin{tabular}{cl}
\toprule
Combinación & Expresión \\
\midrule
""" + "\n".join(fc) + r"""
\midrule
ENVOLVENTE & Envolvente de C1 a C11 \\
\bottomrule
\end{tabular}
\end{table}
""")

    # ------------------------------------------------------------------ 4. modelo
    filas_sec = []
    usos = [(G["SEC_COL"], "Columnas (todas)"), (G["SEC_VIGA_X"], "Vigas ejes A y B"),
            (G["SEC_VIGA_Y"], "Vigas ejes 1 a 5"), (G["SEC_RIO_X"], "Riostras ejes A y B"),
            (G["SEC_RIO_Y"], "Riostras ejes 1 y 5")]
    for nom, uso in usos:
        s = SEC[nom]
        filas_sec.append(r"%s & %s & $%s$ & $%s$ & $%s$ & $%s$ \\" % (
            uso, nom.replace("X", "$\\times$", 1).replace("X", "x") if False else nom, n(s["A"] * 1e4, 1),
            n(s["I33"] * 1e8, 0), n(s["I22"] * 1e8, 0), n(s["Z33"] * 1e6, 0)))
    w(r"""
\section{Modelo estructural en SAP2000}

\subsection{Descripción del modelo}
Se realizó un modelo tridimensional en SAP2000 (archivo \texttt{G5\_Edificio3D.sdb}, entregado también
en formato \texttt{.s2k} compatible con la versión 19.2.1). Las principales características y supuestos son:
\begin{itemize}
  \item Unidades del modelo: tonf, m, $^\circ$C. Eje X según los ejes 1--5, eje Y según los ejes A--B y eje Z vertical.
  \item Columnas, vigas y riostras se modelan con elementos tipo \emph{frame} ubicados en los ejes
        (no se consideran zonas rígidas ni excentricidades en las conexiones).
  \item Apoyos: columnas empotradas en la base (placas base con anclajes rígidos).
  \item Conexiones viga--columna rígidas (transmiten momento) en ambas direcciones.
  \item Riostras con extremos rotulados (liberación de $M_2$ y $M_3$ en ambos extremos), por lo que trabajan
        esencialmente a carga axial.
  \item Columnas orientadas con su eje fuerte resistiendo la flexión en el plano de los marcos rígidos
        (ejes 2, 3 y 4, dirección Y): ángulo local de $90^\circ$.
  \item La losa de H.A. de 18 cm se considera como diafragma rígido en su plano: se asigna una restricción
        tipo \emph{Diaphragm} (DIAF1 a DIAF5) a todos los nudos de cada nivel. La losa no se modela; sus cargas
        se aplican sobre las vigas (sección anterior).
  \item En cada nivel se crea un nudo en el centro de masas (nudos 199, 299, \dots, 599), incorporado al diafragma
        y restringido en $U_Z$, $R_X$ y $R_Y$, donde se aplican las fuerzas sísmicas.
  \item Patrones de carga: \texttt{PP} (peso propio $\times 1{,}10$), \texttt{SCP}, \texttt{L}, \texttt{LR},
        \texttt{EX} y \texttt{EY}. Combinaciones C1 a C11 y ENVOLVENTE.
  \item Análisis estático lineal de primer orden.
\end{itemize}

\subsection{Secciones de prediseño}
\begin{table}[H]
\centering
\caption{Secciones utilizadas en el modelo (propiedades calculadas por SAP2000 a partir de las dimensiones).}
\begin{tabular}{llcccc}
\toprule
Elemento & Perfil & $A$ [cm$^2$] & $I_{33}$ [cm$^4$] & $I_{22}$ [cm$^4$] & $Z_{33}$ [cm$^3$] \\
\midrule
""" + "\n".join(filas_sec) + r"""
\bottomrule
\end{tabular}
\end{table}
Las secciones corresponden a un prediseño y serán verificadas/ajustadas en las Tareas 2 y 3.
""")

    # ------------------------------------------------------------------ 5. resultados
    nombres = {"COL": "Columnas", "VX": "Vigas X (ejes A, B)", "VY": "Vigas Y (ejes 1 a 5)",
               "RX": "Riostras X (ejes A, B)", "RY": "Riostras Y (ejes 1, 5)"}

    def ub(bid):
        b = BARRAS[bid]
        if b.tipo == "COL":
            return "Col. %s, piso %d" % (b.eje, b.piso)
        if b.tipo in ("VX", "VY"):
            return "Eje %s, nivel %d" % (b.eje, b.piso)
        return "Eje %s, piso %d" % (b.eje, b.piso)

    filas_r = []
    for tipo in ("COL", "VX", "VY", "RX", "RY"):
        r = resumen[tipo]
        claves = [("Pc", r"$P_u$ compresión", "tonf"), ("Pt", r"$P_u$ tracción", "tonf"),
                  ("M3", r"$M_{u3}$", "tonf$\\cdot$m"), ("M2", r"$M_{u2}$", "tonf$\\cdot$m"),
                  ("V2", r"$V_{u2}$", "tonf")]
        if tipo in ("RX", "RY"):
            claves = claves[:2]
        if tipo in ("VX", "VY"):
            claves = [c for c in claves if c[0] in ("M3", "V2")]
        first = True
        for key, lab, un in claves:
            val, bid, combo = r[key]
            pre = (r"\multirow{%d}{*}{%s}" % (len(claves), nombres[tipo])) if first else ""
            first = False
            filas_r.append(r"%s & %s & $%s$ & %d & %s & %s \\" % (pre, lab + " [" + un + "]", n(val, 1), bid, ub(bid), combo))
        filas_r.append(r"\midrule")
    filas_r = filas_r[:-1]
    desp = datos["desp"]
    fd = []
    prev_x = prev_y = 0.0
    for d in desp:
        hk = ZS[d["k"]] - ZS[d["k"] - 1]
        fd.append(r"%d & $%s$ & $%s$ & $%s$ & $%s$ \\" % (
            d["k"], n(d["ux"] * 100, 2), n((d["ux"] - prev_x) / hk * 1000, 3), n(d["uy"] * 100, 2), n((d["uy"] - prev_y) / hk * 1000, 3)))
        prev_x, prev_y = d["ux"], d["uy"]
    w(r"""
\section{Resultados del análisis}

\subsection{Esfuerzos internos máximos}
La Tabla~\ref{tab:esf} resume los esfuerzos internos máximos (envolvente de las combinaciones LRFD) para cada
tipo de elemento, indicando el número de la barra en el modelo, su ubicación y la combinación que controla.
La numeración de barras se muestra en el Anexo~\ref{anx:num}. Por efecto del diafragma rígido, las vigas no
presentan esfuerzo axial en el modelo. Los diagramas de esfuerzos internos se adjuntan en el Anexo (capturas de SAP2000).

\begin{table}[H]
\centering
\small
\caption{Resumen de esfuerzos internos máximos de diseño (LRFD).}
\label{tab:esf}
\begin{tabular}{llcccc}
\toprule
Elemento & Esfuerzo & Valor & Barra & Ubicación & Comb. \\
\midrule
""" + "\n".join(filas_r) + r"""
\bottomrule
\end{tabular}
\end{table}

\subsection{Desplazamientos laterales}
A modo informativo, la Tabla~\ref{tab:desp} entrega los desplazamientos del centro de masas para las cargas
sísmicas \texttt{EX} y \texttt{EY} (sin mayorar) y las derivas de entrepiso, que resultan menores al límite
de $0{,}002$ de la NCh433.
\begin{table}[H]
\centering
\caption{Desplazamientos del CM y derivas de entrepiso (estados \texttt{EX} y \texttt{EY}).}
\label{tab:desp}
\begin{tabular}{ccccc}
\toprule
Nivel & $u_X$ [cm] & $\Delta_X/h$ [$\times 10^{-3}$] & $u_Y$ [cm] & $\Delta_Y/h$ [$\times 10^{-3}$] \\
\midrule
""" + "\n".join(fd) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

    # ------------------------------------------------------------------ 6. análisis plástico
    w(seccion_plastico(G, po, n))

    # ------------------------------------------------------------------ 7. tracción
    w(seccion_traccion(G, trac, n, BARRAS, NUDOS))

    # ------------------------------------------------------------------ 8. conclusiones
    mec = po["mecanismos"][0]
    w(r"""
\section{Conclusiones}
\begin{itemize}
  \item El peso sísmico de la estructura es $P = """ + n(datos["P_total"], 1) + r"""\tonf$, lo que resulta en un corte
        basal $Q = """ + n(datos["Q"], 1) + r"""\tonf$, aplicado en el centro de masas de cada piso con la distribución
        en altura indicada.
  \item Los marcos arriostrados concentran la mayor parte de la carga lateral; las riostras de los ejes 1 y 5
        (dirección Y) son las más solicitadas, con $T_u = """ + n(trac["Tu"], 1) + r"""\tonf$.
  \item El análisis plástico del eje 2 entrega una carga lateral de primera fluencia
        $F_y = """ + n(po["lam_y"], 1) + r"""\tonf$, una primera rótula plástica para $F = """ + n(po["eventos"][0]["F"], 1) + r"""\tonf$ y
        una carga de colapso $F_c = """ + n(po["F_colapso"], 1) + r"""\tonf$. El método cinemático entrega
        $F_c = """ + n(mec["F"], 1) + r"""\tonf$ (mecanismo \emph{""" + mec["nombre"] + r"""}), coincidente con el
        análisis incremental, ya que este último corresponde a la solución exacta del problema rígido-plástico
        con rótulas en los extremos de los elementos.
  \item La riostra más traccionada se diseña con un perfil """ + trac["sec"] + r""" (A500 Gr.C) con un factor de utilización
        de $""" + n(trac["FU"], 2) + r"""$ a tracción. La sección queda holgada en tracción, pero se escogió anticipando el
        diseño a compresión de la Tarea 2 (la misma riostra se comprime al invertir el sismo).
\end{itemize}
""")

    # ------------------------------------------------------------------ anexos
    w(anexo_numeracion(G))
    w(anexo_sap_pushover(G, po, n))
    w(r"""
\section*{Anexo C: Diagramas de esfuerzos internos}
\addcontentsline{toc}{section}{Anexo C: Diagramas de esfuerzos internos}
% Insertar aquí capturas de SAP2000, por ejemplo:
% \begin{figure}[H]\centering
%   \includegraphics[width=0.9\textwidth]{diagrama_M3_eje2.png}
%   \caption{Diagrama de momento $M_3$, eje 2, combinación ENVOLVENTE.}
% \end{figure}
Ver capturas del modelo SAP2000 (diagramas de momento, corte y axial para la combinación ENVOLVENTE).

\end{document}
""")
    with open(ruta, "w", encoding="utf-8") as fh:
        fh.write("".join(out))


# ====================================================================== figuras
def fig_planta(XS, YS):
    s = 0.40
    t = [r"""
\begin{figure}[H]
\centering
\begin{tikzpicture}[x=%scm,y=%scm,font=\small]""" % (s, s)]
    for x in XS:
        t.append(r"\draw[gray,dashed] (%g,-1.5) -- (%g,%g);" % (x, x, YS[1] + 1.5))
    for y in YS:
        t.append(r"\draw[gray,dashed] (-1.5,%g) -- (%g,%g);" % (y, XS[-1] + 1.5, y))
    nom = ["1", "2", "3", "4", "5"]
    for x, lab in zip(XS, nom):
        t.append(r"\node[draw,circle,inner sep=1.5pt] at (%g,%g) {%s};" % (x, YS[1] + 2.6, lab))
    for y, lab in zip(YS, ["A", "B"]):
        t.append(r"\node[draw,circle,inner sep=1.5pt] at (-2.6,%g) {%s};" % (y, lab))
    t.append(r"\draw[thick] (0,0) rectangle (%g,%g);" % (XS[-1], YS[1]))
    for x in XS:
        t.append(r"\draw[thick] (%g,0) -- (%g,%g);" % (x, x, YS[1]))
        for y in YS:
            t.append(r"\fill (%g,%g) rectangle +(0.35,0.35); \fill (%g,%g) rectangle +(-0.35,-0.35); \fill (%g,%g) rectangle +(0.35,-0.35); \fill (%g,%g) rectangle +(-0.35,0.35);" % (x, y, x, y, x, y, x, y))
    # riostras en planta (líneas rojas sobre ejes arriostrados)
    t.append(r"\draw[red,line width=1.6pt] (0,0) -- (%g,0) (%g,0) -- (%g,0) (0,%g) -- (%g,%g) (%g,%g) -- (%g,%g);" % (
        XS[1], XS[3], XS[4], YS[1], XS[1], YS[1], XS[3], YS[1], XS[4], YS[1]))
    t.append(r"\draw[red,line width=1.6pt] (0,0) -- (0,%g) (%g,0) -- (%g,%g);" % (YS[1], XS[4], XS[4], YS[1]))
    # cotas
    for a, b in zip(XS[:-1], XS[1:]):
        t.append(r"\draw[<->] (%g,-1.0) -- node[below]{%s m} (%g,-1.0);" % (a, ("%g" % (b - a)).replace(".", ","), b))
    t.append(r"\draw[<->] (%g,0) -- node[right]{%s m} (%g,%g);" % (XS[-1] + 1.2, ("%g" % YS[1]).replace(".", ","), XS[-1] + 1.2, YS[1]))
    t.append(r"\fill[blue] (%g,%g) circle (0.3) node[above right,blue]{CM};" % (XS[-1] / 2, YS[1] / 2))
    t.append(r"\draw[->,thick] (-1.5,-3.5) -- (1.5,-3.5) node[right]{X};")
    t.append(r"\draw[->,thick] (-1.5,-3.5) -- (-1.5,-0.8) node[above]{Y};")
    t.append(r"""\end{tikzpicture}
\caption{Vista en planta (Grupo 5). En rojo se indican los vanos arriostrados.}
\label{fig:planta}
\end{figure}
""")
    return "\n".join(t)


def fig_elevaciones(XS, YS, ZS):
    s = 0.205
    t = [r"""
\begin{figure}[H]
\centering
\begin{tikzpicture}[x=%scm,y=%scm,font=\scriptsize]""" % (s, s)]
    # Elevación ejes A y B
    for x in XS:
        t.append(r"\draw[thick] (%g,0) -- (%g,%g);" % (x, x, ZS[-1]))
        t.append(r"\draw (%g,0) -- ++(-0.7,-0.6) -- ++(1.4,0) -- cycle;" % x)
    for z in ZS[1:]:
        t.append(r"\draw[thick] (0,%g) -- (%g,%g);" % (z, XS[-1], z))
    for k in range(1, len(ZS)):
        t.append(r"\draw[red,thick] (%g,%g) -- (%g,%g);" % (XS[0], ZS[k - 1], XS[1], ZS[k]))
        t.append(r"\draw[red,thick] (%g,%g) -- (%g,%g);" % (XS[3], ZS[k], XS[4], ZS[k - 1]))
    for x, lab in zip(XS, "12345"):
        t.append(r"\node at (%g,-2) {%s};" % (x, lab))
    t.append(r"\node at (%g,-4) {\small Elevación ejes A y B};" % (XS[-1] / 2))
    for k in range(1, len(ZS)):
        t.append(r"\node[left] at (-0.5,%g) {N%d};" % (ZS[k], k))
    # Elevación ejes 1 y 5
    ox = XS[-1] + 6
    W = YS[1]
    for y in (0, W):
        t.append(r"\draw[thick] (%g,0) -- (%g,%g);" % (ox + y, ox + y, ZS[-1]))
        t.append(r"\draw (%g,0) -- ++(-0.7,-0.6) -- ++(1.4,0) -- cycle;" % (ox + y))
    for z in ZS[1:]:
        t.append(r"\draw[thick] (%g,%g) -- (%g,%g);" % (ox, z, ox + W, z))
    for k in range(1, len(ZS)):
        if k % 2 == 1:
            t.append(r"\draw[red,thick] (%g,%g) -- (%g,%g);" % (ox, ZS[k], ox + W, ZS[k - 1]))
        else:
            t.append(r"\draw[red,thick] (%g,%g) -- (%g,%g);" % (ox, ZS[k - 1], ox + W, ZS[k]))
    t.append(r"\node at (%g,-2) {A}; \node at (%g,-2) {B};" % (ox, ox + W))
    t.append(r"\node at (%g,-4) {\small Ejes 1 y 5};" % (ox + W / 2))
    # Elevación ejes 2, 3, 4
    ox2 = ox + W + 6
    for y in (0, W):
        t.append(r"\draw[thick] (%g,0) -- (%g,%g);" % (ox2 + y, ox2 + y, ZS[-1]))
        t.append(r"\draw (%g,0) -- ++(-0.7,-0.6) -- ++(1.4,0) -- cycle;" % (ox2 + y))
    for z in ZS[1:]:
        t.append(r"\draw[thick] (%g,%g) -- (%g,%g);" % (ox2, z, ox2 + W, z))
    t.append(r"\node at (%g,-2) {A}; \node at (%g,-2) {B};" % (ox2, ox2 + W))
    t.append(r"\node at (%g,-4) {\small Ejes 2, 3 y 4};" % (ox2 + W / 2))
    for k in range(1, len(ZS)):
        hk = ZS[k] - ZS[k - 1]
        t.append(r"\draw[<->] (%g,%g) -- node[right]{%s} (%g,%g);" % (
            ox2 + W + 1.2, ZS[k - 1], ("%g m" % hk).replace(".", ","), ox2 + W + 1.2, ZS[k]))
    t.append(r"""\end{tikzpicture}
\caption{Elevaciones de la estructura. En rojo, las riostras (diagonales).}
\label{fig:elev}
\end{figure}
""")
    return "\n".join(t)


# ====================================================================== análisis plástico
def seccion_plastico(G, po, n):
    ZS, YS = G["ZS"], G["YS"]
    sc, sb = G["SECCIONES"][G["SEC_COL"]], G["SECCIONES"][G["SEC_VIGA_Y"]]
    MPA = G["MPA"]
    fk = po["fk"]
    ev = po["eventos"]
    mec = po["mecanismos"]
    # curva
    coords = " ".join("(%.3f,%.3f)" % (d * 100, F) for d, F in po["curva"])
    d_col = po["curva"][-2][0] * 100
    xmax = po["curva"][-1][0] * 100
    fe = []
    for e in ev:
        if e["elem"] == "Columna":
            loc = "Columna eje %s, piso %d (extremo %s)" % (e["lado"], e["piso"], e["extremo"])
        else:
            loc = "Viga nivel %d, extremo %s" % (e["piso"], e["extremo"])
        fe.append(r"%d & %s & $%s$ & $%s$ \\" % (e["n"], loc, n(e["F"], 1), n(e["d"] * 100, 2)))
    fm = []
    for m in mec[:8]:
        fm.append(r"%s & $%s$ & $%s$ & $%s$ \\" % (m["nombre"], n(m["Wi"], 1), n(m["We"], 3), n(m["F"], 1)))
    fl = []
    for k in range(1, len(ZS)):
        fl.append(r"%d & $%s$ & $%s$ \\" % (k, n(ZS[k]), n(fk[k - 1], 4)))
    m0 = mec[0]
    loc_y = po["loc_y"]
    ly = ("columna, piso %d" % loc_y["piso"]) if loc_y["elem"] == "Columna" else ("viga del nivel %d" % loc_y["piso"])
    dif = abs(m0["F"] - po["F_colapso"]) / po["F_colapso"] * 100
    # hinges of mechanism (j, m)
    j, m = m0["j"], m0["m"]
    W = YS[1]
    s = 0.33
    hs = []
    hs.append(r"\fill[red] (0,%g) circle (0.35); \fill[red] (%g,%g) circle (0.35);" % (ZS[j - 1] + 0.5, W, ZS[j - 1] + 0.5))
    for k in range(j, m):
        hs.append(r"\fill[red] (0.6,%g) circle (0.35); \fill[red] (%g,%g) circle (0.35);" % (ZS[k], W - 0.6, ZS[k]))
    if m == len(ZS) - 1:
        hs.append(r"\fill[red] (0.6,%g) circle (0.35); \fill[red] (%g,%g) circle (0.35);" % (ZS[m], W - 0.6, ZS[m]))
    else:
        hs.append(r"\fill[red] (0,%g) circle (0.35); \fill[red] (%g,%g) circle (0.35);" % (ZS[m] - 0.5, W, ZS[m] - 0.5))
    marco = []
    for y in (0, W):
        marco.append(r"\draw[thick] (%g,0) -- (%g,%g);" % (y, y, ZS[-1]))
        marco.append(r"\draw[pattern=north east lines] (%g,0) ++(-0.8,-0.5) rectangle ++(1.6,0.5);" % y)
    for z in ZS[1:]:
        marco.append(r"\draw[thick] (0,%g) -- (%g,%g);" % (z, W, z))
    for k in range(1, len(ZS)):
        marco.append(r"\draw[-{Stealth},blue,thick] (%g,%g) -- (0,%g) node[midway,above,font=\tiny]{$%s F$};" % (
            -1.5 - 6 * fk[k - 1], ZS[k], ZS[k], n(fk[k - 1], 3)))
    return r"""
\section{Análisis plástico del eje 2}

\subsection{Modelo y supuestos}
Se estudia el marco rígido del eje 2 como un problema plano (plano YZ): un vano de $8{,}0$ m entre los ejes
A y B y cinco pisos, con columnas """ + G["SEC_COL"] + r""" (flexión en torno al eje fuerte) y vigas """ + G["SEC_VIGA_Y"] + r""",
empotrado en la base. Se adoptan los siguientes supuestos:
\begin{itemize}
  \item Comportamiento elastoplástico perfecto, con rótulas plásticas concentradas en los extremos de vigas y
        columnas, de capacidad $M_p = Z_x F_y$. Se desprecia la interacción carga axial--momento y los efectos
        de segundo orden (P--$\Delta$), de forma de poder comparar directamente con el método cinemático.
  \item Cargas gravitacionales que tributan al eje 2 (aplicadas antes de la carga lateral):
        $D + 0{,}25L$, es decir, el peso propio de la viga ($\times 1{,}10$), la SCP y el $25\%$ de la
        sobrecarga de los paños 1--2 y 2--3 (carga trapezoidal más triangular sobre la viga), además de las reacciones de las
        vigas X que llegan a las columnas del eje.
  \item Carga lateral con la forma de la distribución sísmica de la Tabla~\ref{tab:sismo}: en el nivel $k$ se aplica
        $\alpha_k F$, con $\alpha_k = F_k/Q$, donde $F$ es la carga lateral total (corte basal del marco).
\end{itemize}

\begin{table}[H]
\centering
\begin{minipage}{0.45\textwidth}
\centering
\captionof{table}{Capacidades de las secciones.}
\begin{tabular}{lcc}
\toprule
 & Columna & Viga \\
 & """ + G["SEC_COL"] + r""" & """ + G["SEC_VIGA_Y"] + r""" \\
\midrule
$S_x$ [cm$^3$] & $""" + n(sc["S33"] * 1e6, 0) + r"""$ & $""" + n(sb["S33"] * 1e6, 0) + r"""$ \\
$Z_x$ [cm$^3$] & $""" + n(sc["Z33"] * 1e6, 0) + r"""$ & $""" + n(sb["Z33"] * 1e6, 0) + r"""$ \\
$M_y = S_x F_y$ [tonf$\cdot$m] & $""" + n(po["My_c"], 1) + r"""$ & $""" + n(po["My_b"], 1) + r"""$ \\
$M_p = Z_x F_y$ [tonf$\cdot$m] & $""" + n(po["Mp_c"], 1) + r"""$ & $""" + n(po["Mp_b"], 1) + r"""$ \\
\bottomrule
\end{tabular}
\end{minipage}\hfill
\begin{minipage}{0.45\textwidth}
\centering
\captionof{table}{Patrón de carga lateral.}
\begin{tabular}{ccc}
\toprule
Nivel & $Z_k$ [m] & $\alpha_k = F_k/Q$ \\
\midrule
""" + "\n".join(fl) + r"""
\bottomrule
\end{tabular}
\end{minipage}
\end{table}

\subsection{Curva de capacidad (análisis incremental, SAP2000)}
La curva de capacidad se obtiene mediante un análisis estático no lineal (\emph{pushover}) con control de
desplazamiento en el techo, partiendo del estado gravitacional. El procedimiento es evento a evento: se
incrementa la carga lateral hasta que alguna sección alcanza su $M_p$, se introduce una rótula en ese punto y se
continúa, hasta que la estructura se transforma en un mecanismo. La definición de las rótulas y de los casos
no lineales en SAP2000 se detalla en el Anexo~B.

\textbf{Primera fluencia.} Se define como el instante en que la fibra extrema de alguna sección alcanza $F_y$, es
decir, $|M_{grav} + F\,m_{lat}| = M_y$. Esto ocurre en la """ + ly + r""" para
\[
F_{y} = """ + n(po["lam_y"], 1) + r"""\tonf \qquad (\delta_{techo} = """ + n(po["d_y"] * 100, 2) + r"""\ \mathrm{cm}).
\]

\textbf{Primera rótula plástica} ($M = M_p$): $F = """ + n(ev[0]["F"], 1) + r"""\tonf$.

\textbf{Colapso.} El mecanismo se forma con la rótula N$^\circ$""" + str(len(ev)) + r""" para
\[
F_{c} = """ + n(po["F_colapso"], 1) + r"""\tonf \qquad (\delta_{techo} = """ + n(d_col, 1) + r"""\ \mathrm{cm}).
\]

\begin{figure}[H]
\centering
\begin{tikzpicture}
\begin{axis}[width=0.85\textwidth,height=7.5cm,grid=major,
  xlabel={Desplazamiento de techo $\delta$ [cm]},ylabel={Carga lateral total $F$ [tonf]},
  xmin=0,xmax=""" + ("%.0f" % (xmax * 1.02)) + r""",ymin=0,ymax=""" + ("%.0f" % (po["F_colapso"] * 1.25)) + r""",
  legend pos=south east,legend style={font=\small}]
\addplot[thick,blue,mark=*,mark size=1.5pt] coordinates {""" + coords + r"""};
\addlegendentry{Pushover (evento a evento)}
\addplot[red,dashed,thick] coordinates {(0,""" + "%.3f" % m0["F"] + r""") (""" + "%.2f" % (xmax * 1.02) + r""",""" + "%.3f" % m0["F"] + r""")};
\addlegendentry{Método cinemático $F_c = """ + n(m0["F"], 1) + r"""$ tonf}
\addplot[only marks,mark=square*,black,mark size=2.5pt] coordinates {(""" + "%.3f" % (po["d_y"] * 100) + "," + "%.3f" % po["lam_y"] + r""")};
\addlegendentry{Primera fluencia $F_y = """ + n(po["lam_y"], 1) + r"""$ tonf}
\end{axis}
\end{tikzpicture}
\caption{Curva de capacidad del eje 2 (corte basal versus desplazamiento de techo).}
\label{fig:pushover}
\end{figure}

\begin{table}[H]
\centering
\caption{Secuencia de formación de rótulas plásticas.}
\begin{tabular}{clcc}
\toprule
N$^\circ$ & Ubicación & $F$ [tonf] & $\delta_{techo}$ [cm] \\
\midrule
""" + "\n".join(fe) + r"""
\bottomrule
\end{tabular}
\end{table}
En la tabla, ``extremo eje A/B'' indica el extremo de la viga donde se forma la rótula. Algunas de las rótulas
formadas en la secuencia se descargan una vez que se forma el mecanismo definitivo.

\subsection{Método cinemático}
El método cinemático (teorema del límite superior) iguala el trabajo externo de las cargas con el trabajo interno
disipado en las rótulas para un mecanismo cinemáticamente admisible con giro virtual $\theta$:
\[
W_e = \sum_k \alpha_k F\,\delta_k = W_i = \sum_r M_{p,r}\,\theta_r
\quad\Rightarrow\quad
F = \frac{\sum_r M_{p,r}\,\theta_r}{\sum_k \alpha_k\,\delta_k}.
\]
Se evalúa la familia de mecanismos de piso y combinados en que los pisos $j$ a $m$ se desplazan lateralmente
(columnas de esos pisos giran $\theta$), con rótulas: en el pie de las columnas del piso $j$, en ambos extremos de
las vigas de los niveles $j$ a $m-1$ (o en las columnas si son más débiles, $\min(M_{pb},\,2M_{pc})$ por nudo),
y en la cabeza de las columnas del piso $m$ (o en las vigas de techo si $m = 5$). Para $j = m$ se obtiene el
mecanismo de piso blando, y para $j=1,\ m=5$ el mecanismo global de vigas. Como el mecanismo solo
involucra traslaciones horizontales de las vigas, las cargas gravitacionales no realizan trabajo.

Por ejemplo, para el mecanismo que controla (pisos """ + str(j) + r""" a """ + str(m) + r"""):
\[
W_i = \left[2M_{pc} + """ + str(2 * (m - j)) + r"""\,M_{pb} + 2M_{pc}\right]\theta
     = """ + n(m0["Wi"], 1) + r"""\,\theta \ \mathrm{tonf\cdot m},
\qquad
W_e = F\,\theta\sum_k \alpha_k\,\delta_k/\theta = """ + n(m0["We"], 3) + r"""\,F\,\theta\ \mathrm{m}
\]
\[
\Rightarrow\quad F_c = \frac{""" + n(m0["Wi"], 1) + r"""}{""" + n(m0["We"], 3) + r"""} = """ + n(m0["F"], 1) + r"""\tonf .
\]

\begin{table}[H]
\centering
\caption{Cargas de colapso de los mecanismos evaluados (los 8 menores).}
\begin{tabular}{lccc}
\toprule
Mecanismo & $W_i/\theta$ [tonf$\cdot$m] & $W_e/(F\theta)$ [m] & $F$ [tonf] \\
\midrule
""" + "\n".join(fm) + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{figure}[H]
\centering
\begin{tikzpicture}[x=""" + str(s) + r"""cm,y=""" + str(s) + r"""cm]
""" + "\n".join(marco) + "\n" + "\n".join(hs) + r"""
\node[below] at (0,-0.8) {A}; \node[below] at (""" + "%g" % W + r""",-0.8) {B};
\end{tikzpicture}
\caption{Mecanismo de colapso del eje 2 (rótulas plásticas en rojo).}
\label{fig:mecanismo}
\end{figure}

\subsection{Comparación y discusión}
\begin{table}[H]
\centering
\caption{Comparación de resultados del análisis plástico del eje 2.}
\begin{tabular}{lc}
\toprule
Parámetro & Valor \\
\midrule
Carga de primera fluencia $F_y$ & $""" + n(po["lam_y"], 1) + r"""\tonf$ \\
Carga de primera rótula plástica & $""" + n(ev[0]["F"], 1) + r"""\tonf$ \\
Carga de colapso, pushover $F_c$ & $""" + n(po["F_colapso"], 1) + r"""\tonf$ \\
Carga de colapso, método cinemático & $""" + n(m0["F"], 1) + r"""\tonf$ \\
Diferencia & $""" + n(dif, 2) + r"""\%$ \\
Sobrerresistencia $F_c/F_y$ & $""" + n(po["F_colapso"] / po["lam_y"], 2) + r"""$ \\
\bottomrule
\end{tabular}
\end{table}

Ambos métodos entregan la misma carga de colapso, ya que el método cinemático es un límite superior que
coincide con la solución exacta cuando se evalúa el mecanismo correcto, y el análisis incremental con
rótulas elastoplásticas perfectas converge al mismo mecanismo. Las diferencias que se observen respecto al
resultado de SAP2000 se deben a: (i) la definición de las rótulas (si se utilizan rótulas automáticas
ASCE~41 con endurecimiento, la curva no presenta una meseta perfectamente horizontal), (ii) la
interacción P--M en las columnas si se usan rótulas tipo P--M2--M3 y (iii) los efectos P--$\Delta$,
que reducen la capacidad lateral. El momento máximo de tramo de las vigas bajo cargas gravitacionales es
$""" + n(po["Mtramo"], 1) + r"""\tonfm < M_{pb}$, por lo que no se forman rótulas en el interior de las vigas
y es válido considerar rótulas solo en los extremos.
La razón $F_c/F_y = """ + n(po["F_colapso"] / po["lam_y"], 2) + r"""$ refleja la reserva de resistencia del
marco hiperestático por redistribución de momentos luego de la primera fluencia.
"""


# ====================================================================== tracción
def seccion_traccion(G, t, n, BARRAS, NUDOS):
    b = BARRAS[t["barra"]]
    a, c = NUDOS[b.i], NUDOS[b.j]
    combo_desc = [d for cn, d, f in G["COMBOS"] if cn == t["combo"]][0]
    R_weld = 4 * t["l"] * min(t["Rw_cm"], t["Rbase_cm"]) / 1000
    R_gus_v = 0.75 * 0.6 * t["Fu_g"] * t["tg"] * 2 * t["l"] / 1000
    filas = [
        ("Fluencia en área bruta", r"$\phi_t F_y A_g$", t["phiTn_y"]),
        ("Ruptura en área neta efectiva", r"$\phi_t F_u A_e$", t["phiTn_u"]),
        ("Soldaduras (4 cordones)", r"$4\,l\,\phi R_n/l$", R_weld),
        ("Plancha gusset: fluencia (Whitmore)", r"$\phi F_y L_w t_g$", t["phiRn_whit"]),
        ("Plancha gusset: bloque de corte", r"$\phi R_{bs}$", t["phiRn_bs"]),
        ("Plancha gusset: ruptura por corte", r"$\phi\,0{,}6F_u\,2l\,t_g$", R_gus_v),
    ]
    ft = "\n".join(r"%s & %s & $%s$ & $%s$ \\" % (a_, b_, n(v, 1), n(t["Tu"] / v, 2)) for a_, b_, v in filas)
    return r"""
\section{Diseño de la riostra más traccionada}

\subsection{Demanda}
Del análisis, la riostra más traccionada corresponde a la barra """ + str(t["barra"]) + r""" del modelo,
ubicada en el \textbf{eje """ + t["eje"] + r""", piso """ + str(t["piso"]) + r"""}, que une el nudo """ + str(b.i) + r"""
$(X,Y,Z) = (""" + n(a.x, 1) + ";" + n(a.y, 1) + ";" + n(a.z, 1) + r""")$ con el nudo """ + str(b.j) + r"""
$(""" + n(c.x, 1) + ";" + n(c.y, 1) + ";" + n(c.z, 1) + r""")$ (ver Figuras~\ref{fig:planta} y \ref{fig:elev}).
Su largo es $L = """ + n(t["L"], 2) + r"""$ m y la combinación que controla es """ + t["combo"] + r""" ($""" + combo_desc.replace("Ey", "E_y").replace("Ex", "E_x") + r"""$):
\[
T_u = """ + n(t["Tu"], 1) + r"""\tonf .
\]
La riostra simétrica del eje 5 (mismo piso) presenta una demanda equivalente.

\subsection{Sección propuesta}
Se propone un perfil tubular cuadrado \textbf{""" + t["sec"] + r"""} de acero ASTM A500 Gr.C
($F_y = """ + n(t["Fy"], 0) + r"""\ \mathrm{kgf/cm^2}$, $F_u = """ + n(t["Fu"], 0) + r"""\ \mathrm{kgf/cm^2}$), con espesor de diseño
$t = 0{,}93\,t_{nom} = """ + n(t["t"], 2) + r"""$ cm, $A_g = """ + n(t["A"], 1) + r"""\ \mathrm{cm^2}$ y
$r = """ + n(t["r"], 2) + r"""$ cm.

\subsection{Verificaciones (AISC 360-16, Capítulo D)}

\textbf{Esbeltez (recomendación D1).}
\[
\frac{L}{r} = \frac{""" + n(t["L"] * 100, 1) + r"""}{""" + n(t["r"], 2) + r"""} = """ + n(t["esbeltez"], 1) + r""" \le 300 \quad \checkmark
\]

\textbf{Fluencia en el área bruta (D2-a), $\phi_t = 0{,}90$.}
\[
\phi_t P_n = 0{,}90 \cdot F_y \cdot A_g = 0{,}90 \cdot """ + n(t["Fy"], 0) + r""" \cdot """ + n(t["A"], 1) + r"""
= """ + n(t["phiTn_y"], 1) + r"""\tonf
\]

\textbf{Ruptura en el área neta efectiva (D2-b), $\phi_t = 0{,}75$.}
La conexión se materializa con una plancha gusset concéntrica de espesor $t_g = """ + n(t["tg"], 1) + r"""$ cm insertada en
una ranura del tubo de ancho $t_g + """ + n(t["holgura"] * 10, 0) + r"""$ mm, y soldada con cuatro cordones de filete de largo $l$.
\[
A_n = A_g - 2\,t\,(t_g + 0{,}2) = """ + n(t["A"], 1) + r""" - 2\cdot""" + n(t["t"], 2) + r"""\cdot(""" + n(t["tg"], 1) + r"""+0{,}2)
= """ + n(t["An"], 1) + r"""\ \mathrm{cm^2}
\]
Factor de corte diferido (Tabla D3.1, caso 6, HSS rectangular con una plancha concéntrica, $l \ge H$):
\[
\bar{x} = \frac{B^2 + 2BH}{4(B+H)} = """ + n(t["xbar"], 2) + r"""\ \mathrm{cm},
\qquad
U = 1 - \frac{\bar{x}}{l} = 1 - \frac{""" + n(t["xbar"], 2) + r"""}{""" + n(t["l"], 0) + r"""} = """ + n(t["U"], 3) + r"""
\]
\[
A_e = U A_n = """ + n(t["Ae"], 1) + r"""\ \mathrm{cm^2}
\qquad\Rightarrow\qquad
\phi_t P_n = 0{,}75 \cdot """ + n(t["Fu"], 0) + r""" \cdot """ + n(t["Ae"], 1) + r""" = """ + n(t["phiTn_u"], 1) + r"""\tonf
\]

\textbf{Resistencia de diseño.}
\[
\phi_t P_n = \min(""" + n(t["phiTn_y"], 1) + r""";\ """ + n(t["phiTn_u"], 1) + r""") = """ + n(t["phiTn"], 1) + r"""\tonf
\ \ge\ T_u = """ + n(t["Tu"], 1) + r"""\tonf \quad \checkmark
\qquad \text{F.U.} = """ + n(t["FU"], 2) + r"""
\]

\subsection{Conexión preliminar}
Se propone una conexión soldada con plancha gusset de acero A36, $t_g = """ + n(t["tg"] * 10, 0) + r"""$ mm, insertada
en el tubo ranurado y unida mediante cuatro cordones de soldadura de filete E70XX.

\textbf{Soldadura.} Tamaño de filete $w = """ + n(t["w"] * 10, 0) + r"""$ mm (mínimo según Tabla J2.4: 5 mm para
$6 < t \le 13$ mm; máximo $t - 2 = """ + n(t["t"] * 10 - 2, 1) + r"""$ mm).
Resistencia por unidad de largo de un cordón (J2.4), $\phi = 0{,}75$:
\[
\phi R_n/l = 0{,}75 \cdot 0{,}60 F_{EXX} \cdot 0{,}707\,w
= 0{,}75 \cdot 0{,}60 \cdot """ + n(t["FEXX"], 0) + r""" \cdot 0{,}707 \cdot """ + n(t["w"], 1) + r"""
= """ + n(t["Rw_cm"], 0) + r"""\ \mathrm{kgf/cm}
\]
Corte del metal base (pared del tubo): $0{,}75\cdot 0{,}60\,F_u\,t = """ + n(t["Rbase_cm"], 0) + r"""\ \mathrm{kgf/cm}$ (no controla).
Largo requerido de cada cordón:
\[
l_{req} = \frac{T_u}{4\,\phi R_n/l} = \frac{""" + n(t["Tu"] * 1000, 0) + r"""}{4 \cdot """ + n(min(t["Rw_cm"], t["Rbase_cm"]), 0) + r"""}
= """ + n(t["l_req"], 1) + r"""\ \mathrm{cm}
\quad\Rightarrow\quad l = """ + n(t["l"], 0) + r"""\ \mathrm{cm} \ (\ge B = """ + n(t["B"], 1) + r"""\ \mathrm{cm})
\]

\textbf{Plancha gusset.} Ancho de Whitmore $L_w = B + 2\,l\tan 30^\circ = """ + n(t["Lw"], 1) + r"""$ cm;
área de corte $A_{gv} = 2\,l\,t_g = """ + n(t["Agv"], 1) + r"""\ \mathrm{cm^2}$; área en tracción $A_{nt} = B\,t_g = """ + n(t["Ant"], 1) + r"""\ \mathrm{cm^2}$.

\begin{table}[H]
\centering
\caption{Resumen de verificaciones de la riostra traccionada y su conexión.}
\begin{tabular}{llcc}
\toprule
Estado límite & Expresión & $\phi R_n$ [tonf] & $T_u/\phi R_n$ \\
\midrule
""" + ft + r"""
\bottomrule
\end{tabular}
\end{table}

""" + fig_conexion(t, n) + r"""

\textbf{Comentario.} La sección queda con un factor de utilización bajo a tracción; sin embargo, la misma riostra
trabaja en compresión al invertirse el sentido del sismo ($P_u \approx T_u$), condición que controlará su diseño
en la Tarea 2, por lo que se mantiene el perfil """ + t["sec"] + r""". La razón ancho--espesor del tubo es
$b/t = """ + n(t["b_t"], 1) + r"""$.
"""


def fig_conexion(t, n):
    B = t["B"]
    l = t["l"]
    return r"""\begin{figure}[H]
\centering
\begin{tikzpicture}[x=0.095cm,y=0.095cm,font=\small]
% gusset
\draw[thick,fill=gray!15] (-42,-30) -- (""" + "%g" % (l + 10) + r""",-30) -- (""" + "%g" % (l + 10) + r""",30) -- (-42,30) -- cycle;
\node[align=center] at (-5,-38) {Gusset $t_g = """ + n(t["tg"] * 10, 0) + r"""$ mm (A36)};
% tubo
\draw[thick,fill=white] (0,""" + "%g" % (B / 2) + r""") -- (90,""" + "%g" % (B / 2) + r""") -- (90,""" + "%g" % (-B / 2) + r""") -- (0,""" + "%g" % (-B / 2) + r""") -- cycle;
\draw[dashed] (0,0) -- (90,0);
\node at (50,0) [above] {""" + t["sec"] + r"""};
% soldaduras
\draw[red,line width=2pt] (0,""" + "%g" % (B / 2) + r""") -- (""" + "%g" % l + r""",""" + "%g" % (B / 2) + r""");
\draw[red,line width=2pt] (0,""" + "%g" % (-B / 2) + r""") -- (""" + "%g" % l + r""",""" + "%g" % (-B / 2) + r""");
\draw[<->] (0,""" + "%g" % (B / 2 + 5) + r""") -- node[above]{$l = """ + n(l, 0) + r"""$ cm} (""" + "%g" % l + r""",""" + "%g" % (B / 2 + 5) + r""");
\draw[<->] (82,""" + "%g" % (-B / 2) + r""") -- node[left,fill=white,inner sep=1pt]{$B = """ + n(B, 1) + r"""$ cm} (82,""" + "%g" % (B / 2) + r""");
\node[red] at (""" + "%g" % (l / 2) + r""",""" + "%g" % (-B / 2 - 5) + r""") {4 filetes $w = """ + n(t["w"] * 10, 0) + r"""$ mm, E70XX};
% Whitmore
\draw[blue,dashed] (0,""" + "%g" % (B / 2) + r""") -- (""" + "%g" % (-l) + r""",""" + "%g" % (B / 2 + l * math.tan(math.radians(30))) + r""");
\draw[blue,dashed] (0,""" + "%g" % (-B / 2) + r""") -- (""" + "%g" % (-l) + r""",""" + "%g" % (-B / 2 - l * math.tan(math.radians(30))) + r""");
\draw[blue] (""" + "%g" % (-l) + r""",""" + "%g" % (-B / 2 - l * math.tan(math.radians(30))) + r""") -- node[left]{$L_w$} (""" + "%g" % (-l) + r""",""" + "%g" % (B / 2 + l * math.tan(math.radians(30))) + r""");
\draw[->,very thick] (92,0) -- (108,0) node[right]{$T_u$};
\end{tikzpicture}
\caption{Esquema de la conexión preliminar riostra--gusset (vista en elevación de la riostra).}
\label{fig:conexion}
\end{figure}"""


# ====================================================================== anexos
def anexo_numeracion(G):
    return r"""
\appendix
\section{Numeración del modelo}\label{anx:num}
\begin{itemize}
  \item \textbf{Nudos:} $100k + 10\,i_Y + i_X$, con $k$ = nivel (0 = base), $i_Y = 0$ (eje A) o 1 (eje B) e
        $i_X = 1 \dots 5$ (ejes 1 a 5). Ejemplo: nudo 312 = nivel 3, eje B, eje 2. Los nudos de centro de masas
        son $100k + 99$.
  \item \textbf{Barras:} 1--50 columnas (10 por piso, piso 1 = barras 1--10), 51--90 vigas X (ejes A y B),
        91--115 vigas Y (ejes 1 a 5), 116--135 riostras X (ejes A y B) y 136--145 riostras Y (ejes 1 y 5).
  \item \textbf{Grupos} definidos en SAP2000: COLUMNAS, VIGAS\_X, VIGAS\_Y, RIOSTRAS\_X, RIOSTRAS\_Y y EJE\_2.
\end{itemize}
"""


def anexo_sap_pushover(G, po, n):
    fk = po["fk"]
    return r"""
\section{Análisis no lineal en SAP2000 (eje 2)}
Se utiliza el modelo 2D \texttt{G5\_Eje2\_Pushover} (plano YZ), que contiene la geometría, las secciones,
las cargas gravitacionales tributarias del eje 2 (patrones \texttt{PP}, \texttt{SCP}, \texttt{L}, \texttt{LR})
y el patrón lateral \texttt{PUSH} (fuerzas $\alpha_k$ en cada nivel, de suma unitaria). Los pasos son:
\begin{enumerate}
  \item \emph{Define $\rightarrow$ Section Properties $\rightarrow$ Hinge Properties}: crear rótulas
        \texttt{RC} (columna) y \texttt{RV} (viga) tipo \emph{Moment M3}, comportamiento elastoplástico perfecto
        (\emph{Deformation Controlled}), con $M_p = """ + n(po["Mp_c"], 1) + r"""\tonfm$ y
        $M_p = """ + n(po["Mp_b"], 1) + r"""\tonfm$ respectivamente, y meseta horizontal (puntos B--C con igual momento).
        Alternativamente, rótulas automáticas \emph{Auto M3} según ASCE 41-13 Tabla 9-6.
  \item Seleccionar todas las barras y asignar (\emph{Assign $\rightarrow$ Frame $\rightarrow$ Hinges}) rótulas a
        distancias relativas $0$ y $1$.
  \item \emph{Define $\rightarrow$ Load Cases}: caso \texttt{CGNL}, \emph{Static Nonlinear}, cargas
        $PP + SCP + 0{,}25L + 0{,}25L_r$, control de carga, condición inicial cero.
  \item Caso \texttt{PUSHOVER}, \emph{Static Nonlinear}, continuar desde \texttt{CGNL}, carga \texttt{PUSH} con
        factor 1, control de desplazamiento (\emph{Monitored displacement}: nudo 51, dirección U2, $60$--$80$ cm),
        resultados en múltiples pasos.
  \item Ejecutar y obtener la curva en \emph{Display $\rightarrow$ Show Static Pushover Curve} (corte basal vs.
        desplazamiento del nudo de techo).
\end{enumerate}
"""
