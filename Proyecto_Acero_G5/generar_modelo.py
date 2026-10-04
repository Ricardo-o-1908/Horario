# -*- coding: utf-8 -*-
"""
Proyecto Diseño en Acero (CIV-336) - USM - Segundo Semestre 2026 - GRUPO 5
==========================================================================

Genera, a partir de los datos del enunciado (L = 10 m, H = 3.0 m), para SAP2000 v27:

  1. G5_Edificio3D<suf>.s2k     -> modelo 3D completo (File > Import > SAP2000 .s2k)
  2. G5_Eje2_Pushover<suf>.s2k  -> modelo 2D del eje 2 para el análisis plástico (pushover)
  3. resultados_G5<suf>.json    -> resultados del análisis lineal propio (verificación del modelo)
  4. Memoria_Tarea1_G5<suf>.tex -> memoria de cálculo Tarea 1 en LaTeX (autocontenida)

Dos configuraciones de perfiles (argumento de línea de comandos):
  AISC -> perfiles laminados W y HSS (ASTM A572 Gr.50 / A500 Gr.C)          <suf> = ""
  CL   -> perfiles chilenos SOLDADOS (no conformados en frío): HN, IN y
          cajones soldados de planchas, acero NCh203 A270ES                  <suf> = "_CL"
Sin argumento se generan ambas.

Además resuelve el modelo con un análisis matricial propio (pórtico 3D con diafragma
rígido) para tener los esfuerzos de diseño, y realiza el pushover del eje 2
(evento a evento, rótulas elastoplásticas perfectas) y el método cinemático.

Unidades internas: tonf, m.

Uso:  python3 generar_modelo.py [AISC|CL]
Requiere: numpy
"""
import json
import sys
import math
import datetime
import numpy as np

# =============================================================================
# 1. PARÁMETROS DEL PROYECTO
# =============================================================================
GRUPO = 5
L = 10.0            # [m]  Tabla 2, grupo 5
H = 3.0             # [m]  Tabla 2, grupo 5
N_PISOS = 5

g = 9.80665
KPA = 1.0 / g       # 1 kPa = 0.10197 tonf/m2

GAMMA_ACERO = 7.85      # tonf/m3   (7850 kgf/m3)
GAMMA_HA = 2.50         # tonf/m3   (2500 kgf/m3)
E_LOSA = 0.18           # m
Q_TABIQUES = 0.100      # tonf/m2   (100 kgf/m2)
F_CONEX = 1.10          # +10 % peso propio por conexiones
Q_USO = 5.0 * KPA       # tonf/m2   (5 kPa)
# NCh1537 - techo transitable. >>> VERIFICAR con la Tabla 4 de la NCh1537 que usan en clases <<<
Q_TECHO_KPA = 2.0
Q_TECHO = Q_TECHO_KPA * KPA
FRAC_SC_SISMO = 0.25    # 25 % de la sobrecarga de uso en el peso sísmico
C_SISMO = 0.18          # Q = 0.18 P

# Configuración de perfiles: "AISC" (W + HSS) o "CL" (perfiles soldados chilenos)
CATALOGO = sys.argv[1].upper() if len(sys.argv) > 1 and sys.argv[1].upper() in ("AISC", "CL") else "AISC"
SUF = "" if CATALOGO == "AISC" else "_CL"
SAP_VERSION = "27.0.0"

# Acero  (nombre SAP, Fy [MPa], Fu [MPa], grado, descripción)
if CATALOGO == "AISC":
    MAT_PERFIL = ("A572Gr50", 345.0, 450.0, "Grade 50", "ASTM A572 Gr.50")
    MAT_RIOSTRA = ("A500GrC", 345.0, 427.0, "Grade C", "ASTM A500 Gr.C")
    MAT_PLANCHA = ("A36", 248.0, 400.0, "Grade 36", "ASTM A36")
else:
    # NCh203.Of2006: A270ES -> Fy = 270 MPa, Fu = 410 MPa (perfiles soldados y planchas)
    MAT_PERFIL = ("A270ES", 270.0, 410.0, "A270ES", "NCh203 A270ES")
    MAT_RIOSTRA = MAT_PERFIL
    MAT_PLANCHA = MAT_PERFIL
FY_W, FU_W = MAT_PERFIL[1], MAT_PERFIL[2]
FY_HSS, FU_HSS = MAT_RIOSTRA[1], MAT_RIOSTRA[2]
E_ACERO = 200000.0      # MPa
NU = 0.30
MPA = 1e6 / (g * 1e3)   # 1 MPa -> tonf/m2  (= 101.97)

# =============================================================================
# 2. GEOMETRÍA
# =============================================================================
XS = [0.0, 0.75 * L, 1.75 * L, 2.75 * L, 3.5 * L]          # ejes 1..5
YS = [0.0, 0.8 * L]                                         # ejes A, B
EJES_X = ["1", "2", "3", "4", "5"]
EJES_Y = ["A", "B"]
HP = [1.1 * H] + [H] * (N_PISOS - 1)                        # alturas de entrepiso
ZS = [0.0]
for h in HP:
    ZS.append(round(ZS[-1] + h, 6))
HTOT = ZS[-1]
XCM = (XS[0] + XS[-1]) / 2.0
YCM = (YS[0] + YS[-1]) / 2.0
AREA_PLANTA = (XS[-1] - XS[0]) * (YS[-1] - YS[0])


# =============================================================================
# 3. SECCIONES (prediseño)
# =============================================================================
IN = 0.0254
# Perfiles W (AISC): d, bf, tf, tw en pulgadas
W_DIMS = {
    "W14X90": (14.0, 14.5, 0.710, 0.440),
    "W14X74": (14.2, 10.1, 0.785, 0.450),
    "W21X62": (21.0, 8.24, 0.615, 0.400),
    "W21X50": (20.8, 6.53, 0.535, 0.380),
    "W18X50": (18.0, 7.50, 0.570, 0.355),
}
# HSS cuadrados (AISC): b, t nominal en pulgadas  (t diseño = 0.93 t nominal)
HSS_DIMS = {
    "HSS10X10X1/2": (10.0, 0.500),
    "HSS8X8X1/2": (8.0, 0.500),
}

# Perfiles soldados chilenos (nomenclatura tipo ICHA: serie, altura [cm] x peso [kgf/m]),
# armados con planchas de espesores comerciales. Dimensiones en mm: d, bf, tf, tw
CL_I_DIMS = {
    "HN40": (400, 400, 20, 10),     # columnas
    "IN50a": (500, 200, 14, 8),     # vigas X
    "IN50b": (500, 200, 16, 8),     # vigas Y
}
# Cajones soldados de 4 planchas (b, t en mm), no conformados en frío
CL_BOX_DIMS = {
    "CAJ250": (250, 12),            # riostras ejes 1 y 5
    "CAJ200": (200, 12),            # riostras ejes A y B
}

if CATALOGO == "AISC":
    SEC_COL = "W14X90"
    SEC_VIGA_X = "W21X50"      # vigas ejes A y B (luces 7.5 y 10 m)
    SEC_VIGA_Y = "W21X62"      # vigas ejes 1 a 5 (luz 8 m)
    SEC_RIO_X = "HSS8X8X1/2"   # riostras ejes A y B
    SEC_RIO_Y = "HSS10X10X1/2"  # riostras ejes 1 y 5


def props_I(d, bf, tf, tw):
    """Propiedades de perfil I a partir de las placas (sin filetes, igual que SAP2000)."""
    hw = d - 2 * tf
    A = 2 * bf * tf + hw * tw
    I33 = (bf * d ** 3 - (bf - tw) * hw ** 3) / 12.0
    I22 = 2 * tf * bf ** 3 / 12.0 + hw * tw ** 3 / 12.0
    J = (2 * bf * tf ** 3 + (d - tf) * tw ** 3) / 3.0
    Z33 = bf * tf * (d - tf) + tw * hw ** 2 / 4.0
    Z22 = tf * bf ** 2 / 2.0 + hw * tw ** 2 / 4.0
    return dict(A=A, I33=I33, I22=I22, J=J, S33=2 * I33 / d, S22=2 * I22 / bf, Z33=Z33, Z22=Z22,
                AS2=d * tw, AS3=5.0 / 3.0 * bf * tf, r33=math.sqrt(I33 / A), r22=math.sqrt(I22 / A))


def props_box(b, h, t):
    """Tubo rectangular de esquinas rectas (como Box/Tube de SAP2000)."""
    A = 2 * t * (b + h) - 4 * t * t
    I33 = (b * h ** 3 - (b - 2 * t) * (h - 2 * t) ** 3) / 12.0
    I22 = (h * b ** 3 - (h - 2 * t) * (b - 2 * t) ** 3) / 12.0
    Am = (b - t) * (h - t)
    J = 4 * Am ** 2 * t / (2 * ((b - t) + (h - t)))
    Z33 = b * h ** 2 / 4.0 - (b - 2 * t) * (h - 2 * t) ** 2 / 4.0
    Z22 = h * b ** 2 / 4.0 - (h - 2 * t) * (b - 2 * t) ** 2 / 4.0
    return dict(A=A, I33=I33, I22=I22, J=J, S33=2 * I33 / h, S22=2 * I22 / b, Z33=Z33, Z22=Z22,
                AS2=2 * h * t, AS3=2 * b * t, r33=math.sqrt(I33 / A), r22=math.sqrt(I22 / A))


def props_hss_aisc(b, tn, n=1500):
    """Propiedades de diseño de HSS cuadrado con esquinas redondeadas (r_ext = 2t), como el manual AISC."""
    t = 0.93 * tn
    ro, ri = 2 * t, t
    xs = (np.arange(n) + 0.5) / n * b - b / 2
    X, Y = np.meshgrid(xs, xs)
    dA = (b / n) ** 2

    def inside(a, r):
        ax, ay, hh = np.abs(X), np.abs(Y), a / 2
        m = (ax <= hh) & (ay <= hh)
        c = (ax > hh - r) & (ay > hh - r)
        return m & (~c | ((ax - (hh - r)) ** 2 + (ay - (hh - r)) ** 2 <= r * r))

    m = inside(b, ro) & ~inside(b - 2 * t, ri)
    A = m.sum() * dA
    I = (Y ** 2 * m).sum() * dA
    Z = (np.abs(Y) * m).sum() * dA
    return dict(b=b, tn=tn, t=t, A=A, I=I, r=math.sqrt(I / A), Z=Z, b_t=(b - 3 * t) / t)


SECCIONES = {}
if CATALOGO == "AISC":
    for nom, (d, bf, tf, tw) in W_DIMS.items():
        p = props_I(d * IN, bf * IN, tf * IN, tw * IN)
        p.update(tipo="I", mat=MAT_PERFIL[0], t3=d * IN, t2=bf * IN, tf=tf * IN, tw=tw * IN, Fy=FY_W, Fu=FU_W)
        SECCIONES[nom] = p
    for nom, (b, tn) in HSS_DIMS.items():
        t = 0.93 * tn * IN
        p = props_box(b * IN, b * IN, t)
        p.update(tipo="BOX", mat=MAT_RIOSTRA[0], t3=b * IN, t2=b * IN, tf=t, tw=t, Fy=FY_HSS, Fu=FU_HSS)
        SECCIONES[nom] = p
else:
    _cl = {}
    for clave, (d, bf, tf, tw) in CL_I_DIMS.items():
        p = props_I(d / 1e3, bf / 1e3, tf / 1e3, tw / 1e3)
        peso = p["A"] * GAMMA_ACERO * 1000          # kgf/m
        nom = "%s%gx%.0f" % (clave[:2], d / 10, peso)
        p.update(tipo="I", mat=MAT_PERFIL[0], t3=d / 1e3, t2=bf / 1e3, tf=tf / 1e3, tw=tw / 1e3, Fy=FY_W, Fu=FU_W,
                 dims_mm=(d, bf, tf, tw), peso=peso)
        SECCIONES[nom] = p
        _cl[clave] = nom
    for clave, (b, t) in CL_BOX_DIMS.items():
        p = props_box(b / 1e3, b / 1e3, t / 1e3)
        peso = p["A"] * GAMMA_ACERO * 1000
        nom = "CAJ%dx%dx%d" % (b, b, t)
        p.update(tipo="BOX", mat=MAT_RIOSTRA[0], t3=b / 1e3, t2=b / 1e3, tf=t / 1e3, tw=t / 1e3, Fy=FY_HSS,
                 Fu=FU_HSS, dims_mm=(b, t), peso=peso)
        SECCIONES[nom] = p
        _cl[clave] = nom
    SEC_COL = _cl["HN40"]
    SEC_VIGA_X = _cl["IN50a"]
    SEC_VIGA_Y = _cl["IN50b"]
    SEC_RIO_X = _cl["CAJ200"]
    SEC_RIO_Y = _cl["CAJ250"]

E = E_ACERO * MPA
G_MOD = E / (2 * (1 + NU))


# =============================================================================
# 4. NUDOS Y ELEMENTOS DEL MODELO 3D
# =============================================================================
class Nudo:
    def __init__(self, id, x, y, z, nivel):
        self.id, self.x, self.y, self.z, self.nivel = id, x, y, z, nivel


class Barra:
    def __init__(self, id, i, j, sec, tipo, eje, piso):
        self.id, self.i, self.j, self.sec, self.tipo, self.eje, self.piso = id, i, j, sec, tipo, eje, piso
        self.cargas = []   # cargas distribuidas: (patron, [(rel_a, w_a, rel_b, w_b), ...]) en tonf/m (gravedad)


NUDOS = {}
NUDO_ID = {}   # (ix, iy, k) -> id
for k in range(N_PISOS + 1):
    for iy, y in enumerate(YS):
        for ix, x in enumerate(XS):
            nid = 100 * k + 10 * iy + ix + 1        # ej.: 301 = nivel 3, eje A, eje 1
            NUDOS[nid] = Nudo(nid, x, y, ZS[k], k)
            NUDO_ID[(ix, iy, k)] = nid
NUDO_CM = {}
for k in range(1, N_PISOS + 1):
    nid = 100 * k + 99
    NUDOS[nid] = Nudo(nid, XCM, YCM, ZS[k], k)
    NUDO_CM[k] = nid

BARRAS = []


def nueva_barra(i, j, sec, tipo, eje, piso):
    b = Barra(len(BARRAS) + 1, i, j, sec, tipo, eje, piso)
    BARRAS.append(b)
    return b


# Columnas
for k in range(1, N_PISOS + 1):
    for iy in range(2):
        for ix in range(5):
            nueva_barra(NUDO_ID[(ix, iy, k - 1)], NUDO_ID[(ix, iy, k)], SEC_COL, "COL",
                        EJES_X[ix] + EJES_Y[iy], k)
# Vigas dirección X (ejes A y B)
for k in range(1, N_PISOS + 1):
    for iy in range(2):
        for ix in range(4):
            nueva_barra(NUDO_ID[(ix, iy, k)], NUDO_ID[(ix + 1, iy, k)], SEC_VIGA_X, "VX",
                        EJES_Y[iy] + " (" + EJES_X[ix] + "-" + EJES_X[ix + 1] + ")", k)
# Vigas dirección Y (ejes 1 a 5)
for k in range(1, N_PISOS + 1):
    for ix in range(5):
        nueva_barra(NUDO_ID[(ix, 0, k)], NUDO_ID[(ix, 1, k)], SEC_VIGA_Y, "VY", EJES_X[ix] + " (A-B)", k)
# Riostras ejes A y B: vano 1-2 "/" y vano 4-5 "\"
for k in range(1, N_PISOS + 1):
    for iy in range(2):
        nueva_barra(NUDO_ID[(0, iy, k - 1)], NUDO_ID[(1, iy, k)], SEC_RIO_X, "RX", EJES_Y[iy] + " (1-2)", k)
        nueva_barra(NUDO_ID[(3, iy, k)], NUDO_ID[(4, iy, k - 1)], SEC_RIO_X, "RX", EJES_Y[iy] + " (4-5)", k)
# Riostras ejes 1 y 5 (zig-zag): pisos impares "\" (A arriba -> B abajo), pares "/" (A abajo -> B arriba)
for k in range(1, N_PISOS + 1):
    for ix in (0, 4):
        if k % 2 == 1:
            nueva_barra(NUDO_ID[(ix, 0, k)], NUDO_ID[(ix, 1, k - 1)], SEC_RIO_Y, "RY", EJES_X[ix], k)
        else:
            nueva_barra(NUDO_ID[(ix, 0, k - 1)], NUDO_ID[(ix, 1, k)], SEC_RIO_Y, "RY", EJES_X[ix], k)


def largo(b):
    a, c = NUDOS[b.i], NUDOS[b.j]
    return math.sqrt((c.x - a.x) ** 2 + (c.y - a.y) ** 2 + (c.z - a.z) ** 2)


def es_riostra(b):
    return b.tipo in ("RX", "RY")


# =============================================================================
# 5. CARGAS DE LOSA (distribución en 2 direcciones, líneas a 45°)
# =============================================================================
Q_SCP = E_LOSA * GAMMA_HA + Q_TABIQUES     # tonf/m2


def perfil_tributario(luz, s, triangular):
    """Ancho tributario a lo largo de una viga (puntos relativos, ancho [m])."""
    if triangular:
        return [(0.0, 0.0), (0.5, s / 2.0), (1.0, 0.0)]
    a = (s / 2.0) / luz
    return [(0.0, 0.0), (a, s / 2.0), (1.0 - a, s / 2.0), (1.0, 0.0)]


def barra_entre(i, j):
    for b in BARRAS:
        if {b.i, b.j} == {i, j}:
            return b
    raise KeyError((i, j))


ANCHO_TRIB = {}   # id barra -> lista de perfiles tributarios (uno por paño)
for k in range(1, N_PISOS + 1):
    for ix in range(4):
        lx = XS[ix + 1] - XS[ix]
        ly = YS[1] - YS[0]
        s = min(lx, ly)
        # vigas X (bordes A y B del paño)
        for iy in range(2):
            b = barra_entre(NUDO_ID[(ix, iy, k)], NUDO_ID[(ix + 1, iy, k)])
            ANCHO_TRIB.setdefault(b.id, []).append(perfil_tributario(lx, s, lx <= ly))
        # vigas Y (bordes ix e ix+1 del paño)
        for jx in (ix, ix + 1):
            b = barra_entre(NUDO_ID[(jx, 0, k)], NUDO_ID[(jx, 1, k)])
            ANCHO_TRIB.setdefault(b.id, []).append(perfil_tributario(ly, s, ly < lx))

for b in BARRAS:
    if b.id not in ANCHO_TRIB:
        continue
    for perfil in ANCHO_TRIB[b.id]:
        segs = [(perfil[n][0], perfil[n][1], perfil[n + 1][0], perfil[n + 1][1]) for n in range(len(perfil) - 1)]
        b.cargas.append(("SCP", [(ra, wa * Q_SCP, rb, wb * Q_SCP) for ra, wa, rb, wb in segs]))
        q_sc = Q_TECHO if b.piso == N_PISOS else Q_USO
        pat = "LR" if b.piso == N_PISOS else "L"
        b.cargas.append((pat, [(ra, wa * q_sc, rb, wb * q_sc) for ra, wa, rb, wb in segs]))

# =============================================================================
# 6. PESO SÍSMICO Y FUERZAS SÍSMICAS
# =============================================================================
peso_acero_nivel = [0.0] * (N_PISOS + 1)
peso_acero_tipo = {}
for b in BARRAS:
    w = GAMMA_ACERO * SECCIONES[b.sec]["A"] * largo(b) * F_CONEX
    peso_acero_tipo[b.tipo] = peso_acero_tipo.get(b.tipo, 0.0) + w
    zi, zj = NUDOS[b.i].nivel, NUDOS[b.j].nivel
    peso_acero_nivel[zi] += w / 2
    peso_acero_nivel[zj] += w / 2

NIVELES = []
for k in range(1, N_PISOS + 1):
    scp = Q_SCP * AREA_PLANTA
    sc = (Q_TECHO if k == N_PISOS else Q_USO) * AREA_PLANTA
    Pk = peso_acero_nivel[k] + scp + FRAC_SC_SISMO * sc
    Ak = math.sqrt(1 - ZS[k - 1] / HTOT) - math.sqrt(1 - ZS[k] / HTOT)
    NIVELES.append(dict(k=k, z=ZS[k], acero=peso_acero_nivel[k], scp=scp, sc=sc,
                        sc25=FRAC_SC_SISMO * sc, P=Pk, A=Ak, AP=Ak * Pk))
P_TOTAL = sum(n["P"] for n in NIVELES)
Q_BASAL = C_SISMO * P_TOTAL
S_AP = sum(n["AP"] for n in NIVELES)
for n in NIVELES:
    n["F"] = n["AP"] / S_AP * Q_BASAL

# =============================================================================
# 7. COMBINACIONES NCh3171 (LRFD)
# =============================================================================
PATRONES = ["PP", "SCP", "L", "LR", "EX", "EY"]
COMBOS = [
    ("C1", "1.4D", {"PP": 1.4, "SCP": 1.4}),
    ("C2", "1.2D + 1.6L + 0.5Lr", {"PP": 1.2, "SCP": 1.2, "L": 1.6, "LR": 0.5}),
    ("C3", "1.2D + 1.6Lr + L", {"PP": 1.2, "SCP": 1.2, "L": 1.0, "LR": 1.6}),
    ("C4", "1.2D + L + 1.4Ex", {"PP": 1.2, "SCP": 1.2, "L": 1.0, "EX": 1.4}),
    ("C5", "1.2D + L - 1.4Ex", {"PP": 1.2, "SCP": 1.2, "L": 1.0, "EX": -1.4}),
    ("C6", "1.2D + L + 1.4Ey", {"PP": 1.2, "SCP": 1.2, "L": 1.0, "EY": 1.4}),
    ("C7", "1.2D + L - 1.4Ey", {"PP": 1.2, "SCP": 1.2, "L": 1.0, "EY": -1.4}),
    ("C8", "0.9D + 1.4Ex", {"PP": 0.9, "SCP": 0.9, "EX": 1.4}),
    ("C9", "0.9D - 1.4Ex", {"PP": 0.9, "SCP": 0.9, "EX": -1.4}),
    ("C10", "0.9D + 1.4Ey", {"PP": 0.9, "SCP": 0.9, "EY": 1.4}),
    ("C11", "0.9D - 1.4Ey", {"PP": 0.9, "SCP": 0.9, "EY": -1.4}),
]


# =============================================================================
# 8. ANÁLISIS MATRICIAL 3D PROPIO (verificación del modelo SAP2000)
# =============================================================================
def ejes_locales(b):
    a, c = NUDOS[b.i], NUDOS[b.j]
    ex = np.array([c.x - a.x, c.y - a.y, c.z - a.z])
    Lb = np.linalg.norm(ex)
    ex = ex / Lb
    if abs(ex[2]) > 0.999:          # vertical: eje 2 = +Y (ángulo 90° en SAP)
        ey = np.array([0.0, 1.0, 0.0])
    else:                           # horizontal/inclinada: eje 2 hacia arriba
        ey = np.array([0.0, 0.0, 1.0]) - ex[2] * ex
        ey = ey / np.linalg.norm(ey)
    ez = np.cross(ex, ey)
    return np.vstack([ex, ey, ez]), Lb


def k_local(b, Lb):
    s = SECCIONES[b.sec]
    A, Iz, Iy, J = s["A"], s["I33"], s["I22"], s["J"]
    k = np.zeros((12, 12))
    if es_riostra(b):   # extremos rotulados M2, M3 -> axial + torsión
        ea = E * A / Lb
        k[0, 0] = k[6, 6] = ea
        k[0, 6] = k[6, 0] = -ea
        gj = G_MOD * J / Lb
        k[3, 3] = k[9, 9] = gj
        k[3, 9] = k[9, 3] = -gj
        return k
    ea, gj = E * A / Lb, G_MOD * J / Lb
    k[0, 0] = k[6, 6] = ea; k[0, 6] = k[6, 0] = -ea
    k[3, 3] = k[9, 9] = gj; k[3, 9] = k[9, 3] = -gj
    a1, a2, a3, a4 = 12 * E * Iz / Lb ** 3, 6 * E * Iz / Lb ** 2, 4 * E * Iz / Lb, 2 * E * Iz / Lb
    for (r, c, v) in [(1, 1, a1), (1, 5, a2), (1, 7, -a1), (1, 11, a2), (5, 5, a3), (5, 7, -a2), (5, 11, a4),
                      (7, 7, a1), (7, 11, -a2), (11, 11, a3)]:
        k[r, c] = k[c, r] = v
    b1, b2, b3, b4 = 12 * E * Iy / Lb ** 3, 6 * E * Iy / Lb ** 2, 4 * E * Iy / Lb, 2 * E * Iy / Lb
    for (r, c, v) in [(2, 2, b1), (2, 4, -b2), (2, 8, -b1), (2, 10, -b2), (4, 4, b3), (4, 8, b2), (4, 10, b4),
                      (8, 8, b1), (8, 10, b2), (10, 10, b3)]:
        k[r, c] = k[c, r] = v
    return k


GAUSS_N = 400


def cargas_locales(b, patron, R, Lb):
    """Devuelve funciones q(s) [tonf/m] en ejes locales (x,y,z), muestreadas en GAUSS_N puntos."""
    s = (np.arange(GAUSS_N) + 0.5) / GAUSS_N * Lb
    qg = np.zeros((GAUSS_N, 3))        # en ejes globales
    if patron == "PP":
        qg[:, 2] -= GAMMA_ACERO * SECCIONES[b.sec]["A"] * F_CONEX
    for pat, segs in b.cargas:
        if pat != patron:
            continue
        rel = s / Lb
        for ra, wa, rb, wb in segs:
            m = (rel >= ra) & (rel < rb)
            qg[m, 2] -= wa + (wb - wa) * (rel[m] - ra) / (rb - ra)
    ql = qg @ R.T                       # componentes locales
    return s, ql


def cargas_equivalentes(b, patron, R, Lb):
    s, ql = cargas_locales(b, patron, R, Lb)
    ds = Lb / GAUSS_N
    xi = s / Lb
    N1, N2 = 1 - xi, xi
    H1, H2 = 1 - 3 * xi ** 2 + 2 * xi ** 3, Lb * (xi - 2 * xi ** 2 + xi ** 3)
    H3, H4 = 3 * xi ** 2 - 2 * xi ** 3, Lb * (-xi ** 2 + xi ** 3)
    f = np.zeros(12)
    if es_riostra(b):   # biarticulada: se reparte como viga simple
        f[0] = (N1 * ql[:, 0]).sum() * ds; f[6] = (N2 * ql[:, 0]).sum() * ds
        f[1] = (N1 * ql[:, 1]).sum() * ds; f[7] = (N2 * ql[:, 1]).sum() * ds
        f[2] = (N1 * ql[:, 2]).sum() * ds; f[8] = (N2 * ql[:, 2]).sum() * ds
        return f, s, ql
    f[0] = (N1 * ql[:, 0]).sum() * ds; f[6] = (N2 * ql[:, 0]).sum() * ds
    f[1] = (H1 * ql[:, 1]).sum() * ds; f[5] = (H2 * ql[:, 1]).sum() * ds
    f[7] = (H3 * ql[:, 1]).sum() * ds; f[11] = (H4 * ql[:, 1]).sum() * ds
    f[2] = (H1 * ql[:, 2]).sum() * ds; f[4] = -(H2 * ql[:, 2]).sum() * ds
    f[8] = (H3 * ql[:, 2]).sum() * ds; f[10] = -(H4 * ql[:, 2]).sum() * ds
    return f, s, ql


def analizar_3d():
    ids = sorted(NUDOS)
    pos = {nid: n for n, nid in enumerate(ids)}
    ndof = 6 * len(ids)
    K = np.zeros((ndof, ndof))
    F = {p: np.zeros(ndof) for p in PATRONES}
    datos = {}
    for b in BARRAS:
        R, Lb = ejes_locales(b)
        T = np.zeros((12, 12))
        for q in range(4):
            T[3 * q:3 * q + 3, 3 * q:3 * q + 3] = R
        kl = k_local(b, Lb)
        kg = T.T @ kl @ T
        dofs = [6 * pos[b.i] + q for q in range(6)] + [6 * pos[b.j] + q for q in range(6)]
        K[np.ix_(dofs, dofs)] += kg
        feq = {}
        for p in ("PP", "SCP", "L", "LR"):
            f, s, ql = cargas_equivalentes(b, p, R, Lb)
            feq[p] = (f, s, ql)
            F[p][dofs] += T.T @ f
        datos[b.id] = (R, Lb, T, kl, dofs, feq)
    for nl in NIVELES:
        nid = NUDO_CM[nl["k"]]
        F["EX"][6 * pos[nid] + 0] += nl["F"]
        F["EY"][6 * pos[nid] + 1] += nl["F"]
    # Restricciones: base empotrada; nudos CM: UZ, RX, RY restringidos
    fijos = set()
    for nid, n in NUDOS.items():
        if n.nivel == 0:
            fijos.update(6 * pos[nid] + q for q in range(6))
        if nid in NUDO_CM.values():
            fijos.update(6 * pos[nid] + q for q in (2, 3, 4))
    # Diafragma rígido: maestro = nudo CM de cada nivel (UX, UY, RZ)
    cols = []           # columnas de T (gdl reducidos)
    mapa = {}           # gdl completo -> lista (col, coef)
    for k in range(1, N_PISOS + 1):
        m = NUDO_CM[k]
        for q in (0, 1, 5):
            mapa[(m, q)] = len(cols)
            cols.append((m, q))
    Tm = []
    idx_red = {}
    for nid in ids:
        for q in range(6):
            gdl = 6 * pos[nid] + q
            if gdl in fijos:
                continue
            n = NUDOS[nid]
            if n.nivel > 0 and q in (0, 1, 5) and nid not in NUDO_CM.values():
                continue
            if nid in NUDO_CM.values() and q in (0, 1, 5):
                continue
            idx_red[gdl] = len(cols)
            cols.append((nid, q))
    nred = len(cols)
    Tr = np.zeros((ndof, nred))
    for gdl, c in idx_red.items():
        Tr[gdl, c] = 1.0
    for nid in ids:
        n = NUDOS[nid]
        if n.nivel == 0:
            continue
        m = NUDO_CM[n.nivel]
        cx, cy, cr = mapa[(m, 0)], mapa[(m, 1)], mapa[(m, 5)]
        dx, dy = n.x - XCM, n.y - YCM
        Tr[6 * pos[nid] + 0, cx] = 1.0
        Tr[6 * pos[nid] + 0, cr] = -dy
        Tr[6 * pos[nid] + 1, cy] = 1.0
        Tr[6 * pos[nid] + 1, cr] = dx
        Tr[6 * pos[nid] + 5, cr] = 1.0
    Kr = Tr.T @ K @ Tr
    U = {}
    for p in PATRONES:
        ur = np.linalg.solve(Kr, Tr.T @ F[p])
        U[p] = Tr @ ur
    # Esfuerzos por barra y patrón, en estaciones
    res = {}
    NST = 21
    for b in BARRAS:
        R, Lb, T, kl, dofs, feq = datos[b.id]
        st = np.linspace(0, Lb, NST)
        out = {}
        for p in PATRONES:
            ul = T @ U[p][dofs]
            fl = kl @ ul
            if p in feq:
                f, s, ql = feq[p]
                fl = fl - f
                ds = Lb / GAUSS_N
            else:
                s = ql = None
            # esfuerzos internos en estaciones (cuerpo libre desde el extremo i)
            N = np.zeros(NST); V2 = np.zeros(NST); V3 = np.zeros(NST)
            M3 = np.zeros(NST); M2 = np.zeros(NST); Tq = np.zeros(NST)
            for n_, x in enumerate(st):
                qx = qy = qz = mqy = mqz = 0.0
                if s is not None:
                    m = s < x
                    qx = ql[m, 0].sum() * ds
                    qy = ql[m, 1].sum() * ds
                    qz = ql[m, 2].sum() * ds
                    mqy = (ql[m, 1] * (x - s[m])).sum() * ds
                    mqz = (ql[m, 2] * (x - s[m])).sum() * ds
                N[n_] = -(fl[0] + qx)
                V2[n_] = fl[1] + qy
                V3[n_] = fl[2] + qz
                M3[n_] = x * fl[1] + mqy - fl[5]
                M2[n_] = -(x * fl[2] + mqz) - fl[4]
                Tq[n_] = -fl[3]
            out[p] = dict(N=N, V2=V2, V3=V3, M3=M3, M2=M2, T=Tq)
        res[b.id] = out
    return U, res, pos


def envolventes(res):
    env = {}
    for b in BARRAS:
        r = res[b.id]
        e = {}
        for cname, cdesc, fac in COMBOS:
            comb = {key: sum(f * r[p][key] for p, f in fac.items()) for key in ("N", "V2", "V3", "M3", "M2", "T")}
            e[cname] = comb
        env[b.id] = e
    return env


def resumen_por_tipo(env):
    out = {}
    for tipo in ("COL", "VX", "VY", "RX", "RY"):
        best = {}
        for b in BARRAS:
            if b.tipo != tipo:
                continue
            for cname, c in env[b.id].items():
                cand = {
                    "Pt": (c["N"].max(), c), "Pc": (-c["N"].min(), c),
                    "M3": (np.abs(c["M3"]).max(), c), "M2": (np.abs(c["M2"]).max(), c),
                    "V2": (np.abs(c["V2"]).max(), c), "V3": (np.abs(c["V3"]).max(), c),
                }
                for key, (val, _) in cand.items():
                    if key not in best or val > best[key][0]:
                        best[key] = (float(val), b.id, cname)
        out[tipo] = best
    return out


# =============================================================================
# 9. ANÁLISIS PLÁSTICO EJE 2 (pórtico 2D, 1 vano, 5 pisos)
# =============================================================================
def cargas_eje2():
    """Carga gravitacional distribuida sobre las vigas del eje 2 (perfil de w(x) por nivel) y
    cargas puntuales en las columnas (reacciones de vigas X que llegan al eje 2)."""
    ix = 1
    vig = {}
    for k in range(1, N_PISOS + 1):
        b = barra_entre(NUDO_ID[(ix, 0, k)], NUDO_ID[(ix, 1, k)])
        vig[k] = b
    return vig


def k2d(EI, EA, Lb, c, s, hi, hj):
    """Rigidez 2D global de barra con rótulas (hi, hj = True si el extremo está rotulado)."""
    kl = np.zeros((6, 6))
    kl[0, 0] = kl[3, 3] = EA / Lb
    kl[0, 3] = kl[3, 0] = -EA / Lb
    if not hi and not hj:
        a = np.array([[12, 6 * Lb, -12, 6 * Lb], [6 * Lb, 4 * Lb ** 2, -6 * Lb, 2 * Lb ** 2],
                      [-12, -6 * Lb, 12, -6 * Lb], [6 * Lb, 2 * Lb ** 2, -6 * Lb, 4 * Lb ** 2]]) * EI / Lb ** 3
    elif hi and not hj:
        a = np.array([[3, 0, -3, 3 * Lb], [0, 0, 0, 0], [-3, 0, 3, -3 * Lb], [3 * Lb, 0, -3 * Lb, 3 * Lb ** 2]]) * EI / Lb ** 3
    elif hj and not hi:
        a = np.array([[3, 3 * Lb, -3, 0], [3 * Lb, 3 * Lb ** 2, -3 * Lb, 0], [-3, -3 * Lb, 3, 0], [0, 0, 0, 0]]) * EI / Lb ** 3
    else:
        a = np.zeros((4, 4))
    idx = [1, 2, 4, 5]
    for r in range(4):
        for q in range(4):
            kl[idx[r], idx[q]] = a[r, q]
    T = np.array([[c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]])
    return kl, T


def pushover_eje2(frac_L=0.25):
    """Pushover evento a evento del eje 2 bajo D + 0.25L (gravitacional) + patrón lateral F_k."""
    ycols = [YS[0], YS[1]]
    sc = SECCIONES[SEC_COL]; sb = SECCIONES[SEC_VIGA_Y]
    Mp_c = sc["Z33"] * FY_W * MPA; Mp_b = sb["Z33"] * FY_W * MPA
    My_c = sc["S33"] * FY_W * MPA; My_b = sb["S33"] * FY_W * MPA
    # nudos 2D
    nod = {}
    for k in range(N_PISOS + 1):
        for c_, y in enumerate(ycols):
            nod[(k, c_)] = (y, ZS[k])
    keys = sorted(nod)
    pos = {key: i for i, key in enumerate(keys)}
    ndof = 3 * len(keys)
    elems = []    # (ni, nj, sec, tipo, piso, w_profile(None))
    for k in range(1, N_PISOS + 1):
        for c_ in range(2):
            elems.append(dict(i=(k - 1, c_), j=(k, c_), sec=SEC_COL, tipo="COL", piso=k))
    vig = cargas_eje2()
    for k in range(1, N_PISOS + 1):
        elems.append(dict(i=(k, 0), j=(k, 1), sec=SEC_VIGA_Y, tipo="VIG", piso=k, b3d=vig[k]))
    # carga gravitacional sobre vigas: PP*1.1 + SCP + frac*(L o LR)
    NG = GAUSS_N
    for e in elems:
        yi, zi = nod[e["i"]]; yj, zj = nod[e["j"]]
        Lb = math.hypot(yj - yi, zj - zi)
        e["L"] = Lb; e["c"] = (yj - yi) / Lb; e["s"] = (zj - zi) / Lb
        s_ = SECCIONES[e["sec"]]
        e["EI"] = E * s_["I33"]; e["EA"] = E * s_["A"]
        e["Mp"] = Mp_c if e["tipo"] == "COL" else Mp_b
        e["My"] = My_c if e["tipo"] == "COL" else My_b
        xs = (np.arange(NG) + 0.5) / NG * Lb
        q = np.zeros(NG)   # carga transversal local (hacia +y local); para vigas = vertical
        if e["tipo"] == "VIG":
            q -= GAMMA_ACERO * s_["A"] * F_CONEX
            for pat, segs in e["b3d"].cargas:
                fac = 1.0 if pat == "SCP" else frac_L
                rel = xs / Lb
                for ra, wa, rb, wb in segs:
                    m = (rel >= ra) & (rel < rb)
                    q[m] -= fac * (wa + (wb - wa) * (rel[m] - ra) / (rb - ra))
        e["xs"] = xs; e["q"] = q
    # cargas puntuales en columnas: reacciones (viga simple) de vigas X de los paños 1-2 y 2-3
    Pcol = {}
    for k in range(1, N_PISOS + 1):
        for c_, iy in enumerate((0, 1)):
            tot = 0.0
            for (a, bb) in ((0, 1), (1, 2)):
                bx = barra_entre(NUDO_ID[(a, iy, k)], NUDO_ID[(bb, iy, k)])
                Lx = largo(bx)
                w = GAMMA_ACERO * SECCIONES[bx.sec]["A"] * F_CONEX * Lx
                for pat, segs in bx.cargas:
                    fac = 1.0 if pat == "SCP" else frac_L
                    for ra, wa, rb, wb in segs:
                        w += fac * (wa + wb) / 2 * (rb - ra) * Lx
                tot += w / 2
            # peso propio columnas (mitad superior e inferior) incluido como puntual
            Pcol[(k, c_)] = tot
    for k in range(1, N_PISOS + 1):
        for c_ in range(2):
            wcol = GAMMA_ACERO * SECCIONES[SEC_COL]["A"] * F_CONEX
            Pcol[(k, c_)] += wcol * HP[k - 1] / 2 + (wcol * HP[k] / 2 if k < N_PISOS else 0.0)

    fk = np.array([n["F"] for n in NIVELES]) / Q_BASAL   # patrón lateral normalizado (suma 1)
    fijos = [3 * pos[(0, c_)] + q for c_ in range(2) for q in range(3)]
    libres = [d for d in range(ndof) if d not in fijos]

    def ensamblar(rot):
        K = np.zeros((ndof, ndof))
        for n_, e in enumerate(elems):
            kl, T = k2d(e["EI"], e["EA"], e["L"], e["c"], e["s"], rot[(n_, 0)], rot[(n_, 1)])
            d = [3 * pos[e["i"]] + q for q in range(3)] + [3 * pos[e["j"]] + q for q in range(3)]
            K[np.ix_(d, d)] += T.T @ kl @ T
            e["_kl"], e["_T"], e["_d"] = kl, T, d
        return K

    def momentos_extremos(u, rot):
        M = {}
        for n_, e in enumerate(elems):
            ul = e["_T"] @ u[e["_d"]]
            fl = e["_kl"] @ ul
            # momento interno (convención: positivo = tracción en cara +y local... se usa en valor con signo)
            M[(n_, 0)] = -fl[2]
            M[(n_, 1)] = fl[5]
        return M

    # --- Etapa gravitacional (elástica) ---
    rot = {(n_, q): False for n_ in range(len(elems)) for q in (0, 1)}
    K = ensamblar(rot)
    Fg = np.zeros(ndof)
    fem = {}
    for n_, e in enumerate(elems):
        Lb = e["L"]; xi = e["xs"] / Lb; ds = Lb / NG
        H1, H2 = 1 - 3 * xi ** 2 + 2 * xi ** 3, Lb * (xi - 2 * xi ** 2 + xi ** 3)
        H3, H4 = 3 * xi ** 2 - 2 * xi ** 3, Lb * (-xi ** 2 + xi ** 3)
        f = np.array([0, (H1 * e["q"]).sum() * ds, (H2 * e["q"]).sum() * ds,
                      0, (H3 * e["q"]).sum() * ds, (H4 * e["q"]).sum() * ds])
        fem[n_] = f
        Fg[e["_d"]] += e["_T"].T @ f
    for (k, c_), P in Pcol.items():
        Fg[3 * pos[(k, c_)] + 1] -= P
    ug = np.zeros(ndof)
    ug[libres] = np.linalg.solve(K[np.ix_(libres, libres)], Fg[libres])
    Mg = {}
    for n_, e in enumerate(elems):
        fl = e["_kl"] @ (e["_T"] @ ug[e["_d"]]) - fem[n_]
        Mg[(n_, 0)] = -fl[2]
        Mg[(n_, 1)] = fl[5]
    # momento máximo de tramo en vigas bajo gravedad (verificación de no rotulación interior)
    Mtramo = 0.0
    for n_, e in enumerate(elems):
        if e["tipo"] != "VIG":
            continue
        fl = e["_kl"] @ (e["_T"] @ ug[e["_d"]]) - fem[n_]
        ds = e["L"] / NG
        for x in np.linspace(0, e["L"], 41):
            m = e["xs"] < x
            Mx = x * fl[1] + (e["q"][m] * (x - e["xs"][m])).sum() * ds - fl[2]
            Mtramo = max(Mtramo, abs(Mx))

    techo = 3 * pos[(N_PISOS, 0)] + 0
    Fl = np.zeros(ndof)
    for k in range(1, N_PISOS + 1):
        Fl[3 * pos[(k, 0)] + 0] += fk[k - 1]     # carga lateral aplicada en el nivel (eje A del pórtico)

    # --- Primera fluencia (M = My = S*Fy) con análisis elástico ---
    ul1 = np.zeros(ndof)
    ul1[libres] = np.linalg.solve(K[np.ix_(libres, libres)], Fl[libres])
    Ml1 = momentos_extremos(ul1, rot)
    lam_y = 1e9; loc_y = None
    for key in Ml1:
        e = elems[key[0]]
        for sgn in (1, -1):
            dm = sgn * Ml1[key]
            if dm > 1e-12:
                lam = (e["My"] - sgn * Mg[key]) / dm
                if lam < lam_y:
                    lam_y, loc_y = lam, key
    d_y = ug[techo] + lam_y * ul1[techo]

    # --- Pushover evento a evento ---
    lam = 0.0
    u = ug.copy()
    M = dict(Mg)
    curva = [(float(ug[techo]), 0.0)]
    eventos = []
    for paso in range(40):
        K = ensamblar(rot)
        Kff = K[np.ix_(libres, libres)]
        if np.linalg.matrix_rank(Kff, tol=1e-6 * np.abs(Kff).max()) < len(libres):
            break
        du = np.zeros(ndof)
        du[libres] = np.linalg.solve(Kff, Fl[libres])
        dM = momentos_extremos(du, rot)
        best = 1e9; bkey = None
        for key, val in dM.items():
            if rot[key] or abs(val) < 1e-10:
                continue
            e = elems[key[0]]
            lim = e["Mp"] if val > 0 else -e["Mp"]
            dl = (lim - M[key]) / val
            if 0 <= dl < best:
                best, bkey = dl, key
        if bkey is None or du[techo] > 10.0:
            break
        lam += best
        u = u + best * du
        for key in M:
            if not rot[key]:
                M[key] += best * dM[key]
        rot[bkey] = True
        e = elems[bkey[0]]
        eventos.append(dict(n=len(eventos) + 1, F=float(lam), d=float(u[techo]),
                            elem=("Columna" if e["tipo"] == "COL" else "Viga"), piso=e["piso"],
                            extremo=("inferior" if e["tipo"] == "COL" and bkey[1] == 0 else
                                     "superior" if e["tipo"] == "COL" else
                                     ("eje A" if bkey[1] == 0 else "eje B")),
                            lado=(["A", "B"][e["i"][1]] if e["tipo"] == "COL" else "")))
        curva.append((float(u[techo]), float(lam)))
    # tramo final (meseta) para graficar
    curva.append((curva[-1][0] * 1.25, curva[-1][1]))
    # --- Método cinemático ---
    zs = ZS
    mecs = []
    for j in range(1, N_PISOS + 1):
        for m in range(j, N_PISOS + 1):
            Wi = 2 * Mp_c   # rótulas al pie de las columnas del piso j
            for k in range(j, m):
                Wi += 2 * min(Mp_b, 2 * Mp_c)
            Wi += 2 * (min(Mp_b, Mp_c) if m == N_PISOS else Mp_c)
            We = 0.0
            for k in range(1, N_PISOS + 1):
                if k < j:
                    dk = 0.0
                elif k <= m:
                    dk = zs[k] - zs[j - 1]
                else:
                    dk = zs[m] - zs[j - 1]
                We += fk[k - 1] * dk
            if j == m:
                nombre = "Piso %d (piso blando)" % j
            elif j == 1 and m == N_PISOS:
                nombre = "Global (vigas + base)"
            else:
                nombre = "Combinado pisos %d a %d" % (j, m)
            mecs.append(dict(j=j, m=m, nombre=nombre, Wi=Wi, We=We, F=Wi / We))
    mecs.sort(key=lambda d: d["F"])
    return dict(Mp_c=Mp_c, Mp_b=Mp_b, My_c=My_c, My_b=My_b, lam_y=lam_y,
                loc_y=dict(elem="Columna" if elems[loc_y[0]]["tipo"] == "COL" else "Viga",
                           piso=elems[loc_y[0]]["piso"], extremo=loc_y[1]),
                d_y=float(d_y), curva=curva, eventos=eventos, F_colapso=lam, mecanismos=mecs,
                fk=fk.tolist(), Mtramo=Mtramo, Pcol=Pcol, elems2d=elems)


# =============================================================================
# 10. DISEÑO A TRACCIÓN DE LA RIOSTRA MÁS TRACCIONADA (AISC 360-16, cap. D y J)
# =============================================================================
def disenar_traccion(resumen):
    cand = []
    for tipo in ("RX", "RY"):
        Tu, bid, combo = resumen[tipo]["Pt"]
        cand.append((Tu, bid, combo, tipo))
    Tu, bid, combo, tipo = max(cand)
    b = [x for x in BARRAS if x.id == bid][0]
    Lb = largo(b)
    if CATALOGO == "AISC":
        bb, tn = HSS_DIMS[b.sec]
        h = props_hss_aisc(bb, tn)                      # esquinas redondeadas, t = 0.93 t_nom
        cm = 2.54
        A = h["A"] * cm ** 2; t = h["t"] * cm; B = bb * cm; r = h["r"] * cm
        b_t = h["b_t"]
    else:
        s_ = SECCIONES[b.sec]                           # cajón soldado: esquinas rectas, t nominal
        A = s_["A"] * 1e4; t = s_["tf"] * 100; B = s_["t2"] * 100; r = s_["r33"] * 100
        b_t = (B - 2 * t) / t
    Fy = FY_HSS * 10.197; Fu = FU_HSS * 10.197          # kgf/cm2
    # Fluencia en área bruta
    phiTn_y = 0.90 * Fy * A / 1000.0                    # tonf
    # Conexión: plancha gusset concéntrica, tubo ranurado, 4 soldaduras de filete E70
    tg = 2.2                                            # cm espesor gusset (7/8")
    holgura = 0.2                                       # cm (ranura = tg + 2 mm)
    An = A - 2 * t * (tg + holgura)
    FEXX = 4920.0                                       # kgf/cm2 (E70XX)
    w = 0.8                                             # cm, filete 8 mm
    w_min = 0.5
    # largo de soldadura requerido (4 cordones)
    Rw_cm = 0.75 * 0.60 * FEXX * 0.707 * w              # kgf/cm de cordón
    # resistencia al corte del metal base (pared del tubo): 0.75*0.6*Fu*t
    Rbase_cm = 0.75 * 0.60 * Fu * t
    Rcm = min(Rw_cm, Rbase_cm)
    l_req = Tu * 1000.0 / (4 * Rcm)
    l = math.ceil(max(l_req, B) / 5.0) * 5.0             # l >= B (para usar la tabla D3.1 caso 6)
    xbar = (B ** 2 + 2 * B * B) / (4 * (B + B))
    U = 1 - xbar / l
    Ae = U * An
    phiTn_u = 0.75 * Fu * Ae / 1000.0
    phiTn = min(phiTn_y, phiTn_u)
    # Gusset: fluencia en sección de Whitmore y bloque de corte
    Lw = B + 2 * l * math.tan(math.radians(30))
    Fy_g = MAT_PLANCHA[1] * 10.197; Fu_g = MAT_PLANCHA[2] * 10.197   # plancha gusset, kgf/cm2
    phiRn_whit = 0.90 * Fy_g * Lw * tg / 1000.0
    Agv = 2 * l * tg; Ant = B * tg
    phiRn_bs = 0.75 * min(0.6 * Fu_g * Agv + Fu_g * Ant, 0.6 * Fy_g * Agv + Fu_g * Ant) / 1000.0
    esbeltez = Lb * 100 / r
    return dict(Tu=Tu, barra=bid, combo=combo, tipo=tipo, eje=b.eje, piso=b.piso, L=Lb, sec=b.sec,
                A=A, t=t, B=B, r=r, Fy=Fy, Fu=Fu, phiTn_y=phiTn_y, An=An, tg=tg, holgura=holgura,
                xbar=xbar, l=l, l_req=l_req, U=U, Ae=Ae, phiTn_u=phiTn_u, phiTn=phiTn,
                FU=Tu / phiTn, w=w, w_min=w_min, Rw_cm=Rw_cm, Rbase_cm=Rbase_cm, FEXX=FEXX,
                Lw=Lw, phiRn_whit=phiRn_whit, phiRn_bs=phiRn_bs, Agv=Agv, Ant=Ant,
                Fy_g=Fy_g, Fu_g=Fu_g, esbeltez=esbeltez, b_t=b_t,
                mat_g=MAT_PLANCHA[4], mat_r=MAT_RIOSTRA[4], catalogo=CATALOGO)


# =============================================================================
# 11. ESCRITURA DE ARCHIVOS SAP2000 (.s2k, formato tablas v19)
# =============================================================================
def fmt(v):
    if isinstance(v, str):
        return '"%s"' % v if (" " in v or "/" in v or "-" in v or "," in v) else v
    if isinstance(v, bool):
        return "Yes" if v else "No"
    if abs(v) < 1e-12:
        return "0"
    return ("%.6g" % v) if abs(v) >= 1e-3 else ("%.6E" % v)


def fila(**kw):
    return "   " + "   ".join("%s=%s" % (k, fmt(v)) for k, v in kw.items())


def tabla(nombre, filas):
    return ['TABLE:  "%s"' % nombre] + filas + [""]


def bloque_comun(lineas, titulo):
    hoy = datetime.datetime.now()
    lineas += ["File %s was saved on %s" % (titulo, hoy.strftime("%m/%d/%y at %H:%M:%S")), ""]
    lineas += tabla("PROGRAM CONTROL", [
        '   ProgramName=SAP2000   Version=%s   ProgLevel=Ultimate   CurrUnits="Tonf, m, C"' % SAP_VERSION])
    lineas += tabla("COORDINATE SYSTEMS", [fila(Name="GLOBAL", Type="Cartesian", X=0, Y=0, Z=0,
                                                AboutZ=0, AboutY=0, AboutX=0)])
    mats = []
    for m in (MAT_PERFIL, MAT_RIOSTRA):
        if m[0] not in [x[0] for x in mats]:
            mats.append((m[0], m[1], m[2], m[3]))
    lineas += tabla("MATERIAL PROPERTIES 01 - GENERAL",
                    [fila(Material=m, Type="Steel", SymType="Isotropic", TempDepend=False,
                          Color="Blue") for m, _, _, gr in mats])
    lineas += tabla("MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES",
                    [fila(Material=m, UnitWeight=GAMMA_ACERO, UnitMass=GAMMA_ACERO / g, E1=E, G12=G_MOD,
                          U12=NU, A1=1.17e-5) for m, _, _, _ in mats])
    lineas += tabla("MATERIAL PROPERTIES 03A - STEEL DATA",
                    [fila(Material=m, Fy=fy * MPA, Fu=fu * MPA, EffFy=1.1 * fy * MPA, EffFu=1.1 * fu * MPA,
                          SSCurveOpt="Simple", SSHysType="Kinematic", SHard=0.02, SMax=0.14, SRup=0.2,
                          FinalSlope=-0.1, CoupModType="Von Mises") for m, fy, fu, _ in mats])


def filas_secciones(nombres):
    f = []
    for nom in nombres:
        s = SECCIONES[nom]
        shape = "I/Wide Flange" if s["tipo"] == "I" else "Box/Tube"
        kw = dict(SectionName=nom, Material=s["mat"], Shape=shape, t3=s["t3"], t2=s["t2"], tf=s["tf"], tw=s["tw"])
        if s["tipo"] == "I":
            kw.update(t2b=s["t2"], tfb=s["tf"])
        kw.update(Area=s["A"], TorsConst=s["J"], I33=s["I33"], I22=s["I22"], I23=0, AS2=s["AS2"], AS3=s["AS3"],
                  S33=s["S33"], S22=s["S22"], Z33=s["Z33"], Z22=s["Z22"], R33=s["r33"], R22=s["r22"],
                  ConcCol=False, ConcBeam=False, Color=("Cyan" if s["tipo"] == "I" else "Red"),
                  FromFile=False, AMod=1, A2Mod=1, A3Mod=1, JMod=1, I2Mod=1, I3Mod=1, MMod=1, WMod=1)
        f.append(fila(**kw))
    return f


def escribir_s2k_3d(ruta):
    lin = []
    bloque_comun(lin, "G5_Edificio3D%s.s2k" % SUF)
    lin += tabla("ACTIVE DEGREES OF FREEDOM", ["   UX=Yes   UY=Yes   UZ=Yes   RX=Yes   RY=Yes   RZ=Yes"])
    g_ = []
    for ix, x in enumerate(XS):
        g_.append(fila(CoordSys="GLOBAL", AxisDir="X", GridID=EJES_X[ix], XRYZCoord=x, LineType="Primary",
                       LineColor="Gray8Dark", Visible=True, BubbleLoc="End", AllVisible=True, BubbleSize=1.25))
    for iy, y in enumerate(YS):
        g_.append(fila(CoordSys="GLOBAL", AxisDir="Y", GridID=EJES_Y[iy], XRYZCoord=y, LineType="Primary",
                       LineColor="Gray8Dark", Visible=True, BubbleLoc="Start"))
    for k, z in enumerate(ZS):
        g_.append(fila(CoordSys="GLOBAL", AxisDir="Z", GridID=("BASE" if k == 0 else "N%d" % k), XRYZCoord=z,
                       LineType="Primary", LineColor="Gray8Dark", Visible=True, BubbleLoc="End"))
    lin += tabla("GRID LINES", g_)
    lin += tabla("JOINT COORDINATES", [fila(Joint=n.id, CoordSys="GLOBAL", CoordType="Cartesian", XorR=n.x, Y=n.y,
                                            Z=n.z, SpecialJt=(n.id in NUDO_CM.values()), GlobalX=n.x, GlobalY=n.y,
                                            GlobalZ=n.z) for n in sorted(NUDOS.values(), key=lambda q: q.id)])
    cf = []
    for b in BARRAS:
        a, c = NUDOS[b.i], NUDOS[b.j]
        cf.append(fila(Frame=b.id, JointI=b.i, JointJ=b.j, IsCurved=False, Length=largo(b),
                       CentroidX=(a.x + c.x) / 2, CentroidY=(a.y + c.y) / 2, CentroidZ=(a.z + c.z) / 2))
    lin += tabla("CONNECTIVITY - FRAME", cf)
    rest = []
    for n in sorted(NUDOS.values(), key=lambda q: q.id):
        if n.nivel == 0:
            rest.append(fila(Joint=n.id, U1=True, U2=True, U3=True, R1=True, R2=True, R3=True))
        elif n.id in NUDO_CM.values():
            rest.append(fila(Joint=n.id, U1=False, U2=False, U3=True, R1=True, R2=True, R3=False))
    lin += tabla("JOINT RESTRAINT ASSIGNMENTS", rest)
    lin += tabla("CONSTRAINT DEFINITIONS - DIAPHRAGM",
                 [fila(Name="DIAF%d" % k, CoordSys="GLOBAL", Axis="Z", MultiLevel=False) for k in range(1, N_PISOS + 1)])
    lin += tabla("JOINT CONSTRAINT ASSIGNMENTS",
                 [fila(Joint=n.id, Constraint="DIAF%d" % n.nivel, Type="Diaphragm")
                  for n in sorted(NUDOS.values(), key=lambda q: q.id) if n.nivel > 0])
    lin += tabla("FRAME SECTION PROPERTIES 01 - GENERAL",
                 filas_secciones([SEC_COL, SEC_VIGA_X, SEC_VIGA_Y, SEC_RIO_X, SEC_RIO_Y]))
    lin += tabla("FRAME SECTION ASSIGNMENTS",
                 [fila(Frame=b.id, SectionType=("I/Wide Flange" if SECCIONES[b.sec]["tipo"] == "I" else "Box/Tube"),
                       AutoSelect="N.A.", AnalSect=b.sec, DesignSect=b.sec, MatProp="Default") for b in BARRAS])
    lin += tabla("FRAME LOCAL AXES ASSIGNMENTS 1 - TYPICAL",
                 [fila(Frame=b.id, Angle=90, AdvanceAxes=False) for b in BARRAS if b.tipo == "COL"])
    lin += tabla("FRAME RELEASE ASSIGNMENTS 1 - GENERAL",
                 [fila(Frame=b.id, PI=False, V2I=False, V3I=False, TI=False, M2I=True, M3I=True,
                       PJ=False, V2J=False, V3J=False, TJ=False, M2J=True, M3J=True, PartialFix=False)
                  for b in BARRAS if es_riostra(b)])
    lin += tabla("LOAD PATTERN DEFINITIONS", [
        fila(LoadPat="PP", DesignType="Dead", SelfWtMult=F_CONEX),
        fila(LoadPat="SCP", DesignType="Super Dead", SelfWtMult=0),
        fila(LoadPat="L", DesignType="Live", SelfWtMult=0),
        fila(LoadPat="LR", DesignType="Roof Live", SelfWtMult=0),
        fila(LoadPat="EX", DesignType="Quake", SelfWtMult=0, AutoLoad="None"),
        fila(LoadPat="EY", DesignType="Quake", SelfWtMult=0, AutoLoad="None")])
    tipos_d = {"PP": "Dead", "SCP": "Super Dead", "L": "Live", "LR": "Roof Live", "EX": "Quake", "EY": "Quake"}
    lin += tabla("LOAD CASE DEFINITIONS", [
        fila(Case=p, Type="LinStatic", InitialCond="Zero", DesTypeOpt="Prog Det", DesignType=tipos_d[p],
             AutoType="None", RunCase=True) for p in PATRONES])
    lin += tabla("CASE - STATIC 1 - LOAD ASSIGNMENTS",
                 [fila(Case=p, LoadType="Load pattern", LoadName=p, LoadSF=1) for p in PATRONES])
    cd = []
    for cname, desc, fac in COMBOS:
        first = True
        for p, f in fac.items():
            if first:
                cd.append(fila(ComboName=cname, ComboType="Linear Add", AutoDesign=False, CaseType="Linear Static",
                               CaseName=p, ScaleFactor=f, SteelDesign="Strength", ConcDesign="None",
                               AlumDesign="None", ColdDesign="None"))
                first = False
            else:
                cd.append(fila(ComboName=cname, CaseType="Linear Static", CaseName=p, ScaleFactor=f))
    first = True
    for cname, _, _ in COMBOS:
        if first:
            cd.append(fila(ComboName="ENVOLVENTE", ComboType="Envelope", AutoDesign=False, CaseType="Response Combo",
                           CaseName=cname, ScaleFactor=1, SteelDesign="None", ConcDesign="None",
                           AlumDesign="None", ColdDesign="None"))
            first = False
        else:
            cd.append(fila(ComboName="ENVOLVENTE", CaseType="Response Combo", CaseName=cname, ScaleFactor=1))
    lin += tabla("COMBINATION DEFINITIONS", cd)
    dl = []
    for b in BARRAS:
        Lb = largo(b)
        for pat, segs in b.cargas:
            for ra, wa, rb, wb in segs:
                dl.append(fila(Frame=b.id, LoadPat=pat, CoordSys="GLOBAL", Type="Force", Dir="Gravity",
                               DistType="RelDist", RelDistA=ra, RelDistB=rb, AbsDistA=ra * Lb, AbsDistB=rb * Lb,
                               FOverLA=wa, FOverLB=wb))
    lin += tabla("FRAME LOADS - DISTRIBUTED", dl)
    jl = []
    for nl in NIVELES:
        nid = NUDO_CM[nl["k"]]
        jl.append(fila(Joint=nid, LoadPat="EX", CoordSys="GLOBAL", F1=nl["F"], F2=0, F3=0, M1=0, M2=0, M3=0))
        jl.append(fila(Joint=nid, LoadPat="EY", CoordSys="GLOBAL", F1=0, F2=nl["F"], F3=0, M1=0, M2=0, M3=0))
    lin += tabla("JOINT LOADS - FORCE", jl)
    lin += ["END TABLE DATA"]
    with open(ruta, "w", newline="\r\n") as fh:
        fh.write("\n".join(lin) + "\n")


def escribir_s2k_eje2(ruta, po):
    """Modelo 2D del eje 2 (plano YZ, X = 7.5 m) para el pushover."""
    lin = []
    bloque_comun(lin, "G5_Eje2_Pushover%s.s2k" % SUF)
    lin += tabla("ACTIVE DEGREES OF FREEDOM", ["   UX=No   UY=Yes   UZ=Yes   RX=Yes   RY=No   RZ=No"])
    g_ = [fila(CoordSys="GLOBAL", AxisDir="X", GridID="2", XRYZCoord=XS[1], LineType="Primary",
               LineColor="Gray8Dark", Visible=True, BubbleLoc="End", AllVisible=True, BubbleSize=1.25)]
    for iy, y in enumerate(YS):
        g_.append(fila(CoordSys="GLOBAL", AxisDir="Y", GridID=EJES_Y[iy], XRYZCoord=y, LineType="Primary",
                       LineColor="Gray8Dark", Visible=True, BubbleLoc="Start"))
    for k, z in enumerate(ZS):
        g_.append(fila(CoordSys="GLOBAL", AxisDir="Z", GridID=("BASE" if k == 0 else "N%d" % k), XRYZCoord=z,
                       LineType="Primary", LineColor="Gray8Dark", Visible=True, BubbleLoc="End"))
    lin += tabla("GRID LINES", g_)
    jn = []
    for k in range(N_PISOS + 1):
        for c_, y in enumerate(YS):
            jn.append(fila(Joint=10 * k + c_ + 1, CoordSys="GLOBAL", CoordType="Cartesian", XorR=XS[1], Y=y, Z=ZS[k],
                           SpecialJt=False, GlobalX=XS[1], GlobalY=y, GlobalZ=ZS[k]))
    lin += tabla("JOINT COORDINATES", jn)
    cf, sa, la = [], [], []
    fid = 0
    vig_id = {}
    for k in range(1, N_PISOS + 1):
        for c_ in range(2):
            fid += 1
            cf.append(fila(Frame=fid, JointI=10 * (k - 1) + c_ + 1, JointJ=10 * k + c_ + 1, IsCurved=False,
                           Length=HP[k - 1], CentroidX=XS[1], CentroidY=YS[c_], CentroidZ=(ZS[k - 1] + ZS[k]) / 2))
            sa.append(fila(Frame=fid, SectionType="I/Wide Flange", AutoSelect="N.A.", AnalSect=SEC_COL,
                           DesignSect=SEC_COL, MatProp="Default"))
            la.append(fila(Frame=fid, Angle=90, AdvanceAxes=False))
    for k in range(1, N_PISOS + 1):
        fid += 1
        vig_id[k] = fid
        cf.append(fila(Frame=fid, JointI=10 * k + 1, JointJ=10 * k + 2, IsCurved=False, Length=YS[1],
                       CentroidX=XS[1], CentroidY=YCM, CentroidZ=ZS[k]))
        sa.append(fila(Frame=fid, SectionType="I/Wide Flange", AutoSelect="N.A.", AnalSect=SEC_VIGA_Y,
                       DesignSect=SEC_VIGA_Y, MatProp="Default"))
    lin += tabla("CONNECTIVITY - FRAME", cf)
    lin += tabla("JOINT RESTRAINT ASSIGNMENTS",
                 [fila(Joint=c_ + 1, U1=True, U2=True, U3=True, R1=True, R2=True, R3=True) for c_ in range(2)])
    lin += tabla("FRAME SECTION PROPERTIES 01 - GENERAL", filas_secciones([SEC_COL, SEC_VIGA_Y]))
    lin += tabla("FRAME SECTION ASSIGNMENTS", sa)
    lin += tabla("FRAME LOCAL AXES ASSIGNMENTS 1 - TYPICAL", la)
    lin += tabla("LOAD PATTERN DEFINITIONS", [
        fila(LoadPat="PP", DesignType="Dead", SelfWtMult=F_CONEX),
        fila(LoadPat="SCP", DesignType="Super Dead", SelfWtMult=0),
        fila(LoadPat="L", DesignType="Live", SelfWtMult=0),
        fila(LoadPat="LR", DesignType="Roof Live", SelfWtMult=0),
        fila(LoadPat="PUSH", DesignType="Other", SelfWtMult=0)])
    lin += tabla("LOAD CASE DEFINITIONS", [
        fila(Case=p, Type="LinStatic", InitialCond="Zero", DesTypeOpt="Prog Det",
             DesignType={"PP": "Dead", "SCP": "Super Dead", "L": "Live", "LR": "Roof Live", "PUSH": "Other"}[p],
             AutoType="None", RunCase=True) for p in ("PP", "SCP", "L", "LR", "PUSH")])
    lin += tabla("CASE - STATIC 1 - LOAD ASSIGNMENTS",
                 [fila(Case=p, LoadType="Load pattern", LoadName=p, LoadSF=1) for p in ("PP", "SCP", "L", "LR", "PUSH")])
    vig = cargas_eje2()
    dl = []
    for k in range(1, N_PISOS + 1):
        for pat, segs in vig[k].cargas:
            for ra, wa, rb, wb in segs:
                dl.append(fila(Frame=vig_id[k], LoadPat=pat, CoordSys="GLOBAL", Type="Force", Dir="Gravity",
                               DistType="RelDist", RelDistA=ra, RelDistB=rb, AbsDistA=ra * YS[1], AbsDistB=rb * YS[1],
                               FOverLA=wa, FOverLB=wb))
    lin += tabla("FRAME LOADS - DISTRIBUTED", dl)
    jl = []
    # Reacciones de vigas X (paños 1-2 y 2-3) que llegan a las columnas del eje 2, separadas por patrón
    for k in range(1, N_PISOS + 1):
        for c_, iy in enumerate((0, 1)):
            porpat = {}
            for (a, bb) in ((0, 1), (1, 2)):
                bx = barra_entre(NUDO_ID[(a, iy, k)], NUDO_ID[(bb, iy, k)])
                Lx = largo(bx)
                porpat["PP"] = porpat.get("PP", 0.0) + GAMMA_ACERO * SECCIONES[bx.sec]["A"] * F_CONEX * Lx / 2
                for pat, segs in bx.cargas:
                    porpat[pat] = porpat.get(pat, 0.0) + sum((wa + wb) / 2 * (rb - ra) * Lx for ra, wa, rb, wb in segs) / 2
            for pat, P in porpat.items():
                jl.append(fila(Joint=10 * k + c_ + 1, LoadPat=pat, CoordSys="GLOBAL", F1=0, F2=0, F3=-P, M1=0, M2=0, M3=0))
    for k in range(1, N_PISOS + 1):
        jl.append(fila(Joint=10 * k + 1, LoadPat="PUSH", CoordSys="GLOBAL", F1=0, F2=po["fk"][k - 1], F3=0,
                       M1=0, M2=0, M3=0))
    lin += tabla("JOINT LOADS - FORCE", jl)
    lin += ["END TABLE DATA"]
    with open(ruta, "w", newline="\r\n") as fh:
        fh.write("\n".join(lin) + "\n")


# =============================================================================
# 12. MAIN
# =============================================================================
if __name__ == "__main__":
    import os
    import subprocess
    import memoria_latex
    if len(sys.argv) == 1:      # sin argumento: generar ambas configuraciones
        for cat in ("AISC", "CL"):
            print("=" * 20, cat, "=" * 20)
            subprocess.run([sys.executable, os.path.abspath(__file__), cat], check=True)
        sys.exit(0)
    aqui = os.path.dirname(os.path.abspath(__file__))
    U, res, pos = analizar_3d()
    env = envolventes(res)
    resumen = resumen_por_tipo(env)
    po = pushover_eje2()
    trac = disenar_traccion(resumen)
    escribir_s2k_3d(os.path.join(aqui, "G5_Edificio3D%s.s2k" % SUF))
    escribir_s2k_eje2(os.path.join(aqui, "G5_Eje2_Pushover%s.s2k" % SUF), po)

    # reacciones verticales totales (chequeo de equilibrio)
    ids = sorted(NUDOS)
    # desplazamientos de CM por piso para EX, EY
    desp = []
    for k in range(1, N_PISOS + 1):
        m = NUDO_CM[k]
        desp.append(dict(k=k, ux=float(U["EX"][6 * pos[m]]), uy=float(U["EY"][6 * pos[m] + 1])))
    datos = dict(niveles=NIVELES, P_total=P_TOTAL, Q=Q_BASAL, resumen=resumen, desp=desp,
                 pushover={k: v for k, v in po.items() if k not in ("Pcol", "elems2d")},
                 traccion=trac, peso_acero_tipo=peso_acero_tipo)
    with open(os.path.join(aqui, "resultados_G5%s.json" % SUF), "w") as fh:
        json.dump(datos, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    memoria_latex.escribir(os.path.join(aqui, "Memoria_Tarea1_G5%s.tex" % SUF), globals(), datos, po, trac, resumen, env)
    print("P total = %.1f tonf ; Q = %.1f tonf" % (P_TOTAL, Q_BASAL))
    for n in NIVELES:
        print("  nivel %d  P=%.1f  A=%.4f  F=%.2f" % (n["k"], n["P"], n["A"], n["F"]))
    for tipo, d in resumen.items():
        print(tipo, {k: ("%.1f" % v[0], v[1], v[2]) for k, v in d.items()})
    print("Pushover: F_1a fluencia = %.1f  F_colapso = %.1f  F_cinem = %.1f (%s)" %
          (po["lam_y"], po["F_colapso"], po["mecanismos"][0]["F"], po["mecanismos"][0]["nombre"]))
    for ev in po["eventos"]:
        print("   ", ev)
    print("Traccion:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in trac.items()})
    print("Desp:", desp)
