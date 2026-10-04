# -*- coding: utf-8 -*-
"""
Crea los modelos del Grupo 5 DIRECTAMENTE en SAP2000 (v27, también v20+) usando la API (OAPI)
y los guarda como .sdb. No depende de importar archivos .s2k.

Requisitos (en el PC con SAP2000 instalado, Windows):
    pip install numpy comtypes
    (SAP2000 debe estar instalado y registrado; si la API no se encuentra, ejecutar
     "RegisterSAP2000.exe" como administrador desde la carpeta de instalación de SAP2000)

Uso (desde la carpeta Proyecto_Acero_G5):
    python crear_modelo_SAP_API.py CL      -> perfiles soldados chilenos (A270ES)
    python crear_modelo_SAP_API.py AISC    -> perfiles W / HSS
    Opcional: agregar  --correr  para ejecutar el análisis al final.

Genera:
    G5_Edificio3D[_CL].sdb      modelo 3D completo
    G5_Eje2_Pushover[_CL].sdb   modelo 2D del eje 2 (casos CGNL y PUSHOVER ya definidos)
"""
import os
import sys

CAT = "CL"
for a in sys.argv[1:]:
    if a.upper() in ("CL", "AISC"):
        CAT = a.upper()
CORRER = "--correr" in sys.argv
sys.argv = [sys.argv[0], CAT]          # generar_modelo lee la configuración desde sys.argv
import generar_modelo as G             # noqa: E402  (geometría, secciones y cargas)

import comtypes.client                 # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))

# Enumeraciones de la OAPI
TON_M_C = 12
MAT_STEEL = 1
LP = {"Dead": 1, "Super Dead": 2, "Live": 3, "Quake": 5, "Other": 8, "Roof Live": 11}
DIR_GRAVITY = 10
EJE_Z = 3


def chk(ret, que):
    """La OAPI devuelve 0 si todo salió bien (a veces dentro de una lista)."""
    r = ret[-1] if isinstance(ret, (list, tuple)) else ret
    if r != 0:
        print("   ADVERTENCIA: %s devolvió %s" % (que, r))
    return ret


def iniciar_sap():
    helper = comtypes.client.CreateObject("SAP2000v1.Helper")
    helper = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    try:
        sap = helper.GetObject("CSI.SAP2000.API.SapObject")    # usa un SAP2000 ya abierto
        print("Conectado a una sesión abierta de SAP2000.")
    except Exception:
        sap = helper.CreateObjectProgID("CSI.SAP2000.API.SapObject")
        sap.ApplicationStart()
        print("SAP2000 iniciado.")
    return sap


def nuevo_modelo(sap):
    m = sap.SapModel
    chk(m.InitializeNewModel(TON_M_C), "InitializeNewModel")
    chk(m.File.NewBlank(), "NewBlank")
    return m


def materiales_y_secciones(m, nombres):
    mats = {}
    for mat in (G.MAT_PERFIL, G.MAT_RIOSTRA):
        mats[mat[0]] = mat
    for nom, fy, fu, _, _ in mats.values():
        chk(m.PropMaterial.SetMaterial(nom, MAT_STEEL), "SetMaterial " + nom)
        chk(m.PropMaterial.SetMPIsotropic(nom, G.E, G.NU, 1.17e-5), "SetMPIsotropic")
        chk(m.PropMaterial.SetWeightAndMass(nom, 1, G.GAMMA_ACERO), "SetWeightAndMass")
        try:
            m.PropMaterial.SetOSteel_1(nom, fy * G.MPA, fu * G.MPA, 1.1 * fy * G.MPA, 1.1 * fu * G.MPA,
                                       1, 1, 0.02, 0.14, 0.20, -0.1)
        except Exception as e:
            print("   (SetOSteel_1 no disponible: %s)" % e)
    for nom in nombres:
        s = G.SECCIONES[nom]
        if s["tipo"] == "I":
            chk(m.PropFrame.SetISection(nom, s["mat"], s["t3"], s["t2"], s["tf"], s["tw"], s["t2"], s["tf"]),
                "SetISection " + nom)
        else:
            try:
                chk(m.PropFrame.SetTube(nom, s["mat"], s["t3"], s["t2"], s["tf"], s["tw"]), "SetTube " + nom)
            except Exception:
                chk(m.PropFrame.SetTube_1(nom, s["mat"], s["t3"], s["t2"], s["tf"], s["tw"], 0.0), "SetTube_1 " + nom)


def patrones(m, lista):
    for nom, tipo, sw in lista:
        chk(m.LoadPatterns.Add(nom, LP[tipo], sw, True), "LoadPattern " + nom)
    try:
        m.LoadPatterns.Delete("DEAD")       # patrón por defecto del modelo en blanco
    except Exception:
        pass


def cargas_distribuidas(m, nombre_sap, barra, Lb, factor_pat=None):
    for pat, segs in barra.cargas:
        for ra, wa, rb, wb in segs:
            chk(m.FrameObj.SetLoadDistributed(nombre_sap, pat, 1, DIR_GRAVITY, ra, rb, float(wa), float(wb),
                                              "Global", True, False), "Carga distribuida")


# =============================================================================
def modelo_3d(sap):
    m = nuevo_modelo(sap)
    print("Creando modelo 3D (%s)..." % CAT)
    secs = [G.SEC_COL, G.SEC_VIGA_X, G.SEC_VIGA_Y, G.SEC_RIO_X, G.SEC_RIO_Y]
    materiales_y_secciones(m, secs)
    # nudos
    for n in sorted(G.NUDOS.values(), key=lambda q: q.id):
        chk(m.PointObj.AddCartesian(n.x, n.y, n.z, "", str(n.id)), "Nudo %d" % n.id)
    # barras
    for b in G.BARRAS:
        chk(m.FrameObj.AddByPoint(str(b.i), str(b.j), "", b.sec, str(b.id)), "Barra %d" % b.id)
        if b.tipo == "COL":
            chk(m.FrameObj.SetLocalAxes(str(b.id), 90.0), "Ejes locales")
        if G.es_riostra(b):
            ii = [False, False, False, False, True, True]
            jj = [False, False, False, False, True, True]
            chk(m.FrameObj.SetReleases(str(b.id), ii, jj, [0.0] * 6, [0.0] * 6), "Liberaciones")
    # apoyos y diafragmas
    for n in G.NUDOS.values():
        if n.nivel == 0:
            chk(m.PointObj.SetRestraint(str(n.id), [True] * 6), "Empotramiento")
    for k in range(1, G.N_PISOS + 1):
        chk(m.ConstraintDef.SetDiaphragm("DIAF%d" % k, EJE_Z), "Diafragma")
    for n in G.NUDOS.values():
        if n.nivel > 0:
            chk(m.PointObj.SetConstraint(str(n.id), "DIAF%d" % n.nivel), "Asignar diafragma")
    for nid in G.NUDO_CM.values():
        chk(m.PointObj.SetRestraint(str(nid), [False, False, True, True, True, False]), "Restricción CM")
    # cargas
    patrones(m, [("PP", "Dead", G.F_CONEX), ("SCP", "Super Dead", 0), ("L", "Live", 0),
                 ("LR", "Roof Live", 0), ("EX", "Quake", 0), ("EY", "Quake", 0)])
    for b in G.BARRAS:
        if b.cargas:
            cargas_distribuidas(m, str(b.id), b, G.largo(b))
    for nl in G.NIVELES:
        nid = str(G.NUDO_CM[nl["k"]])
        chk(m.PointObj.SetLoadForce(nid, "EX", [float(nl["F"]), 0.0, 0.0, 0.0, 0.0, 0.0], False), "Carga EX")
        chk(m.PointObj.SetLoadForce(nid, "EY", [0.0, float(nl["F"]), 0.0, 0.0, 0.0, 0.0], False), "Carga EY")
    # combinaciones NCh3171
    for cname, _, fac in G.COMBOS:
        chk(m.RespCombo.Add(cname, 0), "Combo " + cname)
        for p, f in fac.items():
            chk(m.RespCombo.SetCaseList(cname, 0, p, f), "Combo caso")
    chk(m.RespCombo.Add("ENVOLVENTE", 1), "Envolvente")
    for cname, _, _ in G.COMBOS:
        chk(m.RespCombo.SetCaseList("ENVOLVENTE", 1, cname, 1.0), "Envolvente caso")
    # grupos
    grupos = {"COLUMNAS": "COL", "VIGAS_X": "VX", "VIGAS_Y": "VY", "RIOSTRAS_X": "RX", "RIOSTRAS_Y": "RY"}
    for g, t in grupos.items():
        m.GroupDef.SetGroup(g)
        for b in G.BARRAS:
            if b.tipo == t:
                m.FrameObj.SetGroupAssign(str(b.id), g)
    m.GroupDef.SetGroup("EJE_2")
    for b in G.BARRAS:
        if b.tipo in ("COL", "VY") and b.eje[0] == "2":
            m.FrameObj.SetGroupAssign(str(b.id), "EJE_2")
    m.View.RefreshView(0, False)
    ruta = os.path.join(AQUI, "G5_Edificio3D%s.sdb" % G.SUF)
    chk(m.File.Save(ruta), "Guardar")
    print("Guardado:", ruta)
    if CORRER:
        chk(m.Analyze.RunAnalysis(), "RunAnalysis")
        print("Análisis 3D ejecutado.")


# =============================================================================
def modelo_eje2(sap):
    m = nuevo_modelo(sap)
    print("Creando modelo 2D del eje 2 (%s)..." % CAT)
    try:
        m.Analyze.SetActiveDOF([False, True, True, True, False, False])   # pórtico plano YZ
    except Exception:
        pass
    materiales_y_secciones(m, [G.SEC_COL, G.SEC_VIGA_Y])
    x = G.XS[1]
    for k in range(G.N_PISOS + 1):
        for c, y in enumerate(G.YS):
            chk(m.PointObj.AddCartesian(x, y, G.ZS[k], "", str(10 * k + c + 1)), "Nudo")
    fid = 0
    for k in range(1, G.N_PISOS + 1):
        for c in range(2):
            fid += 1
            chk(m.FrameObj.AddByPoint(str(10 * (k - 1) + c + 1), str(10 * k + c + 1), "", G.SEC_COL, str(fid)), "Columna")
            chk(m.FrameObj.SetLocalAxes(str(fid), 90.0), "Ejes locales")
    vig = G.cargas_eje2()
    for k in range(1, G.N_PISOS + 1):
        fid += 1
        chk(m.FrameObj.AddByPoint(str(10 * k + 1), str(10 * k + 2), "", G.SEC_VIGA_Y, str(fid)), "Viga")
        cargas_distribuidas(m, str(fid), vig[k], G.YS[1])
    for c in range(2):
        chk(m.PointObj.SetRestraint(str(c + 1), [True] * 6), "Empotramiento")
    patrones(m, [("PP", "Dead", G.F_CONEX), ("SCP", "Super Dead", 0), ("L", "Live", 0),
                 ("LR", "Roof Live", 0), ("PUSH", "Other", 0)])
    # reacciones de vigas X que llegan a las columnas del eje 2
    for k in range(1, G.N_PISOS + 1):
        for c, iy in enumerate((0, 1)):
            porpat = {}
            for (a, bb) in ((0, 1), (1, 2)):
                bx = G.barra_entre(G.NUDO_ID[(a, iy, k)], G.NUDO_ID[(bb, iy, k)])
                Lx = G.largo(bx)
                porpat["PP"] = porpat.get("PP", 0.0) + G.GAMMA_ACERO * G.SECCIONES[bx.sec]["A"] * G.F_CONEX * Lx / 2
                for pat, segs in bx.cargas:
                    porpat[pat] = porpat.get(pat, 0.0) + sum((wa + wb) / 2 * (rb - ra) * Lx
                                                             for ra, wa, rb, wb in segs) / 2
            for pat, P in porpat.items():
                chk(m.PointObj.SetLoadForce(str(10 * k + c + 1), pat, [0.0, 0.0, -float(P), 0.0, 0.0, 0.0], False), "Carga columna")
    po_fk = [n["F"] / G.Q_BASAL for n in G.NIVELES]
    for k in range(1, G.N_PISOS + 1):
        chk(m.PointObj.SetLoadForce(str(10 * k + 1), "PUSH", [0.0, float(po_fk[k - 1]), 0.0, 0.0, 0.0, 0.0], False), "Carga PUSH")
    # casos no lineales (CGNL y PUSHOVER)
    try:
        nl = m.LoadCases.StaticNonlinear
        chk(nl.SetCase("CGNL"), "CGNL")
        chk(nl.SetLoads("CGNL", 4, ["Load"] * 4, ["PP", "SCP", "L", "LR"], [1.0, 1.0, 0.25, 0.25]), "CGNL cargas")
        chk(nl.SetCase("PUSHOVER"), "PUSHOVER")
        chk(nl.SetInitialCase("PUSHOVER", "CGNL"), "PUSHOVER inicial")
        chk(nl.SetLoads("PUSHOVER", 1, ["Load"], ["PUSH"], [1.0]), "PUSHOVER cargas")
        # control de desplazamiento: nudo de techo 51, dirección U2, 0.80 m
        chk(nl.SetLoadApplication("PUSHOVER", 2, 2, 0.80, 1, 2, str(10 * G.N_PISOS + 1), ""), "Control despl.")
        chk(nl.SetResultsSaved("PUSHOVER", True, 10, 100, True), "Pasos guardados")
    except Exception as e:
        print("   No se pudieron crear los casos no lineales por API (%s). Crearlos a mano (Anexo B)." % e)
    m.View.RefreshView(0, False)
    ruta = os.path.join(AQUI, "G5_Eje2_Pushover%s.sdb" % G.SUF)
    chk(m.File.Save(ruta), "Guardar")
    print("Guardado:", ruta)
    print("FALTA (a mano): asignar rótulas M3 a todas las barras en 0 y 1 (Assign > Frame > Hinges),")
    print("  Mp columna = %.1f tonf-m, Mp viga = %.1f tonf-m (ver Anexo B de la memoria)." %
          (G.SECCIONES[G.SEC_COL]["Z33"] * G.FY_W * G.MPA, G.SECCIONES[G.SEC_VIGA_Y]["Z33"] * G.FY_W * G.MPA))


if __name__ == "__main__":
    sap = iniciar_sap()
    modelo_3d(sap)
    modelo_eje2(sap)
    print("Listo. Los modelos quedaron guardados como .sdb en", AQUI)
