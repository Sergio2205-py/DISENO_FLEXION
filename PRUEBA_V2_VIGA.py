import math
import pandas as pd

def calculoFlexion(
    b,
    h,
    fc,
    fy,
    Es,
    Ecu,
    phiFlexion,
    acero,
    r
):
    # Cálculo del factor beta1
    if fc <= 280:
        beta1 = 0.85
    elif fc <= 560:
        beta1 = round(1.05 - 0.714 * (fc / 1000), 3)
    else:
        beta1 = 0.65

    d = h - r

    eps_y = fy / Es
    cb = d * Ecu / (Ecu + eps_y)
    ab = beta1 * cb 

    aceroMinimo = 0.7 * ((math.sqrt(fc)) / fy) * b * d
    aceroBalanceado = 0.85 * fc * b * ab / fy
    aceroMaximo = 0.75 * aceroBalanceado

        if acero < aceroBalanceado:
        a = acero * fy / (0.85 * fc * b)
        c = a / beta1
        
        Mn = (acero * fy * (d - a/2)) / (1000 * 100)
        phiMn = phiFlexion * Mn
        
        Cc = (0.85 * fc * b * a)/1000
        T = acero * fy / 1000  # en tonf
    else:
        # Caso: acero no fluye
        A = (0.85 * fc) / (Ecu * Es * (acero / (b * d)))
        B = d
        C = -beta1 * d * d

        a = (-B + math.sqrt(B * B - 4 * A * C)) / (2 * A)
        c = a / beta1
        Mn = (0.85 * fc * a * b * (d - a / 2)) / (1000 * 100)
        phiMn = phiFlexion * Mn
        Cc = (0.85 * fc * b * a)/10**3    # ya estaba aquí
        defAs = Ecu * (d - c) / c
        fs = Es * defAs
        T = acero * fs / 1000

    # Tipo de falla
    if round(acero, 2) < round(aceroBalanceado, 2):
        tipoFalla = "Tracción"
    elif round(acero, 2) > round(aceroBalanceado, 2):
        tipoFalla = "Compresión"
    else:
        tipoFalla = "Balanceada"

    defAs = round(Ecu * (d - c) / c, 6)
    
    resultado = {
        "beta1": beta1,
        "d": d,
        "aceroMinimo_val": aceroMinimo,
        "aceroBalanceado_val": aceroBalanceado,
        "aceroMaximo_val": aceroMaximo,
        "a_val": a,
        "c_val": c,
        "cb_val": cb,
        "cb": f"{cb:.2f} cm",
        "Mn_val": Mn,
        "phiMn_val": phiMn,
        "defAs": defAs,
        "Cc_val": Cc,
        "aceroMinimo": f"{aceroMinimo:.2f} cm²",
        "aceroBalanceado": f"{aceroBalanceado:.2f} cm²",
        "aceroMaximo": f"{aceroMaximo:.2f} cm²",
        "a": f"{a:.2f} cm",
        "c": f"{c:.2f} cm",
        "Mn": f"{Mn:.2f} ton·m",
        "phiMn": f"{phiMn:.2f} ton·m",
        "tipoFalla": tipoFalla,
        "Cc": f"{Cc:.2f} tonf",
        "T_val": T,
        "T": f"{T:.2f} tonf",
        "fs": fs,
    }

    return resultado

def acero_requerido_flexion_simple_formula(
    b, h, r, fc, fy, phi, Mu
):
    """
    Calcula As requerido para sección rectangular
    reforzada simple a partir de Mu.
    """

    d = h - r

    Mn = Mu / phi   # ton-m

    # pasar a kgf-cm
    Mn = Mn * 1000 * 100

    A = fy**2 / (2 * 0.85 * fc * b)
    B = -fy * d
    C = Mn

    discriminante = B**2 - 4*A*C

    if discriminante < 0:
        return None, None

    As1 = (-B - discriminante**0.5) / (2*A)
    As2 = (-B + discriminante**0.5) / (2*A)

    As = min(x for x in [As1, As2] if x > 0)

    a = As * fy / (0.85 * fc * b)

    return round(As,2), round(a,2)

tablaAceros = pd.DataFrame(
    {
        "Diametro": [
            "6mm",
            '1/4"',
            "8mm",
            '3/8"',
            "12mm",
            '1/2"',
            '5/8"',
            '3/4"',
            '1"',
            '1 3/8"',
        ],
        "Área(cm2)": [0.28, 0.32, 0.5, 0.71, 1.13, 1.29, 2, 2.84, 5.1, 10.06],
    }
)
#tablaAceros = tablaAceros.set_index("Diametro")

def areaAs (numero, diametro):
    return numero * tablaAceros.loc[tablaAceros['Diametro'] == diametro, "Área(cm2)"].values[0] 


def ancho_minimo_acero(
    grupos_acero,
    recubrimiento=4.0,
    diam_estribo=1.0,
    sep_min_aci=2.54
):
    """
    Calcula el ancho mínimo necesario para alojar el acero en una capa
    grupos_acero: lista de tuplas (n_barras, diametro_str)
    """

    diametros_cm = []

    for n, diametro_str in grupos_acero:
        if n == 0:
            continue

        area = tablaAceros.loc[
            tablaAceros["Diametro"] == diametro_str,
            "Área(cm2)"
        ].values[0]

        db_cm = (4 * area / 3.1416) ** 0.5
        diametros_cm.append((n, db_cm))

    if not diametros_cm:
        return 0.0, 0

    # Ancho ocupado por barras
    ancho_barras = sum(n * db for n, db in diametros_cm)

    # Separaciones entre barras
    n_total = sum(n for n, _ in diametros_cm)
    ancho_separaciones = (n_total - 1) * sep_min_aci

    # Ancho total requerido
    b_min = (
        ancho_barras
        + ancho_separaciones
        + 2 * recubrimiento
        + 2 * diam_estribo
    )

    return round(b_min, 2), n_total


def calculoFlexionDoble(
    b,
    h,
    fc,
    fy,
    Es,
    Ecu,
    phiFlexion,
    As_trac,
    As_comp,
    r_trac,
    r_comp,
):
    import math

    # -------------------------------
    # PARÁMETROS DEL CONCRETO
    # -------------------------------
    if fc <= 280:
        beta1 = 0.85
    elif fc <= 560:
        beta1 = round(1.05 - 0.714 * (fc / 1000), 3)
    else:
        beta1 = 0.65

  # ---------------------------------
# DEFINICIÓN CORRECTA DE DISTANCIAS
# ---------------------------------
    # d = distancia desde fibra comprimida hasta acero a tracción
    d_trac = h - r_trac
    # d' = distancia desde fibra comprimida hasta acero en compresión
    d_comp = r_comp
    # -------------------------------
    # ECUACIÓN CUADRÁTICA EN c
    # -------------------------------
    A = 0.85 * fc * beta1 * b
    B = As_comp * Es * Ecu - As_trac * fy
    C = -As_comp * Es * Ecu * d_comp

    discriminante = B**2 - 4 * A * C
    c = (-B + math.sqrt(discriminante)) / (2 * A)

    a = beta1 * c
    
    # DEFORMACIONES
    eps_s  = Ecu * (d_trac - c) / c
    eps_sp = Ecu * (c - d_comp) / c

    # ESFUERZOS
    fs = min(Es * eps_s, fy)
    fs_p = Es * eps_sp

    # -------------------------------
    # LIMITAR ESFUERZO EN ACERO COMPRESIÓN
    # -------------------------------
    if abs(fs_p) > fy:
        fs_p = fy * (1 if fs_p > 0 else -1)
    
    # -------------------------------
    # FUERZAS (en tonf)
    # -------------------------------
    T  = As_trac * fs / 1000
    Cc = 0.85 * fc * b * a / 1000
    Cs = As_comp * abs(fs_p) / 1000

    # -------------------------------
    # MOMENTO NOMINAL
    # -------------------------------
    Mn = (Cc * (d_trac - a / 2) + Cs * (d_trac - d_comp)) / (1000 * 100)
    phiMn = phiFlexion * Mn

    # -------------------------------
    # ACEROS DE REFERENCIA
    # -------------------------------
    As_min = 0.7 * (fc ** 0.5) / fy * b * d_trac

    eps_y = fy / Es
    cb = (Ecu / (Ecu + eps_y)) * d_trac

    ab = beta1 * cb

    eps_sp_bal = Ecu * (cb - d_comp) / cb
    fs_p_bal = min(Es * eps_sp_bal, fy)

    Cc_bal = 0.85 * fc * b * ab
    Cs_bal = As_comp * fs_p_bal

    As_bal = (Cc_bal + Cs_bal) / fy
    As_max = 0.75 * As_bal
    # -------------------------------
    # TIPO DE FALLA
    # -------------------------------
    tipoFalla = "Tracción" if eps_s >= 0.005 else "Compresión"

    return {
        "beta1": beta1,
        "d": d_trac,
        "aceroMinimo_val": As_min,
        "aceroBalanceado_val": As_bal,
        "aceroMaximo_val": As_max,
        "a_val": a,
        "c_val": c,
        "phiMn_val": phiMn,
        "defAs": round(eps_s, 5),
        "Cc_val": Cc / 10**3,
        "Mn_val": Mn,
        "cb_val": cb,
        "cb": f"{cb:.2f} cm",
        "aceroMinimo": f"{As_min:.2f} cm²",
        "aceroBalanceado": f"{As_bal:.2f} cm²",
        "aceroMaximo": f"{As_max:.2f} cm²",
        "a": f"{a:.2f} cm",
        "c": f"{c:.2f} cm",
        "phiMn": f"{phiMn:.2f} ton·m",
        "tipoFalla": tipoFalla,
        "T_val": T,
        "T": f"{T:.2f} tonf",
        "Cs_val": Cs,
        "Cs": f"{Cs:.2f} tonf",
        "fs": fs,
        "fs_comp": fs_p,
    }

def sugerir_acero(As_req):
    mejores = []

    for _, fila in tablaAceros.iterrows():
        nombre = fila["Diametro"]
        area = fila["Área(cm2)"]

        for n in range(2, 9):   # de 2 a 8 barras
            As_total = n * area

            if As_total >= As_req:
                exceso = As_total - As_req
                mejores.append((exceso, n, nombre, As_total))

    mejores.sort()

    if mejores:
        _, n, nombre, As_total = mejores[0]
        return f"{n} barras de {nombre} = {As_total:.2f} cm²"

    return "Sin sugerencia"


def diseno_flexion_doble(
    b, h, fc, fy, Es, Ecu, phi,
    Mu, As_comp, d, d_prima
):

    Mu = Mu * 100000  # ton·m → kg·cm

    if fc <= 280:
        beta1 = 0.85
    elif fc <= 560:
        beta1 = 1.05 - 0.714*(fc/1000)
    else:
        beta1 = 0.65

    c = 1.0
    paso = 0.05

    for _ in range(1000):

        a = beta1 * c

        # deformación acero compresión
        eps_s_comp = Ecu * (c - d_prima) / c
        fs_comp = min(Es * eps_s_comp, fy)

        # fuerzas
        Cc = 0.85 * fc * b * a
        Cs = As_comp * fs_comp

        # acero requerido
        As_trac = (Cc + Cs) / fy

        # momento
        Mn = Cc * (d - a/2) + Cs * (d - d_prima)

        if phi * Mn >= Mu:
            break

        c += paso

    return {
        "c": c,
        "a": a,
        "As_trac": As_trac,
        "Mn": Mn / 100000
    }
