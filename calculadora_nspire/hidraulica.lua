-- =====================================================================
--  HIDRAULICA APLICADA CIV-346  --  Solucionador para TI-Nspire CX II CAS
--  Clases 01-12: aguas subterraneas, drenaje y drenaje superficial
--
--  USO
--   * Menu principal: flechas + [enter] (o numero 1-9). [esc] = volver.
--   * Formulario: flechas arriba/abajo para elegir dato, escribir el
--     valor (se aceptan expresiones: 2E-4, 3*10^-5, sqrt(2), pi, ln(3)).
--     [del] borra, [clear] limpia el campo.
--   * ECUACIONES (marcadas "="): deje VACIO el dato que quiere calcular
--     y presione [enter]; el programa despeja numericamente esa variable.
--   * PROCEDIMIENTOS (marcados ">"): complete los datos y [enter].
--     Listas: valores separados por coma  ej. 10,20,35
--   * GRAFICOS: los ejercicios marcados [graf] muestran graficos al
--     presionar [enter] en el resultado. En el grafico: flechas <- ->
--     mueven el cursor (lee x,y), flechas arriba/abajo cambian de serie,
--     [tab] o [enter] pasa al siguiente grafico, [esc] vuelve.
--   * Seccion "GRAFICOS y simuladores": pruebas de gasto/recuperacion,
--     escalonada, curva tipo interactiva, conos de depresion, Hantush.
--   * Seccion "ESQUEMAS": dibujos de alcantarillas (5 casos), tipos de
--     atraviesos, cunetas, sumideros, pozos, drenes, puentes, etc.
--   * Unidades: use un sistema coherente (SI: m, s, m3/s) salvo que el
--     campo indique otra unidad.
-- =====================================================================
platform.apiLevel = "2.2"

-- ---------------------------------------------------------------------
--  Utilidades numericas
-- ---------------------------------------------------------------------
local G = 9.81
local CJ = 2.24          -- constante de Cooper-Jacob usada en clases (2.25 exacto)
local pi, ln, sqrt, exp, abs = math.pi, math.log, math.sqrt, math.exp, math.abs

local function ok(x) return type(x) == "number" and x == x and x ~= 1/0 and x ~= -1/0 end
local function fmt(x)
  if type(x) ~= "number" then return tostring(x) end
  if not ok(x) then return "indef." end
  if x ~= 0 and (abs(x) >= 1e6 or abs(x) < 1e-3) then return string.format("%.4e", x) end
  return string.format("%.5g", x)
end
local function cosh(x) return (exp(x) + exp(-x)) / 2 end
local function acos(x) return math.acos(x) end

-- evaluacion de expresiones escritas por el usuario
local env = { pi = pi, e = exp(1), sqrt = sqrt, ln = ln, log = function(x) return ln(x)/ln(10) end,
  exp = exp, sin = math.sin, cos = math.cos, tan = math.tan, abs = abs,
  asin = math.asin, acos = math.acos, atan = math.atan }
local function evalstr(s)
  s = s:gsub("\226\136\146", "-"):gsub("\225\180\135", "E"):gsub("\207\128", "pi")
       :gsub("\195\151", "*"):gsub("\195\183", "/"):gsub("\226\136\154", "sqrt")
       :gsub("\194\178", "^2"):gsub("\226\132\175", "e")
  if s:match("^%s*$") then return nil end
  local f
  if setfenv then
    f = loadstring("return " .. s)
    if f then setfenv(f, env) end
  else
    f = load("return " .. s, "expr", "t", env)
  end
  if f then
    local okc, v = pcall(f)
    if okc and type(v) == "number" then return v end
  end
  if math.eval then
    local okc, v = pcall(math.eval, s)
    if okc and type(v) == "number" then return v end
  end
  error("valor invalido: " .. s)
end
local function evallist(s)
  local t = {}
  for part in (s .. ","):gmatch("([^,;]*)[,;]") do
    if not part:match("^%s*$") then t[#t + 1] = evalstr(part) end
  end
  return t
end

-- raices de f en un intervalo [a,b] (escaneo + biseccion). Devuelve lista.
local function roots(f, a, b, n, logscale)
  n = n or 600
  local res = {}
  local function xat(i)
    if logscale then return a * (b / a) ^ (i / n) end
    return a + (b - a) * i / n
  end
  local x0, f0 = xat(0), nil
  local okc, v = pcall(f, x0); if okc and ok(v) then f0 = v end
  for i = 1, n do
    local x1 = xat(i)
    local okc1, f1 = pcall(f, x1)
    if not (okc1 and ok(f1)) then f1 = nil end
    if f0 and f1 then
      if f0 == 0 then res[#res + 1] = x0
      elseif f0 * f1 < 0 then
        local lo, hi, flo = x0, x1, f0
        for _ = 1, 200 do
          local mid = logscale and sqrt(lo * hi) or (lo + hi) / 2
          local fm = f(mid)
          if not ok(fm) then break end
          if flo * fm <= 0 then hi = mid else lo, flo = mid, fm end
          if abs(hi - lo) <= 1e-13 * abs(hi) then break end
        end
        local r = (lo + hi) / 2
        -- descarta cambios de signo por discontinuidad
        local fr = f(r)
        if ok(fr) and abs(fr) <= 1e-6 * (abs(f0) + abs(f1)) + 1e-12 then res[#res + 1] = r end
      end
    end
    x0, f0 = x1, f1
  end
  return res
end
-- busca raices en todo el rango positivo (y negativo si se pide)
local function solveAll(f, neg)
  local r = roots(f, 1e-10, 1e10, 1200, true)
  if #r == 0 then
    -- sin raiz positiva: probar cero y valores negativos
    local okz, f0 = pcall(f, 0)
    if okz and ok(f0) and abs(f0) < 1e-12 then r[1] = 0 end
    neg = true
  end
  if neg then
    local rn = roots(function(x) return f(-x) end, 1e-10, 1e10, 1200, true)
    for _, x in ipairs(rn) do r[#r + 1] = -x end
  end
  return r
end
local function largest(t) local m; for _, x in ipairs(t) do if not m or x > m then m = x end end; return m end

-- ---------------------------------------------------------------------
--  Funciones especiales
-- ---------------------------------------------------------------------
-- Funcion de pozo de Theis W(u)
local function W(u)
  if u <= 0 then return 1/0 end
  if u < 1 then
    local s, term = -0.5772156649 - ln(u), 1
    for k = 1, 60 do term = term * (-u) / k; s = s - term / k end
    return s
  end
  local m = u^4 + 8.5733287401*u^3 + 18.059016973*u^2 + 8.6347608925*u + 0.2677737343
  local d = u^4 + 9.5733223454*u^3 + 25.6329561486*u^2 + 21.0996530827*u + 3.9584969228
  return m / (d * u * exp(u))
end
-- Funcion de Hantush W(u, r/B) por integracion numerica
local function Wh(u, rB)
  if rB == 0 then return W(u) end
  local smax = ln((u + 50) / u)
  local n, h = 500, nil
  h = smax / n
  local function g(s) local y = u * exp(s); return exp(-y - rB * rB / (4 * y)) end
  local sum = g(0) + g(smax)
  for i = 1, n - 1 do sum = sum + ((i % 2 == 1) and 4 or 2) * g(i * h) end
  return sum * h / 3
end
-- Kirkham F(x,y)
local function Fk(x, y)
  local s = ln(2 / (pi * x))
  for m = 1, 20000 do
    local c = 2 / (exp(4 * m * pi * y) - 1)
    if c < 1e-14 then break end
    s = s + (math.cos(m * pi * x) - math.cos(m * pi)) * c / m
  end
  return s / pi
end
-- Slichter: 1/K1 en funcion de n (interpolacion lineal)
local SL_n = {0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46}
local SL_v = {0.00187,0.01517,0.01905,0.02356,0.02878,0.03473,0.04154,0.04922,0.05789,0.06776,0.07838}
local function slichter(n)
  if n <= SL_n[1] then return SL_v[1] end
  for i = 2, #SL_n do
    if n <= SL_n[i] then
      return SL_v[i-1] + (SL_v[i] - SL_v[i-1]) * (n - SL_n[i-1]) / (SL_n[i] - SL_n[i-1])
    end
  end
  return SL_v[#SL_v]
end
-- Geometria de secciones: tipo 1 rect(b), 2 trapecio(b,k), 3 triang(k), 4 circular(D)
local function seccion(tp, b, k, D, y)
  if tp == 1 then return b*y, b + 2*y, b
  elseif tp == 2 then return (b + k*y)*y, b + 2*y*sqrt(1 + k*k), b + 2*k*y
  elseif tp == 3 then return k*y*y, 2*y*sqrt(1 + k*k), 2*k*y
  else
    if y >= D then return pi*D*D/4, pi*D, 0 end
    local th = 2 * acos(1 - 2*y/D)
    return (th - math.sin(th)) * D * D / 8, th * D / 2, math.sin(th/2) * D
  end
end
-- regresion lineal y = a*x + b
local function regresion(xs, ys)
  local n, sx, sy, sxx, sxy = #xs, 0, 0, 0, 0
  for i = 1, n do sx = sx + xs[i]; sy = sy + ys[i]; sxx = sxx + xs[i]^2; sxy = sxy + xs[i]*ys[i] end
  local a = (n*sxy - sx*sy) / (n*sxx - sx*sx)
  return a, (sy - a*sx) / n
end

-- ---------------------------------------------------------------------
--  Definicion de ejercicios
--  E(nombre, formula, campos, residuo, extra)  -> ecuacion despejable
--  P(nombre, formula, campos, calculo)         -> procedimiento
--  campo = {id, etiqueta, descripcion, opt=valor_por_defecto, list=true}
-- ---------------------------------------------------------------------
local function E(name, f, fields, res, extra, neg)
  return {kind = "E", name = name, f = f, fields = fields, res = res, extra = extra, neg = neg}
end
local function P(name, f, fields, calc)
  return {kind = "P", name = name, f = f, fields = fields, calc = calc}
end

local CATS = {}
local function cat(name, items) CATS[#CATS + 1] = {name = name, items = items} end

-- ===================== CLASE 01 =====================
cat("C01 Suelos, Darcy y K", {
  E("Ley de Darcy", "Q = K*A*i   (i = dh/L)",
    {{"Q","Q","Caudal [m3/s]"},{"K","K","Conductividad [m/s]"},{"A","A","Area [m2]"},{"i","i","Gradiente dh/L [-]"}},
    function(v) return v.Q - v.K*v.A*v.i end,
    function(v) return {"v Darcy = Q/A = " .. fmt(v.Q/v.A) .. " m/s"} end),
  E("Velocidad real", "vr = v/ne = K*i/ne",
    {{"vr","vr","Velocidad real [m/s]"},{"K","K","Conductividad [m/s]"},{"i","i","Gradiente [-]"},{"ne","ne","Porosidad efectiva"}},
    function(v) return v.vr - v.K*v.i/v.ne end),
  E("Carga hidraulica", "h = z + p/gamma",
    {{"h","h","Carga [m]"},{"z","z","Cota [m]"},{"p","p","Presion [Pa]"},{"gam","gamma","Peso especif. [N/m3] (9810)"}},
    function(v) return v.h - v.z - v.p/v.gam end, nil, true),
  E("K <-> permeab. intrinseca", "K = k*rho*g/mu",
    {{"K","K","Conductividad [m/s]"},{"k","k","Permeab. intrinseca [m2] (1 darcy=9.87E-13 m2)"},
     {"rho","rho","Densidad [kg/m3] (1000)"},{"mu","mu","Viscosidad [Pa s] (1E-3)"}},
    function(v) return v.K - v.k*v.rho*G/v.mu end),
  E("Permeametro carga cte.", "K = 4*L*Q/(pi*(h1-h2)*D^2)",
    {{"K","K","Conductividad"},{"L","L","Largo muestra"},{"Q","Q","Caudal"},{"dh","h1-h2","Perdida de carga"},{"D","D","Diametro muestra"}},
    function(v) return v.K - 4*v.L*v.Q/(pi*v.dh*v.D^2) end),
  E("Permeametro carga var.", "K = dt^2*L/(dc^2*t)*ln(h0/h)",
    {{"K","K","Conductividad"},{"dt","dt","Diametro tubo"},{"dc","dc","Diametro muestra"},{"L","L","Largo muestra"},
     {"t","t","Tiempo"},{"h0","h0","Carga inicial"},{"h","h","Carga final"}},
    function(v) return v.K - v.dt^2*v.L/(v.dc^2*v.t)*ln(v.h0/v.h) end),
  E("Hazen", "K[cm/s] = C*d10^2, d10 en cm",
    {{"K","K","Conductividad [cm/s]"},{"C","C","Coef. (40-150)"},{"d10","d10","Diametro efectivo [cm]"}},
    function(v) return v.K - v.C*v.d10^2 end),
  P("Slichter", "K = (gamma/mu)*10.0219*d10^2/K1  (1/K1 de tabla segun n)",
    {{"n","n","Porosidad (0.26-0.46)"},{"d10","d10","d10 [cm]"},{"gam","gamma","Peso esp. (g/cm2s2 -> 981 cgs)", opt="981"},
     {"mu","mu","Viscosidad (cgs: 0.01 poise)", opt="0.01"}},
    function(v)
      local iK1 = slichter(v.n)
      local K = v.gam/v.mu*10.0219*v.d10^2*iK1
      return {"1/K1 = " .. fmt(iK1), "K = " .. fmt(K) .. " (cm/s si datos cgs)"}
    end),
  E("Kozeny-Carman", "K = (rho*g/mu)*n^3/(1-n)^2*d50^2/180",
    {{"K","K","Conductividad [m/s]"},{"n","n","Porosidad"},{"d50","d50","d50 [m]"},
     {"rho","rho","Densidad (1000)"},{"mu","mu","Viscosidad (1E-3)"}},
    function(v) return v.K - v.rho*G/v.mu*v.n^3/(1-v.n)^2*v.d50^2/180 end),
  E("Metodo de 2 pozos", "K = L^2/(t*dh)",
    {{"K","K","Conductividad"},{"L","L","Distancia entre pozos"},{"t","t","Tiempo medio"},{"dh","dh","Desnivel de carga"}},
    function(v) return v.K - v.L^2/(v.t*v.dh) end),
  E("Trazador en 1 pozo", "ln(C0/C) = 4*v*t/(pi*D) ,  v = K*i",
    {{"C0","C0","Conc. inicial"},{"C","C","Conc. en t"},{"t","t","Tiempo"},{"D","D","Diametro pozo"},
     {"K","K","Conductividad"},{"i","i","Gradiente"}},
    function(v) return ln(v.C0/v.C) - 4*v.K*v.i*v.t/(pi*v.D) end),
  E("Porchet (infiltracion)", "f = R/(2(t2-t1))*ln((2h1+R)/(2h2+R))",
    {{"f","f","Tasa de infiltracion"},{"R","R","Radio del pozo"},{"dt","t2-t1","Intervalo"},{"h1","h1","Altura agua t1"},{"h2","h2","Altura agua t2"}},
    function(v) return v.f - v.R/(2*v.dt)*ln((2*v.h1 + v.R)/(2*v.h2 + v.R)) end),
  P("K equiv. estratos", "Paralelo: Keq=Sum(Ki mi)/Sum(mi) | Perpendic.: Keq=Sum(mi)/Sum(mi/Ki)",
    {{"K","K1,K2,..","Conductividades (lista)", list=true},{"m","m1,m2,..","Espesores (lista)", list=true}},
    function(v)
      local a, b, c = 0, 0, 0
      for i = 1, #v.K do a = a + v.K[i]*v.m[i]; b = b + v.m[i]; c = c + v.m[i]/v.K[i] end
      return {"Flujo PARALELO:", "  Keq = " .. fmt(a/b), "  T = Sum(Ki mi) = " .. fmt(a),
              "Flujo PERPENDICULAR:", "  Keq = " .. fmt(b/c), "Espesor total = " .. fmt(b)}
    end),
  E("Transmisibilidad", "T = K*m",
    {{"T","T","Transmisibilidad [m2/s]"},{"K","K","Conductividad"},{"m","m","Espesor saturado"}},
    function(v) return v.T - v.K*v.m end),
  E("Volumen drenado", "V = S*A*dh",
    {{"V","V","Volumen drenado"},{"S","S","Coef. almacenamiento"},{"A","A","Area"},{"dh","dh","Descenso"}},
    function(v) return v.V - v.S*v.A*v.dh end),
})

-- ===================== CLASE 02 =====================
cat("C02 Pozos permanente", {
  E("Almacenamiento especif.", "Ss = g*rho*(alfa + n*beta) ,  S = Ss*m",
    {{"S","S","Coef. almacenamiento"},{"m","m","Espesor acuifero"},{"rho","rho","Densidad (1000)"},
     {"a","alfa","Compresib. matriz [m2/N]"},{"n","n","Porosidad"},{"b","beta","Compresib. agua (4.4E-10)"}},
    function(v) return v.S - G*v.rho*(v.a + v.n*v.b)*v.m end),
  E("Thiem (confinado)", "H - h = Q/(2*pi*T)*ln(R/r)",
    {{"dh","H-h","Descenso"},{"Q","Q","Caudal"},{"T","T","Transmisibilidad K*m"},{"R","R","Radio de influencia"},{"r","r","Distancia al pozo"}},
    function(v) return v.dh - v.Q/(2*pi*v.T)*ln(v.R/v.r) end),
  E("Thiem 2 piezometros", "h2 - h1 = Q/(2*pi*T)*ln(r2/r1)",
    {{"dh","h2-h1","Diferencia de nivel"},{"Q","Q","Caudal"},{"T","T","Transmisibilidad"},{"r2","r2","Dist. lejana"},{"r1","r1","Dist. cercana"}},
    function(v) return v.dh - v.Q/(2*pi*v.T)*ln(v.r2/v.r1) end),
  E("Dupuit (libre)", "H^2 - h^2 = Q/(pi*K)*ln(R/r)",
    {{"H","H","Espesor saturado inicial"},{"h","h","Nivel a distancia r"},{"Q","Q","Caudal"},{"K","K","Conductividad"},
     {"R","R","Radio de influencia"},{"r","r","Distancia"}},
    function(v) return v.H^2 - v.h^2 - v.Q/(pi*v.K)*ln(v.R/v.r) end,
    function(v) local dh = v.H - v.h; return {"Descenso dh = " .. fmt(dh),
      "dh/H = " .. fmt(dh/v.H) .. ((dh < 0.1*v.H) and " (<0.1H: vale Thiem)" or "")} end),
  E("Libre: correccion Jacob", "dh - dh^2/(2H) = Q/(2*pi*K*H)*ln(R/r)",
    {{"dh","dh","Descenso real"},{"H","H","Espesor inicial"},{"Q","Q","Caudal"},{"K","K","Conductividad"},{"R","R","Radio infl."},{"r","r","Distancia"}},
    function(v) return v.dh - v.dh^2/(2*v.H) - v.Q/(2*pi*v.K*v.H)*ln(v.R/v.r) end),
  E("Flujo esferico", "H - h = Q/(4*pi*K)*(1/r - 1/R)",
    {{"dh","H-h","Descenso"},{"Q","Q","Caudal"},{"K","K","Conductividad"},{"r","r","Distancia"},{"R","R","Radio infl."}},
    function(v) return v.dh - v.Q/(4*pi*v.K)*(1/v.r - 1/v.R) end),
  E("Flujo semi-esferico", "H - h = Q/(2*pi*K)*(1/r - 1/R)",
    {{"dh","H-h","Descenso"},{"Q","Q","Caudal"},{"K","K","Conductividad"},{"r","r","Distancia"},{"R","R","Radio infl."}},
    function(v) return v.dh - v.Q/(2*pi*v.K)*(1/v.r - 1/v.R) end),
  E("Penetr. parcial Girinsky", "h1-h2 = Q/(2piK)*(1/b*ln(1.6b/r2) - asinh(b/r1)/b)  b/m<0.3",
    {{"dh","h1-h2","Diferencia de carga"},{"Q","Q","Caudal"},{"K","K","Conductividad"},{"b","b","Penetracion del pozo"},
     {"r1","r1","Distancia 1"},{"r2","r2","Distancia 2 (radio pozo)"}},
    function(v)
      local ash = ln(v.b/v.r1 + sqrt((v.b/v.r1)^2 + 1))
      return v.dh - v.Q/(2*pi*v.K)*(ln(1.6*v.b/v.r2) - ash)/v.b
    end),
  E("Penetr. parcial Nasberg", "h1-h2 = Q/(2piKm)[ln(r1/r2)+ash(m/r1)-ash(m/r2)-(m/b)(ash(b/r1)-ash(b/r2))]",
    {{"dh","h1-h2","Diferencia de carga"},{"Q","Q","Caudal"},{"K","K","Conductividad"},{"m","m","Espesor acuifero"},
     {"b","b","Penetracion (b/m>=0.3)"},{"r1","r1","Distancia 1"},{"r2","r2","Distancia 2"}},
    function(v)
      local function ash(x) return ln(x + sqrt(x*x + 1)) end
      return v.dh - v.Q/(2*pi*v.K*v.m)*(ln(v.r1/v.r2) + ash(v.m/v.r1) - ash(v.m/v.r2)
             - v.m/v.b*(ash(v.b/v.r1) - ash(v.b/v.r2)))
    end, nil, true),
  E("Radio infl. Sichardt", "R = 3000*dh0*sqrt(K)   (K en m/s, R y dh0 en m)",
    {{"R","R","Radio de influencia [m]"},{"dh","dh0","Descenso en el pozo [m]"},{"K","K","Conductividad [m/s]"}},
    function(v) return v.R - 3000*v.dh*sqrt(v.K) end),
})

-- ===================== CLASE 03 =====================
local function superpos(v, transiente)
  local lines, tot = {}, 0
  local X, Y, Q = v.X, v.Y, v.Q
  for i = 1, #Q do
    local r = sqrt((X[i] - v.xp)^2 + (Y[i] - v.yp)^2)
    if r <= 0 then r = v.r0 or 0.1 end
    local R = v.R
    local ti = transiente and (v.ts and v.ts[i] or 0) or 0
    if transiente then
      local dt = v.t - ti
      R = (dt > 0) and sqrt(CJ*v.T*dt/v.S) or 0
    end
    local d = 0
    if r < R then d = Q[i]/(2*pi*v.T)*ln(R/r) end
    tot = tot + d
    lines[#lines + 1] = "Pozo " .. i .. ": r=" .. fmt(r) .. " R=" .. fmt(R) .. " dh=" .. fmt(d)
  end
  lines[#lines + 1] = "DESCENSO TOTAL = " .. fmt(tot)
  if v.H and v.H > 0 then
    local h2 = v.H^2 - 2*v.H*tot
    lines[#lines + 1] = "Si es LIBRE (linealizado, T=K*H):"
    lines[#lines + 1] = "  h = sqrt(H^2-2H*dhT) = " .. (h2 > 0 and fmt(sqrt(h2)) or "seco")
  end
  return lines
end
cat("C03 Interferencia/imagen", {
  P("Superposicion (x,y)", "dhT = Sum Qi/(2piT)*ln(Ri/ri), solo si ri<R. Q>0 extrae, Q<0 inyecta (imagen)",
    {{"X","x pozos","Coordenadas x (lista)", list=true},{"Y","y pozos","Coordenadas y (lista)", list=true},
     {"Q","Q pozos","Caudales (lista, signo)", list=true},{"T","T","Transmisibilidad"},{"R","R","Radio de influencia"},
     {"xp","x punto","x del punto P"},{"yp","y punto","y del punto P"},{"H","H (libre)","Esp. inicial si es libre (0=conf.)", opt="0"}},
    function(v) return superpos(v, false) end),
  P("Imagen: 1 pozo y borde", "Borde en x=0, pozo en (b,0). Tipo 1=recarga (imagen inyecta), 2=barrera (imagen extrae)",
    {{"tp","Tipo 1/2","1 recarga lineal, 2 barrera impermeable"},{"b","b","Distancia pozo-borde"},{"Q","Q","Caudal"},
     {"T","T","Transmisibilidad"},{"R","R","Radio de influencia"},{"xp","x punto","x de P (lado del pozo, >0)"},{"yp","y punto","y de P", opt="0"}},
    function(v)
      local rr = sqrt((v.xp - v.b)^2 + v.yp^2)
      local ri = sqrt((v.xp + v.b)^2 + v.yp^2)
      local dr = (rr < v.R) and v.Q/(2*pi*v.T)*ln(v.R/rr) or 0
      local di = (ri < v.R) and v.Q/(2*pi*v.T)*ln(v.R/ri) or 0
      local tot = (v.tp == 1) and (dr - di) or (dr + di)
      return {"r real = " .. fmt(rr), "r imagen = " .. fmt(ri), "dh real = " .. fmt(dr), "dh imagen = " .. fmt(di),
        (v.tp == 1 and "Recarga: dh = dh_r - dh_i" or "Barrera: dh = dh_r + dh_i"), "DESCENSO = " .. fmt(tot),
        "Si pozos imagen afectan otro borde, use Superposicion (x,y)."}
    end),
  E("n pozos iguales", "dhT = Q/(2piT)*ln(R^n/(r1*r2*...*rn))",
    {{"dh","dhT","Descenso total"},{"Q","Q","Caudal de cada pozo"},{"T","T","Transmisibilidad"},{"R","R","Radio infl."},
     {"n","n","Numero de pozos"},{"P","r1*r2*..","Producto de distancias"}},
    function(v) return v.dh - v.Q/(2*pi*v.T)*(v.n*ln(v.R) - ln(v.P)) end),
  P("Pozo en varias napas", "Q = 2pi/ln(R/r) * Sum(dhi*Ki*mi)",
    {{"dh","dh1,dh2,..","Descensos en cada napa (lista)", list=true},{"K","K1,K2,..","Conductividades", list=true},
     {"m","m1,m2,..","Espesores", list=true},{"R","R","Radio de influencia"},{"r","r","Radio del pozo"}},
    function(v)
      local s = 0
      for i = 1, #v.K do s = s + v.dh[i]*v.K[i]*v.m[i] end
      local Q = 2*pi/ln(v.R/v.r)*s
      local out = {"Q total = " .. fmt(Q)}
      for i = 1, #v.K do out[#out+1] = "Aporte napa " .. i .. " = " .. fmt(2*pi/ln(v.R/v.r)*v.dh[i]*v.K[i]*v.m[i]) end
      return out
    end),
  E("Napas: condicion inicial", "dhT = dh1*(T1+T2)/T2",
    {{"dhT","dhT","Diferencia entre niveles piezom."},{"dh1","dh1","Descenso napa 1"},{"T1","T1","Transm. napa 1"},{"T2","T2","Transm. napa 2"}},
    function(v) return v.dhT - v.dh1*(v.T1 + v.T2)/v.T2 end,
    function(v) return {"dh2 = dhT - dh1 = " .. fmt(v.dhT - v.dh1)} end),
  E("Napa 1 libre + Q0", "Q = Q0 + 2pi*K2*m2*dh2/ln(R/r0)",
    {{"Q","Q","Caudal total"},{"Q0","Q0","Aporte cte. napa 1"},{"K2","K2","Conductividad napa 2"},{"m2","m2","Espesor napa 2"},
     {"dh2","dh2","Descenso napa 2"},{"R","R","Radio infl."},{"r0","r0","Radio pozo"}},
    function(v) return v.Q - v.Q0 - 2*pi*v.K2*v.m2*v.dh2/ln(v.R/v.r0) end),
})

-- ===================== CLASE 04 =====================
cat("C04 Redes de flujo", {
  E("Caudal red de flujo", "q = K*h*(Nt/Nd)*(a/l)  [por unidad de ancho]",
    {{"q","q","Caudal por metro [m2/s]"},{"K","K","Conductividad (o Keq)"},{"h","h","Perdida de carga total"},
     {"Nt","Nt","N. canales de flujo"},{"Nd","Nd","N. caidas de potencial"},{"al","a/l","Relacion de forma (1=cuadrados)", opt="1"}},
    function(v) return v.q - v.K*v.h*v.Nt/v.Nd*v.al end,
    function(v) return {"Caida por equipotencial dh = h/Nd = " .. fmt(v.h/v.Nd)} end),
  P("Medio anisotropico", "a = sqrt(Kz/Kx), x' = a*x, Keq = sqrt(Kx*Kz)",
    {{"Kx","Kx","Conductividad horizontal"},{"Kz","Kz","Conductividad vertical"},{"x","x","Distancia real a escalar", opt="1"}},
    function(v)
      local a = sqrt(v.Kz/v.Kx)
      return {"a = sqrt(Kz/Kx) = " .. fmt(a), "Keq = " .. fmt(sqrt(v.Kx*v.Kz)), "x' = a*x = " .. fmt(a*v.x)}
    end),
  E("Refraccion (2 medios)", "K1/K2 = tan(alfa)/tan(beta)  (angulos en grados)",
    {{"K1","K1","Conductividad medio 1"},{"K2","K2","Conductividad medio 2"},{"a","alfa","Angulo medio 1 [grados]"},{"b","beta","Angulo medio 2 [grados]"}},
    function(v) return v.K1/v.K2 - math.tan(v.a*pi/180)/math.tan(v.b*pi/180) end),
})

-- ===================== CLASE 05 =====================
cat("C05 Costero/impermanente", {
  E("Ghyben-Herzberg", "hs = hd*gd/(gs-gd)   (gs/gd=1.025 -> hs=40hd)",
    {{"hs","hs","Prof. interfaz bajo nivel mar"},{"hd","hd","Altura agua dulce sobre mar"},
     {"gd","gd","Peso esp. dulce (1000)", opt="1000"},{"gs","gs","Peso esp. salada (1025)", opt="1025"}},
    function(v) return v.hs - v.hd*v.gd/(v.gs - v.gd) end),
  P("Isla circular (recarga)", "hd^2 = N*(Ri^2 - r^2)/(2K(1+delta)), delta=gd/(gs-gd)",
    {{"Ri","Radio isla","Radio de la isla [m]"},{"K","K","Conductividad [m/s]"},{"N","N","Recarga [m/s] (mm/dia /8.64E7)"},
     {"r","r","Dist. al centro (0=centro)", opt="0"},{"x","dist. costa","Distancia desde la costa (opcional)", opt="-1"},
     {"gs","gs/gd","Densidad relativa salada", opt="1.025"}},
    function(v)
      local dl = 1/(v.gs - 1)
      local function hd(r) return sqrt(v.N*(v.Ri^2 - r^2)/(2*v.K*(1 + dl))) end
      local out = {"delta = " .. fmt(dl), "hd(r=" .. fmt(v.r) .. ") = " .. fmt(hd(v.r)) .. " m",
                   "hs = delta*hd = " .. fmt(dl*hd(v.r)) .. " m"}
      if v.x >= 0 then
        local r = v.Ri - v.x
        out[#out+1] = "A " .. fmt(v.x) .. " de la costa: hd=" .. fmt(hd(r)) .. " hs=" .. fmt(dl*hd(r))
      end
      return out
    end),
  P("Funcion de pozo W(u)", "W(u) = -0.5772 - ln u + u - u^2/(2*2!) + ...",
    {{"u","u","Argumento u"}},
    function(v) return {"W(u) = " .. fmt(W(v.u)), "1/u = " .. fmt(1/v.u)} end),
  E("Theis (confinado)", "dh = Q/(4piT)*W(u),  u = r^2*S/(4*T*t)",
    {{"dh","dh","Descenso"},{"Q","Q","Caudal"},{"T","T","Transmisibilidad"},{"S","S","Coef. almacenamiento"},{"r","r","Distancia"},{"t","t","Tiempo"}},
    function(v) return v.dh - v.Q/(4*pi*v.T)*W(v.r^2*v.S/(4*v.T*v.t)) end,
    function(v) local u = v.r^2*v.S/(4*v.T*v.t); return {"u = " .. fmt(u), "W(u) = " .. fmt(W(u)),
      (u < 0.01 and "u<0.01: vale Cooper-Jacob" or "u>=0.01: usar Theis")} end),
  E("Cooper-Jacob", "dh = Q/(4piT)*ln(2.24*T*t/(r^2*S))",
    {{"dh","dh","Descenso"},{"Q","Q","Caudal"},{"T","T","Transmisibilidad"},{"S","S","Coef. almacenamiento"},{"r","r","Distancia"},{"t","t","Tiempo"}},
    function(v) return v.dh - v.Q/(4*pi*v.T)*ln(CJ*v.T*v.t/(v.r^2*v.S)) end,
    function(v) local u = v.r^2*v.S/(4*v.T*v.t); return {"u = " .. fmt(u) .. (u < 0.01 and " (ok)" or " (>0.01 revisar)"),
      "R(t) = " .. fmt(sqrt(CJ*v.T*v.t/v.S))} end),
  E("Radio de infl. R(t)", "R(t) = sqrt(2.24*T*t/S)",
    {{"R","R(t)","Radio de influencia"},{"T","T","Transmisibilidad"},{"t","t","Tiempo"},{"S","S","Coef. almacenamiento"}},
    function(v) return v.R - sqrt(CJ*v.T*v.t/v.S) end),
  E("Libre transiente", "H^2 - h^2 = Q/(pi*K)*ln(R(t)/r)",
    {{"H","H","Espesor inicial"},{"h","h","Nivel en r,t"},{"Q","Q","Caudal"},{"K","K","Conductividad"},
     {"S","S","Coef. almac. (Sy)"},{"t","t","Tiempo"},{"r","r","Distancia"}},
    function(v) local R = sqrt(CJ*v.K*v.H*v.t/v.S); return v.H^2 - v.h^2 - v.Q/(pi*v.K)*ln(R/v.r) end,
    function(v) return {"R(t) = " .. fmt(sqrt(CJ*v.K*v.H*v.t/v.S)) .. " (T=K*H)", "dh = " .. fmt(v.H - v.h)} end),
  E("Tiempo de influencia", "ti = r^2*S/(2.24*T)",
    {{"ti","ti","Tiempo en que el pozo afecta a r"},{"r","r","Distancia"},{"S","S","Coef. almac."},{"T","T","Transmisibilidad"}},
    function(v) return v.ti - v.r^2*v.S/(CJ*v.T) end),
  P("Superposicion transiente", "dhT = Sum Qi/(2piT)*ln(R(t-ti)/ri), R=sqrt(2.24T(t-ti)/S). Imagenes: Q con signo",
    {{"X","x pozos","Coordenadas x (lista)", list=true},{"Y","y pozos","Coordenadas y (lista)", list=true},
     {"Q","Q pozos","Caudales (lista, signo)", list=true},{"ts","t inicio","Tiempos de inicio (lista, 0..)", list=true},
     {"T","T","Transmisibilidad"},{"S","S","Coef. almac."},{"t","t","Instante de calculo"},
     {"xp","x punto","x del punto"},{"yp","y punto","y del punto"},{"H","H (libre)","Esp. inicial si libre (0=conf.)", opt="0"}},
    function(v) return superpos(v, true) end),
  P("Bombeo variable", "dh = 1/(4piT) Sum +-Qi*ln(2.24T(t-ti)/(r^2 S)), t>ti",
    {{"Q","Q1,Q2,..","Cambios de caudal (+ encender, - apagar)", list=true},{"ti","t1,t2,..","Instantes de cada cambio", list=true},
     {"T","T","Transmisibilidad"},{"S","S","Coef. almac."},{"r","r","Distancia"},{"t","t","Instante de calculo"}},
    function(v)
      local s, out = 0, {}
      for i = 1, #v.Q do
        if v.t > v.ti[i] then
          local d = v.Q[i]/(4*pi*v.T)*ln(CJ*v.T*(v.t - v.ti[i])/(v.r^2*v.S))
          s = s + d; out[#out+1] = "Termino " .. i .. " = " .. fmt(d)
        end
      end
      out[#out+1] = "DESCENSO = " .. fmt(s)
      return out
    end),
})

-- ===================== CLASE 06 =====================
cat("C06 Pruebas de bombeo", {
  P("Curva tipo Theis (ajuste)", "T = Q*W(u)/(4pi*dh*),  S = 4*T*u*t*/r^2",
    {{"Q","Q","Caudal"},{"W","W(u)","W(u) del punto elegido"},{"iu","1/u","1/u del punto elegido"},
     {"dh","dh*","Descenso del punto"},{"t","t*","Tiempo del punto"},{"r","r","Distancia pozo obs."}},
    function(v) local T = v.Q*v.W/(4*pi*v.dh); return {"T = " .. fmt(T), "S = " .. fmt(4*T*(1/v.iu)*v.t/v.r^2)} end),
  P("Cooper-Jacob con datos", "Regresion dh vs ln(t): T=Q/(4pi*a), S=2.24*T*t0/r^2",
    {{"t","t1,t2,..","Tiempos (lista)", list=true},{"dh","dh1,dh2,..","Descensos (lista)", list=true},
     {"Q","Q","Caudal"},{"r","r","Distancia (radio pozo si es el mismo)"}},
    function(v)
      local xs = {}
      for i = 1, #v.t do xs[i] = ln(v.t[i]) end
      local a, b = regresion(xs, v.dh)
      local T = v.Q/(4*pi*a)
      local t0 = exp(-b/a)
      return {"Pendiente a (por ln t) = " .. fmt(a), "Pendiente por ciclo log10 = " .. fmt(a*ln(10)),
              "Intercepto = " .. fmt(b), "T = " .. fmt(T), "t0 (dh=0) = " .. fmt(t0), "S = " .. fmt(CJ*T*t0/v.r^2)}
    end),
  P("Recuperacion con datos", "Regresion dh vs ln(t/(t-tf)): T = Q/(4pi*a)",
    {{"t","t (desde inicio)","Tiempos totales (lista, >tf)", list=true},{"dh","dh residual","Descensos residuales (lista)", list=true},
     {"tf","tf","Tiempo en que se detuvo el bombeo"},{"Q","Q","Caudal bombeado"}},
    function(v)
      local xs = {}
      for i = 1, #v.t do xs[i] = ln(v.t[i]/(v.t[i] - v.tf)) end
      local a, b = regresion(xs, v.dh)
      return {"Pendiente a = " .. fmt(a), "Intercepto = " .. fmt(b), "T = " .. fmt(v.Q/(4*pi*a))}
    end),
  E("Hantush-Jacob (semiconf.)", "dh = Q/(4piT)*W(u,r/B), B = sqrt(T*m'/K')",
    {{"dh","dh","Descenso"},{"Q","Q","Caudal"},{"T","T","Transmisibilidad"},{"S","S","Coef. almac."},{"r","r","Distancia"},
     {"t","t","Tiempo"},{"mp","m'","Espesor acuitardo"},{"Kp","K'","Conductividad acuitardo"}},
    function(v)
      local B = sqrt(v.T*v.mp/v.Kp)
      return v.dh - v.Q/(4*pi*v.T)*Wh(v.r^2*v.S/(4*v.T*v.t), v.r/B)
    end,
    function(v) local B = sqrt(v.T*v.mp/v.Kp); local u = v.r^2*v.S/(4*v.T*v.t)
      return {"B = " .. fmt(B), "r/B = " .. fmt(v.r/B), "u = " .. fmt(u), "W(u,r/B) = " .. fmt(Wh(u, v.r/B))} end),
  E("Libre: descenso corregido", "dh' = dh - dh^2/(2H)",
    {{"dhp","dh'","Descenso corregido"},{"dh","dh","Descenso medido"},{"H","H","Espesor saturado inicial"}},
    function(v) return v.dhp - (v.dh - v.dh^2/(2*v.H)) end),
  P("Distancia a borde", "r2 = r1*sqrt(t2/t1), d = (r1+r2)/2",
    {{"r1","r1","Dist. pozo bombeo - obs."},{"t1","t1","Tiempo recta original con dh"},{"t2","t2","Tiempo con desviacion dh"}},
    function(v) local r2 = v.r1*sqrt(v.t2/v.t1); return {"r2 (pozo imagen) = " .. fmt(r2), "d al borde = " .. fmt((v.r1 + r2)/2)} end),
})

-- ===================== CLASE 07 =====================
cat("C07 Drenes permanente", {
  E("Dren interceptor", "H^2 - h0^2 = 2*q*x1/K",
    {{"H","H","Nivel lejos del dren"},{"h0","h0","Altura en el canal"},{"q","q","Caudal por metro [m2/s]"},{"x1","x1","Distancia de influencia"},{"K","K","Conductividad"}},
    function(v) return v.H^2 - v.h0^2 - 2*v.q*v.x1/v.K end),
  E("Caudal dren (largo L)", "Qz = (q - q')*L",
    {{"Qz","Qz","Caudal que lleva el canal"},{"q","q","Caudal aguas arriba [m2/s]"},{"qp","q'","Caudal que pasa [m2/s]"},{"L","L","Largo del dren"}},
    function(v) return v.Qz - (v.q - v.qp)*v.L end),
  E("Donnan (abiertos)", "H^2 - h0^2 = f*D^2/(4K)",
    {{"H","H","Altura max. napa"},{"h0","h0","Altura en el dren"},{"f","f","Recarga [m/s]"},{"D","D","Separacion drenes"},{"K","K","Conductividad"}},
    function(v) return v.H^2 - v.h0^2 - v.f*v.D^2/(4*v.K) end),
  E("Hooghoudt (abiertos)", "(H-h0)(H+h0+2d) = f*D^2/(4K)",
    {{"H","H","Altura max."},{"h0","h0","Altura en el dren"},{"d","d (o d')","Dist. al impermeable"},{"f","f","Recarga"},{"D","D","Separacion"},{"K","K","Conductividad"}},
    function(v) return (v.H - v.h0)*(v.H + v.h0 + 2*v.d) - v.f*v.D^2/(4*v.K) end),
  E("Dagan (cerrados)", "(H+d)^2-d^2 = fD^2/(4K)*(1+alfa*4d/D), alfa=-ln(2(cosh(pi r0/d)-1))/pi",
    {{"H","H","Altura max. sobre drenes"},{"d","d","Dist. dren-impermeable"},{"f","f","Recarga"},{"D","D","Separacion"},{"K","K","Conductividad"},{"r0","r0","Radio del dren"}},
    function(v)
      local al = -ln(2*(cosh(pi*v.r0/v.d) - 1))/pi
      if v.r0/v.d >= 0.3 then al = 0 end
      return (v.H + v.d)^2 - v.d^2 - v.f*v.D^2/(4*v.K)*(1 + al*4*v.d/v.D)
    end,
    function(v) local al = -ln(2*(cosh(pi*v.r0/v.d) - 1))/pi; return {"alfa = " .. fmt(al)} end),
  E("Drenes superficiales", "H = D*f/(K*pi)*ln(D/(5.44*r0))  (espesor infinito)",
    {{"H","H","Altura max."},{"D","D","Separacion"},{"f","f","Recarga"},{"K","K","Conductividad"},{"r0","r0","Radio del dren"}},
    function(v) return v.H - v.D*v.f/(v.K*pi)*ln(v.D/(5.44*v.r0)) end),
  E("Kirkham", "H = D*f/K * F(2r0/D, d/D)  (F por serie)",
    {{"H","H","Altura max."},{"D","D","Separacion"},{"f","f","Recarga"},{"K","K","Conductividad"},{"r0","r0","Radio dren"},{"d","d","Dist. al impermeable"}},
    function(v) return v.H - v.D*v.f/v.K*Fk(2*v.r0/v.D, v.d/v.D) end,
    function(v) return {"F = " .. fmt(Fk(2*v.r0/v.D, v.d/v.D))} end),
  P("Tabla Kirkham F(x,y)", "F(x,y) = 1/pi[ln(2/(pi x)) + Sum 1/m(cos m pi x - cos m pi)(coth 2m pi y - 1)]",
    {{"x","2r0/D","x = 2r0/D"},{"y","d/D","y = d/D"}},
    function(v) return {"F = " .. fmt(Fk(v.x, v.y))} end),
})

-- ===================== CLASE 08 =====================
local function dprime(d, D, r0)
  local chi = pi*r0
  if d/D < 0.25 then return d/(1 + 8/pi*d/D*ln(d/chi)) end
  return d*pi/(2*(ln(D/chi) + 0.18))
end
local function gloverD(K, d, t, S, y0, y) return pi*sqrt(K*d*t/(S*ln(1.16*y0/y))) end
cat("C08 Drenes impermanente", {
  E("Glover-Dumm", "y = 1.16*y0*exp(-alfa*t), alfa = pi^2*K*d/(S*D^2)",
    {{"y","y","Altura final (sobre drenes)"},{"y0","y0","Altura inicial"},{"t","t","Tiempo"},{"K","K","Conductividad"},
     {"d","d","Prof. equivalente"},{"S","S","Porosidad drenable"},{"D","D","Separacion drenes"}},
    function(v) return v.y - 1.16*v.y0*exp(-pi^2*v.K*v.d/(v.S*v.D^2)*v.t) end,
    function(v) local a = pi^2*v.K*v.d/(v.S*v.D^2); return {"alfa = " .. fmt(a), "alfa*t = " .. fmt(a*v.t) .. (a*v.t >= 0.2 and " (ok >=0.2)" or " (<0.2!)")} end),
  P("Hooghoudt d' (Glover)", "d/D<0.25: d'=d/(1+8/pi d/D ln(d/chi)) | si no: d'=d*pi/(2(ln(D/chi)+0.18)), chi=pi*r0",
    {{"d","d","Dist. al impermeable"},{"D","D","Separacion"},{"r0","r0","Radio del dren"}},
    function(v) return {"chi = " .. fmt(pi*v.r0), "d/D = " .. fmt(v.d/v.D), "d' = " .. fmt(dprime(v.d, v.D, v.r0))} end),
  P("Iterar D (Glover+Hooghoudt)", "1) D con d; 2) d' con D; 3) D con d'; repetir",
    {{"K","K","Conductividad"},{"d","d","Dist. al impermeable"},{"S","S","Porosidad drenable"},{"t","t","Tiempo"},
     {"y0","y0","Altura inicial"},{"y","y","Altura final"},{"r0","r0","Radio dren"}},
    function(v)
      local out, D, dp = {}, gloverD(v.K, v.d, v.t, v.S, v.y0, v.y), v.d
      out[1] = "Iter 0: D = " .. fmt(D)
      for i = 1, 30 do
        dp = dprime(v.d, D, v.r0)
        local Dn = gloverD(v.K, dp, v.t, v.S, v.y0, v.y)
        out[#out+1] = "Iter " .. i .. ": d'=" .. fmt(dp) .. " D=" .. fmt(Dn)
        if abs(Dn - D) < 1e-6*D then D = Dn; break end
        D = Dn
      end
      out[#out+1] = "D FINAL = " .. fmt(D)
      return out
    end),
  P("Equilibrio dinamico", "y0,i = yf,i-1 + Ri/S ; yf,i = 1.16 y0,i exp(-alfa*ti) ; Is=int(y>ya)dt",
    {{"R","R mensual","Recargas por mes [mm] (lista)", list=true},{"dias","dias/mes","Dias de cada mes", opt="30"},
     {"K","K [m/dia]","Conductividad [m/dia]"},{"d","d","Prof. equivalente [m]"},{"S","S","Porosidad drenable"},
     {"D","D","Separacion [m]"},{"yi","y inicial","Altura al inicio [m]", opt="0"},{"ya","y raices","Altura limite (H0-Za) [m]", opt="1E9"}},
    function(v)
      local al = pi^2*v.K*v.d/(v.S*v.D^2)
      local out, y, Is = {"alfa = " .. fmt(al) .. " 1/dia"}, v.yi, 0
      local y0first
      for i = 1, #v.R do
        local y0 = y + v.R[i]/1000/v.S
        if i == 1 then y0first = y0 end
        local yf = 1.16*y0*exp(-al*v.dias)
        local n = 60
        for k = 0, n - 1 do
          local yy = 1.16*y0*exp(-al*(k + 0.5)*v.dias/n)
          if yy > v.ya then Is = Is + (yy - v.ya)*1000*v.dias/n end
        end
        out[#out+1] = "Mes " .. i .. ": y0=" .. fmt(y0) .. " yf=" .. fmt(yf)
        y = yf
      end
      out[#out+1] = "yf final = " .. fmt(y) .. ((y <= y0first) and " <= y0 (equilibrio OK)" or " > y0 (NO)")
      out[#out+1] = "Is = " .. fmt(Is) .. " mm*dia" .. ((Is <= 200) and " (OK <=200)" or " (>200 NO)")
      return out
    end),
})

-- ===================== CLASE 09 =====================
local function canal(v)
  local out = {}
  local tp, b, k, D = v.tp, v.b, v.k, v.D
  local function fn(y) local A, Pm = seccion(tp, b, k, D, y); return v.Q*v.n/sqrt(v.S) - A^(5/3)/Pm^(2/3) end
  local function fc(y) local A, _, B = seccion(tp, b, k, D, y); return v.Q^2/G - A^3/B end
  local lim = (tp == 4) and D*0.9999 or 1e4
  local yn = roots(fn, 1e-6, lim, 2000, true)
  local yc = roots(fc, 1e-6, lim, 2000, true)
  local ync, ycc = yn[1], yc[1]
  if #yn == 0 then out[#out+1] = "Sin altura normal (seccion insuficiente)" end
  for i, y in ipairs(yn) do
    local A, Pm, B = seccion(tp, b, k, D, y)
    out[#out+1] = "yn" .. (#yn > 1 and i or "") .. " = " .. fmt(y) .. "  v=" .. fmt(v.Q/A) .. " Fr=" .. fmt(v.Q/A/sqrt(G*A/B))
    out[#out+1] = "   A=" .. fmt(A) .. " P=" .. fmt(Pm) .. " Rh=" .. fmt(A/Pm) .. " B=" .. fmt(B)
  end
  for _, y in ipairs(yc) do out[#out+1] = "yc = " .. fmt(y) end
  if ync and ycc then out[#out+1] = (ync > ycc) and "Regimen SUBCRITICO (rio)" or "Regimen SUPERCRITICO (torrente)" end
  if ync then out[#out+1] = "Revancha 0.15-0.2 h: " .. fmt(0.15*ync) .. " a " .. fmt(0.2*ync) end
  return out
end
cat("C09 Drenaje longitudinal", {
  E("Formula racional", "Q = C*i*A/3.6   (i mm/h, A km2, Q m3/s)",
    {{"Q","Q","Caudal [m3/s]"},{"C","C","Coef. escorrentia"},{"i","i","Intensidad [mm/h]"},{"A","A","Area [km2]"}},
    function(v) return v.Q - v.C*v.i*v.A/3.6 end),
  P("Tiempo de concentracion", "Complete lo que tenga; calcula las formulas posibles (min)",
    {{"L","L [km]","Longitud cauce [km]"},{"S","S [%]","Pendiente [%]", opt="-1"},{"H","H [m]","Desnivel total [m]", opt="-1"},
     {"A","A [km2]","Area [km2]", opt="-1"},{"Hm","Hm [m]","Cota media - salida [m]", opt="-1"},{"CN","CN","Curva numero", opt="-1"}},
    function(v)
      local o = {}
      if v.S > 0 then o[#o+1] = "Normas Esp.: " .. fmt(18*v.L^0.76/v.S^0.19) end
      if v.H > 0 then o[#o+1] = "California: " .. fmt(57*(v.L^3/v.H)^0.385) end
      if v.A > 0 and v.Hm > 0 then o[#o+1] = "Giandotti: " .. fmt(60*(4*sqrt(v.A) + 1.5*v.L)/(0.8*sqrt(v.Hm))) end
      if v.CN > 0 and v.S > 0 then o[#o+1] = "SCS: " .. fmt(3.42*v.L^0.8*(1000/v.CN - 9)^0.7/v.S^0.5) end
      o[#o+1] = "Verificar tc >= 10 min"
      return o
    end),
  P("Tc areas planas", "FAA, Izzard, Morgali-Linsley (min)",
    {{"Ls","Ls [m]","Long. escurrimiento superficial [m]"},{"S","S","Pendiente (m/m)"},{"C","C","Coef. escorrentia", opt="-1"},
     {"i","i [mm/h]","Intensidad", opt="-1"},{"n","n","Manning superficial", opt="-1"}},
    function(v)
      local o = {}
      if v.C > 0 then o[#o+1] = "FAA: " .. fmt(3.26*(1.1 - v.C)*v.Ls^0.5/(100*v.S)^0.33) end
      if v.C > 0 and v.i > 0 then o[#o+1] = "Izzard: " .. fmt(525.28*(0.0000276*v.i + v.C)*v.Ls^0.33/(v.i^0.667*v.S^0.333)) end
      if v.n > 0 and v.i > 0 then o[#o+1] = "Morgali: " .. fmt(7*v.Ls^0.6*v.n^0.6/(v.i^0.4*v.S^0.3)) end
      return o
    end),
  E("Grunsky", "it = i24*sqrt(24/t)  (t en horas)",
    {{"it","it","Intensidad para t"},{"i24","i24","Intensidad 24 h"},{"t","t [h]","Duracion"}},
    function(v) return v.it - v.i24*sqrt(24/v.t) end),
  E("Coef. duracion/frecuencia", "P = K*CDt*CFT*PD10  (K=1.1)",
    {{"P","P","Lluvia T, t [mm]"},{"CD","CDt","Coef. duracion"},{"CF","CFT","Coef. frecuencia"},{"PD","PD10","Lluvia diaria T=10"},{"K","K","Correccion", opt="1.1"}},
    function(v) return v.P - v.K*v.CD*v.CF*v.PD end,
    function(v) return {"Intensidad = P/t (dividir por la duracion en h)"} end),
  E("Bell (t<1 h)", "P = (0.54 t^0.25 - 0.50)(0.21 lnT + 0.52) P1^10  (t min)",
    {{"P","P","Lluvia [mm]"},{"t","t [min]","Duracion"},{"T","T [anos]","Periodo retorno"},{"P1","P1^10","Lluvia 1 h, T=10 [mm]"}},
    function(v) return v.P - (0.54*v.t^0.25 - 0.5)*(0.21*ln(v.T) + 0.52)*v.P1 end,
    function(v) return {"i = P*60/t = " .. fmt(v.P*60/v.t) .. " mm/h"} end),
  P("Canal: Manning yn, yc", "Q*n/sqrt(S) = A^(5/3)/P^(2/3) ; Q^2/g = A^3/B",
    {{"tp","Tipo 1-4","1 rect(b) 2 trapecio(b,k) 3 triang(k) 4 circular(D)"},{"Q","Q","Caudal [m3/s]"},{"n","n","Manning"},
     {"S","S","Pendiente [m/m]"},{"b","b","Ancho basal", opt="0"},{"k","k talud","Talud H:V", opt="0"},{"D","D","Diametro", opt="0"}},
    canal),
  E("Riesgo de falla", "R = 1 - (1 - 1/T)^n",
    {{"R","R","Riesgo (0-1)"},{"T","T","Periodo de retorno"},{"n","n","Vida util"}},
    function(v) return v.R - (1 - (1 - 1/v.T)^v.n) end),
})

-- ===================== CLASE 10 =====================
local function areaFull(tp, b, D) if tp == 4 then return pi*D*D/4, D/4 end; return b*D, b*D/(2*b + 2*D) end
cat("C10 Alcantarillas", {
  E("Rect. control entrada", "D = 1.5/g^(1/3)*(Q/b)^(2/3)  (D=H=Hc)",
    {{"D","D=H","Altura cajon"},{"Q","Q","Caudal T1"},{"b","b","Ancho cajon"}},
    function(v) return v.D - 1.5/G^(1/3)*(v.Q/v.b)^(2/3) end,
    function(v) local hc = (v.Q^2/(G*v.b^2))^(1/3); return {"hc = " .. fmt(hc), "v = " .. fmt(v.Q/(v.b*hc)) .. " (<5 m/s)"} end),
  E("Rect. pendiente critica", "Sc = 2/3*g*H*n^2*((b+4/3H)/(2/3*H*b))^(4/3)",
    {{"Sc","Sc","Pendiente critica"},{"H","H","Altura (=D)"},{"n","n","Manning"},{"b","b","Ancho"}},
    function(v) return v.Sc - 2/3*G*v.H*v.n^2*((v.b + 4/3*v.H)/(2/3*v.H*v.b))^(4/3) end,
    function(v) return {"Si SD >= Sc -> torrente; si no SD = Sc"} end),
  E("Circular control entrada", "Q = 1.425*D^(5/2)  (hc=0.7D)",
    {{"Q","Q","Caudal [m3/s]"},{"D","D","Diametro [m]"}},
    function(v) return v.Q - 1.425*v.D^2.5 end,
    function(v) return {"Sc = 31.14 n^2/D^(1/3): use item siguiente"} end),
  E("Circular pendiente critica", "Sc = 31.14*n^2/D^(1/3)",
    {{"Sc","Sc","Pendiente critica"},{"n","n","Manning"},{"D","D","Diametro"}},
    function(v) return v.Sc - 31.14*v.n^2/v.D^(1/3) end),
  E("Sn circular llena", "Sn = 10.3*Q^2*n^2/D^(16/3)",
    {{"Sn","Sn","Pendiente normal a seccion llena"},{"Q","Q2D","Caudal T2"},{"n","n","Manning"},{"D","D","Diametro"}},
    function(v) return v.Sn - 10.3*v.Q^2*v.n^2/v.D^(16/3) end),
  E("Sn rectangular llena", "Q*n/sqrt(Sn) = A*Rh^(2/3), A=bD, Rh=bD/(2b+2D)",
    {{"Sn","Sn","Pendiente normal"},{"Q","Q2D","Caudal T2"},{"n","n","Manning"},{"b","b","Ancho"},{"D","D","Altura"}},
    function(v) local A, R = v.b*v.D, v.b*v.D/(2*v.b + 2*v.D); return v.Q*v.n/sqrt(v.Sn) - A*R^(2/3) end,
    function(v) return {"Caso 1: SD>Sn | Caso 2: SD=Sn | Caso 3: SD<Sn"} end),
  E("Caso 1: SD>Sn (compuerta)", "H' = cc*D + 1/(2g)*(Q/(b*cc*D))^2",
    {{"Hp","H'","Carga aguas arriba"},{"cc","cc","Coef. contraccion", opt="0.611"},{"D","D","Altura"},{"Q","Q2D","Caudal"},{"b","b","Ancho"}},
    function(v) return v.Hp - v.cc*v.D - (v.Q/(v.b*v.cc*v.D))^2/(2*G) end),
  E("Caso 2: SD=Sn", "H' = D + (ke+1)/(2g)*Q^2/(D*b)^2",
    {{"Hp","H'","Carga aguas arriba"},{"D","D","Altura"},{"ke","ke","Perdida entrada (0.2/0.5/0.9)"},{"Q","Q2D","Caudal"},{"b","b","Ancho"}},
    function(v) return v.Hp - v.D - (v.ke + 1)*v.Q^2/((v.D*v.b)^2*2*G) end),
  P("Caso 3 / control salida", "H' = D + v^2/2g + ke v^2/2g + ks(v-v3)^2/2g + J*L - SD*L, J=v^2 n^2/Rh^(4/3)",
    {{"tp","Tipo 1/4","1 rectangular, 4 circular"},{"Q","Q2D","Caudal"},{"b","b","Ancho (rect)", opt="0"},{"D","D","Altura/diametro"},
     {"n","n","Manning"},{"L","L","Largo alcantarilla"},{"SD","SD","Pendiente de diseno"},{"ke","ke","Coef. entrada"},
     {"ks","ks","Coef. salida (0=caso 3; 0.4)", opt="0"},{"v3","v3","Velocidad aguas abajo", opt="0"}},
    function(v)
      local A, R = areaFull(v.tp, v.b, v.D)
      local vel = v.Q/A
      local J = vel^2*v.n^2/R^(4/3)
      local Hp = v.D + vel^2/(2*G) + v.ke*vel^2/(2*G) + v.ks*(vel - v.v3)^2/(2*G) + J*v.L - v.SD*v.L
      return {"A = " .. fmt(A) .. "  Rh = " .. fmt(R), "v = " .. fmt(vel), "J = " .. fmt(J), "H' = " .. fmt(Hp)}
    end),
})

-- ===================== CLASE 11 =====================
cat("C11 Puentes", {
  P("L minimo (momenta)", "M1 = Q^2/(g b h1)+b h1^2/2 = M2c(L), hc2=(Q^2/(g b L))^(1/3)",
    {{"Q","Q","Caudal"},{"b","b","Ancho del cauce"},{"h1","hn1","Altura normal en el cauce"}},
    function(v)
      local M1 = v.Q^2/(G*v.b*v.h1) + v.b*v.h1^2/2
      local function M2c(L) local hc = (v.Q^2/(G*v.b*L))^(1/3); return v.Q^2/(G*L*hc) + v.b*hc^2/2 end
      local r = roots(function(L) return M2c(L) - M1 end, 1e-4*v.b, v.b, 2000, true)
      local L = largest(r)
      local o = {"M1 = " .. fmt(M1)}
      if L then o[#o+1] = "L minimo = " .. fmt(L); o[#o+1] = "hc2 = " .. fmt((v.Q^2/(G*v.b*L))^(1/3))
      else o[#o+1] = "No hay solucion con L<=b" end
      return o
    end),
  P("Alturas h2 y h3 (energia)", "h2+v2^2/2g = hn1+v1^2/2g+Ls ; h3+v3^2/2g = h2+(1+ke)v2^2/2g",
    {{"Q","Q","Caudal"},{"b","b","Ancho cauce"},{"L","L","Luz del puente"},{"h1","hn1","Altura normal aguas abajo"},
     {"ke","ke","Coef. entrada (0.1-0.3)", opt="0.2"}},
    function(v)
      local v1 = v.Q/(v.h1*v.b)
      local function f2(h2)
        local v2 = v.Q/(h2*v.L)
        local Ls = (v2 - v1)^2/(2*G) - (v.h1 - h2)^2/(2*v.h1)
        return h2 + v2^2/(2*G) - v.h1 - v1^2/(2*G) - Ls
      end
      local h2 = largest(roots(f2, 1e-4, 100*v.h1, 3000, true))
      if not h2 then return {"Sin solucion para h2 (resalto?)"} end
      local v2 = v.Q/(h2*v.L)
      local E2 = h2 + (1 + v.ke)*v2^2/(2*G)
      local h3 = largest(roots(function(h) return h + (v.Q/(h*v.b))^2/(2*G) - E2 end, 1e-4, 100*v.h1, 3000, true))
      local o = {"v1 = " .. fmt(v1), "h2 = " .. fmt(h2), "v2 = " .. fmt(v2) .. " (verificar < vmax)"}
      if h3 then o[#o+1] = "h3 = " .. fmt(h3) .. " (verificar < hmax)"; o[#o+1] = "Remanso h3-hn1 = " .. fmt(h3 - v.h1) end
      return o
    end),
  E("Momentum 3-2 (X2,X3)", "1/X2 - 1/(n X3) = 0.5(X3^2 - X2^2), X=h/hc2, n=b3/b2",
    {{"X3","X3","h3/hc2"},{"X2","X2","h2/hc2 (1-6)"},{"n","n","b3/b2 (1-4.5)"}},
    function(v) return 1/v.X2 - 1/(v.n*v.X3) - 0.5*(v.X3^2 - v.X2^2) end),
  P("Yarnell (cepas)", "dh/h1 = K Fr1^2 (K + 5Fr1^2 - 0.6)(sig + 15 sig^4), sig = nD/b",
    {{"K","K","Forma cepa (0.9-2.5)"},{"Q","Q","Caudal"},{"b","b","Ancho del puente"},{"h1","h1","Altura aguas abajo"},
     {"n","n","N. de cepas"},{"D","D","Ancho de cada cepa"}},
    function(v)
      local Fr2 = (v.Q/(v.b*v.h1))^2/(G*v.h1)
      local sig = v.n*v.D/v.b
      local r = v.K*Fr2*(v.K + 5*Fr2 - 0.6)*(sig + 15*sig^4)
      return {"Fr1^2 = " .. fmt(Fr2), "sigma = " .. fmt(sig) .. "  be = " .. fmt(v.b - v.n*v.D),
              "dh/h1 = " .. fmt(r), "dh = " .. fmt(r*v.h1), "h3 = " .. fmt(v.h1*(1 + r))}
    end),
})

-- ===================== CLASE 12 =====================
cat("C12 Drenaje urbano", {
  E("Cuneta simple", "Q = sqrt(S)/n * b^2 i/2 * (b i/(2(i+1)))^(2/3)",
    {{"Q","Q","Caudal [m3/s]"},{"S","S","Pendiente long. [m/m]"},{"n","n","Manning"},{"b","b","Ancho inundado [m]"},{"i","i","Bombeo [m/m]"}},
    function(v) return v.Q - sqrt(v.S)/v.n*v.b^2*v.i/2*(v.b*v.i/(2*(v.i + 1)))^(2/3) end,
    function(v) return {"y = b*i = " .. fmt(v.b*v.i) .. " (<= 0.15 m)", "v = " .. fmt(v.Q/(v.b^2*v.i/2))} end),
  E("Cuneta compuesta", "Q = 0.315 sqrt(S)/n b^(8/3)[i2+(i1-i2)(w/b)^2]^(5/3)",
    {{"Q","Q","Caudal"},{"S","S","Pendiente long."},{"n","n","Manning"},{"b","b","Ancho inundado (>w)"},
     {"i1","i1","Pendiente zarpa"},{"i2","i2","Bombeo calle"},{"w","w","Ancho zarpa (~0.6)"}},
    function(v) return v.Q - 0.315*sqrt(v.S)/v.n*v.b^(8/3)*(v.i2 + (v.i1 - v.i2)*(v.w/v.b)^2)^(5/3) end),
  P("Sumidero de fondo Qm", "Vertedero: Qm=1.66(Le+2we)h^1.5 | Orificio: Qm=2.66 Ae h^0.5",
    {{"w","w","Ancho reja"},{"L","L","Largo reja"},{"e","e","Espesor barras"},{"nL","nL","N. barras longit."},
     {"nT","nT","N. barras transv."},{"h","h","Altura de agua"}},
    function(v)
      local we, Le = v.w - v.e*v.nL, v.L - v.e*v.nT
      local Ae = we*Le
      local lim = 1.6*Ae/(Le + 2*we)
      local Qm = (v.h < lim) and 1.66*(Le + 2*we)*v.h^1.5 or 2.66*Ae*v.h^0.5
      return {"we = " .. fmt(we) .. "  Le = " .. fmt(Le) .. "  Ae = " .. fmt(Ae), "h limite = " .. fmt(lim),
              (v.h < lim) and "Funciona como VERTEDERO" or "Funciona como ORIFICIO", "Qm = " .. fmt(Qm)}
    end),
  P("Eficiencia sumid. fondo", "etaH = Rf eta0 + Rs(1-eta0); Qs = min(etaH*Q, Qm)",
    {{"Q","Q","Caudal en la cuneta"},{"w","w","Ancho reja"},{"bs","bs","Ancho inundado"},{"v","v","Velocidad"},
     {"v0","v0","Veloc. de salpicadura"},{"i","i","Bombeo"},{"L","L","Largo reja"},{"Qm","Qm","Capacidad max.", opt="1E9"},
     {"F","F","Factor seguridad", opt="1"}},
    function(v)
      local function cl(x) return math.max(0, math.min(1, x)) end
      local e0 = cl(1 - (1 - v.w/v.bs)^2.67)
      local Rs = cl(1/(1 + 0.0828*v.v^1.8/(v.i*v.L^2.3)))
      local Rf = cl(1 - 0.295*(v.v - v.v0))
      local eH = cl(Rf*e0 + Rs*(1 - e0))*v.F
      return {"eta0 = " .. fmt(e0), "Rs = " .. fmt(Rs), "Rf = " .. fmt(Rf), "etaH (con F) = " .. fmt(eH),
              "Qs = " .. fmt(math.min(eH*v.Q, v.Qm))}
    end),
  P("Sumidero lateral", "Qm=1.27 L h^1.5 (h<a) | 2.66 L a h^0.5 ; etaL=1-(1-L/LT)^1.8",
    {{"Q","Q","Caudal en la cuneta"},{"L","L","Largo ventana"},{"a","a","Altura ventana"},{"h","h","Altura de agua"},
     {"S","S","Pendiente long."},{"n","n","Manning"},{"i","i (o i2)","Bombeo"},{"F","F","Factor seguridad", opt="1"}},
    function(v)
      local Qm = (v.h < v.a) and 1.27*v.L*v.h^1.5 or 2.66*v.L*v.a*v.h^0.5
      local LT = 0.81*v.Q^0.42*v.S^0.3*(v.n*v.i)^-0.6
      local eL = (v.h > v.a) and 1 or (1 - (1 - math.min(1, v.L/LT))^1.8)
      eL = eL*v.F
      return {"Qm = " .. fmt(Qm), "LT = " .. fmt(LT), "etaL = " .. fmt(eL), "Qs = " .. fmt(math.min(eL*v.Q, Qm)),
              "Mixto: Qs = min((etaH+etaL)Q, Qm)"}
    end),
  E("Distancia entre sumideros", "L = 3600*F*eta*Qmax/(C*i*bc) , bc = 0.5e+wc o e+wc",
    {{"L","L","Separacion [m]"},{"F","F","Factor seguridad"},{"eta","eta","Eficiencia"},{"Q","Qmax","Caudal max. [m3/s]"},
     {"C","C","Coef. escorrentia"},{"i","i [mm/h]","Intensidad"},{"bc","bc","Ancho aportante [m]"}},
    function(v) return v.L - 3600*v.F*v.eta*v.Q/(v.C*v.i/1000*v.bc) end,
    function(v) return {"Nota: i convertida a m/h (/1000)"} end),
  P("Colector circular", "Manning con h = fraccion*D ; Vc = 0.397 D^(2/3) sqrt(S)/n >= 0.6",
    {{"Q","Q","Caudal de diseno"},{"n","n","Manning"},{"S","S","Pendiente"},{"r","h/D","Llenado de diseno", opt="0.8"},
     {"D","D (opc.)","Diametro existente (0 = calcular)", opt="0"}},
    function(v)
      local o = {}
      local function Qcap(D) local A, Pm = seccion(4, 0, 0, D, v.r*D); return sqrt(v.S)/v.n*A*(A/Pm)^(2/3) end
      local D = v.D
      if D <= 0 then
        D = largest(roots(function(d) return Qcap(d) - v.Q end, 1e-3, 20, 2000, true))
        o[#o+1] = "D requerido = " .. fmt(D) .. " (min 0.3 m)"
      end
      local A, Pm = seccion(4, 0, 0, D, v.r*D)
      o[#o+1] = "Q capacidad (h=" .. fmt(v.r) .. "D) = " .. fmt(Qcap(D))
      o[#o+1] = "v = " .. fmt(v.Q/A)
      o[#o+1] = "Vc autolavado = " .. fmt(0.397*D^(2/3)*sqrt(v.S)/v.n) .. " (>=0.6)"
      local yn = roots(function(y) local a, p = seccion(4, 0, 0, D, y); return v.Q*v.n/sqrt(v.S) - a^(5/3)/p^(2/3) end, 1e-5, D*0.9999, 2000, true)
      if yn[1] then o[#o+1] = "yn real = " .. fmt(yn[1]) .. " (h/D=" .. fmt(yn[1]/D) .. ")" end
      return o
    end),
  P("Cuencas en serie", "C medio, tc de interconexion: tc = max(tc_prev + L/v, tc_i)",
    {{"C","C1,C2,..","Coef. escorrentia (lista)", list=true},{"A","A1,A2,..","Areas [km2] (lista)", list=true},
     {"tc","tc1,tc2,..","Tiempos de conc. propios (lista)", list=true},{"L","L1,L2,..","Largos colectores (lista)", list=true},
     {"v","v1,v2,..","Velocidades en colectores (lista)", list=true}},
    function(v)
      local o, sC, sA, tcI = {}, 0, 0, 0
      for k = 1, #v.C do
        sC = sC + v.C[k]; sA = sA + v.A[k]
        if k == 1 then tcI = v.tc[1] else tcI = math.max(tcI + v.L[k-1]/v.v[k-1], v.tc[k]) end
        o[#o+1] = "Colector " .. k .. ": Cmed=" .. fmt(sC/k) .. " Aacum=" .. fmt(sA) .. " tcI=" .. fmt(tcI)
      end
      o[#o+1] = "QI = Cmed*i(tcI)*Aacum/3.6 ; QD = max(QD ant, Qi, QI)"
      o[#o+1] = "(tt = L/v en las mismas unidades que tc)"
      return o
    end),
  P("Hidrograma SCS", "tp = tLL/2 + 0.6 tc ; tr = 1.67 tp ; tB = 2.67 tp ; Qp = 2 Pef A/tB",
    {{"tLL","tLL [h]","Duracion lluvia efectiva"},{"tc","tc [h]","Tiempo de concentracion"},{"P","Pef [mm]","Lluvia efectiva"},{"A","A [km2]","Area"}},
    function(v)
      local tp = v.tLL/2 + 0.6*v.tc
      local tB = 2.67*tp
      return {"tp = " .. fmt(tp) .. " h", "tr = " .. fmt(1.67*tp) .. " h", "tB = " .. fmt(tB) .. " h",
              "Qp = " .. fmt(2*v.P*v.A/tB/3.6) .. " m3/s"}
    end),
})

-- =====================================================================
--  MOTOR DE GRAFICOS
--  Un grafico = { title, xl, yl, xlog, ylog, yinv, series = { {name, pts,
--   style="line"|"pts"|"both", col} }, vl = { {x, label} }, hl = {...},
--   info = texto o function(P), move = function(P,dx,dy) (interactivo),
--   draw = function(gc, box, P) (esquema dibujado) }
-- =====================================================================
local COL = {{0,90,190},{210,30,30},{0,140,60},{230,120,0},{120,40,160},{0,150,160},{90,90,90}}
local function log10(x) return ln(x)/ln(10) end
local function linspace(a, b, n) local t = {}; for i = 0, n - 1 do t[#t+1] = a + (b - a)*i/(n - 1) end; return t end
local function logspace(a, b, n) local t = {}; for i = 0, n - 1 do t[#t+1] = a*(b/a)^(i/(n - 1)) end; return t end
local function S(name, xs, ys, style, col) local p = {}; for i = 1, #xs do p[i] = {xs[i], ys[i]} end
  return {name = name, pts = p, style = style or "line", col = col} end
local function fn2(x) -- numero corto para ejes
  if x == 0 then return "0" end
  local a = abs(x)
  if a >= 1e4 or a < 1e-2 then return string.format("%.0e", x) end
  return string.format("%.3g", x)
end

local function niceStep(r)
  local s = 10^math.floor(log10(r/4))
  local q = r/4/s
  if q > 5 then s = s*10 elseif q > 2 then s = s*5 elseif q > 1 then s = s*2 end
  return s
end

local function tx(P, x) if P.xlog then if x > 0 then return log10(x) end return nil end return x end
local function ty(P, y) if P.ylog then if y > 0 then return log10(y) end return nil end return y end

local function bounds(P)
  local x0, x1, y0, y1
  local function add(x, y)
    if x and ok(x) then x0 = x0 and math.min(x0, x) or x; x1 = x1 and math.max(x1, x) or x end
    if y and ok(y) then y0 = y0 and math.min(y0, y) or y; y1 = y1 and math.max(y1, y) or y end
  end
  for _, s in ipairs(P.series or {}) do
    for _, p in ipairs(s.pts) do
      local X, Y = tx(P, p[1]), ty(P, p[2])
      if X and Y and ok(X) and ok(Y) then add(X, Y) end
    end
  end
  for _, v in ipairs(P.vl or {}) do add(tx(P, v[1]), nil) end
  for _, h in ipairs(P.hl or {}) do add(nil, ty(P, h[1])) end
  if P.xr then x0, x1 = tx(P, P.xr[1]), tx(P, P.xr[2]) end
  if P.yr then y0, y1 = ty(P, P.yr[1]), ty(P, P.yr[2]) end
  x0, x1, y0, y1 = x0 or 0, x1 or 1, y0 or 0, y1 or 1
  if x1 - x0 < 1e-12 then x0, x1 = x0 - 1, x1 + 1 end
  if y1 - y0 < 1e-12 then y0, y1 = y0 - 1, y1 + 1 end
  local px, py = (x1 - x0)*0.04, (y1 - y0)*0.06
  if not P.xr then x0, x1 = x0 - px, x1 + px end
  if not P.yr then y0, y1 = y0 - py, y1 + py end
  if P.y0 ~= nil and not P.ylog and y0 > 0 and y0 < 0.3*(y1 - y0) then y0 = 0 end
  return x0, x1, y0, y1
end

local function drawPlot(gc, P, bx, by, bw, bh)
  if P.draw then P.draw(gc, bx, by, bw, bh, P); return end
  local L, R, T, B = bx + 34, bx + bw - 6, by + 4, by + bh - 24
  local x0, x1, y0, y1 = bounds(P)
  local function X(x) return L + (x - x0)/(x1 - x0)*(R - L) end
  local function Y(y)
    local f = (y - y0)/(y1 - y0)
    if P.yinv then return T + f*(B - T) end
    return B - f*(B - T)
  end
  gc:setFont("sansserif", "r", 7)
  -- rejilla y marcas
  local function ticks(a, b, islog)
    local t = {}
    if islog then
      for k = math.ceil(a), math.floor(b) do t[#t+1] = {k, "1E" .. k} end
      if #t < 2 then for k = math.floor(a), math.ceil(b) do
        for _, m in ipairs({2, 5}) do local v = k + log10(m); if v > a and v < b then t[#t+1] = {v, fn2(10^v)} end end end end
    else
      local s = niceStep(b - a)
      local v = math.ceil(a/s)*s
      while v <= b + 1e-9*s do t[#t+1] = {v, fn2(abs(v) < 1e-9*s and 0 or v)}; v = v + s end
    end
    return t
  end
  gc:setColorRGB(225, 225, 225)
  for _, t in ipairs(ticks(x0, x1, P.xlog)) do gc:drawLine(X(t[1]), T, X(t[1]), B) end
  for _, t in ipairs(ticks(y0, y1, P.ylog)) do gc:drawLine(L, Y(t[1]), R, Y(t[1])) end
  gc:setColorRGB(0, 0, 0)
  gc:drawRect(L, T, R - L, B - T)
  for _, t in ipairs(ticks(x0, x1, P.xlog)) do
    local w = gc:getStringWidth(t[2]); gc:drawString(t[2], X(t[1]) - w/2, B + 1, "top") end
  for _, t in ipairs(ticks(y0, y1, P.ylog)) do
    local w = gc:getStringWidth(t[2]); gc:drawString(t[2], L - w - 2, Y(t[1]) - 5, "top") end
  if P.xl then local w = gc:getStringWidth(P.xl); gc:drawString(P.xl, R - w, B + 11, "top") end
  if P.yl then gc:drawString(P.yl, bx + 1, B + 11, "top") end
  -- lineas de referencia
  gc:setPen("thin", "dashed")
  for _, v in ipairs(P.vl or {}) do
    local xv = tx(P, v[1])
    if xv and xv >= x0 and xv <= x1 then
      gc:setColorRGB(150, 0, 150); gc:drawLine(X(xv), T, X(xv), B)
      if v[2] then gc:drawString(v[2], X(xv) + 2, T + 1, "top") end
    end
  end
  for _, h in ipairs(P.hl or {}) do
    local yv = ty(P, h[1])
    if yv and yv >= math.min(y0, y1) and yv <= math.max(y0, y1) then
      gc:setColorRGB(150, 0, 150); gc:drawLine(L, Y(yv), R, Y(yv))
      if h[2] then gc:drawString(h[2], L + 2, Y(yv) - 10, "top") end
    end
  end
  gc:setPen("thin", "smooth")
  -- series
  local function inside(px, py) return px >= L - 1 and px <= R + 1 and py >= T - 1 and py <= B + 1 end
  for k, s in ipairs(P.series or {}) do
    local c = s.col or COL[(k - 1) % #COL + 1]
    local st = s.style or "line"
    gc:setColorRGB(c[1], c[2], c[3])
    local prev
    for _, p in ipairs(s.pts) do
      local xv, yv = tx(P, p[1]), ty(P, p[2])
      local cur
      if xv and yv and ok(xv) and ok(yv) then cur = {X(xv), Y(yv)} end
      if cur and st ~= "pts" and prev and inside(prev[1], prev[2]) and inside(cur[1], cur[2]) then
        gc:drawLine(prev[1], prev[2], cur[1], cur[2])
      end
      if cur and st ~= "line" and inside(cur[1], cur[2]) then gc:fillRect(cur[1] - 1, cur[2] - 1, 3, 3) end
      prev = cur
    end
  end
  -- leyenda
  local ly = T + 2
  for k, s in ipairs(P.series or {}) do
    if s.name and s.name ~= "" then
      local c = s.col or COL[(k - 1) % #COL + 1]
      local w = gc:getStringWidth(s.name)
      gc:setColorRGB(255, 255, 255); gc:fillRect(R - w - 14, ly, w + 12, 10)
      gc:setColorRGB(c[1], c[2], c[3]); gc:fillRect(R - w - 12, ly + 3, 6, 4)
      gc:drawString(s.name, R - w - 4, ly - 1, "top")
      ly = ly + 10
    end
  end
  -- cursor
  if not P.move and P.series and P.series[P.cs or 1] then
    local s = P.series[P.cs or 1]
    local p = s.pts[P.ci or 1]
    if p then
      local xv, yv = tx(P, p[1]), ty(P, p[2])
      if xv and yv and ok(xv) and ok(yv) then
        gc:setColorRGB(0, 0, 0)
        local cx, cy = X(xv), Y(yv)
        gc:drawRect(cx - 3, cy - 3, 6, 6)
      end
    end
  end
  gc:setFont("sansserif", "r", 9)
end

local function plotInfo(P)
  if type(P.info) == "function" then return P.info(P) end
  if P.move then return P.info or "" end
  local s = P.series and P.series[P.cs or 1]
  if s and s.pts[P.ci or 1] then
    local p = s.pts[P.ci or 1]
    return (s.name or "") .. " x=" .. fmt(p[1]) .. " y=" .. fmt(p[2])
  end
  return P.info or ""
end

-- ---------------------------------------------------------------------
--  Herramientas de dibujo de esquemas (mundo 100 x 60, y hacia arriba)
-- ---------------------------------------------------------------------
local function canvas(gc, bx, by, bw, bh, WW, HH)
  WW, HH = WW or 100, HH or 60
  local s = math.min(bw/WW, bh/HH)
  local ox, oy = bx + (bw - WW*s)/2, by + (bh - HH*s)/2 + HH*s
  local c = {}
  function c.p(x, y) return ox + x*s, oy - y*s end
  function c.col(r, g, b) gc:setColorRGB(r, g, b) end
  function c.line(x1, y1, x2, y2) local a, b = c.p(x1, y1); local d, e = c.p(x2, y2); gc:drawLine(a, b, d, e) end
  function c.dash(x1, y1, x2, y2) gc:setPen("thin", "dashed"); c.line(x1, y1, x2, y2); gc:setPen("thin", "smooth") end
  function c.poly(pts, fill)
    local t = {}
    for i = 1, #pts, 2 do local a, b = c.p(pts[i], pts[i+1]); t[#t+1] = a; t[#t+1] = b end
    t[#t+1] = t[1]; t[#t+1] = t[2]
    if fill then gc:fillPolygon(t) else gc:drawPolyLine(t) end
  end
  function c.rect(x, y, w, h, fill) c.poly({x, y, x + w, y, x + w, y + h, x, y + h}, fill) end
  function c.circle(x, y, r, fill)
    local a, b = c.p(x - r, y + r)
    if fill then gc:fillArc(a, b, 2*r*s, 2*r*s, 0, 360) else gc:drawArc(a, b, 2*r*s, 2*r*s, 0, 360) end
  end
  function c.text(t, x, y, size)
    gc:setFont("sansserif", "r", size or 7)
    local a, b = c.p(x, y); gc:drawString(t, a, b, "top")
    gc:setFont("sansserif", "r", 9)
  end
  function c.ctext(t, x, y, size)
    gc:setFont("sansserif", "r", size or 7)
    local a, b = c.p(x, y); gc:drawString(t, a - gc:getStringWidth(t)/2, b, "top")
    gc:setFont("sansserif", "r", 9)
  end
  function c.arrow(x1, y1, x2, y2)
    c.line(x1, y1, x2, y2)
    local ang = math.atan2(y2 - y1, x2 - x1)
    for _, d in ipairs({2.6, -2.6}) do c.line(x2, y2, x2 - 2.2*math.cos(ang + d*0.2), y2 - 2.2*math.sin(ang + d*0.2)) end
  end
  function c.dim(x1, y1, x2, y2, label)
    c.col(0, 0, 0); c.arrow(x1, y1, x2, y2); c.arrow(x2, y2, x1, y1)
    if label then c.text(label, (x1 + x2)/2 + 1, (y1 + y2)/2 + 2) end
  end
  function c.tri(x, y) -- simbolo nivel de agua
    c.col(0, 60, 160); c.poly({x - 1.4, y + 2.2, x + 1.4, y + 2.2, x, y}, true)
  end
  function c.hatch(x1, x2, y) -- suelo impermeable
    c.col(90, 90, 90); c.line(x1, y, x2, y)
    local x = x1; while x < x2 do c.line(x, y, x - 1.5, y - 1.5); x = x + 2.5 end
  end
  return c
end
local SUELO, AGUA, AGUA2, GRIS, TIERRA = {200,160,70}, {150,195,240}, {60,120,200}, {185,185,190}, {160,110,50}
local function fillc(c, col) c.col(col[1], col[2], col[3]) end

-- ---------------------------------------------------------------------
--  ESQUEMAS
-- ---------------------------------------------------------------------
local DRAW = {}

-- Perfil de alcantarilla segun caso
function DRAW.alcantarilla(caso)
  return function(gc, bx, by, bw, bh)
    local c = canvas(gc, bx, by, bw, bh, 100, 62)
    local function fondo(x) return 20 - 0.12*x end   -- lecho con pendiente SD
    local xa, xb = 34, 78                               -- entrada / salida
    local D = 12
    -- interior del conducto
    fillc(c, {232, 232, 232})
    c.poly({xa, fondo(xa), xb, fondo(xb), xb, fondo(xb) + D, xa, fondo(xa) + D}, true)
    -- agua segun caso
    fillc(c, AGUA)
    local up, dn
    if caso == 0 then up = 34 elseif caso == 1 or caso == 2 then up = 40 elseif caso == 3 then up = 43 else up = 45 end
    -- aguas arriba
    c.poly({2, fondo(2), xa, fondo(xa), xa, up, 2, up}, true)
    -- dentro de la alcantarilla
    local hin
    if caso == 0 then hin = function(x) return fondo(x) + 6 - 2.5*(x - xa)/(xb - xa) end
    elseif caso == 1 then hin = function(x) return fondo(x) + 0.61*D - 2*(x - xa)/(xb - xa) end
    else hin = function(x) return fondo(x) + D end end
    local pts = {}
    for i = 0, 10 do local x = xa + (xb - xa)*i/10; pts[#pts+1] = x; pts[#pts+1] = fondo(x) end
    for i = 10, 0, -1 do local x = xa + (xb - xa)*i/10; pts[#pts+1] = x; pts[#pts+1] = hin(x) end
    c.poly(pts, true)
    -- transicion en la entrada (caso 0: caida hidraulica)
    if caso == 0 then c.poly({xa - 8, fondo(xa - 8), xa, fondo(xa), xa, hin(xa), xa - 8, up}, true) end
    -- aguas abajo
    if caso == 4 then dn = fondo(xb) + D + 3 else dn = fondo(xb) + (caso == 0 and 3 or 4) end
    c.poly({xb, fondo(xb), 98, fondo(98), 98, dn, xb, dn}, true)
    -- terraplen y conducto
    fillc(c, SUELO)
    c.poly({xa, fondo(xa) + D, xa + 5, 50, xb - 5, 48, xb, fondo(xb) + D}, true)
    c.col(0, 0, 0)
    c.line(xa, fondo(xa) + D, xb, fondo(xb) + D)
    c.line(xa, fondo(xa) + D, xa, fondo(xa) + D + 2)
    c.line(xb, fondo(xb) + D, xb, fondo(xb) + D + 2)
    c.col(90, 90, 90)
    c.line(2, fondo(2), 98, fondo(98)); c.line(2, fondo(2) - 1, 98, fondo(98) - 1)
    c.dash(2, fondo(0) - 5, 98, fondo(0) - 5)
    -- superficie libre
    c.col(0, 60, 160)
    c.line(2, up, (caso == 0) and xa - 8 or xa, up)
    if caso == 0 then c.line(xa - 8, up, xa, hin(xa)) end
    c.line(xb, dn, 98, dn)
    c.tri(10, up); c.tri(90, dn)
    -- cotas
    c.dim(8, fondo(8), 8, up, (caso == 0) and "H=hn" or "H'")
    c.dim(xb - 4, fondo(xb - 4), xb - 4, fondo(xb - 4) + D, "D")
    c.col(0, 0, 0); c.text("SD", 20, fondo(20) - 1.5)
    local tit = ({[0] = "Criterio 1: control de entrada (torrente)", "Caso 1: SD>Sn, actua como compuerta",
      "Caso 2: SD=Sn, seccion llena", "Caso 3: SD<Sn, llena con sobrepresion", "Control de salida (ahogada aguas abajo)"})[caso]
    c.col(0, 69, 137); c.text(tit, 2, 62, 8)
    c.col(0, 0, 0)
    if caso == 0 then c.text("Rio", 14, up + 4); c.text("hc", xa + 1, hin(xa) + 1); c.text("Torrente", 52, fondo(52) + 6)
      c.text("H = 1.5 hc ; D = H", 2, 55)
    elseif caso == 1 then c.text("h2 = cc D (cc=0.611)", 44, fondo(44) + 9); c.text("H' = cc D + v^2/2g", 2, 55)
    elseif caso == 2 then c.text("H' = D + (1+ke) v^2/2g", 2, 55)
    elseif caso == 3 then c.text("H' = D + v^2/2g + Le + J L - SD L", 2, 55)
      c.col(0, 60, 160); for k = 0, 3 do c.arrow(xb, fondo(xb) + 2 + 3*k, xb + 7, fondo(xb) + 3 + 3*k) end
    else c.text("h3 = hn3", xb + 2, dn - 1); c.text("H' = D+v^2/2g+Le+Ls+JL-SD L", 2, 55) end
    c.text("(1)", 18, fondo(18) - 5); c.text("(2)", xb - 2, fondo(xb) - 5)
  end
end

function DRAW.tipos(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 62)
  -- terraplen comun
  local function terraplen(x0)
    fillc(c, SUELO); c.poly({x0, 30, x0 + 30, 30, x0 + 30, 18, x0 + 22, 8, x0 + 8, 8, x0, 18}, true)
    c.col(80, 80, 80); c.line(x0 - 1, 30, x0 + 31, 30)
  end
  -- 1 tubo
  terraplen(2)
  fillc(c, {230,230,230}); c.circle(17, 14, 5, true); c.col(0, 0, 0); c.circle(17, 14, 5)
  fillc(c, AGUA); c.rect(13.5, 9.8, 7, 2.5, true)
  c.dim(12, 4, 22, 4, "D")
  c.col(0, 69, 137); c.ctext("1) Tubo hormigon/HDPE", 17, 44, 7); c.ctext("Q pequeno", 17, 39, 7); c.ctext("libre >= D/2", 17, 35, 7)
  -- 2 cajon multiple
  terraplen(35)
  for k = 0, 1 do fillc(c, {230,230,230}); c.rect(40 + k*10.5, 9, 9, 9, true); c.col(0, 0, 0); c.rect(40 + k*10.5, 9, 9, 9)
    fillc(c, AGUA); c.rect(40.5 + k*10.5, 9.5, 8, 3, true) end
  c.col(0, 69, 137); c.ctext("2) Cajon H.A.", 50, 44, 7); c.ctext("(simple o multiple)", 50, 39, 7); c.ctext("luz <= 6 m", 50, 35, 7)
  -- 3 puente
  fillc(c, SUELO); c.poly({68, 30, 70, 30, 76, 8, 68, 8}, true); c.poly({98, 30, 96, 30, 90, 8, 98, 8}, true)
  fillc(c, AGUA); c.rect(76, 8, 14, 4, true)
  fillc(c, GRIS); c.rect(67, 30, 32, 2, true); c.rect(82.5, 8, 1.5, 22, true)
  c.col(0, 0, 0); c.text("estribo", 66, 5); c.text("cepa", 81, 5)
  c.col(0, 69, 137); c.ctext("3) Puente", 83, 44, 7); c.ctext("Q importante", 83, 39, 7); c.ctext("(rios, quebradas)", 83, 35, 7)
  c.col(0, 0, 0); c.text("Tubos y cajones: muros guia en entrada y salida", 4, 58, 7)
  c.text("Alcantarilla = estructura de luz <= 6 m (MC 3.703.101)", 4, 53, 7)
end

function DRAW.entradas(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 50)
  local function tubo(x0, y0)
    fillc(c, AGUA); c.rect(x0, y0, 24, 8, true); c.col(0, 0, 0); c.line(x0, y0 + 8, x0 + 24, y0 + 8); c.line(x0, y0, x0 + 24, y0)
    c.col(0, 60, 160); c.arrow(x0 + 12, y0 + 4, x0 + 20, y0 + 4)
  end
  -- reentrante
  tubo(4, 20); fillc(c, GRIS); c.rect(8, 28, 2, 8, true); c.rect(8, 12, 2, 8, true)
  c.col(0, 0, 0); c.ctext("Tubo saliente", 16, 42); c.ctext("ke = 0.9", 16, 8, 9)
  -- abocinada
  tubo(40, 20); fillc(c, GRIS); c.poly({36, 36, 42, 36, 42, 30, 40, 28}, true); c.poly({36, 12, 42, 12, 42, 18, 40, 20}, true)
  c.col(0, 0, 0); c.ctext("Abocinada", 52, 42); c.ctext("ke = 0.2", 52, 8, 9)
  -- a ras
  tubo(74, 20); fillc(c, GRIS); c.rect(72, 28, 2, 8, true); c.rect(72, 12, 2, 8, true)
  c.col(0, 0, 0); c.ctext("A ras de muro", 86, 42); c.ctext("ke = 0.5", 86, 8, 9)
end

function DRAW.cuneta(tp, b, i, w, i1)
  return function(gc, bx, by, bw, bh)
    local c = canvas(gc, bx, by, bw, bh, 100, 50)
    local sc = 55/math.max(b, 0.3)
    local x0, y0 = 15, 15
    local ym = (tp == 2) and (w*i1 + (b - w)*i) or b*i
    local vf = math.min(4, 18/math.max(ym*sc, 1e-6))
    local function pt(x, y) return x0 + x*sc, y0 + y*sc*vf end   -- escala vertical exagerada
    -- solera
    fillc(c, GRIS); c.rect(x0 - 5, y0 - 8, 5, 8 + math.max(ym*sc*vf*1.4, 6), true)
    local ymax = (tp == 2) and (w*i1 + (b - w)*i) or b*i
    local surf
    if tp == 2 then
      surf = {0, 0, w, w*i1, 1.5*b, w*i1 + (1.5*b - w)*i}
    else surf = {0, 0, 1.5*b, 1.5*b*i} end
    -- agua
    fillc(c, AGUA)
    local wp = {}
    local ax, ay = pt(0, ymax); wp = {x0, y0, }
    if tp == 2 then local a1, a2 = pt(w, w*i1); local a3, a4 = pt(b, ymax); c.poly({x0, y0, a1, a2, a3, a4, ax, ay}, true)
    else local a3, a4 = pt(b, ymax); c.poly({x0, y0, a3, a4, ax, ay}, true) end
    -- calzada
    fillc(c, {150, 150, 155})
    local q = {}
    for k = 1, #surf, 2 do local a, bb = pt(surf[k], surf[k+1]); q[#q+1] = a; q[#q+1] = bb end
    q[#q+1] = q[#q-1]; q[#q+1] = y0 - 8; q[#q+1] = x0; q[#q+1] = y0 - 8
    c.poly(q, true)
    c.col(0, 60, 160); local a, bb = pt(0, ymax); local d, e = pt(b, ymax); c.line(a, bb, d, e); c.tri(a + 3, bb)
    c.dim(a, bb + 6, d, e + 6, "b = " .. fmt(b))
    c.dim(x0 - 7, y0, x0 - 7, bb, "y=" .. fmt(ymax))
    c.col(0, 0, 0)
    if tp == 2 then local p1, p2 = pt(w, 0); c.dash(p1, y0 - 8, p1, p2 + 8); c.text("w=" .. fmt(w), p1 - 8, y0 - 3)
      c.text("i1=" .. fmt(i1), x0 + 2, y0 - 3) end
    c.text("i=" .. fmt(i), 70, 12)
    c.col(0, 69, 137); c.text(tp == 2 and "Cuneta compuesta (esc. vert. exagerada)" or "Cuneta simple (esc. vert. exagerada)", 2, 49, 8)
    c.col(0, 0, 0); c.text(tp == 2 and "Q=0.315 sqrt(S)/n b^(8/3)[i2+(i1-i2)(w/b)^2]^(5/3)" or "Q=sqrt(S)/n b^2 i/2 (b i/(2(i+1)))^(2/3)", 2, 43, 7)
    c.text("Verificar y = b i <= 15 cm, b <= 1-2 m", 2, 38, 7)
  end
end

function DRAW.contrafoso(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 62)
  fillc(c, TIERRA)
  c.poly({2, 5, 2, 50, 12, 48, 22, 44, 30, 42, 30, 35, 37, 35, 38, 40, 45, 38, 52, 32, 58, 22, 64, 12, 68, 8, 68, 5}, true)
  fillc(c, AGUA); c.rect(30.5, 35, 6, 2.5, true)
  fillc(c, AGUA); c.poly({68, 8, 70, 5, 72, 7}, true)
  fillc(c, GRIS); c.poly({68, 5, 70, 5, 72, 7, 88, 9, 98, 8, 98, 2, 68, 2}, true)
  c.col(0, 60, 160)
  for k = 0, 3 do c.arrow(5 + 6*k, 52 - 2*k, 10 + 6*k, 49 - 2*k) end
  c.col(0, 69, 137)
  c.text("Contrafoso", 26, 55, 8); c.arrow(33, 54, 33, 38)
  c.text("Cuneta", 58, 22, 8); c.arrow(66, 20, 70, 8)
  c.text("Bombeo i2=2-4%", 76, 20, 7); c.arrow(85, 17, 86, 10)
  c.text("Calzada", 80, 4, 7)
  c.col(0, 0, 0)
  c.text("Cuneta: w>=0.5 m ; i1>=8% (w>0.5)", 2, 30, 7)
  c.text("Pend. long. min: 0.12% revest.", 2, 25, 7)
  c.text("0.25% sin revestir", 2, 21, 7)
  c.text("Descarga a quebrada/rio", 2, 15, 7)
  c.text("Contrafoso tierra: n=0.023-0.025", 2, 10, 7)
end

function DRAW.sumideros(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 62)
  local function losa(x0)
    fillc(c, {215, 215, 215}); c.poly({x0, 12, x0 + 24, 12, x0 + 30, 26, x0 + 6, 26}, true)
    fillc(c, GRIS); c.poly({x0 + 6, 26, x0 + 30, 26, x0 + 30, 31, x0 + 6, 31}, true)
    c.col(0, 0, 0); c.poly({x0, 12, x0 + 24, 12, x0 + 30, 26, x0 + 6, 26})
  end
  -- fondo
  losa(2); fillc(c, {70, 70, 70}); c.poly({8, 17, 22, 17, 24.5, 23, 10.5, 23}, true)
  c.col(255, 255, 255); for k = 0, 5 do c.line(10 + 2*k, 17.5, 12.5 + 2*k, 22.5) end
  c.col(0, 0, 0); c.ctext("Horizontal de fondo", 16, 42); c.ctext("S3/S4 SERVIU", 16, 37)
  c.ctext("w=0.4-0.7 L=1-2 m", 16, 8); c.ctext("vertedero u orificio", 16, 4)
  -- lateral
  losa(35); fillc(c, {40, 40, 40}); c.rect(47, 26.5, 10, 2.5, true)
  c.col(0, 0, 0); c.ctext("Lateral (ventana)", 50, 42); c.ctext("S2 SERVIU", 50, 37)
  c.ctext("usar si S < 3%", 50, 8); c.ctext("h<a vert. / h>=a ahog.", 50, 4)
  -- mixto
  losa(68); fillc(c, {70, 70, 70}); c.poly({74, 17, 88, 17, 90.5, 23, 76.5, 23}, true)
  fillc(c, {40, 40, 40}); c.rect(80, 26.5, 10, 2.5, true)
  c.col(0, 0, 0); c.ctext("Mixto", 83, 42); c.ctext("S1+S2", 83, 37)
  c.ctext("+10% vs fondo", 83, 8); c.ctext("hojas: lateral emerg.", 83, 4)
  c.col(0, 69, 137); c.text("Sumideros: captan agua de la calle al colector", 2, 58, 8)
  c.text("Qs = min(eta*F*Q , Qm)", 2, 52, 7)
end

function DRAW.pozo(tipo)
  return function(gc, bx, by, bw, bh)
    local c = canvas(gc, bx, by, bw, bh, 100, 62)
    c.hatch(2, 98, 6)
    local top = (tipo == 1) and 26 or 44
    fillc(c, AGUA)
    local function prof(x)
      local r = abs(x - 50); if r < 3 then r = 3 end
      local d = 14*ln(46/r)/ln(46/3)
      return 44 - d
    end
    local pts = {2, 6, 98, 6}
    if tipo == 1 then c.rect(2, 6, 96, 20, true)
    else
      for x = 98, 2, -2 do pts[#pts+1] = x; pts[#pts+1] = prof(x) end
      c.poly(pts, true)
    end
    if tipo == 1 then fillc(c, SUELO); c.rect(2, 26, 96, 4, true); c.col(0, 0, 0); c.hatch(2, 98, 26) end
    -- pozo
    fillc(c, {255, 255, 255}); c.rect(47, 6, 6, 50, true); c.col(0, 0, 0); c.line(47, 6, 47, 56); c.line(53, 6, 53, 56)
    fillc(c, AGUA); c.rect(47, 6, 6, prof(50) - 6, true)
    -- niveles
    c.col(0, 60, 160); c.dash(2, 44, 98, 44); c.text("Nivel original", 60, 50, 7)
    local prev
    for x = 2, 98, 1 do
      if not (x > 47 and x < 53) then
        local y = prof(x)
        if prev and not (prev[1] < 47 and x > 47) then c.line(prev[1], prev[2], x, y) end
        prev = {x, y}
      end
    end
    c.tri(50, prof(50))
    for _, s in ipairs({-1, 1}) do for k = 0, 2 do local y = 9 + 5*k; c.arrow(50 + s*20, y, 50 + s*5, y) end end
    c.dim(75, 6, 75, prof(75), "h")
    c.dim(90, 6, 90, 44, "H")
    c.col(0, 0, 0); c.arrow(50, 58, 55, 60); c.text("Q", 56, 61)
    c.dim(50, 2, 92, 2); c.text("R", 70, 3)
    if tipo == 1 then c.dim(8, 6, 8, 26, "m"); c.col(0, 69, 137)
      c.text("Confinado: H-h = Q/(2 pi T) ln(R/r)", 2, 61, 7)
    else c.col(0, 69, 137); c.text("Libre: H^2-h^2 = Q/(pi K) ln(R/r)", 2, 61, 7) end
  end
end

function DRAW.drenes(tipo)
  return function(gc, bx, by, bw, bh)
    local c = canvas(gc, bx, by, bw, bh, 100, 62)
    c.hatch(2, 98, 6)
    fillc(c, AGUA)
    local yd = (tipo == 1) and 14 or 18
    local pts = {8, 6, 92, 6}
    for x = 92, 8, -2 do local u = (x - 8)/84; pts[#pts+1] = x; pts[#pts+1] = yd + 22*sqrt(math.max(0, u*(1 - u)))*2*0.5 + 0 end
    c.poly(pts, true)
    if tipo == 1 then
      fillc(c, {255,255,255}); c.rect(4, 6, 4, 44, true); c.rect(92, 6, 4, 44, true)
      fillc(c, AGUA); c.rect(4, 6, 4, yd - 6, true); c.rect(92, 6, 4, yd - 6, true)
      c.col(0, 0, 0); c.line(4, 6, 4, 50); c.line(8, 6, 8, 50); c.line(92, 6, 92, 50); c.line(96, 6, 96, 50)
      c.line(0, 50, 4, 50); c.line(8, 50, 92, 50); c.line(96, 50, 100, 50)
      c.dim(12, 6, 12, yd, "h0")
    else
      c.col(0, 0, 0); c.line(0, 50, 100, 50)
      for _, x in ipairs({8, 92}) do fillc(c, {255,255,255}); c.circle(x, yd, 3, true); c.col(0, 0, 0); c.circle(x, yd, 3) end
      c.dash(2, yd, 98, yd); c.dim(94, 6, 94, yd, "d")
    end
    c.col(0, 60, 160)
    for x = 12, 88, 8 do c.arrow(x, 58, x, 51) end
    c.text("f", 50, 61, 8)
    c.dim(50, yd, 50, yd + 11, "H")
    c.dim(8, 3, 92, 3); c.text("D", 50, 4)
    c.col(0, 69, 137)
    c.text(tipo == 1 and "Drenes abiertos: Donnan / Hooghoudt" or "Drenes cerrados: Dagan / Kirkham / Glover", 12, 46, 7)
  end
end

function DRAW.interceptor(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 50)
  c.hatch(2, 98, 6)
  fillc(c, AGUA); c.poly({2, 6, 36, 6, 36, 14, 2, 13}, true); c.rect(36, 6, 12, 8, true)
  local pts = {48, 6, 98, 6, 98, 32}
  for x = 98, 48, -2 do pts[#pts+1] = x; pts[#pts+1] = 14 + 18*(1 - exp(-(x - 48)/10)) end
  c.poly(pts, true)
  c.col(0, 0, 0); c.line(36, 40, 36, 6); c.line(48, 6, 48, 40); c.line(36, 6, 48, 6)
  c.tri(42, 14); c.dim(40, 6, 40, 14, "h0"); c.dim(92, 6, 92, 32, "H")
  c.col(0, 60, 160); c.arrow(80, 20, 70, 20); c.text("q", 74, 23); c.arrow(28, 9, 18, 9); c.text("q'", 22, 12)
  c.dim(48, 2, 92, 2); c.text("x1", 70, 3)
  c.col(0, 69, 137); c.text("H^2 - h0^2 = 2 q x1/K ;  Qz = (q - q') L", 2, 48, 7)
end

function DRAW.costero(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 62)
  fillc(c, AGUA2); c.poly({2, 4, 98, 4, 98, 40, 60, 40}, true)
  fillc(c, AGUA); c.poly({2, 4, 60, 40, 30, 43, 2, 45}, true)
  fillc(c, {170, 200, 250}); c.poly({60, 40, 75, 34, 98, 30, 98, 40}, true)
  fillc(c, SUELO); c.poly({2, 45, 30, 43, 60, 40, 75, 34, 98, 30, 98, 29, 75, 33, 60, 40.5, 30, 47, 2, 50}, true)
  c.col(0, 0, 0); c.hatch(2, 98, 4)
  c.col(0, 60, 160); c.dash(2, 40, 98, 40); c.tri(90, 40); c.text("Mar", 84, 37)
  c.dim(30, 40, 30, 43, "hd"); c.dim(30, 40, 30, 22, "hs")
  c.col(255, 255, 255); c.text("Agua salada", 70, 15)
  c.col(0, 0, 0); c.text("Agua dulce", 6, 30)
  c.col(0, 69, 137); c.text("hs = hd gd/(gs-gd) ~ 40 hd", 2, 60, 8)
end

function DRAW.puente(gc, bx, by, bw, bh)
  local c = canvas(gc, bx, by, bw, bh, 100, 62)
  -- planta: contraccion por terraplenes
  fillc(c, AGUA); c.rect(2, 33, 96, 23, true)
  fillc(c, SUELO); c.rect(40, 49, 14, 7, true); c.rect(40, 33, 14, 7, true)
  c.col(0, 0, 0); c.line(2, 56, 98, 56); c.line(2, 33, 98, 33)
  c.col(210, 30, 30)
  for _, s in ipairs({{30, "(3)"}, {54, "(2)"}, {66, "(1)"}}) do c.dash(s[1], 31, s[1], 58); c.text(s[2], s[1] - 3, 32.5) end
  c.col(0, 60, 160); c.arrow(94, 43, 80, 43); c.text("Q", 86, 47)
  c.dim(8, 33, 8, 56, "b"); c.dim(47, 40, 47, 49, "L")
  c.col(0, 69, 137); c.text("Opcion 1: terraplen y estribos (planta)", 2, 61, 7)
  -- cepas Yarnell
  fillc(c, AGUA); c.rect(2, 2, 96, 20, true)
  fillc(c, GRIS); for k = 0, 2 do c.rect(45, 4 + 6*k, 10, 3, true) end
  c.col(0, 60, 160); for k = 0, 3 do c.line(10, 3.5 + 6*k, 90, 3.5 + 6*k) end
  c.col(0, 0, 0); c.dim(58, 4, 58, 7, "D"); c.text("Opcion 2: cepas (Yarnell)", 4, 27, 7)
  c.text("sigma = nD/b", 70, 27, 7)
end

function DRAW.seccion(tp, b, k, D, y)
  return function(gc, bx, by, bw, bh)
    local c = canvas(gc, bx, by, bw, bh, 100, 62)
    local Bmax
    if tp == 4 then Bmax = D elseif tp == 3 then Bmax = 2*k*y*1.4 else Bmax = b + 2*k*y*1.4 end
    local Hmax = (tp == 4) and D or y*1.4
    local sc = math.min(80/Bmax, 45/Hmax)
    local cx, y0 = 50, 8
    local function P(x, yy) return cx + x*sc, y0 + yy*sc end
    fillc(c, AGUA)
    if tp == 4 then
      local N = 40
      local pts = {}
      local th0 = math.acos(1 - 2*y/D)
      for i = 0, N do local a = -th0 + 2*th0*i/N; local xx, yy = P((D/2)*math.sin(a), D/2 - (D/2)*math.cos(a)); pts[#pts+1] = xx; pts[#pts+1] = yy end
      c.poly(pts, true)
      c.col(0, 0, 0); local xx, yy = P(0, D/2); c.circle(xx, yy, D/2*sc)
    else
      local kk = (tp == 1) and 0 or k
      local bb = (tp == 3) and 0 or b
      local a1, a2 = P(-bb/2, 0); local a3, a4 = P(bb/2, 0); local a5, a6 = P(bb/2 + kk*y, y); local a7, a8 = P(-bb/2 - kk*y, y)
      c.poly({a1, a2, a3, a4, a5, a6, a7, a8}, true)
      c.col(0, 0, 0)
      local t1, t2 = P(-bb/2 - kk*Hmax, Hmax); local t3, t4 = P(bb/2 + kk*Hmax, Hmax)
      if tp == 1 then t1, t2 = P(-bb/2, Hmax); t3, t4 = P(bb/2, Hmax) end
      c.line(t1, t2, a1, a2); c.line(a1, a2, a3, a4); c.line(a3, a4, t3, t4)
      c.dim(a1, a2 - 3, a3, a4 - 3, (tp ~= 3) and ("b=" .. fmt(b)) or nil)
    end
    local _, ys = P(0, y); local xl, _ = P(-Bmax/2, 0)
    c.tri(cx, ys); c.col(0, 0, 0); local xd = math.min(cx + Bmax/2*sc + 3, 86); c.dim(xd, y0, xd, ys, "y=" .. fmt(y))
    c.col(0, 69, 137)
    local A, Pm, B = seccion(tp, b, k, D, y)
    c.text("A=" .. fmt(A) .. "  P=" .. fmt(Pm) .. "  Rh=" .. fmt(A/Pm) .. "  B=" .. fmt(B), 2, 61, 7)
  end
end
-- =====================================================================
--  GRAFICOS ASOCIADOS A EJERCICIOS Y SIMULADORES
-- =====================================================================
local function theisDh(Q, T, Sx, r, t) if t <= 0 then return 0 end return Q/(4*pi*T)*W(r*r*Sx/(4*T*t)) end

-- datos de prueba de gasto: semilog + ajuste
local function plotGasto(t, dh, a, b, t0)
  local fx = {t0 or t[1], t[#t]*1.5}
  local fy = {a*ln(fx[1]) + b, a*ln(fx[2]) + b}
  return {
    {title = "Gasto cte: dh vs t (semilog)", xl = "t", yl = "dh", xlog = true, yinv = true, y0 = true,
     series = {S("datos", t, dh, "pts"), S("recta Jacob", fx, fy, "line")},
     vl = t0 and {{t0, "t0"}} or nil},
    {title = "Gasto cte: dh vs t", xl = "t", yl = "dh", yinv = true, y0 = true,
     series = {S("datos", t, dh, "both")}},
  }
end

local function plotRecup(t, dh, tf, a, b)
  local xs = {}
  for i = 1, #t do xs[i] = t[i]/(t[i] - tf) end
  local mx = 1
  for _, x in ipairs(xs) do if x > mx then mx = x end end
  return {
    {title = "Recuperacion: s' vs t/t'", xl = "t/t'", yl = "s'", xlog = true, yinv = true, y0 = true,
     series = {S("datos", xs, dh, "pts"), S("ajuste", {1, mx*1.2}, {b, a*ln(mx*1.2) + b}, "line")}},
    {title = "Recuperacion: s' vs t", xl = "t", yl = "s'", yinv = true, y0 = true,
     series = {S("datos", t, dh, "both")}, vl = {{tf, "tf"}}},
  }
end

local ATT = {}   -- graficos que se agregan a ejercicios existentes (por nombre)

ATT["Cooper-Jacob con datos"] = function(v)
  local xs = {}; for i = 1, #v.t do xs[i] = ln(v.t[i]) end
  local a, b = regresion(xs, v.dh)
  return plotGasto(v.t, v.dh, a, b, exp(-b/a))
end
ATT["Recuperacion con datos"] = function(v)
  local xs = {}; for i = 1, #v.t do xs[i] = ln(v.t[i]/(v.t[i] - v.tf)) end
  local a, b = regresion(xs, v.dh)
  return plotRecup(v.t, v.dh, v.tf, a, b)
end
ATT["Bombeo variable"] = function(v)
  local tmax = v.t
  for _, x in ipairs(v.ti) do if x*1.3 > tmax then tmax = x*1.3 end end
  local ts, ds = linspace(tmax/300, tmax, 300), {}
  for k, t in ipairs(ts) do
    local s = 0
    for i = 1, #v.Q do if t > v.ti[i] then s = s + v.Q[i]/(4*pi*v.T)*W(v.r^2*v.S/(4*v.T*(t - v.ti[i]))) end end
    ds[k] = s
  end
  local vl = {}; for i = 1, #v.ti do vl[i] = {v.ti[i], "t" .. i} end
  return {{title = "Descenso con bombeo variable (Theis)", xl = "t", yl = "dh", yinv = true, y0 = true,
    series = {S("dh(t)", ts, ds)}, vl = vl}}
end
ATT["Theis (confinado)"] = function(v)
  local ts, ds = logspace(v.t/100, v.t*10, 120), {}
  for k, t in ipairs(ts) do ds[k] = theisDh(v.Q, v.T, v.S, v.r, t) end
  local rs, hs, hn = linspace(-3*sqrt(CJ*v.T*v.t/v.S), 3*sqrt(CJ*v.T*v.t/v.S), 161), {}, {}
  for k, x in ipairs(rs) do hs[k] = -theisDh(v.Q, v.T, v.S, math.max(abs(x), 0.05), v.t) end
  return {{title = "Theis: dh vs t en r=" .. fmt(v.r), xl = "t", yl = "dh", xlog = true, yinv = true, y0 = true,
            series = {S("dh(t)", ts, ds)}, vl = {{v.t, "t"}}},
          {title = "Cono de depresion en t=" .. fmt(v.t), xl = "r", yl = "-dh", series = {S("nivel", rs, hs)}, hl = {{0, "original"}}}}
end
ATT["Cooper-Jacob"] = ATT["Theis (confinado)"]
ATT["Hantush-Jacob (semiconf.)"] = function(v)
  local B = sqrt(v.T*v.mp/v.Kp)
  local ts, d1, d0 = logspace(v.t/100, v.t*100, 100), {}, {}
  for k, t in ipairs(ts) do local u = v.r^2*v.S/(4*v.T*t); d1[k] = v.Q/(4*pi*v.T)*Wh(u, v.r/B); d0[k] = v.Q/(4*pi*v.T)*W(u) end
  return {{title = "Hantush vs Theis", xl = "t", yl = "dh", xlog = true, yinv = true, y0 = true,
    series = {S("Hantush", ts, d1), S("Theis", ts, d0)}}}
end
ATT["Glover-Dumm"] = function(v)
  local al = pi^2*v.K*v.d/(v.S*v.D^2)
  local ts, ys = linspace(0, 3/al, 120), {}
  for k, t in ipairs(ts) do
    local s = 0
    for n = 1, 60 do s = s + 1/(2*n - 1)*math.sin((2*n - 1)*pi/2)*exp(-(2*n - 1)^2*al*t) end
    ys[k] = 4*v.y0/pi*s
  end
  local xs = linspace(0, v.D, 81)
  local P2 = {title = "Perfil de la napa entre drenes", xl = "x", yl = "y", series = {}, y0 = true}
  for _, f in ipairs({0.05, 0.25, 0.5, 1}) do
    local t = f*(v.t > 0 and v.t or 1/al)
    local yy = {}
    for k, x in ipairs(xs) do
      local s = 0
      for n = 1, 60 do s = s + 1/(2*n - 1)*math.sin((2*n - 1)*pi*x/v.D)*exp(-(2*n - 1)^2*al*t) end
      yy[k] = 4*v.y0/pi*s
    end
    P2.series[#P2.series + 1] = S("t=" .. fmt(t), xs, yy)
  end
  return {{title = "Glover-Dumm: y(D/2) vs t", xl = "t", yl = "y", y0 = true, series = {S("serie", ts, ys)},
           hl = {{v.y, "y obj."}}, vl = {{v.t, "t"}}}, P2}
end
ATT["Equilibrio dinamico"] = function(v)
  local al = pi^2*v.K*v.d/(v.S*v.D^2)
  local ts, ys, y, t0 = {}, {}, v.yi, 0
  for i = 1, #v.R do
    local y0 = y + v.R[i]/1000/v.S
    for k = 0, 20 do local t = v.dias*k/20; ts[#ts+1] = t0 + t; ys[#ys+1] = 1.16*y0*exp(-al*t) end
    y = 1.16*y0*exp(-al*v.dias); t0 = t0 + v.dias
  end
  local P = {title = "Equilibrio dinamico: y(t)", xl = "dias", yl = "y", y0 = true, series = {S("y", ts, ys)}}
  if v.ya < 1e8 then P.hl = {{v.ya, "raices"}} end
  local ms, rs = {}, {}
  for i = 1, #v.R do ms[i] = i; rs[i] = v.R[i] end
  return {P, {title = "Recarga mensual", xl = "mes", yl = "R mm", y0 = true, series = {S("R", ms, rs, "both")}}}
end
ATT["Hidrograma SCS"] = function(v)
  local tp = v.tLL/2 + 0.6*v.tc
  local tB = 2.67*tp
  local Qp = 2*v.P*v.A/tB/3.6
  return {{title = "Hidrograma unitario triangular SCS", xl = "t [h]", yl = "Q", y0 = true,
    series = {S("Q", {0, tp, tB}, {0, Qp, 0}, "both")}, vl = {{tp, "tp"}, {tB, "tB"}}}}
end
ATT["L minimo (momenta)"] = function(v)
  local M1 = v.Q^2/(G*v.b*v.h1) + v.b*v.h1^2/2
  local function M2c(L) local hc = (v.Q^2/(G*v.b*L))^(1/3); return v.Q^2/(G*L*hc) + v.b*hc^2/2 end
  local L = largest(roots(function(L) return M2c(L) - M1 end, 1e-4*v.b, v.b, 2000, true)) or v.b
  local hs = linspace(0.05*v.h1, 2.5*v.h1, 120)
  local m1, m2 = {}, {}
  for k, h in ipairs(hs) do m1[k] = v.Q^2/(G*v.b*h) + v.b*h*h/2; m2[k] = v.Q^2/(G*L*h) + v.b*h*h/2 end
  local P = {title = "Momenta: seccion 1 y 2 (L=" .. fmt(L) .. ")", xl = "M", yl = "h", series = {}}
  P.series[1] = {name = "M1", pts = {}}; P.series[2] = {name = "M2", pts = {}}
  for k = 1, #hs do P.series[1].pts[k] = {m1[k], hs[k]}; P.series[2].pts[k] = {m2[k], hs[k]} end
  P.hl = {{v.h1, "hn1"}}; P.vl = {{M1, "M1"}}
  P.xr = {0, M1*2.5}
  return {P}
end
ATT["Canal: Manning yn, yc"] = function(v)
  local tp, b, k, D = v.tp, v.b, v.k, v.D
  local lim = (tp == 4) and D*0.9999 or 1e4
  local yn = roots(function(y) local A, Pm = seccion(tp, b, k, D, y); return v.Q*v.n/sqrt(v.S) - A^(5/3)/Pm^(2/3) end, 1e-6, lim, 2000, true)[1]
  local yc = roots(function(y) local A, _, B = seccion(tp, b, k, D, y); return v.Q^2/G - A^3/B end, 1e-6, lim, 2000, true)[1]
  local out = {}
  if yn then out[#out+1] = {title = "Seccion con altura normal", draw = DRAW.seccion(tp, b, k, D, yn)} end
  if yc then
    local ymax = (tp == 4) and 0.999*D or math.max(yn or 0, yc)*3
    local ys, Es = linspace(yc*0.35, ymax, 150), {}
    for i, y in ipairs(ys) do local A = seccion(tp, b, k, D, y); Es[i] = y + v.Q^2/(2*G*A*A) end
    local P = {title = "Energia especifica E(y)", xl = "E", yl = "y", series = {{name = "E", pts = {}}}, hl = {{yc, "yc"}}}
    for i = 1, #ys do P.series[1].pts[i] = {Es[i], ys[i]} end
    if yn then P.hl[2] = {yn, "yn"} end
    out[#out+1] = P
  end
  return out
end
ATT["Superposicion (x,y)"] = function(v)
  local pe, pi_ = {}, {}
  for i = 1, #v.Q do if v.Q[i] >= 0 then pe[#pe+1] = {v.X[i], v.Y[i]} else pi_[#pi_+1] = {v.X[i], v.Y[i]} end end
  return {{title = "Planta: pozos (Q>0 azul, Q<0 rojo) y punto", xl = "x", yl = "y",
    series = {{name = "extrae", pts = pe, style = "pts"}, {name = "inyecta", pts = pi_, style = "pts"},
              {name = "P", pts = {{v.xp, v.yp}}, style = "pts", col = {0, 140, 60}}}}}
end
ATT["Imagen: 1 pozo y borde"] = function(v)
  local xs = linspace(0.01, 3*v.b + abs(v.xp), 200)
  local d = {}
  for k, x in ipairs(xs) do
    local rr = math.max(sqrt((x - v.b)^2 + v.yp^2), 0.05)
    local ri = sqrt((x + v.b)^2 + v.yp^2)
    local dr = (rr < v.R) and v.Q/(2*pi*v.T)*ln(v.R/rr) or 0
    local di = (ri < v.R) and v.Q/(2*pi*v.T)*ln(v.R/ri) or 0
    d[k] = (v.tp == 1) and (dr - di) or (dr + di)
  end
  return {{title = (v.tp == 1 and "Recarga lineal" or "Barrera") .. ": dh a lo largo de x (y=" .. fmt(v.yp) .. ")",
    xl = "x (borde en 0)", yl = "dh", yinv = true, series = {S("dh", xs, d)}, vl = {{0, "borde"}, {v.b, "pozo"}}}}
end
ATT["Dupuit (libre)"] = function(v)
  local rs, hs = linspace(-v.R, v.R, 161), {}
  for k, x in ipairs(rs) do local r = math.max(abs(x), 0.05); hs[k] = (r < v.R) and sqrt(math.max(0, v.H^2 - v.Q/(pi*v.K)*ln(v.R/r))) or v.H end
  return {{title = "Cono de depresion (Dupuit)", xl = "r", yl = "h", y0 = true, series = {S("h(r)", rs, hs)}, hl = {{v.H, "H"}}}}
end
ATT["Thiem (confinado)"] = function(v)
  local rs, hs = linspace(-v.R, v.R, 161), {}
  for k, x in ipairs(rs) do local r = math.max(abs(x), 0.05); hs[k] = (r < v.R) and -v.Q/(2*pi*v.T)*ln(v.R/r) or 0 end
  return {{title = "Cono de depresion (Thiem)", xl = "r", yl = "-dh", series = {S("-dh(r)", rs, hs)}, hl = {{0, "original"}}}}
end
ATT["Funcion de pozo W(u)"] = function(v)
  local iu, w = logspace(0.1, 1e4, 120), {}
  for k, x in ipairs(iu) do w[k] = W(1/x) end
  return {{title = "Curva tipo de Theis W(u)", xl = "1/u", yl = "W(u)", xlog = true, ylog = true, yr = {0.01, 10},
    series = {S("W", iu, w)}, vl = {{1/v.u, "u"}}}}
end
ATT["Yarnell (cepas)"] = function() return {{title = "Esquema", draw = DRAW.puente}} end

-- --------------------- Simuladores y graficos nuevos -----------------
local GRAF = {
  P("Simular gasto + recuperacion", "Theis: dh = Q/(4piT)[W(u) - W(u')], u'=r^2S/(4T(t-tb)) para t>tb",
    {{"Q","Q","Caudal"},{"T","T","Transmisibilidad"},{"S","S","Coef. almacenamiento"},{"r","r","Distancia (radio pozo)"},
     {"tb","t bombeo","Duracion del bombeo"},{"tt","t total","Tiempo total (> tb)"}},
    function(v)
      local dtb = theisDh(v.Q, v.T, v.S, v.r, v.tb)
      local res = theisDh(v.Q, v.T, v.S, v.r, v.tt) - theisDh(v.Q, v.T, v.S, v.r, v.tt - v.tb)
      local ds = 2.303*v.Q/(4*pi*v.T)
      return {"dh al fin del bombeo = " .. fmt(dtb), "dh residual en t total = " .. fmt(res),
              "Pendiente por ciclo log = 2.3Q/(4piT) = " .. fmt(ds), "R(tb) = " .. fmt(sqrt(CJ*v.T*v.tb/v.S)),
              "[enter] ver graficos, [tab] cambia grafico"}
    end),
  P("Simular prueba escalonada", "Gasto variable: escalones de igual duracion, superposicion Theis",
    {{"Q","Q1,Q2,..","Caudales de cada escalon", list=true},{"dur","duracion","Duracion de cada escalon"},
     {"T","T","Transmisibilidad"},{"S","S","Coef. almac."},{"r","r","Radio del pozo"}},
    function(v)
      local o = {}
      for i = 1, #v.Q do
        local t = i*v.dur
        local s = 0
        for j = 1, i do local dQ = v.Q[j] - (v.Q[j-1] or 0); s = s + theisDh(dQ, v.T, v.S, v.r, t - (j - 1)*v.dur) end
        o[#o+1] = "Esc." .. i .. " Q=" .. fmt(v.Q[i]) .. " dh=" .. fmt(s) .. " Q/dh=" .. fmt(v.Q[i]/s)
      end
      o[#o+1] = "Prueba gasto cte: usar 80% Qmax = " .. fmt(0.8*v.Q[#v.Q])
      return o
    end),
  P("Prueba escalonada (datos)", "dh = B Q + C Q^2 : regresion dh/Q vs Q",
    {{"Q","Q1,Q2,..","Caudales (lista)", list=true},{"dh","dh1,dh2,..","Descensos estabilizados", list=true}},
    function(v)
      local y = {}; for i = 1, #v.Q do y[i] = v.dh[i]/v.Q[i] end
      local C, B = regresion(v.Q, y)
      local o = {"B (acuifero) = " .. fmt(B), "C (perdidas pozo) = " .. fmt(C)}
      for i = 1, #v.Q do o[#o+1] = "Q=" .. fmt(v.Q[i]) .. " cap.esp.=" .. fmt(v.Q[i]/v.dh[i]) end
      o[#o+1] = "80% Qmax = " .. fmt(0.8*v.Q[#v.Q])
      return o
    end),
  P("Prueba completa (datos)", "Bombeo (Cooper-Jacob) + recuperacion en una sola figura",
    {{"tb","t bombeo","Tiempos de bombeo (lista)", list=true},{"db","dh bombeo","Descensos (lista)", list=true},
     {"tr","t recup.","Tiempos desde inicio (>tf) (lista)", list=true},{"dr","s' recup.","Desc. residuales (lista)", list=true},
     {"tf","tf","Fin del bombeo"},{"Q","Q","Caudal"},{"r","r","Distancia"}},
    function(v)
      local xs = {}; for i = 1, #v.tb do xs[i] = ln(v.tb[i]) end
      local a, b = regresion(xs, v.db)
      local T = v.Q/(4*pi*a); local t0 = exp(-b/a)
      local xr = {}; for i = 1, #v.tr do xr[i] = ln(v.tr[i]/(v.tr[i] - v.tf)) end
      local a2 = regresion(xr, v.dr)
      return {"Bombeo: T = " .. fmt(T) .. "  S = " .. fmt(CJ*T*t0/v.r^2), "Recuperacion: T = " .. fmt(v.Q/(4*pi*a2)),
              "[enter] ver graficos"}
    end),
  P("Curva tipo interactiva", "Mueva la curva W(u) con flechas hasta calzar los datos. Punto de ajuste W=1, 1/u=1",
    {{"t","t1,t2,..","Tiempos (lista)", list=true},{"dh","dh1,dh2,..","Descensos (lista)", list=true},
     {"Q","Q","Caudal"},{"r","r","Distancia al pozo obs."}},
    function(v) return {"[enter] abre el grafico interactivo.", "Flechas: mover curva. Teclas * y / cambian el paso.",
                        "T y S se actualizan en la parte superior."} end),
  P("Cono de depresion", "Perfil del nivel para varios tiempos (Theis) o permanente",
    {{"tipo","Tipo 1-3","1 Thiem conf., 2 Dupuit libre, 3 Theis (tiempos)"},{"Q","Q","Caudal"},{"T","T o K","Transmisibilidad (1,3) o K (2)"},
     {"H","H","Nivel/espesor inicial", opt="0"},{"R","R","Radio de influencia (1,2)", opt="100"},
     {"S","S","Coef. almac. (3)", opt="1E-4"},{"t","t1,t2,..","Tiempos (3)", list=true, opt="3600,86400"}},
    function(v) return {"[enter] ver grafico"} end),
  P("Curvas de Hantush", "W(u, r/B) para varios r/B",
    {{"rb","r/B lista","Valores de r/B", list=true, opt="0.01,0.1,0.5,1,2"}},
    function(v) return {"[enter] ver familia de curvas"} end),
}
local GRAFPLOT = {}
GRAFPLOT["Simular gasto + recuperacion"] = function(v)
  local ts = linspace(v.tt/400, v.tt, 400)
  local ds = {}
  for k, t in ipairs(ts) do
    ds[k] = theisDh(v.Q, v.T, v.S, v.r, t) - ((t > v.tb) and theisDh(v.Q, v.T, v.S, v.r, t - v.tb) or 0)
  end
  local tp, dp = logspace(v.tb/1000, v.tb, 80), {}
  for k, t in ipairs(tp) do dp[k] = theisDh(v.Q, v.T, v.S, v.r, t) end
  local cj = {}; for k, t in ipairs(tp) do cj[k] = v.Q/(4*pi*v.T)*ln(CJ*v.T*t/(v.r^2*v.S)) end
  local xr, dr = {}, {}
  for _, t in ipairs(logspace(v.tb*1.001, v.tt, 80)) do
    xr[#xr+1] = t/(t - v.tb); dr[#dr+1] = theisDh(v.Q, v.T, v.S, v.r, t) - theisDh(v.Q, v.T, v.S, v.r, t - v.tb)
  end
  return {
    {title = "Gasto constante + recuperacion", xl = "t", yl = "dh", yinv = true, y0 = true, series = {S("dh(t)", ts, ds)}, vl = {{v.tb, "tb"}}},
    {title = "Bombeo semilog (Jacob)", xl = "t", yl = "dh", xlog = true, yinv = true, y0 = true,
     series = {S("Theis", tp, dp), S("Jacob", tp, cj)}},
    {title = "Recuperacion: s' vs t/t'", xl = "t/t'", yl = "s'", xlog = true, yinv = true, y0 = true, series = {S("s'", xr, dr)}},
  }
end
GRAFPLOT["Simular prueba escalonada"] = function(v)
  local n = #v.Q
  local ts = linspace(v.dur*n/400, v.dur*n, 400)
  local ds = {}
  for k, t in ipairs(ts) do
    local s = 0
    for j = 1, n do local tj = (j - 1)*v.dur; if t > tj then s = s + theisDh(v.Q[j] - (v.Q[j-1] or 0), v.T, v.S, v.r, t - tj) end end
    ds[k] = s
  end
  local qe, de = {}, {}
  for i = 1, n do
    local t = i*v.dur; local s = 0
    for j = 1, i do s = s + theisDh(v.Q[j] - (v.Q[j-1] or 0), v.T, v.S, v.r, t - (j - 1)*v.dur) end
    qe[i] = v.Q[i]; de[i] = s
  end
  local vl = {}; for i = 1, n - 1 do vl[i] = {i*v.dur} end
  local qq, tq = {}, {}
  for i = 1, n do qq[#qq+1] = v.Q[i]; tq[#tq+1] = (i - 1)*v.dur; qq[#qq+1] = v.Q[i]; tq[#tq+1] = i*v.dur end
  return {
    {title = "Prueba escalonada: dh(t)", xl = "t", yl = "dh", yinv = true, y0 = true, series = {S("dh", ts, ds)}, vl = vl},
    {title = "Caudal Q(t)", xl = "t", yl = "Q", y0 = true, series = {S("Q", tq, qq)}},
    {title = "Curva caracteristica dh vs Q", xl = "Q", yl = "dh", yinv = true, y0 = true, series = {S("dh", qe, de, "both")}},
  }
end
GRAFPLOT["Prueba escalonada (datos)"] = function(v)
  local y = {}; for i = 1, #v.Q do y[i] = v.dh[i]/v.Q[i] end
  local C, B = regresion(v.Q, y)
  local qs = linspace(0, v.Q[#v.Q]*1.3, 60); local fs = {}
  for k, q in ipairs(qs) do fs[k] = B*q + C*q*q end
  return {{title = "dh vs Q", xl = "Q", yl = "dh", yinv = true, y0 = true, series = {S("datos", v.Q, v.dh, "pts"), S("BQ+CQ^2", qs, fs)}},
          {title = "dh/Q vs Q (Jacob)", xl = "Q", yl = "dh/Q", series = {S("datos", v.Q, y, "pts"), S("recta", {0, qs[#qs]}, {B, B + C*qs[#qs]})}}}
end
GRAFPLOT["Prueba completa (datos)"] = function(v)
  local xs = {}; for i = 1, #v.tb do xs[i] = ln(v.tb[i]) end
  local a, b = regresion(xs, v.db)
  local xr = {}; for i = 1, #v.tr do xr[i] = ln(v.tr[i]/(v.tr[i] - v.tf)) end
  local a2, b2 = regresion(xr, v.dr)
  local tt, dd = {}, {}
  for i = 1, #v.tb do tt[#tt+1] = v.tb[i]; dd[#dd+1] = v.db[i] end
  for i = 1, #v.tr do tt[#tt+1] = v.tr[i]; dd[#dd+1] = v.dr[i] end
  local out = {{title = "Prueba completa dh(t)", xl = "t", yl = "dh", yinv = true, y0 = true, series = {S("datos", tt, dd, "both")}, vl = {{v.tf, "tf"}}}}
  for _, p in ipairs(plotGasto(v.tb, v.db, a, b, exp(-b/a))) do out[#out+1] = p end
  out[#out+1] = plotRecup(v.tr, v.dr, v.tf, a2, b2)[1]
  return out
end
GRAFPLOT["Curva tipo interactiva"] = function(v)
  local iu = logspace(0.1, 1e4, 100)
  local P = {title = "Ajuste curva tipo (log-log)", xl = "t", yl = "dh", xlog = true, ylog = true,
             series = {S("datos", v.t, v.dh, "pts"), {name = "W(u)", pts = {}, style = "line"}}}
  -- posicion inicial: alinear centros
  local mt, md = 0, 0
  for i = 1, #v.t do mt = mt + log10(v.t[i]); md = md + log10(v.dh[i]) end
  mt, md = mt/#v.t, md/#v.t
  P.ox, P.oy, P.step = mt - 1.5, md - 0.4, 0.05
  local x0, x1, y0, y1 = 1e30, -1e30, 1e30, -1e30
  for i = 1, #v.t do x0 = math.min(x0, v.t[i]); x1 = math.max(x1, v.t[i]); y0 = math.min(y0, v.dh[i]); y1 = math.max(y1, v.dh[i]) end
  P.xr = {x0/10, x1*10}; P.yr = {y0/10, y1*10}
  local function upd()
    local pts = {}
    for k, x in ipairs(iu) do pts[k] = {x*10^P.ox, W(1/x)*10^P.oy} end
    P.series[2].pts = pts
    local dhs, ts = 10^P.oy, 10^P.ox
    local T = v.Q/(4*pi*dhs)
    P.info = "dh*=" .. fmt(dhs) .. " t*=" .. fmt(ts) .. " T=" .. fmt(T) .. " S=" .. fmt(4*T*ts/v.r^2)
  end
  P.move = function(P_, dx, dy) P.ox = P.ox + dx*P.step; P.oy = P.oy + dy*P.step; upd() end
  upd()
  return {P}
end
GRAFPLOT["Cono de depresion"] = function(v)
  local P = {xl = "r", series = {}}
  if v.tipo == 3 then
    local Rm = sqrt(CJ*v.T*v.t[#v.t]/v.S)
    local rs = linspace(-1.5*Rm, 1.5*Rm, 161)
    P.title = "Cono de depresion (Theis)"; P.yl = "nivel"
    for _, t in ipairs(v.t) do
      local hs = {}
      for k, x in ipairs(rs) do hs[k] = v.H - theisDh(v.Q, v.T, v.S, math.max(abs(x), 0.05), t) end
      P.series[#P.series+1] = S("t=" .. fmt(t), rs, hs)
    end
  else
    local rs = linspace(-1.2*v.R, 1.2*v.R, 161); local hs = {}
    for k, x in ipairs(rs) do
      local r = math.max(abs(x), 0.05)
      if r >= v.R then hs[k] = v.H
      elseif v.tipo == 1 then hs[k] = v.H - v.Q/(2*pi*v.T)*ln(v.R/r)
      else hs[k] = sqrt(math.max(0, v.H^2 - v.Q/(pi*v.T)*ln(v.R/r))) end
    end
    P.title = v.tipo == 1 and "Cono (Thiem, confinado)" or "Cono (Dupuit, libre)"; P.yl = "h"
    P.series[1] = S("h(r)", rs, hs)
  end
  P.hl = {{v.H, "inicial"}}
  return {P}
end
GRAFPLOT["Curvas de Hantush"] = function(v)
  local iu = logspace(0.1, 1e4, 80)
  local P = {title = "Hantush W(u, r/B)", xl = "1/u", yl = "W", xlog = true, ylog = true, yr = {0.01, 10}, series = {}}
  for _, rb in ipairs(v.rb) do
    local w = {}; for k, x in ipairs(iu) do w[k] = Wh(1/x, rb) end
    P.series[#P.series+1] = S("r/B=" .. fmt(rb), iu, w)
  end
  return P and {P}
end
for _, it in ipairs(GRAF) do it.plot = GRAFPLOT[it.name] end
CATS[#CATS + 1] = {name = "GRAFICOS y simuladores", items = GRAF}

-- --------------------- Esquemas -------------------------------------
local function Dg(name, desc, fields, drawf)
  return {kind = "D", name = name, f = desc, fields = fields or {}, drawf = drawf}
end
CATS[#CATS + 1] = {name = "ESQUEMAS (como son)", items = {
  Dg("Tipos de atraviesos", "Tubo, cajon y puente", nil, function() return DRAW.tipos end),
  Dg("Alcantarilla: casos", "Perfil longitudinal segun el caso de diseno",
     {{"caso","Caso 0-4","0 control entrada,1 SD>Sn,2 SD=Sn,3 SD<Sn,4 control salida", opt="0"}},
     function(v) return DRAW.alcantarilla(math.floor(v.caso)) end),
  Dg("Alcantarilla: 5 casos", "Recorra los 5 casos con [tab]", nil, "todos"),
  Dg("Entradas y coef. ke", "Tipos de entrada a la alcantarilla", nil, function() return DRAW.entradas end),
  Dg("Seccion hidraulica", "Dibuja la seccion con agua a altura y",
     {{"tp","Tipo 1-4","1 rect 2 trapecio 3 triangular 4 circular", opt="4"},{"b","b","Ancho basal", opt="1"},
      {"k","k","Talud H:V", opt="1"},{"D","D","Diametro", opt="1"},{"y","y","Altura de agua", opt="0.7"}},
     function(v) return DRAW.seccion(math.floor(v.tp), v.b, v.k, v.D, v.y) end),
  Dg("Esquema cuneta simple", "Seccion a escala",
     {{"b","b","Ancho inundado [m]", opt="1"},{"i","i","Bombeo", opt="0.04"}},
     function(v) return DRAW.cuneta(1, v.b, v.i, 0, 0) end),
  Dg("Esquema cuneta compuesta", "Seccion a escala",
     {{"b","b","Ancho inundado [m]", opt="1.5"},{"i2","i2","Bombeo", opt="0.03"},{"w","w","Ancho zarpa", opt="0.6"},{"i1","i1","Pend. zarpa", opt="0.08"}},
     function(v) return DRAW.cuneta(2, v.b, v.i2, v.w, v.i1) end),
  Dg("Contrafoso y cuneta", "Camino a media ladera (perfil transversal)", nil, function() return DRAW.contrafoso end),
  Dg("Sumideros", "Fondo, lateral y mixto", nil, function() return DRAW.sumideros end),
  Dg("Pozo confinado", "Thiem", nil, function() return DRAW.pozo(1) end),
  Dg("Pozo libre", "Dupuit", nil, function() return DRAW.pozo(2) end),
  Dg("Drenes abiertos", "Zanjas paralelas", nil, function() return DRAW.drenes(1) end),
  Dg("Drenes cerrados", "Tuberias paralelas", nil, function() return DRAW.drenes(2) end),
  Dg("Esquema dren interceptor", "Corte transversal", nil, function() return DRAW.interceptor end),
  Dg("Acuifero costero", "Ghyben-Herzberg", nil, function() return DRAW.costero end),
  Dg("Puente", "Contraccion (planta) y cepas", nil, function() return DRAW.puente end),
}}

for _, c in ipairs(CATS) do for _, it in ipairs(c.items) do if ATT[it.name] then it.plot = ATT[it.name] end end end
-- =====================================================================
--  Resolucion de un ejercicio
-- =====================================================================
local function leer(item, vals)
  local v, blanks = {}, {}
  for _, fd in ipairs(item.fields) do
    local s = vals[fd[1]] or ""
    if s:match("^%s*$") and fd.opt then s = fd.opt end
    if s:match("^%s*$") then blanks[#blanks + 1] = fd
    elseif fd.list then v[fd[1]] = evallist(s)
    else v[fd[1]] = evalstr(s) end
  end
  return v, blanks
end

local function resolver(item, vals)
  local v, blanks = leer(item, vals)
  if item.kind == "P" or item.kind == "D" then
    if #blanks > 0 then error("Falta dato: " .. blanks[1][2]) end
    if item.kind == "D" then return {}, v end
    return item.calc(v), v
  end
  local out = {}
  if #blanks == 0 then
    local r = item.res(v)
    out[#out + 1] = "Todos los datos ingresados."
    out[#out + 1] = "Residuo (lhs-rhs) = " .. fmt(r)
    out[#out + 1] = "Deje UN campo vacio para despejarlo."
  elseif #blanks > 1 then
    error("Deje solo UN dato vacio (hay " .. #blanks .. ")")
  else
    local id = blanks[1][1]
    local rs = solveAll(function(x) v[id] = x; return item.res(v) end, item.neg)
    if #rs == 0 then error("No se encontro solucion para " .. blanks[1][2]) end
    for k, x in ipairs(rs) do
      out[#out + 1] = blanks[1][2] .. (#rs > 1 and (" (sol " .. k .. ")") or "") .. " = " .. fmt(x)
    end
    if #rs > 1 then out[#out + 1] = "(Varias raices: elija la fisica)" end
    v[id] = rs[1]
  end
  if item.extra then
    local okx, ex = pcall(item.extra, v)
    if okx and ex then for _, l in ipairs(ex) do out[#out + 1] = l end end
  end
  return out, v
end

local function graficos(item, v)
  if item.kind == "D" then
    if item.drawf == "todos" then
      local t = {}
      for k = 0, 4 do t[#t+1] = {title = "Alcantarilla caso " .. k, draw = DRAW.alcantarilla(k)} end
      return t
    end
    return {{title = item.name, draw = item.drawf(v)}}
  end
  if item.plot then
    local okp, p = pcall(item.plot, v)
    if okp and p and #p > 0 then return p end
    if not okp then error(p) end
  end
  return nil
end

-- exportar para pruebas fuera de la calculadora
HID = {CATS = CATS, resolver = resolver, graficos = graficos, drawPlot = drawPlot, plotInfo = plotInfo,
       W = W, Wh = Wh, Fk = Fk, evalstr = evalstr}

-- =====================================================================
--  Interfaz grafica (TI-Nspire)
-- =====================================================================
local scr = "main"      -- main | sub | form | res | graf
local ci, ii, fi = 1, 1, 1
local top = 1
local item
local vals = {}
local reslines, rtop = {}, 1
local plots, pk = nil, 1
local grafBack = "res"
local LH = 15

local function inv() platform.window:invalidate() end
local function W_() return platform.window:width() end
local function H_() return platform.window:height() end
local function nvis() return math.floor((H_() - 36) / LH) end

local function wrap(gc, text, width)
  local out, line = {}, ""
  for word in text:gmatch("%S+") do
    local test = (line == "") and word or (line .. " " .. word)
    if gc:getStringWidth(test) > width and line ~= "" then out[#out + 1] = line; line = word
    else line = test end
  end
  if line ~= "" then out[#out + 1] = line end
  return out
end

local function header(gc, title)
  gc:setColorRGB(0, 69, 137)
  gc:fillRect(0, 0, W_(), 18)
  gc:setColorRGB(255, 255, 255)
  gc:setFont("sansserif", "b", 10)
  gc:drawString(title, 4, 1, "top")
  gc:setFont("sansserif", "r", 9)
end
local function footer(gc, text)
  gc:setColorRGB(230, 230, 230)
  gc:fillRect(0, H_() - 14, W_(), 14)
  gc:setColorRGB(60, 60, 60)
  gc:setFont("sansserif", "r", 7)
  gc:drawString(text, 3, H_() - 13, "top")
  gc:setFont("sansserif", "r", 9)
end
local function drawList(gc, lines, sel)
  local n = nvis()
  if sel then
    if sel < top then top = sel end
    if sel > top + n - 1 then top = sel - n + 1 end
  end
  for k = 0, n - 1 do
    local idx = top + k
    local l = lines[idx]
    if not l then break end
    local y = 20 + k * LH
    if idx == sel then
      gc:setColorRGB(247, 176, 0); gc:fillRect(0, y, W_(), LH)
    end
    gc:setColorRGB(0, 0, 0)
    gc:drawString(l, 4, y, "top")
  end
end

local MARCA = {E = "= ", P = "> ", D = "# "}

function on.paint(gc)
  gc:setPen("thin", "smooth")
  gc:setFont("sansserif", "r", 9)
  if scr == "main" then
    header(gc, "CIV-346 Hidraulica Aplicada")
    local l = {}
    for k, c in ipairs(CATS) do l[k] = k .. ". " .. c.name end
    drawList(gc, l, ci)
    footer(gc, "enter: abrir   " .. #CATS .. " secciones")
  elseif scr == "sub" then
    header(gc, CATS[ci].name)
    local l = {}
    for k, it in ipairs(CATS[ci].items) do l[k] = (MARCA[it.kind] or "") .. it.name .. (it.plot and "  [graf]" or "") end
    drawList(gc, l, ii)
    footer(gc, "= despeja  > calcula  # esquema  [graf] tiene grafico")
  elseif scr == "form" then
    header(gc, item.name)
    local fl = wrap(gc, item.f, W_() - 8)
    local y = 20
    gc:setColorRGB(0, 69, 137)
    for k = 1, math.min(#fl, 2) do gc:drawString(fl[k], 4, y, "top"); y = y + 12 end
    local lines = {}
    for k, fd in ipairs(item.fields) do
      local s = vals[fd[1]] or ""
      if s == "" and fd.opt then s = "(" .. fd.opt .. ")" end
      if s == "" then s = "?" end
      lines[k] = fd[2] .. " = " .. s
    end
    local n = math.floor((H_() - y - 30) / LH)
    if fi < top then top = fi end
    if fi > top + n - 1 then top = fi - n + 1 end
    for k = 0, n - 1 do
      local idx = top + k
      if not lines[idx] then break end
      local yy = y + k * LH
      if idx == fi then gc:setColorRGB(247, 176, 0); gc:fillRect(0, yy, W_(), LH) end
      gc:setColorRGB(0, 0, 0)
      gc:drawString(lines[idx] .. ((idx == fi) and "_" or ""), 4, yy, "top")
    end
    gc:setColorRGB(90, 90, 90)
    if item.fields[fi] then gc:drawString(item.fields[fi][3], 4, H_() - 29, "top") end
    footer(gc, item.kind == "E" and "Deje vacio el dato a calcular. enter: resolver" or "enter: calcular/dibujar   esc: volver")
  elseif scr == "res" then
    header(gc, "Resultado: " .. item.name)
    local wl = {}
    for _, l in ipairs(reslines) do for _, w in ipairs(wrap(gc, l, W_() - 8)) do wl[#wl + 1] = w end end
    local n = nvis()
    if rtop > math.max(1, #wl - n + 1) then rtop = math.max(1, #wl - n + 1) end
    for k = 0, n - 1 do
      local l = wl[rtop + k]
      if not l then break end
      gc:setColorRGB(0, 0, 0)
      gc:drawString(l, 4, 20 + k * LH, "top")
    end
    footer(gc, plots and "enter: VER GRAFICO   flechas: desplazar   esc: volver" or "flechas: desplazar   esc/enter: volver")
  elseif scr == "graf" then
    local P = plots[pk]
    header(gc, (P.title or item.name) .. ((#plots > 1) and (" (" .. pk .. "/" .. #plots .. ")") or ""))
    gc:setColorRGB(255, 255, 255); gc:fillRect(0, 18, W_(), H_() - 18)
    local info = plotInfo(P)
    local by = 19
    if info ~= "" and not P.draw then
      gc:setFont("sansserif", "r", 8); gc:setColorRGB(150, 0, 0); gc:drawString(info, 3, 18, "top")
      gc:setFont("sansserif", "r", 9); by = 30
    end
    drawPlot(gc, P, 0, by, W_(), H_() - by - 14)
    if P.move then footer(gc, "flechas: mover curva  * / : paso  tab: sig.  esc: volver")
    elseif P.draw then footer(gc, "tab/enter: siguiente   esc: volver")
    else footer(gc, "<- -> cursor   ^v serie   tab/enter: sig. grafico   esc: volver") end
  end
end

local function openItem()
  item = CATS[ci].items[ii]
  vals = {}
  fi, top = 1, 1
  plots = nil
  if item.kind == "D" and #item.fields == 0 then
    local okc, p = pcall(graficos, item, {})
    if okc and p then plots, pk, scr, grafBack = p, 1, "graf", "sub"; return end
  end
  scr = "form"
end
local function compute()
  plots = nil
  local okc, r, v = pcall(resolver, item, vals)
  if okc then
    reslines = r
    local okg, p = pcall(graficos, item, v)
    if okg then plots = p else reslines[#reslines + 1] = "Grafico no disponible: " .. tostring(p):gsub("^.-:%d+: ", "") end
  else reslines = {"ERROR:", (tostring(r):gsub("^.-:%d+: ", ""))} end
  rtop = 1
  if item.kind == "D" and plots then pk, scr, grafBack = 1, "graf", "form"; return end
  scr = "res"
end

local function curStep(P)
  local s = P.series and P.series[P.cs or 1]
  return s and math.max(1, math.floor(#s.pts/60)) or 1
end
local function grafKey(k)
  local P = plots[pk]
  if P.draw then
    if k == "left" or k == "up" then pk = (pk - 2) % #plots + 1 elseif k == "right" or k == "down" then pk = pk % #plots + 1 end
    return
  end
  if P.move then
    if k == "left" then P.move(P, -1, 0) elseif k == "right" then P.move(P, 1, 0)
    elseif k == "up" then P.move(P, 0, 1) elseif k == "down" then P.move(P, 0, -1) end
    return
  end
  P.cs, P.ci = P.cs or 1, P.ci or 1
  local s = P.series[P.cs]
  if k == "left" then P.ci = math.max(1, P.ci - curStep(P))
  elseif k == "right" then P.ci = math.min(#s.pts, P.ci + curStep(P))
  elseif k == "up" or k == "down" then
    P.cs = (k == "up") and ((P.cs - 2) % #P.series + 1) or (P.cs % #P.series + 1)
    P.ci = math.min(P.ci, #P.series[P.cs].pts)
  end
end

function on.arrowUp()
  if scr == "main" then ci = (ci - 2) % #CATS + 1
  elseif scr == "sub" then ii = (ii - 2) % #CATS[ci].items + 1
  elseif scr == "form" then if #item.fields > 0 then fi = (fi - 2) % #item.fields + 1 end
  elseif scr == "res" then rtop = math.max(1, rtop - 1)
  elseif scr == "graf" then grafKey("up") end
  inv()
end
function on.arrowDown()
  if scr == "main" then ci = ci % #CATS + 1
  elseif scr == "sub" then ii = ii % #CATS[ci].items + 1
  elseif scr == "form" then if #item.fields > 0 then fi = fi % #item.fields + 1 end
  elseif scr == "res" then rtop = rtop + 1
  elseif scr == "graf" then grafKey("down") end
  inv()
end
function on.arrowLeft() if scr == "graf" then grafKey("left"); inv() end end
function on.arrowRight() if scr == "graf" then grafKey("right"); inv() end end
function on.tabKey()
  if scr == "graf" then pk = pk % #plots + 1; inv() else on.arrowDown() end
end
function on.backtabKey()
  if scr == "graf" then pk = (pk - 2) % #plots + 1; inv() else on.arrowUp() end
end
function on.enterKey()
  if scr == "main" then scr = "sub"; ii = 1; top = 1
  elseif scr == "sub" then openItem()
  elseif scr == "form" then compute()
  elseif scr == "res" then
    if plots then scr, pk, grafBack = "graf", 1, "res" else scr = "form"; top = 1 end
  elseif scr == "graf" then pk = pk % #plots + 1 end
  inv()
end
on.returnKey = on.enterKey
function on.escapeKey()
  if scr == "sub" then scr = "main"; top = 1
  elseif scr == "form" then scr = "sub"; top = 1
  elseif scr == "res" then scr = "form"; top = 1
  elseif scr == "graf" then scr = grafBack; top = 1 end
  inv()
end
function on.backspaceKey()
  if scr == "form" and item.fields[fi] then
    local id = item.fields[fi][1]
    local s = vals[id] or ""
    local p = #s
    while p > 0 and s:byte(p) >= 128 and s:byte(p) < 192 do p = p - 1 end
    vals[id] = s:sub(1, p - 1)
    inv()
  elseif scr ~= "main" then on.escapeKey() end
end
function on.clearKey()
  if scr == "form" and item.fields[fi] then vals[item.fields[fi][1]] = ""; inv() end
end
function on.charIn(ch)
  if scr == "form" and item.fields[fi] then
    local id = item.fields[fi][1]
    vals[id] = (vals[id] or "") .. ch
  elseif scr == "graf" then
    local P = plots[pk]
    if P.move then
      if ch == "*" then P.step = P.step*2 elseif ch == "/" then P.step = P.step/2 end
    end
  elseif scr == "main" or scr == "sub" then
    local n = tonumber(ch)
    if n and n >= 1 then
      if scr == "main" and n <= #CATS then ci = n; scr = "sub"; ii = 1; top = 1
      elseif scr == "sub" and n <= #CATS[ci].items then ii = n; openItem() end
    end
  end
  inv()
end
