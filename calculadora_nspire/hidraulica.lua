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
  local n, h = 800, nil
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
--  Resolucion de un ejercicio
-- =====================================================================
local function resolver(item, vals)
  local v, blanks = {}, {}
  for _, fd in ipairs(item.fields) do
    local s = vals[fd[1]] or ""
    if s:match("^%s*$") and fd.opt then s = fd.opt end
    if s:match("^%s*$") then blanks[#blanks + 1] = fd
    elseif fd.list then v[fd[1]] = evallist(s)
    else v[fd[1]] = evalstr(s) end
  end
  if item.kind == "P" then
    if #blanks > 0 then error("Falta dato: " .. blanks[1][2]) end
    return item.calc(v)
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
  return out
end

-- exportar para pruebas fuera de la calculadora
HID = {CATS = CATS, resolver = resolver, W = W, Wh = Wh, Fk = Fk, evalstr = evalstr}

-- =====================================================================
--  Interfaz grafica (TI-Nspire)
-- =====================================================================
local scr = "main"      -- main | sub | form | res
local ci, ii, fi = 1, 1, 1
local top = 1
local item
local vals = {}
local reslines, rtop = {}, 1
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
  gc:fillRect(0, H_() - 16, W_(), 16)
  gc:setColorRGB(60, 60, 60)
  gc:setFont("sansserif", "r", 8)
  gc:drawString(text, 3, H_() - 15, "top")
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

function on.paint(gc)
  gc:setFont("sansserif", "r", 9)
  if scr == "main" then
    header(gc, "CIV-346 Hidraulica Aplicada")
    local l = {}
    for k, c in ipairs(CATS) do l[k] = k .. ". " .. c.name end
    drawList(gc, l, ci)
    footer(gc, "enter: abrir   esc: -   " .. #CATS .. " clases")
  elseif scr == "sub" then
    header(gc, CATS[ci].name)
    local l = {}
    for k, it in ipairs(CATS[ci].items) do l[k] = (it.kind == "E" and "= " or "> ") .. it.name end
    drawList(gc, l, ii)
    footer(gc, "= ecuacion (despeja)  > procedimiento  esc: volver")
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
    gc:drawString(item.fields[fi][3], 4, H_() - 30, "top")
    footer(gc, item.kind == "E" and "Deje vacio el dato a calcular. enter: resolver" or "enter: calcular   esc: volver")
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
    footer(gc, "flechas: desplazar   esc/enter: volver")
  end
end

local function openItem()
  item = CATS[ci].items[ii]
  vals = {}
  fi, top = 1, 1
  scr = "form"
end
local function compute()
  local okc, r = pcall(resolver, item, vals)
  if okc then reslines = r else reslines = {"ERROR:", (tostring(r):gsub("^.-:%d+: ", ""))} end
  rtop = 1
  scr = "res"
end

function on.arrowUp()
  if scr == "main" then ci = (ci - 2) % #CATS + 1
  elseif scr == "sub" then ii = (ii - 2) % #CATS[ci].items + 1
  elseif scr == "form" then fi = (fi - 2) % #item.fields + 1
  elseif scr == "res" then rtop = math.max(1, rtop - 1) end
  inv()
end
function on.arrowDown()
  if scr == "main" then ci = ci % #CATS + 1
  elseif scr == "sub" then ii = ii % #CATS[ci].items + 1
  elseif scr == "form" then fi = fi % #item.fields + 1
  elseif scr == "res" then rtop = rtop + 1 end
  inv()
end
on.tabKey = on.arrowDown
on.backtabKey = on.arrowUp
function on.enterKey()
  if scr == "main" then scr = "sub"; ii = 1; top = 1
  elseif scr == "sub" then openItem()
  elseif scr == "form" then compute()
  elseif scr == "res" then scr = "form"; top = 1 end
  inv()
end
on.returnKey = on.enterKey
function on.escapeKey()
  if scr == "sub" then scr = "main"; top = 1
  elseif scr == "form" then scr = "sub"; top = 1
  elseif scr == "res" then scr = "form"; top = 1 end
  inv()
end
function on.backspaceKey()
  if scr == "form" then
    local id = item.fields[fi][1]
    local s = vals[id] or ""
    -- borrar un caracter UTF-8 completo
    local p = #s
    while p > 0 and s:byte(p) >= 128 and s:byte(p) < 192 do p = p - 1 end
    vals[id] = s:sub(1, p - 1)
    inv()
  elseif scr ~= "main" then on.escapeKey() end
end
function on.clearKey()
  if scr == "form" then vals[item.fields[fi][1]] = ""; inv() end
end
function on.charIn(ch)
  if scr == "form" then
    local id = item.fields[fi][1]
    vals[id] = (vals[id] or "") .. ch
  elseif scr == "main" or scr == "sub" then
    local n = tonumber(ch)
    if n and n >= 1 then
      if scr == "main" and n <= #CATS then ci = n; scr = "sub"; ii = 1; top = 1
      elseif scr == "sub" and n <= #CATS[ci].items then ii = n; openItem() end
    end
  end
  inv()
end
