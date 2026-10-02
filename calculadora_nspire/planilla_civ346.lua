-- =====================================================================
--  PLANILLA CIV-346 (tipo Excel) para TI-Nspire CX II CAS
--  Celdas A1..P99 con numeros, texto ('texto) o formulas (=...).
--  Formulas: operadores + - * / ^, referencias A1, $A$1, rangos A1:A9,
--  funciones generales (SUM, AVG, SLOPE, LN, IF...) y todas las
--  funciones hidraulicas del curso (THIEM, THEIS, JACOB, W, YN, ...).
--  [menu]: plantillas, copiar/pegar, rellenar abajo, serie, buscar
--  objetivo, graficar con regresion, exportar a lista, ayuda.
--  Plantillas "CERT ..." resuelven el certamen de practica completo.
--  (si la tecla menu no responde, escriba ? sobre una celda)
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
    gc:setColorRGB(c[1], c[2], c[3])
    local prev
    for _, p in ipairs(s.pts) do
      local xv, yv = tx(P, p[1]), ty(P, p[2])
      local cur
      if xv and yv and ok(xv) and ok(yv) then cur = {X(xv), Y(yv)} end
      if cur and s.style ~= "pts" and prev and inside(prev[1], prev[2]) and inside(cur[1], cur[2]) then
        gc:drawLine(prev[1], prev[2], cur[1], cur[2])
      end
      if cur and s.style ~= "line" and inside(cur[1], cur[2]) then gc:fillRect(cur[1] - 1, cur[2] - 1, 3, 3) end
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

-- =====================================================================
--  PLANILLA (tipo Excel)
-- =====================================================================
local NCOL, NROW = 16, 99
local COLS = "ABCDEFGHIJKLMNOP"
local function colName(c) return COLS:sub(c, c) end
local function colIdx(ch) return COLS:find(ch:upper(), 1, true) end
local function key(c, r) return colName(c) .. r end

local cells = {}       -- contenido crudo por celda "A1" -> "=..."
local cache = {}       -- valores calculados
local compiled = {}    -- funciones compiladas por texto de formula
local evaluating = {}

-- -------------------- Biblioteca de funciones -----------------------
local function flat(...)
  local t = {}
  for _, a in ipairs({...}) do
    if type(a) == "table" then for _, x in ipairs(a) do if type(x) == "number" then t[#t+1] = x end end
    elseif type(a) == "number" then t[#t+1] = a end
  end
  return t
end
local function pairsXY(ys, xs)
  local X, Y = {}, {}
  for i = 1, math.min(#xs, #ys) do if type(xs[i]) == "number" and type(ys[i]) == "number" then X[#X+1] = xs[i]; Y[#Y+1] = ys[i] end end
  return X, Y
end
local function secc(tp, b, k, D, y) return seccion(math.floor(tp + 0.5), b or 0, k or 0, D or 0, y) end
local function yn(tp, Q, n, S, b, k, D)
  local lim = (math.floor(tp + 0.5) == 4) and D*0.9999 or 1e4
  return roots(function(y) local A, P_ = secc(tp, b, k, D, y); return Q*n/sqrt(S) - A^(5/3)/P_^(2/3) end, 1e-6, lim, 1500, true)[1] or 0/0
end
local function yc(tp, Q, b, k, D)
  local lim = (math.floor(tp + 0.5) == 4) and D*0.9999 or 1e4
  return roots(function(y) local A, _, B = secc(tp, b, k, D, y); return Q*Q/G - A^3/B end, 1e-6, lim, 1500, true)[1] or 0/0
end
local function dprime(d, D, r0)
  local chi = pi*r0
  if d/D < 0.25 then return d/(1 + 8/pi*d/D*ln(d/chi)) end
  return d*pi/(2*(ln(D/chi) + 0.18))
end

-- tuberia circular parcialmente llena (r = y/D)
local function qparc(D, n, S, r) local A, P_ = seccion(4, 0, 0, D, r*D); return sqrt(S)/n*A*(A/P_)^(2/3) end
local DCOMS = {0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0}
local function truthy(c) return c and c ~= 0 end

local FUNC = {
  -- generales
  PI = pi, GRAV = G, E = exp(1),
  SUM = function(...) local s = 0; for _, x in ipairs(flat(...)) do s = s + x end; return s end,
  PROD = function(...) local s = 1; for _, x in ipairs(flat(...)) do s = s*x end; return s end,
  AVG = function(...) local t = flat(...); local s = 0; for _, x in ipairs(t) do s = s + x end; return s/#t end,
  MIN = function(...) local t = flat(...); local m = t[1]; for _, x in ipairs(t) do if x < m then m = x end end; return m end,
  MAX = function(...) local t = flat(...); local m = t[1]; for _, x in ipairs(t) do if x > m then m = x end end; return m end,
  COUNT = function(...) return #flat(...) end,
  SLOPE = function(ys, xs) local X, Y = pairsXY(ys, xs); local a = regresion(X, Y); return a end,
  INTERCEPT = function(ys, xs) local X, Y = pairsXY(ys, xs); local _, b = regresion(X, Y); return b end,
  RSQ = function(ys, xs) local X, Y = pairsXY(ys, xs); local a, b = regresion(X, Y)
    local m = 0; for i = 1, #Y do m = m + Y[i] end; m = m/#Y
    local st, sr = 0, 0; for i = 1, #Y do st = st + (Y[i] - m)^2; sr = sr + (Y[i] - a*X[i] - b)^2 end
    return 1 - sr/st end,
  LN = ln, LOG = function(x) return ln(x)/ln(10) end, LOG10 = function(x) return ln(x)/ln(10) end, EXP = exp, SQRT = sqrt, ABS = abs,
  SIN = math.sin, COS = math.cos, TAN = math.tan, ASIN = math.asin, ACOS = math.acos, ATAN = math.atan,
  SINH = function(x) return (exp(x) - exp(-x))/2 end, COSH = function(x) return (exp(x) + exp(-x))/2 end,
  ASINH = function(x) return ln(x + sqrt(x*x + 1)) end,
  DEG = function(x) return x*180/pi end, RAD = function(x) return x*pi/180 end,
  ROUND = function(x, n) local m = 10^(n or 0); return math.floor(x*m + 0.5)/m end,
  IF = function(c, a, b) if truthy(c) then return a else return b end end,
  AND = function(...) for _, c in ipairs({...}) do if not truthy(c) then return false end end; return true end,
  OR = function(...) for _, c in ipairs({...}) do if truthy(c) then return true end end; return false end,
  NOT = function(c) return not truthy(c) end,
  CEIL = function(x, p) p = p or 1; return math.ceil(x/p - 1e-9)*p end,
  FLOOR = function(x, p) p = p or 1; return math.floor(x/p + 1e-9)*p end,
  -- clase 01
  DARCY = function(K, A, i) return K*A*i end,
  HAZEN = function(C, d10) return C*d10^2 end,
  SLICHTER = function(n) return slichter(n) end,
  KOZENY = function(n, d50, rho, mu) rho = rho or 1000; mu = mu or 1e-3; return rho*G/mu*n^3/(1 - n)^2*d50^2/180 end,
  KCTE = function(L, Q, dh, D) return 4*L*Q/(pi*dh*D^2) end,
  KVAR = function(dt, dc, L, t, h0, h) return dt^2*L/(dc^2*t)*ln(h0/h) end,
  -- clase 02-03
  THIEM = function(Q, T, R, r) if r >= R then return 0 end return Q/(2*pi*T)*ln(R/r) end,
  DUPUIT = function(H, Q, K, R, r) if r >= R then return H end return sqrt(H*H - Q/(pi*K)*ln(R/r)) end,
  JACOBCORR = function(dh, H) return dh - dh^2/(2*H) end,
  ESFERICO = function(Q, K, r, R) return Q/(4*pi*K)*(1/r - 1/R) end,
  SEMIESF = function(Q, K, r, R) return Q/(2*pi*K)*(1/r - 1/R) end,
  SICHARDT = function(dh, K) return 3000*dh*sqrt(K) end,
  DIST = function(x1, y1, x2, y2) return sqrt((x1 - x2)^2 + (y1 - y2)^2) end,
  -- clase 04-05
  GHYBEN = function(hd, gs) gs = gs or 1.025; return hd/(gs - 1) end,
  W = function(u) return W(u) end,
  WH = function(u, rb) return Wh(u, rb) end,
  U = function(r, S, T, t) return r*r*S/(4*T*t) end,
  THEIS = function(Q, T, S, r, t) if t <= 0 then return 0 end return Q/(4*pi*T)*W(r*r*S/(4*T*t)) end,
  JACOB = function(Q, T, S, r, t) if t <= 0 then return 0 end local a = CJ*T*t/(r*r*S); if a <= 1 then return 0 end return Q/(4*pi*T)*ln(a) end,
  RT = function(T, t, S) return sqrt(CJ*T*t/S) end,
  TINF = function(r, S, T) return r*r*S/(CJ*T) end,
  HANTUSH = function(Q, T, S, r, t, B) return Q/(4*pi*T)*Wh(r*r*S/(4*T*t), r/B) end,
  -- clase 07-08
  DONNAN = function(h0, f, D, K) return sqrt(h0*h0 + f*D*D/(4*K)) end,
  HOOGH = function(h0, d, f, D, K) return -d + sqrt((d + h0)^2 + f*D*D/(4*K)) end,
  DAGAN = function(r0, d) if r0/d >= 0.3 then return 0 end return -ln(2*(math.cosh and math.cosh(pi*r0/d) or (exp(pi*r0/d) + exp(-pi*r0/d))/2) - 2)/pi end,
  KIRKF = function(x, y) return Fk(x, y) end,
  DPRIME = dprime,
  GLOVERY = function(y0, K, d, S, D, t) return 1.16*y0*exp(-pi^2*K*d/(S*D*D)*t) end,
  GLOVERD = function(K, d, t, S, y0, y) return pi*sqrt(K*d*t/(S*ln(1.16*y0/y))) end,
  -- clase 09-12
  RACIONAL = function(C, i, A) return C*i*A/3.6 end,
  GRUNSKY = function(i24, t) return i24*sqrt(24/t) end,
  BELL = function(t, T, P1) return (0.54*t^0.25 - 0.5)*(0.21*ln(T) + 0.52)*P1 end,
  RIESGO = function(T, n) return 1 - (1 - 1/T)^n end,
  TRIESGO = function(R, n) return 1/(1 - (1 - R)^(1/n)) end,
  BELLI = function(t, T, P1) return 60*(0.54*t^0.25 - 0.5)*(0.21*ln(T) + 0.52)*P1/t end,
  IDF = function(a, b, c, t) return a/(t + b)^c end,
  HCRECT = function(Q, b) return (Q*Q/(G*b*b))^(1/3) end,
  -- pruebas de bombeo
  CJT = function(Q, a) return Q/(4*pi*a) end,
  CJS = function(T, t0, r) return CJ*T*t0/(r*r) end,
  RIMAGEN = function(r1, t1, t2) return r1*sqrt(t2/t1) end,
  EFIC = function(B, C, Q) return B*Q/(B*Q + C*Q*Q) end,
  TCCALIF = function(L, H) return 57*(L^3/H)^0.385 end,
  TCGIAND = function(A, L, Hm) return 60*(4*sqrt(A) + 1.5*L)/(0.8*sqrt(Hm)) end,
  TCSCS = function(L, CN, S) return 3.42*L^0.8*(1000/CN - 9)^0.7/S^0.5 end,
  TCESP = function(L, S) return 18*L^0.76/S^0.19 end,
  AREA = function(tp, b, k, D, y) local A = secc(tp, b, k, D, y); return A end,
  PERIM = function(tp, b, k, D, y) local _, P_ = secc(tp, b, k, D, y); return P_ end,
  ANCHO = function(tp, b, k, D, y) local _, _, B = secc(tp, b, k, D, y); return B end,
  YN = yn, YC = yc,
  MANNINGQ = function(tp, b, k, D, y, n, S) local A, P_ = secc(tp, b, k, D, y); return sqrt(S)/n*A*(A/P_)^(2/3) end,
  FROUDE = function(tp, b, k, D, y, Q) local A, _, B = secc(tp, b, k, D, y); return Q/A/sqrt(G*A/B) end,
  DRECT = function(Q, b) return 1.5/G^(1/3)*(Q/b)^(2/3) end,
  SCRECT = function(H, n, b) return 2/3*G*H*n^2*((b + 4/3*H)/(2/3*H*b))^(4/3) end,
  QCIRC = function(D) return 1.425*D^2.5 end,
  SCCIRC = function(n, D) return 31.14*n^2/D^(1/3) end,
  SNCIRC = function(Q, n, D) return 10.3*Q^2*n^2/D^(16/3) end,
  SNRECT = function(Q, n, b, D) local A, R = b*D, b*D/(2*b + 2*D); return (Q*n/(A*R^(2/3)))^2 end,
  HCASO1 = function(D, Q, b, cc) cc = cc or 0.611; return cc*D + (Q/(b*cc*D))^2/(2*G) end,
  HCASO2 = function(D, ke, Q, b) return D + (ke + 1)*Q^2/((D*b)^2*2*G) end,
  YARNELL = function(K, Q, b, h1, n, D) local Fr2 = (Q/(b*h1))^2/(G*h1); local s = n*D/b
    return h1*K*Fr2*(K + 5*Fr2 - 0.6)*(s + 15*s^4) end,
  CUNETA = function(S, n, b, i) return sqrt(S)/n*b^2*i/2*(b*i/(2*(i + 1)))^(2/3) end,
  CUNETAC = function(S, n, b, i1, i2, w) return 0.315*sqrt(S)/n*b^(8/3)*(i2 + (i1 - i2)*(w/b)^2)^(5/3) end,
  QVERT = function(Le, we, h) return 1.66*(Le + 2*we)*h^1.5 end,
  QORIF = function(Ae, h) return 2.66*Ae*h^0.5 end,
  QLAT = function(L, a, h) if h < a then return 1.27*L*h^1.5 end return 2.66*L*a*h^0.5 end,
  VAUTOL = function(D, S, n) return 0.397*D^(2/3)*sqrt(S)/n end,
  QPARC = qparc,
  INDEXY = function(xs, ys, x) local r = 0/0; for i = 1, math.min(#xs, #ys) do if xs[i] == x then r = ys[i] end end; return r end,
  DCOM = function(Q, n, S, r) r = r or 0.8; for _, D in ipairs(DCOMS) do if qparc(D, n, S, r) >= Q then return D end end; return 0/0 end,
  YND = function(Q, n, S, D) return yn(4, Q, n, S, 0, 0, D)/D end,
  VCIRC = function(Q, n, S, D) local y = yn(4, Q, n, S, 0, 0, D); local A = seccion(4, 0, 0, D, y); return Q/A end,
  HLIMSUM = function(Ae, Le, we) return 1.6*Ae/(Le + 2*we) end,
  QSUMF = function(Le, we, h) local Ae = Le*we; if h < 1.6*Ae/(Le + 2*we) then return 1.66*(Le + 2*we)*h^1.5 end return 2.66*Ae*h^0.5 end,
  SEPSUM = function(eta, Q, C, i, bc) return 3600*eta*Q/(C*i/1000*bc) end,
}
local FHELP = {
  {"SUM(rango)", "suma"}, {"AVG(rango)", "promedio"}, {"MIN/MAX(rango)", ""}, {"COUNT(rango)", "n datos"},
  {"SLOPE(Y,X)", "pendiente regresion"}, {"INTERCEPT(Y,X)", "intercepto"}, {"RSQ(Y,X)", "R^2"},
  {"LN, LOG, EXP, SQRT, ABS", ""}, {"SIN/COS/TAN (rad), DEG, RAD", ""}, {"IF(c,a,b)", "condicional"}, {"ROUND(x,n)", ""},
  {"PI, GRAV (9.81)", "constantes"},
  {"DARCY(K,A,i)", "Q"}, {"HAZEN(C,d10)", "K cm/s"}, {"SLICHTER(n)", "1/K1"}, {"KOZENY(n,d50)", "K"},
  {"KCTE(L,Q,dh,D)", "permeametro cte"}, {"KVAR(dt,dc,L,t,h0,h)", "permeametro var"},
  {"THIEM(Q,T,R,r)", "descenso conf."}, {"DUPUIT(H,Q,K,R,r)", "h libre"}, {"JACOBCORR(dh,H)", "dh'"},
  {"ESFERICO(Q,K,r,R)", ""}, {"SEMIESF(Q,K,r,R)", ""}, {"SICHARDT(dh,K)", "R"}, {"DIST(x1,y1,x2,y2)", "distancia"},
  {"GHYBEN(hd,gs)", "hs"}, {"W(u)", "funcion de pozo"}, {"WH(u,r/B)", "Hantush"}, {"U(r,S,T,t)", "u"},
  {"THEIS(Q,T,S,r,t)", "descenso"}, {"JACOB(Q,T,S,r,t)", "descenso"}, {"RT(T,t,S)", "R(t)"}, {"TINF(r,S,T)", "t influencia"},
  {"HANTUSH(Q,T,S,r,t,B)", "descenso"},
  {"DONNAN(h0,f,D,K)", "H"}, {"HOOGH(h0,d,f,D,K)", "H"}, {"DAGAN(r0,d)", "alfa"}, {"KIRKF(2r0/D,d/D)", "F"},
  {"DPRIME(d,D,r0)", "d'"}, {"GLOVERY(y0,K,d,S,D,t)", "y"}, {"GLOVERD(K,d,t,S,y0,y)", "D"},
  {"RACIONAL(C,i,A)", "Q m3/s"}, {"GRUNSKY(i24,t)", ""}, {"BELL(t,T,P1)", "P"}, {"RIESGO(T,n)", ""},
  {"TCCALIF(L,H) TCGIAND(A,L,Hm)", ""}, {"TCSCS(L,CN,S) TCESP(L,S)", ""},
  {"AREA/PERIM/ANCHO(tp,b,k,D,y)", "tp 1rect 2trap 3tri 4circ"}, {"YN(tp,Q,n,S,b,k,D)", "altura normal"},
  {"YC(tp,Q,b,k,D)", "altura critica"}, {"MANNINGQ(tp,b,k,D,y,n,S)", "Q"}, {"FROUDE(tp,b,k,D,y,Q)", ""},
  {"DRECT(Q,b)", "D control entrada"}, {"SCRECT(H,n,b)", "Sc"}, {"QCIRC(D)", ""}, {"SCCIRC(n,D)", ""},
  {"SNCIRC(Q,n,D)", ""}, {"SNRECT(Q,n,b,D)", ""}, {"HCASO1(D,Q,b,cc)", "H'"}, {"HCASO2(D,ke,Q,b)", "H'"},
  {"YARNELL(K,Q,b,h1,n,D)", "dh"}, {"CUNETA(S,n,b,i)", "Q"}, {"CUNETAC(S,n,b,i1,i2,w)", "Q"},
  {"QVERT(Le,we,h) QORIF(Ae,h)", "sumidero fondo"}, {"QLAT(L,a,h)", "sumidero lateral"}, {"VAUTOL(D,S,n)", "Vc"},
  {"AND/OR/NOT(c..)", "logicos"}, {"CEIL/FLOOR(x,paso)", "redondeo"}, {"TRIESGO(R,n)", "T para riesgo R"},
  {"BELLI(t,T,P1)", "i mm/h (t min)"}, {"IDF(a,b,c,t)", "a/(t+b)^c"}, {"HCRECT(Q,b)", "hc rect"},
  {"CJT(Q,a)", "T con pend. ln"}, {"CJS(T,t0,r)", "S"}, {"RIMAGEN(r1,t1,t2)", "r pozo imagen"}, {"EFIC(B,C,Q)", "eficiencia"},
  {"QPARC(D,n,S,y/D)", "Q tubo parcial"}, {"DCOM(Q,n,S,y/D)", "D comercial"}, {"YND(Q,n,S,D)", "yn/D"},
  {"VCIRC(Q,n,S,D)", "v normal tubo"}, {"HLIMSUM(Ae,Le,we)", "h vert/orif"}, {"QSUMF(Le,we,h)", "Q sumidero fondo"},
  {"SEPSUM(eta,Q,C,i,bc)", "separacion (i mm/h)"}, {"INDEXY(X,Y,x)", "y donde X=x"},
}

-- -------------------- Evaluacion de celdas --------------------------
local REF = "%f[%w%$](%$?)([A-Pa-p])(%$?)(%d+)%f[^%w]"
local RANGE = "%f[%w%$]%$?([A-Pa-p])%$?(%d+):%$?([A-Pa-p])%$?(%d+)%f[^%w]"

local evalCell
local function rng(c1, r1, c2, r2)
  local t = {}
  for c = math.min(c1, c2), math.max(c1, c2) do
    for r = math.min(r1, r2), math.max(r1, r2) do
      local v = evalCell(key(c, r))
      if type(v) == "number" then t[#t+1] = v else t[#t+1] = false end
    end
  end
  return t
end
local function cellv(c, r)
  local v = evalCell(key(c, r))
  if v == nil then return 0 end
  if v == "" then error("#BLANCO") end
  if type(v) ~= "number" then error("#VALOR " .. key(c, r)) end
  return v
end
local FENV = setmetatable({CELL = cellv, RNG = rng}, {__index = FUNC})

local function compile(raw)
  if compiled[raw] then return compiled[raw] end
  local e = raw:sub(2):upper():gsub(";", ","):gsub("<>", "~="):gsub("%%", "/100")
  e = e:gsub("([^<>~=])=([^=])", "%1==%2")
  e = e:gsub(RANGE, function(a, b, c, d) return "RNG(" .. colIdx(a) .. "," .. b .. "," .. colIdx(c) .. "," .. d .. ")" end)
  e = e:gsub(REF, function(_, c, _, r) return "CELL(" .. colIdx(c) .. "," .. r .. ")" end)
  local f, err
  if setfenv then f, err = loadstring("return " .. e); if f then setfenv(f, FENV) end
  else f, err = load("return " .. e, "f", "t", FENV) end
  if not f then f = function() error("#SINTAXIS") end end
  compiled[raw] = f
  return f
end

function evalCell(k)
  if cache[k] ~= nil then return cache[k] end
  local raw = cells[k]
  if raw == nil or raw == "" then return nil end
  if evaluating[k] then error("#CIRCULAR") end
  local v
  if raw:sub(1, 1) == "=" then
    evaluating[k] = true
    local okc, res = pcall(compile(raw))
    evaluating[k] = nil
    if okc and type(res) == "number" and ok(res) then v = res
    elseif okc and type(res) == "number" then v = "#NUM"
    elseif okc and type(res) == "string" then v = res
    elseif okc then v = "#VALOR"
    else
      local m = tostring(res):match("#%u+") or "#ERR"
      if m == "#BLANCO" then m = "" end
      if m == "#CIRCULAR" then evaluating = {} ; error(res) end
      v = m
    end
  else
    local okn, n = pcall(evalstr, raw)
    if okn and n then v = n else v = raw:gsub("^'", "") end
  end
  cache[k] = v
  return v
end
local function valueOf(k)
  local okc, v = pcall(evalCell, k)
  if okc then return v end
  evaluating = {}
  cache[k] = "#CIRC"
  return "#CIRC"
end
local function recalc() cache = {}; evaluating = {} end
local function setCell(k, raw)
  if raw == nil or raw == "" then cells[k] = nil else cells[k] = raw end
  recalc()
  if document and document.markChanged then pcall(document.markChanged) end
end

-- desplaza referencias relativas (copiar / rellenar)
local function shiftRefs(raw, dr, dc)
  if raw:sub(1, 1) ~= "=" then return raw end
  return (raw:gsub(REF, function(dcol, c, drow, r)
    local ci_, ri_ = colIdx(c), tonumber(r)
    if dcol == "" then ci_ = ci_ + dc end
    if drow == "" then ri_ = ri_ + dr end
    if ci_ < 1 or ci_ > NCOL or ri_ < 1 then return "#REF" end
    return dcol .. colName(ci_) .. drow .. ri_
  end))
end

local function parseRef(s)
  s = (s or ""):upper():gsub("%$", ""):gsub("%s", "")
  local c, r = s:match("^([A-P])(%d+)$")
  if c then return colIdx(c), tonumber(r) end
end
local function parseRange(s)
  s = (s or ""):upper():gsub("%$", ""):gsub("%s", "")
  local c1, r1, c2, r2 = s:match("^([A-P])(%d+):([A-P])(%d+)$")
  if c1 then return colIdx(c1), tonumber(r1), colIdx(c2), tonumber(r2) end
end
local function rangeVals(s)
  local c1, r1, c2, r2 = parseRange(s)
  if not c1 then error("rango invalido: " .. tostring(s)) end
  local t = {}
  for c = c1, c2 do for r = r1, r2 do t[#t+1] = valueOf(key(c, r)) end end
  return t
end

-- -------------------- Plantillas ------------------------------------
-- cada plantilla: lista de {celda, contenido}
local TPL = {}
local function T(name, desc, list) TPL[#TPL+1] = {name = name, desc = desc, list = list} end
local function fillCol(list, col, r1, r2, raw)
  for r = r1, r2 do list[#list+1] = {col .. r, shiftRefs(raw, r - r1, 0)} end
end

-- ---- Plantillas del certamen (datos precargados; cambie la columna de datos) ----
local function put(L, col, r1, vals) for i, v in ipairs(vals) do L[#L+1] = {col .. (r1 + i - 1), tostring(v)} end end
do -- P1 atravieso completo
  local L = {{"A1","'DATOS"},{"C1","'unid."},
    {"A2","'A cuenca"},{"B2","0.9"},{"C2","'km2"},{"A3","'L cauce"},{"B3","1.6"},{"C3","'km"},{"A4","'H desniv"},{"B4","120"},{"C4","'m"},
    {"A5","'C10"},{"B5","=0.24+0.07+0.07+0.07"},{"C5","'suma tabla"},{"A6","'T1 diseno"},{"B6","50"},{"C6","'anos"},
    {"A7","'fC T1"},{"B7","1.2"},{"A8","'T2 verif"},{"B8","100"},{"C8","'anos"},{"A9","'fC T2"},{"B9","1.25"},
    {"A10","'P1^10"},{"B10","28"},{"C10","'mm"},{"A11","'vida util"},{"B11","30"},{"C11","'anos"},
    {"A12","'b cajon"},{"B12","2"},{"C12","'m"},{"A13","'n"},{"B13","0.015"},{"A14","'L cajon"},{"B14","24"},{"C14","'m"},
    {"A15","'SD"},{"B15","0.012"},{"C15","'m/m"},{"A16","'ke"},{"B16","0.5"},{"A17","'rasante"},{"B17","3.6"},{"C17","'m"},
    {"A18","'revancha"},{"B18","0.3"},{"C18","'m"},{"A19","'paso D"},{"B19","0.1"},{"C19","'m"},
    {"D1","'RESULTADO"},{"F1","'unid."},
    {"D2","'a) tc Calif"},{"E2","=TCCALIF(B3,B4)"},{"F2","'min"},{"D3","'tc>=10?"},{"E3","=IF(E2>=10,1,0)"},{"F3","'1=si"},
    {"D4","'b) riesgo"},{"E4","=RIESGO(B6,B11)"},{"F4","'T1 en n"},{"D5","'S util"},{"E5","=B12*E15"},{"F5","'m2 (>1.75?)"},
    {"D6","'c) C T1"},{"E6","=B5*B7"},{"D7","'P T1"},{"E7","=BELL(E2,B6,B10)"},{"F7","'mm"},
    {"D8","'i T1"},{"E8","=60*E7/E2"},{"F8","'mm/h"},{"D9","'Q T1"},{"E9","=RACIONAL(E6,E8,B2)"},{"F9","'m3/s"},
    {"D10","'C T2"},{"E10","=B5*B9"},{"D11","'i T2"},{"E11","=BELLI(E2,B8,B10)"},{"F11","'mm/h"},
    {"D12","'Q T2"},{"E12","=RACIONAL(E10,E11,B2)"},{"F12","'m3/s"},
    {"D13","'d) hc"},{"E13","=HCRECT(E9,B12)"},{"F13","'m"},{"D14","'D=1.5hc"},{"E14","=1.5*E13"},{"F14","'m"},
    {"D15","'D adopt"},{"E15","=CEIL(E14,B19)"},{"F15","'m"},{"D16","'v=Q/bhc"},{"E16","=E9/(B12*E13)"},{"F16","'m/s"},
    {"D17","'v<5 ?"},{"E17","=IF(E16<5,1,0)"},{"D18","'Sc"},{"E18","=SCRECT(E14,B13,B12)"},
    {"D19","'SD>=Sc?"},{"E19","=IF(B15>=E18,1,0)"},{"F19","'torrente"},
    {"D20","'e) Sn Q2"},{"E20","=SNRECT(E12,B13,B12,E15)"},{"D21","'caso"},{"E21","=IF(B15>E20,1,IF(B15<E20,3,2))"},
    {"D22","'H' caso1"},{"E22","=HCASO1(E15,E12,B12)"},{"F22","'m"},{"D23","'H' caso2"},{"E23","=HCASO2(E15,B16,E12,B12)"},{"F23","'m"},
    {"D24","'v llena"},{"E24","=E12/(B12*E15)"},{"D25","'J"},{"E25","=E24^2*B13^2/(B12*E15/(2*B12+2*E15))^(4/3)"},
    {"D26","'H' caso3"},{"E26","=E15+(1+B16)*E24^2/(2*GRAV)+E25*B14-B15*B14"},{"F26","'m"},
    {"D27","'H' usar"},{"E27","=IF(E21=1,E22,IF(E21=2,E23,E26))"},{"F27","'m"},
    {"D28","'H max"},{"E28","=B17-B18"},{"F28","'rasante-rev"},{"D29","'H'<=Hmax?"},{"E29","=IF(E27<=E28,1,0)"},
    {"D30","'D+0.6"},{"E30","=E15+0.6"},{"D31","'H'<=D+0.6?"},{"E31","=IF(E27<=E30,1,0)"}}
  T("CERT P1 Atravieso (cajon)", "a)-e) completo; datos en B", L)
end
do -- P2a prueba de gasto variable
  local L = {{"A1","'Q L/s"},{"B1","'sw m"},{"C1","'Q m3/s"},{"D1","'sw/Q"},
    {"F1","'C pozo"},{"G1","=SLOPE(D2:D20,C2:C20)"},{"H1","'s2/m5"},{"F2","'B acuif"},{"G2","=INTERCEPT(D2:D20,C2:C20)"},{"H2","'s/m2"},
    {"F3","'R^2"},{"G3","=RSQ(D2:D20,C2:C20)"},{"F5","'Q prueba"},{"G5","25"},{"H5","'L/s"},
    {"F6","'BQ"},{"G6","=G2*G5/1000"},{"H6","'m"},{"F7","'CQ^2"},{"G7","=G1*(G5/1000)^2"},{"H7","'m"},
    {"F8","'sw"},{"G8","=G6+G7"},{"H8","'m"},{"F9","'eficiencia"},{"G9","=EFIC(G2,G1,G5/1000)"},
    {"F10","'0.8Qmax"},{"G10","=0.8*MAX(A2:A20)"},{"H10","'L/s"},{"F11","'Qp<=0.8Qm"},{"G11","=IF(G5<=G10*1.03,1,0)"}}
  put(L, "A", 2, {8, 16, 24, 32}); put(L, "B", 2, {1.91, 4.04, 6.41, 9.01})
  fillCol(L, "C", 2, 20, "=IF(A2>0,A2/1000,\"\")")
  fillCol(L, "D", 2, 20, "=IF(A2>0,B2/C2,\"\")")
  T("CERT P2a Gasto variable", "sw=BQ+CQ^2, eficiencia", L)
end
do -- P2b-c-e Cooper-Jacob con tramo elegido + borde
  local L = {{"A1","'t min"},{"B1","'s m"},{"C1","'ln t tramo"},{"D1","'s recta"},{"E1","'desvio"},{"F1","'ln t final"},
    {"H1","'Q m3/s"},{"I1","0.025"},{"H2","'r m"},{"I2","25"},{"H3","'m espesor"},{"I3","18"},{"H4","'t->s"},{"I4","60"},
    {"H5","'tramo ini"},{"I5","10"},{"H6","'tramo fin"},{"I6","120"},{"H7","'final desde"},{"I7","360"},
    {"H9","'a=ds/ln"},{"I9","=SLOPE(B2:B40,C2:C40)"},{"H10","'ds/ciclo"},{"I10","=2.303*I9"},{"J10","'m"},
    {"H11","'T"},{"I11","=CJT(I1,I9)"},{"J11","'m2/s"},{"H12","'K=T/m"},{"I12","=I11/I3"},{"J12","'m/s"},
    {"H13","'t0"},{"I13","=EXP(-INTERCEPT(B2:B40,C2:C40)/I9)"},{"J13","'min"},
    {"H14","'S"},{"I14","=CJS(I11,I13*I4,I2)"},{"H15","'u(t ini)"},{"I15","=U(I2,I14,I11,I5*I4)"},
    {"H16","'t(u=.01)"},{"I16","=I2^2*I14/(0.04*I11)/I4"},{"J16","'min"},{"H17","'R^2"},{"I17","=RSQ(B2:B40,C2:C40)"},
    {"H19","'a final"},{"I19","=SLOPE(B2:B40,F2:F40)"},{"H20","'a fin/a"},{"I20","=I19/I9"},
    {"H21","'borde"},{"I21","=IF(I20>1.2,\"IMPERM\",IF(I20<0.8,\"RECARGA\",\"NINGUNO\"))"},
    {"H22","'t2 ultimo"},{"I22","=MAX(A2:A40)"},{"J22","'min"},{"H23","'ds en t2"},{"I23","=ABS(INDEXY(A2:A40,B2:B40,I22)-(I9*LN(I22)+INTERCEPT(B2:B40,C2:C40)))"},
    {"H24","'t1 (s=ds)"},{"I24","=EXP((I23-INTERCEPT(B2:B40,C2:C40))/I9)"},{"J24","'min"},
    {"H25","'r imagen"},{"I25","=RIMAGEN(I2,I24,I22)"},{"J25","'m"},{"H26","'d borde"},{"I26","=(I2+I25)/2"},{"J26","'m"},
    {"H28","'R(t2)"},{"I28","=RT(I11,I22*I4,I14)"},{"J28","'m"},{"H29","'R>2d ?"},{"I29","=IF(I28>2*I26,1,0)"}}
  put(L, "A", 2, {1, 2, 3, 5, 8, 10, 15, 20, 30, 45, 60, 90, 120, 180, 240, 360, 480, 720, 1080, 1440})
  put(L, "B", 2, {0.65, 0.86, 0.99, 1.16, 1.31, 1.39, 1.52, 1.61, 1.75, 1.89, 1.99, 2.14, 2.26, 2.45, 2.59, 2.81, 2.98, 3.22, 3.46, 3.65})
  fillCol(L, "C", 2, 40, "=IF(AND(A2>=$I$5,A2<=$I$6),LN(A2),\"\")")
  fillCol(L, "D", 2, 40, "=IF(A2>0,$I$9*LN(A2)+INTERCEPT($B$2:$B$40,$C$2:$C$40),\"\")")
  fillCol(L, "E", 2, 40, "=IF(A2>0,B2-D2,\"\")")
  fillCol(L, "F", 2, 40, "=IF(A2>=$I$7,LN(A2),\"\")")
  T("CERT P2b Cooper-Jacob+borde", "elija tramo en I5:I6; t en min", L)
end
do -- P2d recuperacion con t'
  local L = {{"A1","'t' min"},{"B1","'s' m"},{"C1","'t/t'"},{"D1","'ln tramo"},
    {"F1","'tf min"},{"G1","1440"},{"F2","'Q m3/s"},{"G2","0.025"},{"F3","'t/t' min"},{"G3","20"},{"H3","'usar t/t'>=G3"},
    {"F5","'a=ds/ln"},{"G5","=SLOPE(B2:B30,D2:D30)"},{"F6","'ds/ciclo"},{"G6","=2.303*G5"},{"H6","'m"},
    {"F7","'T recup"},{"G7","=CJT(G2,G5)"},{"H7","'m2/s"},{"F8","'R^2"},{"G8","=RSQ(B2:B30,D2:D30)"},
    {"F10","'T bombeo"},{"G10","0.0057"},{"F11","'dif %"},{"G11","=100*(G7-G10)/G10"}}
  put(L, "A", 2, {1, 2, 5, 10, 20, 30, 60, 120, 240, 480})
  put(L, "B", 2, {3.00, 2.78, 2.49, 2.26, 2.04, 1.91, 1.68, 1.43, 1.15, 0.85})
  fillCol(L, "C", 2, 30, "=IF(A2>0,($G$1+A2)/A2,\"\")")
  fillCol(L, "D", 2, 30, "=IF(A2>0,IF(C2>=$G$3,LN(C2),\"\"),\"\")")
  T("CERT P2d Recuperacion (t')", "t' desde que se detiene", L)
end
do -- P3a-c cuneta, sumidero y separacion
  local L = {{"A1","'DATOS"},{"A2","'S long"},{"B2","0.025"},{"A3","'n cuneta"},{"B3","0.016"},{"A4","'b max"},{"B4","1.5"},{"C4","'m"},
    {"A5","'i bombeo"},{"B5","0.03"},{"A6","'C"},{"B6","0.75"},{"A7","'I mm/h"},{"B7","42"},
    {"A8","'w sumid"},{"B8","0.6"},{"C8","'m"},{"A9","'L sumid"},{"B9","1"},{"C9","'m"},{"A10","'e barra"},{"B10","0.02"},{"C10","'m"},
    {"A11","'nL"},{"B11","7"},{"A12","'nT"},{"B12","3"},{"A13","'eta"},{"B13","0.75"},{"A14","'F seg"},{"B14","0.8"},{"C14","'tabla MDU"},
    {"A15","'calzada e"},{"B15","8"},{"C15","'m"},{"A16","'doble 1/0"},{"B16","1"},{"A17","'wc lateral"},{"B17","15"},{"C17","'m"},
    {"D1","'RESULTADO"},{"D2","'a) Qmax"},{"E2","=CUNETA(B2,B3,B4,B5)"},{"F2","'m3/s"},
    {"D3","'y=b i"},{"E3","=B4*B5"},{"F3","'m"},{"D4","'y<=0.15?"},{"E4","=IF(E3<=0.15,1,0)"},
    {"D5","'b) we"},{"E5","=B8-B10*B11"},{"F5","'m"},{"D6","'Le"},{"E6","=B9-B10*B12"},{"F6","'m"},
    {"D7","'Ae"},{"E7","=E5*E6"},{"F7","'m2"},{"D8","'h limite"},{"E8","=HLIMSUM(E7,E6,E5)"},{"F8","'m"},
    {"D9","'modo"},{"E9","=IF(E3<E8,\"VERTEDERO\",\"ORIFICIO\")"},{"D10","'Q sumid"},{"E10","=QSUMF(E6,E5,E3)"},{"F10","'m3/s"},
    {"D11","'c) eta d"},{"E11","=B14*B13"},{"D12","'bc"},{"E12","=IF(B16=1,0.5*B15,B15)+B17"},{"F12","'m"},
    {"D13","'L sep"},{"E13","=SEPSUM(E11,E2,B6,B7,E12)"},{"F13","'m"},{"D14","'Q captado"},{"E14","=E11*E2"},{"F14","'m3/s"},
    {"D15","'Qcap<=Qsum?"},{"E15","=IF(E14<=E10,1,0)"}}
  T("CERT P3a-c Cuneta/sumidero", "capacidad, modo y separacion", L)
end
do -- P3d-e colectores en serie + dimensionamiento
  local L = {{"A1","'C"},{"B1","'A km2"},{"C1","'tc min"},{"D1","'L m"},{"E1","'v m/s"},{"F1","'Cmed"},{"G1","'Aacum"},{"H1","'tcI"},
    {"I1","'i(tcI)"},{"J1","'QiI"},{"K1","'i(tci)"},{"L1","'Qi"},{"M1","'QD"},{"N1","'D adopt"},{"O1","'Q0.8D"},{"P1","'yn/D"},
    {"A9","'IDF a"},{"B9","900"},{"A10","'IDF b"},{"B10","10"},{"A11","'IDF c"},{"B11","0.75"},{"C9","'i=a/(t+b)^c"},
    {"A12","'n tubo"},{"B12","0.013"},{"A13","'S tubo"},{"B13","0.005"},{"A14","'y/D max"},{"B14","0.8"},
    {"D9","'col"},{"E9","'v real"},{"F9","'Vc autol"},{"G9","'Vc>=0.6"},{"H9","'v<=vmax"},{"I9","'resp."},{"J9","'vmax"},{"J10","3.5"}}
  put(L, "A", 2, {0.70, 0.60, 0.80}); put(L, "B", 2, {0.04, 0.06, 0.05}); put(L, "C", 2, {10, 14, 12})
  put(L, "D", 2, {250, 300}); put(L, "E", 2, {1.2, 1.5})
  L[#L+1] = {"F2", "=A2"}; L[#L+1] = {"G2", "=B2"}; L[#L+1] = {"H2", "=C2"}
  fillCol(L, "F", 3, 7, "=IF(A3>0,SUM($A$2:A3)/COUNT($A$2:A3),\"\")")
  fillCol(L, "G", 3, 7, "=IF(A3>0,G2+B3,\"\")")
  fillCol(L, "H", 3, 7, "=IF(A3>0,MAX(H2+D2/E2/60,C3),\"\")")
  fillCol(L, "I", 2, 7, "=IF(A2>0,IDF($B$9,$B$10,$B$11,H2),\"\")")
  fillCol(L, "J", 2, 7, "=IF(A2>0,RACIONAL(F2,I2,G2),\"\")")
  fillCol(L, "K", 2, 7, "=IF(A2>0,IDF($B$9,$B$10,$B$11,C2),\"\")")
  fillCol(L, "L", 2, 7, "=IF(A2>0,RACIONAL(A2,K2,B2),\"\")")
  L[#L+1] = {"M2", "=MAX(J2,L2)"}
  fillCol(L, "M", 3, 7, "=IF(A3>0,MAX(M2,J3,L3),\"\")")
  fillCol(L, "N", 2, 7, "=IF(A2>0,DCOM(M2,$B$12,$B$13,$B$14),\"\")")
  fillCol(L, "O", 2, 7, "=IF(A2>0,QPARC(N2,$B$12,$B$13,$B$14),\"\")")
  fillCol(L, "P", 2, 7, "=IF(A2>0,YND(M2,$B$12,$B$13,N2),\"\")")
  for i = 1, 6 do
    local r, s = 9 + i, 1 + i
    L[#L+1] = {"D" .. r, "=IF(A" .. s .. ">0," .. i .. ",\"\")"}
    L[#L+1] = {"E" .. r, "=IF(A" .. s .. ">0,VCIRC(M" .. s .. ",$B$12,$B$13,N" .. s .. "),\"\")"}
    L[#L+1] = {"F" .. r, "=IF(A" .. s .. ">0,VAUTOL(N" .. s .. ",$B$13,$B$12),\"\")"}
    L[#L+1] = {"G" .. r, "=IF(A" .. s .. ">0,IF(F" .. r .. ">=0.6,1,0),\"\")"}
    L[#L+1] = {"H" .. r, "=IF(A" .. s .. ">0,IF(E" .. r .. "<=$J$10,1,0),\"\")"}
    L[#L+1] = {"I" .. r, "=IF(A" .. s .. ">0,IF(N" .. s .. "<=0.8,\"SERVIU\",\"DOH\"),\"\")"}
  end
  T("CERT P3d-e Colectores", "QD en serie + diametro + v", L)
end
do -- Puente: Yarnell
  local L = {{"A1","'K pila"},{"B1","0.9"},{"C1","'forma (tabla)"},{"A2","'Q"},{"B2","150"},{"C2","'m3/s"},{"A3","'b cauce"},{"B3","40"},{"C3","'m"},
    {"A4","'h aguas ab."},{"B4","3"},{"C4","'m (h3)"},{"A5","'n pilas"},{"B5","3"},{"A6","'ancho pila"},{"B6","1.2"},{"C6","'m"},
    {"A8","'sigma=nD/b"},{"B8","=B5*B6/B3"},{"A9","'Fr3"},{"B9","=B2/(B3*B4)/SQRT(GRAV*B4)"},
    {"A10","'dh Yarnell"},{"B10","=YARNELL(B1,B2,B3,B4,B5,B6)"},{"C10","'m"},{"A11","'h aguas arr"},{"B11","=B4+B10"},{"C11","'m"},
    {"A12","'v pilas"},{"B12","=B2/((B3-B5*B6)*B4)"},{"C12","'m/s"}}
  T("Puente: remanso Yarnell", "sobreelevacion por pilas", L)
end
do -- Riesgo y periodo de retorno
  local L = {{"A1","'T"},{"B1","'n vida"},{"C1","'riesgo"},{"E1","'R objetivo"},{"F1","0.2"},{"E2","'n"},{"F2","30"},{"E3","'T requerido"},{"F3","=TRIESGO(F1,F2)"}}
  put(L, "A", 2, {2, 5, 10, 25, 50, 100, 200, 500}); for i = 2, 9 do L[#L+1] = {"B" .. i, "30"} end
  fillCol(L, "C", 2, 9, "=RIESGO(A2,B2)")
  T("Riesgo / periodo de retorno", "R=1-(1-1/T)^n", L)
end
do -- Prueba de gasto constante (Cooper-Jacob)
  local L = {{"A1","'t"},{"B1","'dh"},{"C1","'ln t"},{"E1","'Q"},{"F1","0.05"},{"E2","'r"},{"F2","20"},
    {"A2","60"},{"B2","0.8"},{"A3","120"},{"B3","0.9"},{"A4","300"},{"B4","1.03"},{"A5","600"},{"B5","1.13"},{"A6","1200"},{"B6","1.24"},
    {"E4","'a pend"},{"F4","=SLOPE(B2:B30,C2:C30)"},{"E5","'b interc"},{"F5","=INTERCEPT(B2:B30,C2:C30)"},
    {"E6","'dh/ciclo"},{"F6","=2.303*F4"},{"E7","'T"},{"F7","=F1/(4*PI*F4)"},{"E8","'t0"},{"F8","=EXP(-F5/F4)"},
    {"E9","'S"},{"F9","=2.24*F7*F8/F2^2"},{"E10","'u(tmin)"},{"F10","=U(F2,F9,F7,A2)"},{"E11","'R^2"},{"F11","=RSQ(B2:B30,C2:C30)"}}
  fillCol(L, "C", 2, 30, "=IF(A2>0,LN(A2),\"\")")
  T("Gasto cte (Cooper-Jacob)", "t y dh en A,B; resultados en F", L)
end
do -- Recuperacion
  local L = {{"A1","'t"},{"B1","'s'"},{"C1","'t/t'"},{"D1","'ln t/t'"},{"F1","'tf"},{"G1","3600"},{"F2","'Q"},{"G2","0.05"},
    {"A2","3700"},{"B2","0.6"},{"A3","3900"},{"B3","0.45"},{"A4","4200"},{"B4","0.35"},{"A5","5000"},{"B5","0.25"},{"A6","7200"},{"B6","0.15"},
    {"F4","'a"},{"G4","=SLOPE(B2:B30,D2:D30)"},{"F5","'T"},{"G5","=G2/(4*PI*G4)"}}
  fillCol(L, "C", 2, 30, "=IF(A2>$G$1,A2/(A2-$G$1),\"\")")
  fillCol(L, "D", 2, 30, "=LN(C2)")
  T("Recuperacion", "t (desde inicio) y s' en A,B", L)
end
do -- Superposicion de pozos / imagenes
  local L = {{"A1","'x"},{"B1","'y"},{"C1","'Q(+/-)"},{"D1","'r"},{"E1","'dh"},{"G1","'T"},{"H1","0.01"},{"G2","'R"},{"H2","500"},
    {"G3","'xp"},{"H3","30"},{"G4","'yp"},{"H4","10"},{"G6","'dhT"},{"H6","=SUM(E2:E20)"},{"G7","'H libre"},{"H7","0"},
    {"G8","'h libre"},{"H8","=IF(H7>0,SQRT(H7^2-2*H7*H6),\"\")"},
    {"A2","0"},{"B2","0"},{"C2","0.02"},{"A3","100"},{"B3","0"},{"C3","-0.02"}}
  fillCol(L, "D", 2, 20, "=DIST(A2,B2,$H$3,$H$4)")
  fillCol(L, "E", 2, 20, "=IF(C2<>0,THIEM(C2,$H$1,$H$2,D2),0)")
  T("Superposicion / imagenes", "pozos en A:C (imagen con signo)", L)
end
do -- Superposicion transiente
  local L = {{"A1","'x"},{"B1","'y"},{"C1","'Q"},{"D1","'t ini"},{"E1","'r"},{"F1","'dh"},{"H1","'T"},{"I1","0.005"},{"H2","'S"},{"I2","1E-4"},
    {"H3","'t"},{"I3","86400"},{"H4","'xp"},{"I4","50"},{"H5","'yp"},{"I5","0"},{"H7","'dhT"},{"I7","=SUM(F2:F20)"},
    {"A2","0"},{"B2","0"},{"C2","0.02"},{"D2","0"}}
  fillCol(L, "E", 2, 20, "=DIST(A2,B2,$I$4,$I$5)")
  fillCol(L, "F", 2, 20, "=IF(C2<>0,JACOB(C2,$I$1,$I$2,E2,$I$3-D2),0)")
  T("Superposicion transiente", "Jacob con tiempos de inicio", L)
end
do -- Bombeo variable
  local L = {{"A1","'dQ"},{"B1","'t i"},{"C1","'dh"},{"E1","'T"},{"F1","0.005"},{"E2","'S"},{"F2","1E-4"},{"E3","'r"},{"F3","0.2"},
    {"E4","'t"},{"F4","7200"},{"E6","'dh total"},{"F6","=SUM(C2:C20)"},{"A2","0.02"},{"B2","0"},{"A3","-0.02"},{"B3","3600"}}
  fillCol(L, "C", 2, 20, "=IF(A2<>0,THEIS(A2,$F$1,$F$2,$F$3,$F$4-B2),0)")
  T("Bombeo variable / recup.", "+Q enciende, -Q apaga", L)
end
do -- Prueba escalonada
  local L = {{"A1","'Q"},{"B1","'dh"},{"C1","'dh/Q"},{"E1","'C perd"},{"F1","=SLOPE(C2:C20,A2:A20)"},{"E2","'B acuif"},{"F2","=INTERCEPT(C2:C20,A2:A20)"},
    {"E3","'80%Qmax"},{"F3","=0.8*MAX(A2:A20)"},{"A2","0.01"},{"B2","2.6"},{"A3","0.02"},{"B3","5.4"},{"A4","0.03"},{"B4","8.3"},{"A5","0.04"},{"B5","11.2"}}
  fillCol(L, "C", 2, 20, "=IF(A2>0,B2/A2,\"\")")
  T("Prueba escalonada", "dh = B Q + C Q^2", L)
end
do -- Estratos
  local L = {{"A1","'K"},{"B1","'m"},{"C1","'K*m"},{"D1","'m/K"},{"F1","'Keq par"},{"G1","=SUM(C2:C20)/SUM(B2:B20)"},
    {"F2","'Keq perp"},{"G2","=SUM(B2:B20)/SUM(D2:D20)"},{"F3","'T"},{"G3","=SUM(C2:C20)"},
    {"A2","1E-4"},{"B2","2"},{"A3","1E-6"},{"B3","1"},{"A4","5E-5"},{"B4","3"}}
  fillCol(L, "C", 2, 20, "=IF(A2>0,A2*B2,0)")
  fillCol(L, "D", 2, 20, "=IF(A2>0,B2/A2,0)")
  T("K equivalente estratos", "K y m por estrato", L)
end
do -- Cono de depresion
  local L = {{"A1","'r"},{"B1","'Thiem dh"},{"C1","'Dupuit h"},{"D1","'Theis dh"},{"F1","'Q"},{"G1","0.02"},{"F2","'T"},{"G2","0.005"},
    {"F3","'K"},{"G3","2E-4"},{"F4","'H"},{"G4","25"},{"F5","'R"},{"G5","300"},{"F6","'S"},{"G6","1E-4"},{"F7","'t"},{"G7","3600"}}
  local rs = {0.2, 1, 5, 10, 20, 50, 100, 200, 300}
  for i, r in ipairs(rs) do L[#L+1] = {"A" .. (i + 1), tostring(r)} end
  fillCol(L, "B", 2, 10, "=THIEM($G$1,$G$2,$G$5,A2)")
  fillCol(L, "C", 2, 10, "=DUPUIT($G$4,$G$1,$G$3,$G$5,A2)")
  fillCol(L, "D", 2, 10, "=THEIS($G$1,$G$2,$G$6,A2,$G$7)")
  T("Cono de depresion", "perfiles vs r", L)
end
do -- Theis dh vs t
  local L = {{"A1","'t"},{"B1","'u"},{"C1","'W(u)"},{"D1","'dh"},{"F1","'Q"},{"G1","0.02"},{"F2","'T"},{"G2","0.005"},{"F3","'S"},{"G3","1E-4"},{"F4","'r"},{"G4","50"}}
  local ts = {60, 120, 300, 600, 1200, 1800, 3600, 7200, 14400, 28800, 86400}
  for i, t in ipairs(ts) do L[#L+1] = {"A" .. (i + 1), tostring(t)} end
  fillCol(L, "B", 2, 12, "=U($G$4,$G$3,$G$2,A2)")
  fillCol(L, "C", 2, 12, "=W(B2)")
  fillCol(L, "D", 2, 12, "=$G$1/(4*PI*$G$2)*C2")
  T("Theis: dh vs t", "tabla u, W(u), dh", L)
end
do -- Equilibrio dinamico
  local L = {{"A1","'mes"},{"B1","'R mm"},{"C1","'y0"},{"D1","'yf"},{"F1","'K m/d"},{"G1","0.5"},{"F2","'d"},{"G2","2"},{"F3","'S"},{"G3","0.05"},
    {"F4","'D"},{"G4","25"},{"F5","'dias"},{"G5","30"},{"F6","'y ini"},{"G6","0.2"},{"F7","'alfa"},{"G7","=PI^2*G1*G2/(G3*G4^2)"},
    {"F9","'yf<=y0 ?"},{"G9","=IF(D13<=C2,1,0)"}}
  local R = {50, 80, 120, 60, 10, 0, 0, 0, 20, 40, 60, 70}
  for i = 1, 12 do L[#L+1] = {"A" .. (i + 1), tostring(i)}; L[#L+1] = {"B" .. (i + 1), tostring(R[i])} end
  L[#L+1] = {"C2", "=$G$6+B2/1000/$G$3"}
  fillCol(L, "C", 3, 13, "=D2+B3/1000/$G$3")
  fillCol(L, "D", 2, 13, "=1.16*C2*EXP(-$G$7*$G$5)")
  T("Equilibrio dinamico", "Glover mes a mes", L)
end
do -- Iteracion Glover + Hooghoudt
  local L = {{"A1","'iter"},{"B1","'d'"},{"C1","'D"},{"E1","'K"},{"F1","0.5"},{"E2","'d"},{"F2","3"},{"E3","'S"},{"F3","0.05"},
    {"E4","'t"},{"F4","3"},{"E5","'y0"},{"F5","1"},{"E6","'y"},{"F6","0.5"},{"E7","'r0"},{"F7","0.1"},
    {"A2","0"},{"B2","=F2"},{"C2","=GLOVERD($F$1,B2,$F$4,$F$3,$F$5,$F$6)"}}
  for i = 3, 12 do L[#L+1] = {"A" .. i, tostring(i - 2)} end
  fillCol(L, "B", 3, 12, "=DPRIME($F$2,C2,$F$7)")
  fillCol(L, "C", 3, 12, "=GLOVERD($F$1,B3,$F$4,$F$3,$F$5,$F$6)")
  T("Iteracion Glover-Hooghoudt", "D converge en columna C", L)
end
do -- Cuencas en serie
  local L = {{"A1","'C"},{"B1","'A km2"},{"C1","'tc"},{"D1","'L"},{"E1","'v"},{"F1","'Cmed"},{"G1","'Aacum"},{"H1","'tcI"},{"I1","'i(tcI)"},{"J1","'QI"},
    {"K1","'Qi"},{"L1","'i(tci)"},{"M1","'QD"},
    {"A2","0.5"},{"B2","0.8"},{"C2","15"},{"D2","600"},{"E2","1.5"},{"A3","0.6"},{"B3","0.5"},{"C3","12"},{"D3","800"},{"E3","1.5"},
    {"A4","0.4"},{"B4","1.2"},{"C4","20"},{"N1","'IDF i=a/(t+b)^c"},{"N2","'a"},{"O2","1500"},{"N3","'b"},{"O3","10"},{"N4","'c"},{"O4","0.8"},
    {"F2","=A2"},{"G2","=B2"},{"H2","=C2"}}
  fillCol(L, "F", 3, 10, "=IF(A3>0,SUM($A$2:A3)/COUNT($A$2:A3),\"\")")
  fillCol(L, "G", 3, 10, "=IF(A3>0,G2+B3,\"\")")
  fillCol(L, "H", 3, 10, "=IF(A3>0,MAX(H2+D2/E2/60,C3),\"\")")
  fillCol(L, "I", 2, 10, "=$O$2/(H2+$O$3)^$O$4")
  fillCol(L, "J", 2, 10, "=RACIONAL(F2,I2,G2)")
  fillCol(L, "L", 2, 10, "=IF(A2>0,$O$2/(C2+$O$3)^$O$4,\"\")")
  fillCol(L, "K", 2, 10, "=IF(A2>0,RACIONAL(A2,L2,B2),\"\")")
  L[#L+1] = {"M2", "=MAX(J2,K2)"}
  fillCol(L, "M", 3, 10, "=MAX(M2,J3,K3)")
  T("Cuencas en serie (racional)", "tt=L/v en min; IDF editable", L)
end
do -- Canal Manning
  local L = {{"A1","'tipo"},{"B1","2"},{"C1","'1rect 2trap 3tri 4circ"},{"A2","'Q"},{"B2","1"},{"A3","'n"},{"B3","0.025"},{"A4","'S"},{"B4","0.001"},
    {"A5","'b"},{"B5","1"},{"A6","'k"},{"B6","1.5"},{"A7","'D"},{"B7","0"},
    {"A9","'yn"},{"B9","=YN(B1,B2,B3,B4,B5,B6,B7)"},{"A10","'yc"},{"B10","=YC(B1,B2,B5,B6,B7)"},
    {"A11","'A"},{"B11","=AREA(B1,B5,B6,B7,B9)"},{"A12","'P"},{"B12","=PERIM(B1,B5,B6,B7,B9)"},{"A13","'Rh"},{"B13","=B11/B12"},
    {"A14","'v"},{"B14","=B2/B11"},{"A15","'Fr"},{"B15","=FROUDE(B1,B5,B6,B7,B9,B2)"},{"A16","'rio?"},{"B16","=IF(B9>B10,1,0)"},
    {"A17","'revancha"},{"B17","=0.2*B9"}}
  T("Canal (Manning yn, yc)", "cambie datos en columna B", L)
end
do -- Alcantarilla
  local L = {{"A1","'Q1 (T1)"},{"B1","3"},{"A2","'Q2 (T2)"},{"B2","4"},{"A3","'b"},{"B3","1.5"},{"A4","'n"},{"B4","0.015"},{"A5","'SD"},{"B5","0.01"},
    {"A6","'ke"},{"B6","0.5"},{"A7","'L"},{"B7","20"},
    {"A9","'D=H"},{"B9","=DRECT(B1,B3)"},{"A10","'hc"},{"B10","=(B1^2/(GRAV*B3^2))^(1/3)"},{"A11","'v entrada"},{"B11","=B1/(B3*B10)"},
    {"A12","'Sc"},{"B12","=SCRECT(B9,B4,B3)"},{"A13","'SD>=Sc?"},{"B13","=IF(B5>=B12,1,0)"},
    {"A15","'Sn (Q2)"},{"B15","=SNRECT(B2,B4,B3,B9)"},{"A16","'caso"},{"B16","=IF(B5>B15,1,IF(B5<B15,3,2))"},
    {"A17","'H' caso1"},{"B17","=HCASO1(B9,B2,B3)"},{"A18","'H' caso2"},{"B18","=HCASO2(B9,B6,B2,B3)"},
    {"A19","'v llena"},{"B19","=B2/(B3*B9)"},{"A20","'J"},{"B20","=B19^2*B4^2/(B3*B9/(2*B3+2*B9))^(4/3)"},
    {"A21","'H' caso3"},{"B21","=B9+(1+B6)*B19^2/(2*GRAV)+B20*B7-B5*B7"}}
  T("Alcantarilla rectangular", "criterio 1 y 2", L)
end
do -- Tiempos de concentracion
  local L = {{"A1","'L km"},{"B1","2.5"},{"A2","'S %"},{"B2","4"},{"A3","'H m"},{"B3","120"},{"A4","'A km2"},{"B4","1.8"},{"A5","'Hm m"},{"B5","60"},
    {"A6","'CN"},{"B6","75"},{"A8","'Esp."},{"B8","=TCESP(B1,B2)"},{"A9","'Calif."},{"B9","=TCCALIF(B1,B3)"},
    {"A10","'Giand."},{"B10","=TCGIAND(B4,B1,B5)"},{"A11","'SCS"},{"B11","=TCSCS(B1,B6,B2)"},{"A12","'prom"},{"B12","=AVG(B8:B11)"}}
  T("Tiempos de concentracion", "tc en minutos", L)
end
do -- Drenes permanente
  local L = {{"A1","'h0"},{"B1","0.3"},{"A2","'d"},{"B2","0.5"},{"A3","'f"},{"B3","5E-8"},{"A4","'D"},{"B4","30"},{"A5","'K"},{"B5","1E-5"},{"A6","'r0"},{"B6","0.1"},
    {"A8","'Donnan H"},{"B8","=DONNAN(B1,B3,B4,B5)"},{"A9","'Hoogh H"},{"B9","=HOOGH(B1,B2,B3,B4,B5)"},{"A10","'alfa Dagan"},{"B10","=DAGAN(B6,B2)"},
    {"A11","'Dagan H"},{"B11","=SQRT(B2^2+B3*B4^2/(4*B5)*(1+B10*4*B2/B4))-B2"},{"A12","'Kirk F"},{"B12","=KIRKF(2*B6/B4,B2/B4)"},
    {"A13","'Kirk H"},{"B13","=B4*B3/B5*B12"}}
  T("Drenes regimen permanente", "compara formulas", L)
end

-- =====================================================================
--  Interfaz de la planilla
-- =====================================================================
local selC, selR = 1, 1
local offC, offR = 1, 1
local editing, buf = false, ""
local mode = "sheet"   -- sheet | dialog | msg | plot | lista
local dlg              -- dialogo activo {title, fields, vals, fi, ok}
local msg = {}         -- lineas del mensaje
local mtop = 1
local lista            -- {title, items, sel, top, pick}
local plots, pk
local clip
local showFormulas = false
local CW, RH, HW = 58, 13, 20

local function inv() platform.window:invalidate() end
local function Wd() return platform.window:width() end
local function Hd() return platform.window:height() end
local function visCols() return math.max(1, math.floor((Wd() - HW)/CW)) end
local function visRows() return math.max(1, math.floor((Hd() - 30 - 14 - 12)/RH)) end

local function cellText(k)
  local raw = cells[k]
  if not raw then return "" end
  if showFormulas and raw:sub(1, 1) == "=" then return raw end
  local v = valueOf(k)
  if type(v) == "number" then
    local s = string.format("%.5g", v)
    if #s > 9 then s = string.format("%.2e", v) end
    return s
  end
  return tostring(v or "")
end

local function showMsg(lines) msg = lines; mtop = 1; mode = "msg" end
local function openDialog(title, fields, okfn)
  dlg = {title = title, fields = fields, vals = {}, fi = 1, ok = okfn}
  for i, f in ipairs(fields) do dlg.vals[i] = f[2] or "" end
  mode = "dialog"
end
local function openList(title, items, pick) lista = {title = title, items = items, sel = 1, top = 1, pick = pick}; mode = "lista" end

local function commitEdit()
  if editing then setCell(key(selC, selR), buf); editing = false; buf = "" end
end
local function moveSel(dc, dr)
  selC = math.max(1, math.min(NCOL, selC + dc))
  selR = math.max(1, math.min(NROW, selR + dr))
  if selC < offC then offC = selC end
  if selC >= offC + visCols() then offC = selC - visCols() + 1 end
  if selR < offR then offR = selR end
  if selR >= offR + visRows() then offR = selR - visRows() + 1 end
end

-- ----------------------- Acciones de menu ---------------------------
local A = {}
function A.nueva() cells = {}; recalc(); selC, selR, offC, offR = 1, 1, 1, 1 end
function A.copiar() clip = {raw = cells[key(selC, selR)], c = selC, r = selR} end
function A.pegar()
  if clip and clip.raw then setCell(key(selC, selR), shiftRefs(clip.raw, selR - clip.r, selC - clip.c)) end
end
function A.borrar() setCell(key(selC, selR), nil) end
function A.rellenar()
  openDialog("Rellenar hacia abajo", {{"Hasta fila", tostring(math.min(NROW, selR + 10))}}, function(v)
    local r2 = math.floor(evalstr(v[1]))
    local raw = cells[key(selC, selR)]
    if not raw then error("La celda esta vacia") end
    for r = selR + 1, math.min(NROW, r2) do cells[key(selC, r)] = shiftRefs(raw, r - selR, 0) end
    recalc()
  end)
end
function A.rellenarSerie()
  openDialog("Serie numerica", {{"Desde", "1"}, {"Paso", "1"}, {"Hasta fila", tostring(selR + 9)}}, function(v)
    local a, s, r2 = evalstr(v[1]), evalstr(v[2]), math.floor(evalstr(v[3]))
    for r = selR, math.min(NROW, r2) do cells[key(selC, r)] = tostring(a + (r - selR)*s) end
    recalc()
  end)
end
function A.objetivo()
  openDialog("Buscar objetivo", {{"Celda objetivo (formula)", key(selC, selR)}, {"Valor deseado", "0"}, {"Celda a cambiar (dato)", ""}},
    function(v)
      local tc, tr = parseRef(v[1]); local vc, vr = parseRef(v[3])
      if not tc or not vc then error("Referencia invalida") end
      local target = evalstr(v[2])
      local kT, kV = key(tc, tr), key(vc, vr)
      local old = cells[kV]
      local x0 = tonumber(old or "") or 1
      local function f(x) cells[kV] = string.format("%.15g", x); recalc(); local y = valueOf(kT)
        if type(y) ~= "number" then return 0/0 end; return y - target end
      local rs = roots(f, 1e-10, 1e10, 1200, true)
      local rn = roots(function(x) return f(-x) end, 1e-10, 1e10, 1200, true)
      for _, x in ipairs(rn) do rs[#rs+1] = -x end
      if #rs == 0 then cells[kV] = old; recalc(); error("No se encontro solucion") end
      local best = rs[1]
      for _, x in ipairs(rs) do if abs(x - x0) < abs(best - x0) then best = x end end
      cells[kV] = string.format("%.10g", best); recalc()
      local out = {kV .. " = " .. fmt(best), kT .. " = " .. fmt(valueOf(kT))}
      if #rs > 1 then out[#out+1] = "Otras soluciones:"; for _, x in ipairs(rs) do if x ~= best then out[#out+1] = "  " .. fmt(x) end end end
      return out
    end)
end
local function tryNum(x) return type(x) == "number" and ok(x) end
function A.graficar()
  openDialog("Graficar", {{"Rango X", "A2:A20"}, {"Rango Y", "B2:B20"}, {"Rango Y2 (opcional)", ""},
      {"Escala X log (0/1)", "0"}, {"Y invertido (0/1)", "0"}, {"Ajuste 0 no,1 lineal,2 ln(x)", "0"}},
    function(v)
      local xs = rangeVals(v[1]); local ys = rangeVals(v[2])
      local P = {title = "Grafico " .. v[2] .. " vs " .. v[1], xl = v[1], yl = v[2], series = {}, xlog = v[4] == "1", yinv = v[5] == "1"}
      local s1 = {name = "Y", pts = {}, style = "both"}
      local X, Y = {}, {}
      for i = 1, math.min(#xs, #ys) do if tryNum(xs[i]) and tryNum(ys[i]) then s1.pts[#s1.pts+1] = {xs[i], ys[i]}; X[#X+1] = xs[i]; Y[#Y+1] = ys[i] end end
      P.series[1] = s1
      if v[3] ~= "" then
        local y2 = rangeVals(v[3]); local s2 = {name = "Y2", pts = {}, style = "both"}
        for i = 1, math.min(#xs, #y2) do if tryNum(xs[i]) and tryNum(y2[i]) then s2.pts[#s2.pts+1] = {xs[i], y2[i]} end end
        P.series[#P.series+1] = s2
      end
      local info
      if v[6] == "1" or v[6] == "2" then
        local XX = {}
        for i = 1, #X do XX[i] = (v[6] == "2") and ln(X[i]) or X[i] end
        local a, b = regresion(XX, Y)
        local fit = {name = "ajuste", pts = {}, style = "line"}
        local mn, mx = math.huge, -math.huge
        for _, x in ipairs(X) do mn = math.min(mn, x); mx = math.max(mx, x) end
        for i = 0, 40 do local x = mn + (mx - mn)*i/40; if v[6] == "2" then x = mn*(mx/mn)^(i/40) end
          fit.pts[#fit.pts+1] = {x, a*((v[6] == "2") and ln(x) or x) + b} end
        P.series[#P.series+1] = fit
        info = "y = " .. fmt(a) .. ((v[6] == "2") and " ln(x) + " or " x + ") .. fmt(b)
        P.info = function(PP) local s = PP.series[PP.cs or 1]; local p = s.pts[PP.ci or 1]
          return info .. (p and ("  (" .. fmt(p[1]) .. ", " .. fmt(p[2]) .. ")") or "") end
      end
      plots, pk = {P}, 1
      mode = "plot"
      return nil
    end)
end
function A.exportar()
  openDialog("Exportar columna a lista TI", {{"Rango", colName(selC) .. "2:" .. colName(selC) .. "20"}, {"Nombre de lista", "col" .. colName(selC):lower()}},
    function(v)
      local t = {}
      for _, x in ipairs(rangeVals(v[1])) do if type(x) == "number" then t[#t+1] = x end end
      if var and var.store then var.store(v[2], t) end
      return {"Lista " .. v[2] .. " guardada (" .. #t .. " datos)."}
    end)
end
function A.plantillas()
  local items = {}
  for i, t in ipairs(TPL) do items[i] = t.name end
  openList("Plantillas", items, function(i)
    local t = TPL[i]
    cells = {}
    for _, p in ipairs(t.list) do cells[p[1]] = p[2] end
    recalc(); selC, selR, offC, offR = 1, 1, 1, 1
    mode = "sheet"
  end)
end
function A.funciones()
  local items = {}
  for i, f in ipairs(FHELP) do items[i] = f[1] .. (f[2] ~= "" and ("  " .. f[2]) or "") end
  openList("Funciones (enter inserta)", items, function(i)
    local name = FHELP[i][1]:match("^([%u%d]+)")
    mode = "sheet"
    if name then
      if not editing then editing = true; buf = cells[key(selC, selR)] or "=" end
      if buf == "" then buf = "=" end
      buf = buf .. name .. (FUNC[name] and type(FUNC[name]) == "function" and "(" or "")
    end
  end)
end
function A.ayuda()
  showMsg({"PLANILLA CIV-346 (tipo Excel)",
    "- Flechas: mover. Escribir: edita la celda.",
    "- enter: confirma (o edita la celda actual).",
    "- Formulas empiezan con = ej: =A2*2+B3",
    "- Rangos: SUM(A2:A10). $ fija: $B$1",
    "- Funciones hidraulicas: menu Funciones",
    "- del: borra caracter / celda. esc: cancela.",
    "- [menu]: plantillas, copiar, pegar, rellenar,",
    "  buscar objetivo, graficar, exportar a lista.",
    "- Texto: empiece con ' (ej: 'Caudal)",
    "- La hoja se guarda con el documento .tns"})
end
function A.verFormulas() showFormulas = not showFormulas end

local MENU = {
  {"Hoja", {{"Ayuda", A.ayuda}, {"Nueva hoja", A.nueva}, {"Ver formulas/valores", A.verFormulas}, {"Exportar a lista TI", A.exportar}}},
  {"Editar", {{"Copiar celda", A.copiar}, {"Pegar celda", A.pegar}, {"Rellenar abajo", A.rellenar},
              {"Serie numerica", A.rellenarSerie}, {"Borrar celda", A.borrar}}},
  {"Herramientas", {{"Buscar objetivo", A.objetivo}, {"Graficar / regresion", A.graficar}, {"Funciones", A.funciones}}},
  {"Plantillas", {{"Elegir plantilla", A.plantillas}}},
}
-- menu propio (tecla menu en la calculadora usa toolpalette)
local function runAction(fn)
  commitEdit()
  local okc, err = pcall(fn)
  if not okc then showMsg({"ERROR:", (tostring(err):gsub("^.-:%d+: ", ""))}) end
  inv()
end
if toolpalette and toolpalette.register then
  local tp = {}
  for _, m in ipairs(MENU) do
    local sub = {m[1]}
    for _, it in ipairs(m[2]) do sub[#sub+1] = {it[1], function() runAction(it[2]) end} end
    tp[#tp+1] = sub
  end
  pcall(toolpalette.register, tp)
end
local function menuPropio()
  local items, acts = {}, {}
  for _, m in ipairs(MENU) do for _, it in ipairs(m[2]) do items[#items+1] = m[1] .. ": " .. it[1]; acts[#acts+1] = it[2] end end
  openList("Menu", items, function(i) mode = "sheet"; runAction(acts[i]) end)
end

-- ----------------------------- Dibujo --------------------------------
local function header(gc, title)
  gc:setColorRGB(0, 69, 137); gc:fillRect(0, 0, Wd(), 16)
  gc:setColorRGB(255, 255, 255); gc:setFont("sansserif", "b", 9); gc:drawString(title, 3, 1, "top")
  gc:setFont("sansserif", "r", 9)
end
local function footer(gc, t)
  gc:setColorRGB(230, 230, 230); gc:fillRect(0, Hd() - 12, Wd(), 12)
  gc:setColorRGB(60, 60, 60); gc:setFont("sansserif", "r", 7); gc:drawString(t, 2, Hd() - 12, "top")
  gc:setFont("sansserif", "r", 9)
end
local function clipText(gc, s, w)
  if gc:getStringWidth(s) <= w then return s end
  while #s > 0 and gc:getStringWidth(s .. "~") > w do s = s:sub(1, -2) end
  return s .. "~"
end

local function drawSheet(gc)
  local k = key(selC, selR)
  -- barra de formula
  gc:setColorRGB(0, 69, 137); gc:fillRect(0, 0, Wd(), 16)
  gc:setColorRGB(255, 255, 255); gc:setFont("sansserif", "b", 9)
  local content = editing and (buf .. "_") or (cells[k] or "")
  gc:drawString(clipText(gc, k .. ": " .. content, Wd() - 6), 3, 1, "top")
  gc:setFont("sansserif", "r", 8)
  gc:setColorRGB(255, 250, 225); gc:fillRect(0, 16, Wd(), 14)
  local v = valueOf(k)
  gc:setColorRGB(120, 0, 0)
  gc:drawString("= " .. ((type(v) == "number") and fmt(v) or tostring(v or "")), 3, 16, "top")
  -- cabeceras
  local y0 = 30
  local nc, nr = visCols(), visRows()
  gc:setColorRGB(215, 222, 235); gc:fillRect(0, y0, Wd(), RH); gc:fillRect(0, y0, HW, RH*(nr + 1))
  gc:setColorRGB(0, 0, 0)
  for i = 0, nc - 1 do
    local c = offC + i
    if c > NCOL then break end
    local x = HW + i*CW
    if c == selC then gc:setColorRGB(247, 176, 0); gc:fillRect(x, y0, CW, RH); gc:setColorRGB(0, 0, 0) end
    gc:drawString(colName(c), x + CW/2 - 3, y0, "top")
  end
  for j = 0, nr - 1 do
    local r = offR + j
    local y = y0 + RH*(j + 1)
    if r == selR then gc:setColorRGB(247, 176, 0); gc:fillRect(0, y, HW, RH); gc:setColorRGB(0, 0, 0) end
    gc:drawString(tostring(r), 2, y, "top")
    for i = 0, nc - 1 do
      local c = offC + i
      if c > NCOL then break end
      local x = HW + i*CW
      local kk = key(c, r)
      if c == selC and r == selR then gc:setColorRGB(255, 235, 170); gc:fillRect(x, y, CW, RH) end
      local raw = cells[kk]
      local t = cellText(kk)
      if t ~= "" then
        local isnum = type(valueOf(kk)) == "number" and not (showFormulas and raw and raw:sub(1, 1) == "=")
        if raw and raw:sub(1, 1) == "=" then gc:setColorRGB(0, 60, 150) else gc:setColorRGB(0, 0, 0) end
        if type(valueOf(kk)) == "string" and tostring(valueOf(kk)):sub(1, 1) == "#" then gc:setColorRGB(200, 0, 0) end
        t = clipText(gc, t, CW - 3)
        local xx = isnum and (x + CW - 2 - gc:getStringWidth(t)) or (x + 2)
        gc:drawString(t, xx, y, "top")
      end
    end
  end
  -- grilla
  gc:setColorRGB(200, 200, 200)
  for i = 0, nc do gc:drawLine(HW + i*CW, y0, HW + i*CW, y0 + RH*(nr + 1)) end
  for j = 0, nr + 1 do gc:drawLine(0, y0 + RH*j, HW + nc*CW, y0 + RH*j) end
  gc:setColorRGB(0, 0, 0)
  gc:drawRect(HW + (selC - offC)*CW, y0 + RH*(selR - offR + 1), CW, RH)
  footer(gc, editing and "enter: confirmar  esc: cancelar  del: borrar" or "menu: opciones  enter: editar  =: formula  del: borrar celda")
end

local function drawDialog(gc)
  header(gc, dlg.title)
  for i, f in ipairs(dlg.fields) do
    local y = 20 + (i - 1)*28
    gc:setColorRGB(60, 60, 60); gc:setFont("sansserif", "r", 8); gc:drawString(f[1], 6, y, "top")
    gc:setFont("sansserif", "r", 9)
    if i == dlg.fi then gc:setColorRGB(255, 235, 170) else gc:setColorRGB(245, 245, 245) end
    gc:fillRect(6, y + 11, Wd() - 12, 15)
    gc:setColorRGB(0, 0, 0); gc:drawRect(6, y + 11, Wd() - 12, 15)
    gc:drawString(dlg.vals[i] .. ((i == dlg.fi) and "_" or ""), 9, y + 11, "top")
  end
  footer(gc, "flechas: campo  enter: aceptar  esc: cancelar")
end

local function drawList(gc)
  header(gc, lista.title)
  local n = math.floor((Hd() - 30)/14)
  if lista.sel < lista.top then lista.top = lista.sel end
  if lista.sel > lista.top + n - 1 then lista.top = lista.sel - n + 1 end
  for i = 0, n - 1 do
    local idx = lista.top + i
    local it = lista.items[idx]
    if not it then break end
    local y = 17 + i*14
    if idx == lista.sel then gc:setColorRGB(247, 176, 0); gc:fillRect(0, y, Wd(), 14) end
    gc:setColorRGB(0, 0, 0); gc:drawString(clipText(gc, it, Wd() - 6), 4, y, "top")
  end
  footer(gc, "enter: elegir  esc: volver")
end

local function drawMsg(gc)
  header(gc, "Mensaje")
  local n = math.floor((Hd() - 30)/14)
  for i = 0, n - 1 do
    local l = msg[mtop + i]
    if not l then break end
    gc:setColorRGB(0, 0, 0); gc:drawString(l, 4, 18 + i*14, "top")
  end
  footer(gc, "enter/esc: volver")
end

function on.paint(gc)
  gc:setPen("thin", "smooth")
  gc:setFont("sansserif", "r", 9)
  if mode == "sheet" then drawSheet(gc)
  elseif mode == "dialog" then drawDialog(gc)
  elseif mode == "lista" then drawList(gc)
  elseif mode == "msg" then drawMsg(gc)
  elseif mode == "plot" then
    local P = plots[pk]
    header(gc, P.title or "Grafico")
    gc:setColorRGB(255, 255, 255); gc:fillRect(0, 16, Wd(), Hd() - 16)
    local info = plotInfo(P)
    gc:setFont("sansserif", "r", 8); gc:setColorRGB(150, 0, 0); gc:drawString(info, 3, 16, "top"); gc:setFont("sansserif", "r", 9)
    drawPlot(gc, P, 0, 28, Wd(), Hd() - 28 - 12)
    footer(gc, "<- -> cursor  ^v serie  esc: volver")
  end
end

-- ----------------------------- Teclado -------------------------------
local function dlgAccept()
  local okc, res = pcall(dlg.ok, dlg.vals)
  if not okc then showMsg({"ERROR:", (tostring(res):gsub("^.-:%d+: ", ""))})
  elseif type(res) == "table" then showMsg(res)
  elseif mode == "dialog" then mode = "sheet" end
end

function on.arrowUp()
  if mode == "sheet" then commitEdit(); moveSel(0, -1)
  elseif mode == "dialog" then dlg.fi = (dlg.fi - 2) % #dlg.fields + 1
  elseif mode == "lista" then lista.sel = (lista.sel - 2) % #lista.items + 1
  elseif mode == "msg" then mtop = math.max(1, mtop - 1)
  elseif mode == "plot" then local P = plots[pk]; P.cs = ((P.cs or 1) - 2) % #P.series + 1; P.ci = 1 end
  inv()
end
function on.arrowDown()
  if mode == "sheet" then commitEdit(); moveSel(0, 1)
  elseif mode == "dialog" then dlg.fi = dlg.fi % #dlg.fields + 1
  elseif mode == "lista" then lista.sel = lista.sel % #lista.items + 1
  elseif mode == "msg" then mtop = mtop + 1
  elseif mode == "plot" then local P = plots[pk]; P.cs = (P.cs or 1) % #P.series + 1; P.ci = 1 end
  inv()
end
function on.arrowLeft()
  if mode == "sheet" then commitEdit(); moveSel(-1, 0)
  elseif mode == "plot" then local P = plots[pk]; P.ci = math.max(1, (P.ci or 1) - 1) end
  inv()
end
function on.arrowRight()
  if mode == "sheet" then commitEdit(); moveSel(1, 0)
  elseif mode == "plot" then local P = plots[pk]; P.ci = math.min(#P.series[P.cs or 1].pts, (P.ci or 1) + 1) end
  inv()
end
function on.tabKey() if mode == "sheet" then commitEdit(); moveSel(1, 0) elseif mode == "dialog" then on.arrowDown() end; inv() end
function on.enterKey()
  if mode == "sheet" then
    if editing then commitEdit(); moveSel(0, 1) else editing = true; buf = cells[key(selC, selR)] or "" end
  elseif mode == "dialog" then dlgAccept()
  elseif mode == "lista" then local l = lista; l.pick(l.sel)
  elseif mode == "msg" or mode == "plot" then mode = "sheet" end
  inv()
end
on.returnKey = on.enterKey
function on.escapeKey()
  if mode == "sheet" then editing = false; buf = "" else mode = "sheet" end
  inv()
end
function on.backspaceKey()
  if mode == "sheet" then
    if editing then
      local p = #buf
      while p > 0 and buf:byte(p) >= 128 and buf:byte(p) < 192 do p = p - 1 end
      buf = buf:sub(1, p - 1)
    else A.borrar() end
  elseif mode == "dialog" then local s = dlg.vals[dlg.fi]; dlg.vals[dlg.fi] = s:sub(1, -2)
  else mode = "sheet" end
  inv()
end
function on.clearKey()
  if mode == "sheet" then if editing then buf = "" else A.borrar() end
  elseif mode == "dialog" then dlg.vals[dlg.fi] = "" end
  inv()
end
local MAPCH = {["\226\136\146"] = "-", ["\225\180\135"] = "E", ["\207\128"] = "PI", ["\195\151"] = "*", ["\195\183"] = "/",
  ["\226\136\154"] = "SQRT(", ["\194\178"] = "^2", ["\226\132\175"] = "EXP(1)"}
function on.charIn(ch)
  ch = MAPCH[ch] or ch
  if mode == "sheet" then
    if not editing and ch == "?" then menuPropio(); inv(); return end
    if not editing then editing = true; buf = "" end
    buf = buf .. ch
  elseif mode == "dialog" then dlg.vals[dlg.fi] = dlg.vals[dlg.fi] .. ch
  elseif mode == "lista" then local n = tonumber(ch); if n and n >= 1 and n <= #lista.items then lista.sel = n end
  end
  inv()
end
function on.contextMenu() menuPropio(); inv() end   -- ctrl+menu

-- Persistencia en el documento
function on.save() return {cells = cells} end
function on.restore(st) if type(st) == "table" and type(st.cells) == "table" then cells = st.cells; recalc() end end

-- Hoja inicial: ayuda rapida
cells = {A1 = "'CIV-346 PLANILLA", A2 = "'menu: plantillas", A3 = "'ej: =THIEM(0.05,0.01,500,10)", A4 = "=THIEM(0.05,0.01,500,10)"}

PLAN = {cells = function() return cells end, set = setCell, val = valueOf, A = A, TPL = TPL, FUNC = FUNC, shift = shiftRefs,
  dialog = function() return dlg end, mode = function() return mode end, lista = function() return lista end, msg = function() return msg end}
