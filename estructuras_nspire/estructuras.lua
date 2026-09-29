platform.apiLevel = "2.4"
--[[==========================================================================
  ESTRUCTURAS 2D  -  Analisis matricial de marcos y enrejados
  Para TI-Nspire CX II / CX II CAS (Lua, Script Editor)

  * Marcos (porticos, vigas continuas) hiperestaticos de cualquier grado
      - Articulaciones internas (liberacion de momento en extremos de barra)
      - Apoyos: empotrado, articulado, rodillos, guias, inclinados,
        resortes y asentamientos
      - Cargas: nodales, puntuales, momentos, distribuidas uniformes,
        triangulares, trapezoidales (parciales), en ejes locales o
        globales (o proyectadas), temperatura y errores de fabricacion
  * Enrejados (armaduras) isostaticos o hiperestaticos
  * Las cargas pueden ser numeros o letras (P, w, 2*P+5, w*L/2 ...).
    El resultado se entrega en forma literal por superposicion:
       M = 12.5 + 3.2*P - 0.5*w
  * Diagramas de Momento, Corte, Axial y Deformada

  Convenciones:
    Ejes globales X derecha, Y arriba, giro antihorario (+).
    Eje local x de nudo i a nudo j; y local a 90 grados antihorario.
    N (+) traccion. V (+) segun convencion viga (izq. sube).
    M (+) traccion abajo (fibra -y local). Diagrama M dibujado del lado
    traccionado.
============================================================================]]

local abs, sqrt, floor, max, min = math.abs, math.sqrt, math.floor, math.max, math.min
local PI = math.pi

local function fail(msg) error(msg, 0) end

------------------------------------------------------------------------------
-- 1. EXPRESIONES LINEALES CON LETRAS
--    Un valor "Lin" es una tabla {[""]=constante, P=coef, w=coef, ...}
------------------------------------------------------------------------------
local function lin_const(c) return {[""] = c} end

local function lin_isconst(a)
  for k, v in pairs(a) do if k ~= "" and v ~= 0 then return false end end
  return true
end

local function lin_add(a, b, sb)
  sb = sb or 1
  local r = {}
  for k, v in pairs(a) do r[k] = v end
  for k, v in pairs(b) do r[k] = (r[k] or 0) + sb * v end
  return r
end

local function lin_scale(a, s)
  local r = {}
  for k, v in pairs(a) do r[k] = v * s end
  return r
end

local function lin_copy(a) return lin_scale(a, 1) end

local FUNCS = {
  sqrt = function(x) if x < 0 then fail("sqrt de negativo") end return sqrt(x) end,
  abs = abs,
  sin = function(x) return math.sin(x * PI / 180) end,  -- grados
  cos = function(x) return math.cos(x * PI / 180) end,
  tan = function(x) return math.tan(x * PI / 180) end,
  ln = function(x) return math.log(x) end,
  exp = math.exp,
}

local REPL = {
  {"\226\136\146", "-"},  -- signo negativo de la TI
  {"\195\151", "*"},      -- x de multiplicar
  {"\194\183", "*"},      -- punto medio
  {"\226\128\162", "*"},  -- bullet
  {"\195\183", "/"},      -- division
  {"\225\180\135", "e"},  -- tecla EE
  {"\207\128", " pi "},   -- pi
  {"\226\136\154", " sqrt "}, -- raiz
  {"\194\178", "^2"},     -- superindice 2
  {"\194\179", "^3"},
  {"\194\176", ""},       -- grados
}

local function normalize(s)
  for _, r in ipairs(REPL) do
    local a, b = r[1], r[2]
    local i = 1
    while true do
      local p = s:find(a, i, true)
      if not p then break end
      s = s:sub(1, p - 1) .. b .. s:sub(p + #a)
      i = p + #b
    end
  end
  return s
end

local function tokenize(s)
  s = normalize(s)
  local toks, i, n = {}, 1, #s
  while i <= n do
    local ch = s:sub(i, i)
    if ch:match("%s") then
      i = i + 1
    elseif ch:match("[%d%.]") then
      local num = s:match("^%d*%.?%d*", i)
      local e = s:match("^[eE][%+%-]?%d+", i + #num)
      if e then num = num .. e end
      local v = tonumber(num)
      if not v then fail("numero invalido '" .. num .. "'") end
      toks[#toks + 1] = {t = "n", v = v}
      i = i + #num
    elseif ch:match("[%a_]") then
      local id = s:match("^[%a_][%w_]*", i)
      toks[#toks + 1] = {t = "i", v = id}
      i = i + #id
    elseif ch:match("[%+%-%*/%^%(%)]") then
      toks[#toks + 1] = {t = ch}
      i = i + 1
    else
      fail("caracter invalido '" .. ch .. "'")
    end
  end
  toks[#toks + 1] = {t = "eof"}
  return toks
end

-- resolve(nombre) -> Lin
local function parseLin(str, resolve)
  local toks = tokenize(str or "")
  local p = 1
  local function peek() return toks[p] end
  local function nxt() p = p + 1; return toks[p - 1] end
  local expr, unary, power

  local function atom()
    local t = nxt()
    if t.t == "n" then
      return lin_const(t.v)
    elseif t.t == "(" then
      local v = expr()
      if nxt().t ~= ")" then fail("falta ')'") end
      return v
    elseif t.t == "i" then
      local f = FUNCS[t.v:lower()]
      if f then
        local arg
        if peek().t == "(" then
          nxt(); arg = expr()
          if nxt().t ~= ")" then fail("falta ')'") end
        else
          arg = power()
        end
        if not lin_isconst(arg) then fail(t.v .. "() de una letra no es lineal") end
        return lin_const(f(arg[""] or 0))
      end
      if t.v == "pi" then return lin_const(PI) end
      return lin_copy(resolve(t.v))
    elseif t.t == "eof" then
      fail("expresion incompleta")
    else
      fail("simbolo inesperado '" .. t.t .. "'")
    end
  end

  power = function()
    local b = atom()
    if peek().t == "^" then
      nxt()
      local e = unary()
      if not lin_isconst(e) then fail("exponente con letra") end
      local ev = e[""] or 0
      if lin_isconst(b) then return lin_const((b[""] or 0) ^ ev) end
      if ev == 1 then return b end
      if ev == 0 then return lin_const(1) end
      fail("potencia de una letra (no lineal)")
    end
    return b
  end

  unary = function()
    local t = peek()
    if t.t == "-" then nxt(); return lin_scale(unary(), -1) end
    if t.t == "+" then nxt(); return unary() end
    return power()
  end

  local function mul(a, b)
    if lin_isconst(a) then return lin_scale(b, a[""] or 0) end
    if lin_isconst(b) then return lin_scale(a, b[""] or 0) end
    fail("producto de letras (no lineal)")
  end

  local function term()
    local v = unary()
    while true do
      local t = peek()
      if t.t == "*" then
        nxt(); v = mul(v, unary())
      elseif t.t == "/" then
        nxt()
        local d = unary()
        if not lin_isconst(d) then fail("division por una letra (no lineal)") end
        if (d[""] or 0) == 0 then fail("division por cero") end
        v = lin_scale(v, 1 / d[""])
      elseif t.t == "n" or t.t == "i" or t.t == "(" then
        v = mul(v, unary())            -- multiplicacion implicita: 2P, 3(w+1)
      else
        break
      end
    end
    return v
  end

  expr = function()
    local v = term()
    while true do
      local t = peek()
      if t.t == "+" then nxt(); v = lin_add(v, term())
      elseif t.t == "-" then nxt(); v = lin_add(v, term(), -1)
      else break end
    end
    return v
  end

  if peek().t == "eof" then return lin_const(0) end
  local v = expr()
  if peek().t ~= "eof" then fail("sobra texto en la expresion") end
  return v
end

local function evalConst(str, resolve)
  local v = parseLin(str, resolve)
  if not lin_isconst(v) then fail("se requiere un numero (sin letras)") end
  return v[""] or 0
end

local function isBlank(s) return s == nil or not tostring(s):match("%S") end

-- formato numerico
local function fmt(v, dig)
  if v ~= v then return "NaN" end
  local a = abs(v)
  if a < 1e-12 then return "0" end
  local s
  if a >= 1e6 or a < 1e-3 then
    s = string.format("%." .. (dig or 4) .. "e", v)
    s = s:gsub("%.?0+e", "e")
  else
    s = string.format("%." .. (dig or 4) .. "f", v)
    s = s:gsub("0+$", "")
    s = s:gsub("%.$", "")
  end
  if s == "-0" then s = "0" end
  return s
end

local function fmtShort(v)
  if abs(v) < 1e-12 then return "0" end
  local s = string.format("%.3g", v)
  return s
end

local function linKeys(L)
  local ks = {}
  for k, v in pairs(L) do if k ~= "" and v ~= 0 then ks[#ks + 1] = k end end
  table.sort(ks)
  return ks
end

local function linStr(L)
  local parts = {}
  local c = L[""] or 0
  if c ~= 0 then parts[#parts + 1] = {c, ""} end
  for _, k in ipairs(linKeys(L)) do parts[#parts + 1] = {L[k], k} end
  if #parts == 0 then return "0" end
  local out = ""
  for i, pr in ipairs(parts) do
    local v, k = pr[1], pr[2]
    local sgn = v < 0 and "-" or "+"
    local a = abs(v)
    local body
    if k == "" then
      body = fmt(a)
    else
      local fa = fmt(a)
      if fa == "1" then body = k
      elseif #k == 1 then body = fa .. k
      else body = fa .. "*" .. k end
    end
    if i == 1 then
      out = (sgn == "-" and "-" or "") .. body
    else
      out = out .. " " .. sgn .. " " .. body
    end
  end
  return out
end

------------------------------------------------------------------------------
-- 2. MODELO
------------------------------------------------------------------------------
local TIPOS_CARGA = {"Puntual", "Momento", "Distribuida", "Temperatura", "Error fabric."}
local DIRS = {"Perpend. local", "Axial local", "Global X", "Global Y",
              "Gravedad (-Y)", "Global X proy.", "Global Y proy.", "Grav. proy."}
local SUP_FRAME = {"Empotrado", "Articulado", "Rodillo libre X", "Rodillo libre Y",
                   "Guia libre X", "Guia libre Y", "Solo giro fijo", "Personalizado"}
local SUP_FRAME_R = {{true, true, true}, {true, true, false}, {false, true, false},
                     {true, false, false}, {false, true, true}, {true, false, true},
                     {false, false, true}}
local SUP_TRUSS = {"Articulado", "Rodillo libre X", "Rodillo libre Y", "Personalizado"}
local SUP_TRUSS_R = {{true, true, false}, {false, true, false}, {true, false, false}}

local function newModel(kind)
  local m = {kind = kind, nodes = {}, secs = {}, mems = {}, sups = {},
             nl = {}, ml = {}, vars = {}}
  if kind == "truss" then
    m.secs[1] = {E = "1", A = "1", I = "1", al = "1e-5", h = "1"}
  else
    m.secs[1] = {E = "1", A = "1e5", I = "1", al = "1e-5", h = "0.5"}
  end
  return m
end

------------------------------------------------------------------------------
-- 3. MOTOR DE CALCULO (METODO DE RIGIDEZ)
------------------------------------------------------------------------------
local Eng = {}

local function zeros(n, m)
  local A = {}
  for i = 1, n do
    local r = {}
    for j = 1, m do r[j] = 0 end
    A[i] = r
  end
  return A
end

-- tabla de variables definidas por el usuario
local function varTable(model)
  local vt = {}
  for _, v in ipairs(model.vars or {}) do
    local rec = {keep = v.keep}
    if not isBlank(v.val) then
      local ok, val = pcall(evalConst, v.val, function(n)
        local r = vt[n]
        if r and r.value then return lin_const(r.value) end
        fail("variable '" .. n .. "' sin valor")
      end)
      if not ok then fail("Variable " .. v.name .. ": " .. tostring(val)) end
      rec.value = val
    end
    vt[v.name] = rec
  end
  return vt
end

local function makeGres(vt)
  return function(n)
    local r = vt[n]
    if r and r.value then return lin_const(r.value) end
    fail("letra '" .. n .. "' sin valor numerico (definala en Variables)")
  end
end

local function ev(str, res, ctx, default)
  if isBlank(str) then
    if default ~= nil then return default end
    fail(ctx .. ": falta un valor")
  end
  local ok, v = pcall(evalConst, str, res)
  if not ok then fail(ctx .. ": " .. tostring(v)) end
  return v
end

function Eng.geometry(model, vt)
  vt = vt or varTable(model)
  local gres = makeGres(vt)
  local G = {nodes = {}, mems = {}, gres = gres, vt = vt}
  for i, nd in ipairs(model.nodes) do
    G.nodes[i] = {x = ev(nd.x, gres, "Nudo " .. i .. " X"),
                  y = ev(nd.y, gres, "Nudo " .. i .. " Y")}
  end
  local secs = {}
  for i, s in ipairs(model.secs) do
    local ctx = "Seccion " .. i
    secs[i] = {E = ev(s.E, gres, ctx .. " E"), A = ev(s.A, gres, ctx .. " A"),
               I = ev(s.I, gres, ctx .. " I", 1), al = ev(s.al, gres, ctx .. " alfa", 0),
               h = ev(s.h, gres, ctx .. " h", 1)}
    if secs[i].E <= 0 or secs[i].A <= 0 then fail(ctx .. ": E y A deben ser > 0") end
    if model.kind ~= "truss" and secs[i].I <= 0 then fail(ctx .. ": I debe ser > 0") end
    if secs[i].I <= 0 then secs[i].I = 1 end
    if secs[i].h <= 0 then secs[i].h = 1 end
  end
  for m, mb in ipairs(model.mems) do
    local ni, nj = G.nodes[mb.i], G.nodes[mb.j]
    if not ni or not nj then fail("Barra " .. m .. ": nudo inexistente") end
    if mb.i == mb.j then fail("Barra " .. m .. ": nudos iguales") end
    local s = secs[mb.sec or 1]
    if not s then fail("Barra " .. m .. ": seccion inexistente") end
    local dx, dy = nj.x - ni.x, nj.y - ni.y
    local L = sqrt(dx * dx + dy * dy)
    if L < 1e-12 then fail("Barra " .. m .. ": longitud cero") end
    local truss = model.kind == "truss"
    G.mems[m] = {i = mb.i, j = mb.j, L = L, c = dx / L, s = dy / L,
                 E = s.E, A = s.A, I = s.I, al = s.al, h = s.h,
                 ri = truss or mb.ri or false, rj = truss or mb.rj or false}
  end
  return G
end

local function kLocal(E, A, I, L)
  local a = E * A / L
  local b = 12 * E * I / L ^ 3
  local c = 6 * E * I / L ^ 2
  local d = 4 * E * I / L
  local e = 2 * E * I / L
  return {
    { a, 0, 0, -a, 0, 0},
    { 0, b, c, 0, -b, c},
    { 0, c, d, 0, -c, e},
    {-a, 0, 0, a, 0, 0},
    { 0, -b, -c, 0, b, -c},
    { 0, c, e, 0, -c, d}}
end

local function toLocal(v, c, s)
  return {c * v[1] + s * v[2], -s * v[1] + c * v[2], v[3],
          c * v[4] + s * v[5], -s * v[4] + c * v[5], v[6]}
end

local function toGlobal(v, c, s)
  return {c * v[1] - s * v[2], s * v[1] + c * v[2], v[3],
          c * v[4] - s * v[5], s * v[4] + c * v[5], v[6]}
end

-- kg = T' k T
local function kGlobal(k, c, s)
  local T = zeros(6, 6)
  for b = 0, 3, 3 do
    T[b + 1][b + 1] = c;  T[b + 1][b + 2] = s
    T[b + 2][b + 1] = -s; T[b + 2][b + 2] = c
    T[b + 3][b + 3] = 1
  end
  local kT = zeros(6, 6)
  for i = 1, 6 do for j = 1, 6 do
    local sum = 0
    for p = 1, 6 do sum = sum + k[i][p] * T[p][j] end
    kT[i][j] = sum
  end end
  local kg = zeros(6, 6)
  for i = 1, 6 do for j = 1, 6 do
    local sum = 0
    for p = 1, 6 do sum = sum + T[p][i] * kT[p][j] end
    kg[i][j] = sum
  end end
  return kg
end

-- condensacion estatica de extremos articulados
local function condense(k, rel)
  if #rel == 0 then return k, nil end
  local inv
  if #rel == 1 then
    local r = rel[1]
    if abs(k[r][r]) < 1e-300 then fail("rigidez nula al condensar") end
    inv = {{1 / k[r][r]}}
  else
    local a, b = k[rel[1]][rel[1]], k[rel[1]][rel[2]]
    local c, d = k[rel[2]][rel[1]], k[rel[2]][rel[2]]
    local det = a * d - b * c
    if abs(det) < 1e-300 then fail("rigidez nula al condensar") end
    inv = {{d / det, -b / det}, {-c / det, a / det}}
  end
  local isr = {}
  for _, r in ipairs(rel) do isr[r] = true end
  local kc = zeros(6, 6)
  for i = 1, 6 do
    if not isr[i] then
      for j = 1, 6 do
        if not isr[j] then
          local v = k[i][j]
          for p = 1, #rel do for q = 1, #rel do
            v = v - k[i][rel[p]] * inv[p][q] * k[rel[q]][j]
          end end
          kc[i][j] = v
        end
      end
    end
  end
  return kc, inv
end

local function shapes(x, L)
  local xi = x / L
  return 1 - xi, 1 - 3 * xi ^ 2 + 2 * xi ^ 3, L * (xi - 2 * xi ^ 2 + xi ^ 3),
         xi, 3 * xi ^ 2 - 2 * xi ^ 3, L * (-xi ^ 2 + xi ^ 3)
end

local function dshapes(x, L)
  local xi = x / L
  return (-6 * xi + 6 * xi ^ 2) / L, 1 - 4 * xi + 3 * xi ^ 2,
         (6 * xi - 6 * xi ^ 2) / L, -2 * xi + 3 * xi ^ 2
end

local GP = {{-0.7745966692414834, 5 / 9}, {0, 8 / 9}, {0.7745966692414834, 5 / 9}}

-- fuerzas de empotramiento perfecto (locales) de una carga ya evaluada
local function fixedEnd(d, mb)
  local L = mb.L
  local F = {0, 0, 0, 0, 0, 0}
  if d.t == 1 then
    local n1, n2, n3, n4, n5, n6 = shapes(d.a, L)
    F = {-d.px * n1, -d.py * n2, -d.py * n3, -d.px * n4, -d.py * n5, -d.py * n6}
  elseif d.t == 2 then
    local b2, b3, b5, b6 = dshapes(d.a, L)
    F = {0, -d.m * b2, -d.m * b3, 0, -d.m * b5, -d.m * b6}
  elseif d.t == 3 then
    local h = (d.b - d.a) / 2
    local mid = (d.a + d.b) / 2
    for _, g in ipairs(GP) do
      local s = mid + g[1] * h
      local t = (s - d.a) / (d.b - d.a)
      local px = d.px1 + (d.px2 - d.px1) * t
      local py = d.py1 + (d.py2 - d.py1) * t
      local wgt = g[2] * h
      local n1, n2, n3, n4, n5, n6 = shapes(s, L)
      F[1] = F[1] - wgt * px * n1
      F[2] = F[2] - wgt * py * n2
      F[3] = F[3] - wgt * py * n3
      F[4] = F[4] - wgt * px * n4
      F[5] = F[5] - wgt * py * n5
      F[6] = F[6] - wgt * py * n6
    end
  elseif d.t == 4 then
    local nT = mb.E * mb.A * mb.al * d.T
    local mT = mb.E * mb.I * mb.al * d.G / mb.h
    F = {nT, 0, mT, -nT, 0, -mT}
  elseif d.t == 5 then
    local n = mb.E * mb.A * d.dL / mb.L
    F = {n, 0, 0, -n, 0, 0}
  end
  return F
end

-- componentes locales unitarias de una direccion de carga
local function dirLocal(dir, c, s, proj)
  if dir == 1 then return 0, 1 end
  if dir == 2 then return 1, 0 end
  local gx, gy, f = 0, 0, 1
  if dir == 3 or dir == 6 then
    gx = 1; if dir == 6 and proj then f = abs(s) end
  elseif dir == 4 or dir == 7 then
    gy = 1; if dir == 7 and proj then f = abs(c) end
  else
    gy = -1; if dir == 8 and proj then f = abs(c) end
  end
  return f * (gx * c + gy * s), f * (-gx * s + gy * c)
end
Eng.dirLocal = dirLocal

-- Gauss con pivoteo parcial y varios lados derechos
local function solveMulti(A, B, n, m, tol)
  for col = 1, n do
    local piv, best = col, abs(A[col][col])
    for r = col + 1, n do
      local v = abs(A[r][col])
      if v > best then best = v; piv = r end
    end
    if best <= tol then return nil, col end
    if piv ~= col then
      A[piv], A[col] = A[col], A[piv]
      B[piv], B[col] = B[col], B[piv]
    end
    local pr, pb = A[col], B[col]
    local d = pr[col]
    for r = col + 1, n do
      local row = A[r]
      local f = row[col] / d
      if f ~= 0 then
        for c = col, n do row[c] = row[c] - f * pr[c] end
        local br = B[r]
        for k = 1, m do br[k] = br[k] - f * pb[k] end
      end
    end
  end
  local X = zeros(n, m)
  for k = 1, m do
    for i = n, 1, -1 do
      local s = B[i][k]
      local row = A[i]
      for j = i + 1, n do s = s - row[j] * X[j][k] end
      X[i][k] = s / row[i]
    end
  end
  return X
end

function Eng.solve(model)
  local vt = varTable(model)
  local G = Eng.geometry(model, vt)
  local nodes, mems, gres = G.nodes, G.mems, G.gres
  local nn = #nodes
  local truss = model.kind == "truss"
  if nn < 2 then fail("Defina al menos 2 nudos") end
  if #mems < 1 then fail("Defina al menos 1 barra") end
  if #model.sups < 1 then fail("Defina los apoyos") end

  local lres = function(n)
    local r = vt[n]
    if r and r.value and not r.keep then return lin_const(r.value) end
    local t = {}; t[n] = 1
    return t
  end
  local function evL(str, ctx)
    if isBlank(str) then return lin_const(0) end
    local ok, v = pcall(parseLin, str, lres)
    if not ok then fail(ctx .. ": " .. tostring(v)) end
    return v
  end

  -- cargas nodales
  local NL = {}
  for i, l in ipairs(model.nl) do
    if not nodes[l.node] then fail("Carga nodal " .. i .. ": nudo inexistente") end
    local ctx = "Carga nodal " .. i
    NL[#NL + 1] = {node = l.node, fx = evL(l.fx, ctx), fy = evL(l.fy, ctx),
                   m = truss and lin_const(0) or evL(l.m, ctx)}
  end

  -- cargas en barras
  local ML = {}
  for m in ipairs(mems) do mems[m].crit = {} end
  for i, l in ipairs(model.ml) do
    local mb = mems[l.mem]
    if not mb then fail("Carga en barra " .. i .. ": barra inexistente") end
    local ctx = "Carga barra " .. i
    local pres = function(n) if n == "L" then return lin_const(mb.L) end return gres(n) end
    local rec = {mem = l.mem, t = l.t, dir = l.dir or 1}
    rec.v1 = evL(l.v1, ctx)
    if l.t == 3 and isBlank(l.v2) then rec.v2 = lin_copy(rec.v1) else rec.v2 = evL(l.v2, ctx) end
    local L = mb.L
    if l.t == 1 or l.t == 2 then
      rec.a = ev(l.a, pres, ctx .. " (a)", L / 2)
      if rec.a < -1e-9 * L or rec.a > L * (1 + 1e-9) then fail(ctx .. ": a fuera de la barra (0..L)") end
      rec.a = max(0, min(L, rec.a))
      if rec.a > 0 and rec.a < L then table.insert(mb.crit, rec.a) end
      rec.ux, rec.uy = dirLocal(rec.dir, mb.c, mb.s, false)
    elseif l.t == 3 then
      rec.a = ev(l.a, pres, ctx .. " (a)", 0)
      rec.b = ev(l.b, pres, ctx .. " (b)", L)
      rec.a = max(0, rec.a); rec.b = min(L, rec.b)
      if rec.b - rec.a <= 1e-12 * L then fail(ctx .. ": se requiere 0 <= a < b <= L") end
      if rec.a > 0 then table.insert(mb.crit, rec.a) end
      if rec.b < L then table.insert(mb.crit, rec.b) end
      rec.ux, rec.uy = dirLocal(rec.dir, mb.c, mb.s, true)
    end
    ML[#ML + 1] = rec
  end

  -- apoyos
  local sup = {}
  local lin0 = lin_const(0)
  for i, s in ipairs(model.sups) do
    if not nodes[s.node] then fail("Apoyo " .. i .. ": nudo inexistente") end
    local ctx = "Apoyo " .. i
    local S = sup[s.node] or {rx = false, ry = false, rz = false, kx = 0, ky = 0, kr = 0,
                              dx = lin0, dy = lin0, dr = lin0, ang = 0}
    S.rx = S.rx or s.rx or false
    S.ry = S.ry or s.ry or false
    S.rz = S.rz or s.rz or false
    S.kx = S.kx + ev(s.kx, gres, ctx .. " kx", 0)
    S.ky = S.ky + ev(s.ky, gres, ctx .. " ky", 0)
    S.kr = S.kr + ev(s.kr, gres, ctx .. " kgiro", 0)
    S.dx = lin_add(S.dx, evL(s.dx, ctx))
    S.dy = lin_add(S.dy, evL(s.dy, ctx))
    S.dr = lin_add(S.dr, evL(s.dr, ctx))
    local a = ev(s.ang, gres, ctx .. " angulo", 0)
    if a ~= 0 then S.ang = a end
    if truss then S.rz = false; S.kr = 0; S.dr = lin0 end
    S.ca, S.sa = math.cos(S.ang * PI / 180), math.sin(S.ang * PI / 180)
    sup[s.node] = S
  end
  for m, mb in ipairs(mems) do table.sort(mb.crit) end

  -- letras presentes en las cargas
  local symset = {}
  local function addSyms(L) for k, v in pairs(L) do if k ~= "" and v ~= 0 then symset[k] = true end end end
  for _, l in ipairs(NL) do addSyms(l.fx); addSyms(l.fy); addSyms(l.m) end
  for _, l in ipairs(ML) do addSyms(l.v1); addSyms(l.v2) end
  for _, S in pairs(sup) do addSyms(S.dx); addSyms(S.dy); addSyms(S.dr) end
  local syms = {}
  for k in pairs(symset) do syms[#syms + 1] = k end
  table.sort(syms)
  local ns = #syms
  local function lv(L, k)
    if k == 0 then return L[""] or 0 end
    return L[syms[k]] or 0
  end

  -- ensamblaje de la matriz de rigidez
  local nd = 3 * nn
  local K = zeros(nd, nd)
  for m, mb in ipairs(mems) do
    mb.k = kLocal(mb.E, mb.A, mb.I, mb.L)
    mb.rel = {}
    if mb.ri then mb.rel[#mb.rel + 1] = 3 end
    if mb.rj then mb.rel[#mb.rel + 1] = 6 end
    mb.kc, mb.inv = condense(mb.k, mb.rel)
    local kg = kGlobal(mb.kc, mb.c, mb.s)
    local i, j = mb.i, mb.j
    mb.dofs = {3 * i - 2, 3 * i - 1, 3 * i, 3 * j - 2, 3 * j - 1, 3 * j}
    for a = 1, 6 do
      local ra = K[mb.dofs[a]]
      for b = 1, 6 do ra[mb.dofs[b]] = ra[mb.dofs[b]] + kg[a][b] end
    end
  end

  -- rotacion de nudos con apoyo inclinado: K' = R' K R
  local rotNodes = {}
  for n, S in pairs(sup) do if S.ang ~= 0 then rotNodes[#rotNodes + 1] = n end end
  local function rotRows(M, n, c, s, ncol)
    local p, q = 3 * n - 2, 3 * n - 1
    for j = 1, ncol do
      local a, b = M[p][j], M[q][j]
      M[p][j] = c * a + s * b
      M[q][j] = -s * a + c * b
    end
  end
  for _, n in ipairs(rotNodes) do
    local S = sup[n]
    rotRows(K, n, S.ca, S.sa, nd)
    local p, q = 3 * n - 2, 3 * n - 1
    for i = 1, nd do
      local a, b = K[i][p], K[i][q]
      K[i][p] = S.ca * a + S.sa * b
      K[i][q] = -S.sa * a + S.ca * b
    end
  end

  -- resortes y restricciones
  local known = {}   -- dof -> Lin del desplazamiento impuesto
  for n, S in pairs(sup) do
    local b = 3 * n - 3
    K[b + 1][b + 1] = K[b + 1][b + 1] + S.kx
    K[b + 2][b + 2] = K[b + 2][b + 2] + S.ky
    K[b + 3][b + 3] = K[b + 3][b + 3] + S.kr
    if S.rx then known[b + 1] = S.dx end
    if S.ry then known[b + 2] = S.dy end
    if S.rz then known[b + 3] = S.dr end
  end
  local maxd = 0
  for i = 1, nd do maxd = max(maxd, abs(K[i][i])) end
  if maxd == 0 then fail("Matriz de rigidez nula") end
  local warn = {}
  local autofix = {}
  for i = 1, nd do
    if not known[i] and abs(K[i][i]) <= 1e-12 * maxd then
      local n = floor((i - 1) / 3) + 1
      if i % 3 == 0 then
        known[i] = lin0; autofix[n] = true
      else
        fail("Nudo " .. n .. " sin rigidez en " .. (i % 3 == 1 and "X" or "Y") ..
             " (nudo suelto o mecanismo)")
      end
    end
  end
  for _, l in ipairs(NL) do
    if autofix[l.node] and not (lin_isconst(l.m) and (l.m[""] or 0) == 0) then
      warn[#warn + 1] = "Momento en nudo " .. l.node .. " (todo articulado) se ignora"
    end
  end

  local free, fidx = {}, {}
  for i = 1, nd do if not known[i] then free[#free + 1] = i; fidx[i] = #free end end
  local nf = #free

  -- vectores de carga por caso (caso 0 = constantes, caso k = coef. de letra k)
  local cases = {}
  local B = zeros(nf, ns + 1)
  for k = 0, ns do
    local C = {F = {}, D = {}, L = {}, R = {}, d = {}, F0 = {}, P = {}}
    cases[k] = C
    local Fv = {}
    for i = 1, nd do Fv[i] = 0 end
    for n = 1, nn do C.P[n] = {0, 0, 0} end
    for _, l in ipairs(NL) do
      local P = C.P[l.node]
      P[1] = P[1] + lv(l.fx, k); P[2] = P[2] + lv(l.fy, k); P[3] = P[3] + lv(l.m, k)
    end
    for n = 1, nn do
      Fv[3 * n - 2] = C.P[n][1]; Fv[3 * n - 1] = C.P[n][2]; Fv[3 * n] = C.P[n][3]
    end
    for m in ipairs(mems) do C.L[m] = {} end
    for _, l in ipairs(ML) do
      local mb = mems[l.mem]
      local d
      local v1, v2 = lv(l.v1, k), lv(l.v2, k)
      if l.t == 1 then d = {t = 1, a = l.a, px = v1 * l.ux, py = v1 * l.uy}
      elseif l.t == 2 then d = {t = 2, a = l.a, m = v1}
      elseif l.t == 3 then d = {t = 3, a = l.a, b = l.b, px1 = v1 * l.ux, px2 = v2 * l.ux,
                                py1 = v1 * l.uy, py2 = v2 * l.uy}
      elseif l.t == 4 then d = {t = 4, T = v1, G = v2}
      else d = {t = 5, dL = v1} end
      table.insert(C.L[l.mem], d)
    end
    for m, mb in ipairs(mems) do
      local F0 = {0, 0, 0, 0, 0, 0}
      for _, d in ipairs(C.L[m]) do
        local f = fixedEnd(d, mb)
        for a = 1, 6 do F0[a] = F0[a] + f[a] end
      end
      C.F0[m] = F0
      local F0c = F0
      if mb.inv then
        F0c = {}
        local isr = {}
        for _, r in ipairs(mb.rel) do isr[r] = true end
        for a = 1, 6 do
          if isr[a] then F0c[a] = 0 else
            local v = F0[a]
            for p = 1, #mb.rel do for q = 1, #mb.rel do
              v = v - mb.k[a][mb.rel[p]] * mb.inv[p][q] * F0[mb.rel[q]]
            end end
            F0c[a] = v
          end
        end
      end
      local fg = toGlobal(F0c, mb.c, mb.s)
      for a = 1, 6 do Fv[mb.dofs[a]] = Fv[mb.dofs[a]] - fg[a] end
    end
    for _, n in ipairs(rotNodes) do
      local S = sup[n]
      local p, q = 3 * n - 2, 3 * n - 1
      local a, b = Fv[p], Fv[q]
      Fv[p] = S.ca * a + S.sa * b
      Fv[q] = -S.sa * a + S.ca * b
    end
    local dk = {}
    for i, L in pairs(known) do dk[i] = lv(L, k) end
    C.dk = dk
    for fi, i in ipairs(free) do
      local v = Fv[i]
      local row = K[i]
      for j, dv in pairs(dk) do if dv ~= 0 then v = v - row[j] * dv end end
      B[fi][k + 1] = v
    end
  end

  -- resolver
  local X = {}
  if nf > 0 then
    local A = zeros(nf, nf)
    for a = 1, nf do
      local ra, Ka = A[a], K[free[a]]
      for b = 1, nf do ra[b] = Ka[free[b]] end
    end
    local col
    X, col = solveMulti(A, B, nf, ns + 1, 1e-11 * maxd)
    if not X then
      local dof = free[col]
      local n = floor((dof - 1) / 3) + 1
      local comp = ({"X", "Y", "giro"})[(dof - 1) % 3 + 1]
      fail("ESTRUCTURA INESTABLE (mecanismo). Revise apoyos y articulaciones cerca del nudo " ..
           n .. " (direccion " .. comp .. ").")
    end
  end

  -- post-proceso
  local supNodes = {}
  for n in pairs(sup) do supNodes[#supNodes + 1] = n end
  table.sort(supNodes)
  for k = 0, ns do
    local C = cases[k]
    local d = {}
    for i = 1, nd do
      if fidx[i] then d[i] = X[fidx[i]][k + 1] else d[i] = C.dk[i] end
    end
    for _, n in ipairs(rotNodes) do
      local S = sup[n]
      local p, q = 3 * n - 2, 3 * n - 1
      local a, b = d[p], d[q]
      d[p] = S.ca * a - S.sa * b
      d[q] = S.sa * a + S.ca * b
    end
    C.d = d
    local fs, ds = 0, 0
    for i = 1, nd do ds = max(ds, abs(d[i])) end
    local Rn = {}
    for n = 1, nn do Rn[n] = {-C.P[n][1], -C.P[n][2], -C.P[n][3]} end
    for m, mb in ipairs(mems) do
      local dg = {}
      for a = 1, 6 do dg[a] = d[mb.dofs[a]] end
      local dl = toLocal(dg, mb.c, mb.s)
      local F0 = C.F0[m]
      if mb.inv then
        local isr = {}
        for _, r in ipairs(mb.rel) do isr[r] = true end
        local tmp = {}
        for q, r in ipairs(mb.rel) do
          local v = F0[r]
          for j = 1, 6 do if not isr[j] then v = v + mb.k[r][j] * dl[j] end end
          tmp[q] = v
        end
        for p, r in ipairs(mb.rel) do
          local v = 0
          for q = 1, #mb.rel do v = v - mb.inv[p][q] * tmp[q] end
          dl[r] = v
        end
      end
      local F = {}
      for a = 1, 6 do
        local v = F0[a]
        for b = 1, 6 do v = v + mb.k[a][b] * dl[b] end
        F[a] = v
        fs = max(fs, abs(v))
      end
      if mb.ri then F[3] = 0 end
      if mb.rj then F[6] = 0 end
      C.F[m] = F
      C.D[m] = dl
      local fg = toGlobal(F, mb.c, mb.s)
      for a = 1, 3 do
        Rn[mb.i][a] = Rn[mb.i][a] + fg[a]
        Rn[mb.j][a] = Rn[mb.j][a] + fg[a + 3]
      end
    end
    for _, n in ipairs(supNodes) do
      C.R[n] = Rn[n]
      for a = 1, 3 do fs = max(fs, abs(Rn[n][a])) end
    end
    C.fs, C.ds = fs, ds
  end

  -- valores numericos de las letras (para graficos)
  local w = {[0] = 1}
  local vals = {}
  for k, s in ipairs(syms) do
    local r = vt[s]
    local v = (r and r.value) or 1
    w[k] = v; vals[k] = v
    if not (r and r.value) then
      warn[#warn + 1] = "Letra " .. s .. " sin valor: se usa " .. s .. "=1 en graficos"
    end
  end

  return {kind = model.kind, nodes = nodes, mems = mems, syms = syms, ns = ns,
          cases = cases, sup = sup, supNodes = supNodes, w = w, vals = vals,
          warn = warn, NL = NL, ML = ML, autofix = autofix}
end

-- esfuerzos internos de un caso en la seccion x
local function internalCase(C, m, x, right)
  local F = C.F[m]
  local N, V, M = -F[1], F[2], -F[3] + F[2] * x
  for _, d in ipairs(C.L[m]) do
    if d.t == 1 then
      if d.a < x or (right and d.a == x) then
        N = N - d.px; V = V + d.py; M = M + d.py * (x - d.a)
      end
    elseif d.t == 2 then
      if d.a < x or (right and d.a == x) then M = M - d.m end
    elseif d.t == 3 then
      local c = min(x, d.b)
      if c > d.a then
        local u, X = c - d.a, x - d.a
        local len = d.b - d.a
        local kx = (d.px2 - d.px1) / len
        local ky = (d.py2 - d.py1) / len
        N = N - (d.px1 * u + kx * u * u / 2)
        V = V + (d.py1 * u + ky * u * u / 2)
        M = M + d.py1 * (X * u - u * u / 2) + ky * (X * u * u / 2 - u ^ 3 / 3)
      end
    end
  end
  return N, V, M
end

function Eng.internal(res, m, x, right, w)
  w = w or res.w
  local N, V, M = 0, 0, 0
  for k = 0, res.ns do
    local wk = w[k]
    if wk ~= 0 then
      local n, v, mm = internalCase(res.cases[k], m, x, right)
      N = N + wk * n; V = V + wk * v; M = M + wk * mm
    end
  end
  return N, V, M
end

local function mkLin(res, f, sk)
  local L = {}
  for k = 0, res.ns do
    local v = f(k)
    local sc = res.cases[k][sk] or 0
    if abs(v) > 1e-9 * sc then
      L[k == 0 and "" or res.syms[k]] = v
    end
  end
  return L
end

function Eng.internalLin(res, m, x, right)
  local cache = {}
  for k = 0, res.ns do cache[k] = {internalCase(res.cases[k], m, x, right)} end
  return mkLin(res, function(k) return cache[k][1] end, "fs"),
         mkLin(res, function(k) return cache[k][2] end, "fs"),
         mkLin(res, function(k) return cache[k][3] end, "fs")
end

function Eng.linNum(res, L)
  local v = L[""] or 0
  for k, s in ipairs(res.syms) do v = v + (L[s] or 0) * res.w[k] end
  return v
end

function Eng.dispLin(res, n, a)
  return mkLin(res, function(k) return res.cases[k].d[3 * n - 3 + a] end, "ds")
end

function Eng.reacLin(res, n, a)
  return mkLin(res, function(k) return res.cases[k].R[n][a] end, "fs")
end

-- puntos de muestreo (x, lado derecho) incluyendo discontinuidades
function Eng.samples(res, m, n)
  local mb = res.mems[m]
  local L = mb.L
  local pts = {}
  for i = 0, n do pts[#pts + 1] = {L * i / n, false} end
  for _, a in ipairs(mb.crit) do
    pts[#pts + 1] = {a, false}; pts[#pts + 1] = {a, true}
  end
  table.sort(pts, function(p, q)
    if p[1] ~= q[1] then return p[1] < q[1] end
    return (not p[2]) and q[2]
  end)
  return pts
end

-- deformada numerica: devuelve lista {x, ux_global, uy_global}
function Eng.deflect(res, m, nseg)
  local mb = res.mems[m]
  local L = mb.L
  local w = res.w
  local dl = {0, 0, 0, 0, 0, 0}
  local k0 = 0
  for k = 0, res.ns do
    local D = res.cases[k].D[m]
    for a = 1, 6 do dl[a] = dl[a] + w[k] * D[a] end
    for _, d in ipairs(res.cases[k].L[m]) do
      if d.t == 4 then k0 = k0 + w[k] * mb.al * d.G / mb.h end
    end
  end
  local EI = mb.E * mb.I
  local h = L / nseg
  local vs, xs = {dl[2]}, {0}
  local th, v = dl[3], dl[2]
  local _, _, Mp = Eng.internal(res, m, 0, true)
  local kp = Mp / EI + k0
  if res.kind == "truss" and #res.cases[0].L[m] == 0 then kp = 0 end
  for i = 1, nseg do
    local x = i * h
    local _, _, M = Eng.internal(res, m, x, false)
    local kc = M / EI + k0
    local thn = th + (kp + kc) / 2 * h
    v = v + (th + thn) / 2 * h
    th = thn; kp = kc
    vs[#vs + 1] = v; xs[#xs + 1] = x
  end
  local err = dl[5] - v
  local out = {}
  for i = 1, #xs do
    local x = xs[i]
    local vv = vs[i] + err * x / L
    local uu = dl[1] + (dl[4] - dl[1]) * x / L
    out[i] = {x, mb.c * uu - mb.s * vv, mb.s * uu + mb.c * vv}
  end
  return out
end

-- extremos de N, V, M en una barra (numericos)
function Eng.extremes(res, m)
  local pts = Eng.samples(res, m, 120)
  local e = {}
  for q = 1, 3 do e[q] = {mx = -1e300, xmx = 0, mn = 1e300, xmn = 0} end
  for _, p in ipairs(pts) do
    local r = {Eng.internal(res, m, p[1], p[2])}
    for q = 1, 3 do
      if r[q] > e[q].mx then e[q].mx = r[q]; e[q].xmx = p[1] end
      if r[q] < e[q].mn then e[q].mn = r[q]; e[q].xmn = p[1] end
    end
  end
  return e
end

-- resultantes de las cargas aplicadas (para verificar equilibrio) caso k
function Eng.appliedTotals(res, k)
  local C = res.cases[k]
  local fx, fy, mo = 0, 0, 0
  for n, P in ipairs(C.P) do
    local nd = res.nodes[n]
    fx = fx + P[1]; fy = fy + P[2]
    mo = mo + P[3] + nd.x * P[2] - nd.y * P[1]
  end
  for m, mb in ipairs(res.mems) do
    local ni = res.nodes[mb.i]
    local function addF(s, px, py)
      local gx, gy = mb.c * px - mb.s * py, mb.s * px + mb.c * py
      local x, y = ni.x + s * mb.c, ni.y + s * mb.s
      fx = fx + gx; fy = fy + gy
      mo = mo + x * gy - y * gx
    end
    for _, d in ipairs(C.L[m]) do
      if d.t == 1 then addF(d.a, d.px, d.py)
      elseif d.t == 2 then mo = mo + d.m
      elseif d.t == 3 then
        local h = (d.b - d.a) / 2
        local mid = (d.a + d.b) / 2
        for _, g in ipairs(GP) do
          local s = mid + g[1] * h
          local t = (s - d.a) / (d.b - d.a)
          addF(s, g[2] * h * (d.px1 + (d.px2 - d.px1) * t), g[2] * h * (d.py1 + (d.py2 - d.py1) * t))
        end
      end
    end
  end
  local rx, ry, rm = 0, 0, 0
  for _, n in ipairs(res.supNodes) do
    local R = C.R[n]
    local nd = res.nodes[n]
    rx = rx + R[1]; ry = ry + R[2]
    rm = rm + R[3] + nd.x * R[2] - nd.y * R[1]
  end
  return fx, fy, mo, rx, ry, rm
end

------------------------------------------------------------------------------
-- 4. EJEMPLOS
------------------------------------------------------------------------------
local Examples = {}

Examples[1] = {"Viga continua 3 apoyos (w, P)", function()
  local m = newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = "6", y = "0"}, {x = "10", y = "0"}}
  m.mems = {{i = 1, j = 2, sec = 1}, {i = 2, j = 3, sec = 1}}
  m.sups = {{node = 1, tipo = 2, rx = true, ry = true, rz = false},
            {node = 2, tipo = 3, rx = false, ry = true, rz = false},
            {node = 3, tipo = 1, rx = true, ry = true, rz = true}}
  m.ml = {{mem = 1, t = 3, dir = 5, v1 = "w", v2 = "", a = "0", b = "L"},
          {mem = 2, t = 1, dir = 5, v1 = "P", a = "L/2"}}
  m.vars = {{name = "w", val = "10", keep = true}, {name = "P", val = "20", keep = true}}
  return m
end}

Examples[2] = {"Portico empotrado-articulado", function()
  local m = newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = "0", y = "4"}, {x = "6", y = "4"}, {x = "6", y = "0"}}
  m.mems = {{i = 1, j = 2, sec = 1}, {i = 2, j = 3, sec = 1}, {i = 3, j = 4, sec = 1}}
  m.sups = {{node = 1, tipo = 1, rx = true, ry = true, rz = true},
            {node = 4, tipo = 2, rx = true, ry = true, rz = false}}
  m.nl = {{node = 2, fx = "P", fy = "", m = ""}}
  m.ml = {{mem = 2, t = 3, dir = 5, v1 = "w", v2 = "", a = "0", b = "L"}}
  return m
end}

Examples[3] = {"Portico 2 aguas con rotula", function()
  local m = newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = "0", y = "5"}, {x = "5", y = "7"},
             {x = "10", y = "5"}, {x = "10", y = "0"}}
  m.mems = {{i = 1, j = 2, sec = 1}, {i = 2, j = 3, sec = 1, rj = true},
            {i = 3, j = 4, sec = 1}, {i = 4, j = 5, sec = 1}}
  m.sups = {{node = 1, tipo = 1, rx = true, ry = true, rz = true},
            {node = 5, tipo = 1, rx = true, ry = true, rz = true}}
  m.ml = {{mem = 2, t = 3, dir = 8, v1 = "q", v2 = "", a = "0", b = "L"},
          {mem = 3, t = 3, dir = 8, v1 = "q", v2 = "", a = "0", b = "L"},
          {mem = 1, t = 3, dir = 3, v1 = "0", v2 = "h", a = "0", b = "L"}}
  return m
end}

Examples[4] = {"Enrejado hiperestatico", function()
  local m = newModel("truss")
  m.nodes = {{x = "0", y = "0"}, {x = "3", y = "0"}, {x = "6", y = "0"},
             {x = "0", y = "3"}, {x = "3", y = "3"}, {x = "6", y = "3"}}
  m.mems = {{i = 1, j = 2, sec = 1}, {i = 2, j = 3, sec = 1}, {i = 4, j = 5, sec = 1},
            {i = 5, j = 6, sec = 1}, {i = 1, j = 4, sec = 1}, {i = 2, j = 5, sec = 1},
            {i = 3, j = 6, sec = 1}, {i = 1, j = 5, sec = 1}, {i = 2, j = 4, sec = 1},
            {i = 2, j = 6, sec = 1}, {i = 3, j = 5, sec = 1}}
  m.sups = {{node = 1, tipo = 1, rx = true, ry = true},
            {node = 3, tipo = 1, rx = true, ry = true}}
  m.nl = {{node = 5, fx = "", fy = "-P"}, {node = 4, fx = "H", fy = ""}}
  return m
end}

Examples[5] = {"Armadura isostatica (Pratt)", function()
  local m = newModel("truss")
  m.nodes = {{x = "0", y = "0"}, {x = "4", y = "0"}, {x = "8", y = "0"}, {x = "12", y = "0"},
             {x = "4", y = "3"}, {x = "8", y = "3"}}
  m.mems = {{i = 1, j = 2}, {i = 2, j = 3}, {i = 3, j = 4}, {i = 1, j = 5}, {i = 5, j = 6},
            {i = 6, j = 4}, {i = 2, j = 5}, {i = 3, j = 6}, {i = 2, j = 6}}
  for _, mb in ipairs(m.mems) do mb.sec = 1 end
  m.sups = {{node = 1, tipo = 1, rx = true, ry = true},
            {node = 4, tipo = 2, rx = false, ry = true}}
  m.nl = {{node = 2, fx = "", fy = "-P"}, {node = 3, fx = "", fy = "-P"}}
  return m
end}

------------------------------------------------------------------------------
-- 5. INTERFAZ GRAFICA
------------------------------------------------------------------------------
local W, H = 318, 212
local stack = {}
local App = {model = nil, res = nil}

local function inval()
  if platform and platform.window then platform.window:invalidate() end
end
local function top() return stack[#stack] end
local function push(s) stack[#stack + 1] = s; inval() end
local function pop() if #stack > 1 then stack[#stack] = nil end; inval() end

local C = {
  title = {14, 38, 72}, title2 = {32, 84, 150}, titleTx = {255, 255, 255},
  accent = {242, 138, 36}, bg = {243, 246, 250}, card = {255, 255, 255},
  border = {200, 209, 222}, grid = {228, 234, 243}, sel = {218, 232, 252}, selBar = {32, 104, 214},
  tx = {22, 28, 38}, dim = {110, 118, 130}, err = {200, 30, 30}, foot = {232, 236, 242},
  ok = {20, 140, 70}, row = {248, 250, 253},
  M = {210, 36, 36}, V = {20, 140, 45}, N = {36, 84, 214}, D = {150, 50, 170},
  mem = {45, 52, 64}, light = {196, 204, 216}, load = {232, 118, 0}, sup = {0, 118, 128},
}
local function col(gc, c) gc:setColorRGB(c[1], c[2], c[3]) end
local function font(gc, sz, st) gc:setFont("sansserif", st or "r", sz) end

-- degradado vertical por franjas
local function grad(gc, x, y, w, h, c1, c2)
  local step = 2
  for i = 0, h - 1, step do
    local t = h > 1 and i / (h - 1) or 0
    gc:setColorRGB(floor(c1[1] + (c2[1] - c1[1]) * t), floor(c1[2] + (c2[2] - c1[2]) * t),
                   floor(c1[3] + (c2[3] - c1[3]) * t))
    gc:fillRect(x, y + i, w, min(step, h - i))
  end
end

-- rectangulo redondeado relleno
local function fillRound(gc, x, y, w, h, r)
  r = min(r, floor(h / 2), floor(w / 2))
  if r < 1 then gc:fillRect(x, y, w, h); return end
  gc:fillRect(x + r, y, w - 2 * r, h)
  gc:fillRect(x, y + r, r, h - 2 * r)
  gc:fillRect(x + w - r, y + r, r, h - 2 * r)
  local d = 2 * r
  gc:fillArc(x, y, d, d, 90, 90)
  gc:fillArc(x + w - d - 1, y, d, d, 0, 90)
  gc:fillArc(x, y + h - d - 1, d, d, 180, 90)
  gc:fillArc(x + w - d - 1, y + h - d - 1, d, d, 270, 90)
end

-- contorno redondeado
local function strokeRound(gc, x, y, w, h, r)
  r = min(r, floor(h / 2), floor(w / 2))
  local d = 2 * r
  gc:drawLine(x + r, y, x + w - r, y)
  gc:drawLine(x + r, y + h, x + w - r, y + h)
  gc:drawLine(x, y + r, x, y + h - r)
  gc:drawLine(x + w, y + r, x + w, y + h - r)
  if r >= 1 then
    gc:drawArc(x, y, d, d, 90, 90)
    gc:drawArc(x + w - d, y, d, d, 0, 90)
    gc:drawArc(x, y + h - d, d, d, 180, 90)
    gc:drawArc(x + w - d, y + h - d, d, d, 270, 90)
  end
end

-- fondo tipo plano (cuadricula)
local function gridBg(gc, x, y, w, h, step)
  step = step or 12
  col(gc, C.card); gc:fillRect(x, y, w, h)
  col(gc, C.grid)
  for gx = x + step, x + w - 1, step do gc:fillRect(gx, y, 1, h) end
  for gy = y + step, y + h - 1, step do gc:fillRect(x, gy, w, 1) end
end

-- logo: pequena armadura
local function drawLogo(gc, x, y, s, c)
  col(gc, c)
  gc:setPen("thin", "smooth")
  local p = {{0, 1}, {0.5, 0}, {1, 1}}
  gc:drawLine(x, y + s, x + s, y + s)
  gc:drawLine(x, y + s, x + s / 2, y)
  gc:drawLine(x + s / 2, y, x + s, y + s)
  gc:drawLine(x + s / 4, y + s / 2, x + s / 2, y + s)
  gc:drawLine(x + s / 2, y + s, x + 3 * s / 4, y + s / 2)
  gc:drawLine(x + s / 4, y + s / 2, x + 3 * s / 4, y + s / 2)
end

local function drawHeader(gc, title, right)
  grad(gc, 0, 0, W, 17, C.title2, C.title)
  col(gc, C.accent); gc:fillRect(0, 17, W, 1)
  drawLogo(gc, 4, 3, 11, C.accent)
  col(gc, C.titleTx); font(gc, 10, "b")
  gc:drawString(title, 20, 1, "top")
  if right then
    font(gc, 9)
    local w = gc:getStringWidth(right)
    col(gc, {190, 208, 235})
    gc:drawString(right, W - w - 5, 2, "top")
  end
end

-- pie con "chips" de teclas: "enter:abrir|esc:volver"
local function drawFooter(gc, text, err)
  grad(gc, 0, H - 14, W, 14, {240, 243, 248}, {222, 228, 237})
  col(gc, C.border); gc:fillRect(0, H - 14, W, 1)
  text = text or ""
  font(gc, 7)
  if err or not text:find(":") then
    col(gc, err and C.err or C.dim)
    gc:drawString(text, 4, H - 13, "top")
    return
  end
  local x = 3
  for chunk in text:gmatch("[^|]+") do
    local k, d = chunk:match("^([^:]*):(.*)$")
    if not k then k, d = "", chunk end
    local kw = gc:getStringWidth(k) + 6
    local dw = gc:getStringWidth(d)
    if x + kw + dw + 4 > W then break end
    if k ~= "" then
      col(gc, C.title2); fillRound(gc, x, H - 12, kw, 11, 3)
      col(gc, C.titleTx); gc:drawString(k, x + 3, H - 13, "top")
      x = x + kw + 2
    end
    col(gc, C.tx); gc:drawString(d, x, H - 13, "top")
    x = x + dw + 7
  end
end

-- iconos de 12x12 para los menus
local function drawIcon(gc, name, x, y, c)
  col(gc, c)
  gc:setPen("thin", "smooth")
  if name == "node" then
    gc:fillArc(x + 3, y + 3, 7, 7, 0, 360)
    gc:drawLine(x, y + 6, x + 12, y + 6); gc:drawLine(x + 6, y, x + 6, y + 12)
  elseif name == "section" then
    gc:fillRect(x + 1, y + 1, 11, 3); gc:fillRect(x + 5, y + 3, 3, 7); gc:fillRect(x + 1, y + 9, 11, 3)
  elseif name == "bar" then
    gc:setPen("medium", "smooth"); gc:drawLine(x + 2, y + 10, x + 10, y + 2); gc:setPen("thin", "smooth")
    gc:fillArc(x, y + 8, 4, 4, 0, 360); gc:fillArc(x + 8, y, 4, 4, 0, 360)
  elseif name == "support" then
    gc:fillPolygon({x + 6, y + 1, x + 1, y + 9, x + 11, y + 9, x + 6, y + 1})
    gc:drawLine(x, y + 11, x + 12, y + 11)
  elseif name == "nload" then
    gc:drawLine(x + 6, y, x + 6, y + 9)
    gc:fillPolygon({x + 2, y + 6, x + 10, y + 6, x + 6, y + 11, x + 2, y + 6})
  elseif name == "mload" then
    gc:drawLine(x, y + 1, x + 12, y + 1)
    for i = 0, 2 do
      local ax = x + 1 + i * 5
      gc:drawLine(ax, y + 1, ax, y + 8)
      gc:fillPolygon({ax - 2, y + 7, ax + 2, y + 7, ax, y + 10, ax - 2, y + 7})
    end
    gc:drawLine(x, y + 11, x + 12, y + 11)
  elseif name == "var" then
    font(gc, 9, "bi"); gc:drawString("x", x + 2, y - 3, "top")
  elseif name == "eye" then
    gc:drawArc(x, y + 2, 12, 8, 0, 360); gc:fillArc(x + 4, y + 4, 4, 4, 0, 360)
  elseif name == "play" then
    gc:fillPolygon({x + 2, y, x + 12, y + 6, x + 2, y + 12, x + 2, y})
  elseif name == "frame" then
    gc:setPen("medium", "smooth")
    gc:drawLine(x + 1, y + 12, x + 1, y + 3); gc:drawLine(x + 1, y + 3, x + 11, y + 3)
    gc:drawLine(x + 11, y + 3, x + 11, y + 12)
    gc:setPen("thin", "smooth")
  elseif name == "truss" then
    drawLogo(gc, x, y, 11, c)
  elseif name == "go" then
    gc:fillRect(x, y + 5, 7, 3)
    gc:fillPolygon({x + 6, y + 1, x + 12, y + 6, x + 6, y + 11, x + 6, y + 1})
  elseif name == "book" then
    gc:drawRect(x + 1, y + 1, 10, 10); gc:fillRect(x + 1, y + 1, 3, 11)
    gc:drawLine(x + 6, y + 4, x + 9, y + 4); gc:drawLine(x + 6, y + 7, x + 9, y + 7)
  elseif name == "help" then
    gc:drawArc(x, y, 12, 12, 0, 360); font(gc, 7, "b"); gc:drawString("?", x + 4, y - 1, "top")
  elseif name == "info" then
    gc:drawArc(x, y, 12, 12, 0, 360); font(gc, 7, "b"); gc:drawString("i", x + 5, y - 1, "top")
  elseif name == "diagram" then
    gc:drawLine(x, y + 3, x + 12, y + 3)
    gc:drawPolyLine({x, y + 3, x + 3, y + 9, x + 6, y + 11, x + 9, y + 9, x + 12, y + 3})
  elseif name == "reaction" then
    gc:drawLine(x + 6, y + 2, x + 6, y + 11)
    gc:fillPolygon({x + 2, y + 5, x + 10, y + 5, x + 6, y, x + 2, y + 5})
    gc:drawLine(x, y + 11, x + 12, y + 11)
  elseif name == "disp" then
    gc:drawPolyLine({x, y + 10, x + 4, y + 4, x + 8, y + 8, x + 12, y + 1})
  elseif name == "forces" then
    gc:drawLine(x, y + 6, x + 12, y + 6)
    gc:fillPolygon({x, y + 6, x + 4, y + 2, x + 4, y + 10, x, y + 6})
    gc:fillPolygon({x + 12, y + 6, x + 8, y + 2, x + 8, y + 10, x + 12, y + 6})
  elseif name == "table" then
    gc:drawRect(x, y + 1, 12, 10); gc:drawLine(x, y + 4, x + 12, y + 4)
    gc:drawLine(x + 4, y + 1, x + 4, y + 11); gc:drawLine(x + 8, y + 1, x + 8, y + 11)
  elseif name == "check" then
    gc:setPen("medium", "smooth")
    gc:drawLine(x + 1, y + 6, x + 5, y + 10); gc:drawLine(x + 5, y + 10, x + 11, y + 2)
    gc:setPen("thin", "smooth")
  elseif name == "warn" then
    gc:fillPolygon({x + 6, y, x + 12, y + 11, x, y + 11, x + 6, y})
    col(gc, C.card); gc:fillRect(x + 5, y + 4, 2, 4); gc:fillRect(x + 5, y + 9, 2, 1)
  end
end

local function fitText(gc, s, wmax)
  if gc:getStringWidth(s) <= wmax then return s end
  while #s > 1 and gc:getStringWidth(s .. "..") > wmax do s = s:sub(1, -2) end
  return s .. ".."
end

local function modelChanged() App.res = nil end

-- pantalla dividida: panel izquierdo de datos y vista previa a la derecha
local SPLIT = 150
local drawPreviewPane, cyclePreview

local function msgBox(title, lines) end  -- se define mas abajo

------------------------- Menu ------------------------------------------------
local Menu = {}; Menu.__index = Menu
function Menu.new(title, items, footer)
  return setmetatable({title = title, items = items, sel = 1, off = 0,
                       footer = footer or "enter:abrir|1-9:elegir|esc:volver"}, Menu)
end
function Menu:getItems() return type(self.items) == "function" and self.items() or self.items end
function Menu:paint(gc)
  local items = self:getItems()
  drawHeader(gc, type(self.title) == "function" and self.title() or self.title)
  if self.split then drawPreviewPane(gc, self, H - 15) end
  local W = self.split and SPLIT or W
  col(gc, C.bg); gc:fillRect(0, 18, W, H - 32)
  local rh, y0 = self.split and 17 or 19, 21
  local vis = floor((H - 16 - y0) / rh)
  if self.sel > #items then self.sel = #items end
  if self.sel < 1 then self.sel = 1 end
  if self.sel > self.off + vis then self.off = self.sel - vis end
  if self.sel <= self.off then self.off = self.sel - 1 end
  for r = 1, vis do
    local i = self.off + r
    local it = items[i]
    if not it then break end
    local y = y0 + (r - 1) * rh
    local selected = (i == self.sel)
    if selected then
      col(gc, C.sel); fillRound(gc, 3, y - 1, W - 6, rh - 2, 4)
      col(gc, C.selBar); fillRound(gc, 3, y - 1, 3, rh - 2, 1)
    end
    -- numero en circulo
    local num = i <= 9 and tostring(i) or (i == 10 and "0" or "")
    local cx = 9
    if num ~= "" then
      col(gc, selected and C.selBar or C.light)
      gc:fillArc(cx, y + (rh - 15) / 2, 12, 12, 0, 360)
      col(gc, C.titleTx); font(gc, 7, "b")
      gc:drawString(num, cx + 6 - gc:getStringWidth(num) / 2, y + (rh - 15) / 2 - 1, "top")
    end
    local tx = cx + 16
    if it.icon then
      drawIcon(gc, it.icon, tx, y + (rh - 15) / 2, selected and C.selBar or C.title2)
      tx = tx + 17
    end
    col(gc, it.dim and C.dim or C.tx)
    font(gc, self.split and 9 or 10, selected and "b" or "r")
    local lab = type(it.label) == "function" and it.label() or it.label
    gc:drawString(fitText(gc, lab, W - tx - 6), tx, y + (rh - 17) / 2, "top")
  end
  drawFooter(gc, self.footer)
end
function Menu:arrow(k)
  local n = #self:getItems()
  if k == "up" then self.sel = self.sel > 1 and self.sel - 1 or n
  elseif k == "down" then self.sel = self.sel < n and self.sel + 1 or 1 end
end
function Menu:enter()
  local it = self:getItems()[self.sel]
  if it and it.action then it.action() end
end
function Menu:char(ch)
  if self.split and (ch == "p" or ch == "P") then cyclePreview(); return end
  local n = tonumber(ch)
  local items = self:getItems()
  if n then
    if n == 0 then n = 10 end
    if n >= 1 and n <= #items then self.sel = n; self:enter() end
  end
end
function Menu:esc() pop() end
function Menu:click(x, y)
  if self.split and x > SPLIT then return end
  local i = self.off + floor((y - 21) / (self.split and 17 or 19)) + 1
  if i >= 1 and i <= #self:getItems() then
    if i == self.sel then self:enter() else self.sel = i end
  end
end

------------------------- Lista de elementos ---------------------------------
local List = {}; List.__index = List
function List.new(o)
  o.sel = 1; o.off = 0
  return setmetatable(o, List)
end
function List:paint(gc)
  local items = self.items()
  drawHeader(gc, self.title, #items .. " elem.")
  drawPreviewPane(gc, self, H - 15)
  local W = SPLIT
  col(gc, C.bg); gc:fillRect(0, 18, W, H - 32)
  local rh, y0 = 15, 21
  local n = #items + 1
  local vis = floor((H - 16 - y0) / rh)
  if self.sel > n then self.sel = n end
  if self.sel > self.off + vis then self.off = self.sel - vis end
  if self.sel <= self.off then self.off = self.sel - 1 end
  for r = 1, vis do
    local i = self.off + r
    if i > n then break end
    local y = y0 + (r - 1) * rh
    if i == self.sel then
      col(gc, C.sel); fillRound(gc, 2, y - 1, W - 4, rh - 1, 4)
      col(gc, C.selBar); fillRound(gc, 2, y - 1, 3, rh - 1, 1)
    elseif i <= #items and i % 2 == 0 then
      col(gc, C.row); gc:fillRect(2, y - 1, W - 4, rh - 1)
    end
    if i <= #items then
      col(gc, C.tx); font(gc, 9, i == self.sel and "b" or "r")
      gc:drawString(fitText(gc, items[i], W - 12), 8, y, "top")
    else
      col(gc, C.accent); gc:fillArc(7, y + 1, 11, 11, 0, 360)
      col(gc, C.titleTx); gc:fillRect(10, y + 6, 5, 1); gc:fillRect(12, y + 4, 1, 5)
      col(gc, C.title2); font(gc, 9, "b")
      gc:drawString("Agregar nuevo", 22, y, "top")
    end
  end
  drawFooter(gc, self.msg or "enter:editar|+:nuevo|del:borrar|p:vista|esc:volver", self.msg ~= nil)
  self.msg = nil
end
function List:arrow(k)
  local n = #self.items() + 1
  if k == "up" then self.sel = self.sel > 1 and self.sel - 1 or n
  elseif k == "down" then self.sel = self.sel < n and self.sel + 1 or 1 end
end
function List:enter()
  if self.sel > #self.items() then self.add() else self.edit(self.sel) end
end
function List:char(ch)
  if ch == "+" then self.add() elseif ch == "p" or ch == "P" then cyclePreview() end
end
function List:back()
  if self.sel <= #self.items() then
    local e = self.del(self.sel)
    if e then self.msg = e end
    modelChanged()
  end
end
function List:esc() pop() end
function List:click(x, y)
  if x > SPLIT then return end
  local i = self.off + floor((y - 21) / 15) + 1
  if i >= 1 and i <= #self.items() + 1 then
    if i == self.sel then self:enter() else self.sel = i end
  end
end

------------------------- Formulario ------------------------------------------
local function utf8chars(s)
  local t = {}
  for ch in tostring(s or ""):gmatch("[%z\1-\127\194-\244][\128-\191]*") do t[#t + 1] = ch end
  return t
end

local Form = {}; Form.__index = Form
-- fields: {k=, label= (str|fn), t="txt"|"opt"|"chk", opts=, help=, show=fn(vals), change=fn(vals)}
function Form.new(title, fields, vals, ok)
  local f = setmetatable({title = title, fields = fields, vals = vals, ok = ok,
                          sel = 1, off = 0, cur = nil, fresh = true, split = App.model ~= nil}, Form)
  return f
end
function Form:visible()
  local v = {}
  for _, fd in ipairs(self.fields) do
    if not fd.show or fd.show(self.vals) then v[#v + 1] = fd end
  end
  return v
end
function Form:curField()
  local v = self:visible()
  return v[self.sel], v
end
function Form:resetCursor()
  local fd = self:curField()
  if fd and fd.t == "txt" then self.cur = #utf8chars(self.vals[fd.k]) end
  self.fresh = true   -- el primer caracter tecleado reemplaza el contenido
end
function Form:paint(gc)
  local vis = self:visible()
  local n = #vis + 1
  if self.sel > n then self.sel = n end
  if not self.cur then self:resetCursor() end
  drawHeader(gc, self.title)
  local realW = W
  local split = self.split
  if split then drawPreviewPane(gc, self, H - 30) end
  local W = split and SPLIT or W
  col(gc, C.bg); gc:fillRect(0, 18, W, H - 47)
  local rh, y0 = split and 26 or 17, 21
  self.rh = rh
  local nv = floor((H - 30 - y0) / rh)
  if self.sel > self.off + nv then self.off = self.sel - nv end
  if self.sel <= self.off then self.off = self.sel - 1 end
  local lx, vx, vy, bh = 6, 124, 0, 14
  if split then vx, vy = 5, 11 end
  for r = 1, nv do
    local i = self.off + r
    if i > n then break end
    local y = y0 + (r - 1) * rh
    local selected = (i == self.sel)
    if selected and i <= #vis then
      col(gc, C.sel); fillRound(gc, 2, y - 2, W - 4, rh, 4)
      col(gc, C.selBar); fillRound(gc, 2, y - 2, 3, rh, 1)
    end
    if i <= #vis then
      local fd = vis[i]
      local lab = type(fd.label) == "function" and fd.label(self.vals) or fd.label
      font(gc, split and 7 or 9, selected and "b" or "r"); col(gc, selected and C.selBar or C.dim)
      if not split then col(gc, C.tx) end
      gc:drawString(fitText(gc, lab, split and (W - 10) or (vx - lx - 4)), lx + (split and 2 or 0),
                    y + (split and -2 or 0), "top")
      local v = self.vals[fd.k]
      y = y + vy
      local bw = W - vx - 6
      if fd.t == "txt" then
        col(gc, C.card); fillRound(gc, vx, y, bw, bh, 3)
        if selected then
          col(gc, C.selBar); strokeRound(gc, vx, y, bw, bh, 3)
          strokeRound(gc, vx + 1, y + 1, bw - 2, bh - 2, 2)
        else
          col(gc, C.border); strokeRound(gc, vx, y, bw, bh, 3)
        end
        col(gc, C.tx); font(gc, 10)
        local chars = utf8chars(v)
        local wmax = bw - 8
        local first = 1
        local cur = selected and self.cur or #chars
        local function seg(a, b) return table.concat(chars, "", a, b) end
        while first < cur and gc:getStringWidth(seg(first, cur)) > wmax do first = first + 1 end
        local shown = seg(first, #chars)
        gc:drawString(fitText(gc, shown, wmax + 4), vx + 4, y - 1, "top")
        if selected then
          local cx = vx + 4 + gc:getStringWidth(seg(first, cur))
          col(gc, C.accent); gc:fillRect(cx, y + 2, 2, bh - 4)
        end
      elseif fd.t == "opt" then
        col(gc, selected and C.selBar or C.title2); fillRound(gc, vx, y, bw, bh, 7)
        col(gc, C.titleTx); font(gc, 9, "b")
        local txt = fitText(gc, fd.opts[v] or "?", bw - 26)
        gc:drawString(txt, vx + bw / 2 - gc:getStringWidth(txt) / 2, y, "top")
        gc:fillPolygon({vx + 5, y + 7, vx + 10, y + 3, vx + 10, y + 11, vx + 5, y + 7})
        gc:fillPolygon({vx + bw - 5, y + 7, vx + bw - 10, y + 3, vx + bw - 10, y + 11, vx + bw - 5, y + 7})
      elseif fd.t == "chk" then
        -- interruptor
        col(gc, v and C.ok or C.light); fillRound(gc, vx, y + 1, 24, 12, 6)
        col(gc, C.card); gc:fillArc(v and (vx + 13) or (vx + 1), y + 2, 10, 10, 0, 360)
        col(gc, C.tx); font(gc, 9, "b")
        gc:drawString(v and "Si" or "No", vx + 30, y, "top")
      end
    else
      local bw = split and (W - 20) or 110
      local bx = (W - bw) / 2
      local by = y + (split and 4 or 1)
      col(gc, selected and C.accent or C.title2); fillRound(gc, bx, by, bw, 16, 8)
      col(gc, C.titleTx); font(gc, 10, "b")
      gc:drawString("ACEPTAR", bx + bw / 2 - gc:getStringWidth("ACEPTAR") / 2, by, "top")
    end
  end
  -- ayuda
  local fd = vis[self.sel]
  local help = self.err or (fd and fd.help) or "enter en ACEPTAR para guardar"
  if fd and not self.err then
    if fd.t == "opt" then help = (fd.help and fd.help .. "  " or "") .. "<- -> cambia"
    elseif fd.t == "chk" then help = "<- -> o espacio cambia" end
  end
  col(gc, self.err and {253, 232, 232} or {252, 246, 226}); gc:fillRect(0, H - 29, realW, 15)
  col(gc, self.err and C.err or C.accent); gc:fillRect(0, H - 29, 3, 15)
  col(gc, self.err and C.err or C.tx); font(gc, 7, self.err and "b" or "r")
  gc:drawString(fitText(gc, help, realW - 10), 7, H - 27, "top")
  drawFooter(gc, "enter:siguiente|<- ->:editar|esc:cancelar")
end
function Form:move(d)
  local n = #self:visible() + 1
  self.sel = self.sel + d
  if self.sel < 1 then self.sel = n end
  if self.sel > n then self.sel = 1 end
  self:resetCursor()
end
function Form:arrow(k)
  local fd = self:curField()
  if k == "up" then self:move(-1); return end
  if k == "down" then self:move(1); return end
  if not fd then return end
  local d = (k == "right") and 1 or -1
  if fd.t == "txt" then
    local n = #utf8chars(self.vals[fd.k])
    self.cur = max(0, min(n, (self.cur or n) + d))
    self.fresh = false
  elseif fd.t == "opt" then
    local v = (self.vals[fd.k] or 1) + d
    if v < 1 then v = #fd.opts end
    if v > #fd.opts then v = 1 end
    self.vals[fd.k] = v
    if fd.change then fd.change(self.vals) end
  elseif fd.t == "chk" then
    self.vals[fd.k] = not self.vals[fd.k]
    if fd.change then fd.change(self.vals) end
  end
end
function Form:char(ch)
  self.err = nil
  local fd = self:curField()
  if not fd then return end
  if fd.t == "txt" then
    if self.fresh then self.vals[fd.k] = ""; self.cur = 0; self.fresh = false end
    local chars = utf8chars(self.vals[fd.k])
    local c = self.cur or #chars
    table.insert(chars, c + 1, ch)
    self.vals[fd.k] = table.concat(chars)
    self.cur = c + #utf8chars(ch)
  elseif fd.t == "chk" then
    self.vals[fd.k] = not self.vals[fd.k]
    if fd.change then fd.change(self.vals) end
  end
end
function Form:back()
  local fd = self:curField()
  if fd and fd.t == "txt" then
    self.fresh = false
    local chars = utf8chars(self.vals[fd.k])
    local c = self.cur or #chars
    if c > 0 then
      table.remove(chars, c)
      self.vals[fd.k] = table.concat(chars)
      self.cur = c - 1
    end
  end
end
function Form:clear()
  local fd = self:curField()
  if fd and fd.t == "txt" then self.vals[fd.k] = ""; self.cur = 0 end
end
function Form:enter()
  local fd, vis = self:curField()
  if self.sel > #vis then
    local e = self.ok(self.vals)
    if e then self.err = e else modelChanged(); pop() end
  else
    self:move(1)
  end
end
function Form:tab() self:move(1) end
function Form:esc() pop() end
function Form:click(x, y)
  if self.split and x > SPLIT then return end
  local i = self.off + floor((y - 20) / (self.rh or 16)) + 1
  if i >= 1 and i <= #self:visible() + 1 then
    if i == self.sel and i > #self:visible() then self:enter() else self.sel = i; self:resetCursor() end
  end
end

------------------------- Visor de texto --------------------------------------
local TextView = {}; TextView.__index = TextView
function TextView.new(title, lines, footer)
  return setmetatable({title = title, raw = lines, off = 0, footer = footer}, TextView)
end
function TextView:paint(gc)
  font(gc, 9)
  if not self.lines then
    self.lines = {}
    for _, ln in ipairs(self.raw) do
      local s, style = ln, nil
      if type(ln) == "table" then s, style = ln[1], ln[2] end
      local ind = ""
      if gc:getStringWidth(s) <= W - 24 then
        self.lines[#self.lines + 1] = {s, style}
      else
        local cur = ""
        for word in s:gmatch("%S+%s*") do
          if gc:getStringWidth(cur .. word) > W - 24 and cur ~= "" then
            self.lines[#self.lines + 1] = {cur, style}
            cur = "   " .. word
          else
            cur = cur .. word
          end
        end
        if cur ~= "" then self.lines[#self.lines + 1] = {cur, style} end
      end
    end
  end
  drawHeader(gc, self.title)
  col(gc, C.bg); gc:fillRect(0, 18, W, H - 32)
  col(gc, C.card); fillRound(gc, 3, 20, W - 6, H - 37, 5)
  col(gc, C.border); strokeRound(gc, 3, 20, W - 7, H - 38, 5)
  local rh, y0 = 13, 22
  self.vis = floor((H - 19 - y0) / rh)
  local maxoff = max(0, #self.lines - self.vis)
  if self.off > maxoff then self.off = maxoff end
  for r = 1, self.vis do
    local ln = self.lines[self.off + r]
    if not ln then break end
    local st = ln[2]
    if st == "h" then
      col(gc, C.accent); gc:fillRect(6, y0 + (r - 1) * rh + 2, 2, 9)
      col(gc, C.title2); font(gc, 9, "b")
    elseif st == "e" then col(gc, C.err); font(gc, 9, "b")
    elseif st == "d" then col(gc, C.dim); font(gc, 9)
    else col(gc, C.tx); font(gc, 9) end
    gc:drawString(ln[1], st == "h" and 11 or 8, y0 + (r - 1) * rh, "top")
  end
  if #self.lines > self.vis then
    local bh = H - 19 - y0
    local th = max(10, bh * self.vis / #self.lines)
    local ty = y0 + (bh - th) * (maxoff > 0 and self.off / maxoff or 0)
    col(gc, C.light); fillRound(gc, W - 9, ty, 3, th, 1)
  end
  drawFooter(gc, self.footer or "^v:desplazar|<- ->:pagina|esc:volver")
end
function TextView:arrow(k)
  local vis = self.vis or 10
  if k == "up" then self.off = max(0, self.off - 1)
  elseif k == "down" then self.off = self.off + 1
  elseif k == "left" then self.off = max(0, self.off - vis)
  elseif k == "right" then self.off = self.off + vis end
end
function TextView:esc() pop() end
function TextView:enter() pop() end

msgBox = function(title, lines)
  if type(lines) == "string" then lines = {lines} end
  push(TextView.new(title, lines, "enter:cerrar|esc:cerrar|^v:desplazar"))
end

------------------------- Dibujo de la estructura ----------------------------
local function makeView(nodes, x0, y0, w, h, pad)
  local xmin, xmax, ymin, ymax = 1e300, -1e300, 1e300, -1e300
  for _, n in ipairs(nodes) do
    xmin = min(xmin, n.x); xmax = max(xmax, n.x)
    ymin = min(ymin, n.y); ymax = max(ymax, n.y)
  end
  if xmin > xmax then xmin, xmax, ymin, ymax = 0, 1, 0, 1 end
  local dx, dy = xmax - xmin, ymax - ymin
  local sx = dx > 1e-12 and (w - 2 * pad) / dx or 1e300
  local sy = dy > 1e-12 and (h - 2 * pad) / dy or 1e300
  local sc = min(sx, sy)
  if sc > 1e299 then sc = 1 end
  local cx = x0 + w / 2 - (xmin + xmax) / 2 * sc
  local cy = y0 + h / 2 + (ymin + ymax) / 2 * sc
  return function(x, y) return cx + x * sc, cy - y * sc end, sc
end

local function drawArrow(gc, x1, y1, x2, y2)
  gc:drawLine(x1, y1, x2, y2)
  local dx, dy = x2 - x1, y2 - y1
  local L = sqrt(dx * dx + dy * dy)
  if L < 1 then return end
  dx, dy = dx / L, dy / L
  local hx, hy = 5, 3
  gc:drawLine(x2, y2, x2 - hx * dx + hy * dy, y2 - hx * dy - hy * dx)
  gc:drawLine(x2, y2, x2 - hx * dx - hy * dy, y2 - hx * dy + hy * dx)
end

local function drawSupport(gc, px, py, S)
  local a = (S.ang or 0) * PI / 180
  local ca, sa = math.cos(a), math.sin(a)
  local function P(u, v, rot90)
    if rot90 then u, v = -v, u end   -- glifo girado 90 grados (restriccion en X)
    return px + u * ca + v * sa, py - u * sa + v * ca
  end
  local function line(u1, v1, u2, v2, r)
    local x1, y1 = P(u1, v1, r); local x2, y2 = P(u2, v2, r)
    gc:drawLine(x1, y1, x2, y2)
  end
  col(gc, C.sup)
  local rx, ry, rz = S.rx, S.ry, S.rz
  if rx and ry and rz then
    line(-9, 0, 9, 0)
    for i = 0, 4 do line(-8 + i * 4, 0, -11 + i * 4, 5) end
  elseif rx and ry then
    line(0, 0, -6, 9); line(0, 0, 6, 9); line(-9, 9, 9, 9)
    for i = 0, 4 do line(-7 + i * 4, 9, -10 + i * 4, 13) end
  elseif rz and (rx or ry) then
    local r = rx and not ry
    line(-7, 0, 7, 0, r); line(-7, 0, -7, -4, r); line(7, 0, 7, -4, r)
    line(-9, 4, 9, 4, r)
    for i = 0, 4 do line(-7 + i * 4, 4, -10 + i * 4, 8, r) end
  elseif ry or rx then
    local r = rx and not ry
    line(0, 0, -6, 8, r); line(0, 0, 6, 8, r); line(-6, 8, 6, 8, r)
    line(-9, 12, 9, 12, r)
  elseif rz then
    gc:drawRect(px - 4, py - 4, 8, 8)
  end
  if (S.kx or 0) ~= 0 or (S.ky or 0) ~= 0 or (S.kr or 0) ~= 0 then
    font(gc, 7); gc:drawString("k", px + 5, py + 2, "top")
  end
end

local function supFlags(model, G)
  local sup = {}
  for _, s in ipairs(model.sups) do
    local S = sup[s.node] or {}
    S.rx = S.rx or s.rx; S.ry = S.ry or s.ry
    S.rz = (S.rz or s.rz) and model.kind ~= "truss"
    local ok, a = pcall(ev, s.ang, G.gres, "", 0)
    S.ang = ok and a or 0
    local okk, kx = pcall(ev, s.kx, G.gres, "", 0)
    S.kx = okk and kx or 0
    sup[s.node] = S
  end
  return sup
end

local function drawFrame(gc, G, tr, sup, opts)
  opts = opts or {}
  gc:setPen("medium", "smooth")
  for m, mb in ipairs(G.mems) do
    local x1, y1 = tr(G.nodes[mb.i].x, G.nodes[mb.i].y)
    local x2, y2 = tr(G.nodes[mb.j].x, G.nodes[mb.j].y)
    col(gc, opts.sel == m and C.load or (opts.memColor or C.mem))
    gc:drawLine(x1, y1, x2, y2)
  end
  gc:setPen("thin", "smooth")
  -- articulaciones
  if G.kind ~= "truss" then
    for _, mb in ipairs(G.mems) do
      local ends = {}
      if mb.ri then ends[#ends + 1] = {mb.i, 1} end
      if mb.rj then ends[#ends + 1] = {mb.j, -1} end
      for _, e in ipairs(ends) do
        local n = G.nodes[e[1]]
        local x, y = tr(n.x, n.y)
        local cx, cy = x + e[2] * 6 * mb.c, y - e[2] * 6 * mb.s
        col(gc, {255, 255, 255}); gc:fillArc(cx - 3, cy - 3, 6, 6, 0, 360)
        col(gc, C.mem); gc:drawArc(cx - 3, cy - 3, 6, 6, 0, 360)
      end
    end
  end
  for n, nd in ipairs(G.nodes) do
    local x, y = tr(nd.x, nd.y)
    if sup[n] then drawSupport(gc, x, y, sup[n]) end
    col(gc, C.tx)
    if G.kind == "truss" then
      col(gc, {255, 255, 255}); gc:fillArc(x - 2, y - 2, 5, 5, 0, 360)
      col(gc, C.tx); gc:drawArc(x - 2, y - 2, 5, 5, 0, 360)
    else
      gc:fillRect(x - 1, y - 1, 3, 3)
    end
  end
  if opts.labels then
    font(gc, 7)
    for n, nd in ipairs(G.nodes) do
      local x, y = tr(nd.x, nd.y)
      col(gc, C.N); gc:drawString(tostring(n), x + 3, y - 12, "top")
    end
    for m, mb in ipairs(G.mems) do
      local x1, y1 = tr(G.nodes[mb.i].x, G.nodes[mb.i].y)
      local x2, y2 = tr(G.nodes[mb.j].x, G.nodes[mb.j].y)
      col(gc, C.M); gc:drawString("b" .. m, (x1 + x2) / 2 - 4 - 8 * mb.s, (y1 + y2) / 2 - 6 - 8 * mb.c, "top")
    end
  end
end

------------------------- Ver estructura -------------------------------------
local function drawLoads(gc, model, G, tr)
  local vt = G.vt
  local lres = function(n)
    local r = vt[n]
    if r and r.value then return lin_const(r.value) end
    return lin_const(1)
  end
  local function num(s) local ok, v = pcall(evalConst, s or "", lres); return ok and v or 0 end
  col(gc, C.load); font(gc, 7)
  for _, l in ipairs(model.nl) do
    local nd = G.nodes[l.node]
    if nd then
      local x, y = tr(nd.x, nd.y)
      local fx, fy = num(l.fx), num(l.fy)
      if fx ~= 0 then
        local sg = fx > 0 and 1 or -1
        drawArrow(gc, x - sg * 22, y, x - sg * 3, y)
        gc:drawString(l.fx, x - sg * 22 - (sg > 0 and 12 or -2), y - 10, "top")
      end
      if fy ~= 0 then
        local sg = fy > 0 and 1 or -1
        drawArrow(gc, x, y + sg * 22, x, y + sg * 3)
        gc:drawString(l.fy, x + 3, y + sg * 18 - 5, "top")
      end
      if model.kind ~= "truss" and num(l.m) ~= 0 then
        gc:drawArc(x - 9, y - 9, 18, 18, 30, 240)
        gc:drawString(l.m, x + 8, y + 3, "top")
      end
    end
  end
  for _, l in ipairs(model.ml) do
    local mb = G.mems[l.mem]
    if mb then
      local ni = G.nodes[mb.i]
      local pres = function(n) if n == "L" then return lin_const(mb.L) end return lres(n) end
      local function pos(s, d) if isBlank(s) then return d end local ok, v = pcall(evalConst, s, pres); return ok and v or d end
      if l.t == 1 or l.t == 3 then
        local ux, uy = dirLocal(l.dir or 1, mb.c, mb.s, false)
        local gx, gy = mb.c * ux - mb.s * uy, mb.s * ux + mb.c * uy
        local v1 = num(l.v1)
        local v2 = (l.t == 3 and not isBlank(l.v2)) and num(l.v2) or v1
        local a = pos(l.a, l.t == 1 and mb.L / 2 or 0)
        local b = l.t == 3 and pos(l.b, mb.L) or a
        local nar = l.t == 1 and 1 or 6
        local vm = max(abs(v1), abs(v2), 1e-12)
        for q = 0, nar - 1 do
          local t = nar == 1 and 0 or q / (nar - 1)
          local sx = a + (b - a) * t
          local val = v1 + (v2 - v1) * t
          local len = (l.t == 1 and 22 or 16) * abs(val) / vm
          if len > 2 then
            local sg = val > 0 and 1 or -1
            local x, y = tr(ni.x + sx * mb.c, ni.y + sx * mb.s)
            drawArrow(gc, x - sg * gx * len, y + sg * gy * len, x, y)
          end
        end
        local sx = (a + b) / 2
        local x, y = tr(ni.x + sx * mb.c, ni.y + sx * mb.s)
        local txt = l.v1 .. ((l.t == 3 and not isBlank(l.v2)) and (".." .. l.v2) or "")
        local sg = v1 >= 0 and 1 or -1
        gc:drawString(txt, x - sg * gx * 26 + 2, y + sg * gy * 26 - 6, "top")
      elseif l.t == 2 then
        local a = pos(l.a, mb.L / 2)
        local x, y = tr(ni.x + a * mb.c, ni.y + a * mb.s)
        gc:drawArc(x - 8, y - 8, 16, 16, 30, 240)
        gc:drawString(l.v1, x + 7, y + 2, "top")
      else
        local x, y = tr(ni.x + mb.L / 2 * mb.c, ni.y + mb.L / 2 * mb.s)
        gc:drawString(l.t == 4 and "T" or "e", x + 3, y + 1, "top")
      end
    end
  end
end

local StructView = {}; StructView.__index = StructView
function StructView.new()
  local s = setmetatable({labels = true, loads = true}, StructView)
  local ok, G = pcall(Eng.geometry, App.model)
  if not ok then msgBox("Error", tostring(G)); return nil end
  G.kind = App.model.kind
  s.G = G
  s.sup = supFlags(App.model, G)
  return s
end
function StructView:paint(gc)
  local G = self.G
  drawHeader(gc, "Estructura", #G.nodes .. " nudos, " .. #G.mems .. " barras")
  gridBg(gc, 0, 18, W, H - 32, 14)
  if #G.nodes == 0 then drawFooter(gc, "Sin nudos"); return end
  local tr, sc = makeView(G.nodes, 0, 18, W, H - 32, 32)
  drawFrame(gc, G, tr, self.sup, {labels = self.labels})
  if self.loads then drawLoads(gc, App.model, G, tr) end
  drawFooter(gc, "tab:numeros|c:cargas|esc:volver")
end
function StructView:tab() self.labels = not self.labels end
function StructView:char(ch) if ch == "c" or ch == "C" then self.loads = not self.loads end end
function StructView:esc() pop() end
function StructView:enter() pop() end

------------------------- Diagramas ------------------------------------------
local QN = {"N", "V", "M", "D"}
local QNAME = {N = "Axial N", V = "Corte V", M = "Momento M", D = "Deformada"}
local Diagram = {}; Diagram.__index = Diagram
function Diagram.new(q, res)
  res = res or App.res
  local d = setmetatable({q = q or (res.kind == "truss" and "N" or "M"), sel = 1, res = res,
                          zoom = 1, labels = res.kind == "truss"}, Diagram)
  d.G = {nodes = res.nodes, mems = res.mems, kind = res.kind}
  d.sup = {}
  for n, S in pairs(res.sup) do d.sup[n] = S end
  d.data = {}
  for m in ipairs(res.mems) do
    local pts = Eng.samples(res, m, 24)
    local vals = {}
    for i, p in ipairs(pts) do
      local N, V, M = Eng.internal(res, m, p[1], p[2])
      vals[i] = {p[1], N, V, M}
    end
    d.data[m] = {pts = vals, def = Eng.deflect(res, m, 24), ext = Eng.extremes(res, m)}
  end
  return d
end
function Diagram:paint(gc)
  local res = self.res
  local q = self.q
  local qi = ({N = 2, V = 3, M = 4})[q]
  local mb = res.mems[self.sel]
  drawHeader(gc, "Diagrama " .. QNAME[q], "barra " .. self.sel .. "/" .. #res.mems)
  gridBg(gc, 0, 18, W, H - 47, 14)
  self:drawBody(gc, 0, 18, W, H - 46, 30)
  -- informacion
  col(gc, {236, 242, 252}); gc:fillRect(0, H - 29, W, 15)
  col(gc, C[q]); gc:fillRect(0, H - 29, 3, 15)
  font(gc, 7, "b"); col(gc, C.tx)
  local info
  if q == "D" then
    local D = {0, 0, 0, 0, 0, 0}
    for kk = 0, res.ns do
      for a = 1, 6 do D[a] = D[a] + res.w[kk] * res.cases[kk].D[self.sel][a] end
    end
    info = string.format("b%d (%d-%d)  giro i=%s  giro j=%s", self.sel, mb.i, mb.j, fmtShort(D[3]), fmtShort(D[6]))
  else
    local e = self.data[self.sel].ext[qi - 1]
    info = string.format("b%d (%d-%d) L=%s  max=%s @%s  min=%s @%s", self.sel, mb.i, mb.j, fmt(mb.L, 3),
                         fmtShort(e.mx), fmt(e.xmx, 2), fmtShort(e.mn), fmt(e.xmn, 2))
  end
  gc:drawString(fitText(gc, info, W - 10), 7, H - 27, "top")
  drawFooter(gc, "<- ->:barra|m v n d:tipo|^v:escala|tab:valores|enter:tabla")
end

-- dibuja el diagrama dentro del rectangulo (x0, y0, w, h)
function Diagram:drawBody(gc, x0, y0, w, h, pad)
  local res = self.res
  local q = self.q
  local qi = ({N = 2, V = 3, M = 4})[q]
  local tr, sc = makeView(res.nodes, x0, y0, w, h, pad)
  gc:setPen("thin", "smooth")
  drawFrame(gc, self.G, tr, self.sup, {memColor = C.light, labels = false})
  -- escala
  local mx = 0
  for m in ipairs(res.mems) do
    if q == "D" then
      for _, p in ipairs(self.data[m].def) do mx = max(mx, sqrt(p[2] ^ 2 + p[3] ^ 2)) end
    else
      for _, p in ipairs(self.data[m].pts) do mx = max(mx, abs(p[qi])) end
    end
  end
  local colr = C[q]
  local maxpix = 0.16 * min(w, h) * self.zoom
  local k = mx > 1e-12 and maxpix / mx or 0
  if q == "D" then k = mx > 1e-12 and (0.12 * min(w, h) * self.zoom) / mx or 0 end
  local trussStyle = res.kind == "truss" and q == "N" and not self.offsetMode
  for m, mbm in ipairs(res.mems) do
    local ni = res.nodes[mbm.i]
    local thick = (m == self.sel)
    gc:setPen(thick and "medium" or "thin", "smooth")
    if trussStyle then
      -- enrejado: barra coloreada (azul traccion, rojo compresion) + valor
      local pts = self.data[m].pts
      local v = pts[math.ceil(#pts / 2)][2]
      local nj = res.nodes[mbm.j]
      local x1, y1 = tr(ni.x, ni.y)
      local x2, y2 = tr(nj.x, nj.y)
      local tol = 1e-9 * max(mx, 1e-30)
      if abs(v) <= tol then col(gc, C.light) elseif v > 0 then col(gc, C.N) else col(gc, C.M) end
      gc:setPen(thick and "thick" or "medium", "smooth")
      gc:drawLine(x1, y1, x2, y2)
      if thick then col(gc, C.load); gc:drawLine(x1, y1, x2, y2) end
      if self.labels or thick then
        font(gc, 7)
        local s = fmtShort(v)
        local tw = gc:getStringWidth(s)
        local f = 0.5 + ((m % 3) - 1) * 0.12
        local cx, cy = x1 + (x2 - x1) * f, y1 + (y2 - y1) * f
        col(gc, {255, 255, 255}); gc:fillRect(cx - tw / 2 - 1, cy - 5, tw + 2, 10)
        col(gc, C.tx); gc:drawString(s, cx - tw / 2, cy - 6, "top")
      end
    elseif q == "D" then
      col(gc, colr)
      local poly = {}
      for _, p in ipairs(self.data[m].def) do
        local x, y = tr(ni.x + p[1] * mbm.c + k * p[2] / sc * 1, ni.y + p[1] * mbm.s + k * p[3] / sc)
        poly[#poly + 1] = x; poly[#poly + 1] = y
      end
      if #poly >= 4 then gc:drawPolyLine(poly) end
    else
      local sgn = (q == "M") and -1 or 1
      local nx, ny = -mbm.s, mbm.c
      local poly = {}
      local pts = self.data[m].pts
      for i, p in ipairs(pts) do
        local bx, by = tr(ni.x + p[1] * mbm.c, ni.y + p[1] * mbm.s)
        local off = sgn * p[qi] * k
        local x, y = bx + off * nx, by - off * ny
        poly[#poly + 1] = x; poly[#poly + 1] = y
        if i % 2 == 1 or i == #pts then
          col(gc, thick and {250, 190, 150} or {215, 215, 235})
          gc:drawLine(bx, by, x, y)
        end
      end
      col(gc, colr)
      local x1, y1 = tr(ni.x, ni.y)
      local nj = res.nodes[mbm.j]
      local x2, y2 = tr(nj.x, nj.y)
      if #poly >= 2 then
        gc:drawLine(x1, y1, poly[1], poly[2])
        gc:drawLine(x2, y2, poly[#poly - 1], poly[#poly])
      end
      if #poly >= 4 then gc:drawPolyLine(poly) end
      -- etiquetas
      if self.labels or thick then
        font(gc, 7); col(gc, C.tx)
        local e = self.data[m].ext[qi - 1]
        local marks = {}
        local function mark(x, v)
          if abs(v) < 1e-9 * max(mx, 1e-30) then return end
          for _, o in ipairs(marks) do if abs(o - x) < 0.08 * mbm.L then return end end
          marks[#marks + 1] = x
          local bx, by = tr(ni.x + x * mbm.c, ni.y + x * mbm.s)
          local off = sgn * v * k
          local s = fmtShort(v)
          local tw = gc:getStringWidth(s)
          local px, py = bx + off * nx, by - off * ny
          local dx = (sgn * v >= 0) and 1 or -1
          gc:drawString(s, px + dx * 3 * nx - (nx < -0.3 and tw or 0) - (abs(nx) <= 0.3 and tw / 2 or 0),
                        py - dx * 3 * ny - ((dx * ny) > 0.3 and 10 or 0), "top")
        end
        local p0, pL = pts[1], pts[#pts]
        if res.kind == "truss" and #res.cases[0].L[m] == 0 then
          mark(mbm.L / 2, pts[math.ceil(#pts / 2)][qi])
        else
          mark(0, p0[qi]); mark(mbm.L, pL[qi])
          mark(e.xmx, e.mx); mark(e.xmn, e.mn)
        end
      end
    end
  end
  gc:setPen("thin", "smooth")
  if trussStyle then
    for n, nd in ipairs(res.nodes) do
      local x, y = tr(nd.x, nd.y)
      col(gc, {255, 255, 255}); gc:fillArc(x - 2, y - 2, 5, 5, 0, 360)
      col(gc, C.tx); gc:drawArc(x - 2, y - 2, 5, 5, 0, 360)
    end
    font(gc, 7); col(gc, C.N); gc:drawString("traccion", x0 + 4, y0 + 2, "top")
    col(gc, C.M); gc:drawString("compresion", x0 + 48, y0 + 2, "top")
  end
  return mx
end
function Diagram:arrow(k)
  local n = #self.res.mems
  if k == "left" then self.sel = self.sel > 1 and self.sel - 1 or n
  elseif k == "right" then self.sel = self.sel < n and self.sel + 1 or 1
  elseif k == "up" then self.zoom = self.zoom * 1.3
  elseif k == "down" then self.zoom = self.zoom / 1.3 end
end
function Diagram:char(ch)
  local c = ch:upper()
  if c == "M" or c == "V" or c == "N" or c == "D" then self.q = c; return end
  if c == "O" then self.offsetMode = not self.offsetMode; return end
  if ch == "+" then self.zoom = self.zoom * 1.3 elseif ch == "-" then self.zoom = self.zoom / 1.3 end
  local n = tonumber(ch)
  if n and n >= 1 and n <= 4 then self.q = QN[n] end
end
function Diagram:tab() self.labels = not self.labels end
function Diagram:enter() App.memberTable(self.sel) end
function Diagram:esc() pop() end

------------------------- Vista previa en vivo -------------------------------
local PV = {mode = 1}
local PV_MODES = {"E", "M", "V", "N", "D"}
local PV_NAMES = {E = "Estructura y cargas", M = "Momento M", V = "Corte V", N = "Axial N", D = "Deformada"}

local function serialize(v, out)
  if type(v) == "table" then
    local keys = {}
    for k in pairs(v) do keys[#keys + 1] = k end
    table.sort(keys, function(a, b) return tostring(a) < tostring(b) end)
    out[#out + 1] = "{"
    for _, k in ipairs(keys) do
      out[#out + 1] = tostring(k) .. "="
      serialize(v[k], out)
      out[#out + 1] = ";"
    end
    out[#out + 1] = "}"
  else
    out[#out + 1] = tostring(v)
  end
end

local function deepcopy(v)
  if type(v) ~= "table" then return v end
  local r = {}
  for k, x in pairs(v) do r[k] = deepcopy(x) end
  return r
end

cyclePreview = function() PV.mode = PV.mode % #PV_MODES + 1 end

-- modelo a mostrar: incluye lo que se esta escribiendo en el formulario
local function previewModel(scr)
  local model = App.model
  if scr.ok and scr.vals then
    local copy = deepcopy(model)
    App.model = copy
    local ok, e = pcall(scr.ok, deepcopy(scr.vals))
    App.model = model
    if ok and e == nil then return copy end
  end
  return model
end

-- grado de hiperestaticidad
local function ghText(model, G)
  local n, m = #G.nodes, #G.mems
  local truss = model.kind == "truss"
  local sup = {}
  for _, s in ipairs(model.sups) do
    local S = sup[s.node] or {}
    S.rx = S.rx or s.rx or not isBlank(s.kx)
    S.ry = S.ry or s.ry or not isBlank(s.ky)
    S.rz = (S.rz or s.rz or not isBlank(s.kr)) and not truss
    sup[s.node] = S
  end
  local r = 0
  for _, S in pairs(sup) do r = r + (S.rx and 1 or 0) + (S.ry and 1 or 0) + (S.rz and 1 or 0) end
  local gh
  if truss then
    gh = m + r - 2 * n
  else
    local rel, ends, relEnds = 0, {}, {}
    for _, mb in ipairs(G.mems) do
      for _, e in ipairs({{mb.i, mb.ri}, {mb.j, mb.rj}}) do
        ends[e[1]] = (ends[e[1]] or 0) + 1
        if e[2] then rel = rel + 1; relEnds[e[1]] = (relEnds[e[1]] or 0) + 1 end
      end
    end
    local auto = 0
    for nd, c in pairs(ends) do
      if relEnds[nd] == c and not (sup[nd] and sup[nd].rz) then auto = auto + 1 end
    end
    gh = 3 * m + r - 3 * n - rel + auto
  end
  if gh == 0 then return "Estable: ISOSTATICA" end
  return "Estable: HIPERESTATICA grado " .. gh
end

local function previewState(model)
  local o = {}
  serialize(model, o)
  local sig = table.concat(o)
  if PV.sig == sig then return PV.st end
  local st = {}
  local ok, G = pcall(Eng.geometry, model)
  if not ok then
    st.err = tostring(G)
  else
    G.kind = model.kind
    st.G = G
    st.sup = supFlags(model, G)
    if #G.nodes == 0 then st.msg = "Agregue nudos"
    elseif #G.mems == 0 then st.msg = "Agregue barras"
    elseif #model.sups == 0 then st.msg = "Agregue apoyos"
    else
      local ok2, res = pcall(Eng.solve, model)
      if ok2 then
        st.res = res
        st.msg = ghText(model, G)
        st.good = true
      else
        st.err = tostring(res)
      end
    end
  end
  PV.sig, PV.st = sig, st
  return st
end

local function wrapLines(gc, text, wmax, maxl)
  local lines, cur = {}, ""
  for word in text:gmatch("%S+") do
    local t = cur == "" and word or (cur .. " " .. word)
    if gc:getStringWidth(t) > wmax and cur ~= "" then
      lines[#lines + 1] = cur; cur = word
    else
      cur = t
    end
  end
  if cur ~= "" then lines[#lines + 1] = cur end
  while #lines > maxl do table.remove(lines) end
  return lines
end

drawPreviewPane = function(gc, scr, ybot)
  local x0 = SPLIT + 1
  local w = W - x0
  local y0 = 18
  gridBg(gc, x0, y0, w, ybot - y0, 12)
  col(gc, C.border); gc:fillRect(SPLIT, y0, 1, ybot - y0)
  col(gc, C.selBar); gc:fillRect(SPLIT + 1, y0, 1, ybot - y0)
  if not App.model then return end
  local model = previewModel(scr)
  local st = previewState(model)
  local mode = PV_MODES[PV.mode]
  font(gc, 7, "b")
  local tag = PV_NAMES[mode]
  local tw = gc:getStringWidth(tag) + 22
  col(gc, C.title2); fillRound(gc, x0 + 3, y0 + 2, tw, 11, 5)
  col(gc, C.accent); gc:fillArc(x0 + 5, y0 + 3, 9, 9, 0, 360)
  col(gc, C.titleTx); gc:drawString("p", x0 + 7, y0 + 1, "top")
  gc:drawString(tag, x0 + 17, y0 + 1, "top")
  -- elemento resaltado
  local hlNode, hlMem
  local hl = scr.hl
  if scr.key and scr.sel then hl = {scr.key, scr.sel} end
  if hl then
    local key, idx = hl[1], hl[2]
    local it = model[key] and model[key][idx]
    if it then
      if key == "nodes" then hlNode = idx
      elseif key == "mems" then hlMem = idx
      elseif key == "sups" or key == "nl" then hlNode = it.node
      elseif key == "ml" then hlMem = it.mem end
    end
  end
  local top, bot = y0 + 14, ybot - 14
  if st.G and #st.G.nodes > 0 then
    local tr
    if mode ~= "E" and st.res then
      if not st.dg then st.dg = Diagram.new(nil, st.res) end
      st.dg.q = mode
      st.dg.sel = hlMem or 0
      st.dg.labels = false
      st.dg:drawBody(gc, x0, top, w, bot - top, 18)
      tr = makeView(st.G.nodes, x0, top, w, bot - top, 18)
      if mode ~= "D" then
        -- valores extremos globales
        local qi = ({N = 1, V = 2, M = 3})[mode]
        local mx, mn = -1e300, 1e300
        for m in ipairs(st.res.mems) do
          local e = st.dg.data[m].ext[qi]
          mx = max(mx, e.mx); mn = min(mn, e.mn)
        end
        font(gc, 7, "b"); col(gc, C[mode])
        gc:drawString(fitText(gc, "max " .. fmtShort(mx) .. "   min " .. fmtShort(mn), w - 6), x0 + 4, bot - 13, "top")
      end
    else
      tr = makeView(st.G.nodes, x0, top, w, bot - top, 20)
      drawFrame(gc, st.G, tr, st.sup, {labels = true, sel = hlMem})
      drawLoads(gc, model, st.G, tr)
    end
    if hlNode and st.G.nodes[hlNode] then
      local x, y = tr(st.G.nodes[hlNode].x, st.G.nodes[hlNode].y)
      col(gc, C.load); gc:setPen("medium", "smooth")
      gc:drawArc(x - 6, y - 6, 12, 12, 0, 360)
      gc:setPen("thin", "smooth")
    end
  end
  -- estado
  font(gc, 7, "b")
  if st.err then
    local ls = wrapLines(gc, st.err, w - 14, 3)
    local hh = 11 * #ls + 3
    col(gc, {253, 232, 232}); fillRound(gc, x0 + 3, ybot - hh - 2, w - 6, hh, 4)
    col(gc, C.err); gc:fillRect(x0 + 3, ybot - hh - 2, 3, hh)
    for i, l in ipairs(ls) do gc:drawString(l, x0 + 9, ybot - hh - 2 + 11 * (i - 1), "top") end
  else
    local msg = fitText(gc, st.msg or "", w - 26)
    local mw = gc:getStringWidth(msg) + 20
    col(gc, st.good and {222, 245, 230} or {235, 238, 243}); fillRound(gc, x0 + 3, ybot - 14, mw, 12, 6)
    col(gc, st.good and C.ok or C.dim); gc:fillArc(x0 + 6, ybot - 12, 8, 8, 0, 360)
    gc:drawString(msg, x0 + 17, ybot - 15, "top")
  end
end

------------------------- Resultados en texto --------------------------------
local function valsLine(res)
  if res.ns == 0 then return nil end
  local t = {}
  for k, s in ipairs(res.syms) do t[#t + 1] = s .. "=" .. fmt(res.vals[k]) end
  return "Valores numericos: " .. table.concat(t, ", ")
end

local function numTag(res, L)
  if res.ns == 0 then return "" end
  return "  [" .. fmt(Eng.linNum(res, L)) .. "]"
end

function App.memberTable(m)
  local res = App.res
  local mb = res.mems[m]
  local truss = res.kind == "truss"
  local lines = {}
  local function add(s, st) lines[#lines + 1] = {s, st} end
  local ang = math.atan2(mb.s, mb.c) * 180 / PI
  add(string.format("Barra %d: nudo %d -> %d", m, mb.i, mb.j), "h")
  add("L = " .. fmt(mb.L) .. "   angulo = " .. fmt(ang, 2) .. " grados")
  add("EA=" .. fmt(mb.E * mb.A) .. (truss and "" or ("  EI=" .. fmt(mb.E * mb.I))) ..
      ((mb.ri or mb.rj) and not truss and ("  artic.: " .. (mb.ri and "i " or "") .. (mb.rj and "j" or "")) or ""))
  add("ESFUERZOS (literal):", "h")
  local function row(lbl, x, right)
    local N, V, M = Eng.internalLin(res, m, x, right)
    add(lbl .. "  N = " .. linStr(N) .. numTag(res, N))
    if not truss or #res.cases[0].L[m] > 0 then
      add("        V = " .. linStr(V) .. numTag(res, V))
      add("        M = " .. linStr(M) .. numTag(res, M))
    end
  end
  row("x=0  ", 0, true)
  if not truss or #res.cases[0].L[m] > 0 then row("x=L/2", mb.L / 2, false) end
  row("x=L  ", mb.L, false)
  local vl = valsLine(res)
  if vl then add(vl, "d") end
  add("TABLA NUMERICA  x | N | V | M", "h")
  local pts = Eng.samples(res, m, 10)
  local lastx
  for _, p in ipairs(pts) do
    local N, V, M = Eng.internal(res, m, p[1], p[2])
    add(string.format("%s%s | %s | %s | %s", fmt(p[1], 3), (lastx == p[1]) and "+" or "",
                      fmt(N, 3), fmt(V, 3), fmt(M, 3)))
    lastx = p[1]
  end
  local e = Eng.extremes(res, m)
  add("EXTREMOS", "h")
  local nm = {"N", "V", "M"}
  for q = 1, 3 do
    add(nm[q] .. "max = " .. fmt(e[q].mx) .. " en x=" .. fmt(e[q].xmx, 3) ..
        "   " .. nm[q] .. "min = " .. fmt(e[q].mn) .. " en x=" .. fmt(e[q].xmn, 3))
  end
  -- deformada
  local df = Eng.deflect(res, m, 40)
  local dmax, xd = 0, 0
  for _, p in ipairs(df) do
    local ux, uy = p[2], p[3]
    local v = -mb.s * ux + mb.c * uy
    if abs(v) > abs(dmax) then dmax = v; xd = p[1] end
  end
  add("Flecha max (y local) = " .. fmt(dmax) .. " en x=" .. fmt(xd, 3))
  push(TextView.new("Barra " .. m, lines))
end

function App.showEndForces()
  local res = App.res
  local lines = {}
  local function add(s, st) lines[#lines + 1] = {s, st} end
  local truss = res.kind == "truss"
  if truss then add("FUERZAS AXIALES (+ traccion)", "h")
  else add("ESFUERZOS EN EXTREMOS (N, V, M)", "h") end
  local vl = valsLine(res)
  if vl then add(vl, "d") end
  for m, mb in ipairs(res.mems) do
    if truss and #res.cases[0].L[m] == 0 then
      local N = Eng.internalLin(res, m, mb.L / 2, false)
      local nv = Eng.linNum(res, N)
      local tag = abs(nv) < 1e-9 and "(0)" or (nv > 0 and "(T)" or "(C)")
      add(string.format("b%d (%d-%d): N = %s%s %s", m, mb.i, mb.j, linStr(N), numTag(res, N), tag))
    else
      add(string.format("Barra %d (%d-%d)  L=%s", m, mb.i, mb.j, fmt(mb.L)), "h")
      for _, e in ipairs({{"i", 0, true}, {"j", mb.L, false}}) do
        local N, V, M = Eng.internalLin(res, m, e[2], e[3])
        add(" " .. e[1] .. ": N=" .. linStr(N) .. numTag(res, N))
        add("    V=" .. linStr(V) .. numTag(res, V))
        add("    M=" .. linStr(M) .. numTag(res, M))
      end
    end
  end
  push(TextView.new(truss and "Fuerzas axiales" or "Esfuerzos", lines))
end

function App.showReactions()
  local res = App.res
  local lines = {}
  local function add(s, st) lines[#lines + 1] = {s, st} end
  add("REACCIONES (ejes globales)", "h")
  local vl = valsLine(res)
  if vl then add(vl, "d") end
  local nm = {"Rx", "Ry", "Mz"}
  for _, n in ipairs(res.supNodes) do
    local S = res.sup[n]
    add("Nudo " .. n .. ":", "h")
    for a = 1, (res.kind == "truss" and 2 or 3) do
      local L = Eng.reacLin(res, n, a)
      add("  " .. nm[a] .. " = " .. linStr(L) .. numTag(res, L))
    end
    if S.ang ~= 0 then
      local Rx, Ry = Eng.reacLin(res, n, 1), Eng.reacLin(res, n, 2)
      local R1 = lin_add(lin_scale(Rx, S.ca), lin_scale(Ry, S.sa))
      local R2 = lin_add(lin_scale(Rx, -S.sa), lin_scale(Ry, S.ca))
      add("  (ejes girados " .. fmt(S.ang) .. " grados)", "d")
      add("  R1 = " .. linStr(R1) .. numTag(res, R1))
      add("  R2 = " .. linStr(R2) .. numTag(res, R2))
    end
  end
  push(TextView.new("Reacciones", lines))
end

function App.showDisplacements()
  local res = App.res
  local lines = {}
  local function add(s, st) lines[#lines + 1] = {s, st} end
  add("DESPLAZAMIENTOS NODALES", "h")
  add("(dx, dy globales; giro rad antihorario)", "d")
  local vl = valsLine(res)
  if vl then add(vl, "d") end
  local nm = {"dx", "dy", "giro"}
  for n = 1, #res.nodes do
    add("Nudo " .. n .. ":", "h")
    for a = 1, 3 do
      if not (a == 3 and (res.kind == "truss" or res.autofix[n])) then
        local L = Eng.dispLin(res, n, a)
        add("  " .. nm[a] .. " = " .. linStr(L) .. numTag(res, L))
      end
    end
  end
  push(TextView.new("Desplazamientos", lines))
end

function App.showEquilibrium()
  local res = App.res
  local lines = {}
  local function add(s, st) lines[#lines + 1] = {s, st} end
  add("EQUILIBRIO GLOBAL (cargas + reacciones)", "h")
  local tot = {}
  for k = 0, res.ns do tot[k] = {Eng.appliedTotals(res, k)} end
  local function L(i) return mkLin(res, function(k) return tot[k][i] end, "fs") end
  local function S(i) return mkLin(res, function(k) return tot[k][i] + tot[k][i + 3] end, "fs") end
  add("Cargas:  Fx=" .. linStr(L(1)))
  add("         Fy=" .. linStr(L(2)))
  add("         M0=" .. linStr(L(3)))
  add("Reacc.:  Rx=" .. linStr(L(4)))
  add("         Ry=" .. linStr(L(5)))
  add("         M0=" .. linStr(L(6)))
  add("Suma (debe ser 0):", "h")
  add("  SFx=" .. linStr(S(1)) .. "  SFy=" .. linStr(S(2)) .. "  SM=" .. linStr(S(3)))
  add("M0 = momento respecto al origen (0,0)", "d")
  push(TextView.new("Equilibrio", lines))
end

------------------------- Editores de datos ----------------------------------
local function M() return App.model end
local function isInt(s) local n = tonumber(s); return n and n == floor(n) and n or nil end
local function gresNow()
  local ok, vt = pcall(varTable, M())
  if not ok then return nil, vt end
  return makeGres(vt)
end
local function checkNum(s, label, allowEmpty)
  if isBlank(s) then
    if allowEmpty then return nil end
    return label .. ": falta valor"
  end
  local g, e = gresNow()
  if not g then return e end
  local ok, err = pcall(evalConst, s, g)
  if not ok then return label .. ": " .. tostring(err) end
end
local function checkLoad(s, label)
  if isBlank(s) then return nil end
  local ok, err = pcall(parseLin, s, function(n) local t = {}; t[n] = 1; return t end)
  if not ok then return label .. ": " .. tostring(err) end
end

local HELP_NUM = "numero o expresion (ej: 3, 2.5, 4*sqrt(2), a/2)"
local HELP_LOAD = "numero o letras: 10, P, -w, 2*P+5, w*3/2"

-- NUDOS
local function nodeStr(i)
  local n = M().nodes[i]
  return "N" .. i .. ":  X = " .. n.x .. "    Y = " .. n.y
end
local function editNode(i)
  local n = M().nodes[i] or {x = "0", y = "0"}
  local vals = {x = n.x, y = n.y}
  push(Form.new(i and ("Nudo " .. i) or ("Nuevo nudo " .. (#M().nodes + 1)), {
    {k = "x", label = "Coordenada X", t = "txt", help = HELP_NUM},
    {k = "y", label = "Coordenada Y", t = "txt", help = HELP_NUM},
  }, vals, function(v)
    local e = checkNum(v.x, "X") or checkNum(v.y, "Y")
    if e then return e end
    if i then M().nodes[i] = {x = v.x, y = v.y}
    else M().nodes[#M().nodes + 1] = {x = v.x, y = v.y} end
  end))
end
local function delNode(i)
  local m = M()
  table.remove(m.nodes, i)
  local function fix(n) if n > i then return n - 1 end return n end
  -- barras que usan el nudo
  local k = 1
  while k <= #m.mems do
    local mb = m.mems[k]
    if mb.i == i or mb.j == i then
      table.remove(m.mems, k)
      local q = 1
      while q <= #m.ml do
        if m.ml[q].mem == k then table.remove(m.ml, q)
        else if m.ml[q].mem > k then m.ml[q].mem = m.ml[q].mem - 1 end; q = q + 1 end
      end
    else
      mb.i = fix(mb.i); mb.j = fix(mb.j); k = k + 1
    end
  end
  for _, lst in ipairs({m.sups, m.nl}) do
    local q = 1
    while q <= #lst do
      if lst[q].node == i then table.remove(lst, q)
      else lst[q].node = fix(lst[q].node); q = q + 1 end
    end
  end
end

-- SECCIONES
local function secStr(i)
  local s = M().secs[i]
  if M().kind == "truss" then return "S" .. i .. ": E=" .. s.E .. "  A=" .. s.A .. "  a=" .. (s.al or "") end
  return "S" .. i .. ": E=" .. s.E .. " A=" .. s.A .. " I=" .. s.I .. " a=" .. (s.al or "") .. " h=" .. (s.h or "")
end
local function editSec(i)
  local s = M().secs[i] or {E = "1", A = M().kind == "truss" and "1" or "1e5", I = "1", al = "1e-5", h = "0.5"}
  local vals = {E = s.E, A = s.A, I = s.I, al = s.al or "", h = s.h or ""}
  local fr = M().kind ~= "truss"
  push(Form.new(i and ("Seccion " .. i) or "Nueva seccion", {
    {k = "E", label = "E (modulo elast.)", t = "txt", help = "use E=1 e I=1 para resultados en 1/EI"},
    {k = "A", label = "A (area)", t = "txt", help = fr and "A grande (1e5) = barra axialmente rigida" or HELP_NUM},
    {k = "I", label = "I (inercia)", t = "txt", help = HELP_NUM, show = function() return fr end},
    {k = "al", label = "alfa (coef. term.)", t = "txt", help = "solo para cargas de temperatura"},
    {k = "h", label = "h (altura secc.)", t = "txt", help = "solo para gradiente de temperatura",
     show = function() return fr end},
  }, vals, function(v)
    local e = checkNum(v.E, "E") or checkNum(v.A, "A") or (fr and checkNum(v.I, "I")) or
              checkNum(v.al, "alfa", true) or (fr and checkNum(v.h, "h", true))
    if e then return e end
    local rec = {E = v.E, A = v.A, I = fr and v.I or "1", al = v.al, h = v.h}
    if i then M().secs[i] = rec else M().secs[#M().secs + 1] = rec end
  end))
end
local function delSec(i)
  local m = M()
  if #m.secs <= 1 then return "Debe existir al menos una seccion" end
  table.remove(m.secs, i)
  for _, mb in ipairs(m.mems) do
    if mb.sec == i then mb.sec = 1 elseif mb.sec > i then mb.sec = mb.sec - 1 end
  end
end

-- BARRAS
local function memStr(k)
  local mb = M().mems[k]
  local s = "B" .. k .. ": " .. mb.i .. " -> " .. mb.j .. "   secc." .. (mb.sec or 1)
  if M().kind ~= "truss" and (mb.ri or mb.rj) then
    s = s .. "   rotula:" .. (mb.ri and " i" or "") .. (mb.rj and " j" or "")
  end
  return s
end
local function editMem(k)
  local mb = M().mems[k] or {i = #M().nodes - 1, j = #M().nodes, sec = 1}
  local vals = {i = tostring(mb.i or ""), j = tostring(mb.j or ""), sec = tostring(mb.sec or 1),
                ri = mb.ri or false, rj = mb.rj or false}
  local fr = M().kind ~= "truss"
  push(Form.new(k and ("Barra " .. k) or ("Nueva barra " .. (#M().mems + 1)), {
    {k = "i", label = "Nudo inicial i", t = "txt", help = "numero del nudo"},
    {k = "j", label = "Nudo final j", t = "txt", help = "numero del nudo"},
    {k = "sec", label = "Seccion No.", t = "txt", help = "numero de seccion/material"},
    {k = "ri", label = "Rotula en i", t = "chk", show = function() return fr end},
    {k = "rj", label = "Rotula en j", t = "chk", show = function() return fr end},
  }, vals, function(v)
    local i, j, s = isInt(v.i), isInt(v.j), isInt(v.sec)
    local nn = #M().nodes
    if not i or i < 1 or i > nn then return "Nudo i invalido (1.." .. nn .. ")" end
    if not j or j < 1 or j > nn then return "Nudo j invalido (1.." .. nn .. ")" end
    if i == j then return "i y j deben ser distintos" end
    if not s or s < 1 or s > #M().secs then return "Seccion invalida (1.." .. #M().secs .. ")" end
    local rec = {i = i, j = j, sec = s, ri = fr and v.ri or false, rj = fr and v.rj or false}
    if k then M().mems[k] = rec else M().mems[#M().mems + 1] = rec end
  end))
end
local function delMem(k)
  local m = M()
  table.remove(m.mems, k)
  local q = 1
  while q <= #m.ml do
    if m.ml[q].mem == k then table.remove(m.ml, q)
    else if m.ml[q].mem > k then m.ml[q].mem = m.ml[q].mem - 1 end; q = q + 1 end
  end
end

-- APOYOS
local function supStr(k)
  local s = M().sups[k]
  local opts = M().kind == "truss" and SUP_TRUSS or SUP_FRAME
  local t = "A" .. k .. ": nudo " .. s.node .. "  " .. (opts[s.tipo or #opts] or "?")
  local r = (s.rx and "X" or "") .. (s.ry and "Y" or "") .. (s.rz and M().kind ~= "truss" and "G" or "")
  t = t .. " [" .. r .. "]"
  if not isBlank(s.ang) and tonumber(s.ang) ~= 0 then t = t .. " ang=" .. s.ang end
  if not isBlank(s.kx) or not isBlank(s.ky) or not isBlank(s.kr) then t = t .. " resorte" end
  if not isBlank(s.dx) or not isBlank(s.dy) or not isBlank(s.dr) then t = t .. " asent." end
  return t
end
local function editSup(k)
  local fr = M().kind ~= "truss"
  local opts = fr and SUP_FRAME or SUP_TRUSS
  local pres = fr and SUP_FRAME_R or SUP_TRUSS_R
  local s = M().sups[k] or {node = "", tipo = 1, rx = true, ry = true, rz = fr}
  local vals = {node = tostring(s.node), tipo = s.tipo or #opts, rx = s.rx or false, ry = s.ry or false,
                rz = s.rz or false, kx = s.kx or "", ky = s.ky or "", kr = s.kr or "",
                dx = s.dx or "", dy = s.dy or "", dr = s.dr or "", ang = s.ang or ""}
  local function preset(v)
    local p = pres[v.tipo]
    if p then v.rx, v.ry, v.rz = p[1], p[2], p[3] end
  end
  local function custom(v) v.tipo = #opts end
  push(Form.new(k and ("Apoyo " .. k) or "Nuevo apoyo", {
    {k = "node", label = "Nudo", t = "txt", help = "numero del nudo apoyado"},
    {k = "tipo", label = "Tipo", t = "opt", opts = opts, change = preset},
    {k = "rx", label = "Restringe X (1)", t = "chk", change = custom},
    {k = "ry", label = "Restringe Y (2)", t = "chk", change = custom},
    {k = "rz", label = "Restringe giro", t = "chk", change = custom, show = function() return fr end},
    {k = "ang", label = "Angulo apoyo (gr)", t = "txt", help = "gira ejes 1-2 del apoyo (antihorario). 0 = normal"},
    {k = "kx", label = "Resorte k1 (X)", t = "txt", help = "rigidez resorte (vacio = sin resorte)"},
    {k = "ky", label = "Resorte k2 (Y)", t = "txt", help = "rigidez resorte (vacio = sin resorte)"},
    {k = "kr", label = "Resorte giro", t = "txt", help = "rigidez rotacional", show = function() return fr end},
    {k = "dx", label = "Asentam. d1 (X)", t = "txt", help = "desplaz. impuesto (requiere restriccion) " .. HELP_LOAD},
    {k = "dy", label = "Asentam. d2 (Y)", t = "txt", help = "desplaz. impuesto, ej: -0.01 o -d"},
    {k = "dr", label = "Giro impuesto", t = "txt", help = "rotacion impuesta (rad)", show = function() return fr end},
  }, vals, function(v)
    local n = isInt(v.node)
    if not n or n < 1 or n > #M().nodes then return "Nudo invalido (1.." .. #M().nodes .. ")" end
    local e = checkNum(v.ang, "angulo", true) or checkNum(v.kx, "k1", true) or checkNum(v.ky, "k2", true) or
              checkNum(v.kr, "kgiro", true) or checkLoad(v.dx, "d1") or checkLoad(v.dy, "d2") or checkLoad(v.dr, "giro")
    if e then return e end
    if not (v.rx or v.ry or v.rz) and isBlank(v.kx) and isBlank(v.ky) and isBlank(v.kr) then
      return "El apoyo no restringe nada"
    end
    local rec = {node = n, tipo = v.tipo, rx = v.rx, ry = v.ry, rz = fr and v.rz or false,
                 kx = v.kx, ky = v.ky, kr = v.kr, dx = v.dx, dy = v.dy, dr = v.dr, ang = v.ang}
    if k then M().sups[k] = rec else M().sups[#M().sups + 1] = rec end
  end))
end

-- CARGAS NODALES
local function nlStr(k)
  local l = M().nl[k]
  local s = "Q" .. k .. ": nudo " .. l.node
  if not isBlank(l.fx) then s = s .. "  Fx=" .. l.fx end
  if not isBlank(l.fy) then s = s .. "  Fy=" .. l.fy end
  if M().kind ~= "truss" and not isBlank(l.m) then s = s .. "  M=" .. l.m end
  return s
end
local function editNL(k)
  local fr = M().kind ~= "truss"
  local l = M().nl[k] or {node = "", fx = "", fy = "", m = ""}
  local vals = {node = tostring(l.node), fx = l.fx or "", fy = l.fy or "", m = l.m or ""}
  push(Form.new(k and ("Carga nodal " .. k) or "Nueva carga nodal", {
    {k = "node", label = "Nudo", t = "txt", help = "numero del nudo"},
    {k = "fx", label = "Fx (+ derecha)", t = "txt", help = HELP_LOAD},
    {k = "fy", label = "Fy (+ arriba)", t = "txt", help = HELP_LOAD .. "  (abajo: -P)"},
    {k = "m", label = "Momento (+ antihor.)", t = "txt", help = HELP_LOAD, show = function() return fr end},
  }, vals, function(v)
    local n = isInt(v.node)
    if not n or n < 1 or n > #M().nodes then return "Nudo invalido (1.." .. #M().nodes .. ")" end
    local e = checkLoad(v.fx, "Fx") or checkLoad(v.fy, "Fy") or checkLoad(v.m, "M")
    if e then return e end
    local rec = {node = n, fx = v.fx, fy = v.fy, m = fr and v.m or ""}
    if k then M().nl[k] = rec else M().nl[#M().nl + 1] = rec end
  end))
end

-- CARGAS EN BARRAS
local DIRS_SHORT = {"perp", "axial", "GX", "GY", "grav", "GXproy", "GYproy", "grav.proy"}
local function mlStr(k)
  local l = M().ml[k]
  local s = "C" .. k .. ": b" .. l.mem .. " " .. TIPOS_CARGA[l.t]
  if l.t == 1 then s = s .. " " .. l.v1 .. " " .. DIRS_SHORT[l.dir] .. " a=" .. (isBlank(l.a) and "L/2" or l.a)
  elseif l.t == 2 then s = s .. " " .. l.v1 .. " a=" .. (isBlank(l.a) and "L/2" or l.a)
  elseif l.t == 3 then
    s = s .. " " .. l.v1 .. (isBlank(l.v2) and "" or (".." .. l.v2)) .. " " .. DIRS_SHORT[l.dir]
    if not (isBlank(l.a) or l.a == "0") or not (isBlank(l.b) or l.b == "L") then
      s = s .. " [" .. (isBlank(l.a) and "0" or l.a) .. "," .. (isBlank(l.b) and "L" or l.b) .. "]"
    end
  elseif l.t == 4 then s = s .. " dT=" .. l.v1 .. (isBlank(l.v2) and "" or (" grad=" .. l.v2))
  else s = s .. " dL=" .. l.v1 end
  return s
end
local function editML(k)
  local l = M().ml[k] or {mem = "", t = 3, dir = 5, v1 = "", v2 = "", a = "", b = ""}
  local vals = {mem = tostring(l.mem), t = l.t, dir = l.dir or 5, v1 = l.v1 or "", v2 = l.v2 or "",
                a = l.a or "", b = l.b or ""}
  local fr = M().kind ~= "truss"
  local tipos = fr and TIPOS_CARGA or TIPOS_CARGA
  local function T(n) return function(v) return v.t == n end end
  push(Form.new(k and ("Carga en barra " .. k) or "Nueva carga en barra", {
    {k = "mem", label = "Barra No.", t = "txt", help = "numero de la barra"},
    {k = "t", label = "Tipo", t = "opt", opts = tipos},
    {k = "dir", label = "Direccion", t = "opt", opts = DIRS,
     help = "proy.: por metro de proyeccion horizontal/vertical",
     show = function(v) return v.t == 1 or v.t == 3 end},
    {k = "v1", t = "txt", help = HELP_LOAD, label = function(v)
      return ({"Valor P", "Momento M (+antih.)", "w inicial (en a)", "dT uniforme", "dL (alargamiento +)"})[v.t] end},
    {k = "v2", t = "txt", show = function(v) return v.t == 3 or v.t == 4 end,
     help = "vacio = igual a w inicial (uniforme); 0 = triangular",
     label = function(v) return v.t == 3 and "w final (en b)" or "dT inf - dT sup" end},
    {k = "a", t = "txt", show = function(v) return v.t <= 3 end,
     help = "distancia desde nudo i. Puede usar L (ej: L/3). vacio: " ,
     label = function(v) return v.t == 3 and "Desde a (vacio=0)" or "Posicion a (vacio=L/2)" end},
    {k = "b", t = "txt", show = T(3), label = "Hasta b (vacio=L)", help = "distancia desde nudo i. Puede usar L"},
  }, vals, function(v)
    local mm = isInt(v.mem)
    if not mm or mm < 1 or mm > #M().mems then return "Barra invalida (1.." .. #M().mems .. ")" end
    if isBlank(v.v1) then return "Falta el valor de la carga" end
    local e = checkLoad(v.v1, "valor") or checkLoad(v.v2, "valor 2")
    if e then return e end
    local rec = {mem = mm, t = v.t, dir = v.dir, v1 = v.v1, v2 = (v.t == 3 or v.t == 4) and v.v2 or "",
                 a = v.t <= 3 and v.a or "", b = v.t == 3 and v.b or ""}
    if k then M().ml[k] = rec else M().ml[#M().ml + 1] = rec end
  end))
end

-- VARIABLES
local RESERVED = {L = true, pi = true, sqrt = true, sin = true, cos = true, tan = true, abs = true, ln = true, exp = true}
local function varStr(k)
  local v = M().vars[k]
  local s = v.name .. " = " .. (isBlank(v.val) and "(sin valor)" or v.val)
  if v.keep then s = s .. "   [letra en cargas]" else s = s .. "   [se reemplaza]" end
  return s
end
local function editVar(k)
  local v0 = M().vars[k] or {name = "", val = "", keep = true}
  local vals = {name = v0.name, val = v0.val or "", keep = v0.keep}
  push(Form.new(k and ("Variable " .. v0.name) or "Nueva variable", {
    {k = "name", label = "Nombre (letra)", t = "txt", help = "ej: P, w, q1, H  (L esta reservada)"},
    {k = "val", label = "Valor numerico", t = "txt", help = "para graficos / geometria. Vacio = 1 en graficos"},
    {k = "keep", label = "Letra en resultados", t = "chk"},
  }, vals, function(v)
    local nm = v.name:match("^%s*([%a_][%w_]*)%s*$")
    if not nm then return "Nombre invalido (letras/numeros, empieza con letra)" end
    if RESERVED[nm] then return "Nombre reservado: " .. nm end
    for q, o in ipairs(M().vars) do if o.name == nm and q ~= k then return "Ya existe " .. nm end end
    if not isBlank(v.val) then
      local ok, e = pcall(evalConst, v.val, function(n)
        for _, o in ipairs(M().vars) do
          if o.name == n and o.name ~= nm and not isBlank(o.val) then
            return lin_const(evalConst(o.val, function() fail("anidado") end))
          end
        end
        fail("'" .. n .. "' sin valor")
      end)
      if not ok then return "Valor: " .. tostring(e) end
    end
    local rec = {name = nm, val = v.val, keep = v.keep}
    if k then M().vars[k] = rec else M().vars[#M().vars + 1] = rec end
  end))
end
local function delVar(k) table.remove(M().vars, k) end

local function listScreen(title, key, strf, editf, delf)
  local function mark(i)
    local t = top()
    if t and t.ok then t.hl = {key, i or (#M()[key] + 1)} end
  end
  return List.new{title = title, key = key,
    items = function()
      local t = {}
      for i = 1, #M()[key] do t[i] = strf(i) end
      return t
    end,
    add = function() editf(nil); mark(nil) end,
    edit = function(i) editf(i); mark(i) end,
    del = delf or function(i) table.remove(M()[key], i) end}
end

------------------------- Menus principales ----------------------------------
local ensureSolved

function App.results()
  if not ensureSolved() then return end
  local res = App.res
  local truss = res.kind == "truss"
  local items = {
    {label = truss and "Grafico fuerzas axiales" or "Diagramas M, V, N, deformada", icon = "diagram",
     action = function() if ensureSolved() then push(Diagram.new()) end end},
    {label = "Reacciones", icon = "reaction", action = function() if ensureSolved() then App.showReactions() end end},
    {label = "Desplazamientos", icon = "disp", action = function() if ensureSolved() then App.showDisplacements() end end},
    {label = truss and "Fuerzas axiales (lista)" or "Esfuerzos en extremos", icon = "forces",
     action = function() if ensureSolved() then App.showEndForces() end end},
    {label = "Tabla por barra", icon = "table", action = function()
      if not ensureSolved() then return end
      local its = {}
      for m, mb in ipairs(App.res.mems) do
        its[#its + 1] = {label = "Barra " .. m .. " (" .. mb.i .. "-" .. mb.j .. ")",
                         icon = "bar", action = function() App.memberTable(m) end}
      end
      push(Menu.new("Elija barra", its))
    end},
    {label = "Verificar equilibrio", icon = "check", action = function() if ensureSolved() then App.showEquilibrium() end end},
  }
  if #res.warn > 0 then
    items[#items + 1] = {label = "Avisos (" .. #res.warn .. ")", icon = "warn", action = function() msgBox("Avisos", res.warn) end}
  end
  local lt = res.ns > 0 and (" letras: " .. table.concat(res.syms, ",")) or ""
  push(Menu.new("Resultados" .. lt, items))
end

ensureSolved = function()
  if App.res then return true end
  if not App.model then msgBox("Aviso", "No hay modelo. Cree uno nuevo."); return false end
  local ok, r = pcall(Eng.solve, App.model)
  if not ok then
    msgBox("No se pudo resolver", {{"ERROR:", "e"}, tostring(r)})
    return false
  end
  App.res = r
  return true
end

function App.modelMenu()
  local fr = function() return M().kind ~= "truss" end
  local items = {
    {icon = "node", label = function() return "Nudos (" .. #M().nodes .. ")" end,
     action = function() push(listScreen("Nudos", "nodes", nodeStr, editNode, delNode)) end},
    {icon = "section", label = function() return "Secciones / material (" .. #M().secs .. ")" end,
     action = function() push(listScreen("Secciones", "secs", secStr, editSec, delSec)) end},
    {icon = "bar", label = function() return "Barras (" .. #M().mems .. ")" end,
     action = function() push(listScreen("Barras", "mems", memStr, editMem, delMem)) end},
    {icon = "support", label = function() return "Apoyos (" .. #M().sups .. ")" end,
     action = function() push(listScreen("Apoyos", "sups", supStr, editSup)) end},
    {icon = "nload", label = function() return "Cargas en nudos (" .. #M().nl .. ")" end,
     action = function() push(listScreen("Cargas en nudos", "nl", nlStr, editNL)) end},
    {icon = "mload", label = function() return "Cargas en barras (" .. #M().ml .. ")" end,
     action = function() push(listScreen("Cargas en barras", "ml", mlStr, editML)) end},
    {icon = "var", label = function() return "Variables / letras (" .. #M().vars .. ")" end,
     action = function() push(listScreen("Variables", "vars", varStr, editVar, delVar)) end},
    {label = "Ver estructura", icon = "eye", action = function()
      local s = StructView.new(); if s then push(s) end end},
    {label = "RESOLVER y ver resultados", icon = "play", action = function() App.results() end},
  }
  local mm = Menu.new(function() return fr() and "Marco / Portico / Viga" or "Enrejado / Armadura" end, items,
                      "enter:abrir|p:cambia vista|esc:volver")
  mm.split = true
  push(mm)
end

local function help()
  msgBox("Ayuda", {
    {"ESTRUCTURAS 2D - metodo de rigidez", "h"},
    "1) Cree nudos (X,Y). 2) Secciones (E,A,I). 3) Barras (nudo i -> j). 4) Apoyos. 5) Cargas. 6) Resolver.",
    {"UNIDADES", "h"},
    "Use unidades consistentes (kN, m). Con E=1 e I=1 los desplazamientos quedan multiplicados por EI.",
    "A grande (1e5) aproxima barras axialmente rigidas.",
    {"CARGAS CON LETRAS", "h"},
    "Escriba 10, P, -w, 2*P+5, w*L/2 ... El resultado sale como combinacion lineal: M = 3.5*P - 12*w.",
    "En Variables puede dar valor a una letra: si 'Letra en resultados' = Si, se mantiene literal (el valor solo se usa en graficos); si = No, se reemplaza por el numero.",
    "Letras con valor tambien sirven para geometria (ej: X = a/2).",
    {"CONVENCIONES", "h"},
    "Global: X derecha, Y arriba, giro antihorario +. Gravedad: use direccion 'Gravedad (-Y)' con valor positivo.",
    "Local: x de i a j, y a 90 grados antihorario. 'Perpend. local' + es hacia +y local.",
    "N + traccion. M + tracciona la fibra inferior (lado -y local); el diagrama M se dibuja del lado traccionado.",
    "Cargas distribuidas: w inicial en a y w final en b (trapecio/triangulo). Posiciones desde el nudo i; puede usar L.",
    "Temperatura: dT uniforme y gradiente (T inferior - T superior) con alfa y h de la seccion.",
    "Error de fabricacion: dL + si la barra es mas larga.",
    "Rotula en i/j: libera el momento en ese extremo (articulacion interna).",
    "Apoyo inclinado: el angulo gira los ejes 1-2 del apoyo. Rodillo inclinado = restringe solo el eje 2.",
    {"PANTALLA DIVIDIDA", "h"},
    "Al editar datos, la mitad derecha muestra la estructura en vivo (incluso mientras escribe), resalta el elemento seleccionado e indica si es isostatica, hiperestatica (grado) o inestable.",
    "Tecla p (en menus y listas): cambia la vista entre Estructura, M, V, N y Deformada, que se recalculan solos.",
    {"TECLAS", "h"},
    "Menus: flechas, enter, numeros. Listas: + nuevo, del borrar. Formularios: flechas, enter siguiente campo, ACEPTAR guarda.",
    "Diagramas: <- -> barra, m/v/n/d tipo, flechas arriba/abajo escala, tab etiquetas, enter tabla de la barra.",
    "Tecla menu: acceso rapido a todas las opciones.",
  })
end

local function newModelConfirm(kind)
  local function go()
    App.model = newModel(kind); App.res = nil
    while #stack > 1 do stack[#stack] = nil end
    App.modelMenu()
  end
  if App.model and #App.model.nodes > 0 then
    push(Menu.new("Borrar el modelo actual?", {
      {label = "Si, crear nuevo", action = function() pop(); go() end},
      {label = "No, volver", action = function() pop() end}}))
  else
    go()
  end
end

local function examplesMenu()
  local its = {}
  for i, ex in ipairs(Examples) do
    its[#its + 1] = {label = ex[1], icon = (i >= 4) and "truss" or "frame", action = function()
      App.model = ex[2](); App.res = nil
      while #stack > 1 do stack[#stack] = nil end
      App.modelMenu()
    end}
  end
  push(Menu.new("Ejemplos", its))
end

------------------------- Portada ---------------------------------------------
local Splash = {}; Splash.__index = Splash
function Splash.new()
  local sp = setmetatable({t = 0}, Splash)
  local ok, res = pcall(Eng.solve, Examples[3][2]())
  if ok then
    sp.dg = Diagram.new("M", res)
    sp.dg.sel = 0; sp.dg.labels = false; sp.dg.zoom = 0
  end
  return sp
end
function Splash:tick()
  self.t = self.t + 1
  if self.t == 32 and timer then pcall(timer.start, 0.5) end
end
function Splash:paint(gc)
  if not self.started then
    self.started = true
    if timer then pcall(timer.start, 0.05) end
  end
  local t = self.t
  grad(gc, 0, 0, W, H, C.title2, C.title)
  -- cuadricula tenue
  gc:setColorRGB(38, 76, 128)
  for x = 8, W, 16 do gc:fillRect(x, 0, 1, H) end
  for y = 8, H, 16 do gc:fillRect(0, y, W, 1) end
  -- titulo
  font(gc, 24, "b")
  local title = "ESTRUCTURAS 2D"
  local tw = gc:getStringWidth(title)
  local tx = (W - tw) / 2 + 12
  drawLogo(gc, tx - 30, 12, 22, C.accent)
  col(gc, C.titleTx); font(gc, 24, "b")
  gc:drawString(title, tx, 4, "top")
  local lw = floor(min(1, t / 12) * (tw + 30))
  col(gc, C.accent); gc:fillRect((W - lw) / 2, 38, lw, 2)
  font(gc, 9); gc:setColorRGB(190, 210, 238)
  local sub = "Analisis matricial de marcos y enrejados"
  gc:drawString(sub, (W - gc:getStringWidth(sub)) / 2, 42, "top")
  -- tarjeta con un portico real resuelto
  local cx, cy, cw, ch = 60, 60, W - 120, 104
  gc:setColorRGB(8, 24, 48); fillRound(gc, cx + 2, cy + 3, cw, ch, 8)
  col(gc, C.card); fillRound(gc, cx, cy, cw, ch, 8)
  col(gc, C.grid)
  for x = cx + 12, cx + cw - 6, 12 do gc:fillRect(x, cy + 4, 1, ch - 8) end
  for y = cy + 12, cy + ch - 6, 12 do gc:fillRect(cx + 4, y, cw - 8, 1) end
  if self.dg then
    self.dg.zoom = min(1, max(0, (t - 6) / 16)) * 1.1
    if self.dg.zoom > 0 then
      self.dg:drawBody(gc, cx + 4, cy + 12, cw - 8, ch - 16, 14)
    else
      local tr = makeView(self.dg.res.nodes, cx + 4, cy + 12, cw - 8, ch - 16, 14)
      drawFrame(gc, self.dg.G, tr, self.dg.sup, {})
    end
  end
  col(gc, C.M); fillRound(gc, cx + 6, cy + 5, 58, 11, 5)
  col(gc, C.titleTx); font(gc, 7, "b"); gc:drawString("Momento M", cx + 11, cy + 3, "top")
  -- rasgos
  font(gc, 7, "b"); gc:setColorRGB(190, 210, 238)
  local feats = {"Hiperestaticos", "Enrejados", "Cargas con letras", "N  V  M"}
  local fx = 8
  for i, f in ipairs(feats) do
    local fw = gc:getStringWidth(f) + 10
    fx = fx + fw + 4
  end
  fx = (W - (fx - 8 - 4)) / 2
  for _, f in ipairs(feats) do
    local fw = gc:getStringWidth(f) + 10
    gc:setColorRGB(28, 64, 112); fillRound(gc, fx, 170, fw, 12, 6)
    gc:setColorRGB(200, 218, 242); gc:drawString(f, fx + 5, 169, "top")
    fx = fx + fw + 4
  end
  -- invitacion
  if t < 32 or t % 2 == 0 then
    col(gc, C.accent); font(gc, 10, "b")
    local m = "Presione ENTER para comenzar"
    gc:drawString(m, (W - gc:getStringWidth(m)) / 2, 185, "top")
  end
  font(gc, 7); gc:setColorRGB(130, 158, 196)
  local v = "TI-Nspire CX II  -  Metodo de rigidez  -  v2.0"
  gc:drawString(v, (W - gc:getStringWidth(v)) / 2, 200, "top")
end
function Splash:leave()
  if timer then pcall(timer.stop) end
  pop()
end
Splash.enter = Splash.leave
Splash.esc = Splash.leave
Splash.tab = Splash.leave
Splash.click = Splash.leave
Splash.arrow = Splash.leave
Splash.char = Splash.leave
Splash.back = Splash.leave

local mainMenu = Menu.new("ESTRUCTURAS 2D", {
  {label = "Nuevo marco / portico / viga", icon = "frame", action = function() newModelConfirm("frame") end},
  {label = "Nuevo enrejado / armadura", icon = "truss", action = function() newModelConfirm("truss") end},
  {icon = "go", label = function()
     if not App.model then return "Continuar (sin modelo)" end
     return "Continuar modelo actual (" .. (App.model.kind == "truss" and "enrejado" or "marco") .. ")"
   end, action = function()
     if App.model then App.modelMenu() else msgBox("Aviso", "No hay modelo. Cree uno nuevo o cargue un ejemplo.") end
   end},
  {label = "Ejemplos", icon = "book", action = examplesMenu},
  {label = "Ayuda", icon = "help", action = help},
  {label = "Acerca de / portada", icon = "info", action = function() push(Splash.new()) end},
}, "enter:abrir|1-6:elegir|menu:accesos rapidos")
function mainMenu:esc() end
stack[1] = mainMenu
stack[2] = Splash.new()

------------------------- Eventos TI-Nspire ----------------------------------
local function showError(e)
  msgBox("Error", {{"Error interno:", "e"}, tostring(e)})
end

local function safe(f)
  local ok, e = pcall(f)
  if not ok then pcall(showError, e) end
  inval()
end

local function dispatch(name, ...)
  local args = {...}
  safe(function()
    local s = top()
    if s[name] then s[name](s, unpack(args)) end
  end)
end

function on.paint(gc)
  local ok, e = pcall(function()
    col(gc, {255, 255, 255}); gc:fillRect(0, 0, W, H)
    top():paint(gc)
  end)
  if not ok then
    col(gc, C.err); font(gc, 9)
    gc:drawString("Error: " .. tostring(e), 2, 30, "top")
    stack = {mainMenu}
  end
end
function on.resize(w, h) W, H = w, h end
function on.timer()
  local s = top()
  if s and s.tick then s:tick() elseif timer then timer.stop() end
  inval()
end
function on.arrowKey(k) dispatch("arrow", k) end
function on.enterKey() dispatch("enter") end
function on.returnKey() dispatch("enter") end
function on.escapeKey() dispatch("esc") end
function on.tabKey() dispatch("tab") end
function on.backtabKey() dispatch("arrow", "up") end
function on.backspaceKey() dispatch("back") end
function on.deleteKey() dispatch("back") end
function on.clearKey() dispatch("clear") end
function on.charIn(ch) dispatch("char", ch) end
function on.mouseDown(x, y) dispatch("click", x, y) end
function on.save() return {model = App.model} end
function on.restore(st)
  if type(st) == "table" and type(st.model) == "table" then App.model = st.model end
end
function on.construction()
  if platform.window then W, H = platform.window:width(), platform.window:height() end
end

-- menu de herramientas (tecla menu)
local function jump(f) return function() safe(function()
  if not App.model then msgBox("Aviso", "Primero cree un modelo"); return end
  f()
end) end end
local function openList(title, key, s, e, d) return jump(function() push(listScreen(title, key, s, e, d)) end) end
if toolpalette then
  pcall(function()
    toolpalette.register({
      {"Archivo",
        {"Nuevo marco", function() safe(function() newModelConfirm("frame") end) end},
        {"Nuevo enrejado", function() safe(function() newModelConfirm("truss") end) end},
        {"Ejemplos", function() safe(examplesMenu) end},
        {"Menu principal", function() safe(function() while #stack > 1 do stack[#stack] = nil end end) end},
        {"Ayuda", function() safe(help) end}},
      {"Datos",
        {"Nudos", openList("Nudos", "nodes", nodeStr, editNode, delNode)},
        {"Secciones", openList("Secciones", "secs", secStr, editSec, delSec)},
        {"Barras", openList("Barras", "mems", memStr, editMem, delMem)},
        {"Apoyos", openList("Apoyos", "sups", supStr, editSup)},
        {"Cargas en nudos", openList("Cargas en nudos", "nl", nlStr, editNL)},
        {"Cargas en barras", openList("Cargas en barras", "ml", mlStr, editML)},
        {"Variables", openList("Variables", "vars", varStr, editVar, delVar)}},
      {"Calcular",
        {"Resolver", jump(App.results)},
        {"Ver estructura", jump(function() local s = StructView.new(); if s then push(s) end end)}},
      {"Resultados",
        {"Diagramas", jump(function() if ensureSolved() then push(Diagram.new()) end end)},
        {"Reacciones", jump(function() if ensureSolved() then App.showReactions() end end)},
        {"Desplazamientos", jump(function() if ensureSolved() then App.showDisplacements() end end)},
        {"Esfuerzos", jump(function() if ensureSolved() then App.showEndForces() end end)},
        {"Equilibrio", jump(function() if ensureSolved() then App.showEquilibrium() end end)}},
    })
  end)
end

-- acceso para pruebas fuera de la calculadora
if _TESTING then
  _TESTING.Eng = Eng; _TESTING.App = App; _TESTING.parseLin = parseLin; _TESTING.linStr = linStr
  _TESTING.newModel = newModel; _TESTING.Examples = Examples; _TESTING.stack = function() return stack end
  _TESTING.Diagram = Diagram; _TESTING.StructView = StructView; _TESTING.Splash = Splash
end
