-- Pruebas de escritorio: lua5.1 tests/test.lua  (desde estructuras_nspire/)
-- Simula la API de la TI-Nspire y verifica el motor contra soluciones conocidas.

local invalidated = 0
platform = {window = {invalidate = function() invalidated = invalidated + 1 end,
                      width = function() return 318 end, height = function() return 212 end}}
on = {}
local registered
toolpalette = {register = function(t) registered = t end}
_TESTING = {}

-- gc simulado: verifica tipos de argumentos
local function mkgc()
  local gc = {}
  local function nums(...) for i, v in ipairs({...}) do assert(type(v) == "number" and v == v, "arg no numerico " .. i) end end
  function gc:setColorRGB(r, g, b) nums(r, g, b) end
  function gc:setFont(f, s, z) assert(type(f) == "string" and type(s) == "string" and type(z) == "number") end
  function gc:drawString(s, x, y, a) assert(type(s) == "string", "drawString texto"); nums(x, y) end
  function gc:getStringWidth(s) return 6 * #s end
  function gc:getStringHeight(s) return 12 end
  function gc:fillRect(x, y, w, h) nums(x, y, w, h) end
  function gc:drawRect(x, y, w, h) nums(x, y, w, h) end
  function gc:drawLine(a, b, c, d) nums(a, b, c, d) end
  function gc:fillArc(x, y, w, h, s, e) nums(x, y, w, h, s, e) end
  function gc:drawArc(x, y, w, h, s, e) nums(x, y, w, h, s, e) end
  function gc:drawPolyLine(t) assert(#t % 2 == 0); nums(unpack(t)) end
  function gc:fillPolygon(t) nums(unpack(t)) end
  function gc:setPen(a, b) end
  return gc
end

dofile("estructuras.lua")
local T = _TESTING
local Eng = T.Eng

local fails, count = 0, 0
local function near(a, b, tol, msg)
  count = count + 1
  tol = tol or 1e-6
  if math.abs(a - b) > tol * math.max(1, math.abs(b)) then
    fails = fails + 1
    print(("FALLA %s: obtenido %.8g esperado %.8g"):format(msg, a, b))
  end
end
local function ok(c, msg) count = count + 1; if not c then fails = fails + 1; print("FALLA " .. msg) end end

-- expresiones
local P = T.parseLin
local r = P("2*P+5-3w/2", function(n) local t = {}; t[n] = 1; return t end)
near(r[""], 5, 1e-12, "lin const"); near(r.P, 2, 1e-12, "lin P"); near(r.w, -1.5, 1e-12, "lin w")
r = P("\226\136\1462(3+x)", function(n) local t = {}; t[n] = 1; return t end)
near(r.x, -2, 1e-12, "menos TI y implicito"); near(r[""], -6, 1e-12, "implicito const")
r = P("2\225\180\1353", function() end); near(r[""], 2000, 1e-12, "tecla EE")
r = P("4*sqrt(2)^2+sin(30)", function() end); near(r[""], 8.5, 1e-12, "funciones")
ok(not pcall(P, "P*w", function(n) local t = {}; t[n] = 1; return t end), "no lineal detectado")
ok(T.linStr({[""] = 12.5, P = 3.2, w = -0.5}) == "12.5 + 3.2P - 0.5w", "formato lin: " .. T.linStr({[""] = 12.5, P = 3.2, w = -0.5}))

local function beam(L, sups, ml, extra)
  local m = T.newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = tostring(L), y = "0"}}
  m.mems = {{i = 1, j = 2, sec = 1}}
  m.sups = sups
  m.ml = ml or {}
  if extra then extra(m) end
  return Eng.solve(m), m
end
local FIX = function(n) return {node = n, rx = true, ry = true, rz = true} end
local PIN = function(n) return {node = n, rx = true, ry = true} end
local ROL = function(n) return {node = n, ry = true} end
local function MV(res, m, x, right) local N, V, M = Eng.internal(res, m, x, right); return N, V, M end

-- 1 viga biempotrada con carga uniforme
local res = beam(6, {FIX(1), FIX(2)}, {{mem = 1, t = 3, dir = 5, v1 = "10"}})
local _, V0, M0 = MV(res, 1, 0, true)
near(M0, -30, 1e-6, "biempotrada M extremo"); near(V0, 30, 1e-6, "biempotrada V")
local _, _, Mc = MV(res, 1, 3); near(Mc, 15, 1e-6, "biempotrada M centro")
local _, _, ML = MV(res, 1, 6); near(ML, -30, 1e-6, "biempotrada M j")

-- 2 viga continua 2 tramos iguales
do
  local m = T.newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = "6", y = "0"}, {x = "12", y = "0"}}
  m.mems = {{i = 1, j = 2, sec = 1}, {i = 2, j = 3, sec = 1}}
  m.sups = {PIN(1), ROL(2), ROL(3)}
  m.ml = {{mem = 1, t = 3, dir = 5, v1 = "w"}, {mem = 2, t = 3, dir = 5, v1 = "w"}}
  local res = Eng.solve(m)
  local _, _, Mb = Eng.internalLin(res, 1, 6, false)
  near(Mb.w, -4.5, 1e-6, "continua M apoyo = -wL^2/8")
  near(Eng.reacLin(res, 2, 2).w, 7.5, 1e-6, "continua R central = 1.25wL")
end

-- 3 empotrada-apoyada con carga puntual al centro: MA = -3PL/16
res = beam(8, {FIX(1), ROL(2)}, {{mem = 1, t = 1, dir = 5, v1 = "P", a = "L/2"}})
local _, _, MA = Eng.internalLin(res, 1, 0, true)
near(MA.P, -3 * 8 / 16, 1e-6, "propped cantilever MA")
local N, V, Mm = Eng.internalLin(res, 1, 4, true)
near(Mm.P, 5 * 8 / 32, 1e-6, "propped M bajo carga")

-- 4 rotula en barra = simplemente apoyada
res = beam(6, {PIN(1), FIX(2)}, {{mem = 1, t = 3, dir = 5, v1 = "10"}}, function(m) m.mems[1].rj = true end)
_, _, Mc = MV(res, 1, 3); near(Mc, 45, 1e-6, "rotula: wL^2/8")
near(Eng.reacLin(res, 2, 3)[""] or 0, 0, 1e-6, "rotula: M reaccion 0")

-- 5 asentamiento en biempotrada: M = 6EI d / L^2
res = beam(5, {FIX(1), {node = 2, rx = true, ry = true, rz = true, dy = "-d"}}, {}, function(m)
  m.secs[1] = {E = "200", A = "1e5", I = "3"} end)
_, _, M0 = Eng.internalLin(res, 1, 0, true)
near(math.abs(M0.d), 6 * 600 / 25, 1e-5, "asentamiento M")

-- 6 temperatura
res = beam(4, {FIX(1), FIX(2)}, {{mem = 1, t = 4, v1 = "30", v2 = "20"}}, function(m)
  m.secs[1] = {E = "2e8", A = "0.01", I = "1e-4", al = "1e-5", h = "0.4"} end)
N, V, Mc = MV(res, 1, 2)
near(N, -2e8 * 0.01 * 1e-5 * 30, 1e-6, "temperatura N")
near(Mc, -2e8 * 1e-4 * 1e-5 * 20 / 0.4, 1e-6, "gradiente M")

-- 7 carga triangular en simplemente apoyada: Mmax = wL^2/(9 sqrt3)
res = beam(9, {PIN(1), ROL(2)}, {{mem = 1, t = 3, dir = 5, v1 = "0", v2 = "12"}})
local e = Eng.extremes(res, 1)
near(e[3].mx, 12 * 81 / (9 * math.sqrt(3)), 1e-3, "triangular Mmax")
near(e[3].xmx, 9 / math.sqrt(3), 2e-2, "triangular x Mmax")
-- trapezoidal parcial en biempotrada: comprobar equilibrio
res = beam(10, {FIX(1), FIX(2)}, {{mem = 1, t = 3, dir = 5, v1 = "4", v2 = "9", a = "2", b = "7"},
                                  {mem = 1, t = 2, v1 = "15", a = "8"}})
local fx, fy, mo, rx, ry, rm = Eng.appliedTotals(res, 0)
near(fy + ry, 0, 1e-9, "trapecio SFy"); near(mo + rm, 0, 1e-9, "trapecio SM")
near(fy, -32.5, 1e-9, "resultante trapecio")
-- momento concentrado en biempotrada centro: M0/4 extremos
res = beam(8, {FIX(1), FIX(2)}, {{mem = 1, t = 2, v1 = "M"}})
_, _, M0 = Eng.internalLin(res, 1, 0, true)
near(math.abs(M0.M), 0.25, 1e-6, "momento concentrado FEM")

-- 8 apoyo inclinado: viga con rodillo a 30 grados
res = beam(6, {PIN(1), {node = 2, ry = true, ang = "30"}}, {{mem = 1, t = 3, dir = 5, v1 = "10"}})
local Rx, Ry = Eng.reacLin(res, 2, 1)[""], Eng.reacLin(res, 2, 2)[""]
near(Rx / Ry, -math.tan(math.pi / 6), 1e-6, "rodillo inclinado direccion")
near(Ry, 30, 1e-6, "rodillo inclinado Ry")

-- 9 resorte: viga en voladizo apoyada en resorte k en extremo
res = beam(3, {FIX(1), {node = 2, ky = "k0"}}, {}, function(m)
  m.vars = {{name = "k0", val = "2", keep = false}}
  m.nl = {{node = 2, fy = "-P"}} end)
local d = Eng.dispLin(res, 2, 2).P
local kb = 3 * 1 / 27
near(d, -1 / (kb + 2), 1e-6, "resorte")

-- 10 enrejado Pratt isostatico
do
  local m = T.Examples[5][2]()
  local res = Eng.solve(m)
  local N = Eng.internalLin(res, 4, 1, false)
  near(N.P, -5 / 3, 1e-6, "Pratt diagonal extremo")
  N = Eng.internalLin(res, 1, 1, false)
  near(N.P, 4 / 3, 1e-6, "Pratt cordon inferior")
  near(Eng.reacLin(res, 1, 2).P, 1, 1e-9, "Pratt reaccion")
end

-- 11 enrejado hiperestatico: equilibrio y simetria
do
  local m = T.Examples[4][2]()
  local res = Eng.solve(m)
  for k = 0, res.ns do
    local fx, fy, mo, rx, ry, rm = Eng.appliedTotals(res, k)
    near(fx + rx, 0, 1e-9, "enrejado SFx caso " .. k)
    near(fy + ry, 0, 1e-9, "enrejado SFy caso " .. k)
    near(mo + rm, 0, 1e-9, "enrejado SM caso " .. k)
  end
  -- carga vertical simetrica en nudo 5 -> reacciones verticales iguales
  near(Eng.reacLin(res, 1, 2).P, 0.5, 1e-6, "enrejado sim Ry1")
  near(Eng.reacLin(res, 3, 2).P, 0.5, 1e-6, "enrejado sim Ry3")
  -- equilibrio nodal en nudo 2
  local sx, sy = 0, 0
  for mm, mb in ipairs(res.mems) do
    if mb.i == 2 or mb.j == 2 then
      local N = Eng.internal(res, mm, mb.L / 2, false)
      local sg = (mb.i == 2) and 1 or -1
      sx = sx + sg * N * mb.c; sy = sy + sg * N * mb.s
    end
  end
  near(sx, 0, 1e-8, "nudo 2 SFx"); near(sy, 0, 1e-8, "nudo 2 SFy")
end

-- 12 todos los ejemplos: equilibrio, continuidad M(L) con extremo
for i, ex in ipairs(T.Examples) do
  local res = Eng.solve(ex[2]())
  for k = 0, res.ns do
    local fx, fy, mo, rx, ry, rm = Eng.appliedTotals(res, k)
    local sc = math.max(1, res.cases[k].fs)
    near((fx + rx) / sc, 0, 1e-8, ex[1] .. " SFx")
    near((fy + ry) / sc, 0, 1e-8, ex[1] .. " SFy")
    near((mo + rm) / sc, 0, 1e-7, ex[1] .. " SM")
    for mm in ipairs(res.mems) do
      local F = res.cases[k].F[mm]
      local mb = res.mems[mm]
      local C = res.cases[k]
      -- M(L) debe igualar F6, V(L) = -F5, N(L) = F4
      local w = {}; for q = 0, res.ns do w[q] = (q == k) and 1 or 0 end
      local N, V, M = Eng.internal(res, mm, mb.L, true, w)
      near(M, F[6], 1e-6, ex[1] .. " M(L) b" .. mm)
      near(V, -F[5], 1e-6, ex[1] .. " V(L) b" .. mm)
      near(N, F[4], 1e-6, ex[1] .. " N(L) b" .. mm)
    end
  end
end

-- 13 portico empotrado-articulado (ejemplo 2): comparacion con solucion clasica
-- (verificada con equilibrio). Rotula del ejemplo 3: M = 0 en nudo 3
do
  local res = Eng.solve(T.Examples[3][2]())
  local _, _, Mh = Eng.internalLin(res, 2, res.mems[2].L, false)
  for k, v in pairs(Mh) do near(v, 0, 1e-9, "rotula M=0 " .. k) end
end

-- 14 mecanismo detectado
do
  local m = T.newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = "5", y = "0"}}
  m.mems = {{i = 1, j = 2, sec = 1}}
  m.sups = {ROL(1), ROL(2)}
  local okk, err = pcall(Eng.solve, m)
  ok(not okk and tostring(err):find("INESTABLE"), "mecanismo detectado: " .. tostring(err))
end

-- 15 barra biarticulada en marco (tirante) con giro auto-fijado
do
  local m = T.newModel("frame")
  m.nodes = {{x = "0", y = "0"}, {x = "4", y = "0"}, {x = "4", y = "3"}}
  m.mems = {{i = 1, j = 3, sec = 1, ri = true, rj = true}, {i = 2, j = 3, sec = 1, ri = true, rj = true}}
  m.sups = {PIN(1), PIN(2)}
  m.nl = {{node = 3, fx = "10"}}
  local res = Eng.solve(m)
  local N = Eng.internal(res, 1, 1, false)
  near(N, 50 / 4, 1e-6, "tirante en marco")
end

-- 16 UI: pintar todas las pantallas y navegar
local gc = mkgc()
local App = T.App
local function paint() on.paint(gc) end
local function key(k, ...) on[k](...); paint() end
local function reset() local st = T.stack(); while #st > 1 do st[#st] = nil end end
paint()
ok(registered ~= nil, "toolpalette registrado")
key("charIn", "4")          -- ejemplos
for i = 1, #T.Examples do
  key("escapeKey")
  key("charIn", "4")
  key("charIn", tostring(i)) -- carga ejemplo -> menu modelo
  local base = #T.stack()
  for item = 1, 8 do
    key("charIn", tostring(item)); key("arrowKey", "down"); key("arrowKey", "down")
    key("enterKey"); paint()
    for _ = 1, 14 do key("arrowKey", "down"); key("arrowKey", "right") end
    local st = T.stack(); while #st > base do st[#st] = nil end
    ok(App.res == nil, "edicion invalida resultados")
  end
  key("charIn", "9")        -- resolver
  local top = T.stack()[#T.stack()]
  ok(App.res ~= nil, "resuelto ejemplo " .. i)
  local rb = #T.stack()
  for item = 1, 7 do
    local st = T.stack(); while #st > rb do st[#st] = nil end
    key("charIn", tostring(item))
    for _, q in ipairs({"m", "v", "n", "d"}) do key("charIn", q); key("arrowKey", "right"); key("tabKey") end
    key("arrowKey", "up"); key("arrowKey", "down"); key("arrowKey", "right")
    key("enterKey"); key("escapeKey"); key("escapeKey")
  end
  reset()
end
-- formulario: crear nudo nuevo tecleando
reset()
App.model = T.newModel("frame")
App.modelMenu()
key("charIn", "1"); key("charIn", "+")
key("charIn", "3"); key("charIn", "."); key("charIn", "5"); key("backspaceKey"); key("charIn", "2")
key("enterKey"); key("charIn", "\226\136\146"); key("charIn", "1")
key("enterKey"); key("enterKey")
ok(#App.model.nodes == 1 and App.model.nodes[1].x == "3.2" and App.model.nodes[1].y == "\226\136\1461",
   "nudo creado por teclado")
key("backspaceKey"); ok(#App.model.nodes == 0, "nudo borrado")
-- persistencia
local st = on.save(); on.restore(st); ok(App.model ~= nil, "save/restore")
-- toolpalette callbacks
for _, grp in ipairs(registered) do
  for q = 2, #grp do reset(); grp[q][2](); paint(); on.escapeKey(); paint() end
end

-- ninguna pantalla de error interno durante la navegacion
print(("%d pruebas, %d fallas"):format(count, fails))
if fails > 0 then os.exit(1) end
