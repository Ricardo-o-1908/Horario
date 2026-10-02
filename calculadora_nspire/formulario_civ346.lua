-- =====================================================================
--  FORMULARIO CIV-346  --  Solo lectura (TI-Nspire CX II CAS)
--  Motor de tipografia matematica propio: fracciones, raices,
--  sub/superindices, letras griegas, cursivas y parentesis escalados.
--
--  Teclas:  flechas arriba/abajo = desplazar
--           flechas izq/der      = pagina anterior / siguiente
--           tab / shift+tab      = capitulo siguiente / anterior
--           + / -                = zoom          enter/esc = indice
-- =====================================================================
platform.apiLevel = "2.2"

-- ------------------------------------------------------------- simbolos
local SYM = {
  cdots = "\194\183\194\183\194\183", ldots = "...",
  alpha = "α", beta = "β", gamma = "γ", delta = "δ", Delta = "Δ", epsilon = "ε", eta = "η",
  theta = "θ", lambda = "λ", Lambda = "Λ", mu = "μ", nu = "ν", pi = "π", rho = "ρ", sigma = "σ",
  Sigma = "Σ", tau = "τ", phi = "φ", Phi = "Φ", psi = "ψ", Psi = "Ψ", omega = "ω", Omega = "Ω",
  chi = "χ", xi = "ξ", nabla = "∇", partial = "∂", cdot = "·", times = "×", le = "≤", ge = "≥",
  approx = "≈", sim = "~", to = "→", Rightarrow = "⇒", infty = "∞", pm = "±", neq = "≠",
  sum = "Σ", int = "∫", deg = "°", lt = "<", gt = ">", vee = "∨", wedge = "∧", prime = "'",
}
local FUN = {ln = 1, log = 1, exp = 1, sin = 1, cos = 1, tan = 1, sen = 1, cosh = 1, coth = 1, senh = 1,
  asinh = 1, max = 1, min = 1, tg = 1, arcsenh = 1, acos = 1}
local GAP = {[","] = 0.17, [";"] = 0.28, quad = 1.0, qquad = 2.0}

-- ------------------------------------------------------------- parser
local function parse(s)
  local i, n = 1, #s
  local seq
  local function arg()
    while s:sub(i, i) == " " do i = i + 1 end
    local c = s:sub(i, i)
    if c == "{" then i = i + 1; return seq("}") end
    if c == "\\" then
      local j = i + 1
      local name = s:match("^%a+", j)
      if name then i = j + #name; return {{k = "t", s = SYM[name] or name, it = false}} end
    end
    -- un caracter UTF-8
    local ch = s:match("^[%z\1-\127\194-\244][\128-\191]*", i)
    i = i + #ch
    return {{k = "t", s = ch, it = ch:match("^%a$") ~= nil}}
  end
  function seq(stop)
    local out = {}
    while i <= n do
      local c = s:sub(i, i)
      if stop and c == stop then i = i + 1; return out end
      if c == "{" then i = i + 1; out[#out + 1] = {k = "g", a = seq("}")}
      elseif c == " " then i = i + 1; out[#out + 1] = {k = "sp", w = 0.25}
      elseif c == "^" or c == "_" then
        i = i + 1
        local base = table.remove(out) or {k = "t", s = "", it = false}
        if base.k ~= "scr" then base = {k = "scr", b = base} end
        if c == "^" then base.sup = arg() else base.sub = arg() end
        out[#out + 1] = base
      elseif c == "\\" then
        local name = s:match("^%a+", i + 1)
        if not name then
          local sym = s:sub(i + 1, i + 1); i = i + 2
          if GAP[sym] then out[#out + 1] = {k = "sp", w = GAP[sym]}
          else out[#out + 1] = {k = "t", s = sym, it = false} end
        else
          i = i + 1 + #name
          if name == "frac" then local a = arg(); local b = arg(); out[#out + 1] = {k = "frac", a = a, b = b}
          elseif name == "sqrt" then
            local idx
            if s:sub(i, i) == "[" then local j = s:find("]", i, true); idx = s:sub(i + 1, j - 1); i = j + 1 end
            out[#out + 1] = {k = "sqrt", a = arg(), idx = idx}
          elseif name == "p" then out[#out + 1] = {k = "par", a = arg(), l = "(", r = ")"}
          elseif name == "b" then out[#out + 1] = {k = "par", a = arg(), l = "[", r = "]"}
          elseif name == "abs" then out[#out + 1] = {k = "par", a = arg(), l = "|", r = "|"}
          elseif name == "text" or name == "mathrm" then
            local a = arg()
            for _, nd in ipairs(a) do if nd.k == "t" then nd.it = false elseif nd.k == "sp" then nd.w = 0.3 end end
            out[#out + 1] = {k = "g", a = a}
          elseif name == "bar" then out[#out + 1] = {k = "bar", a = arg()}
          elseif name == "vec" then out[#out + 1] = {k = "g", a = arg()}
          elseif GAP[name] then out[#out + 1] = {k = "sp", w = GAP[name]}
          elseif FUN[name] then out[#out + 1] = {k = "t", s = name, it = false, fn = true}
          else out[#out + 1] = {k = "t", s = SYM[name] or name, it = false, op = (name == "cdot" or name == "times" or name == "le" or name == "ge" or name == "approx" or name == "Rightarrow" or name == "to" or name == "pm")} end
        end
      else
        local ch = s:match("^[%z\1-\127\194-\244][\128-\191]*", i)
        i = i + #ch
        local isop = ch:match("^[=+%-<>]$") ~= nil
        out[#out + 1] = {k = "t", s = (ch == "-") and "\226\136\146" or ch, it = ch:match("^%a$") ~= nil, op = isop}
      end
    end
    return out
  end
  return seq(nil)
end

-- ------------------------------------------------------------- layout
local ZOOM = 2
local SIZES = {{9, 7, 7}, {11, 9, 7}, {12, 9, 7}, {16, 11, 9}}
local gcm   -- contexto grafico activo
local function fsz(L) local t = SIZES[ZOOM]; return t[math.min(L + 1, 3)] end
local function setf(L, it, bold) gcm:setFont("serif", bold and "b" or (it and "i" or "r"), fsz(L)) end

local measure
local function mseq(seq, L)
  local w, a, d = 0, 0, 0
  for _, nd in ipairs(seq) do
    local nw, na, nde = measure(nd, L)
    nd.mw, nd.ma, nd.md = nw, na, nde
    w = w + nw; a = math.max(a, na); d = math.max(d, nde)
  end
  return w, a, d
end
function measure(nd, L)
  local sz = fsz(L)
  if nd.k == "t" then
    setf(L, nd.it)
    local w = gcm:getStringWidth(nd.s)
    if nd.op then w = w + 0.45*sz end
    if nd.fn then w = w + 0.15*sz end
    if nd.it then w = w + 0.06*sz end
    return w, 0.74*sz, 0.24*sz
  elseif nd.k == "sp" then return nd.w*sz, 0, 0
  elseif nd.k == "g" then return mseq(nd.a, L)
  elseif nd.k == "frac" then
    local Ln = (L >= 1) and L + 1 or L
    if nd.nested then Ln = L + 1 end
    for _, c in ipairs(nd.a) do if c.k == "frac" then c.nested = true end end
    for _, c in ipairs(nd.b) do if c.k == "frac" then c.nested = true end end
    local wa, aa, da = mseq(nd.a, Ln); local wb, ab, db = mseq(nd.b, Ln)
    nd.wa, nd.aa, nd.da, nd.wb, nd.ab, nd.db, nd.Ln = wa, aa, da, wb, ab, db, Ln
    local axis = 0.28*sz
    local w = math.max(wa, wb) + 0.4*sz
    return w, axis + 2 + da + aa, -axis + 3 + ab + db
  elseif nd.k == "sqrt" then
    local wa, aa, da = mseq(nd.a, L)
    nd.wa, nd.aa, nd.da = wa, aa, da
    return wa + 0.75*sz, aa + 3, da + 1
  elseif nd.k == "par" then
    local wa, aa, da = mseq(nd.a, L)
    nd.wa, nd.aa, nd.da = wa, aa, da
    local big = (aa + da) > 1.25*sz
    nd.big = big
    local pw = big and 0.42*sz or 0.36*sz
    nd.pw = pw
    return wa + 2*pw, aa + (big and 1 or 0), da + (big and 1 or 0)
  elseif nd.k == "bar" then
    local wa, aa, da = mseq(nd.a, L)
    return wa, aa + 3, da
  elseif nd.k == "scr" then
    local bw, ba, bd = measure(nd.b, L); nd.b.mw, nd.b.ma, nd.b.md = bw, ba, bd
    local sw, sa, sd, uw, ua, ud = 0, 0, 0, 0, 0, 0
    if nd.sup then sw, sa, sd = mseq(nd.sup, L + 1) end
    if nd.sub then uw, ua, ud = mseq(nd.sub, L + 1) end
    nd.sw, nd.uw = sw, uw
    nd.up = math.max(0.42*sz, ba - 0.35*sz)
    nd.dn = 0.22*sz + (nd.sup and 0.08*sz or 0)
    return bw + math.max(sw, uw) + 0.5, math.max(ba, nd.up + sa), math.max(bd, nd.dn + ud)
  end
  return 0, 0, 0
end

local drawn
local function dseq(seq, x, y, L)
  for _, nd in ipairs(seq) do drawn(nd, x, y, L); x = x + nd.mw end
end
function drawn(nd, x, y, L)
  local sz = fsz(L)
  if nd.k == "t" then
    setf(L, nd.it)
    local ox = nd.op and 0.22*sz or 0
    gcm:drawString(nd.s, x + ox, y, "baseline")
  elseif nd.k == "g" then dseq(nd.a, x, y, L)
  elseif nd.k == "frac" then
    local axis = 0.28*sz
    local yb = y - axis
    gcm:fillRect(x + 1, yb, nd.mw - 2, 1)
    dseq(nd.a, x + (nd.mw - nd.wa)/2, yb - 2 - nd.da, nd.Ln)
    dseq(nd.b, x + (nd.mw - nd.wb)/2, yb + 3 + nd.ab, nd.Ln)
  elseif nd.k == "sqrt" then
    local top = y - nd.aa - 2
    local bot = y + nd.da
    local x0 = x + 1
    gcm:drawLine(x0, y - 0.3*sz, x0 + 0.2*sz, y - 0.4*sz)
    gcm:drawLine(x0 + 0.2*sz, y - 0.4*sz, x0 + 0.38*sz, bot)
    gcm:drawLine(x0 + 0.38*sz, bot, x0 + 0.62*sz, top)
    gcm:drawLine(x0 + 0.62*sz, top, x + nd.mw, top)
    if nd.idx then setf(L + 1, false); gcm:drawString(nd.idx, x, y - 0.45*sz, "baseline") end
    dseq(nd.a, x + 0.7*sz, y, L)
  elseif nd.k == "par" then
    local top, bot = y - nd.aa, y + nd.da
    local h = bot - top
    if not nd.big then
      setf(L, false)
      gcm:drawString(nd.l, x, y, "baseline")
      gcm:drawString(nd.r, x + nd.pw + nd.wa, y, "baseline")
    else
      local xl, xr = x + 1, x + nd.pw + nd.wa + nd.pw - 2
      if nd.l == "(" then
        gcm:drawArc(xl, top, nd.pw*1.6, h, 100, 160)
        gcm:drawArc(xr - nd.pw*1.6 + 2, top, nd.pw*1.6, h, -80, 160)
      elseif nd.l == "[" then
        gcm:drawLine(xl + 2, top, xl + 2, bot); gcm:drawLine(xl + 2, top, xl + nd.pw - 1, top); gcm:drawLine(xl + 2, bot, xl + nd.pw - 1, bot)
        gcm:drawLine(xr, top, xr, bot); gcm:drawLine(xr, top, xr - nd.pw + 3, top); gcm:drawLine(xr, bot, xr - nd.pw + 3, bot)
      else
        gcm:drawLine(xl + 2, top, xl + 2, bot); gcm:drawLine(xr, top, xr, bot)
      end
    end
    dseq(nd.a, x + nd.pw, y, L)
  elseif nd.k == "bar" then
    gcm:fillRect(x, y - nd.ma + 1, nd.mw, 1)
    dseq(nd.a, x, y, L)
  elseif nd.k == "scr" then
    drawn(nd.b, x, y, L)
    local xs = x + nd.b.mw
    local ic = (nd.b.k == "t" and nd.b.it) and 0.08*sz or 0
    if nd.sup then dseq(nd.sup, xs, y - nd.up, L + 1) end
    if nd.sub then dseq(nd.sub, xs - ic, y + nd.dn, L + 1) end
  end
end

local function mathBox(src)
  local tree = parse(src)
  local w, a, d = mseq(tree, 0)
  return {tree = tree, w = w, a = a, d = d}
end
local function drawMath(box, x, y) dseq(box.tree, x, y, 0) end

MATH = {parse = parse, box = mathBox, draw = drawMath, setgc = function(g) gcm = g end, zoom = function(z) if z then ZOOM = z end; return ZOOM end, NZ = #SIZES}

-- =====================================================================
--  CONTENIDO DEL FORMULARIO
-- =====================================================================
DOC = {}
local function CAP(num, title, blocks) DOC[#DOC + 1] = {num = num, title = title, blocks = blocks} end
local function h(t) return {"h", t} end
local function p(t) return {"p", t} end
local function b(t) return {"b", t} end
local function f(t, n) return {"f", t, n} end
local function n(t) return {"n", t} end
local function tab(hd, rows) return {"t", hd, rows} end

-- ------------------------------------------------------------ 01
CAP("C01", "Suelos, Darcy y K", {
  h("Porosidad"),
  f([[n = \frac{V_{vac}}{V_{T}} \qquad n_{e} = \frac{V_{drenado}}{V_{T}}]], "n_e: fracción que participa en el flujo"),
  tab({"Suelo", "n [%]", "ne [%]"}, {{"Arcillas", "50", "4"}, {"Limos", "48", "8"}, {"Arenas", "40", "33"}, {"Gravas", "25", "22"}, {"Rocas", "1-8", "0.5-4"}}),
  h("Formaciones"),
  tab({"Tipo", "Almacena", "Transmite", "Ejemplo"}, {{"Acuífero", "alta", "alta", "gravas, arenas"}, {"Acuitardo", "alta", "baja", "limos"},
    {"Acuicludo", "alta", "nula", "arcillas"}, {"Acuífugo", "nula", "nula", "rocas"}}),
  h("Carga hidráulica"),
  f([[h = z + \frac{p}{\gamma} + \frac{v^{2}}{2g} \approx z + \frac{p}{\gamma}]], "v muy pequeña en flujo subterráneo"),
  h("Ley de Darcy (Re < 4)"),
  f([[Q = -K A \frac{dh}{dl} = K A i]]),
  f([[v = \frac{Q}{A} = K i \qquad v_{r} = \frac{v}{n_{e}}]], "v: velocidad de Darcy; vr: velocidad real"),
  f([[v = -K \nabla h \qquad K_{x} = K_{y} = K_{z} \Rightarrow \text{isotrópico}]]),
  f([[K = k \frac{\rho g}{\mu} \qquad 1 \text{darcy} = 9.87 \times 10^{-9} cm^{2}]]),
  h("Métodos indirectos"),
  f([[K = C d_{10}^{2} \qquad [cm/s], d_{10} [cm] ]], "Hazen: arenas con d10 = 0.1–0.3 mm"),
  tab({"Arena", "C"}, {{"muy fina / con finos", "40-80"}, {"media bien distribuida", "80-120"}, {"gruesa mal distribuida", "80-120"}, {"gruesa limpia", "120-150"}}),
  f([[K = \frac{\gamma}{\mu} \frac{10.0219 d_{10}^{2}}{K_{1}}]], "Slichter (1/K1 de tabla según n)"),
  tab({"n", "1/K1", "n", "1/K1"}, {{"0.26", "0.00187", "0.38", "0.04154"}, {"0.28", "0.01517", "0.40", "0.04922"}, {"0.30", "0.01905", "0.42", "0.05789"},
    {"0.32", "0.02356", "0.44", "0.06776"}, {"0.34", "0.02878", "0.46", "0.07838"}, {"0.36", "0.03473", "", ""}}),
  f([[K = \frac{\rho g}{\mu} \frac{n^{3}}{(1-n)^{2}} \frac{d_{50}^{2}}{180}]], "Kozeny–Carman"),
  h("Métodos directos"),
  f([[K = \frac{4 L Q}{\pi (h_{1} - h_{2}) D^{2}}]], "Permeámetro de carga constante (h1−h2 < 0.5 L)"),
  f([[K = \frac{d_{t}^{2} L}{d_{c}^{2} t} \ln \p{\frac{h_{0}}{h}}]], "Permeámetro de carga variable (cohesivos)"),
  f([[f = \frac{R}{2(t_{2} - t_{1})} \ln \p{\frac{2h_{1} + R}{2h_{2} + R}}]], "Método Porchet"),
  f([[K \approx \frac{L^{2}}{\bar{t} \Delta h}]], "Dos pozos de observación"),
  f([[\ln \p{\frac{C_{0}}{C}} = \frac{4 v}{\pi D} t , \quad v = K i]], "Trazador en un pozo"),
  h("Estratificación"),
  f([[K_{eq} = \frac{\Sigma K_{i} m_{i}}{\Sigma m_{i}}]], "Flujo paralelo a los estratos"),
  f([[K_{eq} = \frac{\Sigma m_{i}}{\Sigma \frac{m_{i}}{K_{i}}}]], "Flujo perpendicular a los estratos"),
  h("Propiedades del acuífero"),
  f([[T = K m \qquad V_{dren} = S A \Delta h]]),
})

-- ------------------------------------------------------------ 02
CAP("C02", "Pozos en régimen permanente", {
  h("Almacenamiento"),
  f([[S_{s} = g \rho (\alpha + \eta \beta) \qquad S = S_{s} m]], "Confinado: 1E−4 < S < 0.005"),
  f([[S = S_{y} = n_{e}]], "Libre: 0.05 < S < 0.2"),
  h("Ecuación general (confinado)"),
  f([[\frac{\partial}{\partial x}\p{K_{x}\frac{\partial h}{\partial x}} + \frac{\partial}{\partial y}\p{K_{y}\frac{\partial h}{\partial y}} + \frac{\partial}{\partial z}\p{K_{z}\frac{\partial h}{\partial z}} = S_{s}\frac{\partial h}{\partial t}]]),
  f([[\nabla^{2} h = 0]], "Permanente, homogéneo, isotrópico (Laplace)"),
  f([[\frac{1}{r}\frac{\partial h}{\partial r} + \frac{\partial^{2} h}{\partial r^{2}} + \frac{1}{r^{2}}\frac{\partial^{2} h}{\partial \theta^{2}} + \frac{\partial^{2} h}{\partial z^{2}} = 0]], "Laplace en cilíndricas"),
  h("Supuestos"),
  b("Acuífero homogéneo, isótropo, infinito; nivel inicial horizontal."),
  b("Flujo radial; radio de influencia R; sin pérdidas en el pozo."),
  h("Pozo de penetración total"),
  f([[\Delta h = H - h = \frac{Q}{2 \pi T} \ln \p{\frac{R}{r}}]], "Thiem – acuífero confinado"),
  f([[H^{2} - h^{2} = \frac{Q}{\pi K} \ln \p{\frac{R}{r}}]], "Dupuit – acuífero libre (Dupuit–Forchheimer)"),
  f([[\Delta h = \frac{Q}{2 \pi K H} \ln \frac{R}{r}]], "Libre con Δh < 0.1 H"),
  f([[\Delta h' = \Delta h - \frac{\Delta h^{2}}{2H} = \frac{Q}{2 \pi K H} \ln \frac{R}{r}]], "Corrección de Jacob (Δh < 0.5 H)"),
  h("Acuífero libre de gran potencia"),
  f([[\Delta h = \frac{Q}{4 \pi K} \p{\frac{1}{r} - \frac{1}{R}}]], "Escurrimiento esférico"),
  f([[\Delta h = \frac{Q}{2 \pi K} \p{\frac{1}{r} - \frac{1}{R}}]], "Pozo somero: semi-esférico"),
  h("Penetración parcial (confinado)"),
  f([[h_{1} - h_{2} = \frac{Q}{2 \pi K} \p{\frac{1}{b} \ln \frac{1.6 b}{r_{2}} - \frac{1}{b}\text{asinh}\frac{b}{r_{1}}}]], "Girinsky (b/m < 0.3)"),
  f([[h_{1} - h_{2} = \frac{Q}{2 \pi K m}\b{\ln\frac{r_{1}}{r_{2}} + \text{asinh}\frac{m}{r_{1}} - \text{asinh}\frac{m}{r_{2}} - \frac{m}{b}\p{\text{asinh}\frac{b}{r_{1}} - \text{asinh}\frac{b}{r_{2}}}}]], "Nasberg (b/m ≥ 0.3)"),
  h("Radio de influencia"),
  f([[R = 3000 \Delta h_{0} \sqrt{K}]], "Sichardt: R [m], Δh0 [m], K [m/s]"),
  tab({"Material", "D [mm]", "R [m]"}, {{"Gravas", ">10", "1500"}, {"Gravillas", "2-10", "500-1500"}, {"Gravilla-arcilla", "1-2", "400-500"},
    {"Arena gruesa", "0.5-1", "200-400"}, {"Arena media", "0.25-0.5", "100-200"}, {"Arena fina", "0.1-0.25", "50-100"}, {"Arena muy fina", "0.05-0.1", "10-50"}}),
})

-- ------------------------------------------------------------ 03
CAP("C03", "Interferencia e imágenes", {
  h("Superposición"),
  f([[\Delta h_{T} = \frac{1}{2 \pi T} \Sigma Q_{i} \ln \p{\frac{R_{i}}{r_{i}}}]], "Sólo pozos con ri < Ri. ¡Libre: linealizar!"),
  f([[\Delta h_{T} = \frac{Q}{2 \pi T} \ln \p{\frac{R^{n}}{r_{1} r_{2} \cdots r_{n}}}]], "n pozos con igual Q y R"),
  h("Recarga lineal (imagen inyecta −Q)"),
  f([[\Delta h = \frac{Q}{2 \pi T} \ln \frac{r_{i}}{r_{r}}]], "rr < R y ri < R"),
  f([[\Delta h = \frac{Q}{2 \pi T} \ln \frac{R}{r_{r}}]], "rr < R ≤ ri ;  Δh = 0 si ambos ≥ R"),
  h("Barrera impermeable (imagen extrae +Q)"),
  f([[\Delta h = \frac{Q}{2 \pi T} \ln \frac{R^{2}}{r_{i} r_{r}}]], "rr < R y ri < R"),
  f([[\Delta h = \frac{Q}{2 \pi T} \ln \frac{R}{r_{r}}]], "rr < R ≤ ri"),
  n("Esquina de recarga: 3 imágenes (−Q, −Q, +Q). Recarga + barrera: verificar 2b + a < R y 2a + b < R (imágenes de imágenes)."),
  h("Pozo en varias napas"),
  f([[Q = \frac{2 \pi}{\ln \p{\frac{R}{r}}} \Sigma \Delta h_{i} K_{i} m_{i}]]),
  f([[\Delta h_{T} = \Delta h_{1} \p{\frac{T_{1} + T_{2}}{T_{2}}}]], "Condición inicial: una napa alimenta a la otra"),
  f([[Q = Q_{0} + \frac{2 \pi K_{2} m_{2} \Delta h_{2}}{\ln \p{\frac{R}{r_{0}}}}]], "N.E. en el estrato impermeable de la napa 2"),
})

-- ------------------------------------------------------------ 04
CAP("C04", "Redes de flujo y analogías", {
  h("Funciones potencial y de flujo"),
  f([[\phi = K h \qquad \nabla^{2} \phi = 0 \qquad \nabla^{2} \psi = 0]]),
  f([[v_{x} = \frac{\partial \phi}{\partial x} = -\frac{\partial \psi}{\partial z} \qquad v_{z} = \frac{\partial \phi}{\partial z} = \frac{\partial \psi}{\partial x}]], "Cauchy–Riemann"),
  f([[m_{\psi} \cdot m_{\phi} = -1]], "Familias ortogonales"),
  h("Caudal"),
  f([[\Delta q = K \frac{h}{N_{d}} \frac{a}{l} \qquad q = K h \frac{N_{t}}{N_{d}} \frac{a}{l}]], "Nt: canales; Nd: caídas; a/l = 1 cuadrados"),
  h("Condiciones de borde"),
  b("Barrera impermeable y eje de simetría: líneas de flujo ψ."),
  b("Superficie de agua en reposo: equipotencial φ."),
  b("Superficie libre en movimiento: línea de flujo con p/γ = 0 (Δφ = KΔz)."),
  h("Medios heterogéneo y anisotrópico"),
  f([[\frac{K_{1}}{K_{2}} = \frac{\tan \alpha}{\tan \beta}]], "Refracción"),
  f([[x' = a x , \quad a = \sqrt{\frac{K_{z}}{K_{x}}} , \quad K_{eq} = \sqrt{K_{x} K_{z}}]], "Transformación del dominio"),
  h("Analogías"),
  f([[\bar{v} = \frac{\gamma b^{2}}{12 \mu} \frac{dh}{dx}]], "Hele–Shaw (placas paralelas)"),
  f([[i = \frac{1}{R} V \Leftrightarrow v = K \frac{\Delta h}{L}]], "Eléctrica (Ohm)"),
})

-- ------------------------------------------------------------ 05
CAP("C05", "Costero y régimen impermanente", {
  h("Intrusión salina (Ghyben–Herzberg)"),
  f([[h_{s} = h_{d} \p{\frac{\gamma_{d}}{\gamma_{s} - \gamma_{d}}} \approx 40 h_{d}]]),
  f([[h_{d}^{2} = \frac{N (R^{2} - r^{2})}{2 K (1 + \delta)} , \quad \delta = \frac{\gamma_{d}}{\gamma_{s} - \gamma_{d}}]], "Isla circular con recarga N"),
  h("Theis (1935) – confinado"),
  f([[\frac{\partial^{2} h}{\partial r^{2}} + \frac{1}{r}\frac{\partial h}{\partial r} = \frac{S}{T}\frac{\partial h}{\partial t}]]),
  f([[\Delta h = \frac{Q}{4 \pi T} W(u) , \qquad u = \frac{r^{2} S}{4 T t}]]),
  f([[W(u) = -0.5772 - \ln u + u - \frac{u^{2}}{2 \cdot 2!} + \frac{u^{3}}{3 \cdot 3!} - \cdots]]),
  tab({"u", "W(u)", "u", "W(u)"}, {{"1E-4", "8.63", "0.05", "2.47"}, {"5E-4", "7.02", "0.1", "1.82"}, {"1E-3", "6.33", "0.2", "1.22"},
    {"5E-3", "4.73", "0.5", "0.56"}, {"0.01", "4.04", "1", "0.219"}, {"0.02", "3.35", "2", "0.049"}}),
  h("Cooper–Jacob (u < 0.01)"),
  f([[\Delta h = \frac{Q}{4 \pi T} \ln \p{\frac{2.24 T t}{S r^{2}}}]], "Error: u = 0.2 → 20%; 0.14 → 10%; 0.09 → 5%"),
  f([[R(t) = \sqrt{\frac{2.24 T t}{S}} \qquad \Delta h = \frac{Q}{2 \pi T} \ln \frac{R(t)}{r}]]),
  f([[H^{2} - h^{2} = \frac{Q}{\pi K} \ln \p{\frac{R(t)}{r}}]], "Libre transiente"),
  h("Interferencia transiente"),
  f([[t_{i} = \frac{r_{i}^{2} S}{2.24 T}]], "Tiempo de influencia: el pozo afecta si t > ti"),
  f([[\Delta h_{A} = \Sigma \frac{Q_{i}}{4 \pi T} \ln \p{\frac{2.24 T t}{r_{i}^{2} S}}]]),
  n("Imágenes: si tr < t < ti sólo actúa el real. Si t > ti: recarga Δh = Q/(2πT) ln(ri/rr); barrera Δh = Q/(2πT) ln(R(t)²/(rr ri))."),
  h("Bombeo variable"),
  f([[\Delta h = \frac{1}{4 \pi T} \Sigma \pm Q_{i} \ln \p{\frac{2.24 T}{r^{2} S}(t - t_{i})}]], "t ≥ ti ; +Q enciende, −Q apaga"),
})

-- ------------------------------------------------------------ 06
CAP("C06", "Pruebas de bombeo", {
  h("Tipos de prueba"),
  b("Gasto variable (escalonada): define Qmax y la cota de la bomba."),
  b("Gasto constante: 24 h, 180 min estabilizado, con 80% de Qmax."),
  b("Recuperación: desde tf se detiene el bombeo."),
  h("Curva tipo de Theis"),
  f([[T = \frac{Q}{4 \pi \Delta h^{*}} W(u) \qquad S = \frac{4 T u t^{*}}{r^{2}}]], "Superponer log-log sin rotar; leer W, 1/u, Δh*, t*"),
  f([[W(u) = \frac{u^{4} + 8.573u^{3} + 18.059u^{2} + 8.635u + 0.268}{u e^{u}(u^{4} + 9.573u^{3} + 25.633u^{2} + 21.097u + 3.958)}]], "Aproximación para u ≥ 1"),
  h("Cooper–Jacob (semilog)"),
  f([[T = \frac{Q}{4 \pi a} = \frac{2.303 Q}{4 \pi \Delta s_{ciclo}} \qquad S = \frac{2.24 T t_{0}}{r^{2}}]], "a: pendiente vs ln t; t0: corte con Δh = 0"),
  f([[\Delta h' = \Delta h - \frac{\Delta h^{2}}{2H}]], "Acuífero libre: corregir descensos"),
  h("Otros métodos"),
  f([[\Delta h = \frac{Q}{4 \pi T} W(u, \alpha) , \quad \alpha = S \p{\frac{r_{0}}{r}}^{2}]], "Papadopulos–Cooper (pozos de gran diámetro)"),
  f([[\Delta h = \frac{Q}{4 \pi T} W\p{u, \frac{r}{B}} , \quad B = \sqrt{\frac{T m'}{K'}}]], "Hantush–Jacob (semiconfinado)"),
  h("Recuperación"),
  f([[s' = \frac{Q}{4 \pi T} \ln \p{\frac{t}{t - t_{f}}} \Rightarrow T = \frac{2.303 Q}{4 \pi \Delta s'}]]),
  h("Bordes"),
  b("Barrera impermeable: la pendiente se duplica."),
  b("Recarga: la curva se estabiliza (pendiente nula)."),
  f([[r_{2} = r_{1} \sqrt{\frac{t_{2}}{t_{1}}} \qquad d = \frac{r_{1} + r_{2}}{2}]], "Igual desviación que el descenso original"),
  h("Pérdidas en el pozo"),
  f([[s_{w} = B Q + C Q^{2} \qquad \frac{s_{w}}{Q} = B + C Q]]),
})

-- ------------------------------------------------------------ 07
CAP("C07", "Drenes – régimen permanente", {
  h("Dren interceptor"),
  f([[H^{2} - h_{0}^{2} = \frac{2 q}{K} x_{1} \qquad Q_{z} = (q - q') L]]),
  h("Drenes paralelos abiertos"),
  f([[H^{2} - h_{0}^{2} = \frac{f D^{2}}{4 K}]], "Donnan: drenes abarcan toda la napa"),
  f([[(H - h_{0})(H + h_{0} + 2d) = \frac{f D^{2}}{4 K}]], "Hooghoudt: d < 0.6 m (si no, usar d')"),
  h("Drenes cerrados"),
  f([[(H + d)^{2} - d^{2} = \frac{f D^{2}}{4 K} \p{1 + \alpha \frac{4 d}{D}}]], "Dagan (1964)"),
  f([[\alpha = -\frac{1}{\pi} \ln \b{2 \p{\cosh \p{\frac{\pi r_{0}}{d}} - 1}}]], "r0/d ≥ 0.3 ⇒ α → 0"),
  tab({"r0/d", "α", "r0/d", "α"}, {{"0.01", "2.16", "0.08", "0.86"}, {"0.02", "1.75", "0.10", "0.75"}, {"0.04", "1.31", "0.20", "0.29"}, {"0.06", "1.05", "0.30", "0.01"}}),
  f([[H = \frac{D f}{K \pi} \ln \p{\frac{D}{5.44 r_{0}}}]], "Superficiales, espesor infinito"),
  f([[H = \frac{D f}{K} F \p{\frac{2 r_{0}}{D}, \frac{d}{D}}]], "Kirkham (1958)"),
  tab({"d/D", "0.0025", "0.005", "0.01", "0.02", "0.04"}, {{"0.01", "12.79", "12.57", "12.33", "12.03", "11.52"}, {"0.02", "6.761", "6.541", "6.318", "6.077", "5.771"},
    {"0.04", "3.864", "3.643", "3.421", "3.195", "2.954"}, {"0.08", "2.522", "2.301", "2.080", "1.858", "1.633"}, {"0.16", "1.961", "1.741", "1.520", "1.299", "1.077"},
    {"0.32", "1.787", "1.566", "1.345", "1.125", "0.904"}, {"≥0.64", "1.764", "1.543", "1.323", "1.102", "0.811"}}),
  n("Columnas de la tabla de Kirkham: 2r0/D."),
  h("Infiltración típica"),
  tab({"Suelo", "f [mm/h]"}, {{"Arena", ">30"}, {"Arena limosa", "20-30"}, {"Limo", "10-20"}, {"Arcilla limosa", "5-10"}, {"Arcilla", "1-5"}}),
})

-- ------------------------------------------------------------ 08
CAP("C08", "Drenes – régimen impermanente", {
  h("Ecuación de continuidad"),
  f([[h \frac{\partial^{2} h}{\partial x^{2}} + \p{\frac{\partial h}{\partial x}}^{2} = \frac{S}{K} \frac{\partial h}{\partial t}]]),
  n("Moody (1964): gráficos y/y0 y qD/(Kπy0) vs α = K H0 t/(S D²), con m = y0/H0 (d ≠ 0)."),
  h("Glover–Dumm (1954)"),
  f([[y = h - d \qquad y_{0} = H_{0} - d \qquad \alpha = \frac{\pi^{2} K d}{S D^{2}}]]),
  f([[y = 1.16 y_{0} e^{-\alpha t}]], "x = D/2 y αt ≥ 0.2"),
  f([[D = \pi \sqrt{\frac{K d t}{S \ln \p{\frac{1.16 y_{0}}{y}}}}]]),
  h("Corrección de Hooghoudt (χ = π r0)"),
  f([[\frac{d'}{d} = \frac{1}{1 + \frac{8}{\pi} \frac{d}{D} \ln \p{\frac{d}{\chi}}}]], "d/D < 0.25"),
  f([[\frac{d'}{d} = \frac{\pi}{2 \p{\ln \p{\frac{D}{\chi}} + 0.18}}]], "d/D ≥ 0.25"),
  b("Iterar: D con d → d' con D → D con d' → hasta converger."),
  h("Equilibrio dinámico"),
  f([[y_{0,i} = y_{f,i-1} + \frac{R_{i}}{S} \qquad y_{f,i} = 1.16 y_{0,i} e^{-\alpha t_{i}}]]),
  f([[y_{f} \le y_{0} \qquad I_{s} = \int (y > y_{a}) dt \le 200 mm \cdot dia]]),
})

-- ------------------------------------------------------------ 09
CAP("C09", "Drenaje longitudinal", {
  h("Cuneta y contrafoso"),
  b("Ancho w ≥ 0.5 m; i1 ≥ 8% (w > 0.5 m), máx. 30% (w = 0.5 m)."),
  b("Pendiente longitudinal mínima: 0.12% revestida; 0.25% sin revestir."),
  b("Bombeo de la calzada i2 ≈ 2–4%."),
  h("Caudal (A < 25 km², td > tc)"),
  f([[Q_{max} = \frac{C \bar{i}(t_{c}) A}{3.6}]], "Q [m³/s]; i [mm/h]; A [km²]; verificar tc ≥ 10 min"),
  tab({"Autor", "tc [min]"}, {{"Normas Esp.", "18 L^0.76/S^0.19"}, {"California", "57 (L³/H)^0.385"}, {"Giandotti", "60(4√A+1.5L)/(0.8√Hm)"},
    {"SCS", "3.42 L^0.8 (1000/CN−9)^0.7/S^0.5"}}),
  n("L [km], S [%], H y Hm [m], A [km²]."),
  h("Precipitación de diseño"),
  f([[i_{t} = i_{24} \sqrt{\frac{24}{t}}]], "Grunsky [mm/h]"),
  f([[P_{t}^{T} = K \cdot CD_{t} \cdot CF_{T} \cdot P_{D}^{10} \qquad K = 1.1]]),
  f([[P_{t}^{T} = (0.54 t^{0.25} - 0.50)(0.21 \ln T + 0.52) P_{1}^{10}]], "Bell (t < 1 h, t en min)"),
  h("Coeficiente de escorrentía"),
  tab({"Terreno", "C"}, {{"Adoquín", "0.50-0.70"}, {"Asfalto", "0.70-0.95"}, {"Concreto", "0.80-0.95"}, {"Arenoso con vegetación", "0.15-0.20"},
    {"Arcilloso con pasto", "0.25-0.65"}, {"Cultivo", "0.20-0.40"}}),
  tab({"Factor", "Extremo", "Alto", "Normal", "Bajo"}, {{"Relieve", ".28-.35", ".20-.28", ".14-.20", ".08-.14"}, {"Infiltración", ".12-.16", ".08-.12", ".06-.08", ".04-.06"},
    {"Cobertura", ".12-.16", ".08-.12", ".06-.08", ".04-.06"}, {"Almacenam.", ".10-.12", ".08-.10", ".06-.08", ".04-.06"}}),
  n("Tabla para T = 10 años (sumar factores). T=25: ×1.10; T=50: ×1.20; T=100: ×1.25."),
  h("Altura de agua y velocidades"),
  f([[\frac{Q n}{\sqrt{S}} = \frac{\Omega^{5/3}}{\Psi^{2/3}}]], "Contrafoso de tierra n = 0.023–0.025"),
  tab({"Vmax [m/s]", "Suelo"}, {{"0.6", "Limo"}, {"0.9", "Arena"}, {"1.2", "Arcilla"}, {"3.0", "Asfalto"}, {"4.5", "Hormigón"}}),
  f([[v_{min} = 0.4 \sim 0.5 m/s \qquad R = 0.15 \sim 0.2 h]], "Smín: tierra 0.25%, asfalto 0.15%; R: revancha"),
  b("Contrafoso trapecial; verificar θ < φ (ángulo de fricción)."),
  h("Periodo de retorno"),
  f([[R = 1 - \b{1 - \frac{1}{T}}^{n}]]),
  tab({"Obra", "T dis.", "T ver.", "n", "R dis."}, {{"Carreteras", "10", "25", "10", "65%"}, {"Caminos", "5", "10", "5", "67%"}}),
  b("Bajadas: Ø ≥ 200 mm si H/V ≥ 4; gradería si H/V ≤ 4:1; v = 5 m/s."),
})

-- ------------------------------------------------------------ 10
CAP("C10", "Alcantarillas", {
  h("Generalidades"),
  b("Luz ≤ 6 m (MC 3.703.101); mayor = puente. Libre entre tubos D/2."),
  b("Seguir línea de flujo y pendiente del cauce; muros guía."),
  b("Horacio Mery: operar cercano al crítico, Fr 0.7–0.9."),
  h("Criterio 1: control de entrada (torrente)"),
  f([[H_{c} = \frac{3}{2} h_{c} \qquad h_{c} = \sqrt[3]{\frac{Q^{2}}{g b^{2}}}]], "Rectangular; se impone D = H = Hc"),
  f([[D = \frac{3}{2} \frac{1}{g^{1/3}} \p{\frac{Q}{b}}^{2/3}]]),
  f([[S_{c} = \frac{2}{3} g H n^{2} \p{\frac{b + \frac{4}{3}H}{\frac{2}{3} H b}}^{4/3}]], "SD ≥ Sc ⇒ torrente; si no, SD = Sc"),
  f([[Q = 1.425 D^{5/2} \qquad S_{c} = \frac{31.14 n^{2}}{D^{1/3}}]], "Circular (hc = 0.7 D)"),
  b("Verificar v < 5 m/s y diseñar con D > H."),
  h("Criterio 2: verificación (T2)"),
  f([[S_{n} = \frac{10.3 Q_{2D}^{2} n^{2}}{D^{16/3}} \qquad \frac{Q_{2D} n}{\sqrt{S_{n}}} = \Omega R_{H}^{2/3}]], "Circular | Rectangular (hn = D)"),
  f([[H' = c_{c} D + \frac{1}{2g} \p{\frac{Q_{2D}}{b c_{c} D}}^{2}]], "Caso 1: SD > Sn (compuerta), cc ≈ 0.611"),
  f([[H' = D + \frac{1}{2g}(k_{e} + 1) \frac{Q_{2D}^{2}}{(D b)^{2}}]], "Caso 2: SD = Sn"),
  f([[H' = D + \frac{v_{2}^{2}}{2g} + \Lambda_{e} + J L - S_{D} L]], "Caso 3: SD < Sn"),
  f([[H' = D + \frac{v_{2}^{2}}{2g} + \Lambda_{e} + \Lambda_{s} + J L - S_{D} L]], "Control de salida"),
  f([[\Lambda_{e} = k_{e} \frac{v_{2}^{2}}{2g} \quad \Lambda_{s} = k_{s} \frac{(v_{2} - v_{3})^{2}}{2g} \quad J = \frac{v_{2}^{2} n^{2}}{R_{H}^{4/3}}]], "ks = 0.22–0.44 (0.4)"),
  tab({"Entrada", "ke"}, {{"Tubo saliente", "0.9"}, {"A ras de muro", "0.5"}, {"Abocinada", "0.2"}}),
  h("Periodos de retorno (MC V3)"),
  tab({"Obra", "Ruta", "Dis.", "Ver.", "n"}, {{"Puentes", "Carretera", "200", "300", "50"}, {"", "Camino", "100", "150", "50"},
    {"Alc. S>1.75 m²", "Carretera", "100", "150", "50"}, {"", "Camino", "50", "100", "30"}, {"Alc. S<1.75 m²", "Carretera", "50", "100", "50"}, {"", "Camino", "25", "50", "30"}}),
  h("Carga de diseño He"),
  tab({"Cauce", "Tubos", "Cajones", "Losas"}, {{"Canales", "D", "H", "H−0.10"}, {"Diseño natural", "D+0.3", "H+0.3", "H−0.10"}, {"Verif. natural", "D+0.6", "H+0.6", "H"}}),
  n("He máximo: cota exterior del SAP − 0.3 m."),
  h("Secciones"),
  tab({"Sección", "A", "P", "B"}, {{"Rect.", "b y", "b+2y", "b"}, {"Trapecio", "(b+ky)y", "b+2y√(1+k²)", "b+2ky"}, {"Triáng.", "k y²", "2y√(1+k²)", "2ky"},
    {"Circular", "(θ−senθ)D²/8", "θD/2", "D sen(θ/2)"}}),
  f([[\theta = 2 \text{arccos} \p{1 - \frac{2y}{D}}]]),
})

-- ------------------------------------------------------------ 11
CAP("C11", "Puentes", {
  h("Recomendaciones"),
  b("No generar resaltos; no peraltar aguas arriba; evitar socavación."),
  b("Opción 1: terraplén y estribos (fusible). Opción 2: cepas."),
  h("Opción 1: largo mínimo L (momenta)"),
  f([[M_{1} = \frac{Q^{2}}{g (b h_{1})} + \frac{h_{1}}{2} (b h_{1})]]),
  f([[M_{2} = \frac{Q^{2}}{g (L h_{2})} + \frac{h_{2}}{2} (b h_{2})]]),
  f([[h_{c2,M} = \sqrt[3]{\frac{Q^{2}}{g b L}}]], "Condición: M1 = M2c ⇒ L ≥ Lmin"),
  h("Energía"),
  f([[h_{2} + \frac{v_{2}^{2}}{2g} = h_{n1} + \frac{v_{1}^{2}}{2g} + \Lambda_{s}]], "Entre (1) y (2)"),
  f([[\Lambda_{s} = \frac{(v_{2} - v_{1})^{2}}{2g} - \frac{(h_{n1} - h_{2})^{2}}{2 h_{n1}}]], "Borda"),
  f([[h_{3} + \frac{v_{3}^{2}}{2g} = h_{2} + \frac{v_{2}^{2}}{2g} + \Lambda_{e} , \quad \Lambda_{e} = k_{e} \frac{v_{2}^{2}}{2g}]], "ke = 0.1–0.3; verificar h3 < hmax y v2 < vmax"),
  f([[\frac{1}{X_{2}} - \frac{1}{n X_{3}} = \frac{1}{2}(X_{3}^{2} - X_{2}^{2})]], "Momentum 3–2: X = h/hc2, n = b3/b2"),
  h("Opción 2: cepas (Yarnell)"),
  f([[\frac{h_{3} - h_{1}}{h_{1}} = K Fr_{1}^{2} (K + 5 Fr_{1}^{2} - 0.6)(\sigma + 15 \sigma^{4})]], "σ = 1 − be/b ; be = b − nD"),
  tab({"Forma de la cepa", "K"}, {{"Nariz y cola semicircular", "0.90"}, {"Doble cilindro con diafragma", "0.95"}, {"Doble cilindro sin diafragma", "1.05"},
    {"Triangular 90°", "1.05"}, {"Cuadrada", "1.25"}, {"Trestle bent 10 pilotes", "2.50"}}),
})

-- ------------------------------------------------------------ 12
CAP("C12", "Drenaje urbano", {
  h("Cunetas"),
  f([[Q = \frac{\sqrt{S}}{n} \frac{b^{2} i}{2} \p{\frac{b i}{2(i + 1)}}^{2/3}]], "Simple: y = b i ≤ 15 cm, b ≤ 1–2 m, i = 2–4%"),
  f([[Q = 0.315 \frac{\sqrt{S}}{n} b^{8/3} \b{i_{2} + (i_{1} - i_{2})\p{\frac{w}{b}}^{2}}^{5/3}]], "Compuesta: w ≈ 0.6 m; i1 = 2–10%"),
  h("Sumidero horizontal de fondo"),
  f([[Q_{m} = 1.66 (L_{e} + 2 w_{e}) h^{1.5}]], "Vertedero: h < 1.6 Ae/(Le + 2we)"),
  f([[Q_{m} = 2.66 A_{e} h^{0.5}]], "Orificio: h ≥ 1.6 Ae/(Le + 2we)"),
  f([[w_{e} = w - e n_{L} \quad L_{e} = L - e n_{T} \quad A_{e} = L_{e} w_{e}]], "w = 0.4–0.7 m; L = 1–2 m"),
  f([[\eta_{H} = R_{f} \eta_{0} + R_{s}(1 - \eta_{0}) \qquad \eta_{0} = 1 - \p{1 - \frac{w}{b_{s}}}^{2.67}]]),
  f([[R_{s} = \frac{1}{1 + \frac{0.0828 v^{1.8}}{i L^{2.3}}} \qquad R_{f} = 1 - 0.295(v - v_{0})]]),
  h("Sumidero lateral (S < 3%)"),
  f([[Q_{m} = 1.27 L h^{1.5} (h < a) \qquad Q_{m} = 2.66 L a h^{0.5} (h \ge a)]]),
  f([[\eta_{L} = 1 - \p{1 - \frac{L}{L_{T}}}^{1.8} \qquad L_{T} = 0.81 Q^{0.42} S^{0.3} (n i)^{-0.6}]], "ηL = 1 si h > a; mixto: (ηH + ηL)Q"),
  h("Factor de seguridad (ηd = F η)"),
  tab({"Caudal", "Mantención", "F"}, {{"< 30 L/s", "buena", "1.0"}, {"", "mala", "0.8"}, {"> 30 L/s", "buena", "0.7"}, {"", "mala", "0.6"}}),
  h("Separación entre sumideros"),
  f([[L = \frac{3600 (F \eta) Q_{max}}{C i b_{c}}]], "i (Tr 2 años, 30 min); bc = 0.5e + wc (doble bombeo) o e + wc"),
  b("Tormentas menores (Tr 2–10): ancho recomendado, tc ≥ 30 min."),
  b("Tormentas mayores (Tr 50–100): bajo solera y h < 15 cm."),
  h("Colectores"),
  b("Ø 300 mm a 2 m; h ≤ 0.8D; D > 0.8 m primario (DOH), ≤ 0.8 m SERVIU."),
  f([[V_{c} = 0.397 \frac{D^{2/3} \sqrt{S}}{n} \ge 0.6 m/s]], "Autolavado"),
  tab({"Situación", "h", "θ", "A", "Rh"}, {{"Llena", "D", "360°", "0.785D²", "0.250D"}, {"Q máx", "0.95D", "308°", "0.771D²", "0.287D"}, {"Diseño", "0.80D", "254°", "0.674D²", "0.304D"}}),
  tab({"Tubería", "vmax [m/s]"}, {{"Concreto ≤ 45 cm", "3.0"}, {"Concreto armado ≥ 61 cm", "3.5"}, {"Fibrocemento / PVC / HDPE", "5.0"},
    {"Acero corrugado", "4.5"}, {"Acero con shotcrete", "3.0"}, {"PRFV", "3.5"}}),
  b("Profundidad: aguas lluvias 1.8 m; agua potable 1–1.5 m. Cámaras cada 50–130 m y trampas de arena."),
  h("Caudales en red (racional)"),
  f([[Q_{D2} = \text{max}(Q_{D1}, Q_{2}, Q_{2I}) \qquad Q_{2I} = \frac{\bar{C}_{1,2} \bar{i}(t_{c2I})(A_{1} + A_{2})}{3.6}]]),
  f([[t_{c2I} = \text{max}(t_{c1} + t_{t1}, t_{c2}) \qquad t_{t1} = \frac{L_{1}}{v_{1}}]]),
  h("Hidrograma SCS"),
  f([[t_{p} = \frac{t_{LL}}{2} + 0.6 t_{c} \quad t_{r} = 1.67 t_{p} \quad t_{B} = 2.67 t_{p}]]),
  f([[Q_{p} = \frac{2 P_{ef} A}{2.67 \p{\frac{t_{LL}}{2} + 0.6 t_{c}}}]]),
})

-- ------------------------------------------------------------ constantes
CAP("REF", "Constantes y unidades", {
  h("Constantes"),
  f([[g = 9.81 m/s^{2} \qquad \gamma_{w} = 9810 N/m^{3} \qquad \mu = 10^{-3} Pa \cdot s]]),
  h("Conversiones"),
  tab({"De", "A", "Factor"}, {{"cm/s", "m/s", "×0.01"}, {"m/día", "m/s", "÷86400"}, {"mm/día", "m/s", "÷8.64E7"}, {"L/s", "m³/s", "÷1000"},
    {"km²", "m²", "×1E6"}, {"ha", "m²", "×1E4"}, {"darcy", "m²", "9.87E−13"}, {"log10", "ln", "ln = 2.303 log"}}),
  h("Racional con otras unidades"),
  f([[Q[m^{3}/s] = \frac{C i[mm/h] A[km^{2}]}{3.6} = \frac{C i[mm/h] A[ha]}{360}]]),
})

-- =====================================================================
--  VISOR (solo lectura)
--  Bloques:  {"h", titulo}   {"p", texto}   {"b", vineta}
--            {"f", formula [, nota]}   {"n", nota}   {"t", {encabezado}, {filas...}}
-- =====================================================================
local C = {
  azul = {0, 69, 137}, azul2 = {28, 96, 170}, amar = {247, 176, 0}, amarc = {255, 243, 205},
  celes = {226, 238, 250}, celesb = {110, 160, 215}, gris = {100, 100, 100}, grisc = {238, 238, 238},
  rojo = {200, 25, 35}, negro = {20, 20, 20}, blanco = {255, 255, 255},
}
local scr = "indice"
local cap, isel, itop = 1, 1, 1
local scroll = 0
local layout, layoutKey
local function W_() return platform.window:width() end
local function H_() return platform.window:height() end
local function col(gc, c) gc:setColorRGB(c[1], c[2], c[3]) end
local function inv() platform.window:invalidate() end
local HEAD = 22

local function wrapText(gc, text, width, font, style, size)
  gc:setFont(font, style, size)
  local out, line = {}, ""
  for word in text:gmatch("%S+") do
    local t = (line == "") and word or (line .. " " .. word)
    if gc:getStringWidth(t) > width and line ~= "" then out[#out + 1] = line; line = word else line = t end
  end
  if line ~= "" then out[#out + 1] = line end
  return out
end
local function roundRect(gc, x, y, w, h, r, fill)
  if fill then
    gc:fillRect(x + r, y, w - 2*r, h); gc:fillRect(x, y + r, w, h - 2*r)
    gc:fillArc(x, y, 2*r, 2*r, 90, 90); gc:fillArc(x + w - 2*r, y, 2*r, 2*r, 0, 90)
    gc:fillArc(x, y + h - 2*r, 2*r, 2*r, 180, 90); gc:fillArc(x + w - 2*r, y + h - 2*r, 2*r, 2*r, 270, 90)
  else
    gc:drawLine(x + r, y, x + w - r, y); gc:drawLine(x + r, y + h, x + w - r, y + h)
    gc:drawLine(x, y + r, x, y + h - r); gc:drawLine(x + w, y + r, x + w, y + h - r)
    gc:drawArc(x, y, 2*r, 2*r, 90, 90); gc:drawArc(x + w - 2*r, y, 2*r, 2*r, 0, 90)
    gc:drawArc(x, y + h - 2*r, 2*r, 2*r, 180, 90); gc:drawArc(x + w - 2*r, y + h - 2*r, 2*r, 2*r, 270, 90)
  end
end

-- ------------------------------------------- construccion del layout
-- cada elemento: {y, h, draw = function(gc, y0)}
local function buildLayout(gc)
  local key = cap .. ":" .. MATH.zoom() .. ":" .. W_()
  if layoutKey == key then return layout end
  MATH.setgc(gc)
  local Z = MATH.zoom()
  local base = ({9, 10, 11, 12})[Z]
  local L, y = {}, 6
  local M = 8
  local width = W_() - 2*M - 6
  local doc = DOC[cap]
  for _, bl in ipairs(doc.blocks) do
    local k = bl[1]
    if k == "h" then
      local lines = wrapText(gc, bl[2], width, "sansserif", "b", base + 1)
      local h = #lines*(base + 6) + 8
      local yy = y + 4
      L[#L + 1] = {y = yy, h = h, draw = function(g, y0)
        col(g, C.azul); g:setFont("sansserif", "b", base + 1)
        for j, l in ipairs(lines) do g:drawString(l, M, y0 + (j - 1)*(base + 6), "top") end
        col(g, C.amar); g:fillRect(M, y0 + #lines*(base + 6) + 1, math.min(width, 60), 2)
      end}
      y = yy + h
    elseif k == "p" or k == "b" then
      local ind = (k == "b") and 12 or 0
      local lines = wrapText(gc, bl[2], width - ind, "sansserif", "r", base)
      local lh = base + 5
      L[#L + 1] = {y = y, h = #lines*lh + 3, draw = function(g, y0)
        g:setFont("sansserif", "r", base)
        if k == "b" then col(g, C.rojo); g:fillArc(M + 2, y0 + base/2 - 1, 5, 5, 0, 360) end
        col(g, C.negro)
        for j, l in ipairs(lines) do g:drawString(l, M + ind, y0 + (j - 1)*lh, "top") end
      end}
      y = y + #lines*lh + 3
    elseif k == "f" then
      local z0 = MATH.zoom()
      local box = MATH.box(bl[2])
      while box.w + 10 > width and MATH.zoom() > 1 do MATH.zoom(MATH.zoom() - 1); box = MATH.box(bl[2]) end
      local zb = MATH.zoom()
      MATH.zoom(z0)
      local bw = box.w + 16
      local h = box.a + box.d + 12
      local nl = bl[3] and wrapText(gc, bl[3], width - 12, "sansserif", "i", base - 1) or {}
      local hn = #nl*(base + 3)
      local tot = h + hn + (hn > 0 and 4 or 0)
      local yy = y + 3
      L[#L + 1] = {y = yy, h = tot, w = bw, draw = function(g, y0, xoff)
        local x0 = M
        col(g, C.amarc); roundRect(g, x0, y0, width + 4, tot, 5, true)
        col(g, C.amar); roundRect(g, x0, y0, width + 4, tot, 5, false)
        col(g, C.negro); MATH.setgc(g)
        local xm = math.max(x0 + 4, x0 + (width + 4 - box.w)/2)
        local zc = MATH.zoom(); MATH.zoom(zb)
        MATH.draw(box, xm, y0 + 6 + box.a)
        MATH.zoom(zc)
        if hn > 0 then
          col(g, C.gris); g:setFont("sansserif", "i", base - 1)
          for j, l in ipairs(nl) do g:drawString(l, x0 + 8, y0 + h + (j - 1)*(base + 3), "top") end
        end
      end, wide = box.w + 8 >= width, bw = box.w}
      y = yy + tot + 4
    elseif k == "n" then
      local lines = wrapText(gc, bl[2], width - 14, "sansserif", "r", base - 1)
      local lh = base + 3
      local h = #lines*lh + 8
      L[#L + 1] = {y = y + 2, h = h, draw = function(g, y0)
        col(g, C.celes); roundRect(g, M, y0, width + 4, h, 4, true)
        col(g, C.celesb); g:fillRect(M, y0, 3, h)
        col(g, C.negro); g:setFont("sansserif", "r", base - 1)
        for j, l in ipairs(lines) do g:drawString(l, M + 9, y0 + 4 + (j - 1)*lh, "top") end
      end}
      y = y + h + 6
    elseif k == "t" then
      local hdr, rows = bl[2], bl[3]
      local ncol = #hdr
      local fs = base - 1
      gc:setFont("sansserif", "r", fs)
      local cw = {}
      for c = 1, ncol do
        gc:setFont("sansserif", "b", fs); cw[c] = gc:getStringWidth(hdr[c]) + 8
        gc:setFont("sansserif", "r", fs)
        for _, r in ipairs(rows) do cw[c] = math.max(cw[c], gc:getStringWidth(r[c] or "") + 8) end
      end
      local tw = 0; for c = 1, ncol do tw = tw + cw[c] end
      if tw > width + 4 then local f = (width + 4)/tw; for c = 1, ncol do cw[c] = cw[c]*f end; tw = width + 4 end
      local rh = fs + 6
      local h = rh*(#rows + 1) + 2
      L[#L + 1] = {y = y + 3, h = h, draw = function(g, y0)
        local x0 = M + (width + 4 - tw)/2
        col(g, C.azul2); g:fillRect(x0, y0, tw, rh)
        local xx = x0
        g:setFont("sansserif", "b", fs); col(g, C.blanco)
        for c = 1, ncol do g:drawString(hdr[c], xx + 4, y0 + 2, "top"); xx = xx + cw[c] end
        g:setFont("sansserif", "r", fs)
        for r, row in ipairs(rows) do
          local yr = y0 + rh*r
          if r % 2 == 0 then col(g, C.grisc); g:fillRect(x0, yr, tw, rh) end
          col(g, C.negro); xx = x0
          for c = 1, ncol do g:drawString(row[c] or "", xx + 4, yr + 2, "top"); xx = xx + cw[c] end
        end
        col(g, C.gris); g:drawRect(x0, y0, tw, rh*(#rows + 1))
      end}
      y = y + h + 8
    end
  end
  layout = {items = L, total = y + 10}
  layoutKey = key
  return layout
end

-- ---------------------------------------------------------------- dibujo
local function header(gc, title, sub)
  col(gc, C.azul); gc:fillRect(0, 0, W_(), HEAD)
  col(gc, C.amar); gc:fillRect(0, HEAD, W_(), 2)
  col(gc, C.blanco); gc:setFont("sansserif", "b", 11)
  gc:drawString(title, 6, 3, "top")
  if sub then gc:setFont("sansserif", "r", 8); local w = gc:getStringWidth(sub); gc:drawString(sub, W_() - w - 6, 6, "top") end
end

local function drawIndex(gc)
  header(gc, "Formulario CIV-346", "solo lectura")
  local lh = 20
  local n = math.floor((H_() - HEAD - 8)/lh)
  if isel < itop then itop = isel end
  if isel > itop + n - 1 then itop = isel - n + 1 end
  for k = 0, n - 1 do
    local i = itop + k
    local d = DOC[i]
    if not d then break end
    local y = HEAD + 5 + k*lh
    if i == isel then col(gc, C.amarc); roundRect(gc, 4, y, W_() - 8, lh - 2, 4, true); col(gc, C.amar); roundRect(gc, 4, y, W_() - 8, lh - 2, 4, false) end
    col(gc, C.azul); gc:setFont("sansserif", "b", 10); gc:drawString(d.num or "", 10, y + 2, "top")
    col(gc, C.negro); gc:setFont("sansserif", "r", 10); gc:drawString(d.title, 42, y + 2, "top")
  end
  -- barra de desplazamiento
  if #DOC > n then
    local h = (H_() - HEAD - 8)*n/#DOC
    local yb = HEAD + 5 + (H_() - HEAD - 8 - h)*(itop - 1)/math.max(1, #DOC - n)
    col(gc, C.grisc); gc:fillRect(W_() - 4, HEAD + 4, 3, H_() - HEAD - 8)
    col(gc, C.azul2); gc:fillRect(W_() - 4, yb, 3, h)
  end
end

local function drawPage(gc)
  local d = DOC[cap]
  header(gc, (d.num and (d.num .. "  ") or "") .. d.title, cap .. "/" .. #DOC)
  local lay = buildLayout(gc)
  local viewH = H_() - HEAD - 2
  local maxs = math.max(0, lay.total - viewH)
  if scroll > maxs then scroll = maxs end
  if scroll < 0 then scroll = 0 end
  gc:clipRect("set", 0, HEAD + 2, W_(), viewH)
  for _, it in ipairs(lay.items) do
    local y0 = HEAD + 2 + it.y - scroll
    if y0 + it.h >= HEAD and y0 <= H_() then it.draw(gc, y0) end
  end
  gc:clipRect("reset")
  if maxs > 0 then
    local h = math.max(12, viewH*viewH/lay.total)
    local yb = HEAD + 2 + (viewH - h)*scroll/maxs
    col(gc, C.grisc); gc:fillRect(W_() - 4, HEAD + 2, 3, viewH)
    col(gc, C.azul2); gc:fillRect(W_() - 4, yb, 3, h)
  end
end

function on.paint(gc)
  gc:setPen("thin", "smooth")
  col(gc, C.blanco); gc:fillRect(0, 0, W_(), H_())
  if scr == "indice" then drawIndex(gc) else drawPage(gc) end
end

-- ---------------------------------------------------------------- teclas
local function openCap(i) cap = math.max(1, math.min(#DOC, i)); scroll = 0; scr = "pagina" end
function on.arrowUp() if scr == "indice" then isel = (isel - 2) % #DOC + 1 else scroll = scroll - 24 end; inv() end
function on.arrowDown() if scr == "indice" then isel = isel % #DOC + 1 else scroll = scroll + 24 end; inv() end
function on.arrowLeft() if scr == "pagina" then scroll = scroll - (H_() - HEAD - 30) end; inv() end
function on.arrowRight() if scr == "pagina" then scroll = scroll + (H_() - HEAD - 30) end; inv() end
function on.tabKey() if scr == "pagina" then openCap(cap % #DOC + 1) else isel = isel % #DOC + 1 end; inv() end
function on.backtabKey() if scr == "pagina" then openCap((cap - 2) % #DOC + 1) else isel = (isel - 2) % #DOC + 1 end; inv() end
function on.enterKey() if scr == "indice" then openCap(isel) else scr = "indice"; isel = cap end; inv() end
on.returnKey = on.enterKey
function on.escapeKey() scr = "indice"; isel = cap; inv() end
on.backspaceKey = on.escapeKey
function on.charIn(ch)
  if ch == "+" then MATH.zoom(math.min(MATH.NZ, MATH.zoom() + 1)); layoutKey = nil
  elseif ch == "-" or ch == "\226\136\146" then MATH.zoom(math.max(1, MATH.zoom() - 1)); layoutKey = nil
  elseif scr == "indice" and tonumber(ch) then local n = tonumber(ch); if n == 0 then n = 10 end; if DOC[n] then isel = n; openCap(n) end end
  inv()
end
function on.resize() layoutKey = nil end
if toolpalette and toolpalette.register then
  local caps = {"Capitulos"}
  for i, d in ipairs(DOC) do caps[#caps + 1] = {(d.num and (d.num .. " ") or "") .. d.title, function() openCap(i); inv() end} end
  pcall(toolpalette.register, {caps, {"Zoom", {"Aumentar (+)", function() on.charIn("+") end}, {"Reducir (-)", function() on.charIn("-") end}}})
end
VISOR = {build = buildLayout, open = openCap, paint = on.paint}
