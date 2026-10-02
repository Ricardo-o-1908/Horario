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
