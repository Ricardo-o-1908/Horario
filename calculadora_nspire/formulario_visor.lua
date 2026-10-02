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
