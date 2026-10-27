-- diff-ref-struct Web · the manuals' tables in the PDF: each column as wide as its longest cell asks (a short label
-- narrow, a description wide), instead of the pipe table's equal dashes; the widths sum to the text width.
local function len(cell)
  local s = pandoc.utils.stringify(cell.contents or cell)
  return utf8.len(s) or #s
end

function Table(tbl)
  local n = #tbl.colspecs
  if n == 0 then return nil end
  local want = {}
  for i = 1, n do want[i] = 6 end
  local function measure(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= n then want[i] = math.max(want[i], math.min(len(cell), 70)) end
      end
    end
  end
  measure(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do measure(body.body) end
  local total = 0
  for i = 1, n do total = total + want[i] end
  for i = 1, n do tbl.colspecs[i] = { tbl.colspecs[i][1], 0.97 * want[i] / total } end
  return tbl
end
