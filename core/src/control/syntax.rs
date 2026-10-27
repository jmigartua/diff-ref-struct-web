//! The syntax of the control file (L066): a strict subset of TOML 1.0, so that any TOML reader reads our files and the
//! editor's TOML grammar colours them, parsed here into a document that keeps every comment, the user's order and the
//! exact place of every key and value.
//!
//! The subset: comments (`#` to the end of the line); tables `[a.b]` and arrays of tables `[[a.b]]`; `key = value`
//! with bare keys (`A–Z a–z 0–9 _ -`) or quoted keys; values: basic strings `"…"` (escapes `\" \\ \n \t \r \uXXXX`),
//! literal strings `'…'`, integers, floats (`1.5`, `-2e-3`, `inf`, `nan`), booleans, arrays (over several lines,
//! with comments and a trailing comma) and inline tables `{ k = v, … }`. Not used: dotted keys in a key/value line,
//! multi-line strings, dates; each is reported where it occurs.
//!
//! Every error has its line, column and byte span and a plain message; the parser goes on after an error (to the next
//! line), so that the editor shows them all at once.

use std::fmt;

/// A place in the text: bytes `start..end`, and the line and column (both from 1, the column in characters) of `start`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Span {
    pub start: usize,
    pub end: usize,
    pub line: usize,
    pub column: usize,
}

/// An error or a warning, where it is: a stable code (`syntax`, `unknown-name`, `value`, `missing`; `check` maps them,
/// L076) and, for a name to replace, the replacement text for the span.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Diagnostic {
    pub span: Span,
    pub message: String,
    pub warning: bool,
    pub code: &'static str,
    pub replacement: Option<String>,
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "line {}, column {}: {}{}", self.span.line, self.span.column, if self.warning { "warning: " } else { "" }, self.message)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    String(String),
    Integer(i64),
    Float(f64),
    Bool(bool),
    Array(Vec<Item>),
    Table(Vec<Entry>),
}

impl Value {
    pub fn kind(&self) -> &'static str {
        match self {
            Value::String(_) => "a string",
            Value::Integer(_) => "an integer",
            Value::Float(_) => "a number",
            Value::Bool(_) => "true or false",
            Value::Array(_) => "an array",
            Value::Table(_) => "a table",
        }
    }

    /// The same value, spans aside.
    pub fn same(&self, other: &Value) -> bool {
        match (self, other) {
            (Value::Array(a), Value::Array(b)) => a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.value.same(&y.value)),
            (Value::Table(a), Value::Table(b)) => a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.key == y.key && x.value.value.same(&y.value.value)),
            (Value::Float(a), Value::Float(b)) => a.to_bits() == b.to_bits(),
            (a, b) => a == b,
        }
    }
}

/// A value and the place of its text.
#[derive(Clone, Debug, PartialEq)]
pub struct Item {
    pub value: Value,
    pub span: Span,
}

/// `key = value`: the key, its place, the value, and the bytes of the whole line (from its start to after its newline;
/// several lines for a multi-line array).
#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub key: String,
    pub key_span: Span,
    pub value: Item,
    pub lines: (usize, usize),
}

/// A table of the document: the root (no header), `[path]` or one element of `[[path]]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Table {
    pub path: Vec<String>,
    pub array: bool,
    /// The header's place (None for the root).
    pub header: Option<Span>,
    /// Bytes of the header line.
    pub header_lines: (usize, usize),
    pub entries: Vec<Entry>,
}

impl Table {
    pub fn get(&self, key: &str) -> Option<&Entry> {
        self.entries.iter().find(|e| e.key == key)
    }

    /// The last byte of the table: after its last entry, or after its header.
    pub fn end(&self) -> usize {
        self.entries.last().map_or(self.header_lines.1, |e| e.lines.1)
    }

    pub fn name(&self) -> String {
        if self.array {
            format!("[[{}]]", self.path.join("."))
        } else if self.path.is_empty() {
            "the top of the file".into()
        } else {
            format!("[{}]", self.path.join("."))
        }
    }
}

/// A parsed control file: the text, its tables in order (the root first), and what is wrong with it.
#[derive(Clone, Debug, PartialEq)]
pub struct Document {
    pub source: String,
    pub tables: Vec<Table>,
    pub diagnostics: Vec<Diagnostic>,
}

/// An entry of the outline: a table, its title and the lines it spans (for the editor's outline and folding).
#[derive(Clone, Debug, PartialEq)]
pub struct OutlineItem {
    pub title: String,
    pub path: String,
    pub first_line: usize,
    pub last_line: usize,
}

impl Document {
    pub fn has_errors(&self) -> bool {
        self.diagnostics.iter().any(|d| !d.warning)
    }

    /// The line (from 1) of a byte.
    pub fn line_of(&self, byte: usize) -> usize {
        self.source[..byte.min(self.source.len())].matches('\n').count() + 1
    }

    /// Tables and the lines they span; a table named by its `name` key when it has one.
    pub fn outline(&self) -> Vec<OutlineItem> {
        let mut out = Vec::new();
        for (k, t) in self.tables.iter().enumerate() {
            let Some(h) = t.header else { continue };
            let next = self.tables.get(k + 1).and_then(|n| n.header).map_or(self.source.len(), |s| s.start);
            let last_byte = t.end().max(h.end).min(next).saturating_sub(1);
            let label = t.get("name").or_else(|| t.get("label")).and_then(|e| match &e.value.value {
                Value::String(s) => Some(s.clone()),
                _ => None,
            });
            out.push(OutlineItem {
                title: match label {
                    Some(l) => format!("{} {l}", t.path.join(".")),
                    None => t.path.join("."),
                },
                path: t.path.join("."),
                first_line: h.line,
                last_line: self.line_of(last_byte),
            });
        }
        out
    }
}

struct Parser<'a> {
    s: &'a str,
    b: &'a [u8],
    pos: usize,
    diags: Vec<Diagnostic>,
}

type PResult<T> = Result<T, Diagnostic>;

fn bare_key_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'-'
}

impl<'a> Parser<'a> {
    fn span_at(&self, start: usize, end: usize) -> Span {
        // Line and column of `start`, which is on the current line or before it.
        let line_start = self.s[..start].rfind('\n').map_or(0, |i| i + 1);
        let line = self.s[..start].matches('\n').count() + 1;
        Span { start, end, line, column: self.s[line_start..start].chars().count() + 1 }
    }

    fn err<T>(&self, start: usize, end: usize, message: impl Into<String>) -> PResult<T> {
        Err(Diagnostic { span: self.span_at(start, end.max(start)), message: message.into(), warning: false, code: "syntax", replacement: None })
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.pos).copied()
    }

    fn bump(&mut self) {
        self.pos += 1;
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t')) {
            self.pos += 1;
        }
    }

    fn skip_comment(&mut self) {
        if self.peek() == Some(b'#') {
            while !matches!(self.peek(), None | Some(b'\n')) {
                self.pos += 1;
            }
        }
    }

    /// Blanks, comments and newlines (inside arrays).
    fn skip_all(&mut self) {
        loop {
            self.skip_ws();
            self.skip_comment();
            match self.peek() {
                Some(b'\n') => self.bump(),
                Some(b'\r') if self.b.get(self.pos + 1) == Some(&b'\n') => self.pos += 1,
                _ => return,
            }
        }
    }

    /// The end of a line: blanks, a comment, then a newline or the end of the text.
    fn end_of_line(&mut self) -> PResult<()> {
        self.skip_ws();
        self.skip_comment();
        if self.peek() == Some(b'\r') && self.b.get(self.pos + 1) == Some(&b'\n') {
            self.pos += 1;
        }
        match self.peek() {
            None => Ok(()),
            Some(b'\n') => {
                self.bump();
                Ok(())
            }
            Some(_) => {
                let end = self.s[self.pos..].find('\n').map_or(self.s.len(), |i| self.pos + i);
                self.err(self.pos, end, format!("unexpected \"{}\" after the value; a comment starts with #", self.s[self.pos..end].trim_end()))
            }
        }
    }

    /// Skip to the start of the next line (after an error).
    fn recover(&mut self) {
        while let Some(c) = self.peek() {
            self.bump();
            if c == b'\n' {
                return;
            }
        }
    }

    fn key(&mut self) -> PResult<(String, Span)> {
        let start = self.pos;
        match self.peek() {
            Some(b'"') => {
                let s = self.basic_string()?;
                Ok((s, self.span_at(start, self.pos)))
            }
            Some(b'\'') => {
                let s = self.literal_string()?;
                Ok((s, self.span_at(start, self.pos)))
            }
            Some(c) if bare_key_char(c) => {
                while self.peek().is_some_and(bare_key_char) {
                    self.pos += 1;
                }
                Ok((self.s[start..self.pos].to_string(), self.span_at(start, self.pos)))
            }
            _ => {
                let end = self.s[start..].find(['\n', '=']).map_or(self.s.len(), |i| start + i);
                self.err(start, end, "a key is expected here (letters, digits, _ and -, or a quoted key)")
            }
        }
    }

    fn basic_string(&mut self) -> PResult<String> {
        let start = self.pos;
        if self.s[self.pos..].starts_with("\"\"\"") {
            return self.err(start, start + 3, "multi-line strings are not used in this format");
        }
        self.pos += 1;
        let mut out = String::new();
        loop {
            let Some(c) = self.s[self.pos..].chars().next() else { return self.err(start, self.pos, "this string is not closed (a \" is missing)") };
            match c {
                '"' => {
                    self.pos += 1;
                    return Ok(out);
                }
                '\n' => return self.err(start, self.pos, "this string is not closed on its line (a \" is missing)"),
                '\\' => {
                    let e = self.s[self.pos + 1..].chars().next();
                    let (ch, len) = match e {
                        Some('"') => ('"', 2),
                        Some('\\') => ('\\', 2),
                        Some('n') => ('\n', 2),
                        Some('t') => ('\t', 2),
                        Some('r') => ('\r', 2),
                        Some('u') => {
                            let hex = self.s.get(self.pos + 2..self.pos + 6).unwrap_or("");
                            match u32::from_str_radix(hex, 16).ok().and_then(char::from_u32) {
                                Some(c) => (c, 6),
                                None => return self.err(self.pos, self.pos + 6, "\\u must be followed by four hexadecimal digits"),
                            }
                        }
                        _ => return self.err(self.pos, self.pos + 2, "unknown escape; use \\\" \\\\ \\n \\t \\r or \\uXXXX"),
                    };
                    out.push(ch);
                    self.pos += len;
                }
                c => {
                    out.push(c);
                    self.pos += c.len_utf8();
                }
            }
        }
    }

    fn literal_string(&mut self) -> PResult<String> {
        let start = self.pos;
        if self.s[self.pos..].starts_with("'''") {
            return self.err(start, start + 3, "multi-line strings are not used in this format");
        }
        self.pos += 1;
        match self.s[self.pos..].find(['\'', '\n']) {
            Some(i) if self.b[self.pos + i] == b'\'' => {
                let v = self.s[self.pos..self.pos + i].to_string();
                self.pos += i + 1;
                Ok(v)
            }
            _ => self.err(start, self.pos, "this string is not closed on its line (a ' is missing)"),
        }
    }

    fn number_or_word(&mut self) -> PResult<Value> {
        let start = self.pos;
        while self.peek().is_some_and(|c| c.is_ascii_alphanumeric() || matches!(c, b'+' | b'-' | b'.' | b'_')) {
            self.pos += 1;
        }
        let t = &self.s[start..self.pos];
        if t.is_empty() {
            let c = self.s[start..].chars().next().map_or("the end of the line".to_string(), |c| format!("\"{c}\""));
            return self.err(start, start + 1, format!("a value is expected here, not {c}"));
        }
        match t {
            "true" => return Ok(Value::Bool(true)),
            "false" => return Ok(Value::Bool(false)),
            "inf" | "+inf" => return Ok(Value::Float(f64::INFINITY)),
            "-inf" => return Ok(Value::Float(f64::NEG_INFINITY)),
            "nan" | "+nan" | "-nan" => return Ok(Value::Float(f64::NAN)),
            _ => {}
        }
        let digits = t.replace('_', "");
        let float = digits.contains(['.', 'e', 'E']);
        let bad = |p: &Parser| p.err(start, p.pos, format!("\"{t}\" is not a number, true, false or a quoted string"));
        if !digits.trim_start_matches(['+', '-']).starts_with(|c: char| c.is_ascii_digit()) || digits.ends_with('.') || digits.contains(".e") || digits.contains(".E") {
            if digits.starts_with('.') || digits.ends_with('.') {
                return self.err(start, self.pos, format!("\"{t}\": a number needs digits on both sides of its point (0.5, 2.0)"));
            }
            return bad(self);
        }
        if float {
            digits.parse::<f64>().map(Value::Float).or_else(|_| bad(self))
        } else {
            digits.parse::<i64>().map(Value::Integer).or_else(|_| bad(self))
        }
    }

    fn value(&mut self) -> PResult<Item> {
        let start = self.pos;
        let v = match self.peek() {
            Some(b'"') => Value::String(self.basic_string()?),
            Some(b'\'') => Value::String(self.literal_string()?),
            Some(b'[') => {
                self.pos += 1;
                let mut items = Vec::new();
                loop {
                    let before = self.pos;
                    self.skip_all();
                    self.open_array(before, start)?;
                    match self.peek() {
                        Some(b']') => {
                            self.pos += 1;
                            break;
                        }
                        None => return self.err(start, self.pos, "this array is not closed (a ] is missing)"),
                        _ => {}
                    }
                    items.push(self.value()?);
                    let before = self.pos;
                    self.skip_all();
                    self.open_array(before, start)?;
                    match self.peek() {
                        Some(b',') => self.pos += 1,
                        Some(b']') => {
                            self.pos += 1;
                            break;
                        }
                        _ => return self.err(self.pos, self.pos + 1, "a comma or a ] is expected after an element of an array"),
                    }
                }
                Value::Array(items)
            }
            Some(b'{') => {
                self.pos += 1;
                let mut entries: Vec<Entry> = Vec::new();
                self.skip_ws();
                if self.peek() == Some(b'}') {
                    self.pos += 1;
                } else {
                    loop {
                        self.skip_ws();
                        let ls = self.pos;
                        let (key, key_span) = self.key()?;
                        self.skip_ws();
                        if self.peek() != Some(b'=') {
                            return self.err(self.pos, self.pos + 1, format!("an = is expected after the key \"{key}\""));
                        }
                        self.pos += 1;
                        self.skip_ws();
                        let value = self.value()?;
                        if entries.iter().any(|e| e.key == key) {
                            return Err(Diagnostic { span: key_span, message: format!("\"{key}\" is given twice in this table"), warning: false, code: "syntax", replacement: None });
                        }
                        entries.push(Entry { key, key_span, value, lines: (ls, self.pos) });
                        self.skip_ws();
                        match self.peek() {
                            Some(b',') => self.pos += 1,
                            Some(b'}') => {
                                self.pos += 1;
                                break;
                            }
                            _ => return self.err(self.pos, self.pos + 1, "a comma or a } is expected in an inline table (on one line)"),
                        }
                    }
                }
                Value::Table(entries)
            }
            _ => self.number_or_word()?,
        };
        Ok(Item { value: v, span: self.span_at(start, self.pos) })
    }

    /// After the blanks between `before` and here: if they crossed a line and this line holds `key =` or a table
    /// header, the array that opened at `start` was left open; the error is placed there, and the parser stops at
    /// this line, which is read for itself.
    fn open_array(&mut self, before: usize, start: usize) -> PResult<()> {
        if !self.s[before..self.pos].contains('\n') {
            return Ok(());
        }
        let rest = &self.s[self.pos..self.s[self.pos..].find('\n').map_or(self.s.len(), |i| self.pos + i)];
        let keyval = rest.find('=').is_some_and(|i| !rest[..i].contains([',', ']', '[', '"', '{']));
        if keyval || (rest.starts_with('[') && rest.trim_end().ends_with(']') && !rest.contains(',')) {
            let line = self.span_at(start, start).line;
            let here = self.s[..self.pos].rfind('\n').map_or(0, |i| i + 1);
            let e = self.err(start, start + 1, format!("the array opened on line {line} is not closed (a ] is missing)"));
            self.pos = here;
            return e;
        }
        Ok(())
    }

    fn header(&mut self) -> PResult<(Vec<String>, bool, Span)> {
        let start = self.pos;
        let array = self.s[self.pos..].starts_with("[[");
        self.pos += if array { 2 } else { 1 };
        let mut path = Vec::new();
        loop {
            self.skip_ws();
            let (k, _) = self.key()?;
            path.push(k);
            self.skip_ws();
            match self.peek() {
                Some(b'.') => self.pos += 1,
                Some(b']') => break,
                _ => return self.err(self.pos, self.pos + 1, format!("a . or {} is expected in a table header", if array { "]]" } else { "]" })),
            }
        }
        let close = if array { "]]" } else { "]" };
        if !self.s[self.pos..].starts_with(close) {
            return self.err(start, self.pos + 1, format!("this header opens with {} and must close with {close}", if array { "[[" } else { "[" }));
        }
        self.pos += close.len();
        Ok((path, array, self.span_at(start, self.pos)))
    }
}

/// Parse a control file's text. The document is returned whatever the errors, with what could be read.
pub fn parse(text: &str) -> Document {
    let mut p = Parser { s: text, b: text.as_bytes(), pos: 0, diags: Vec::new() };
    let mut tables = vec![Table { path: Vec::new(), array: false, header: None, header_lines: (0, 0), entries: Vec::new() }];
    // A byte-order mark is skipped.
    if text.starts_with('\u{feff}') {
        p.pos = 3;
    }
    while p.pos < text.len() {
        let line_begin = p.pos;
        p.skip_ws();
        match p.peek() {
            None => break,
            Some(b'#' | b'\n' | b'\r') => {
                if let Err(d) = p.end_of_line() {
                    p.diags.push(d);
                    p.recover();
                }
                continue;
            }
            Some(b'[') => match p.header() {
                Ok((path, array, span)) => {
                    if !array && tables.iter().any(|t| !t.array && t.path == path && t.header.is_some()) {
                        p.diags.push(Diagnostic { span, message: format!("the table [{}] is given twice", path.join(".")), warning: false, code: "syntax", replacement: None });
                    }
                    if let Err(d) = p.end_of_line() {
                        p.diags.push(d);
                        p.recover();
                    }
                    tables.push(Table { path, array, header: Some(span), header_lines: (line_begin, p.pos), entries: Vec::new() });
                }
                Err(d) => {
                    p.diags.push(d);
                    p.recover();
                }
            },
            Some(_) => {
                let r = (|| -> PResult<Entry> {
                    let (key, key_span) = p.key()?;
                    p.skip_ws();
                    if p.peek() == Some(b'.') {
                        return p.err(key_span.start, p.pos + 1, format!("dotted keys are not used in this format: write the table [{key}…] and the key in it"));
                    }
                    if p.peek() != Some(b'=') {
                        return p.err(p.pos, p.pos + 1, format!("an = is expected after the key \"{key}\""));
                    }
                    p.pos += 1;
                    p.skip_ws();
                    if matches!(p.peek(), None | Some(b'\n' | b'#')) {
                        return p.err(p.pos, p.pos, format!("\"{key}\" has no value"));
                    }
                    let value = p.value()?;
                    p.end_of_line()?;
                    Ok(Entry { key, key_span, value, lines: (line_begin, p.pos) })
                })();
                match r {
                    Ok(e) => {
                        let t = tables.last_mut().unwrap();
                        if let Some(old) = t.entries.iter().find(|x| x.key == e.key) {
                            p.diags.push(Diagnostic { span: e.key_span, message: format!("\"{}\" is given twice in {} (first on line {})", e.key, t.name(), old.key_span.line), warning: false, code: "syntax", replacement: None });
                        } else {
                            t.entries.push(e);
                        }
                    }
                    Err(d) => {
                        p.diags.push(d);
                        // Go on at the next line; an array left open stops at the start of the line that is not
                        // part of it, which is read for itself.
                        let fresh = p.pos > line_begin && p.b.get(p.pos.wrapping_sub(1)) == Some(&b'\n');
                        if !fresh {
                            p.recover();
                        }
                    }
                }
            }
        }
    }
    Document { source: text.to_string(), tables, diagnostics: p.diags }
}

/// A float written so that it reads back to the same bits, as a TOML float (always with a point or an exponent).
pub fn float(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let a = v.abs();
    if a != 0.0 && !(1e-5..1e15).contains(&a) {
        return format!("{v:e}");
    }
    let s = format!("{v}");
    if s.contains('.') {
        s
    } else {
        format!("{s}.0")
    }
}

/// A string written as a TOML basic string.
pub fn string(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04X}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_entries_and_spans() {
        let t = "# head\nformat = \"diff-ref-struct 1\"\n\n[[phase]]\nname = 'Y2O3'  # the phase\ncell_a = 10.6046\nxyz = [0.25,\n  0.5, # half\n  -1e-3,\n]\n[[phase.atom]]\nok = true\n";
        let d = parse(t);
        assert!(d.diagnostics.is_empty(), "{:?}", d.diagnostics);
        assert_eq!(d.tables.len(), 3);
        let ph = &d.tables[1];
        assert_eq!(ph.path, ["phase"]);
        let a = ph.get("cell_a").unwrap();
        assert_eq!(a.value.value, Value::Float(10.6046));
        assert_eq!((a.value.span.line, a.value.span.column), (6, 10));
        assert_eq!(&t[a.value.span.start..a.value.span.end], "10.6046");
        let x = ph.get("xyz").unwrap();
        assert!(matches!(&x.value.value, Value::Array(v) if v.len() == 3));
        assert_eq!(d.outline()[0].title, "phase Y2O3");
    }

    #[test]
    fn errors_with_places_and_recovery() {
        let t = "a = 1\nb = .5\nc = \"open\nd = [1, 2\ne.f = 3\ng = 1 2\na = 2\n";
        let d = parse(t);
        let m: Vec<(usize, &str)> = d.diagnostics.iter().map(|x| (x.span.line, x.message.as_str())).collect();
        assert!(m.iter().any(|(l, s)| *l == 2 && s.contains("digits on both sides")), "{m:?}");
        assert!(m.iter().any(|(l, s)| *l == 3 && s.contains("not closed")), "{m:?}");
        assert!(m.iter().any(|(l, s)| *l == 4 && s.contains("opened on line 4 is not closed")), "{m:?}");
        assert!(m.iter().any(|(l, s)| *l == 5 && s.contains("dotted keys")), "{m:?}");
        assert!(m.iter().any(|(l, s)| *l == 6 && s.contains("after the value")), "{m:?}");
        assert!(m.iter().any(|(l, s)| *l == 7 && s.contains("given twice")), "{m:?}");
    }

    #[test]
    fn floats_read_back_to_the_bit() {
        for v in [0.1, 1.0, -2.5e-7, 1e300, 123456.789, 5e-324, 0.0, -0.0, 6.02214076e23] {
            let d = parse(&format!("x = {}\n", float(v)));
            assert!(matches!(d.tables[0].entries[0].value.value, Value::Float(f) if f.to_bits() == v.to_bits()), "{v}");
        }
    }
}
