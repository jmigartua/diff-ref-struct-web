//! A small XML reader for the data formats written in XML (PANalytical `.xrdml`, the `RawData*.xml` of Bruker
//! `.brml`, Rigaku's measurement conditions): elements, attributes and text, with namespaces dropped from names
//! ("ns:tag" is "tag"). Comments, processing instructions, CDATA and the five predefined entities and character
//! references are handled; a DTD is skipped, not read.

use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct XmlError(pub String);

impl fmt::Display for XmlError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&self.0)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Element {
    pub name: String,
    pub attributes: Vec<(String, String)>,
    pub children: Vec<Element>,
    pub text: String,
}

impl Element {
    pub fn attribute(&self, name: &str) -> Option<&str> {
        self.attributes.iter().find(|(k, _)| k == name).map(|(_, v)| v.as_str())
    }

    /// The first child with this name.
    pub fn child(&self, name: &str) -> Option<&Element> {
        self.children.iter().find(|c| c.name == name)
    }

    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Element> + 'a {
        self.children.iter().filter(move |c| c.name == name)
    }

    /// The first element with this name in document order below this one (this one included).
    pub fn find(&self, name: &str) -> Option<&Element> {
        if self.name == name {
            return Some(self);
        }
        self.children.iter().find_map(|c| c.find(name))
    }

    /// Every element with this name below this one, in document order.
    pub fn find_all<'a>(&'a self, name: &str, out: &mut Vec<&'a Element>) {
        if self.name == name {
            out.push(self);
        }
        for c in &self.children {
            c.find_all(name, out);
        }
    }

    /// The text of the first descendant with this name, trimmed.
    pub fn text_of(&self, name: &str) -> Option<&str> {
        self.find(name).map(|e| e.text.trim())
    }
}

fn local(name: &str) -> String {
    name.rsplit(':').next().unwrap_or(name).to_string()
}

fn unescape(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        rest = &rest[i..];
        let Some(j) = rest.find(';') else { break };
        let ent = &rest[1..j];
        let c = match ent {
            "lt" => Some('<'),
            "gt" => Some('>'),
            "amp" => Some('&'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ if ent.starts_with("#x") => u32::from_str_radix(&ent[2..], 16).ok().and_then(char::from_u32),
            _ if ent.starts_with('#') => ent[1..].parse().ok().and_then(char::from_u32),
            _ => None,
        };
        match c {
            Some(c) => out.push(c),
            None => out.push_str(&rest[..=j]),
        }
        rest = &rest[j + 1..];
    }
    out.push_str(rest);
    out
}

/// The root element of a document.
pub fn parse(text: &str) -> Result<Element, XmlError> {
    let t = text.strip_prefix('\u{feff}').unwrap_or(text);
    let err = |m: &str, at: usize| XmlError(format!("XML: {m} at byte {at}"));
    let mut stack: Vec<Element> = Vec::new();
    let mut root: Option<Element> = None;
    let mut i = 0;
    let b = t.as_bytes();
    while i < b.len() {
        if b[i] != b'<' {
            let j = t[i..].find('<').map_or(b.len(), |k| i + k);
            if let Some(top) = stack.last_mut() {
                top.text.push_str(&unescape(&t[i..j]));
            }
            i = j;
            continue;
        }
        let rest = &t[i..];
        if rest.starts_with("<?") {
            i += rest.find("?>").ok_or_else(|| err("unclosed processing instruction", i))? + 2;
        } else if rest.starts_with("<!--") {
            i += rest.find("-->").ok_or_else(|| err("unclosed comment", i))? + 3;
        } else if rest.starts_with("<![CDATA[") {
            let j = rest.find("]]>").ok_or_else(|| err("unclosed CDATA", i))?;
            if let Some(top) = stack.last_mut() {
                top.text.push_str(&rest[9..j]);
            }
            i += j + 3;
        } else if rest.starts_with("<!") {
            // A DOCTYPE, with or without an internal subset.
            let mut depth = 0;
            let mut k = 0;
            for (n, c) in rest.char_indices() {
                match c {
                    '[' => depth += 1,
                    ']' => depth -= 1,
                    '>' if depth == 0 => {
                        k = n;
                        break;
                    }
                    _ => {}
                }
            }
            if k == 0 {
                return Err(err("unclosed declaration", i));
            }
            i += k + 1;
        } else if let Some(r) = rest.strip_prefix("</") {
            let j = r.find('>').ok_or_else(|| err("unclosed end tag", i))?;
            let name = local(r[..j].trim());
            let e = stack.pop().ok_or_else(|| err("end tag with no start", i))?;
            if e.name != name {
                return Err(err(&format!("</{name}> closes <{}>", e.name), i));
            }
            match stack.last_mut() {
                Some(p) => p.children.push(e),
                None => root = Some(e),
            }
            i += 2 + j + 1;
        } else {
            // A start tag: name, attributes (quoted with ' or "), and an optional '/'.
            let mut k = 1;
            let mut quote: Option<u8> = None;
            while i + k < b.len() {
                let c = b[i + k];
                match quote {
                    Some(q) if c == q => quote = None,
                    Some(_) => {}
                    None if c == b'"' || c == b'\'' => quote = Some(c),
                    None if c == b'>' => break,
                    None => {}
                }
                k += 1;
            }
            if i + k >= b.len() {
                return Err(err("unclosed start tag", i));
            }
            let inner = &t[i + 1..i + k];
            let empty = inner.ends_with('/');
            let inner = inner.trim_end_matches('/');
            let name_end = inner.find(|c: char| c.is_whitespace()).unwrap_or(inner.len());
            let mut e = Element { name: local(&inner[..name_end]), ..Default::default() };
            let mut a = inner[name_end..].trim_start();
            while let Some(eq) = a.find('=') {
                let key = a[..eq].trim();
                let v = a[eq + 1..].trim_start();
                let q = v.chars().next().ok_or_else(|| err("attribute without value", i))?;
                if q != '"' && q != '\'' {
                    return Err(err("unquoted attribute", i));
                }
                let end = v[1..].find(q).ok_or_else(|| err("unclosed attribute", i))?;
                e.attributes.push((local(key), unescape(&v[1..1 + end])));
                a = v[end + 2..].trim_start();
            }
            if empty {
                match stack.last_mut() {
                    Some(p) => p.children.push(e),
                    None => root = Some(e),
                }
            } else {
                stack.push(e);
            }
            i += k + 1;
        }
    }
    if !stack.is_empty() {
        return Err(XmlError(format!("XML: <{}> is not closed", stack.last().map_or("", |e| e.name.as_str()))));
    }
    root.ok_or(XmlError("XML: no root element".into()))
}

#[cfg(test)]
mod tests {
    #[test]
    fn elements_attributes_text() {
        let r = super::parse("<?xml version=\"1.0\"?>\n<!-- c --><a:root xmlns:a=\"u\"><b k='v &amp; w'>1 2<![CDATA[ <3>]]></b><c/></a:root>").unwrap();
        assert_eq!(r.name, "root");
        let b = r.child("b").unwrap();
        assert_eq!(b.attribute("k"), Some("v & w"));
        assert_eq!(b.text, "1 2 <3>");
        assert!(r.find("c").is_some());
    }
}
