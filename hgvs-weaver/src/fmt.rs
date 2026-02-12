use std::fmt;
use crate::structs::*;
use crate::edits::{AaEdit, NaEdit};

impl fmt::Display for SequenceVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SequenceVariant::Genomic(v) => write!(f, "{}", v),
            SequenceVariant::Coding(v) => write!(f, "{}", v),
            SequenceVariant::Protein(v) => write!(f, "{}", v),
            SequenceVariant::Mitochondrial(v) => write!(f, "{}", v),
            SequenceVariant::NonCoding(v) => write!(f, "{}", v),
            SequenceVariant::Rna(v) => write!(f, "{}", v),
        }
    }
}

impl fmt::Display for GVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}:g.{}", self.ac, self.gene.as_ref().map(|g| format!("({})", g)).unwrap_or_default(), self.posedit)
    }
}

impl fmt::Display for CVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}:c.{}", self.ac, self.gene.as_ref().map(|g| format!("({})", g)).unwrap_or_default(), self.posedit)
    }
}

impl AaEdit {
    pub fn format_simple(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AaEdit::Del { .. } => write!(f, "del"),
            AaEdit::Ins { alt, .. } => write!(f, "ins{}", alt),
            AaEdit::Dup { .. } => write!(f, "dup"),
            AaEdit::DelIns { alt, .. } => write!(f, "delins{}", alt),
            _ => write!(f, "{}", self),
        }
    }
}

impl<I: fmt::Display, E: fmt::Display> PosEdit<I, E> {
    pub fn format_simple(&self, f: &mut fmt::Formatter) -> fmt::Result
    where E: SimpleFormatter {
        if self.predicted { write!(f, "(")?; }
        if let Some(pos) = &self.pos {
            write!(f, "{}", pos)?;
        }
        self.edit.fmt_simple(f)?;
        if self.uncertain { write!(f, "?")?; }
        if self.predicted { write!(f, ")")?; }
        Ok(())
    }
}

pub trait SimpleFormatter {
    fn fmt_simple(&self, f: &mut fmt::Formatter) -> fmt::Result;
}

impl SimpleFormatter for AaEdit {
    fn fmt_simple(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.format_simple(f)
    }
}

impl SimpleFormatter for NaEdit {
    fn fmt_simple(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self)
    }
}

impl fmt::Display for PVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}:p.", self.ac, self.gene.as_ref().map(|g| format!("({})", g)).unwrap_or_default())?;
        self.posedit.format_simple(f)
    }
}

impl fmt::Display for MVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}:m.{}", self.ac, self.gene.as_ref().map(|g| format!("({})", g)).unwrap_or_default(), self.posedit)
    }
}

impl fmt::Display for NVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}:n.{}", self.ac, self.gene.as_ref().map(|g| format!("({})", g)).unwrap_or_default(), self.posedit)
    }
}

impl fmt::Display for RVariant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}:r.{}", self.ac, self.gene.as_ref().map(|g| format!("({})", g)).unwrap_or_default(), self.posedit)
    }
}

impl<I: fmt::Display, E: fmt::Display> fmt::Display for PosEdit<I, E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.predicted { write!(f, "(")?; }
        if let Some(pos) = &self.pos {
            write!(f, "{}", pos)?;
        }
        write!(f, "{}", self.edit)?;
        if self.uncertain { write!(f, "?")?; }
        if self.predicted { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for SimpleInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.uncertain { write!(f, "(")?; }
        write!(f, "{}", self.start)?;
        if let Some(end) = &self.end {
            write!(f, "_{}", end)?;
        }
        if self.uncertain { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for SimplePosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.uncertain { write!(f, "(")?; }
        write!(f, "{}", self.base.0)?;
        if let Some(end) = self.end {
            write!(f, "_{}", end.0)?;
        }
        if self.uncertain { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for BaseOffsetInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.uncertain { write!(f, "(")?; }
        write!(f, "{}", self.start)?;
        if let Some(end) = &self.end {
            write!(f, "_{}", end)?;
        }
        if self.uncertain { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for BaseOffsetPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.uncertain { write!(f, "(")?; }
        if self.anchor == Anchor::CdsEnd { write!(f, "*")?; }
        write!(f, "{}", self.base.0)?;
        if let Some(offset) = &self.offset {
            if offset.0 >= 0 { write!(f, "+")?; }
            write!(f, "{}", offset.0)?;
        }
        if self.uncertain { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for AaInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.uncertain { write!(f, "(")?; }
        write!(f, "{}", self.start)?;
        if let Some(end) = &self.end {
            write!(f, "_{}", end)?;
        }
        if self.uncertain { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for AAPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.uncertain { write!(f, "(")?; }
        write!(f, "{}{}", self.aa, self.base.0)?;
        if self.uncertain { write!(f, ")")?; }
        Ok(())
    }
}

impl fmt::Display for NaEdit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NaEdit::RefAlt { ref_, alt, uncertain } => {
                match (ref_, alt) {
                    (Some(r), Some(a)) => {
                        if r == a { write!(f, "=")?; }
                        else if r.len() == 1 && a.len() == 1 { write!(f, "{}>{}", r, a)?; }
                        else { write!(f, "del{}ins{}", r, a)?; }
                    },
                    (Some(r), None) => { write!(f, "del{}", r)?; },
                    (None, Some(a)) => { write!(f, "ins{}", a)?; },
                    (None, None) => { write!(f, "=")?; }
                }
                if *uncertain { write!(f, "?")?; }
                Ok(())
            }
            NaEdit::Del { ref_, uncertain } => {
                write!(f, "del")?;
                if let Some(r) = ref_ { write!(f, "{}", r)?; }
                if *uncertain { write!(f, "?")?; }
                Ok(())
            }
            NaEdit::Ins { alt, uncertain } => {
                write!(f, "ins")?;
                if let Some(a) = alt { write!(f, "{}", a)?; }
                if *uncertain { write!(f, "?")?; }
                Ok(())
            }
            NaEdit::Dup { ref_, uncertain } => {
                write!(f, "dup")?;
                if let Some(r) = ref_ { write!(f, "{}", r)?; }
                if *uncertain { write!(f, "?")?; }
                Ok(())
            }
            NaEdit::Inv { ref_, uncertain } => {
                write!(f, "inv")?;
                if let Some(r) = ref_ { write!(f, "{}", r)?; }
                if *uncertain { write!(f, "?")?; }
                Ok(())
            }
            NaEdit::Repeat { ref_, min, max, .. } => {
                if let Some(r) = ref_ { write!(f, "{}", r)?; }
                if min == max { write!(f, "[{}]", min)?; }
                else { write!(f, "[{}_{}]", min, max)?; }
                Ok(())
            }
            NaEdit::None => write!(f, "="),
            _ => write!(f, "unknown_na_edit"),
        }
    }
}

impl fmt::Display for AaEdit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AaEdit::Subst { alt, .. } => {
                let a3 = if alt == "*" { "Ter" } else { alt };
                write!(f, "{}", a3)
            }
            AaEdit::Del { ref_, .. } => {
                write!(f, "del")?;
                if !ref_.is_empty() { write!(f, "{}", ref_)?; }
                Ok(())
            }
            AaEdit::Ins { alt, .. } => {
                write!(f, "ins{}", alt)
            }
            AaEdit::DelIns { ref_, alt, .. } => {
                if !ref_.is_empty() { write!(f, "del{}ins{}", ref_, alt) }
                else { write!(f, "delins{}", alt) }
            }
            AaEdit::Dup { .. } => write!(f, "dup"),
            AaEdit::Identity { .. } => write!(f, "="),
            AaEdit::RefAlt { ref_, alt, .. } => {
                if let (Some(r), Some(a)) = (ref_, alt) {
                    if r == a { write!(f, "=") }
                    else if r.is_empty() { write!(f, "ins{}", a) }
                    else if a.is_empty() { write!(f, "del{}", r) }
                    else { write!(f, "del{}ins{}", r, a) }
                } else if let Some(a) = alt {
                    let a3 = if a == "*" { "Ter" } else { a };
                    write!(f, "delins{}", a3)
                } else {
                    write!(f, "del")
                }
            }
            AaEdit::Fs { alt, term, length, .. } => {
                let a3 = if alt == "*" { "Ter" } else { alt };
                if a3.is_empty() {
                    write!(f, "fs")?;
                } else {
                    write!(f, "{}fs{}", a3, term.as_deref().unwrap_or(""))?;
                }
                if let Some(l) = length { write!(f, "{}", l)?; }
                Ok(())
            }
            AaEdit::Ext { alt, length, .. } => {
                let a3 = if alt == "*" { "Ter" } else { alt };
                write!(f, "ext*{}", a3)?;
                if let Some(l) = length { write!(f, "{}", l)?; }
                Ok(())
            }
            AaEdit::Special { value, .. } => write!(f, "{}", value),
            _ => write!(f, "unknown_aa_edit"),
        }
    }
}

pub fn aa1_to_aa3(c: char) -> &'static str {
    match c {
        'A' => "Ala", 'C' => "Cys", 'D' => "Asp", 'E' => "Glu", 'F' => "Phe",
        'G' => "Gly", 'H' => "His", 'I' => "Ile", 'K' => "Lys", 'L' => "Leu",
        'M' => "Met", 'N' => "Asn", 'P' => "Pro", 'Q' => "Gln", 'R' => "Arg",
        'S' => "Ser", 'T' => "Thr", 'V' => "Val", 'W' => "Trp", 'Y' => "Tyr",
        '*' => "Ter", _ => "Xaa",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::HgvsTranscriptPos;

    #[test]
    fn test_non_coding_with_asterisk_no_panic() {
        // Construct a non-coding variant with a CdsEnd anchor.
        // Documentation says n. shouldn't use '*', but we want to ensure no panic if it happens.
        let pos = BaseOffsetPosition {
            base: HgvsTranscriptPos(1),
            offset: None,
            anchor: Anchor::CdsEnd,
            uncertain: false,
        };
        let var = NVariant {
            ac: "NR_0001.1".to_string(),
            gene: None,
            posedit: PosEdit {
                pos: Some(BaseOffsetInterval { start: pos, end: None, uncertain: false }),
                edit: NaEdit::None,
                uncertain: false,
                predicted: false,
            },
        };
        let formatted = format!("{}", var);
        assert_eq!(formatted, "NR_0001.1:n.*1=");
    }
}
