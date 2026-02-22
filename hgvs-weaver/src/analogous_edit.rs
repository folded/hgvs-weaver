use crate::error::HgvsError;
use crate::structs::{AaEdit, NaEdit};
use crate::utils::{decompose_aa, normalize_aa};
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ResidueToken {
    Known(String),
    Unknown(i32),
    Any,      // Reference position
    Wildcard, // Matches any sequence of length 0 or more (for shorthand frameshifts)
}

impl ResidueToken {
    pub fn is_known(&self) -> bool {
        matches!(self, Self::Known(_))
    }

    pub fn unwrap_known(&self) -> String {
        match self {
            Self::Known(s) => s.clone(),
            Self::Unknown(p) => format!("Ref({})", p),
            Self::Any => "?".to_string(),
            Self::Wildcard => ".*".to_string(),
        }
    }

    pub fn normalized_symbol(&self) -> Option<String> {
        match self {
            Self::Known(s) => Some(normalize_aa(s)),
            Self::Unknown(_) => None,
            Self::Any => None,
            Self::Wildcard => None,
        }
    }
}

#[derive(Debug)]
pub struct SparseReference {
    pub data: BTreeMap<i32, ResidueToken>,
}

impl SparseReference {
    pub fn new() -> Self {
        Self {
            data: BTreeMap::new(),
        }
    }

    pub fn set(&mut self, pos: i32, symbol: String) -> Result<(), HgvsError> {
        let symbol_clean = normalize_aa(&symbol);
        if let Some(existing) = self.data.get(&pos) {
            match existing {
                ResidueToken::Known(ex) => {
                    let existing_clean = normalize_aa(ex);
                    if existing_clean != symbol_clean {
                        return Err(HgvsError::Other(format!(
                            "Consistency error at position {}: {} != {}",
                            pos, ex, symbol
                        )));
                    }
                }
                ResidueToken::Unknown(_) | ResidueToken::Any | ResidueToken::Wildcard => {}
            }
        }
        self.data.insert(pos, ResidueToken::Known(symbol));
        Ok(())
    }

    pub fn merge(&mut self, other: &Self) -> Result<(), HgvsError> {
        for (pos, token) in &other.data {
            if let ResidueToken::Known(symbol) = token {
                self.set(*pos, symbol.clone())?;
            }
        }
        Ok(())
    }

    pub fn project_range(&self, start: i32, end: i32) -> ProjectedSequence {
        let tokens = (start..=end)
            .map(|p| {
                self.data
                    .get(&p)
                    .cloned()
                    .unwrap_or(ResidueToken::Unknown(p))
            })
            .collect();
        ProjectedSequence(tokens)
    }

    pub fn get_range(&self) -> Vec<(i32, ResidueToken)> {
        if self.data.is_empty() {
            return Vec::new();
        }
        let min_p = *self.data.keys().next().unwrap();
        let max_p = *self.data.keys().last().unwrap();
        (min_p..=max_p)
            .map(|p| {
                (
                    p,
                    self.data
                        .get(&p)
                        .cloned()
                        .unwrap_or(ResidueToken::Unknown(p)),
                )
            })
            .collect()
    }
}

#[derive(Debug)]
pub struct ProjectedSequence(pub Vec<ResidueToken>);

impl ProjectedSequence {
    pub fn trim_at_stop(self) -> Self {
        let mut trimmed = Vec::new();
        for token in self.0 {
            if let Some(s) = token.normalized_symbol() {
                if s.contains('*') {
                    break;
                }
            }
            trimmed.push(token.clone());
        }
        ProjectedSequence(trimmed)
    }

    pub fn is_analogous_to(&self, other: &Self) -> bool {
        reconcile_projections(&self.0, &other.0)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

pub fn apply_aa_edit_to_sparse(
    edit: &AaEdit,
    start: i32,
    end: i32,
    sref: &SparseReference,
) -> ProjectedSequence {
    let mut res = Vec::new();

    match edit {
        AaEdit::Subst { alt, .. } => {
            if let Ok(aas) = decompose_aa(alt) {
                res.extend(aas.into_iter().map(|r| ResidueToken::Known(r.to_string())));
            }
        }
        AaEdit::Del { .. } => {}
        AaEdit::Ins { alt, .. } => {
            if let Ok(aas) = decompose_aa(alt) {
                res.extend(aas.into_iter().map(|r| ResidueToken::Known(r.to_string())));
            }
        }
        AaEdit::DelIns { alt, .. } => {
            if let Ok(aas) = decompose_aa(alt) {
                res.extend(aas.into_iter().map(|r| ResidueToken::Known(r.to_string())));
            }
        }
        AaEdit::RefAlt { alt, .. } => {
            if let Some(alt) = alt {
                if let Ok(aas) = decompose_aa(alt) {
                    res.extend(aas.into_iter().map(|r| ResidueToken::Known(r.to_string())));
                }
            }
        }
        AaEdit::Dup { .. } => {
            res.extend(sref.project_range(start, end).0);
            res.extend(sref.project_range(start, end).0);
        }
        AaEdit::Repeat { max, .. } => {
            for _ in 0..*max {
                res.extend(sref.project_range(start, end).0);
            }
        }
        AaEdit::Ext { alt, .. } => {
            if let Ok(aas) = decompose_aa(alt) {
                res.extend(aas.into_iter().map(|r| ResidueToken::Known(r.to_string())));
            }
        }
        AaEdit::Fs {
            alt, length, term, ..
        } => {
            let alt_tokens = if let Ok(aas) = decompose_aa(alt) {
                aas.into_iter()
                    .map(|r| ResidueToken::Known(r.to_string()))
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };

            let term_token = term
                .as_ref()
                .map(|t| ResidueToken::Known(normalize_aa(t).to_string()));

            let total_len = length.as_ref().and_then(|l| l.parse::<usize>().ok());

            if let Some(tlen) = total_len {
                // Explicit length: Alt + Gap + Term
                let alt_len = alt_tokens.len();
                let term_len = if term_token.is_some() { 1 } else { 0 };

                res.extend(alt_tokens);

                if tlen > alt_len + term_len {
                    let gap_len = tlen - alt_len - term_len;
                    for _ in 0..gap_len {
                        res.push(ResidueToken::Any);
                    }
                }

                if let Some(t) = term_token {
                    res.push(t);
                }
            } else {
                // Fallback: Just Alt
                res.extend(alt_tokens);
                // Mark as open-ended so it can match explicit full-form variants
                res.push(ResidueToken::Wildcard);
            }
        }
        AaEdit::Identity { .. } => {
            res.extend(sref.project_range(start, end).0);
        }
        AaEdit::Special { value, .. } if value == "=" => {
            res.extend(sref.project_range(start, end).0);
        }
        _ => {}
    }

    ProjectedSequence(res)
}

pub fn project_aa_variant(
    edit: &AaEdit,
    edit_start: i32,
    mut edit_end: i32,
    view_start: i32,
    view_end: i32,
    sref: &SparseReference,
) -> ProjectedSequence {
    let original_edit_end = edit_end;
    if let AaEdit::Repeat { .. } = edit {
        let unit = sref.project_range(edit_start, edit_end).0;
        let unit_len = unit.len() as i32;
        if unit_len > 0 {
            loop {
                let next_start = edit_end + 1;
                let next_end = edit_end + unit_len;
                let next_unit = sref.project_range(next_start, next_end).0;

                let mut matched = next_unit.len() == unit.len();
                for (a, b) in next_unit.iter().zip(unit.iter()) {
                    if a.normalized_symbol() != b.normalized_symbol()
                        || a.normalized_symbol().is_none()
                    {
                        matched = false;
                        break;
                    }
                }

                if matched {
                    edit_end = next_end;
                } else {
                    break;
                }
            }
        }
    }

    let mut res = Vec::new();
    let is_insert_point = matches!(edit, AaEdit::Ins { .. });

    let (prefix_end, suffix_start) = if is_insert_point {
        if edit_start == edit_end {
            (edit_start, edit_end + 1)
        } else {
            (edit_start, edit_end)
        }
    } else {
        (edit_start - 1, edit_end + 1)
    };

    if view_start <= prefix_end {
        res.extend(sref.project_range(view_start, prefix_end).0);
    }

    res.extend(apply_aa_edit_to_sparse(edit, edit_start, original_edit_end, sref).0);

    if suffix_start <= view_end {
        res.extend(sref.project_range(suffix_start, view_end).0);
    }

    ProjectedSequence(res)
}

pub fn apply_na_edit_to_sparse(
    edit: &NaEdit,
    start: i32,
    end: i32,
    sref: &SparseReference,
) -> ProjectedSequence {
    let mut res = Vec::new();
    match edit {
        NaEdit::RefAlt { alt, .. } => {
            if let Some(alt) = alt {
                for c in alt.chars() {
                    res.push(ResidueToken::Known(c.to_string()));
                }
            }
        }
        NaEdit::Del { .. } => {}
        NaEdit::Ins { alt, .. } => {
            if let Some(alt) = alt {
                for c in alt.chars() {
                    res.push(ResidueToken::Known(c.to_string()));
                }
            }
        }
        NaEdit::Dup { .. } => {
            res.extend(sref.project_range(start, end).0);
            res.extend(sref.project_range(start, end).0);
        }
        NaEdit::Inv { .. } => {
            let ref_tokens = sref.project_range(start, end).0;
            // Iterate in reverse and complement
            for t in ref_tokens.iter().rev() {
                match t {
                    ResidueToken::Known(s) => {
                        // normalized_symbol returns uppercase 1-letter code
                        // We need to complement it.
                        // utils::complement_dna_char works on char.
                        // Assuming s is single char?
                        let rc: String = s.chars().map(crate::utils::complement_dna_char).collect();
                        res.push(ResidueToken::Known(rc));
                    }
                    ResidueToken::Unknown(_) => res.push(t.clone()),
                    ResidueToken::Any => res.push(ResidueToken::Any),
                    ResidueToken::Wildcard => res.push(ResidueToken::Wildcard),
                }
            }
        }
        NaEdit::Repeat { max, .. } => {
            for _ in 0..*max {
                res.extend(sref.project_range(start, end).0);
            }
        }
        _ => {}
    }
    ProjectedSequence(res)
}

pub fn project_na_variant(
    edit: &NaEdit,
    edit_start: i32,
    mut edit_end: i32,
    view_start: i32,
    view_end: i32,
    sref: &SparseReference,
) -> ProjectedSequence {
    let original_edit_end = edit_end;
    if let NaEdit::Repeat { .. } = edit {
        let unit = sref.project_range(edit_start, edit_end).0;
        let unit_len = unit.len() as i32;
        if unit_len > 0 {
            loop {
                let next_start = edit_end + 1;
                let next_end = edit_end + unit_len;
                let next_unit = sref.project_range(next_start, next_end).0;

                let mut matched = next_unit.len() == unit.len();
                for (a, b) in next_unit.iter().zip(unit.iter()) {
                    if a.is_known()
                        && b.is_known()
                        && a.unwrap_known().to_uppercase() == b.unwrap_known().to_uppercase()
                    {
                        continue;
                    } else {
                        matched = false;
                        break;
                    }
                }

                if matched {
                    edit_end = next_end;
                } else {
                    break;
                }
            }
        }
    }

    let mut res = Vec::new();
    let is_insert_point = matches!(edit, NaEdit::Ins { .. });

    let (prefix_end, suffix_start) = if is_insert_point {
        if edit_start == edit_end {
            (edit_start, edit_end + 1)
        } else {
            (edit_start, edit_end)
        }
    } else {
        (edit_start - 1, edit_end + 1)
    };

    if view_start <= prefix_end {
        res.extend(sref.project_range(view_start, prefix_end).0);
    }

    res.extend(apply_na_edit_to_sparse(edit, edit_start, original_edit_end, sref).0);

    if suffix_start <= view_end {
        res.extend(sref.project_range(suffix_start, view_end).0);
    }

    ProjectedSequence(res)
}

struct UnificationEnv {
    aliases: HashMap<i32, ResidueToken>,
}

impl UnificationEnv {
    fn new() -> Self {
        Self {
            aliases: HashMap::new(),
        }
    }

    fn resolve(&self, t: &ResidueToken) -> ResidueToken {
        match t {
            ResidueToken::Known(_) | ResidueToken::Any | ResidueToken::Wildcard => t.clone(),
            ResidueToken::Unknown(p) => {
                let mut curr_p = *p;
                let mut visited = HashSet::new();
                visited.insert(curr_p);
                while let Some(next) = self.aliases.get(&curr_p) {
                    match next {
                        ResidueToken::Known(_) | ResidueToken::Any | ResidueToken::Wildcard => {
                            return next.clone()
                        }
                        ResidueToken::Unknown(next_p) => {
                            if visited.contains(next_p) {
                                break;
                            }
                            curr_p = *next_p;
                            visited.insert(curr_p);
                        }
                    }
                }
                ResidueToken::Unknown(curr_p)
            }
        }
    }

    fn unify(&mut self, t1: &ResidueToken, t2: &ResidueToken) -> bool {
        let r1 = self.resolve(t1);
        let r2 = self.resolve(t2);

        match (r1, r2) {
            (ResidueToken::Any, _) | (_, ResidueToken::Any) => true,
            (ResidueToken::Wildcard, _) | (_, ResidueToken::Wildcard) => true,
            (ResidueToken::Known(k1), ResidueToken::Known(k2)) => {
                let d1 = decompose_aa(&k1);
                let d2 = decompose_aa(&k2);
                if let (Ok(r1), Ok(r2)) = (d1, d2) {
                    r1 == r2
                } else {
                    normalize_aa(&k1) == normalize_aa(&k2)
                }
            }
            (ResidueToken::Unknown(p1), ResidueToken::Known(k2)) => {
                self.aliases.insert(p1, ResidueToken::Known(k2));
                true
            }
            (ResidueToken::Known(k1), ResidueToken::Unknown(p2)) => {
                self.aliases.insert(p2, ResidueToken::Known(k1));
                true
            }
            (ResidueToken::Unknown(p1), ResidueToken::Unknown(p2)) => {
                if p1 != p2 {
                    self.aliases.insert(p1, ResidueToken::Unknown(p2));
                }
                true
            }
        }
    }
}

pub fn reconcile_projections(v1: &[ResidueToken], v2: &[ResidueToken]) -> bool {
    let mut env = UnificationEnv::new();

    let mut i = 0;
    while i < v1.len() && i < v2.len() {
        let t1 = &v1[i];
        let t2 = &v2[i];

        if matches!(t1, ResidueToken::Wildcard) || matches!(t2, ResidueToken::Wildcard) {
            return true; // Matches the rest
        }

        if !env.unify(t1, t2) {
            return false;
        }
        i += 1;
    }

    // If one is longer, it must be because the other ended with a Wildcard
    // or they are same length.
    if i < v1.len() || i < v2.len() {
        // Zip only goes to the shortest. We check if the last processed was Wildcard.
        // Actually the loop handles it. If we reached here without a Wildcard,
        // then they must have same length.
        return v1.len() == v2.len();
    }

    // Second pass: verify consistency for the common part
    for j in 0..i {
        let r1 = env.resolve(&v1[j]);
        let r2 = env.resolve(&v2[j]);

        match (r1, r2) {
            (ResidueToken::Any, _) | (_, ResidueToken::Any) => {}
            (ResidueToken::Wildcard, _) | (_, ResidueToken::Wildcard) => break,
            (ResidueToken::Known(k1), ResidueToken::Known(k2)) => {
                let d1 = decompose_aa(&k1);
                let d2 = decompose_aa(&k2);
                if let (Ok(r1), Ok(r2)) = (d1, d2) {
                    if r1 != r2 {
                        return false;
                    }
                } else {
                    if normalize_aa(&k1) != normalize_aa(&k2) {
                        return false;
                    }
                }
            }
            (ResidueToken::Unknown(_), ResidueToken::Known(_))
            | (ResidueToken::Known(_), ResidueToken::Unknown(_)) => {
                return false;
            }
            (ResidueToken::Unknown(_), ResidueToken::Unknown(_)) => {}
        }
    }

    true
}
