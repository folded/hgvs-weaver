use super::Rule;
use crate::error::HgvsError;
use crate::structs::*;
use pest::iterators::Pair;

pub fn parse_g_posedit(pair: Pair<Rule>) -> Result<PosEdit<SimpleInterval, NaEdit>, HgvsError> {
    let mut inner = pair.into_inner();
    let pos = parse_simple_interval(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?,
    )?;
    let edit = parse_na_edit(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing edit".into()))?,
    )?;
    Ok(PosEdit {
        pos: Some(pos),
        edit,
        uncertain: false,
        predicted: false,
    })
}

pub fn parse_c_posedit(pair: Pair<Rule>) -> Result<PosEdit<BaseOffsetInterval, NaEdit>, HgvsError> {
    let mut inner = pair.into_inner();
    let pos = parse_base_offset_interval(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?,
    )?;
    let edit = parse_na_edit(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing edit".into()))?,
    )?;
    Ok(PosEdit {
        pos: Some(pos),
        edit,
        uncertain: false,
        predicted: false,
    })
}

pub fn parse_p_posedit(pair: Pair<Rule>) -> Result<PosEdit<AaInterval, AaEdit>, HgvsError> {
    let mut predicted = false;
    let s = pair.as_str();
    if s.starts_with('(') && s.ends_with(')') {
        predicted = true;
    }

    let mut inner = pair.into_inner();
    let inner_pair = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty p_posedit".into()))?;
    if inner_pair.as_rule() == Rule::p_posedit_special {
        let special = inner_pair.as_str();
        let edit = AaEdit::Special {
            value: special.replace(['(', ')'], ""),
            uncertain: false,
        };
        return Ok(PosEdit {
            pos: None,
            edit,
            uncertain: false,
            predicted,
        });
    }

    let mut pos = None;
    let mut edit = AaEdit::None;

    if inner_pair.as_rule() == Rule::p_interval {
        pos = Some(parse_aa_interval(inner_pair)?);
        if let Some(e) = inner.next() {
            edit = parse_pro_edit(e)?;
        }
    } else if inner_pair.as_rule() == Rule::pro_edit {
        edit = parse_pro_edit(inner_pair)?;
    }

    Ok(PosEdit {
        pos,
        edit,
        uncertain: false,
        predicted,
    })
}

pub fn parse_simple_interval(pair: Pair<Rule>) -> Result<SimpleInterval, HgvsError> {
    let s = pair.as_str();
    let mut uncertain = false;
    if s.starts_with('(') && s.ends_with(')') && !s.contains('_') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty interval".into()))?;
    match p.as_rule() {
        Rule::def_g_interval | Rule::def_m_interval => {
            let mut parts = p.into_inner();
            let start = parse_simple_pos(
                parts
                    .next()
                    .ok_or_else(|| HgvsError::PestError("Missing start position".into()))?,
            )?;
            let end = parts.next().map(parse_simple_pos).transpose()?;
            Ok(SimpleInterval {
                start,
                end,
                uncertain,
            })
        }
        Rule::uncertain_g_interval => {
            let mut start = None;
            let mut end = None;
            for sub in p.into_inner() {
                if sub.as_rule() == Rule::def_g_interval {
                    let mut parts = sub.into_inner();
                    let s =
                        parse_simple_pos(parts.next().ok_or_else(|| {
                            HgvsError::PestError("Missing start position".into())
                        })?)?;
                    let e = parts.next().map(parse_simple_pos).transpose()?;

                    let pos = SimplePosition {
                        base: s.base,
                        end: e.map(|x| x.base),
                        uncertain: true,
                    };

                    if start.is_none() {
                        start = Some(pos);
                    } else {
                        end = Some(pos);
                    }
                }
            }
            Ok(SimpleInterval {
                start: start.ok_or_else(|| {
                    HgvsError::PestError("Missing start position in uncertain interval".into())
                })?,
                end,
                uncertain: false,
            })
        }
        _ => Err(HgvsError::PestError(format!(
            "Unexpected interval rule: {:?}",
            p.as_rule()
        ))),
    }
}

pub fn parse_simple_pos(pair: Pair<Rule>) -> Result<SimplePosition, HgvsError> {
    let s = pair.as_str();
    if s == "?" {
        return Ok(SimplePosition {
            base: HgvsGenomicPos(0),
            end: None,
            uncertain: true,
        });
    }
    let hgvs_base = s
        .parse::<i32>()
        .map_err(|_| HgvsError::PestError("Invalid position".into()))?;
    Ok(SimplePosition {
        base: HgvsGenomicPos(hgvs_base),
        end: None,
        uncertain: false,
    })
}

pub fn parse_n_posedit(pair: Pair<Rule>) -> Result<PosEdit<BaseOffsetInterval, NaEdit>, HgvsError> {
    let mut inner = pair.into_inner();
    let pos = parse_base_offset_interval_n(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?,
    )?;
    let edit = parse_na_edit(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing edit".into()))?,
    )?;
    Ok(PosEdit {
        pos: Some(pos),
        edit,
        uncertain: false,
        predicted: false,
    })
}

pub fn parse_base_offset_interval_n(pair: Pair<Rule>) -> Result<BaseOffsetInterval, HgvsError> {
    let mut uncertain = false;
    let s = pair.as_str();
    if s.starts_with('(') && s.ends_with(')') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty base offset interval".into()))?;
    let mut p_inner = p.into_inner();

    let start = parse_base_offset_pos_with_default(
        p_inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing start position".into()))?,
        Anchor::TranscriptStart,
    )?;
    let end = p_inner
        .next()
        .map(|p| parse_base_offset_pos_with_default(p, Anchor::TranscriptStart))
        .transpose()?;
    Ok(BaseOffsetInterval {
        start,
        end,
        uncertain,
    })
}

pub fn parse_base_offset_interval(pair: Pair<Rule>) -> Result<BaseOffsetInterval, HgvsError> {
    let mut uncertain = false;
    let s = pair.as_str();
    if s.starts_with('(') && s.ends_with(')') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty base offset interval".into()))?;
    let mut p_inner = p.into_inner();

    let start = parse_base_offset_pos(
        p_inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing start position".into()))?,
    )?;
    let end = p_inner.next().map(parse_base_offset_pos).transpose()?;
    Ok(BaseOffsetInterval {
        start,
        end,
        uncertain,
    })
}

pub fn parse_base_offset_pos(pair: Pair<Rule>) -> Result<BaseOffsetPosition, HgvsError> {
    parse_base_offset_pos_with_default(pair, Anchor::CdsStart)
}

pub fn parse_base_offset_pos_with_default(
    pair: Pair<Rule>,
    default_anchor: Anchor,
) -> Result<BaseOffsetPosition, HgvsError> {
    let mut anchor = default_anchor;
    if pair.as_str().starts_with('*') {
        anchor = Anchor::CdsEnd;
    }

    let mut hgvs_base = 0;
    let mut hgvs_offset: Option<IntronicOffset> = None;

    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::num | Rule::base => {
                hgvs_base = p.as_str().parse().unwrap_or(0);
            }
            Rule::offset => {
                let off_str = p.as_str();
                if !off_str.is_empty() {
                    hgvs_offset = Some(IntronicOffset(
                        off_str.replace('+', "").parse().unwrap_or(0),
                    ));
                }
            }
            _ => {}
        }
    }

    Ok(BaseOffsetPosition {
        base: HgvsTranscriptPos(hgvs_base),
        offset: hgvs_offset,
        anchor,
        uncertain: false,
    })
}

pub fn parse_aa_interval(pair: Pair<Rule>) -> Result<AaInterval, HgvsError> {
    let s = pair.as_str();
    let mut uncertain = false;
    if s.starts_with('(') && s.ends_with(')') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty AA interval".into()))?;
    let mut p_inner = p.into_inner();

    let start = parse_aa_pos(
        p_inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing start AA position".into()))?,
    )?;
    let end = p_inner.next().map(parse_aa_pos).transpose()?;
    Ok(AaInterval {
        start,
        end,
        uncertain,
    })
}

pub fn parse_aa_pos(pair: Pair<Rule>) -> Result<AAPosition, HgvsError> {
    let mut aa = String::new();
    let mut pos = 0;

    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aa13 | Rule::term13 | Rule::aa3 | Rule::aa1 | Rule::term3 | Rule::term1 => {
                aa = p.as_str().to_string();
            }
            Rule::num => {
                pos = p.as_str().parse().unwrap_or(0);
            }
            _ => {}
        }
    }

    Ok(AAPosition {
        base: HgvsProteinPos(pos),
        aa,
        uncertain: false,
    })
}

pub fn parse_na_edit(pair: Pair<Rule>) -> Result<NaEdit, HgvsError> {
    let mut inner = pair.into_inner();
    let inner_feat = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty na_edit".into()))?;
    match inner_feat.as_rule() {
        Rule::dna_subst | Rule::rna_subst => {
            let mut parts = inner_feat.into_inner();
            let ref_ = parts.next().map(|p: Pair<Rule>| p.as_str().to_string());
            let alt = parts.next().map(|p: Pair<Rule>| p.as_str().to_string());
            Ok(NaEdit::RefAlt {
                ref_,
                alt,
                uncertain: false,
            })
        }
        Rule::dna_del | Rule::rna_del => {
            let ref_ = inner_feat
                .into_inner()
                .next()
                .map(|p: Pair<Rule>| p.as_str().to_string());
            Ok(NaEdit::Del {
                ref_,
                uncertain: false,
            })
        }
        Rule::dna_ins | Rule::rna_ins => {
            let alt = inner_feat
                .into_inner()
                .next()
                .map(|p: Pair<Rule>| p.as_str().to_string());
            Ok(NaEdit::Ins {
                alt,
                uncertain: false,
            })
        }
        Rule::dna_delins | Rule::rna_delins => {
            let mut parts = inner_feat.into_inner();
            let first = parts
                .next()
                .map(|p: Pair<Rule>| p.as_str().to_string())
                .unwrap_or_default();
            let second = parts.next().map(|p: Pair<Rule>| p.as_str().to_string());
            if second.is_none() {
                Ok(NaEdit::RefAlt {
                    ref_: Some("".to_string()),
                    alt: Some(first),
                    uncertain: false,
                })
            } else {
                Ok(NaEdit::RefAlt {
                    ref_: Some(first),
                    alt: second,
                    uncertain: false,
                })
            }
        }
        Rule::dna_dup | Rule::rna_dup => {
            let ref_ = inner_feat
                .into_inner()
                .next()
                .map(|p: Pair<Rule>| p.as_str().to_string());
            Ok(NaEdit::Dup {
                ref_,
                uncertain: false,
            })
        }
        Rule::dna_inv | Rule::rna_inv => {
            let ref_ = inner_feat
                .into_inner()
                .next()
                .map(|p: Pair<Rule>| p.as_str().to_string());
            Ok(NaEdit::Inv {
                ref_,
                uncertain: false,
            })
        }
        Rule::dna_ident | Rule::rna_ident => {
            let mut inner = inner_feat.into_inner();
            let mut ref_ = None;
            if let Some(p) = inner.next() {
                if p.as_rule() == Rule::dna || p.as_rule() == Rule::rna {
                    ref_ = Some(p.as_str().to_string());
                }
            }
            Ok(NaEdit::RefAlt {
                ref_: ref_.clone(),
                alt: ref_.clone(),
                uncertain: false,
            })
        }
        Rule::dna_repeat | Rule::rna_repeat => {
            let inner = inner_feat.into_inner();
            let mut ref_ = None;
            let mut first = 0;
            let mut second = None;

            for p in inner {
                match p.as_rule() {
                    Rule::dna | Rule::rna => ref_ = Some(p.as_str().to_string()),
                    Rule::num => {
                        if first == 0 {
                            first = p.as_str().parse().unwrap_or(0);
                        } else {
                            second = Some(p.as_str().parse().unwrap_or(0));
                        }
                    }
                    _ => {}
                }
            }
            let max = second.unwrap_or(first);
            Ok(NaEdit::Repeat {
                ref_,
                min: first,
                max,
                uncertain: false,
            })
        }
        Rule::dna_con | Rule::rna_con => {
            // Placeholder/Generic for now as struct support is minimal
            Ok(NaEdit::None)
        }
        Rule::dna_copy => {
            let mut inner = inner_feat.into_inner();
            let copy = inner
                .next()
                .map(|p| p.as_str().parse().unwrap_or(0))
                .unwrap_or(0);
            Ok(NaEdit::NACopy {
                copy,
                uncertain: false,
            })
        }
        _ => Ok(NaEdit::None),
    }
}

pub fn parse_pro_edit(pair: Pair<Rule>) -> Result<AaEdit, HgvsError> {
    let inner = pair
        .into_inner()
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty pro_edit".into()))?;
    match inner.as_rule() {
        Rule::pro_ident => Ok(AaEdit::Identity { uncertain: false }),
        Rule::pro_subst => Ok(AaEdit::Subst {
            ref_: "".into(),
            alt: inner.as_str().to_string(),
            uncertain: false,
        }),
        Rule::pro_del => Ok(AaEdit::Del {
            ref_: "".into(),
            uncertain: false,
        }),
        Rule::pro_ins => Ok(AaEdit::Ins {
            alt: inner
                .into_inner()
                .next()
                .map(|p| p.as_str().to_string())
                .unwrap_or_default(),
            uncertain: false,
        }),
        Rule::pro_dup => Ok(AaEdit::Dup {
            ref_: None,
            uncertain: false,
        }),
        Rule::pro_delins => Ok(AaEdit::DelIns {
            ref_: "".into(),
            alt: inner
                .into_inner()
                .next()
                .map(|p| p.as_str().to_string())
                .unwrap_or_default(),
            uncertain: false,
        }),
        Rule::pro_fs => {
            let mut alt = String::new();
            let mut term = None;
            let mut length = None;
            for p in inner.into_inner() {
                match p.as_rule() {
                    Rule::aat13 => alt = p.as_str().to_string(),
                    Rule::fs => {
                        let mut fs_inner = p.into_inner();
                        if let Some(aa_fs) = fs_inner.next() {
                            let mut p_inner = aa_fs.into_inner();
                            term = p_inner.next().map(|t| t.as_str().to_string());
                            length = p_inner.next().map(|l| l.as_str().to_string());
                        }
                    }
                    _ => {}
                }
            }
            Ok(AaEdit::Fs {
                ref_: "".into(),
                alt,
                term,
                length,
                uncertain: false,
            })
        }
        Rule::pro_ext => {
            let ref_ = String::new();
            let mut alt = String::new();
            let mut aaterm = None;
            let mut length = None;
            for p in inner.into_inner() {
                match p.as_rule() {
                    Rule::aat13
                    | Rule::aat3
                    | Rule::aat1
                    | Rule::aa3
                    | Rule::aa1
                    | Rule::term3
                    | Rule::term1 => alt = p.as_str().to_string(),
                    Rule::ext => {
                        let mut ext_inner = p.into_inner();
                        if let Some(aa_ext) = ext_inner.next() {
                            let mut p_inner = aa_ext.into_inner();
                            aaterm = p_inner.next().map(|t| t.as_str().to_string());
                            length = p_inner.next().map(|l| l.as_str().to_string());
                        }
                    }
                    _ => {}
                }
            }
            Ok(AaEdit::Ext {
                ref_,
                alt,
                aaterm,
                length,
                uncertain: false,
            })
        }
        Rule::pro_repeat => {
            let inner = inner.into_inner();
            let mut ref_ = None;
            let mut first = 0;
            let mut second = None;

            for p in inner {
                match p.as_rule() {
                    Rule::aat13_seq => ref_ = Some(p.as_str().to_string()),
                    Rule::num => {
                        if first == 0 {
                            first = p.as_str().parse().unwrap_or(0);
                        } else {
                            second = Some(p.as_str().parse().unwrap_or(0));
                        }
                    }
                    _ => {}
                }
            }
            let max = second.unwrap_or(first);
            Ok(AaEdit::Repeat {
                ref_,
                min: first,
                max,
                uncertain: false,
            })
        }
        _ => Ok(AaEdit::None),
    }
}

#[cfg(test)]
mod tests {
    use crate::parse_hgvs_variant;

    #[test]
    fn test_unsupported_circular_no_panic() {
        // 'o.' is mentioned in docs but not supported in SequenceVariant.
        // Ensure it doesn't panic, just returns an error.
        let result = parse_hgvs_variant("NC_000001.11:o.123A>G");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_repeats() {
        let v_c = parse_hgvs_variant("NM_001291285.3:c.7035TGGAAC[3]").unwrap();
        match v_c {
            crate::coords::SequenceVariant::Coding(v) => match v.posedit.edit {
                crate::edits::NaEdit::Repeat { ref_, min, max, .. } => {
                    assert_eq!(ref_, Some("TGGAAC".to_string()));
                    assert_eq!(min, 3);
                    assert_eq!(max, 3);
                }
                _ => panic!("Expected Repeat edit"),
            },
            _ => panic!("Expected Coding variant"),
        }

        let v_p = parse_hgvs_variant("NP_001278214.1:p.2346GT[3]").unwrap();
        match v_p {
            crate::coords::SequenceVariant::Protein(v) => match v.posedit.edit {
                crate::edits::AaEdit::Repeat { ref_, min, max, .. } => {
                    assert_eq!(ref_, Some("GT".to_string()));
                    assert_eq!(min, 3);
                    assert_eq!(max, 3);
                }
                _ => panic!("Expected Repeat edit for protein"),
            },
            _ => panic!("Expected Protein variant"),
        }
    }

    #[test]
    fn test_parse_extension() {
        let v_p = parse_hgvs_variant("NP_001116078.1:p.Ter312Argext*5").unwrap();
        match v_p {
            crate::coords::SequenceVariant::Protein(v) => match v.posedit.edit {
                crate::edits::AaEdit::Ext {
                    ref_,
                    alt,
                    aaterm,
                    length,
                    ..
                } => {
                    assert_eq!(ref_, "");
                    assert_eq!(alt, "Arg");
                    assert_eq!(aaterm, Some("*".to_string()));
                    assert_eq!(length, Some("5".to_string()));
                }
                _ => panic!("Expected Extension edit"),
            },
            _ => panic!("Expected Protein variant"),
        }
    }
}
