use crate::error::HgvsError;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq)]
pub struct CigarOp {
    pub op: char,
    pub len: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cigar {
    pub ops: Vec<CigarOp>,
}

impl Cigar {
    pub fn ref_len(&self) -> i32 {
        self.ops.iter()
            .filter(|op| "=MXDN".contains(op.op))
            .map(|op| op.len)
            .sum()
    }

    pub fn tgt_len(&self) -> i32 {
        self.ops.iter()
            .filter(|op| "=MXI".contains(op.op))
            .map(|op| op.len)
            .sum()
    }
}

impl FromStr for Cigar {
    type Err = HgvsError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ops = Vec::new();
        let mut current_num = String::new();
        for c in s.chars() {
            if c.is_ascii_digit() {
                current_num.push(c);
            } else {
                let len = if current_num.is_empty() {
                    1
                } else {
                    current_num.parse::<i32>().map_err(|_| {
                        HgvsError::Other(format!("Invalid number in CIGAR string: {}", current_num))
                    })?
                };
                current_num.clear();
                ops.push(CigarOp { op: c, len });
            }
        }
        Ok(Cigar { ops })
    }
}

impl std::fmt::Display for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for op in &self.ops {
            if op.len != 1 {
                write!(f, "{}", op.len)?;
            }
            write!(f, "{}", op.op)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct CigarMapper {
    pub cigar: Cigar,
    pub ref_pos: Vec<i32>,
    pub tgt_pos: Vec<i32>,
}

impl CigarMapper {
    pub fn new(cigar_str: &str) -> Result<Self, HgvsError> {
        let cigar = Cigar::from_str(cigar_str)?;
        let mut ref_pos = Vec::with_capacity(cigar.ops.len() + 1);
        let mut tgt_pos = Vec::with_capacity(cigar.ops.len() + 1);

        let mut ref_cur = 0;
        let mut tgt_cur = 0;

        for op in &cigar.ops {
            ref_pos.push(ref_cur);
            tgt_pos.push(tgt_cur);

            if "=MXIN".contains(op.op) {
                ref_cur += op.len;
            }
            if "=MXD".contains(op.op) {
                tgt_cur += op.len;
            }
        }
        ref_pos.push(ref_cur);
        tgt_pos.push(tgt_cur);

        Ok(CigarMapper {
            cigar,
            ref_pos,
            tgt_pos,
        })
    }

    pub fn ref_len(&self) -> i32 {
        *self.ref_pos.last().unwrap()
    }

    pub fn tgt_len(&self) -> i32 {
        *self.tgt_pos.last().unwrap()
    }

    pub fn map_ref_to_tgt(&self, pos: i32, end_strategy: &str, strict_bounds: bool) -> Result<(i32, i32, char), HgvsError> {
        self.map_internal(&self.ref_pos, &self.tgt_pos, pos, end_strategy, strict_bounds)
    }

    pub fn map_tgt_to_ref(&self, pos: i32, end_strategy: &str, strict_bounds: bool) -> Result<(i32, i32, char), HgvsError> {
        self.map_internal(&self.tgt_pos, &self.ref_pos, pos, end_strategy, strict_bounds)
    }

    fn map_internal(
        &self,
        from_pos: &[i32],
        to_pos: &[i32],
        pos: i32,
        end_strategy: &str,
        strict_bounds: bool,
    ) -> Result<(i32, i32, char), HgvsError> {
        let last_pos = *from_pos.last().unwrap();
        if strict_bounds && (pos < 0 || pos > last_pos) {
            return Err(HgvsError::Other("Position is beyond the bounds of sequence".to_string()));
        }

        let mut pos_i = 0;
        for i in 0..self.cigar.ops.len() {
            if pos < from_pos[i+1] {
                pos_i = i;
                break;
            }
            pos_i = i;
        }

        let op_ref = &self.cigar.ops[pos_i];
        let op = op_ref.op;
        match op {
            '=' | 'M' | 'X' => {
                let mapped_pos = to_pos[pos_i] + (pos - from_pos[pos_i]);
                Ok((mapped_pos, 0, op))
            }
            'D' | 'I' => {
                let mut mapped_pos = to_pos[pos_i];
                if end_strategy == "start" {
                    mapped_pos -= 1;
                }
                Ok((mapped_pos, 0, op))
            }
            'N' => {
                if pos - from_pos[pos_i] < from_pos[pos_i + 1] - pos {
                    let mapped_pos = to_pos[pos_i] - 1;
                    let mapped_pos_offset = pos - from_pos[pos_i] + 1;
                    Ok((mapped_pos, mapped_pos_offset, op))
                } else {
                    let mapped_pos = to_pos[pos_i];
                    let mapped_pos_offset = -(from_pos[pos_i + 1] - pos);
                    Ok((mapped_pos, mapped_pos_offset, op))
                }
            }
            _ => Err(HgvsError::Other(format!("Unsupported CIGAR op: {}", op))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigarmapper() {
        let cigar_str = "3=2N=X=3N=I=D=";
        let cm = CigarMapper::new(cigar_str).unwrap();

        assert_eq!(cm.ref_len(), 15);
        assert_eq!(cm.tgt_len(), 10);

        // ref to tgt
        assert_eq!(cm.map_ref_to_tgt(0, "start", true).unwrap(), (0, 0, '='));
        assert_eq!(cm.map_ref_to_tgt(3, "start", true).unwrap(), (2, 1, 'N'));
        assert_eq!(cm.map_ref_to_tgt(4, "start", true).unwrap(), (3, -1, 'N'));
        assert_eq!(cm.map_ref_to_tgt(5, "start", true).unwrap(), (3, 0, '='));
        assert_eq!(cm.map_ref_to_tgt(6, "start", true).unwrap(), (4, 0, 'X'));
        assert_eq!(cm.map_ref_to_tgt(12, "start", true).unwrap(), (6, 0, 'I'));
        assert_eq!(cm.map_ref_to_tgt(12, "end", true).unwrap(), (7, 0, 'I'));
        assert_eq!(cm.map_ref_to_tgt(14, "start", true).unwrap(), (9, 0, '='));

        // tgt to ref
        assert_eq!(cm.map_tgt_to_ref(0, "start", true).unwrap(), (0, 0, '='));
        assert_eq!(cm.map_tgt_to_ref(3, "start", true).unwrap(), (5, 0, '='));
        assert_eq!(cm.map_tgt_to_ref(4, "start", true).unwrap(), (6, 0, 'X'));
        assert_eq!(cm.map_tgt_to_ref(8, "start", true).unwrap(), (13, 0, 'D'));
        assert_eq!(cm.map_tgt_to_ref(8, "end", true).unwrap(), (14, 0, 'D'));
    }

    #[test]
    fn test_cigarmapper_strict_bounds() {
        let cigar_str = "3=2N=X=3N=I=D=";
        let cm = CigarMapper::new(cigar_str).unwrap();

        assert!(cm.map_ref_to_tgt(-1, "start", true).is_err());
        assert!(cm.map_ref_to_tgt(16, "start", true).is_err());

        assert_eq!(cm.map_ref_to_tgt(0, "start", true).unwrap(), (0, 0, '='));
        assert_eq!(cm.map_ref_to_tgt(-1, "start", false).unwrap(), (-1, 0, '='));
        assert_eq!(cm.map_ref_to_tgt(15, "start", true).unwrap(), (10, 0, '='));
        assert_eq!(cm.map_ref_to_tgt(14, "start", false).unwrap(), (9, 0, '='));
    }

    #[test]
    fn test_cigar_parsing() {
        let cigar_str = "3=2N=X=3N=I=D=";
        let cigar = Cigar::from_str(cigar_str).unwrap();
        assert_eq!(cigar.ops.len(), 11);
        assert_eq!(cigar.ops[0], CigarOp { op: '=', len: 3 });
        assert_eq!(cigar.to_string(), cigar_str);
    }
}
