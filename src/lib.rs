use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3_stub_gen::{define_stub_info_gatherer, derive::*};
use ::hgvs_weaver::{SequenceVariant, Variant as VariantTrait, VariantMapper, DataProvider, TranscriptSearch, Transcript, IdentifierKind, HgvsError};
use serde_json;

#[gen_stub_pyclass]
#[pyclass(name = "Variant")]
#[doc = "Represents a parsed HGVS variant.\n\nProvides access to the variant's accession, gene symbol, and coordinate type.\nVariants can be formatted back to HGVS strings or converted to JSON/dict representations."]
#[derive(Clone)]
pub struct PyVariant {
    pub inner: SequenceVariant,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyVariant {
    #[getter]
    #[doc = "The primary accession of the variant (e.g., 'NM_000051.3')."]
    fn ac(&self) -> String {
        self.inner.ac().to_string()
    }

    #[getter]
    #[doc = "The gene symbol associated with the variant, if available."]
    fn gene(&self) -> Option<String> {
        self.inner.gene().map(|s| s.to_string())
    }

    #[getter]
    #[doc = "The coordinate type of the variant ('g', 'c', 'p', etc.)."]
    fn coordinate_type(&self) -> String {
        self.inner.coordinate_type().to_string()
    }

    #[doc = "Formats the variant back into a standard HGVS string."]
    fn format(&self) -> String {
        self.inner.to_string()
    }

    #[doc = "Returns a JSON string representation of the internal variant structure."]
    fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(&self.inner)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    #[doc = "Returns a dictionary representation of the internal variant structure."]
    fn to_dict(&self, py: Python) -> PyResult<Py<PyAny>> {
        let json_str = self.to_json()?;
        let json_mod = py.import("json")?;
        let dict = json_mod.call_method1("loads", (json_str,))?;
        Ok(dict.unbind())
    }

    fn __str__(&self) -> String {
        self.format()
    }

    fn __repr__(&self) -> String {
        format!("<weaver.Variant {}>", self.format())
    }

    #[doc = "Validates the variant's reference sequence against the provided DataProvider.\n\nReturns True if the reference sequence matches, False otherwise.\nMay raise ValueError if coordinates are out of bounds."]
    fn validate(&self, _py: Python, provider: Py<PyAny>) -> PyResult<bool> {
        let bridge = PyDataProviderBridge { provider };
        let result = match &self.inner {
            SequenceVariant::Genomic(v) => self.validate_genomic(v, &bridge),
            SequenceVariant::Coding(v) => self.validate_coding(v, &bridge),
            _ => Err(HgvsError::UnsupportedOperation("Validation not implemented for this variant type".into())),
        };

        match result {
            Ok(is_valid) => Ok(is_valid),
            Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }
}

impl PyVariant {
    fn validate_genomic(&self, v: &::hgvs_weaver::GVariant, bridge: &PyDataProviderBridge) -> Result<bool, HgvsError> {
        let pos = v.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start_0 = pos.start.base.to_index();
        let end_0 = pos.end.as_ref().map_or(start_0 + 1, |e| e.base.to_index() + 1);

        let ref_seq = bridge.get_seq(&v.ac, start_0.0, end_0.0, IdentifierKind::Genomic)?;

        match &v.posedit.edit {
            ::hgvs_weaver::edits::NaEdit::RefAlt { ref_: Some(r), .. } => {
                if r.is_empty() || r.chars().all(|c| c.is_ascii_digit()) { return Ok(true); }
                Ok(r == &ref_seq)
            }
            _ => Ok(true),
        }
    }

    fn validate_coding(&self, v: &::hgvs_weaver::CVariant, bridge: &PyDataProviderBridge) -> Result<bool, HgvsError> {
        let transcript = bridge.get_transcript(&v.ac, None)?;

        let pos = v.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;

        let ref_seq = bridge.get_seq(&v.ac, 0, -1, IdentifierKind::Transcript)?;

        if pos.start.offset.is_some() || pos.end.as_ref().and_then(|e| e.offset).is_some() {
            return Ok(true);
        }

        let tm = ::hgvs_weaver::transcript_mapper::TranscriptMapper::new(transcript)?;
        let n_start = tm.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let n_end = if let Some(e) = &pos.end {
            tm.c_to_n(e.base.to_index(), e.anchor)?
        } else {
            n_start
        };

        let start_idx = n_start.0 as usize;
        let end_idx = (n_end.0 + 1) as usize;
        if start_idx >= ref_seq.len() || end_idx > ref_seq.len() {
            return Err(HgvsError::ValidationError("Transcript sequence too short".into()));
        }
        let sub_seq = &ref_seq[start_idx .. end_idx];

        match &v.posedit.edit {
            ::hgvs_weaver::edits::NaEdit::RefAlt { ref_: Some(r), .. } => {
                if r.is_empty() || r.chars().all(|c| c.is_ascii_digit()) { return Ok(true); }
                Ok(r == sub_seq)
            }
            _ => Ok(true),
        }
    }
}

#[gen_stub_pyfunction]
#[pyfunction]
#[doc = "Parses an HGVS string into a Variant object.\n\nSupported types include genomic (g.), coding cDNA (c.), non-coding (n.),\nmitochondrial (m.), and protein (p.) variants.\n\nArgs:\n    input: The HGVS string to parse.\n\nReturns:\n    A Variant object.\n\nRaises:\n    ValueError: If the HGVS string is malformed or unsupported."]
fn parse(input: &str) -> PyResult<PyVariant> {
    match ::hgvs_weaver::parse_hgvs_variant(input) {
        Ok(inner) => Ok(PyVariant { inner }),
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e.to_string())),
    }
}

// --- Mapper and DataProvider Bridge ---

pub struct PyDataProviderBridge {
    provider: Py<PyAny>,
}

impl DataProvider for PyDataProviderBridge {
    fn get_seq(&self, ac: &str, start: i32, end: i32, kind: IdentifierKind) -> Result<String, HgvsError> {
        Python::attach(|py| {
            let kind_str = match kind {
                IdentifierKind::Genomic => "g",
                IdentifierKind::Transcript => "c",
                IdentifierKind::Protein => "p",
            };
            let res = self.provider.bind(py).call_method1("get_seq", (ac, start, end, kind_str))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<String>().map_err(|e| HgvsError::DataProviderError(e.to_string()))
        })
    }

    fn get_transcript(&self, transcript_ac: &str, reference_ac: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
        Python::attach(|py| {
            let res = self.provider.bind(py).call_method1("get_transcript", (transcript_ac, reference_ac))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;

            let json_str: String = match res.call_method0("to_json") {
                Ok(val) => val.extract::<String>().map_err(|e| HgvsError::DataProviderError(e.to_string()))?,
                Err(_) => {
                    let dict = res.cast::<PyDict>().map_err(|e| HgvsError::DataProviderError(format!("Failed to extract transcript: {}", e)))?;
                    let json_mod = py.import("json").map_err(|e| HgvsError::Other(e.to_string()))?;
                    let s = json_mod.call_method1("dumps", (dict,)).map_err(|e| HgvsError::Other(e.to_string()))?;
                    s.extract::<String>().map_err(|e| HgvsError::Other(e.to_string()))?
                }
            };

            let td: ::hgvs_weaver::data::TranscriptData = serde_json::from_str(&json_str).map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let boxed: Box<dyn Transcript> = Box::new(td);
            Ok(boxed)
        })
    }

    fn get_symbol_accessions(&self, symbol: &str, source_kind: IdentifierKind, target_kind: IdentifierKind) -> Result<Vec<String>, HgvsError> {
        Python::attach(|py| {
            let sk = match source_kind {
                IdentifierKind::Genomic => "g",
                IdentifierKind::Transcript => "c",
                IdentifierKind::Protein => "p",
            };
            let tk = match target_kind {
                IdentifierKind::Genomic => "g",
                IdentifierKind::Transcript => "c",
                IdentifierKind::Protein => "p",
            };
            let res = self.provider.bind(py).call_method1("get_symbol_accessions", (symbol, sk, tk))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<Vec<String>>().map_err(|e| HgvsError::DataProviderError(e.to_string()))
        })
    }
}

pub struct PyTranscriptSearchBridge {
    searcher: Py<PyAny>,
}

impl TranscriptSearch for PyTranscriptSearchBridge {
    fn get_transcripts_for_region(&self, chrom: &str, start: i32, end: i32) -> Result<Vec<String>, HgvsError> {
        Python::attach(|py| {
            let res = self.searcher.bind(py).call_method1("get_transcripts_for_region", (chrom, start, end))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<Vec<String>>().map_err(|e| HgvsError::DataProviderError(e.to_string()))
        })
    }
}

#[gen_stub_pyclass]
#[pyclass(name = "VariantMapper")]
#[doc = "High-level variant mapping engine.\n\nCoordinates mapping between different reference sequences (e.g., g. to c.)\nand projects cDNA variants onto protein sequences (c. to p.).\nRequires a DataProvider to retrieve transcript and sequence information."]
pub struct PyVariantMapper {
    pub bridge: std::sync::Arc<PyDataProviderBridge>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyVariantMapper {
    #[new]
    #[doc = "Creates a new VariantMapper with the given DataProvider."]
    fn new(provider: Py<PyAny>) -> Self {
        PyVariantMapper {
            bridge: std::sync::Arc::new(PyDataProviderBridge { provider }),
        }
    }

    #[pyo3(signature = (var_g, transcript_ac))]
    #[doc = "Maps a genomic variant (g.) to a coding cDNA variant (c.) for a specific transcript.\n\nArgs:\n    var_g: The genomic Variant to map.\n    transcript_ac: The accession of the target transcript.\n\nReturns:\n    A new Variant object in 'c.' coordinates."]
    fn g_to_c(&self, _py: Python, var_g: &PyVariant, transcript_ac: &str) -> PyResult<PyVariant> {
        if let SequenceVariant::Genomic(v) = &var_g.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper.g_to_c(v, transcript_ac).map_err(|e: HgvsError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            Ok(PyVariant { inner: SequenceVariant::Coding(res) })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err("Expected a genomic variant (g.)"))
        }
    }

    #[pyo3(signature = (var_g, searcher))]
    #[doc = "Maps a genomic variant (g.) to all overlapping transcripts discovered via the searcher.\n\nArgs:\n    var_g: The genomic Variant to map.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    A list of Variant objects in 'c.' coordinates."]
    fn g_to_c_all(&self, _py: Python, var_g: &PyVariant, searcher: Py<PyAny>) -> PyResult<Vec<PyVariant>> {
        if let SequenceVariant::Genomic(v) = &var_g.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let bridge = PyTranscriptSearchBridge { searcher };
            let res = mapper.g_to_c_all(v, &bridge).map_err(|e: HgvsError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            Ok(res.into_iter().map(|inner| PyVariant { inner: SequenceVariant::Coding(inner) }).collect())
        } else {
            Err(pyo3::exceptions::PyValueError::new_err("Expected a genomic variant (g.)"))
        }
    }

    #[pyo3(signature = (var_c, reference_ac))]
    #[doc = "Maps a coding cDNA variant (c.) to a genomic variant (g.) on a specific reference.\n\nArgs:\n    var_c: The coding Variant to map.\n    reference_ac: The accession of the target genomic reference.\n\nReturns:\n    A new Variant object in 'g.' coordinates."]
    fn c_to_g(&self, _py: Python, var_c: &PyVariant, reference_ac: &str) -> PyResult<PyVariant> {
        if let SequenceVariant::Coding(v) = &var_c.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper.c_to_g(v, reference_ac).map_err(|e: HgvsError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            Ok(PyVariant { inner: SequenceVariant::Genomic(res) })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err("Expected a coding variant (c.)"))
        }
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Projects a coding cDNA variant (c.) to its protein consequence (p.).\n\nArgs:\n    var_c: The coding Variant to project.\n    protein_ac: Optional protein accession. If not provided, it will be retrieved from the DataProvider.\n\nReturns:\n    A new Variant object in 'p.' coordinates."]
    fn c_to_p(&self, _py: Python, var_c: &PyVariant, protein_ac: Option<&str>) -> PyResult<PyVariant> {
        if let SequenceVariant::Coding(v) = &var_c.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper.c_to_p(v, protein_ac).map_err(|e: HgvsError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            Ok(PyVariant { inner: SequenceVariant::Protein(res) })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err("Expected a coding variant (c.)"))
        }
    }

    #[pyo3(signature = (var))]
    #[doc = "Normalizes a variant by shifting it to its 3'-most position.\n\nNormalization is performed in the coordinate space of the input variant.\n\nArgs:\n    var: The Variant object to normalize.\n\nReturns:\n    A new normalized Variant object."]
    fn normalize_variant(&self, _py: Python, var: &PyVariant) -> PyResult<PyVariant> {
        let mapper = VariantMapper::new(self.bridge.as_ref());
        let res = mapper.normalize_variant(var.inner.clone())
            .map_err(|e: HgvsError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        Ok(PyVariant { inner: res })
    }

    #[pyo3(signature = (var1, var2, searcher))]
    #[doc = "Determines if two variants are biologically equivalent.\n\nHandles normalization, cross-coordinate mapping (g. vs c.), and gene symbol expansion.\n\nArgs:\n    var1: The first Variant object.\n    var2: The second Variant object.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    True if the variants are equivalent, False otherwise."]
    fn equivalent(&self, _py: Python, var1: &PyVariant, var2: &PyVariant, searcher: Py<PyAny>) -> PyResult<bool> {
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let equiv = ::hgvs_weaver::equivalence::VariantEquivalence::new(self.bridge.as_ref(), &bridge_searcher);
        equiv.equivalent(&var1.inner, &var2.inner)
            .map_err(|e: HgvsError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }
}

#[pymodule]
fn _weaver(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_class::<PyVariant>()?;
    m.add_class::<PyVariantMapper>()?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);
