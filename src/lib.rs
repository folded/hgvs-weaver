use ::hgvs_weaver::{
    DataProvider, HgvsError, IdentifierKind, SequenceVariant, Transcript, TranscriptSearch,
    Variant as VariantTrait, VariantMapper,
};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::{Bound, PyErr};
use pyo3_stub_gen::{define_stub_info_gatherer, derive::*};
use serde_json;

pyo3::create_exception!(
    _weaver,
    TranscriptMismatchError,
    pyo3::exceptions::PyValueError
);

fn map_hgvs_error(e: HgvsError) -> PyErr {
    match e {
        HgvsError::TranscriptMismatch {
            expected,
            found,
            start,
            end,
        } => TranscriptMismatchError::new_err(format!(
            "expected {}, found {} at transcript indices {}..{}",
            expected, found, start, end
        )),
        e => pyo3::exceptions::PyValueError::new_err(e.to_string()),
    }
}

#[gen_stub_pyclass_enum]
#[pyclass(name = "IdentifierType", module = "weaver._weaver")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PyIdentifierType {
    GenomicAccession,
    TranscriptAccession,
    ProteinAccession,
    GeneSymbol,
    Unknown,
}

impl From<PyIdentifierType> for ::hgvs_weaver::data::IdentifierType {
    fn from(it: PyIdentifierType) -> Self {
        match it {
            PyIdentifierType::GenomicAccession => {
                ::hgvs_weaver::data::IdentifierType::GenomicAccession
            }
            PyIdentifierType::TranscriptAccession => {
                ::hgvs_weaver::data::IdentifierType::TranscriptAccession
            }
            PyIdentifierType::ProteinAccession => {
                ::hgvs_weaver::data::IdentifierType::ProteinAccession
            }
            PyIdentifierType::GeneSymbol => ::hgvs_weaver::data::IdentifierType::GeneSymbol,
            PyIdentifierType::Unknown => ::hgvs_weaver::data::IdentifierType::Unknown,
        }
    }
}

impl From<::hgvs_weaver::data::IdentifierType> for PyIdentifierType {
    fn from(it: ::hgvs_weaver::data::IdentifierType) -> Self {
        match it {
            ::hgvs_weaver::data::IdentifierType::GenomicAccession => {
                PyIdentifierType::GenomicAccession
            }
            ::hgvs_weaver::data::IdentifierType::TranscriptAccession => {
                PyIdentifierType::TranscriptAccession
            }
            ::hgvs_weaver::data::IdentifierType::ProteinAccession => {
                PyIdentifierType::ProteinAccession
            }
            ::hgvs_weaver::data::IdentifierType::GeneSymbol => PyIdentifierType::GeneSymbol,
            ::hgvs_weaver::data::IdentifierType::Unknown => PyIdentifierType::Unknown,
        }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyIdentifierType {
    fn __repr__(&self) -> String {
        format!("IdentifierType.{:?}", self)
    }
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
    fn __hash__(&self) -> u64 {
        let mut s = std::collections::hash_map::DefaultHasher::new();
        std::hash::Hash::hash(self, &mut s);
        std::hash::Hasher::finish(&s)
    }
}

#[gen_stub_pyclass_enum]
#[pyclass(name = "EquivalenceLevel", module = "weaver._weaver")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PyEquivalenceLevel {
    Identity,
    Analogous,
    Different,
    Unknown,
}

impl From<::hgvs_weaver::equivalence::EquivalenceLevel> for PyEquivalenceLevel {
    fn from(el: ::hgvs_weaver::equivalence::EquivalenceLevel) -> Self {
        match el {
            ::hgvs_weaver::equivalence::EquivalenceLevel::Identity => Self::Identity,
            ::hgvs_weaver::equivalence::EquivalenceLevel::Analogous => Self::Analogous,
            ::hgvs_weaver::equivalence::EquivalenceLevel::Different => Self::Different,
            ::hgvs_weaver::equivalence::EquivalenceLevel::Unknown => Self::Unknown,
        }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyEquivalenceLevel {
    fn __repr__(&self) -> String {
        format!("EquivalenceLevel.{:?}", self)
    }
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
    fn __hash__(&self) -> u64 {
        let mut s = std::collections::hash_map::DefaultHasher::new();
        std::hash::Hash::hash(self, &mut s);
        std::hash::Hasher::finish(&s)
    }
}

#[gen_stub_pyclass]
#[pyclass(name = "Variant", module = "weaver._weaver")]
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
            _ => Err(HgvsError::UnsupportedOperation(
                "Validation not implemented for this variant type".into(),
            )),
        };
        match result {
            Ok(is_valid) => Ok(is_valid),
            Err(e) => Err(map_hgvs_error(e)),
        }
    }

    #[doc = "Converts the variant to an SPDI string representation."]
    fn to_spdi(&self, _py: Python, provider: Py<PyAny>) -> PyResult<String> {
        let bridge = PyDataProviderBridge { provider };
        self.inner.to_spdi(&bridge).map_err(map_hgvs_error)
    }
}

impl PyVariant {
    fn validate_genomic(
        &self,
        v: &::hgvs_weaver::GVariant,
        bridge: &PyDataProviderBridge,
    ) -> Result<bool, HgvsError> {
        let pos = v
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start_0 = pos.start.base.to_index();
        let end_0 = pos
            .end
            .as_ref()
            .map_or(start_0 + 1, |e| e.base.to_index() + 1);

        let ref_seq = bridge.get_seq(
            &v.ac,
            start_0.0,
            end_0.0,
            IdentifierKind::Genomic.into_identifier_type(),
        )?;

        match &v.posedit.edit {
            ::hgvs_weaver::edits::NaEdit::RefAlt { ref_: Some(r), .. } => {
                if r.is_empty() || r.chars().all(|c| c.is_ascii_digit()) {
                    return Ok(true);
                }
                Ok(r == &ref_seq)
            }
            _ => Ok(true),
        }
    }

    fn validate_coding(
        &self,
        v: &::hgvs_weaver::CVariant,
        bridge: &PyDataProviderBridge,
    ) -> Result<bool, HgvsError> {
        let transcript = bridge.get_transcript(&v.ac, None)?;

        let pos = v
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;

        let ref_seq = bridge.get_seq(
            &v.ac,
            0,
            -1,
            IdentifierKind::Transcript.into_identifier_type(),
        )?;

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
            return Err(HgvsError::ValidationError(
                "Transcript sequence too short".into(),
            ));
        }
        let sub_seq = &ref_seq[start_idx..end_idx];

        match &v.posedit.edit {
            ::hgvs_weaver::edits::NaEdit::RefAlt { ref_: Some(r), .. } => {
                if r.is_empty() || r.chars().all(|c| c.is_ascii_digit()) {
                    return Ok(true);
                }
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
    fn get_transcript(
        &self,
        transcript_ac: &str,
        reference_ac: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        Python::attach(|py| {
            let res = self
                .provider
                .bind(py)
                .call_method1("get_transcript", (transcript_ac, reference_ac))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let dict = res.cast_into::<PyDict>().map_err(|e| {
                HgvsError::DataProviderError(format!("Failed to cast to PyDict: {}", e))
            })?;
            let json_mod = py
                .import("json")
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let json_str: String = json_mod
                .call_method1("dumps", (dict,))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?
                .extract::<String>()
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let data: ::hgvs_weaver::data::TranscriptData = serde_json::from_str(&json_str)
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            Ok(Box::new(data) as Box<dyn Transcript>)
        })
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: i32,
        kind: ::hgvs_weaver::data::IdentifierType,
    ) -> Result<String, HgvsError> {
        Python::attach(|py| {
            let py_kind: PyIdentifierType = kind.into();
            let res = self
                .provider
                .bind(py)
                .call_method1("get_seq", (ac, start, end, py_kind))
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<String>()
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))
        })
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        source_kind: IdentifierKind,
        target_kind: IdentifierKind,
    ) -> Result<Vec<(::hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
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
            let res = self
                .provider
                .bind(py)
                .call_method1("get_symbol_accessions", (symbol, sk, tk))
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))?;
            let raw_list: Vec<(Bound<'_, PyAny>, String)> = res
                .extract::<Vec<(Bound<'_, PyAny>, String)>>()
                .map_err(|e: PyErr| {
                    HgvsError::DataProviderError(format!(
                        "Failed to extract symbol accessions: {}",
                        e
                    ))
                })?;

            let mut result = Vec::new();
            for (type_any, ac) in raw_list {
                let it = if let Ok(s) = type_any.extract::<String>() {
                    match s.as_str() {
                        "genomic_accession" => {
                            ::hgvs_weaver::data::IdentifierType::GenomicAccession
                        }
                        "transcript_accession" => {
                            ::hgvs_weaver::data::IdentifierType::TranscriptAccession
                        }
                        "protein_accession" => {
                            ::hgvs_weaver::data::IdentifierType::ProteinAccession
                        }
                        "gene_symbol" => ::hgvs_weaver::data::IdentifierType::GeneSymbol,
                        _ => ::hgvs_weaver::data::IdentifierType::Unknown,
                    }
                } else if let Ok(py_it) = type_any.extract::<PyIdentifierType>() {
                    py_it.into()
                } else {
                    ::hgvs_weaver::data::IdentifierType::Unknown
                };
                result.push((it, ac));
            }
            Ok(result)
        })
    }

    fn get_identifier_type(
        &self,
        identifier: &str,
    ) -> Result<::hgvs_weaver::data::IdentifierType, HgvsError> {
        Python::attach(|py| {
            let res = self
                .provider
                .bind(py)
                .call_method1("get_identifier_type", (identifier,))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;

            // Try to extract as String first (for backward compatibility or simpler mocks)
            if let Ok(s) = res.extract::<String>() {
                return Ok(match s.as_str() {
                    "genomic_accession" => ::hgvs_weaver::data::IdentifierType::GenomicAccession,
                    "transcript_accession" => {
                        ::hgvs_weaver::data::IdentifierType::TranscriptAccession
                    }
                    "protein_accession" => ::hgvs_weaver::data::IdentifierType::ProteinAccession,
                    "gene_symbol" => ::hgvs_weaver::data::IdentifierType::GeneSymbol,
                    _ => ::hgvs_weaver::data::IdentifierType::Unknown,
                });
            }

            // Otherwise try to extract as the enum type
            let py_it: PyIdentifierType = res.extract::<PyIdentifierType>().map_err(|e| {
                HgvsError::DataProviderError(format!("Failed to extract IdentifierType: {}", e))
            })?;
            Ok(py_it.into())
        })
    }

    fn c_to_g(
        &self,
        transcript_ac: &str,
        pos: ::hgvs_weaver::structs::TranscriptPos,
        offset: ::hgvs_weaver::structs::IntronicOffset,
    ) -> Result<(String, ::hgvs_weaver::structs::GenomicPos), HgvsError> {
        Python::attach(|py| {
            let res = self
                .provider
                .bind(py)
                .call_method1("c_to_g", (transcript_ac, pos.0, offset.0))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let (ac, g_pos): (String, i32) = res
                .extract::<(String, i32)>()
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            Ok((ac, ::hgvs_weaver::structs::GenomicPos(g_pos)))
        })
    }
}

pub struct PyTranscriptSearchBridge {
    searcher: Py<PyAny>,
}

impl TranscriptSearch for PyTranscriptSearchBridge {
    fn get_transcripts_for_region(
        &self,
        chrom: &str,
        start: i32,
        end: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Python::attach(|py| {
            let res = self
                .searcher
                .bind(py)
                .call_method1("get_transcripts_for_region", (chrom, start, end))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<Vec<String>>()
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))
        })
    }
}

#[gen_stub_pyclass]
#[pyclass(name = "VariantMapper", module = "weaver._weaver")]
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
    fn g_to_c(&self, _py: Python, var_g: &PyVariant, transcript_ac: String) -> PyResult<PyVariant> {
        if let SequenceVariant::Genomic(v) = &var_g.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper.g_to_c(v, &transcript_ac).map_err(map_hgvs_error)?;
            Ok(PyVariant {
                inner: SequenceVariant::Coding(res),
            })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "Expected a genomic variant (g.)",
            ))
        }
    }

    #[pyo3(signature = (var_g, searcher))]
    #[doc = "Maps a genomic variant (g.) to all overlapping transcripts discovered via the searcher.\n\nArgs:\n    var_g: The genomic Variant to map.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    A list of Variant objects in 'c.' coordinates."]
    fn g_to_c_all(
        &self,
        _py: Python,
        var_g: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<Vec<PyVariant>> {
        if let SequenceVariant::Genomic(v) = &var_g.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let bridge_searcher = PyTranscriptSearchBridge { searcher };
            let res = mapper
                .g_to_c_all(v, &bridge_searcher)
                .map_err(map_hgvs_error)?;
            Ok(res
                .into_iter()
                .map(|v| PyVariant {
                    inner: SequenceVariant::Coding(v),
                })
                .collect())
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "Expected a genomic variant (g.)",
            ))
        }
    }

    #[pyo3(signature = (var_c, reference_ac = None))]
    #[doc = "Maps a coding cDNA variant (c.) to a genomic variant (g.).\n\nArgs:\n    var_c: The coding Variant to map.\n    reference_ac: Optional chromosomal accession. If not provided, the primary chromosome for the transcript will be used.\n\nReturns:\n    A new Variant object in 'g.' coordinates."]
    fn c_to_g(
        &self,
        _py: Python,
        var_c: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        if let SequenceVariant::Coding(v) = &var_c.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper
                .c_to_g(v, reference_ac.as_deref())
                .map_err(map_hgvs_error)?;
            Ok(PyVariant {
                inner: SequenceVariant::Genomic(res),
            })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "Expected a coding variant (c.)",
            ))
        }
    }

    #[pyo3(signature = (var_n, reference_ac = None))]
    #[doc = "Maps a non-coding cDNA variant (n.) to a genomic variant (g.).\n\nArgs:\n    var_n: The non-coding Variant to map.\n    reference_ac: Optional chromosomal accession.\n\nReturns:\n    A new Variant object in 'g.' coordinates."]
    fn n_to_g(
        &self,
        _py: Python,
        var_n: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        if let SequenceVariant::NonCoding(v) = &var_n.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper
                .n_to_g(v, reference_ac.as_deref())
                .map_err(map_hgvs_error)?;
            Ok(PyVariant {
                inner: SequenceVariant::Genomic(res),
            })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "Expected a non-coding variant (n.)",
            ))
        }
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Projects a coding cDNA variant (c.) to its protein consequence (p.).\n\nArgs:\n    var_c: The coding Variant to project.\n    protein_ac: Optional protein accession. If not provided, it will be retrieved from the DataProvider.\n\nReturns:\n    A new Variant object in 'p.' coordinates."]
    fn c_to_p(
        &self,
        _py: Python,
        var_c: &PyVariant,
        protein_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        if let SequenceVariant::Coding(v) = &var_c.inner {
            let mapper = VariantMapper::new(self.bridge.as_ref());
            let res = mapper
                .c_to_p(v, protein_ac.as_deref())
                .map_err(map_hgvs_error)?;
            Ok(PyVariant {
                inner: SequenceVariant::Protein(res),
            })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "Expected a coding variant (c.)",
            ))
        }
    }

    #[pyo3(signature = (var))]
    #[doc = "Normalizes a variant by shifting it to its 3'-most position.\n\nNormalization is performed in the coordinate space of the input variant.\n\nArgs:\n    var: The Variant object to normalize.\n\nReturns:\n    A new normalized Variant object."]
    fn normalize_variant(&self, _py: Python, var: &PyVariant) -> PyResult<PyVariant> {
        let mapper = VariantMapper::new(self.bridge.as_ref());
        let res = mapper
            .normalize_variant(var.inner.clone())
            .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner: res })
    }

    #[pyo3(signature = (var1, var2, searcher))]
    #[doc = "Determines if two variants are biologically equivalent.\n\nHandles normalization, cross-coordinate mapping (g. vs c.), and gene symbol expansion.\n\nArgs:\n    var1: The first Variant object.\n    var2: The second Variant object.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    True if the variants are equivalent, False otherwise."]
    fn equivalent(
        &self,
        _py: Python,
        var1: &PyVariant,
        var2: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<bool> {
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let equiv = ::hgvs_weaver::equivalence::VariantEquivalence::new(
            self.bridge.as_ref(),
            &bridge_searcher,
        );
        equiv
            .equivalent(&var1.inner, &var2.inner)
            .map_err(map_hgvs_error)
    }

    #[pyo3(signature = (var1, var2, searcher))]
    #[doc = "Determines the granular equivalence level of two variants.\n\nArgs:\n    var1: The first Variant object.\n    var2: The second Variant object.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    An EquivalenceLevel enum value."]
    fn equivalent_level(
        &self,
        _py: Python,
        var1: &PyVariant,
        var2: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<PyEquivalenceLevel> {
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let equiv = ::hgvs_weaver::equivalence::VariantEquivalence::new(
            self.bridge.as_ref(),
            &bridge_searcher,
        );
        let res = equiv
            .equivalent_level(&var1.inner, &var2.inner)
            .map_err(map_hgvs_error)?;
        Ok(res.into())
    }

    #[pyo3(signature = (var, unambiguous = false))]
    #[doc = "Converts a variant to a SPDI string format.\n\nArgs:\n    var: The Variant object to convert.\n    unambiguous: If True, expands the variant range to cover the entire ambiguous region of a repeat or homopolymer. Default is False."]
    fn to_spdi(&self, _py: Python, var: &PyVariant, unambiguous: bool) -> PyResult<String> {
        let mapper = VariantMapper::new(self.bridge.as_ref());
        mapper
            .to_spdi(&var.inner, unambiguous)
            .map_err(map_hgvs_error)
    }

    #[pyo3(signature = (var))]
    #[doc = "Converts a variant to an unambiguous SPDI string format.\n\nThis format is independent of specific shifting conventions (like 3' or 5' shifting)\nby expanding the variant range to cover the entire ambiguous region of a repeat or homopolymer."]
    fn to_spdi_unambiguous(&self, _py: Python, var: &PyVariant) -> PyResult<String> {
        let mapper = VariantMapper::new(self.bridge.as_ref());
        mapper.to_spdi(&var.inner, true).map_err(map_hgvs_error)
    }
}

#[pymodule]
fn _weaver(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_class::<PyVariant>()?;
    m.add_class::<PyVariantMapper>()?;
    m.add_class::<PyIdentifierType>()?;
    m.add_class::<PyEquivalenceLevel>()?;
    m.add(
        "TranscriptMismatchError",
        m.py().get_type::<TranscriptMismatchError>(),
    )?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);
