use pyo3_stub_gen::Result;

fn main() -> Result<()> {
    let stub = _weaver::stub_info()?;
    stub.generate()?;
    Ok(())
}
