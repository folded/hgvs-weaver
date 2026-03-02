use pyo3_stub_gen::Result;
use std::fs;
use std::process::Command;

fn main() -> Result<()> {
    let stub = _weaver::stub_info()?;
    stub.generate()?;

    let pyi_path = "weaver/_weaver.pyi";
    if let Ok(content) = fs::read_to_string(pyi_path) {
        let mut lines: Vec<&str> = content.lines().collect();
        lines.retain(|line| !line.contains("def __repr__") && !line.contains("def __str__"));
        let mut new_content = lines.join("\n");
        new_content.push('\n');

        while let Some(start) = new_content.find("typing.Optional[") {
            let prefix_len = "typing.Optional[".len();
            if let Some(end_offset) = new_content[start + prefix_len..].find(']') {
                let end = start + prefix_len + end_offset;
                let inner = new_content[start + prefix_len..end].to_string();
                let replacement = format!("{} | None", inner);
                new_content.replace_range(start..=end, &replacement);
            } else {
                break;
            }
        }

        new_content = new_content.replace("builtins.str", "str");
        new_content = new_content.replace("builtins.int", "int");
        new_content = new_content.replace("builtins.bool", "bool");
        new_content = new_content.replace("builtins.float", "float");
        new_content = new_content.replace("builtins.list", "list");
        new_content = new_content.replace("builtins.dict", "dict");

        new_content = new_content.replace("= None", "= ...");
        new_content = new_content.replace("= False", "= ...");
        new_content = new_content.replace("= True", "= ...");
        new_content = new_content.replace("= StartCodonConvention.Specific", "= ...");

        new_content = new_content.replace(
            "def __eq__(self, other: EquivalenceLevel) -> bool",
            "def __eq__(self, other: object) -> bool",
        );
        new_content = new_content.replace(
            "def __eq__(self, other: IdentifierType) -> bool",
            "def __eq__(self, other: object) -> bool",
        );
        new_content = new_content.replace(
            "def __eq__(self, other: StartCodonConvention) -> bool",
            "def __eq__(self, other: object) -> bool",
        );

        new_content = new_content.replace("import builtins\n", "");

        fs::write(pyi_path, new_content).unwrap();
    }

    let _ = Command::new("uv")
        .args(&["run", "ruff", "format", pyi_path])
        .status();

    Ok(())
}
