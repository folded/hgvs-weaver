"""Parity tests using toy data."""

import typing

import hgvs.dataproviders.interface
import hgvs.parser
import hgvs.variantmapper


class ToyDataProvider(hgvs.dataproviders.interface.Interface):
    """Mock DataProvider for parity testing."""

    def __init__(self) -> None:
        self.url = "toy://local"
        self.required_version = "1.1"
        super().__init__()
        self.sequences = {
            "NC_TOY.1": "A" * 100 + "ATG" + "G" * 100 + "TAG" + "A" * 100,
            "NM_PLUS.1": "ATGGGGCCCAAA",
            "NM_MINUS.1": "TTTGGGCCCCAT",
        }
        self.transcripts: dict[str, typing.Any] = {
            "NM_PLUS.1": {
                "hgnc": "TOY1",
                "cds_start_i": 0,
                "cds_end_i": 12,
                "strand": 1,
                "alt_ac": "NC_TOY.1",
                "alt_aln_method": "splign",
                "exons": [
                    {
                        "tx_start_i": 0,
                        "tx_end_i": 12,
                        "alt_start_i": 100,
                        "alt_end_i": 112,
                        "alt_strand": 1,
                        "cigar": "12=",
                    },
                ],
            },
            "NM_MINUS.1": {
                "hgnc": "TOY2",
                "cds_start_i": 0,
                "cds_end_i": 12,
                "strand": -1,
                "alt_ac": "NC_TOY.1",
                "alt_aln_method": "splign",
                "exons": [
                    {
                        "tx_start_i": 0,
                        "tx_end_i": 12,
                        "alt_start_i": 230,
                        "alt_end_i": 242,
                        "alt_strand": -1,
                        "cigar": "12=",
                    },
                ],
            },
        }

    def get_seq(self, ac: str, start_i: int, end_i: int) -> str:
        """Returns sequence for accession."""
        return str(self.sequences[ac][start_i:end_i])

    def get_tx_info(self, tx_ac: str, _alt_ac: str | None, _alt_aln_method: str | None) -> typing.Any:
        """Returns transcript info."""
        return self.transcripts[tx_ac]

    def get_tx_exons(self, tx_ac: str, _alt_ac: str | None, _alt_aln_method: str | None) -> list[dict[str, typing.Any]]:
        """Returns transcript exons."""
        return list(self.transcripts[tx_ac]["exons"])

    def get_tx_identity_info(self, tx_ac: str) -> dict[str, typing.Any]:
        """Returns transcript identity info."""
        tx = self.transcripts[tx_ac]
        return {
            "hgnc": tx["hgnc"],
            "lengths": [12],
            "tx_ac": tx_ac,
            "alt_acs": [tx["alt_ac"]],
            "cds_start_i": tx["cds_start_i"],
            "cds_end_i": tx["cds_end_i"],
        }

    def get_pro_ac_for_tx_ac(self, tx_ac: str) -> str:
        """Returns protein accession."""
        return tx_ac.replace("NM", "NP")

    def data_version(self) -> str:
        return "1.1"

    def schema_version(self) -> str:
        return "1.1"


def test_toy_parity() -> None:
    """Tests mapping parity using toy data."""
    hp = hgvs.parser.Parser()
    hdp = ToyDataProvider()
    vm = hgvs.variantmapper.VariantMapper(hdp)

    # Plus Strand: g.101A>T -> c.1A>T
    var_g = hp.parse_hgvs_variant("NC_TOY.1:g.101A>T")
    var_c = vm.g_to_c(var_g, "NM_PLUS.1")
    print(f"PLUS g_to_c: {var_g} -> {var_c}")
    assert str(var_c) == "NM_PLUS.1:c.1A>T"

    # Plus Strand: c.1A>T -> p.(Met1Leu)
    var_p = vm.c_to_p(var_c)
    print(f"PLUS c_to_p: {var_c} -> {var_p}")
    assert str(var_p) == "NP_PLUS.1:p.(Met1Leu)"

    # Minus Strand: g.236T>G -> c.-4A>C
    var_g_minus = hp.parse_hgvs_variant("NC_TOY.1:g.236T>G")
    var_c_minus = vm.g_to_c(var_g_minus, "NM_MINUS.1")
    print(f"MINUS g_to_c: {var_g_minus} -> {var_c_minus}")
    assert str(var_c_minus) == "NM_MINUS.1:c.-4A>C"

    # Minus Strand: c.-4A>C -> g.236T>G
    var_g_back = vm.c_to_g(var_c_minus, "NC_TOY.1")
    print(f"MINUS c_to_g: {var_c_minus} -> {var_g_back}")
    assert str(var_g_back) == "NC_TOY.1:g.236T>G"

    print("All Parity Tests PASSED!")


if __name__ == "__main__":
    test_toy_parity()
