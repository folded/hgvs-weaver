### Validation Results (100,000 variants)

Summary of results comparing `weaver` and `ref-hgvs` against ClinVar ground truth:

| Implementation | Protein Match | SPDI Match  | Parse Errors |
| :------------- | :-----------: | :---------: | :----------: |
| weaver         |  **93.869%**  | **98.183%** | **1** |
| ref-hgvs       |  93.352%  | 94.039% | 394 |

RefSeq Data Mismatches: 0 (0.0%)

#### Protein Translation Agreement

|                     | ref-hgvs Match | ref-hgvs Mismatch |
| :------------------ | :------------: | :---------------: |
| **weaver Match**    |     93,345     |     524     |
| **weaver Mismatch** |     7     |     6,124     |

#### SPDI Mapping Agreement

|                     | ref-hgvs Match | ref-hgvs Mismatch |
| :------------------ | :------------: | :---------------: |
| **weaver Match**    |     93,936     |     4,247     |
| **weaver Mismatch** |     103     |     1,714     |
