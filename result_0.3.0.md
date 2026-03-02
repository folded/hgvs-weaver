### Validation Results (100,000 variants)

Summary of results comparing `weaver` and `ref-hgvs` against ClinVar ground truth:

| Implementation | Protein Match | SPDI Match  | Parse Errors |
| :------------- | :-----------: | :---------: | :----------: |
| weaver         |  **93.849%**  | **97.968%** | **1** |
| ref-hgvs       |  93.337%  | 94.022% | 394 |

RefSeq Data Mismatches: 0 (0.0%)

#### Protein Translation Agreement

|                     | ref-hgvs Match | ref-hgvs Mismatch |
| :------------------ | :------------: | :---------------: |
| **weaver Match**    |     93,330     |     519     |
| **weaver Mismatch** |     7     |     6,144     |

#### SPDI Mapping Agreement

|                     | ref-hgvs Match | ref-hgvs Mismatch |
| :------------------ | :------------: | :---------------: |
| **weaver Match**    |     93,919     |     4,049     |
| **weaver Mismatch** |     103     |     1,929     |
