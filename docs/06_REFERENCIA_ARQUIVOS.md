# Referência de cenários e arquivos

[← Início](../README.md) · [Roteiro](../AULA_HANDS_ON.md)

## Regiões simuladas

| Evento | Região GRCh38 |
|---|---|
| SNV em BRCA2 | `chr13:32,314,077-32,401,268` |
| SNV em CYP21A2 | `chr6:32,037,415-32,042,644` |
| Alelo complexo em TMEM67 | `chr8:93,753,844-93,819,121` |
| Deleção de 13 bp em CYP1B1 | `chr2:38,066,509-38,077,151` |
| Variantes em cis/trans em BTD | `chr3:15,600,361-15,723,516` |
| Deleção de 214 bp em KLC2 | `chr11:66,256,087-66,268,860` |
| Duplicação tandem em LAMA2 | `chr6:128,882,138-129,517,566` |
| CN-LOH | Cromossomo 7 |

As coordenadas exatas dos casos 01–08 também estão em `output/gene_windows.grch38.bed`.

## Manifesto resumido

| BAM | Cenário | Variante/evento | Genótipo/fase | Evidência didática |
|---|---|---|---|---|
| `00_control_reference.grch38` | Controle | Todas as janelas | `0|0` | ~30×, CN=2 e sem eventos sintéticos |
| `00_all_scenarios.grch38` | Casos combinados | Múltiplos loci | Preservada | Oito RG/SM; recomendado para visão conjunta no IGV |
| `01_SNV_good_BRCA2` | SNV boa | `chr13:32319080 T>G` | `0|1` | MAPQ 60 e BQ variável |
| `02_SNV_low_quality_CYP21A2` | Região homóloga | `chr6:32038610 A>T` | `0|1` | MAPQ 0–20 e alinhamentos secundários |
| `03_complex_MNV_hom_cis_chr8` | Alelo complexo | `chr8:93797350 TG>T`; `chr8:93797352 C>T` | `1|1`, cis | Componentes nos dois haplótipos |
| `04_small_del_hom_CYP1B1` | Deleção pequena | `chr2:38071278-38071290 del13` | `1|1` | CIGAR `13D` e soft-clips |
| `05_BTD_trans_pair_view` | Variantes em trans | `chr3:15644917 C>T`; `chr3:15645186 G>C` | `1|0`; `0|1` | Pares informativos e tags HP |
| `06_BTD_cis_pair_view` | Variantes em cis | `chr3:15644857 T>A`; `chr3:15645186 G>C` | `1|0`; `1|0` | Pares informativos e tags HP |
| `07_deletion_214bp_het` | Deleção heterozigótica | `chr11:66257087-66257300 del214` | `0|1` | SNVs ~50% fora e ~100% dentro |
| `08_tandem_dup_251kb_het` | Duplicação tandem | `chr6:129049898-129300892 dup` | `0|1` assumido | Cobertura ~45×, reads split e BAF ~66% |
| `09_chr7_copy_neutral_LOH` | CN-LOH | Cromossomo 7 | CN=2, LOH | Cobertura ~30× e BAF sem banda em 0,5 |

## Controle negativo

`00_control_reference.grch38.bam` reúne as sete janelas gênicas em um único sample/read group, com dois haplótipos idênticos à referência. Não contém SNVs sintéticas, deleções, duplicações, soft-clips intencionais nem alinhamentos secundários ou suplementares.

## BAM combinado

`00_all_scenarios.grch38.bam` preserva oito read groups e oito valores `SM`. No IGV, os reads podem ser agrupados ou coloridos por cenário. Para ferramentas que aceitam apenas um valor `SM`, use os BAMs individuais ou normalize deliberadamente os read groups.

## Convenções de coordenadas do chr7

- `.roh.tsv`: coordenadas **1-based inclusivas**.
- `.roh.bed`: coordenadas **0-based half-open**.

Essa diferença é esperada e segue as convenções de cada formato.

