# Variantes pontuais e indels

[← Preparação](00_PREPARACAO.md) · [Roteiro](../AULA_HANDS_ON.md) · [Próximo: fase cis/trans →](02_FASE_CIS_TRANS.md)

Compare cada cenário com `00_control_reference.grch38.bam` e mantenha `expected_variants.grch38.vcf` carregado.

## Caso 1 — BRCA2: SNV de boa qualidade

- BAM: `01_SNV_good_BRCA2.grch38.bam`
- Região: `chr13:32319050-32319110`
- Variante: `chr13:32319080 T>G`
- Genótipo: `0|1`

Observe o alelo alternativo bem sustentado, MAPQ elevado e qualidades de base consistentes.

**Checkpoint:** a proporção de reads alternativos é compatível com uma variante heterozigótica?

## Caso 2 — CYP21A2: região homóloga

- BAM: `02_SNV_low_quality_CYP21A2.grch38.bam`
- Região: `chr6:32038570-32038650`
- Variante: `chr6:32038610 A>T`
- Genótipo: `0|1`

Compare MAPQ, qualidade das bases e alinhamentos secundários com BRCA2. O contexto homólogo próximo a **CYP21A1P** reduz a confiança de mapeamento.

**Checkpoint:** a presença do alelo alternativo é suficiente para classificá-lo como confiável? Justifique usando pelo menos duas métricas.

## Caso 3 — TMEM67: alelo complexo

- BAM: `03_complex_MNV_hom_cis_chr8.grch38.bam`
- Região: `chr8:93797310-93797390`
- Eventos: `chr8:93797350 TG>T` e `chr8:93797352 C>T`
- Genótipo/fase: `1|1`, em cis

O evento combina uma deleção de 1 bp e uma SNV adjacente. Inspecione os reads em conjunto, em vez de interpretar cada linha do VCF isoladamente.

**Checkpoint:** como os reads demonstram que os dois componentes pertencem ao mesmo alelo complexo?

## Caso 4 — CYP1B1: deleção de 13 bp

- BAM: `04_small_del_hom_CYP1B1.grch38.bam`
- Região: `chr2:38071230-38071330`
- Evento: `chr2:38071278-38071290 del13`
- Genótipo: `1|1`

Ative **Show soft-clipped bases**. Observe o CIGAR `13D` e os soft-clips intencionais próximos ao evento.

**Checkpoint:** por que alinhadores diferentes podem representar a mesma deleção de maneiras visuais um pouco diferentes?

## Comparação final

| Caso | Alelo alternativo | MAPQ/BQ | CIGAR/clip | Principal cautela |
|---|---|---|---|---|
| BRCA2 | Bem sustentado | Altos | Match | VAF e qualidade |
| CYP21A2 | Presente | Reduzidos | Alinhamentos secundários | Homologia/pseudogene |
| TMEM67 | Dois componentes | Consistentes | Deleção + mismatch | Interpretar o alelo completo |
| CYP1B1 | Deleção | Consistentes | `13D` + soft-clips | Representação pelo alinhador |

