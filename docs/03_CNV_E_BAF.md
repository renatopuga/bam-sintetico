# CNV, cobertura, VAF e BAF

[← Fase cis/trans](02_FASE_CIS_TRANS.md) · [Roteiro](../AULA_HANDS_ON.md) · [Próximo: CN-LOH →](04_CHR7_CN_LOH.md)

## Caso 1 — KLC2: deleção heterozigótica de 214 bp

- BAM: `07_deletion_214bp_het.grch38.bam`
- Evento: `chr11:66257087-66257300 del214`
- Genótipo: `0|1`

### O que foi simulado

- Fora da deleção, seis SNVs em HP1 apresentam VAF próxima de **50%**.
- Dentro da deleção, três SNVs pertencem ao alelo intacto HP2.
- Como HP1 está ausente no intervalo deletado, as SNVs internas apresentam VAF próxima de **100%**, reproduzindo hemizigose aparente.

### Atividade

1. Compare o BAM de KLC2 com o controle.
2. Observe a queda de cobertura no intervalo.
3. Inspecione SNVs fora e dentro da deleção.
4. Procure evidências do breakpoint.

> [!NOTE]
> VAF próxima de 100% não implica necessariamente homozigose. A perda da outra cópia também produz esse padrão.

## Caso 2 — LAMA2: duplicação tandem heterozigótica

- BAM: `08_tandem_dup_251kb_het.grch38.bam`
- Evento: `chr6:129049898-129300892 dup`
- Região inicial: `chr6:129048000-129302500`
- Genótipo: `0|1` assumido

### O que foi simulado

- cobertura próxima de 30× fora da duplicação;
- cobertura próxima de 45× dentro da duplicação;
- reads split na junção tandem;
- 16 SNVs no haplótipo duplicado HP1;
- uma cópia de referência em HP2 e duas cópias alternativas em HP1;
- BAF esperada de `2/3 = 0,6667`.

### Tracks auxiliares

- `08_LAMA2_duplication_BAF.tsv`
- `08_LAMA2_duplication_BAF.bedgraph`

Neste projeto:

```text
BAF = ALT_COUNT / (REF_COUNT + ALT_COUNT)
```

O valor não é dobrado para o intervalo de 0 a 0,5.

### Atividade

1. Visualize toda a região da duplicação.
2. Compare a cobertura com o controle.
3. Amplie cada breakpoint.
4. Procure reads suplementares/split na junção.
5. Carregue o bedGraph e compare BAF com cobertura.

## Integração das evidências

| Evento | CN esperado | Cobertura relativa | VAF/BAF didática | Evidência estrutural |
|---|---:|---:|---:|---|
| Diploide normal | 2 | 1,0 | ~0,5 em sítio heterozigoto | Ausente |
| Deleção heterozigótica | 1 | ~0,5 | ~1,0 no alelo remanescente | Breakpoint/queda de cobertura |
| Duplicação heterozigótica | 3 | ~1,5 | ~0,67 no haplótipo duplicado | Junção tandem/reads split |

## Checkpoint

Não conclua uma CNV usando apenas uma evidência. Para cada caso, combine ao menos cobertura, padrão alélico e alinhamentos próximos ao breakpoint.

