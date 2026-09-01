# Aula hands-on — interpretação de BAMs sintéticos no IGV

Este é o roteiro principal do aluno. Abra apenas os documentos complementares indicados em cada etapa.

## Objetivos

Ao final da atividade, você deverá conseguir:

- comparar um BAM de cenário com um controle negativo;
- reconhecer SNVs, indels, soft-clips e reads split;
- diferenciar variantes em cis e trans;
- relacionar cobertura, VAF/BAF e número de cópias;
- distinguir uma deleção de uma perda de heterozigosidade copy-neutral.

## Antes de começar

Confirme que você possui:

- IGV com o genoma **Human (GRCh38/hg38)** selecionado;
- arquivos `.bam` acompanhados dos respectivos índices `.bam.bai`;
- `expected_variants.grch38.vcf`;
- tracks auxiliares `.bed`, `.bedgraph` e `.vcf` do cromossomo 7.

Se precisar configurar o ambiente, consulte [Preparação do IGV](docs/00_PREPARACAO.md).

## Sequência da aula

Tempo sugerido: **75–90 minutos**.

| Bloco | Tempo | Atividade | Evidência principal |
|---:|---:|---|---|
| 1 | 10 min | Controle e navegação | Cobertura ~30× e ausência de variantes |
| 2 | 20 min | Variantes pequenas | Bases alternativas, CIGAR, MAPQ e soft-clips |
| 3 | 15 min | Fase em BTD | Coocorrência das variantes nos reads/pares |
| 4 | 20 min | Deleção e duplicação | Cobertura, breakpoint e VAF/BAF |
| 5 | 15 min | CN-LOH no chr7 | Cobertura mantida e perda da banda em 0,5 |
| 6 | 10 min | Discussão | Síntese e limitações |

## Bloco 1 — controle e navegação

1. Carregue `Google Drive 00_control_reference.grch38.bam`.
2. Carregue `output/expected_variants.grch38.vcf`.
3. Navegue para `chr13:32319050-32319110`.
4. Ative **View as pairs**.
5. Clique em alguns reads e localize MAPQ, CIGAR, orientação e qualidade das bases.

### Verificação

- [ ] A cobertura é próxima de 30×.
- [ ] O controle não sustenta o alelo alternativo esperado em BRCA2.
- [ ] Os reads aparecem como pares em orientação compatível com biblioteca Illumina.

## Bloco 2 — variantes pequenas

Siga [Variantes pontuais e indels](docs/01_VARIANTES_PONTUAIS_E_INDELS.md) nesta ordem:

1. **BRCA2** — SNV de boa qualidade;
2. **CYP21A2** — SNV em região homóloga;
3. **TMEM67** — alelo complexo;
4. **CYP1B1** — deleção de 13 bp e soft-clips.

### Perguntas para registrar

1. Qual evidência diferencia a SNV de BRCA2 do sinal em CYP21A2?
2. O alelo de TMEM67 pode ser descrito adequadamente como uma única SNV?
3. Como o CIGAR e os soft-clips ajudam a interpretar CYP1B1?

## Bloco 3 — fase em cis e trans

Abra os BAMs de BTD e siga [Fase cis/trans](docs/02_FASE_CIS_TRANS.md).

### Perguntas para registrar

1. Em qual cenário as duas variantes aparecem no mesmo haplótipo?
2. Em qual cenário cada variante aparece em um haplótipo diferente?
3. Por que a VAF isolada não resolve a fase?

## Bloco 4 — CNV, VAF e BAF

Siga [CNV e BAF](docs/03_CNV_E_BAF.md):

1. compare controle e **KLC2**;
2. compare controle e **LAMA2**;
3. carregue a track de BAF de LAMA2.

### Perguntas para registrar

1. Por que uma SNV dentro da deleção de KLC2 pode apresentar VAF próxima de 100% sem ser homozigota?
2. Qual razão de cobertura é esperada para três cópias em relação ao controle diploide?
3. Por que as SNVs no haplótipo duplicado de LAMA2 apresentam BAF próxima de 2/3?

## Bloco 5 — CN-LOH no cromossomo 7

Siga [CN-LOH no cromossomo 7](docs/04_CHR7_CN_LOH.md).

### Perguntas para registrar

1. Qual evidência indica perda de heterozigosidade?
2. Qual evidência indica que não houve deleção?
3. Por que um ROH isolado não é suficiente para concluir CN-LOH?

## Síntese final

Complete a tabela sem consultar os documentos:

| Evento | Cobertura | Padrão alélico | Evidência decisiva |
|---|---|---|---|
| SNV heterozigótica |  |  |  |
| Deleção heterozigótica |  |  |  |
| Duplicação heterozigótica |  |  |  |
| CN-LOH |  |  |  |

## Critério de conclusão

Você concluiu a prática quando consegue justificar cada interpretação usando ao menos duas fontes de evidência, por exemplo:

- cobertura + VAF/BAF;
- alinhamento + qualidade;
- CIGAR + soft-clip;
- fase nos reads + tags de haplótipo.

> [!IMPORTANT]
> A inspeção visual é uma etapa de investigação e controle de qualidade. Ela não substitui validação analítica, ferramentas especializadas nem interpretação clínica.

