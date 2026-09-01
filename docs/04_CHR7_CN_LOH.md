# Cromossomo 7 com perda de heterozigosidade copy-neutral

[← CNV e BAF](03_CNV_E_BAF.md) · [Roteiro](../AULA_HANDS_ON.md) · [Próximo: reprodução e validação →](05_REPRODUCAO_VALIDACAO.md)

Conjunto sintético para treinamento na identificação de **perda de heterozigosidade sem perda do número de cópias**, também chamada **copy-neutral LOH** ou **CN-LOH**.

> **Ideia central:** o cromossomo 7 continua com duas cópias (`CN = 2`) e cobertura próxima de 30×, mas a banda heterozigótica de BAF em 0,5 desaparece.

## Objetivo da atividade

Distinguir CN-LOH de uma deleção heterozigótica usando conjuntamente cobertura, BAF, genótipos e intervalo de ROH.

## Arquivos a carregar

| Ordem | Arquivo | Finalidade |
|---:|---|---|
| 1 | `09_chr7_copy_neutral_LOH.grch38.bam` | Cobertura e alinhamentos do cromossomo 7 |
| 2 | `09_chr7_copy_neutral_LOH.roh.bed` | Extensão do ROH/CN-LOH |
| 3 | `09_chr7_copy_neutral_LOH.markers.vcf` | Marcadores e genótipos esperados |
| 4 | `09_chr7_copy_neutral_LOH.BAF.bedgraph` | Track quantitativa de BAF |

O índice `09_chr7_copy_neutral_LOH.grch38.bam.bai` deve permanecer no mesmo diretório do BAM.

## Passo a passo no IGV

1. Selecione **Human (GRCh38/hg38)**.
2. Carregue os quatro arquivos na ordem indicada.
3. Visualize primeiro o cromossomo 7 completo.
4. Observe a estabilidade da cobertura.
5. Examine a distribuição da track de BAF.
6. Amplie marcadores `HOM_REF` e `HOM_ALT` no VCF.
7. Verifique se há queda persistente de cobertura ou breakpoint de deleção.

## Resultado esperado

| Evidência | Resultado | Interpretação |
|---|---|---|
| Cobertura | Próxima de 30× | Número de cópias mantido (`CN = 2`) |
| BAF | Bandas próximas de 0,02 e 0,98 | Predomínio de genótipos homozigotos |
| Banda em 0,5 | Ausente | Perda de heterozigosidade |
| VCF | Marcadores `HOM_REF` ou `HOM_ALT` | Padrão alélico simulado |
| CIGAR/read split de perda | Ausente | Sem evidência de deleção estrutural |

## CN-LOH versus deleção

| Padrão | Número de cópias | Cobertura | Bandas de BAF |
|---|---:|---:|---|
| Região diploide normal | 2 | Mantida | 0, 0,5 e 1 |
| **CN-LOH** | **2** | **Mantida** | **0 e 1; sem 0,5** |
| Deleção heterozigótica | 1 | Reduzida | 0 e 1; sem 0,5 |

> [!IMPORTANT]
> Um ROH isolado não comprova CN-LOH. É necessário demonstrar que o número de cópias permanece igual a dois.

## Como a BAF foi simulada

Os marcadores alternam entre `HOM_REF`, com BAF próxima de 0, e `HOM_ALT`, com BAF próxima de 1. Um ruído de 2% evita bandas artificialmente perfeitas.

```text
BAF = ALT_COUNT / (REF_COUNT + ALT_COUNT)
```

Marcadores individuais podem se afastar das médias por amostragem binomial em profundidade próxima de 30×. Interprete o padrão de vários marcadores.

## Checklist do aluno

- [ ] A cobertura permanece aproximadamente estável.
- [ ] O BED cobre o intervalo esperado.
- [ ] A BAF forma bandas próximas de 0 e 1.
- [ ] Não existe banda persistente em 0,5.
- [ ] Os marcadores ampliados são essencialmente homozigotos.
- [ ] Não há breakpoint ou perda de cobertura compatível com deleção.

## Arquivos de apoio

| Arquivo | Conteúdo |
|---|---|
| `09_chr7_copy_neutral_LOH.roh.tsv` | Chamada em coordenadas 1-based inclusivas |
| `09_chr7_copy_neutral_LOH.BAF.tsv` | Alelos, genótipos, contagens, profundidade e BAF |
| `09_chr7_copy_neutral_LOH.validation.tsv` | Métricas observadas na validação |

As regiões ambíguas e os gaps da montagem GRCh38 não recebem reads. Pequenas interrupções nesses locais não devem ser interpretadas como deleções.

## Limitações

- O cenário é idealizado e ocupa todo o cromossomo 7.
- Pureza tumoral, mosaicismo e heterogeneidade clonal não são modelados de forma abrangente.
- Um BAM isolado não diferencia CN-LOH somática de mecanismos constitucionais, como dissomia uniparental.
- Ferramentas reais podem variar conforme algoritmo, versão, parâmetros e densidade de marcadores.

