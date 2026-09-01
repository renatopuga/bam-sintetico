# Fase de variantes em cis e trans

[← Variantes pequenas](01_VARIANTES_PONTUAIS_E_INDELS.md) · [Roteiro](../AULA_HANDS_ON.md) · [Próximo: CNV e BAF →](03_CNV_E_BAF.md)

## Conceito

- **Em cis:** as duas variantes estão na mesma cópia/haplótipo.
- **Em trans:** cada variante está em uma cópia/haplótipo diferente.

A VAF de cada variante, isoladamente, não determina a fase. É necessário observar reads ou pares que conectem as posições, ou utilizar informação de haplótipo.

## Caso BTD em trans

- BAM: `05_BTD_trans_pair_view.grch38.bam`
- Região: `chr3:15644870-15645220`
- Variantes: `chr3:15644917 C>T` e `chr3:15645186 G>C`
- Fase: `1|0` e `0|1`

### Atividade

1. Ative **View as pairs**.
2. Amplie a região até visualizar os dois sítios.
3. Clique nos reads/pares informativos.
4. Compare as tags `HP` e `PS`.
5. Identifique o haplótipo que sustenta cada variante.

## Caso BTD em cis

- BAM: `06_BTD_cis_pair_view.grch38.bam`
- Região: `chr3:15644820-15645220`
- Variantes: `chr3:15644857 T>A` e `chr3:15645186 G>C`
- Fase: `1|0` e `1|0`

Repita a atividade anterior e procure pares que conectem os alelos alternativos no mesmo haplótipo.

## Quadro comparativo

| Pergunta | Cis | Trans |
|---|---|---|
| As variantes pertencem ao mesmo haplótipo? | Sim | Não |
| Um read/par informativo pode sustentar os dois alelos alternativos? | Sim, quando abrange/conecta os sítios | Não; cada alelo alternativo acompanha um haplótipo diferente |
| A VAF individual resolve a fase? | Não | Não |

## Checkpoint

Explique a diferença entre os dois BAMs sem usar apenas o VCF. Sua resposta deve citar a coocorrência dos alelos nos reads/pares ou as tags de haplótipo.

