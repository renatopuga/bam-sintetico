# Preparação do IGV

[← Início](../README.md) · [Roteiro da aula](../AULA_HANDS_ON.md)

## 1. Organize os arquivos

Mantenha cada BAM no mesmo diretório de seu índice:

```text
output/
├── 00_control_reference.grch38.bam
├── 00_control_reference.grch38.bam.bai
├── expected_variants.grch38.vcf
├── gene_windows.grch38.bed
└── ...
```

O IGV precisa do índice para navegar rapidamente por um BAM ordenado por coordenada.

## 2. Selecione a referência

1. Abra o IGV.
2. Selecione **Genomes → Human (GRCh38/hg38)**.
3. Confirme que os nomes dos contigs incluem o prefixo `chr`, por exemplo `chr7` e `chr13`.

## 3. Carregue as tracks básicas

1. Carregue `00_control_reference.grch38.bam`.
2. Carregue o BAM do cenário estudado.
3. Carregue `expected_variants.grch38.vcf`.
4. Quando indicado no roteiro, adicione o `.bed` ou `.bedgraph` correspondente.

## 4. Ajustes úteis

- **View as pairs**: mostra o pareamento dos reads.
- **Color alignments by → read group**: separa cenários no BAM combinado.
- **Show soft-clipped bases**: facilita a inspeção de CYP1B1.
- Tags `HP` e `PS`: auxiliam a inspeção dos casos em cis/trans.
- Ajuste manualmente a escala da track quantitativa quando comparar BAF.

## 5. Glossário mínimo

| Termo | Significado na atividade |
|---|---|
| BAM | Reads alinhados à referência genômica |
| Cobertura | Número de reads que sustentam uma posição |
| VAF | Fração dos reads com o alelo alternativo |
| BAF | Fração do alelo alternativo usada para avaliar desequilíbrio alélico |
| MAPQ | Confiança do posicionamento do read |
| CIGAR | Representação de matches, inserções, deleções e clips |
| Soft-clip | Parte sequenciada do read que não foi alinhada naquele local |
| Read split | Read representado em mais de um alinhamento, útil em junções estruturais |

## 6. Teste rápido

Navegue para `chr13:32319050-32319110`. Se o controle apresenta cobertura e o VCF marca a posição esperada, o ambiente está pronto.

### Problemas frequentes

| Sintoma | Verificação |
|---|---|
| BAM não abre | Confirme a presença do `.bam.bai` e a ordenação por coordenada |
| Região vazia | Confirme GRCh38/hg38 e o prefixo `chr` |
| VCF não aparece | Verifique se as coordenadas e a montagem são GRCh38 |
| BAF parece uma linha única | Ajuste a altura e a escala da track bedGraph |

