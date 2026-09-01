# Formatos utilizados no IGV

[← Início](../README.md) · [Preparação](00_PREPARACAO.md)

Durante esta aula, os formatos essenciais são:

| Categoria | Formato | Uso na aula | Índice/auxiliar |
|---|---|---|---|
| Alinhamento | BAM (`.bam`) | Cobertura, reads, pares, CIGAR, indels e orientação | `.bam.bai` ou `.bai` |
| Variantes | VCF (`.vcf`, `.vcf.gz`) | SNVs, indels, genótipos e eventos esperados | `.idx`, `.tbi` ou `.csi`, conforme o arquivo |
| Regiões | BED (`.bed`) | Janelas, ROH, deleções e duplicações | Opcional para arquivos pequenos |
| Sinal quantitativo | bedGraph (`.bedgraph`) | BAF por posição | Pode ser convertido para bigWig |
| Genoma | FASTA (`.fa`, `.fasta`) | Referência personalizada | `.fai` |
| Sessão | XML (`.xml`) | Estado de uma sessão do IGV | Arquivos referenciados devem permanecer acessíveis |

## Outros formatos comuns

| Categoria | Formato | Finalidade |
|---|---|---|
| Alinhamento | CRAM (`.cram`) | Alinhamentos comprimidos; normalmente dependem da referência |
| Alinhamento | SAM (`.sam`) | Representação textual, maior e mais lenta que BAM |
| Variante somática | MAF (`.maf`) | Lista de mutações e anotações, comum em câncer |
| Variante estrutural | BEDPE (`.bedpe`) | Pares de intervalos e breakpoints |
| Anotação gênica | GFF/GFF3/GTF | Genes, transcritos, éxons e CDS |
| Anotação compacta | bigBed (`.bb`) | Intervalos binários e indexados |
| Sinal quantitativo | WIG, bigWig, TDF | Cobertura e outros sinais contínuos |
| Copy number | SEG/CBS, CN | Segmentos de número de cópias ou log2 ratio |
| LOH | LOH (`.loh`) | Regiões de perda/retenção de heterozigosidade |
| Associação | GWAS/PLINK | Resultados de associação e p-values |
| Genoma | 2bit, JSON | Sequência compacta ou definição de genoma personalizado |

## Regras práticas

- BAM e CRAM devem estar ordenados e indexados.
- A montagem e os nomes dos contigs devem ser compatíveis entre todas as tracks.
- Para arquivos grandes, prefira formatos comprimidos/indexados como bgzip + Tabix, bigBed ou bigWig.
- Um arquivo de sessão XML não incorpora os dados; ele registra referências aos arquivos carregados.

