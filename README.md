# BAMs sintéticos de variantes — GRCh38 / Illumina WGS

Conjunto de BAMs sintéticos para **treinamento em interpretação de variantes no IGV**. Os exemplos reproduzem sinais visuais de SNVs, alelos complexos, deleções, variantes em cis/trans e duplicação.

> [!CAUTION]
> **Uso exclusivamente educacional.**

## Visão geral em 30 segundos

- Referência: **GRCh38/hg38**.
- Perfil simulado: **Illumina WGS paired-end de 150 bp (PE150), cerca de 30×**.
- Conteúdo: **1 BAM controle**, **8 BAMs individuais de cenário** e **1 BAM combinado**.
- Abrangência: somente as **janelas gênicas de interesse**, e não o genoma completo.
- Visualização recomendada: **IGV**, com o VCF esperado carregado junto aos BAMs.
- Reprodutibilidade: todos os arquivos podem ser recriados com `generate_bams.py`.

O caso de BTD é apresentado em dois BAMs — um em **cis** e outro em **trans** — por isso há oito BAMs individuais distribuídos em sete janelas gênicas.

## Download 

> [!TIP]
> [BAMs no Google Drive](https://drive.google.com/drive/folders/1sbXhIkcIHAKPLDYidC0eNIGgh8qrQbuX?usp=sharing)



## Glossário rápido

| Termo | Significado neste projeto |
|---|---|
| **BAM** | Arquivo que armazena os reads já alinhados à referência genômica |
| **Cobertura** | Número de reads que sustentam uma posição; por exemplo, 30× corresponde a cerca de 30 observações |
| **VAF** | Fração dos reads que contém o alelo alternativo em uma variante |
| **BAF** | Fração do alelo alternativo usada para acompanhar desequilíbrio alélico e número de cópias |
| **MAPQ** | Confiança de que o read foi alinhado na posição correta |
| **CIGAR** | Representação de matches, inserções, deleções e clips no alinhamento |
| **Soft-clip** | Trecho do read que foi sequenciado, mas não alinhou naquele local |
| **Read split** | Read representado em mais de um alinhamento, útil para indicar uma junção estrutural |

## O que podemos observar

| Cenário | Gene | O que procurar no IGV |
|---|---|---|
| SNV de boa qualidade | **BRCA2** | Alelo alternativo bem sustentado, bases com alta qualidade e MAPQ elevado |
| SNV em região homóloga | **CYP21A2** | MAPQ reduzido, menor qualidade no alelo alternativo e alinhamentos secundários próximos ao pseudogene **CYP21A1P** |
| Alelo complexo | **TMEM67** | Deleção de 1 bp associada a uma SNV adjacente |
| Deleção de 13 bp | **CYP1B1** | Soft-clips intencionais próximos ao evento |
| Duas variantes em cis | **BTD** | As duas variantes no mesmo haplótipo/read pair |
| Duas variantes em trans | **BTD** | Cada variante em um haplótipo diferente |
| Deleção heterozigótica de 214 bp | **KLC2** | Queda de cobertura, evidência de breakpoint e SNVs em heterozigose/hemizigose aparente |
| Duplicação tandem heterozigótica | **LAMA2** | Aumento de cobertura, reads split na junção e SNVs com fração alélica próxima de 2/3 |

## Início rápido no IGV

1. Abra o IGV e selecione o genoma **Human (GRCh38/hg38)**.
2. Carregue `output/00_control_reference.grch38.bam`.
3. Carregue o BAM do cenário que deseja estudar.
4. Carregue `output/expected_variants.grch38.vcf` para visualizar as variantes esperadas.
5. Compare as tracks de **controle** e **cenário**, observando cobertura, bases alternativas, CIGAR, pares discordantes e reads split.

Para começar, sugerimos esta ordem:

1. **BRCA2** — exemplo mais simples de SNV de boa qualidade;
2. **CYP1B1** — deleção pequena e soft-clips;
3. **BTD** — comparação entre cis e trans;
4. **KLC2** — deleção acompanhada de alteração da fração alélica;
5. **LAMA2** — integração entre cobertura, breakpoint e BAF.

## Regiões simuladas

Cada janela cobre continuamente o gene envolvido — éxons e íntrons — com aproximadamente 1 kb de flanco.

| Evento | Região GRCh38 |
|---|---|
| SNV em BRCA2 | `chr13:32,314,077-32,401,268` |
| SNV em CYP21A2 | `chr6:32,037,415-32,042,644` |
| Alelo complexo em TMEM67 | `chr8:93,753,844-93,819,121` |
| Deleção de 13 bp em CYP1B1 | `chr2:38,066,509-38,077,151` |
| Variantes em cis/trans em BTD | `chr3:15,600,361-15,723,516` |
| Deleção de 214 bp em KLC2 e flanco | `chr11:66,256,087-66,268,860` |
| Duplicação tandem em LAMA2 | `chr6:128,882,138-129,517,566` |

> [!NOTE]  
> As coordenadas exatas usadas na simulação também estão em `output/gene_windows.grch38.bed`.

## Dois exemplos para integrar cobertura e fração alélica

### KLC2: deleção e hemizigose aparente

O cenário contém uma deleção heterozigótica de 214 bp:

- **fora da deleção**, seis SNVs estão em HP1 e apresentam VAF esperado próximo de **50%**;
- **dentro da deleção**, três SNVs estão no alelo intacto HP2;
- como HP1 está ausente no intervalo deletado, as SNVs internas aparecem com VAF próximo de **100%**, reproduzindo **hemizigose aparente**.

Esse exemplo mostra por que uma fração alélica próxima de 100% nem sempre significa homozigose: a perda da outra cópia também pode produzir esse padrão.

### LAMA2: duplicação e BAF próxima de 2/3

Na região duplicada, foram inseridas 16 SNVs em HP1, o haplótipo duplicado. A simulação contém:

- uma cópia de referência em HP2;
- duas cópias alternativas em HP1;
- fração alélica esperada de `2/3 = 0,6667`.

Arquivos auxiliares:

- `output/08_LAMA2_duplication_BAF.tsv` — alelos REF/ALT, contagens, profundidade e BAF;
- `output/08_LAMA2_duplication_BAF.bedgraph` — BAF por posição, pronto para carregar no IGV.

> [!NOTE]  
> Neste projeto, **BAF** é calculado como `ALT_COUNT / (REF_COUNT + ALT_COUNT)`. O valor não é dobrado para o intervalo de 0 a 0,5.

## Arquivos principais

| Arquivo | Finalidade |
|---|---|
| `output/manifest.tsv` | Relação entre BAM, variante, genótipo, fase e evidência sintética |
| `output/expected_variants.grch38.vcf` | Genótipos e variantes esperados |
| `output/gene_windows.grch38.bed` | Coordenadas exatas das janelas simuladas |
| `output/00_control_reference.grch38.bam` | Controle negativo com sequência de referência |
| `output/00_all_scenarios.grch38.bam` | Todos os cenários reunidos em uma única track |
| `output/*.bam.bai` | Índices necessários para abrir os BAMs no IGV |

### Controle negativo

`00_control_reference.grch38.bam` reúne as sete janelas gênicas em um único sample/read group. Ele apresenta:

- aproximadamente 30× de cobertura;
- geometria Illumina WGS igual à linha de base dos cenários;
- dois haplótipos idênticos à referência GRCh38;
- ausência de SNVs sintéticas, deleções, duplicações, soft-clips intencionais e alinhamentos secundários ou suplementares.

> [!NOTE]  
> Em KLC2, o controle representa **CN = 2**, sem deleção. Em LAMA2, representa **CN = 2**, sem ganho de cobertura.

### BAM combinado

`00_all_scenarios.grch38.bam` reúne todos os loci em uma única track e preserva oito read groups e oito valores `SM`, um por cenário. Isso facilita o treinamento no IGV, pois os reads podem ser agrupados ou coloridos por cenário.

Alguns variant callers aceitam apenas um valor `SM`. Para esses programas, use os BAMs individuais ou reescreva deliberadamente todos os valores `SM` para um único nome.

## Configurações úteis no IGV

- Ative **View as pairs** para examinar o pareamento.
- Use **Color alignments by → read group** no BAM combinado.
- Inspecione as tags `HP` e `PS` nos exemplos de fase cis/trans.
- Ative **Show soft-clipped bases** para a deleção pequena em CYP1B1.
- Em LAMA2, visualize primeiro toda a região `chr6:129048000-129302500`; depois amplie cada breakpoint.
- Carregue o arquivo `.bedgraph` de LAMA2 para comparar BAF e profundidade na mesma região.

## Código

### Requisitos

- Python 3;
- FASTA GRCh38 com os contigs `chr2`, `chr3`, `chr6`, `chr8`, `chr11` e `chr13`.

### Execução

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python generate_bams.py --outdir output --reference /caminho/hg38.fa
```

Se `--reference` for omitido, o script baixa somente os segmentos necessários das sequências RefSeq que compõem GRCh38 e os armazena em `output/reference_cache`.

## Validar a geração

```bash
python validate_bams.py output
```

O validador verifica:

- integridade dos BAMs e de seus índices;
- cobertura contínua das janelas no controle;
- cobertura do controle próxima de 30×;
- concordância do controle com GRCh38 nos sítios didáticos;
- alelos e genótipos esperados;
- fase cis/trans e tags de haplótipo;
- CIGARs das deleções;
- VAFs heterozigóticas e de hemizigose aparente;
- BAF de LAMA2 próxima de 2/3;
- razão de cobertura e ausência de picos artificiais nos breakpoints;
- presença dos reads suplementares esperados.

<details>
<summary><strong>Detalhes técnicos da simulação</strong></summary>

### Perfil da biblioteca

- reads paired-end de 150 bp;
- orientação FR;
- insert size com distribuição aproximadamente normal truncada, `N(350,45)`, limitada a 250–500 bp;
- inícios distribuídos por toda a janela;
- MAPQ 60 em regiões não ambíguas;
- perfil de qualidade predominantemente Q30–Q40.

Esse desenho mantém cobertura interna contínua e evita a formação de duas ilhas artificiais de R1/R2.

### Decisões de modelagem

- Genótipos não especificados foram definidos como heterozigotos.
- Em CYP21A2, a baixa qualidade inclui MAPQ reduzido, menor qualidade de base no alelo alternativo e alinhamentos secundários próximos a CYP21A1P. O efeito não foi criado apenas pela alteração artificial do campo QUAL.
- O evento de TMEM67, por vezes chamado de MNV, é estritamente um alelo complexo formado por uma deleção de 1 bp e uma SNV adjacente. Os dois componentes estão em homozigose e coocorrem nos reads.
- A duplicação ClinVar 543888 possui limites imprecisos. Para visualização, foram adotados como breakpoints os limites internos reportados em GRCh38: `chr6:129049898-129300892`.
- A duplicação de LAMA2 foi modelada como tandem e heterozigótica, com aproximadamente 30× fora, 45× dentro e reads split na junção. A terceira cópia percorre circularmente o breakpoint tandem para manter a cobertura estável e evitar picos artificiais nas bordas.
- Os soft-clips de CYP1B1 são intencionais e didáticos. Um alinhador real pode representar a mesma deleção principalmente como `13D`, dependendo dos parâmetros e do contexto local.

</details>

## Limitações

- Os BAMs representam **somente regiões locais**, embora imitem características de uma biblioteca WGS.
- As evidências foram construídas para serem didáticas e visualmente reconhecíveis.
- O comportamento exato de um alinhador ou variant caller real pode variar conforme algoritmo, versão e parâmetros.
- Os arquivos não representam a diversidade completa de artefatos técnicos e biológicos encontrados em amostras clínicas.

## Uso recomendado

Este material pode ser usado em aulas, workshops e treinamentos para:

- reconhecer evidências de diferentes classes de variantes;
- relacionar cobertura, VAF/BAF e número de cópias;
- comparar variantes em cis e trans;
- discutir regiões homólogas e limitações de mapeamento;
- praticar inspeção visual antes da interpretação clínica.

---

**Dados exclusivamente sintéticos. Não usar para diagnóstico, validação clínica ou tomada de decisão sobre pacientes.**


## Anexos

**BAM e Eventos**

| **bam**                             | **scenario**             | **grch38_variant**                     | **genotype_phase**                  | **source**             | **synthetic_evidence**                                                                           |
| ----------------------------------- | ------------------------ | -------------------------------------- | ----------------------------------- | ---------------------- | ------------------------------------------------------------------------------------------------ |
| **00_control_reference.grch38** | controle negativo GRCh38 | todas as janelas do BED                | 0\|0 em todos os sítios             | GRCh38                 | um RG/SM; PE150; ~30x; CN=2; sem variantes, indels, soft-clips ou reads suplementares            |
| **00_all_scenarios.grch38**     | todos os oito cenários   | múltiplos loci                         | preserva genótipos/fase individuais | arquivos 01-08         | BAM combinado com oito RG/SM; recomendado para IGV                                               |
| **01_SNV_good_BRCA2**           | SNV boa qualidade        | chr13:32319080 T>G                     | 0\|1                                | ClinVar 266991         | BRCA2 completo; PE150; ~30x; MAPQ60; BQ variável                                                 |
| **02_SNV_low_quality_CYP21A2**  | SNV em região homóloga   | chr6:32038610 A>T                      | 0\|1                                | ClinVar 12183          | CYP21A2 completo; MAPQ0-20 no gene; alinhamentos secundários CYP21A1P                            |
| **03_complex_MNV_hom_cis_chr8** | alelo complexo/MNV       | chr8:93797350 TG>T + chr8:93797352 C>T | 1\|1 / cis                          | gnomAD r4              | TMEM67 completo; ambos os componentes nos dois haplótipos                                        |
| **04_small_del_hom_CYP1B1**     | deleção pequena          | chr2:38071278-38071290 del13           | 1\|1                                | ClinVar 282564         | CYP1B1 completo; PE150; CIGAR 13D e soft-clips intencionais                                      |
| **05_BTD_trans_pair_view**      | variantes em trans       | chr3:15644917 C>T / chr3:15645186 G>C  | 1\|0 / 0\|1                         | ClinVar 2230099 / 1900 | BTD completo; insert N(350,45); pares de fase adicionais e tags HP                               |
| **06_BTD_cis_pair_view**        | variantes em cis         | chr3:15644857 T>A / chr3:15645186 G>C  | 1\|0 / 1\|0                         | ClinVar 2203317 / 1900 | BTD completo; insert N(350,45); pares de fase adicionais e tags HP                               |
| **07_deletion_214bp_het**       | deleção heterozigótica   | chr11:66257087-66257300 del214         | 0\|1                                | ClinVar 1684657        | KLC2 completo; SNVs flanqueadores ~50% e SNVs internos hemizigóticos ~100%                       |
| **08_tandem_dup_251kb_het**     | duplicação tandem        | chr6:129049898-129300892 dup           | 0\|1 assumido                       | ClinVar 543888         | LAMA2 completo; terceira cópia circularizada; SNVs no haplótipo duplicado ~66%; TSV/bedGraph BAF |

## Formatos Aceitos pelo IGV

| Categoria                | Formato/extensão                                                      | Descrição                                                                                                                                                                 | Arquivo auxiliar/índice                                                              |
| ------------------------ | --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| **Alinhamento**          | **BAM (`.bam`)**                                                      | Formato binário de alinhamentos. Mostra reads, cobertura, pares, qualidade, mismatches, indels e orientação. É o formato recomendado para dados Illumina/WGS/WES/RNA-seq. | **Obrigatório:** `.bam.bai` ou `.bai`. O BAM deve estar ordenado por coordenada.     |
| **Alinhamento**          | **CRAM (`.cram`)**                                                    | Versão comprimida dos alinhamentos, geralmente dependente da sequência de referência.                                                                                     | **Obrigatório:** `.cram.crai` ou `.crai`; referência compatível.                     |
| **Alinhamento**          | **SAM (`.sam`)**                                                      | Versão textual do BAM. Aceito, mas ocupa muito espaço e tem desempenho inferior.                                                                                          | Deve estar ordenado e indexado com índice IGV `.sai`.                                |
| **Conjunto de BAMs**     | **BAM list (`.bam.list`)**                                            | Arquivo-texto contendo caminhos ou URLs de vários BAMs que serão exibidos como uma única track.                                                                           | Cada BAM precisa de seu respectivo `.bai`.                                           |
| **Variantes**            | **VCF (`.vcf`, `.vcf.gz`)**                                           | Representa SNVs, indels, genótipos e variantes estruturais. O IGV suporta VCF versão 4.                                                                                   | Recomendado: `.idx` para VCF simples ou `.tbi`/`.csi` para VCF comprimido com bgzip. |
| **Mutações somáticas**   | **MAF (`.maf`)**                                                      | Mutation Annotation Format, comum em estudos de câncer e dados do TCGA/GDC. Lista variantes e suas anotações.                                                             | Normalmente não exige índice.                                                        |
| **Regiões/anotações**    | **BED (`.bed`)**                                                      | Define intervalos genômicos, genes, éxons, regiões-alvo, CNVs, ROH, deleções e duplicações.                                                                               | Para arquivos muito grandes, recomenda-se bgzip + Tabix ou conversão para bigBed.    |
| **Regiões pareadas/SV**  | **BEDPE (`.bedpe`)**                                                  | Representa pares de intervalos, interações genômicas e breakpoints de variantes estruturais. Pode ser mostrado como arcos ou blocos.                                      | Geralmente não obrigatório.                                                          |
| **Anotação gênica**      | **GFF2 (`.gff`)**                                                     | Anota genes e outras características genômicas.                                                                                                                           | Opcionalmente comprimido e indexado para arquivos grandes.                           |
| **Anotação gênica**      | **GFF3 (`.gff3`)**                                                    | Formato hierárquico para genes, transcritos, éxons, CDS e outras features.                                                                                                | Opcionalmente `.gz` + `.tbi`.                                                        |
| **Anotação gênica**      | **GTF (`.gtf`)**                                                      | Muito usado para anotações de genes e transcritos em RNA-seq, como GENCODE e Ensembl.                                                                                     | Opcionalmente `.gz` + `.tbi`.                                                        |
| **Anotação gênica**      | **genePred / genePredExt / refGene / refFlat**                        | Formatos tabulares do UCSC para genes, transcritos, éxons e CDS.                                                                                                          | Pode ser comprimido e indexado.                                                      |
| **Anotação compactada**  | **bigBed (`.bb`, `.bigBed`)**                                         | Versão binária e indexada do BED. Indicada para grandes conjuntos de intervalos.                                                                                          | Índice já incorporado ao arquivo.                                                    |
| **Sinal quantitativo**   | **bedGraph (`.bedgraph`)**                                            | Valores contínuos associados a intervalos. Útil para cobertura, BAF, metilação, scores e expressão.                                                                       | Pode ser convertido para bigWig para melhor desempenho.                              |
| **Sinal quantitativo**   | **WIG (`.wig`)**                                                      | Representa sinais contínuos em posições ou intervalos genômicos.                                                                                                          | Recomenda-se converter para bigWig ou TDF.                                           |
| **Sinal quantitativo**   | **bigWig (`.bw`, `.bigWig`)**                                         | Formato binário e indexado para cobertura e outros sinais contínuos. Muito eficiente para arquivos grandes.                                                               | Índice incorporado.                                                                  |
| **Sinal otimizado**      | **TDF (`.tdf`)**                                                      | Formato binário criado pelo `igvtools`, otimizado para visualização rápida de dados quantitativos.                                                                        | Índice incorporado.                                                                  |
| **Copy number**          | **SEG/CBS (`.seg`, `.cbs`)**                                          | Segmentos cromossômicos com valores numéricos, como log₂ ratio, copy number ou BAF segmentado.                                                                            | Não costuma exigir índice.                                                           |
| **Copy number**          | **CN (`.cn`)**                                                        | Formato tabular para número de cópias ou log₂ tumor/normal.                                                                                                               | Não costuma exigir índice.                                                           |
| **LOH**                  | **LOH (`.loh`)**                                                      | Representa perda de heterozigosidade. Usa valores entre retenção e perda de heterozigosidade.                                                                             | Não costuma exigir índice.                                                           |
| **Copy number agregado** | **GISTIC (`.gistic`)**                                                | Regiões recorrentes de amplificação ou deleção, com score e q-value.                                                                                                      | Não costuma exigir índice.                                                           |
| **ChIP-seq/ATAC-seq**    | **narrowPeak (`.narrowPeak`)**                                        | Regiões estreitas de enriquecimento, como picos de fatores de transcrição ou ATAC-seq.                                                                                    | Pode ser convertido para bigNarrowPeak.                                              |
| **ChIP-seq**             | **broadPeak (`.broadPeak`)**                                          | Regiões amplas de enriquecimento, como determinadas marcas de histonas.                                                                                                   | Pode ser convertido para bigBed.                                                     |
| **GWAS**                 | **GWAS/PLINK (`.gwas`, `.assoc`, `.qassoc`, `.linear`, `.logistic`)** | Resultados de associação contendo cromossomo, posição, identificador da variante e p-value.                                                                               | Normalmente não obrigatório.                                                         |
| **Alinhamento BLAT**     | **PSL (`.psl`)**                                                      | Resultados de alinhamentos produzidos por BLAT ou ferramentas compatíveis.                                                                                                | Geralmente não obrigatório.                                                          |
| **Dados tabulares**      | **IGV (`.igv`)**                                                      | Formato tabular genérico com cromossomo, início, fim e valores para uma ou mais amostras.                                                                                 | Deve estar ordenado; arquivos grandes podem ser convertidos para TDF.                |
| **Genoma de referência** | **FASTA (`.fa`, `.fasta`)**                                           | Sequência do genoma de referência usada para exibir bases e interpretar alinhamentos.                                                                                     | **Obrigatório/recomendado:** `.fai`; o IGV pode tentar criá-lo.                      |
| **Genoma de referência** | **2bit (`.2bit`)**                                                    | Representação binária compactada da sequência de referência.                                                                                                              | Índice incorporado.                                                                  |
| **Definição de genoma**  | **JSON (`.json`)**                                                    | Define um genoma personalizado, incluindo FASTA/2bit, índice, citobandas, aliases e tracks de anotação. Substitui o antigo `.genome`.                                     | Depende dos arquivos referenciados no JSON.                                          |
| **Citobandas**           | **Cytoband (`.txt`, `.txt.gz`)**                                      | Define as bandas citogenéticas e o ideograma dos cromossomos.                                                                                                             | Não costuma exigir índice.                                                           |
| **Sessão IGV**           | **Session (`.xml`)**                                                  | Salva o genoma selecionado, arquivos carregados, regiões e configurações visuais.                                                                                         | Os arquivos referenciados precisam continuar acessíveis.                             |
