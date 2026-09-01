# BAMs sintéticos de variantes — GRCh38 / Illumina WGS

Material educacional para treinamento prático de interpretação de variantes no IGV. Os dados reproduzem sinais visuais de SNVs, alelos complexos, fase em cis/trans, deleções, duplicação e perda de heterozigosidade copy-neutral.

> [!CAUTION]
> Dados exclusivamente sintéticos. Não usar para diagnóstico, validação clínica ou tomada de decisão sobre pacientes.

## Comece por aqui

Se você é aluno, siga o [roteiro da aula hands-on](AULA_HANDS_ON.md). Ele indica a ordem dos casos, os arquivos que devem ser carregados e as perguntas de interpretação.

Se você prepara o ambiente ou ministra a aula, consulte primeiro a [preparação do IGV](docs/00_PREPARACAO.md) e a [reprodução e validação dos dados](docs/05_REPRODUCAO_VALIDACAO.md).

## Visão geral

- Referência: **GRCh38/hg38**.
- Perfil simulado: **Illumina WGS paired-end 2 × 150 bp**, cobertura de aproximadamente **30×**.
- Conteúdo principal: **1 BAM controle**, **8 BAMs locais**, **1 BAM combinado** e **1 cenário de CN-LOH no cromossomo 7**.
- Abrangência dos cenários 01–08: janelas gênicas de interesse, e não o genoma completo.
- Visualização recomendada: **IGV**, com o VCF e as tracks auxiliares correspondentes.
- Reprodutibilidade: geração por `generate_bams.py` e `generate_chr7_loh.py`.

## Trilha de aprendizagem

| Etapa | Tema | Casos | Documento |
|---:|---|---|---|
| 0 | Preparar o IGV e entender as tracks | Controle | [Preparação](docs/00_PREPARACAO.md) |
| 1 | Reconhecer SNVs e indels | BRCA2, CYP21A2, TMEM67 e CYP1B1 | [Variantes pequenas](docs/01_VARIANTES_PONTUAIS_E_INDELS.md) |
| 2 | Comparar fase em cis e trans | BTD | [Fase cis/trans](docs/02_FASE_CIS_TRANS.md) |
| 3 | Integrar cobertura, VAF e BAF | KLC2 e LAMA2 | [CNV e BAF](docs/03_CNV_E_BAF.md) |
| 4 | Distinguir CN-LOH de deleção | Cromossomo 7 | [CN-LOH](docs/04_CHR7_CN_LOH.md) |

## Mapa dos arquivos

| Arquivo | Finalidade |
|---|---|
| `output/00_control_reference.grch38.bam` | Controle negativo, CN = 2 e sem variantes sintéticas |
| `output/00_all_scenarios.grch38.bam` | Casos 01–08 reunidos em uma única track |
| `output/01_*.bam` a `output/08_*.bam` | BAM individual de cada cenário |
| `output/expected_variants.grch38.vcf` | Variantes e genótipos esperados nos casos 01–08 |
| `output/gene_windows.grch38.bed` | Janelas simuladas |
| `output/09_chr7_copy_neutral_LOH.*` | BAM, BAF, ROH, VCF e validação do CN-LOH |
| `output/manifest.tsv` | Relação completa entre arquivos, variantes e evidências |

Consulte a [referência de cenários e arquivos](docs/06_REFERENCIA_ARQUIVOS.md) para coordenadas e detalhes de cada BAM. Para extensões aceitas pelo IGV, veja [formatos do IGV](docs/07_FORMATOS_IGV.md).

## Download

[BAMs no Google Drive](https://drive.google.com/drive/folders/1sbXhIkcIHAKPLDYidC0eNIGgh8qrQbuX?usp=sharing)

## Reproduzir e validar

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python generate_bams.py --outdir output --reference /caminho/hg38.fa
python validate_bams.py output
```

As instruções completas, incluindo o cromossomo 7, estão em [reprodução e validação](docs/05_REPRODUCAO_VALIDACAO.md).

