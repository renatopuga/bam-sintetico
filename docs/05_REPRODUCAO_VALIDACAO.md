# Reprodução e validação

[← CN-LOH](04_CHR7_CN_LOH.md) · [Início](../README.md) · [Roteiro](../AULA_HANDS_ON.md)

## Requisitos

- Python 3;
- FASTA GRCh38 com `chr2`, `chr3`, `chr6`, `chr8`, `chr11` e `chr13` para os casos locais;
- referência do cromossomo 7 para o cenário de CN-LOH.

## Ambiente

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Gerar os cenários 01–08

```bash
python generate_bams.py \
  --outdir output \
  --reference /caminho/hg38.fa
```

Se `--reference` for omitido, o script baixa apenas os segmentos necessários e os armazena em `output/reference_cache`.

## Gerar CN-LOH no cromossomo 7

```bash
python3 generate_chr7_loh.py \
  --reference reference/NC_000007.14.fa \
  --outdir output \
  --depth 30
```

## Validar

```bash
python validate_bams.py output
```

Consulte também:

```text
output/09_chr7_copy_neutral_LOH.validation.tsv
```

## O que o validador verifica

- integridade dos BAMs e índices;
- cobertura contínua das janelas no controle;
- cobertura do controle próxima de 30×;
- concordância com GRCh38 nos sítios didáticos;
- alelos, genótipos e fase esperados;
- CIGAR das deleções;
- VAF heterozigótica e hemizigose aparente;
- BAF de LAMA2 próxima de 2/3;
- razão de cobertura e ausência de picos artificiais nos breakpoints;
- reads suplementares esperados;
- profundidade e distribuição da BAF no cromossomo 7.

## Perfil técnico da simulação

- paired-end 2 × 150 bp;
- orientação FR;
- insert `N(350,45)` truncado em 250–500 bp;
- MAPQ 60 em regiões não ambíguas;
- qualidades predominantemente Q30–Q40 nos casos locais e Q28–Q40 no chr7;
- inícios distribuídos ao longo das janelas;
- ausência de reads nos gaps e bases ambíguas do chr7.

## Decisões de modelagem

- Genótipos não especificados foram definidos como heterozigotos.
- CYP21A2 inclui MAPQ reduzido, menor qualidade do alelo alternativo e alinhamentos secundários próximos a CYP21A1P.
- TMEM67 é um alelo complexo composto por deleção de 1 bp e SNV adjacente.
- A duplicação de LAMA2 usa os limites internos GRCh38 `chr6:129049898-129300892` e foi modelada como tandem heterozigótica.
- A terceira cópia de LAMA2 percorre circularmente a junção para manter cobertura estável e evitar picos artificiais.
- Os soft-clips de CYP1B1 são intencionais e didáticos.

## Limitações gerais

- Os casos 01–08 representam apenas regiões locais, embora imitem WGS.
- As evidências foram construídas para serem reconhecíveis em aula.
- Alinhadores e variant callers reais podem representar os eventos de maneira diferente.
- O material não cobre toda a diversidade de artefatos técnicos e biológicos de amostras clínicas.

