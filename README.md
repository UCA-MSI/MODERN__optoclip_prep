# MODERN

optiCLIP Pipeline for Paired-End iCLIP-seq Data

## Arborescence de travail

```text
.
├── draft_code/
│   ├── README.md
│   └── raw_bash/
└── snakemake_wf/
    ├── README.md
    ├── Snakefile
    ├── config/
    │   └── config.yaml
    ├── profiles/
    │   └── local/
    │       └── config.yaml
    ├── resources/
    ├── results/
    └── workflow/
        ├── envs/
        ├── rules/
        └── scripts/
```

- `draft_code/` : dépôt du code Bash brut avant migration.
- `snakemake_wf/` : structure dédiée au développement du workflow Snakemake.
