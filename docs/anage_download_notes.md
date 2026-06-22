# AnAge Download Notes

The first raw dataset was downloaded from the AnAge website on 2026-06-10.

- Homepage: https://genomics.senescence.info/species/
- Download URL used by script: https://genomics.senescence.info/species/dataset.zip
- Local raw file: `data/raw/anage/anage_data.zip`
- Normalized raw TSV: `data/interim/anage_raw.tsv`

The project should treat AnAge as the canonical first-pass source for maximum lifespan, body mass, and sexual maturity, followed by manual auditing of outliers and missing values.

Manual supplements should be used only when they include a source name or URL
and a note describing why the value is comparable to the AnAge field being
filled. Congeners should not be used as substitutes for species-level lifespan
or mass unless the downstream analysis explicitly switches to genus-level
imputation.
