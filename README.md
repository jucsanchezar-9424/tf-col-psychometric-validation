# Data availability

The individual-level dataset used in this study is **not shared** in this repository.

The data contain sensitive clinical and sociodemographic information from oncology patients and caregivers. Access and reuse are restricted under the approvals granted by the participating institutions' research ethics committees, which retain custody of the underlying clinical records. For this reason, only the analysis code and aggregate outputs (tables and figures) are made publicly available.

## To run the code

The scripts expect a single Excel file at:

```
data/tfcol_dataset.xlsx
```

with the variables described in `codebook.csv` (85 variables; the 20 TF-Col items are `imf01_IF`–`imf10_IF` and `imps01_Ips`–`imps10_Ips`, scored 0–4). The scripts read the first sheet and apply `janitor::clean_names()` internally.

Researchers who wish to discuss data access for legitimate scientific purposes may contact the corresponding author; any sharing would be subject to the relevant ethics approvals and institutional data-protection requirements.

