# GAGER: Gene regulatory network Assisted Gene Expression Restoration

**GAGER** an algorithm designed to compare GRNs under two different conditions and identify specific genes whose manipulation could shift gene expression from a source condition, such as diseased tissue, towards a target condition, such as healthy tissue.

---

## Contents
- **`Codes/`**: Contains all the necessary codes to run the analysis.
- **`Codes/Datasets/`**: [Download Datasets](https://drive.google.com/drive/folders/17P5WTmnLN7GFJXgmu8fK4aquMZexycOO?usp=sharing) and put those under this folder
- **`Codes/Networks/`**: [Download Networks](https://drive.google.com/drive/folders/1dgUsLAm5XRUvAbX7Zq59Z0uDDMqwQ1L8?usp=sharing) and put those under this folder

## Sample Run

When you run any of the Python files in Codes folder, the output will appear as follows:

```bash
$ python microglia_analysis.py
['Foxm1', 'Rn45s', 'Ptchd1', 'Esco2', 'Dmrtb1', 'Myb', 'Bhlhe40']
```

## Run General Pipeline

You can run the general pipeline using the command below. Ensure that your gene expression matrices are CSV files with genes in the columns and cells in the rows. If your files are in other formats, such as TSV, or if genes are in the rows and cells are in the columns, update `Codes/gager.py` by editing lines 14–17 accordingly.

```bash
$ python gager.py Datasets/heart_data_new/gene_expression_matrix_healthy.csv Datasets/heart_data_new/gene_expression_matrix_group1.csv Networks/heart_control_byscenic.csv Networks/heart_group1_byscenic.csv
```



