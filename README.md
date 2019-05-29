# detect_plasmid

---

An example preparation workflow is presented in `refseq_prepare.sh`

Workflow:
1. `refseq_cds_extractor.py` - extract the nucleotide coding sequences and taxonomy
2. `refseq_cds_filter.py` - filter the extracted cds file by minimum sequence length
3. `refseq_cds_balance.py` - concatenate the matched and not matched taxa files and balance the number of sequences in each class
4. `refseq_cds_savemat.py` - format the coding sequences tagged by taxa class and save into a `.mat` file for use with PyTorch, also partitions train, validation and test

