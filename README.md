## Assessing the replication / novelty of TenK10K findings in IBDverse
### Bradley April 2026
Downloaded summary data from Albert > `data` \\
Clone repos
```
git clone https://github.com/BradleyH017/pi1_pairs-Tenk10K.git
git clone https://github.com/andersonlab/IBDVerse-sc-eQTL-code.git
```

### 1. Perform pi1 replication
Follow the instruction of https://github.com/BradleyH017/pi1_pairs-Tenk10K to calculate the replication rate (pi1) of each set of eQTLs from TenK10K cell types in IBDverse cell-types. \\
This produces the output file: `pi1_pairs-Tenk10K/results/pi1_all.txt`

### 2. Summarise and plot
Follow code in `scripts/summarise.r` to plot the pi1 results and assess the novelty of effector gene novelty in Tenk10K vs IBDverse.