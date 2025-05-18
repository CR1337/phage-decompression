# phage-decompression

## Example usage

- Decompress all genomes in the directory `GenomesDB`:
```bash
python3 phagetools.py -op decompress -n '.*' GenomesDB
```

- Plot the genome `AB002632` in the directory `GenomesDB`:
```bash
python3 phagetools.py -op plot -n 'AB002632' GenomesDB
```

- Analyze all genomes in the directory `GenomesDB` and put the results in the directory `analysis`:
```bash
python3 phagetools.py -op analyze -n '.*' -o analysis GenomesD
```
