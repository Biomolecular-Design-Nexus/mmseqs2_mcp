#!/bin/bash

mamba activate /home/xux/Desktop/ProteinMCP/ProteinMCP/mcp-servers/mmseqs2_mcp/env

# Download UniRef100 fasta file
# wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
# gunzip uniref100.fasta.gz

UNIREF100_DIR="/mnt/data/data_repository/bio-seq/protein/uniref100"

# Create MMseqs2 DB from fasta file
# mmseqs createdb $UNIREF100_DIR/uniref100.fasta $UNIREF100_DIR/uniref100.fasta.db

# Pad database for GPU search
# mmseqs makepaddedseqdb $UNIREF100_DIR/uniref100.fasta.db $UNIREF100_DIR/uniref100.fasta.db_padded

# Search for MSA given a query fasta file

# Step 1: Create query database
mmseqs createdb examples/DHFR.fasta examples/DHFR_db 

# Step 2: Search (creates result DB instead of m8)
mmseqs search examples/DHFR_db $UNIREF100_DIR/uniref100.fasta.db_padded examples/DHFR_result_db examples/tmp --gpu 1 --threads 64
# Increase the number of hits with -s (sensitivity) parameter, --num-iterations, -e (e-value), --max-seqs parameters to achieve a matching performance with Jackhmmer
# mmseqs search examples/DHFR_db $UNIREF100_DIR/uniref100.fasta.db_padded examples/DHFR_result_db examples/tmp --gpu 1 --threads 64 -s 7.5 --num-iterations 10 -e 0.001 --max-seqs 100000

# Step 3: Convert result to MSA
mmseqs result2msa examples/DHFR_db  $UNIREF100_DIR/uniref100.fasta.db_padded examples/DHFR_result_db examples/DHFR_msa_db --msa-format-mode 6

# Step 4: Extract to a3m file
mmseqs unpackdb examples/DHFR_msa_db examples/DHFR_msa --unpack-suffix .a3m

# Step 5: After unpackdb, concatenate if needed
cat examples/DHFR_msa/*.a3m >examples/DHFR.a3m
