#!/bin/bash
NPROC=1

make > make.log
rm timings_block_*.txt 2>/dev/null
blockSizes=(1 2 4 8 16 32 64 128)
nzVals=(1 2 4 8 16 32 64 128 256)
for block in ${blockSizes[@]}
do
    for nz in ${nzVals[@]}
    do
	mpirun -np ${NPROC} ./sdi_contiguous -q -q mesh:nz=${nz} mesh:nx=128 contig:maxBlockSize=${block} | tail -n 1 | tee -a timings_block_${block}.txt
    done
done
