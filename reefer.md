# Reefer pipeline

## Produce a set of repeat annotations for the genome.

```
#!/bin/bash

# Invoked by:
#
# CONTIGS=<contigs> sbatch censor-hum.q
#

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-00:00
#SBATCH --mem=32GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@domain.net

#SBATCH --array=0-125 #result of `ls *.fa | wc -l' less one.

module load wu-blast censor bioperl

FILES=($(ls $CONTIGS/*.fa))
censor -bprm cpus=8 -lib hum $CONTIGS/$(basename ${FILES[$SLURM_ARRAY_TASK_ID]})
```

Convert to RepeatMasker GFF format.
```
REPLIBPATH=/path/to/biolib/humrep.ref # -lib hum
cat *map | map2gff -lib $REPLIBPATH >${GENOME}.rm.gff
```

## Run reefer

```bash
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=20G

#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@domain.net

#SBATCH --array=0-837 # result of `ls *fasta | wc -l' less one.

FILES=($(ls $READPATH))
reefer	-reads ${READPATH}/${FILES[$SLURM_ARRAY_TASK_ID]} \
     	-reference $REF \
     	-suff ${REF}.sa \
     	-procs 8 \
     	-blasr $BLASR \
     	-min 250 \
     	-err ${FILES[$SLURM_ARRAY_TASK_ID]}.log
```

## Obtain sequence discordances

```
head -n 3 $(ls -1 results/*gff | head -n 1) >${REEFERRESULTS}.gff
grep --no-filename -v '^#' results/*gff >>${REEFERRESULTS}.gff
```

This is horribly inefficient, but easy to express. A better approach is to do each invocation of wring on a distinct pair of gff/sam and concatenate the results.
```
wring results/*sam <${REEFERRESULTS}.gff >${REEFERRESULTS}.fasta
```

better
```
parallel 'wring {.}.blasr.sam <{} >{/.}.reefer.fasta' ::: results/*gff
cat *.reefer.fasta >${REEFERRESULTS}.fasta && parallel 'rm -f {/.}.reefer.fasta' ::: results/*gff
```
or
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=20:00
#SBATCH --mem=4G

#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@domain.net

#SBATCH --array=0-722 # result of `ls *fasta | wc -l' less one.

FILES=($(ls ${RESULTS}/*.gff))
file=$(basename ${FILES[$SLURM_ARRAY_TASK_ID]} .gff)
wring	${RESULTS}/${file}.blasr.sam \
		< ${RESULTS}/${file}.gff \
		> ${file}.reefer.fasta
```

## Perform name mangling

```
mangle <${REEFERRESULTS}.fasta >${REEFERRESULTS}-mangled.fasta
```

## Perform sequence bundling

Get event sequences in file batches censor can handle.
```
bundle -in ${REEFERRESULTS}-mangled.fasta
```

## Perform censor analysis of sequence discordances

Generic case:
```bash
#!/bin/bash

# Invoked by:
#
# EVENTSEQPATH=<eventseqs> sbatch censor-hum.q
#

#SBATCH -p cpuq
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-00:00
#SBATCH --mem=32GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@domain.net

#SBATCH --array=0-8 #result of `ls *.fa | wc -l' less one.

module load wu-blast censor bioperl

FILES=($(ls $EVENTSEQPATH/*.fa))
censor -bprm cpus=8 -lib hum $EVENTSEQPATH/$(basename ${FILES[$SLURM_ARRAY_TASK_ID]})
```

L1 case:
```bash
#!/bin/bash

# Invoked by:
#
# EVENTSEQPATH=<eventseqs> REPLIBPATH=<lib> sbatch censor-arb.q
#

#SBATCH -p cpuq
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=2-00:00
#SBATCH --mem=32GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@domain.net

#SBATCH --array=0-8 #result of `ls *.fa | wc -l' less one.

module load wu-blast censor bioperl

FILES=($(ls $EVENTSEQPATH/*.fa))
censor -bprm cpus=16 -lib ${REPLIBPATH} $EVENTSEQPATH/$(basename ${FILES[$SLURM_ARRAY_TASK_ID]})
```

## Filter events for size

General case:
```
awk '{if ($11 > 80 && $12 > 90) print $0;}' *map > ${REEFERRESULTS}-mangled.all.map
```

L1 case:
```
awk '{if ($11 > 90) print $0;}' *map > ${REEFERRESULTS}-mangled.L1.map
```

## Perform name unmangling

```
mangle -unmangle ${REEFERRESULTS}-mangled.$TYPE.map <${REEFERRESULTS}-mangled.fasta > ${REEFERRESULTS}-unmangled.$TYPE.map
```

## Convert to GFF

map2gff is available by invoking `go get github.com/kortschak/quilt/map2gff`.

General case:
```
map2gff -lib $REPLIBPATH <${REEFERRESULTS}-unmangled.all.map >${REEFERRESULTS}-unmangled.all.gff
```

L1 case (only allow 3' abutting elements):
```
# Only allow repeats that are within 100bp of 3' end.
map2gff -lib $REPLIBPATH -class L1 <${REEFERRESULTS}-unmangled.L1.map | awk '{if ($14 < 100) print $0}' >${REEFERRESULTS}-unmangled.L1.gff
```

## Integrate results

```
rinse -in ${REEFERRESULTS}-unmangled.$TYPE.gff -ref ${GENOME}.rm.gff -map ${REEFERRESULTS}.gff -contigs ${GENOME}.mfa >${REEFERRESULTS}-unmangled-rinsed.$TYPE.gff
```

## Annotate repeat locations on reference

Provides read provenance and accounts unique events.

```
press -in ${REEFERRESULTS}-unmangled-rinsed.$TYPE.gff -ref ${REEFERRESULTS}.gff -gff ${REEFERRESULTS}-unmangled-rinsed.pressed.$TYPE.gff
```

## Find Target Site Duplications

```
catch -in ${REEFERRESULTS}-unmangled-rinsed.pressed.$TYPE.gff ${READPATH}/* > ${REEFERRESULTS}-unmangled-rinsed.pressed.$TYPE.tsd.gff
```
