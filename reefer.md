# Reefer pipeline

## Run reefer

```bash
#!/bin/bash
#SBATCH -p cpuq
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
head -n 3 $(ls -1 | head -n 1) >${REEFERRESULTS}.gff
grep --no-filename results/*gff >>${REEFERRESULTS}.gff
```

This is horribly inefficient, but easy to express. A better approach is to do each invocation of wring on a distinct pair of gff/sam and concatenate the results.
```
wring results/*sam <${REEFERRESULTS}.gff >LCYE01.reefer.fasta
```

better
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=20:00
#SBATCH --mem=4G

#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=dan.kortschak@adelaide.edu.au

#SBATCH --array=0-837 # result of `ls *fasta | wc -l' less one.

FILES=($(ls $READPATH/*.gff))
file=$(basename ${FILES[$SLURM_ARRAY_TASK_ID]} .gff)
wring	${file}.blasr.sam \
     	< ${file}.gff \
     	> ${file}.reefer.fasta
```

## Perform name mangling

```
mangle <LCYE01.reefer.fasta >LCYE01.reefer-mangled.fasta
```

## Perform sequence bundling

```
bundle -in LCYE01.reefer-mangled.fasta
```

## Perform censor analysis of sequence discordances

Generic case:
```bash
#!/bin/bash

# Invoked by:
#
# READPATH=<read> sbatch censor-hum.q
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

FILES=($(ls $READPATH/*.fa))
censor -bprm cpus=8 -lib hum $READPATH/$(basename ${FILES[$SLURM_ARRAY_TASK_ID]})
```

L1 case:
```bash
#!/bin/bash

# Invoked by:
#
# READPATH=<read> REPLIBPATH=<lib> sbatch censor-arb.q
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

FILES=($(ls $READPATH/*.fa))
censor -bprm cpus=16 -lib ${REPLIBPATH} $READPATH/$(basename ${FILES[$SLURM_ARRAY_TASK_ID]})
```

## Filter events for size

General case:
```
awk '{if ($11 > 80 && $12 > 90) print $0;}' *map > LCYE01.reefer-mangled.all.map
```

L1 case:
```
awk '{if ($11 > 90) print $0;}' *map > LCYE01.reefer-mangled.L1.map
```

## Perform name unmangling

```
mangle -unmangle LCYE01.reefer-mangled.$TYPE.map <LCYE01.reefer-mangled.fasta > LCYE01.reefer-unmangled.$TYPE.map
```

## Convert to GFF

map2gff is available by invoking `go get github.com/kortschak/quilt/map2gff`.

General case:
```
map2gff -lib $REPLIBPATH <LCYE01.reefer-unmangled.all.map >LCYE01.reefer-unmangled.all.gff
```

L1 case (only allow 3' abutting elements):
```
map2gff -lib $REPLIBPATH -class L1 <LCYE01.reefer-unmangled.L1.map | awk '{if ($14 < 100) print $0}' >LCYE01.reefer-unmangled.L1.gff
```

## Integrate results

```
rinse -in LCYE01.reefer-unmangled.$TYPE.gff -ref LYCE01.rm.gff -map LCYE01.reefer.gff -contigs LCYE01.mfa >LCYE01.reefer-unmangled-rinsed.$TYPE.gff
```

## Anntotate repeat locations on reference

Provides read provenance and accounts unique events.

```
press -in LCYE01.reefer-unmangled-rinsed.$TYPE.gff -ref ${REEFERRESULTS}.gff -gff LCYE01.reefer-unmangled-rinsed.pressed.$TYPE.gff
```

## Find Target Site Duplications

```
catch -in LCYE01.reefer-unmangled-rinsed.pressed.$TYPE.gff ${READPATH}/* > LCYE01.reefer-unmangled-rinsed.pressed.$TYPE.tsd.gff
```