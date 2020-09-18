# Summary

Generate inversions with invertFREGENE, analyze with recombClust and compute accuracy.

# Generate simulated inversions 

Define scenarios: 5 lenghts and 9 inversion frequencies, with 100 simulations in each scenario.

```{bash}
lengths=(50000 100000 250000 500000 1000000)
freqs=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
sim=($(seq 1 1 100))
```

# Create the simulation folders

```{bash}
for l in ${lengths[@]}
  do
  for f in ${freqs[@]}
    do
    echo l_$l.f_$f
    mkdir l_$l.f_$f
    done
  done
```

# Run invertFREGENE

```{bash}
for t in $(seq 1 5 100) 
  do 
  echo block $t
  (
  for s in $(seq $t 1 $(($t+4)))
  do  
  for f in ${freqs[@]}
    do
    for l in ${lengths[@]}
      do
      ## Make blocks of 20 simulation to parallelize
      echo l_$l.f_$f.$s
      # Generate base population
      invertFREGENE -i in.xml -p par.xml -gn 10000 -sd $s -recombsd $s -o ./l_$l.f_$f/l_$l.f_$f.$s.rin.xml -s -freq > ./l_$l.f_$f/l_$l.f_$f.$s.log.txt

      # Generate inversion
      invertFREGENE -i ./l_$l.f_$f/l_$l.f_$f.$s.rin.xml -StopFreqOfInv $f -p par.xml -gn 500000 -sd $s \
      -recombsd $s -StartOfInv 750000 -o ./l_$l.f_$f/l_$l.f_$f.$s.rin_inv.xml -EndOfInv $(($l+750000)) -MaxFreqOfLostInv 0.1 -s -freq > ./l_$l.f_$f/l_$l.f_$f.$s.log_inv.txt
      #Generate samples
      SAMPLE -noshuffle -i ./l_$l.f_$f/l_$l.f_$f.$s.rin_inv.xml -oh ./l_$l.f_$f/l_$l.f_$f.$s.haplotypes -og ./l_$l.f_$f/l_$l.f_$f.$s.genotypes -sd $s -controls 2000 > ./l_$l.f_$f/l_$l.f_$f.$s.log_sample.txt

      # Generate additional files to run models
      sed -e '1,/Position of inverted chromosomes/d' ./l_$l.f_$f/l_$l.f_$f.$s.rin_inv.xml_InversionSummary.txt > ./l_$l.f_$f/l_$l.f_$f.$s.log_ids.txt
      Rscript convertinvertFREGENE2plink.R ./l_$l.f_$f/l_$l.f_$f.$s
      done 
    done 
  done
  ) & 
done
```

# Run recombClust

```{bash}
for t in $(seq 1 1 100) 
  do  
  for f in ${freqs[@]}
    do
    for l in ${lengths[@]}
      do
      ## Make blocks of 20 simulation to parallelize
      echo l_$l.f_$f.$s
      Rscript runLDModelSimulation.R ./l_$l.f_$f/l_$l.f_$f.$s $l 2> log.txt
      done 
    done 
done
```

# Evaluate accuracy 

```{r}
library(parallel)
library(ggplot2)

lengths <- c("50000", "100000", "250000", "500000", "1000000")
freqs <- seq(0.1, 0.9, 0.1)
```

Compute accuracy of recombClust

```{r}
res <- lapply(lengths, function(l)
  lapply(freqs, function(f){
    message(paste0("L: ", l, " f: ", f))
  
    mclapply(1:100, function(i) {
      load(paste0("l_", l, ".f_", f, "/l_", l, ".f_", f,".", i, ".modelres.Rdata"))
      max(mean(ids == (class == 2)), mean(ids == (class == 1)))
    }, mc.cores = 20)
    })
  )

InversionAges <- lapply(lengths, function(l)
  lapply(freqs, function(f){
    message(paste0("L: ", l, " f: ", f))
    
    sapply(1:100, function(i) {
      as.numeric(system2("grep", paste0("Inversion l_", l, ".f_", f, "/l_", l, ".f_", f,".", i, ".rin_inv.xml_InversionSummary.txt | tr -dc '0-9'"), stdout = TRUE))
    })
  })
)


resdf <- data.frame(acc = as.numeric(unlist(res)),
                    lengths = rep(lengths, each = length(freqs)*100),
                    freqs = rep(freqs, each = 100),
                    age = unlist(InversionAges))
save(resdf, file = "simulationRes.Rdata")
```
