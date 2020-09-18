/*
 * Nextflow script to simulate inversions with invertFREGENE
 */

params.initialxml = ""
params.paramFile = ""
params.lens = 0
params.sim = 100
params.freqs = 0.5

paramFile = file(params.paramFile)
initialxml = file(params.initialxml)

// Create channels with simulations
sims = Channel.from( 1.. params.sim )
lens = Channel.from ( params.lens )
freqs = Channel.from ( params.freqs )
//freqs = Channel.from ( ([0.05, (1..9)*.div(10), 0.95]).flatten() )

simulation = lens.combine(freqs).combine(sims)

// Create base population
process runBasePopulation {

  label 'simulation'

  input:
  file (ini) from initialxml
  file (pars) from paramFile
  set len, freq, sim from simulation

  output:
  set len, freq, sim, file("popini.xml") into inipop

  script:
  """
  invertFREGENE -i $ini -p $pars -gn 10000 -sd $sim -recombsd $sim -o popini.xml -s -freq
  """
}

process runInversion {

  label 'simulation'
  errorStrategy 'ignore'

  publishDir "results/", pattern: '*.txt', mode: 'copy'

  input:
  file (pars) from paramFile
  set len, freq, sim, file(inipop) from inipop

  output:
  set len, freq, sim, file("l_${len}.f_${freq}.s_${sim}") into invFiles
  file("l_${len}.f_${freq}.s_${sim}_InversionSummary.txt") into invSum

  script:
  """
  invertFREGENE -i $inipop -StopFreqOfInv $freq -p $pars -gn 50000 -sd $sim \
    -recombsd $sim -StartOfInv 750000 -o l_${len}.f_${freq}.s_${sim} -EndOfInv \$(($len + 750000)) \
    -sub 1 2000 -subout 4000 -MaxFreqOfLostInv 0.1 -s -freq 
  """
}

process outputResults {

  label 'simulation'

  publishDir "results/", mode: 'copy'

  input:
  set len, freq, sim, invxml from invFiles
  file(invSum) from invSum

  output:
  set len, freq, sim, file ("l_${len}.f_${freq}.s_${sim}.haplotypes_0.dat.gz") into haplotypes
  file ("l_${len}.f_${freq}.s_${sim}.log_ids.txt") into invClass

  script:
  """
  SAMPLE -noshuffle -i $invxml -oh l_${len}.f_${freq}.s_${sim}.haplotypes -og l_${len}.f_${freq}.s_${sim}.genotypes \
    -sd $sim -controls 2000

  sed -e '1,/Position of inverted chromosomes/d' $invSum > l_${len}.f_${freq}.s_${sim}.log_ids.txt
  gzip l_${len}.f_${freq}.s_${sim}.haplotypes_0.dat
  """
}

process runRecombClust {

  label 'recombClust'

  publishDir "results/recombClustClassifications/", mode: 'copy'

  input:
  set len, freq, sim, file (haplos) from haplotypes

  output:
  file("l_${len}.f_${freq}.s_${sim}.recombClass.Rdata") into recombClass

  """
  getrecombClustClass.R $haplos l_${len}.f_${freq}.s_${sim} $len
  """

}
