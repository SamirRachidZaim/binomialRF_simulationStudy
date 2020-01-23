#!/bin/bash
counter =0
for dimX in  10000 20000 
  do
      echo $dimX
	for numGenesSeeded in 100 250 500
	  do 
		echo numGenesSeeded
	((counter++))
        echo "working on "$counter
	echo "#!/bin/bash
  #PBS  -N binomialRF"$counter"
  #PBS  -m bea
  #PBS  -M samirrachidzaim@email.arizona.edu
  #PBS  -W group_list=yves
  #PBS  -q 'standard'
  #PBS -l select=1:ncpus=28:mem=180gb
  #PBS -l walltime=100:00:00
  #PBS -l cput=1000:00:00

  module load R/3.5.1
  export R_LIBS=/rsgrps/yves/samir/R_libs 

  cd /rsgrps/yves/samir/binomialRF_study/code

  Rscript simulation_main_effects.R $dimX $numGenesSeeded
  " > "rf."$counter".sh"
                chmod a+rx "rf."$counter".sh"
		  qsub < "rf."$counter".sh"
	  done
done
