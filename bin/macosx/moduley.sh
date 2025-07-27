#!/bin/bash
### change this to !/bin/zsh for zsh 

# 
# moduley.sh 
# 
# this script makes lists of CModules and Boundary Sets 
# 
# generic shell script
# this should work for any shell environment (zsh, bash, etc)
# 
# 

# RUN fLPS using a range of parameters designed 
# to coverage diverse target lengths and data set coverages 
./fLPS2/fLPS2 -o long -m 7  -M 11  -t 5.2e-05  -O $2.moduley  $1 

./fLPS2/fLPS2 -o long -m 12 -M 16  -t 5.4e-06  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 10 -M 20  -t 1.8e-05  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 9  -M 30  -t 6.9e-04  -O $2.moduley  $1 

./fLPS2/fLPS2 -o long -m 21 -M 25  -t 6.2e-09  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 23 -M 33  -t 4.1e-07  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 20 -M 50  -t 4.7e-05  -O $2.moduley  $1 

./fLPS2/fLPS2 -o long -m 32 -M 36  -t 7.9e-14  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 38 -M 48  -t 7.3e-10  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 34 -M 74  -t 5.2e-07  -O $2.moduley  $1 

./fLPS2/fLPS2 -o long -m 68 -M 78  -t 4.1e-18  -O $2.moduley  $1 
./fLPS2/fLPS2 -o long -m 78 -M 128 -t 2.1e-11  -O $2.moduley  $1 


# CONCATENATE fLPS OUTPUT 
ls -1a | grep $2.moduley | grep fLPS.out | awk '{print "cat\t" $1;}' | sh | sort -k1,1 -k8,8g | uniq | grep -v WHOLE > $2.all.moduley.temp.txt 

# MAKE CMODULES & BOUNDARY_SETS
./moduley    $2.all.moduley.temp.txt > $2.moduley.out.txt
./moduley -f $2.all.moduley.temp.txt > $2.Cmodules.out.fasta 

rm -f $2*moduley*fLPS.out $2.all.moduley.temp.txt 

