#!/bin/bash
### change this to !/bin/zsh for zsh 

# 
# bandy.sh 
# 
# this script picks out compositional banding
# 
# generic shell script
# this should work for any shell environment (zsh, bash, etc)
# 
# 

# RUN fLPS 
./fLPS2/fLPS2 -dbo long -m 3  -M 5  -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 3  -M 10 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 3  -M 15 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 3  -M 20 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 5  -M 10 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 5  -M 15 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 5  -M 20 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 10 -M 15 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 10 -M 20 -t 1e-04  -O $2  $1 
./fLPS2/fLPS2 -dbo long -m 15 -M 15 -t 1e-04  -O $2  $1  
./fLPS2/fLPS2 -dbo long -m 15 -M 20 -t 1e-04  -O $2  $1  
./fLPS2/fLPS2 -dbo long -m 20 -M 20 -t 1e-04  -O $2  $1  

# CONCATENATE fLPS OUTPUT AND RUN BANDY 
ls -1a | grep $2 | grep fLPS.out | awk '{print "./bandy " $1;}' | sh | sort -k 2,3 | uniq > $2.all.bandy.temp.txt

# PICK BANDING 
./pickbanding $2.all.bandy.temp.txt > $2.bandy.out.txt 

rm -f $2*.fLPS.out $2.all.bandy.temp.txt 

