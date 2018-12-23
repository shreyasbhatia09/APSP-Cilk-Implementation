#!/bin/bash
export  CILK_NWORKERS=16
export executable_PATH="./"

N=8192
p=16
./tidp $N 128 $p >> tidp_time.txt
./rdp $N 128 $p >> rdp_time.txt
./idp $N 128 $p >> idp_time.txt

./tidp_papi $N 128 $p >> tidp_cachemiss.txt
./rdp_papi $N 128 $p >> rdp_cachemiss.txt
./idp_papi $N 128 $p >> idp_cachemiss.txt

