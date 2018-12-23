#!/bin/bash
export  CILK_NWORKERS=16
export executable_PATH="./"

#---------------------------------------------------------------------------------
#commands for I-DP

#get runtimes for I-DP

numactl -C 0-1 ./idp 4096 256  2 >> idp1.1.txt
numactl -C 0-1 ./idp 4096 256  2 >> idp2.1.txt &  numactl -C 2-3 ./idp 4096 256  2 >>idp2.2.txt
numactl -C 0-1 ./idp 4096 256  2 >> idp3.1.txt &  numactl -C 2-3 ./idp 4096 256  2 >>idp3.2.txt &  numactl -C 4-5 ./idp 4096 256  2  >>idp3.3.txt
numactl -C 0-1 ./idp 4096 256  2 >> idp4.1.txt &  numactl -C 2-3 ./idp 4096 256  2 >>idp4.2.txt &  numactl -C 4-5 ./idp 4096 256  2  >>idp4.3.txt & numactl -C 6-7 ./idp 4096 256  2 >>idp4.4.txt

#get cachemisses for I-DP
numactl -C 0-1 ./idp_papi 4096 256  2 >> idp1.1.txt
numactl -C 0-1 ./idp_papi 4096 256  2 >> idp2.1.txt &  numactl -C 2-3 ./idp_papi 4096 256  2 >>idp2.2.txt
numactl -C 0-1 ./idp_papi 4096 256  2 >> idp3.1.txt &  numactl -C 2-3 ./idp_papi 4096 256  2  >>idp3.2.txt &  numactl -C 4-5 ./idp_papi 4096 256  2  >>idp3.3.txt
numactl -C 0-1 ./idp_papi 4096 256  2 >> idp4.1.txt &  numactl -C 2-3 ./idp_papi 4096 256  2 >>idp4.2.txt &  numactl -C 4-5 ./idp_papi 4096 256  2  >>idp4.3.txt & numactl -C 6-7 ./idp_papi 4096 256  2 >>idp4.4.txt

#get energy values for I-DP
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./idp 4096 256  2 >> idpENERGY1.1.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./idp 4096 256  2  >> idpENERGY2.1.csv & sudo likwid-perfctr  -c 2,3 -g ENERGY -M 1 -O ./idp 4096 256  2 >>idpENERGY2.2.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./idp 4096 256  2  >> idpENERGY3.1.csv & sudo likwid-perfctr  -c 2,3 -g ENERGY -M 1 -O  ./idp 4096 256  2 >>idpENERGY3.2.csv & sudo likwid-perfctr  -c 4-5 -g ENERGY -M 1 -O ./idp 4096 256  2 >>idpENERGY3.3.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./idp 4096 256  2  >> idpENERGY4.1.csv & sudo likwid-perfctr  -c 2,3 -g ENERGY -M 1 -O ./idp 4096 256  2 >>idpENERGY4.2.csv & sudo likwid-perfctr  -c 4-5 -g ENERGY -M 1 -O ./idp 4096 256  2 >>idpENERGY4.3.csv & sudo likwid-perfctr  -c 6-7 -g ENERGY -M 1 -O ./idp 4096 256  2 >>idpENERGY4.4.csv
#---------------------------------------------------------------------------------

#copmmands for R-DP

#get runtimes for R-DP
numactl -C 0-1 ./rdp 4096 256  2 >> rdp1.1.txt
numactl -C 0-1 ./rdp 4096 256  2 >> rdp2.1.txt &  numactl -C 2-3 ./rdp 4096 256  2 >>rdp2.2.txt
numactl -C 0-1 ./rdp 4096 256  2 >> rdp3.1.txt &  numactl -C 2-3 ./rdp 4096 256  2 >>rdp3.2.txt &  numactl -C 4-5 ./rdp 4096 256  2  >>rdp3.3.txt
numactl -C 0-1 ./rdp 4096 256  2 >> rdp4.1.txt &  numactl -C 2-3 ./rdp 4096 256  2 >>rdp4.2.txt &  numactl -C 4-5 ./rdp 4096 256  2  >>rdp4.3.txt & numactl -C 6-7 ./rdp 4096 256  2 >>rdp4.4.txt
 
#get cachemisses for R-DP
numactl -C 0-1 ./rdp_papi 4096 256  2 >> rdp1.1.txt
numactl -C 0-1 ./rdp_papi 4096 256  2 >> rdp2.1.txt &  numactl -C 2-3 ./rdp_papi 4096 256  2 >>rdp2.2.txt
numactl -C 0-1 ./rdp_papi 4096 256  2 >> rdp3.1.txt &  numactl -C 2-3 ./rdp_papi 4096 256  2  >>rdp3.2.txt &  numactl -C 4-5 ./rdp_papi 4096 256  2  >>rdp3.3.txt
numactl -C 0-1 ./rdp_papi 4096 256  2 >> rdp4.1.txt &  numactl -C 2-3 ./rdp_papi 4096 256  2 >>rdp4.2.txt &  numactl -C 4-5 ./rdp_papi 4096 256  2  >>rdp4.3.txt & numactl -C 6-7 ./rdp_papi 4096 256  2 >>rdp4.4.txt


#get energy values for R-DP
sudo likwid-perfctr  -c 0-1 -g ENERGY -M 1 -O ./rdp 4096 256  2 >> rdp1.1.csv
sudo likwid-perfctr  -c 0-1 -g ENERGY -M 1 -O ./rdp 4096 256  2 >> rdp2.1.csv &  sudo likwid-perfctr  -c 2-3 -g ENERGY -M 1 -O ./rdp 4096 256  2 >>rdp2.2.csv
sudo likwid-perfctr  -c 0-1 -g ENERGY -M 1 -O  ./rdp 4096 256  2  >> rdp3.1.csv & sudo likwid-perfctr  -c 2-3 -g ENERGY -M 1 -O ./rdp 4096 256  2 >>rdp3.2.csv &  sudo likwid-perfctr  -c 4-5 -g ENERGY -M 1 -O ./rdp 4096 256  2 >>rdp3.3.csv
sudo likwid-perfctr  -c 0-1 -g ENERGY -M 1 -O  ./rdp 4096 256  2  >> rdp4.1.csv & sudo likwid-perfctr  -c 2-3 -g ENERGY -M 1 -O ./rdp 4096 256  2 >>rdp4.2.csv &  sudo likwid-perfctr  -c 4-5 -g ENERGY -M 1 -O ./rdp 4096 256  2 >>rdp4.3.csv &  sudo likwid-perfctr  -c 6-7 -g ENERGY -M 1 -O ./rdp 4096 256  2 >>rdp4.4.csv



#---------------------------------------------------------------------------------

#copmmands for T-I-DP

#get runtimes for T-I-DP
numactl -C 0-1 ./tidp 4096 1024  2 >> tidp1.1.txt
numactl -C 0-1 ./tidp 4096 1024  2 >> tidp2.1.txt &  numactl -C 2-3 ./tidp 4096 1024  2 >> tidp2.2.txt
numactl -C 0-1 ./tidp 4096 1024  2 >> tidp3.1.txt &  numactl -C 2-3 ./tidp 4096 1024  2  >> tidp3.2.txt &  numactl -C 4-5 ./tidp 4096 1024  2  >> tidp3.3.txt
numactl -C 0-1 ./tidp 4096 1024  2 >> tidp4.1.txt &  numactl -C 2-3 ./tidp 4096 1024  2 >> tidp4.2.txt &  numactl -C 4-5 ./tidp 4096 1024  2  >> tidp4.3.txt & numactl -C 6-7 ./tidp 4096 1024  2  >> tidp4.4.txt


#get cachemisses for T-I-DP
numactl -C 0-1 ./tidp_papi 4096 1024  2 >> tidp1.1.txt
numactl -C 0-1 ./tidp_papi 4096 1024  2 >> tidp2.1.txt &  numactl -C 2-3 ./tidp_papi 4096 1024  2 >> tidp2.2.txt
numactl -C 0-1 ./tidp_papi 4096 1024  2 >> tidp3.1.txt &  numactl -C 2-3 ./tidp_papi 4096 1024  2  >> tidp3.2.txt &  numactl -C 4-5 ./tidp_papi 4096 1024  2  >> tidp3.3.txt
numactl -C 0-1 ./tidp_papi 4096 1024  2 >> tidp4.1.txt &  numactl -C 2-3 ./tidp_papi 4096 1024  2 >> tidp4.2.txt &  numactl -C 4-5 ./tidp_papi 4096 1024  2  >> tidp4.3.txt & numactl -C 6-7 ./tidp_papi 4096 1024  2 >> tidp4.4.txt


#get energy values for T-I-DP
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./tidp 4096 1024  2  >> tidpENERGY1.1.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./tidp 4096 1024  2  >> tidpENERGY2.1.csv & sudo likwid-perfctr  -c 2,3 -g ENERGY -M 1 -O ./tidp 4096 1024  2 >>tidpENERGY2.2.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O  ./tidp 4096 1024  2  >> tidpENERGY3.1.csv & sudo likwid-perfctr  -c 2,3 -g ENERGY -M 1 -O ./tidp 4096 1024  2 >>tidpENERGY3.2.csv & sudo likwid-perfctr  -c 4-5 -g ENERGY -M 1 -O  ./tidp 4096 1024  2 >>tidpENERGY3.3.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g ENERGY -M 1 -O ./tidp 4096 1024  2  >> tidpENERGY4.1.csv & sudo likwid-perfctr  -c 2,3 -g ENERGY -M 1 -O ./tidp 4096 1024  2 >>tidpENERGY4.2.csv & sudo likwid-perfctr  -c 4-5 -g ENERGY -M 1 -O ./tidp 4096 1024  2 >>tidpENERGY4.3.csv & sudo likwid-perfctr  -c 6-7 -g ENERGY -M 1 -O ./tidp 4096 1024  2 >>tidpENERGY4.4.csv

