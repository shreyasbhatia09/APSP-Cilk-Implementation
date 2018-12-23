#!/bin/bash
export  CILK_NWORKERS=16
export executable_PATH="./"

#---------------------------------------------------------------------------------
#commands for I-DP

#get MEM values for I-DP
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./idp 4096 256 2 >> idpMEM1.1.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./idp 4096 256 2  >> idpMEM2.1.csv & sudo likwid-perfctr  -c 2,3 -g MEM -M 1 -O ./idp 4096 256 2 >>idpMEM2.2.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./idp 4096 256 2  >> idpMEM3.1.csv & sudo likwid-perfctr  -c 2,3 -g MEM -M 1 -O  ./idp 4096 256 2 >>idpMEM3.2.csv & sudo likwid-perfctr  -c 4-5 -g MEM -M 1 -O ./idp 4096 256 2 >>idpMEM3.3.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./idp 4096 256 2  >> idpMEM4.1.csv & sudo likwid-perfctr  -c 2,3 -g MEM -M 1 -O ./idp 4096 256 2 >>idpMEM4.2.csv & sudo likwid-perfctr  -c 4-5 -g MEM -M 1 -O ./idp 4096 256 2 >>idpMEM4.3.csv & sudo likwid-perfctr  -c 6-7 -g MEM -M 1 -O ./idp 4096 256 2 >>idpMEM4.4.csv
#---------------------------------------------------------------------------------

#copmmands for R-DP

#get MEM values for R-DP
sudo likwid-perfctr  -c 0-1 -g MEM -M 1 -O ./rdp 4096 256   2 >> rdp1.1.csv
sudo likwid-perfctr  -c 0-1 -g MEM -M 1 -O ./rdp 4096 256   2 >> rdp2.1.csv &  sudo likwid-perfctr  -c 2-3 -g MEM -M 1 -O ./rdp 4096 256   2 >>rdp2.2.csv
sudo likwid-perfctr  -c 0-1 -g MEM -M 1 -O  ./rdp 4096 256   2  >> rdp3.1.csv & sudo likwid-perfctr  -c 2-3 -g MEM -M 1 -O ./rdp 4096 256   2 >>rdp3.2.csv &  sudo likwid-perfctr  -c 4-5 -g MEM -M 1 -O ./rdp 4096 256   2 >>rdp3.3.csv
sudo likwid-perfctr  -c 0-1 -g MEM -M 1 -O  ./rdp 4096 256   2  >> rdp4.1.csv & sudo likwid-perfctr  -c 2-3 -g MEM -M 1 -O ./rdp 4096 256   2 >>rdp4.2.csv &  sudo likwid-perfctr  -c 4-5 -g MEM -M 1 -O ./rdp 4096 256   2 >>rdp4.3.csv &  sudo likwid-perfctr  -c 6-7 -g MEM -M 1 -O ./rdp 4096 256   2 >>rdp4.4.csv



#---------------------------------------------------------------------------------

#copmmands for T-I-DP

#get MEM values for T-I-DP
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./tidp 4096 1024  2  >> tidpMEM1.1.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./tidp 4096 1024  2  >> tidpMEM2.1.csv & sudo likwid-perfctr  -c 2,3 -g MEM -M 1 -O ./tidp 4096 1024  2 >>tidpMEM2.2.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O  ./tidp 4096 1024  2  >> tidpMEM3.1.csv & sudo likwid-perfctr  -c 2,3 -g MEM -M 1 -O ./tidp 4096 1024  2 >>tidpMEM3.2.csv & sudo likwid-perfctr  -c 4-5 -g MEM -M 1 -O  ./tidp 4096 1024  2 >>tidpMEM3.3.csv
sudo CILK_NWORKERS=2  likwid-perfctr  -c 0,1 -g MEM -M 1 -O ./tidp 4096 1024  2  >> tidpMEM4.1.csv & sudo likwid-perfctr  -c 2,3 -g MEM -M 1 -O ./tidp 4096 1024  2 >>tidpMEM4.2.csv & sudo likwid-perfctr  -c 4-5 -g MEM -M 1 -O ./tidp 4096 1024  2 >>tidpMEM4.3.csv & sudo likwid-perfctr  -c 6-7 -g MEM -M 1 -O ./tidp 4096 1024  2 >>tidpMEM4.4.csv

