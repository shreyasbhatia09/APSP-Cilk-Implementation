PAPI=/opt/apps/papi/5.3.0/x86_64
make:
	icc -O3 -o rdp R-DP.cpp -DCO -AVX -xhost -ipo -parallel -I/usr/include/x86_64-linux-gnu/c++/4.8 
	icc -O3 -o tidp Tiled-I-DP.cpp -AVX -xhost -ipo -parallel -Dunroll -I/usr/include/x86_64-linux-gnu/c++/4.8 
	icc -O3 -o idp R-DP.cpp  -DLOOPDP -AVX -xhost -ipo -parallel -I/usr/include/x86_64-linux-gnu/c++/4.8
	icc -O3 -o idp_papi R-DP.cpp  -DLOOPDP -AVX -xhost -ip -parallel -DUSE_PAPI -lpapi -I/usr/include/x86_64-linux-gnu/c++/4.8 -I $(PAPI)/include -L $(PAPI)/lib
	icc -O3 -o rdp_papi R-DP.cpp -DCO -AVX -xhost -DUSE_PAPI -lpapi -ip -I/usr/include/x86_64-linux-gnu/c++/4.8 -I $(PAPI)/include -L $(PAPI)/lib
	icc -O3 -o tidp_papi Tiled-I-DP.cpp -AVX -xhost -ip -parallel -DUSE_PAPI -lpapi -unroll -I/usr/include/x86_64-linux-gnu/c++/4.8 -I $(PAPI)/include -L $(PAPI)/lib
	