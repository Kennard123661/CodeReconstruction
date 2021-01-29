rm *.o
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LCS2.o LCS2.cpp -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o EditDistance.o EditDistance.cpp -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Clone.o Clone.cpp -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Cluster2.o Cluster2.cpp -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LongestPath.o LongestPath.cpp -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o CommonSubstring2.o CommonSubstring2.cpp  -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -std=c++17 -O3 -g3 -Wall -c -fmessage-length=0 -o DNA.o DNA.cpp -fopenmp -DBOOST_NO_CXX11_SCOPED_ENUMS
g++ -o DNA *.o -fopenmp
