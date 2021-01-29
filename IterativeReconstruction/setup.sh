rm *.o
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LCS2.o LCS2.cpp -fopenmp
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o EditDistance.o EditDistance.cpp -fopenmp
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Clone.o Clone.cpp -fopenmp
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Cluster2.o Cluster2.cpp -fopenmp
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LongestPath.o LongestPath.cpp -fopenmp
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o CommonSubstring2.o CommonSubstring2.cpp  -fopenmp
g++ -std=c++17 -O3 -g3 -Wall -c -fmessage-length=0 -o DNA.o DNA.cpp -fopenmp
g++ -o DNA *.o -fopenmp
