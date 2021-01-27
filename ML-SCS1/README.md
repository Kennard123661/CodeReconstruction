# Reconstruction Algorithms

Source code of the ML-SCS reconstruction algorithm. This algorithm works for small cluster of the deletion channel. 
 

## Compilation

  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o SCS2.o SCS2.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o SCSN.o SCSN.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Strings.o Strings.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o MultiString.o MultiString.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o NChooseK.o NChooseK.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o KSCS.o KSCS.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Tests.o Tests.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o IndexVector.o IndexVector.cpp
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Deletions.o Deletions.cpp

  g++ -o Del *.o


## Usage

./Del >results.txt

The output file presents Levenshtein distance histogram. The i-th entry of the histogram shows the number of clusters that their estimated sequence has edit distance of i errors from the original strings. 


## 

## License
TBA
