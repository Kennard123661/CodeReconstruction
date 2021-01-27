# Reconstruction Algorithms

Source code of the BMA reconstruction algorithm. This algorithm works for clusters of the deletion channel. 
 

## Compilation

  g++ -std=c++0x -o BMA BMA.cpp




## Usage

./BMA >results.txt

The output file presents Levenshtein distance histogram. The i-th entry of the histogram shows the number of clusters that their estimated sequence has edit distance of i errors from the original strings. 


## 

## License
TBA
