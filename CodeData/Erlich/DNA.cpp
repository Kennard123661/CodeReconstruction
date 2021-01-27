//============================================================================
// Name        : DNA.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <chrono>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "Clone.hpp"
#include "Cluster2.hpp"
#include "LongestPath.hpp"
#include "EditDistance.hpp"
using namespace std;

void TestFixAll(int testNum, int strandLen, int cloneNum, int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, const int maxReps, const double delProb, const double insProb,
		const double subProb) {
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);

	int cumTotalFinalGuessEditDist = 0, roundFinalGuessEditDist = 0;
	int cumFinalGuessSubstitutions = 0, cumFinalGuessInsertions = 0, cumFinalGuessDeletions = 0;
	map<int, int> editDistanceHist;

	for (int i = 0; i < testNum; i++) {
		Cluster2 cluster(strandLen, cloneNum, delProb, insProb, subProb, generator);

//		cluster.TestFixAll(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS, roundAvgCloneEditDist,
//				roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum, round1stHalfBingoNum, round2ndHalfBingoNum,roundMinED, roundMaxED,
//				roundFinalGuessEditDist, finalGuessNum, finalStage,subPriority, delPriority, insPriority, generator, minEqual, maxReps);

//				cluster.TestFixAllClonesLast(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS,
//				roundAvgCloneEditDist, roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum,round1stHalfBingoNum, round2ndHalfBingoNum,
//				roundMinED, roundMaxED, roundFinalGuessEditDist, finalGuessNum, finalStage, subPriority, delPriority,
//				insPriority, generator, minEqual, maxReps);

//		cluster.TestFixAllClonesLastInitial(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS,
//				roundAvgCloneEditDist, roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum,round1stHalfBingoNum, round2ndHalfBingoNum,
//				roundMinED, roundMaxED, roundFinalGuessEditDist, finalGuessNum, finalStage, subPriority, delPriority,
//				insPriority, generator, minEqual, maxReps);

//		cluster.TestFixAllByErrorType(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS,
//				roundAvgCloneEditDist, roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum,round1stHalfBingoNum, round2ndHalfBingoNum,
//				roundMinED, roundMaxED, roundFinalGuessEditDist, finalGuessNum, finalStage, subPriority, delPriority,
//				insPriority, generator, minEqual, maxReps, convReps);

		string finalGuess = cluster.TestBest(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority,
				insPriority, generator, maxReps);
		roundFinalGuessEditDist = ComputeEditDistanceNum(cluster.Original(), finalGuess);
		editDistanceHist[roundFinalGuessEditDist]++;
		cumTotalFinalGuessEditDist += roundFinalGuessEditDist;

		vector<LetterOps> result = ComputeEditDistancePriority(finalGuess, cluster.Original(), 0, generator);
		map<string, double> countOperations = CountOperations(result);
		assert(countOperations["I"] + countOperations["D"] + countOperations["R"] == roundFinalGuessEditDist);

		cumFinalGuessSubstitutions += countOperations["R"];
		cumFinalGuessInsertions += countOperations["I"];
		cumFinalGuessDeletions += countOperations["D"];
	}
	map<int, int>::reverse_iterator rit = editDistanceHist.rbegin(); // points to last element in map
	int highestED = rit->first;
	int cumDist = 0;
	cout << "Edit distance hist:" << endl;
	for (int i = 0; i <= highestED; i++) {
		cumDist += editDistanceHist[i];
		cout << i << "\t" << cumDist << endl;
	}
	cout << "Avg. guess substitutions:\t" << 1000 * (double) cumFinalGuessSubstitutions / (testNum * strandLen) << endl;
	cout << "Avg. guess deletions:\t" << 1000 * (double) cumFinalGuessDeletions / (testNum * strandLen) << endl;
	cout << "Avg. guess insertions:\t" << 1000 * (double) cumFinalGuessInsertions / (testNum * strandLen) << endl;
	cout << "Avg. guess edit dist:\t" << 1000 * (double) cumTotalFinalGuessEditDist / (testNum * strandLen) << endl;
}

void TestOriginalRetention(int testNum, int strandLen, int cloneNum, int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority) {
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);

	int cumTotalFinalGuessEditDist = 0, roundFinalGuessEditDist = 0;

	for (int i = 0; i < testNum; i++) {
		Cluster2 cluster(strandLen, cloneNum, 0.1, 0.1, 0.1, generator);

		cluster.TestOriginalRetention(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority, insPriority,
				generator);

		cumTotalFinalGuessEditDist += roundFinalGuessEditDist;

	}

	cout << "Avg. original after fix edit dist:\t\t"
			<< 1000 * (double) cumTotalFinalGuessEditDist / (testNum * strandLen) << endl;
}

void TestStats(int testNum, int strandLen, int cloneNum, int subPatternLen, int delPatternLen, double insThreshold,
		const int subPriority, const int delPriority, const int insPriority) {
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);
	for (int i = 0; i < testNum; i++) {
		Cluster2 cluster(strandLen, cloneNum, 0.1, 0.1, 0.1, generator);
		//cluster.Stats(0, subPatternLen, delPatternLen, insThreshold, subPriority, delPriority, insPriority, generator);
	}
}

void GetCase(ifstream& input, string& original, vector<string>& copies, const int maxCopies) {
	string line;
	original.clear();
	copies.clear();
	getline(input, line);
	if (line.empty()) {
		return;
	}
	// first line is original
	original = line;

	// second line "*****" dump
	getline(input, line);

	// each line a copy until 2 empty lines.
	int endCase = 0;
	int copyIndex = 0;
	while (getline(input, line)) {
		if (line.empty()) {
			endCase++;
		}
		else {
			if (copyIndex < maxCopies) {
				copies.push_back(line);
				copyIndex++;
			}
		}
		if (endCase == 2) {
			break;
		}
	}
}

void TestFromFile(const string& inputFilename, int testNum, int strandLen, int maxCopies, int delPatternLen,
		const int subPriority, const int delPriority, const int insPriority, const int maxReps) {
	ifstream input;
	input.open(inputFilename.c_str());
	if (!input.is_open()) {
		cout << "Failed opening input file!" << endl;
		return;
	}

	string original;
	vector<string> copies;
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);

	int cumTotalFinalGuessEditDist = 0, roundFinalGuessEditDist = 0;
	int cumFinalGuessSubstitutions = 0, cumFinalGuessInsertions = 0, cumFinalGuessDeletions = 0;
	map<int, int> editDistanceHist;

	int countFiltered = 0;
	for (int i = 1; i <= testNum; i++) {
		GetCase(input, original, copies, maxCopies);
//		set<int>::iterator it = filter.find(i);
//		if (it != filter.end()) {
//			cout << "filtered case #" << i << endl;
//			countFiltered++;
//			continue;
//		}
		if (original.empty()) {
			break;
		}
		Cluster2 cluster(original, copies);

		string finalGuess;
		if (copies.size() == 1) { // if only 1 copy return copy.
			finalGuess = copies[0];
		}
		else {
			finalGuess = cluster.TestBest(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority, insPriority,
					generator, maxReps);
		}
		roundFinalGuessEditDist = ComputeEditDistanceNum(cluster.Original(), finalGuess);
		editDistanceHist[roundFinalGuessEditDist]++;
		cumTotalFinalGuessEditDist += roundFinalGuessEditDist;

		vector<LetterOps> result = ComputeEditDistancePriority(finalGuess, cluster.Original(), 0, generator);
		map<string, double> countOperations = CountOperations(result);
		assert(countOperations["I"] + countOperations["D"] + countOperations["R"] == roundFinalGuessEditDist);

		cumFinalGuessSubstitutions += countOperations["R"];
		cumFinalGuessInsertions += countOperations["I"];
		cumFinalGuessDeletions += countOperations["D"];
	}
	map<int, int>::reverse_iterator rit = editDistanceHist.rbegin(); // points to last element in map
	int highestED = rit->first;
	int cumDist = 0;
	cout << "Edit distance hist:" << endl;
	for (int i = 0; i <= highestED; i++) {
		cumDist += editDistanceHist[i];
		cout << i << "\t" << cumDist << endl;
	}
	cout << "Number of filtered:\t" << countFiltered << endl;
	cout << "Avg. guess substitutions:\t"
			<< 1000 * (double) cumFinalGuessSubstitutions / ((testNum - countFiltered) * strandLen) << endl;
	cout << "Avg. guess deletions:\t"
			<< 1000 * (double) cumFinalGuessDeletions / ((testNum - countFiltered) * strandLen) << endl;
	cout << "Avg. guess insertions:\t"
			<< 1000 * (double) cumFinalGuessInsertions / ((testNum - countFiltered) * strandLen) << endl;
	cout << "Avg. guess edit dist:\t"
			<< 1000 * (double) cumTotalFinalGuessEditDist / ((testNum - countFiltered) * strandLen) << endl;

	input.close();
}

void CasesFromFile(const string& inputFilename, const int maxCopies) {
	ifstream input;
	input.open(inputFilename.c_str());
	if (!input.is_open()) {
		cout << "Failed opening input file!" << endl;
		return;
	}
	int countCases = 0;
	int minCopiesNum = 1000, maxCopiesNum = 0, currentCopiesNum, cumCopiesNum = 0;
	string original;
	vector<string> copies;
	while (1) {
		GetCase(input, original, copies, maxCopies);
		if (original.empty()) {
			break;
		}
		countCases++;
		currentCopiesNum = copies.size();
		cumCopiesNum += currentCopiesNum;
		if (currentCopiesNum < minCopiesNum) {
			minCopiesNum = currentCopiesNum;
		}
		if (currentCopiesNum > maxCopiesNum) {
			maxCopiesNum = currentCopiesNum;
		}
		if (original.size() != 117 or currentCopiesNum < 6) {
			cout << countCases << "\t" << currentCopiesNum << "\t" << original.size() << endl;
		}

	}
	cout << "min copies num: " << minCopiesNum << endl;
	cout << "max copies num: " << maxCopiesNum << endl;
	cout << "average copies num: " << (double) cumCopiesNum / countCases << endl;
	cout << "cases num: " << countCases << endl;
	input.close();
}

double CalcMed(vector<int> scores) {
	size_t size = scores.size();

	if (size == 0) {
		return 0;  // Undefined, really.
	}
	else {
		sort(scores.begin(), scores.end());
		if (size % 2 == 0) {
			return (scores[size / 2 - 1] + scores[size / 2]) / 2;
		}
		else {
			return scores[size / 2];
		}
	}
}

double AVGEditDistance(const string& original, const vector<string>& copies, int& minED, int& maxED, double& medED) {
	double cumEditDistance = 0;
	int currentED;
	minED = 1000;
	maxED = 0;
	vector<int> eds;
	for (unsigned i = 0; i < copies.size(); i++) {
		currentED = ComputeEditDistanceNum(original, copies[i]);
		eds.push_back(currentED);
		cumEditDistance += currentED;
		if (currentED > maxED) {
			maxED = currentED;
		}
		if (currentED < minED) {
			minED = currentED;
		}
	}
	medED = CalcMed(eds);
	return cumEditDistance / copies.size();
}

void CasesStats(const string& inputFilename, const int maxCopies, const double EDThresh, const set<int>& filter) {
	ifstream input;
	input.open(inputFilename.c_str());
	if (!input.is_open()) {
		cout << "Failed opening input file!" << endl;
		return;
	}
	int countCases = 0;
	int currentCopiesNum, cumCopiesNum = 0;
	int minED = 0, maxED = 0;
	double medED;
	string original;
	vector<string> copies;
	while (1) {
		GetCase(input, original, copies, maxCopies);
		if (original.empty()) {
			break;
		}
		countCases++;
		currentCopiesNum = copies.size();
		cumCopiesNum += currentCopiesNum;
		double currentAVGED = AVGEditDistance(original, copies, minED, maxED, medED);

//		if (medED > EDThresh or (currentCopiesNum < 10 and currentAVGED > 10) or currentCopiesNum == 1) {
//		if (currentAVGED > 1.5 or currentCopiesNum < 10) {
		set<int>::iterator it = filter.find(countCases);
		if (it != filter.end()) {
			cout << countCases << "\t" << currentCopiesNum << "\t" << minED << "\t" << maxED << "\t" << currentAVGED
					<< "\t" << medED << endl;
		}

	}
//	cout << "min copies num: " << minCopiesNum << endl;
//	cout << "max copies num: " << maxCopiesNum << endl;
//	cout << "average copies num: " << (double) cumCopiesNum / countCases << endl;
//	cout << "cases num: " << countCases << endl;
	input.close();
}

// TODO:	decide fix order by copy len. too long -> prioritize fix inserts, too short -> prioritize fix deletions

int main() {

	clock_t begin = clock();
	int maxCopies = 25;
//	double EDThresh = 10;
//	set<int> filter = { 4490, 4756, 4850, 4879, 4896, 4929, 4937, 4940, 4950, 4955, 4969, 4971, 4973, 4977, 4978, 4979,
//			4981, 4983, 4985, 4986, 4987 };
//	CasesStats("evyaR.txt", maxCopies, EDThresh, filter); // second batch: strand length 117. casesNum 4989

	int testNum = 30;
//	int strandLen = 117;

//	int cloneNum = 2;
//	double delProb = 0.05;
//	double insProb = delProb;
//	double subProb = delProb;

	int delPatternLen = 3;

	int subPriority = 0;
	int delPriority = 0;
	int insPriority = 0;

	int maxReps = 2;
//
//
//	TestFromFile("evyaB.txt", testNum, 117, maxCopies, delPatternLen, subPriority, delPriority, insPriority,
//			maxReps);
	TestFromFile("evyaA.txt", testNum, 152, maxCopies, delPatternLen, subPriority, delPriority, insPriority,
				maxReps);

//	TestFixAll(testNum, strandLen, cloneNum, delPatternLen, subPriority, delPriority, insPriority, maxReps, delProb,
//			insProb, subProb);


	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << endl;
	cout << "Time elapsed: " << (int) elapsed_secs << "\tseconds" << endl;
	return 0;
}
