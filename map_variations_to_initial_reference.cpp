/*
 *  fix_mutations.cpp
 *  
 *	fixes the supplied mutation in the given FASTA file
 *
 *  Created by Aziz Mithani on 15/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "Utilities.h"

using namespace std;

const int BLOCK_SIZE = 80;

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

void MapPositions(string Variation_File, map < string, int* >& PositionMap) {
	
	// open file streams
	ifstream ifs ( Variation_File.c_str() );
	ofstream ofs ( (Variation_File + ".map").c_str() );
	
	// go through the list of snp file to save snps each referecne (chromosome)
	string str = "";
	string reference = "";
	string currentReference = ""; // reference (chromosome) being currently read 
	int* RefPositionMap;
	while (ifs.good()) {
		
		// read the variation
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		reference = entries[0];
		int position = atoi(entries[1].c_str());
		
		if (reference.compare(currentReference) != 0) { // new reference
			RefPositionMap = PositionMap[reference];
			currentReference = reference;
		}
		
		ofs << str << "\t" << *(RefPositionMap + position) << endl;		
		
	} // end while
	
	ifs.close();
	ofs.close();
	
	
}

void CalculateRefLengths(string FASTA_File, map < string, int >& RefLength, map< string, int* >& PositionMap) {
	
	ifstream ifs( FASTA_File.c_str() );
	string str = "";
	string reference = "";
	int length = 0;
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str[0] == '>') { // we are reading a new reference (chromosome)
			if (reference.length() > 0) {
				int *RefPositionMap = new int[length + 1];
				for(int i = 0; i <= length; i++) {
					*(RefPositionMap+i) = i;
				}
				
				RefLength[reference] = length;
				PositionMap[reference] = RefPositionMap;
			}
			
			reference = ExtractReferenceName(str);
			length = 0;
		} else {
			length += str.length();			
		}
	} // end while
	
	// process the last reference
	int *RefPositionMap = new int[length + 1];
	for(int i = 0; i <= length; i++) {
		*(RefPositionMap+i) = i;
	}
	
	RefLength[reference] = length;
	PositionMap[reference] = RefPositionMap;
	
	ifs.close();
	
}

void ReadPositionMapFile(string Map_File, map< string, int* >& PositionMap) {
	
	ifstream ifs( Map_File.c_str() );
	
	string reference = "";
	string previous_reference = "";
	int *RefPositionMap;
	string str = "";
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else {
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[0];
			
			if (reference.compare(previous_reference) != 0) {
				if (previous_reference.length() > 0) {
					PositionMap[previous_reference] = RefPositionMap;
				}
				
				RefPositionMap = PositionMap[reference];
				
				previous_reference = reference;
			}
			
			int CurrentPosition = atoi(entries[1].c_str());
			int InitialPosition = atoi(entries[2].c_str());
			
			*(RefPositionMap+CurrentPosition) = InitialPosition;
			
		}
	} // end while
	PositionMap[reference] = RefPositionMap;
	
	ifs.close();
	
}

int main (int argc, char** argv) {
	
	if (argc < 5) {
		cout <<  "Files not specified" << endl;
		exit(0);
	}
	string FASTA_File (argv[1]);
	string SNP_File (argv[2]);
	string Deletion_File (argv[3]);
	string Insertion_File (argv[4]);
	

	map < string, int > RefLength;
	map < string, int* > PositionMap;
	CalculateRefLengths(FASTA_File, RefLength, PositionMap);
	
	string Map_File = FASTA_File + ".map";
	if (FileExists(Map_File)) {
		ReadPositionMapFile(Map_File, PositionMap);
	}
	
	if (SNP_File.length() > 0) {
		MapPositions(SNP_File, PositionMap);
	}
	if (Deletion_File.length() > 0) {
		MapPositions(Deletion_File, PositionMap);
	}
	if (Insertion_File.length() > 0) {
		MapPositions(Insertion_File, PositionMap);
	}
	
	
	return 0;
}
