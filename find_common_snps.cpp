/*
 *
 *  Created by Aziz Mithani on 29/07/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include "Utilities.h"

using namespace std;

typedef pair<long, char> SNP;

map <char, string> NucleotideMap;

void InitialiseSNPList (string Header_File, map < string, char* >& SNPList, map < string, char* >& RefBaseList, map < string, long >& RefLength) {
	// read the header and initialise the variable to hold SNP List
	string str="";
	// open the sam file stream
	ifstream ifs ( Header_File.c_str() );
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			long length = atol(strLength.c_str());
			
			// memory initialisation
			char *RefSNPList = new char[length];
			char *RefRefBaseList = new char[length];
			for(long i = 0; i < length; i++) {
				*(RefSNPList+i) = ' ';		
				*(RefRefBaseList+i) = ' ';		
			}
			SNPList[reference] = RefSNPList;
			RefBaseList[reference] = RefRefBaseList;
			
			RefLength[reference] = length;
		} else {
			break;
		}
		
	}
	ifs.close();
	
}

void InitialiseSNPList (map < string, long >& RefLength, map < string, char* >& SNPList) {
	
	map < string, long >::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		
		string reference = itRefLength->first;
		long length = itRefLength->second;
		
		// memory initialisation
		char *RefSNPList = new char[length];
		for(long i = 0; i < length; i++) {
			*(RefSNPList+i) = ' ';		
		}
		SNPList[reference] = RefSNPList;
		
	}
	
}

void ReadSNPList (string SNP_File, map < string, char* >& SNPList, map < string, char* >& RefBaseList) {
	
	// open the snp file stream
	ifstream ifs_snp ( SNP_File.c_str() );
	
	int offset = -1;
	string strSNP = "";
	string reference;
	string PreviousReference = "";
	char* RefSNPList;
	char* RefRefBaseList;
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		reference = entries[0];
		long pos = atol(entries[1].c_str());
		string ConsensusBase = entries[3+offset];
		string ReferenceBase = entries[2+offset];
		
		if (PreviousReference.compare(reference) != 0) {
			map< string, char* >::iterator itSNPList = SNPList.find(reference);
			map< string, char* >::iterator itRefBaseList = RefBaseList.find(reference);
			// check if it is empty
			if (itSNPList !=  SNPList.end()) { // if not then get it
				RefSNPList = itSNPList->second;
			}					
			// check if it is empty
			if (itRefBaseList !=  RefBaseList.end()) { // if not then get it
				RefRefBaseList = itRefBaseList->second;
			}		
			PreviousReference = reference;
		}
		
		*(RefSNPList + pos - 1) = ConsensusBase[0];
		*(RefRefBaseList + pos - 1) = ReferenceBase[0];
		
	}
	// close the SNP file
	ifs_snp.close();
	
}
void ReadSNPList (string SNP_File, map < string, char* >& SNPList) {
	
	// open the snp file stream
	ifstream ifs_snp ( SNP_File.c_str() );
	
	int offset = -1;
	string strSNP = "";
	string reference;
	string PreviousReference = "";
	char* RefSNPList;
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		reference = entries[0];
		long pos = atol(entries[1].c_str());
		string ConsensusBase = entries[3+offset];
		//string ReferenceBase = entries[2+offset];
		
		if (PreviousReference.compare(reference) != 0) {
			map< string, char* >::iterator itSNPList = SNPList.find(reference);
			// check if it is empty
			if (itSNPList !=  SNPList.end()) { // if not then get it
				RefSNPList = itSNPList->second;
			}		
			
			PreviousReference = reference;
		}
		
		*(RefSNPList + pos - 1) = ConsensusBase[0];
		
	}
	// close the SNP file
	ifs_snp.close();
	
}

bool FindCommonSNPs (map < string, long >& RefLength, map<string, char* >& RefBaseList, map < string, char* >& SNPList1, map < string, char* >& SNPList2) {

	map < string, long >::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		string reference = itRefLength->first;
		long length = itRefLength->second;
		
		char* RefSNPList1 = SNPList1[reference];
		char* RefSNPList2 = SNPList2[reference];
		char* RefRefBaseList = RefBaseList[reference];
		
		for (int i = 0; i < length; i++) {
			char RefBase = *(RefRefBaseList + i);
			char Base1 = *(RefSNPList1 + i);
			char Base2 = *(RefSNPList2 + i);

			if (Base1 != ' ' && Base2 != ' ') {
				string strBase1 = NucleotideMap[Base1];
				string strBase2 = NucleotideMap[Base2];
				
				if (strBase1.length() > 1 && strBase2.length() > 1) { // both are ambiguous SNPs
					if (Base1 == Base2) {
						cout << reference << "\t" << i + 1 << "\t" << RefBase << "\t" << Base1 << endl;
					}
				} else if (strBase1.length() > strBase2.length()) {
					bool baseIncluded = false;
					for (int c = 0; c < strBase2.length(); c++) {
						if (strBase1.find_first_of(strBase2[c]) != string::npos) {
							baseIncluded = true;
						} else {
							baseIncluded = false;
							break;
						}
					}
					if (baseIncluded) {
						cout << reference << "\t" << i + 1 << "\t" << RefBase << "\t" << Base2 << endl;
					}
				} else if (strBase2.length() > strBase1.length()) {
					bool baseIncluded = false;
					for (int c = 0; c < strBase1.length(); c++) {
						if (strBase2.find_first_of(strBase1[c]) != string::npos) {
							baseIncluded = true;
						} else {
							baseIncluded = false;
							break;
						}
					}
					if (baseIncluded) {
						cout << reference << "\t" << i + 1 << "\t" << RefBase << "\t" << Base1 << endl;
					}
				} else if (Base1 == Base2) {
					cout << reference << "\t" << i + 1 << "\t" << RefBase << "\t" << Base1 << endl;
				}// end if strBase1.length() > 1 && strBase2.length() > 1
			} // end if Base 1 or Base 2 != ' '
		} // end for all positions
	} // end for all references

}

bool FindCommonSNPPositions (map < string, long >& RefLength, map<string, char* >& RefBaseList, map < string, char* >& SNPList1, map < string, char* >& SNPList2) {
	
	map < string, long >::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		string reference = itRefLength->first;
		long length = itRefLength->second;
		
		char* RefSNPList1 = SNPList1[reference];
		char* RefSNPList2 = SNPList2[reference];
		char* RefRefBaseList = RefBaseList[reference];
		
		for (int i = 0; i < length; i++) {
			char RefBase = *(RefRefBaseList + i);
			char Base1 = *(RefSNPList1 + i);
			char Base2 = *(RefSNPList2 + i);
			
			if (Base1 != ' ' && Base2 != ' ') {
				cout << reference << "\t" << i + 1 << "\t" << RefBase << "\t" << Base1 << "\t" << Base2 << endl;
			}
		} // end for all positions
	} // end for all references
	
}


int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	string Header_File (argv[1]);
	string SNP_File_1 (argv[2]);
	string SNP_File_2 (argv[3]);
	bool PositionsOnly = false;
	if (argc > 4) {
		if (atoi(argv[4]) == 1) {
			PositionsOnly = true;
		}
	}
	
	NucleotideMap['A'] = "A";
	NucleotideMap['C'] = "C";
	NucleotideMap['G'] = "G";
	NucleotideMap['T'] = "T";
	NucleotideMap['R'] = "AG";
	NucleotideMap['Y'] = "CT";
	NucleotideMap['M'] = "AC";
	NucleotideMap['K'] = "GT";
	NucleotideMap['S'] = "CG";
	NucleotideMap['W'] = "AT";
	NucleotideMap['H'] = "ACT";
	NucleotideMap['B'] = "CGT";
	NucleotideMap['D'] = "AGT";
	NucleotideMap['V'] = "ACG";
	NucleotideMap['N'] = "ACGT";
	
	// variable to hold SNP and other related lists 
	map < string, char* > SNPList1;
	map < string, char* > SNPList2;
	map < string, long > RefLength;
	map < string, char* > RefBaseList;

	// initialise SNP lists
	InitialiseSNPList(Header_File, SNPList1, RefBaseList, RefLength);	
	InitialiseSNPList(RefLength, SNPList2);
	
	// Read SNPs
	ReadSNPList(SNP_File_1, SNPList1, RefBaseList);
	ReadSNPList(SNP_File_2, SNPList2);
	
	if (PositionsOnly) {
		FindCommonSNPPositions(RefLength, RefBaseList, SNPList1, SNPList2);
	} else {
		FindCommonSNPs(RefLength, RefBaseList, SNPList1, SNPList2);
	}

	
	return 0;
}
