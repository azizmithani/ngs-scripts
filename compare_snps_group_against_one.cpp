/*
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
#include <set>
#include <map>
#include "Utilities.h"

using namespace std;

map <char, string> NucleotideMap;

void InitialiseSNPList(string Header_File, map < string, char** >& SNPList, map < string, int >& RefLength) {

	// open the header file stream
	ifstream ifs ( Header_File.c_str() );
	
	// read the header file and initialise the variable to hold SNPs
	string str="";
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(str, entries, "\t ");
		
			string reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			int length = atoi(strLength.c_str());
						
			// memory initialisation
			char **RefSNPList;
			RefSNPList = new char* [length];
			for(int i = 0; i < length; i++) {
				*(RefSNPList+i) = new char[2]; // first char : source SNP and second char: SNP for a line belonging to the group
				for (int j = 0; j < 2; j++) {
					(*(RefSNPList+i))[j] = ' ';
				}
			}
			
			RefLength[reference] = length;
			SNPList[reference] = RefSNPList;
		}

	}
}

void ResetSNPList(map < string, char** >& SNPList, map < string, int >& RefLength) {
	// reset the second letter (SNP for the line belonging to the group)

	map < string, int >::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		// get the SNP list for this reference
		char **RefSNPList = SNPList[itRefLength->first];
		// reset the second letter for the whole genome, i.e the SNP for the line belonging to the group
		for(int i = 0; i < itRefLength->second; i++) {
				(*(RefSNPList+i))[1] = ' ';
		}
		
		// put the reset list back
		SNPList[itRefLength->first] = RefSNPList;
	}
}

void ReadSNPList(string SNP_File, map < string, char** >& SNPList, bool isGroupSNP) {
	// open the file streams
	ifstream ifs_snp ( SNP_File.c_str() );
	
	int idx = 0;
	if (isGroupSNP) {
		idx = 1;
	}
	
	int offset = -1;
	string strSNP = "";
	string reference = "";
	string currentReference = ""; // reference (chromosome) being currently read 
	char **RefSNPList;
	// go through the list of snp file to save snps each referecne (chromosome)
	while (ifs_snp.good()) {
		
		// read the snp 
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
		
		if (reference.compare(currentReference) != 0) { // new reference
			if (currentReference.length() > 0) {
				SNPList[currentReference] = RefSNPList;
				//ofs << currentReference << "\t" << RefSNPList.size() << endl;
			}
			RefSNPList = SNPList[reference];
			currentReference = reference;
		}
		
		int SNPPos = atoi(entries[1].c_str()) - 1; // 0-based indexing
		string SNPBase = entries[3+offset];
		(*(RefSNPList+SNPPos))[idx] = SNPBase[0]; // Get the Consensus (SNP) base
		
	} // end while
	SNPList[currentReference] = RefSNPList;
	//	ofs << currentReference << "\t" << RefSNPList.size() << endl;
	
	ifs_snp.close();
	
	
}

void CheckSNPList(map < string, char** >& SNPList, map < string, int >& RefLength, string Accession, ofstream& ofs) {
	// check for common and unique SNPs
	
	// variable to count the occurrences
	int AccessionSNPCount = 0;
	int SourceSNPCount = 0;
	int CommonSharedCount = 0;
	int CommonDifferentCount = 0;
	
	map < string, int >::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		string reference = itRefLength->first;
		// get the SNP list for this reference
		char **RefSNPList = SNPList[reference];
		// go through each position
		for (int i = 0; i < itRefLength->second; i++) {

			char SourceSNP = (*(RefSNPList+i))[0];
			char GroupSNP = (*(RefSNPList+i))[1];
			if (SourceSNP == ' ' && GroupSNP == ' ') { // no SNP
				continue;
			} else if (SourceSNP == ' ') { // SNP unique to the Group Member
				ofs << Accession << "\t" << reference << "\t" << i + 1 << "\t" << "Accession" << endl;
				AccessionSNPCount++;
			} else if (GroupSNP == ' ') { // SNP unique to the Source
				ofs << Accession << "\t" << reference << "\t" << i + 1 << "\t" << "Source" << endl;
				SourceSNPCount++;
			} else {
				string strBase = NucleotideMap[SourceSNP];
				if (strBase.find_first_of(GroupSNP) == string::npos) { 
					ofs << Accession << "\t" << reference << "\t" << i + 1 << "\t" << "Common-Different" << endl;
					CommonDifferentCount++;
				} else {
					ofs << Accession << "\t" << reference << "\t" << i + 1 << "\t" << "Common-Shared" << endl;
					CommonSharedCount++;
				}
			}

		} // end for each position
		
	} // end for each reference
	
	cout << Accession << "\t" << SourceSNPCount << "\t" << AccessionSNPCount << "\t" << CommonSharedCount << "\t" << CommonDifferentCount << endl;
}


int main (int argc, char** argv) {
	
	if (argc < 5) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	string Header_File (argv[1]);
	string Source_SNP_File (argv[2]);
	string Group_SNP_List_File (argv[3]);
	string OUT_File (argv[4]);
	
	// Nucleotide mapping
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
	NucleotideMap['N'] = "N";
	
	map < string, char** > SNPList;
	map < string, int > RefLength;
	// initialise the SNP list (contains 2 columns - 1. Source SNP, 2. Group-Member SNP)
	InitialiseSNPList(Header_File, SNPList, RefLength);
	// read SNP list of the source
	ReadSNPList(Source_SNP_File, SNPList, false);
	
	// read the file containing the filenames for SNP lists for the group
	map <string, string> Group_SNP_List;
	ReadFile2(Group_SNP_List_File, Group_SNP_List);

	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		
	cout << "Accession\tSource SNPs\tAccession SNPs\tCommon-Shared SNPs\tCommon-Different SNPs" << endl;
	
	map <string, string>::iterator itGroup_SNP_List = Group_SNP_List.begin();
	for (; itGroup_SNP_List != Group_SNP_List.end(); itGroup_SNP_List++) {
		string Accession = itGroup_SNP_List->first;
		string Filename =  itGroup_SNP_List->second;

		// Reset the SNP list
		ResetSNPList(SNPList, RefLength);
		// read SNP list of the current member of the group
		ReadSNPList(Filename, SNPList, true);
		
		// check for common and unique SNPs
		CheckSNPList (SNPList, RefLength, Accession, ofs);
		
	}		
	
	// close the output file
	ofs.close(); 

	return 0;
}

