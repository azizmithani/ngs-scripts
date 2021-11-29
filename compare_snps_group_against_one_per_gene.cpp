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

void ReadGFFFile(string GFF_File, map < string, map < pair<int, int>, string > >& GeneIDs, map < string, map < pair<int, int>, int* > >& GeneSNPCount, int nAccessions) {
	
	int nCol = nAccessions + 1;
	ifstream ifs_gff ( GFF_File.c_str() );
	string strGFF = "";
	
	string reference = "";
	string previous_reference = "";
	map < pair<int, int>, string > RefGeneIDs;
	map < pair<int, int>, int* > RefGeneSNPCount;
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		}
		if (strGFF[0] == '#') { // ignore header
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		int start = atoi(entries[3].c_str());
		int end = atoi(entries[4].c_str());
		string ID = entries[8].substr(entries[8].find("ID=",0) + 3, entries[8].find_first_of(';') - 3);
		
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GeneIDs[previous_reference] = RefGeneIDs;
				GeneSNPCount[previous_reference] = RefGeneSNPCount;
			}
			
			RefGeneSNPCount = GeneSNPCount[reference];
			RefGeneIDs = GeneIDs[reference];

			previous_reference = reference;
		}
		
		RefGeneIDs[pair<int, int>(start,end)] = ID;
		int* Counts = new int[nCol];
		for (int i = 0; i < nCol; i++) {
			Counts[i] = 0;
		}
		RefGeneSNPCount[pair<int, int>(start,end)] = Counts;
		//cout << ID << endl;
	} // end while
	// update the last reference
	GeneIDs[reference] = RefGeneIDs;
	GeneSNPCount[reference] = RefGeneSNPCount;
	
}

void InitialiseSNPList(string Header_File, map < string, char** >& SNPList, map < string, int >& RefLength, int nAccessions) {

	int nCol = nAccessions + 1;
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
				*(RefSNPList+i) = new char[nCol]; // first char : source SNP and remaining chars: SNP for lines belonging to the group
				for (int j = 0; j < nCol; j++) {
					(*(RefSNPList+i))[j] = ' ';
				}
			}
			
			RefLength[reference] = length;
			SNPList[reference] = RefSNPList;
		}
	}
}

void ReadSNPList(string SNP_File, map < string, char** >& SNPList, int idx) {
	// open the file streams
	ifstream ifs_snp ( SNP_File.c_str() );
		
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

void CheckSNPList(map < string, char** >& SNPList, map < string, int >& RefLength, map <string, string>& Group_SNP_List, map < string, map < pair<int, int>, string > >& GeneIDs, map < string, map < pair<int, int>, int* > >& GeneSNPCount, ofstream& ofs) {

	// number of accessions
	int nAccessions = Group_SNP_List.size();
	
	ofs << "Reference\tGene ID\tStart\tEnd\tSource";
	map <string, string>::iterator itGroup_SNP_List = Group_SNP_List.begin();
	for (; itGroup_SNP_List != Group_SNP_List.end(); itGroup_SNP_List++) {
		ofs << "\t" << itGroup_SNP_List->first;
	}
	ofs << endl;
	
	// check SNPs
	map < string, map < pair<int, int>, int* > >::iterator itGeneSNPCount = GeneSNPCount.begin();
	for (; itGeneSNPCount != GeneSNPCount.end(); itGeneSNPCount++) {

		string reference = itGeneSNPCount->first;
		// get the SNP list for this reference
		char **RefSNPList = SNPList[reference];

		// get the GeneSNPCount list for this reference
		map < pair<int, int>, int* > RefGeneSNPCount = itGeneSNPCount->second;
		map < pair<int, int>, int* >::iterator itRefGeneSNPCount = RefGeneSNPCount.begin();
		// and the itertor for the GeneIds for this reference
		map < pair<int, int>, string >::iterator itRefGeneIDs = GeneIDs[reference].begin();
		for (; itRefGeneSNPCount != RefGeneSNPCount.end(); itRefGeneSNPCount++, itRefGeneIDs++) {
			int start = itRefGeneSNPCount->first.first;
			int end = itRefGeneSNPCount->first.second;
			
			for (int i = start - 1; i < end; i++) {
				
				char SourceSNP = (*(RefSNPList+i))[0];
				if (SourceSNP == ' ') { // SNP not present in the source, ignore this SNP
					continue;
				} else {
					int AbsentAccCount = 0;
					for (int idx = 1; idx <= nAccessions; idx++) {
						if ((*(RefSNPList+i))[idx] == ' ') { // No SNP in this accession
							// increment the count of accessions not having this SNP
							AbsentAccCount++;
						}
					}
					if (AbsentAccCount == nAccessions) { // SNP not present in any of the accessions, ignore this SNP
						continue;
					}
				}
				
				for (int idx = 0; idx <= nAccessions; idx++) {
					if ((*(RefSNPList+i))[idx] != ' ') {					
						itRefGeneSNPCount->second[idx]++;
					}
				}
			}
			
			ofs << reference << "\t" << itRefGeneIDs->second << "\t" << start << "\t" << end;
			for (int idx = 0; idx <= nAccessions; idx++) {
				ofs << "\t" << itRefGeneSNPCount->second[idx];
			}
			ofs << endl;
		}
		
	}
	
/*		
	cout << "SNPs common to all members of the group but absent in the Source: " << ExtraAndCommonSNPs << endl;
	cout << "SNPs common to all members of the group and shared with the Source: " << SharedAndCommonSNPs << endl;
	cout << "SNPs present in the source but absent in all members of the group: " << MissingInAllAccessionsSNPs << endl;
	cout << "Total Positions with SNP in the Source and at least one members of the group: " << TotalPositionsSNPInSourceAndAtleastOneAccession << endl; 
	cout << "Total Positions with SNP not in the Source but at least one members of the group: " << TotalPositionsNoSNPInSource << endl; 
	cout << "Total Positions with SNP in the Source: " << TotalPositionsSNPInSource << endl; 
	cout << endl;
	
	map <string, string>::iterator itGroup_SNP_List = Group_SNP_List.begin();
	for (; itGroup_SNP_List != Group_SNP_List.end(); itGroup_SNP_List++) {
		cout << "\t" << itGroup_SNP_List->first;
	}
	cout << endl << "SNPs absent in the source but present in the accession";
	for (int i = 0; i < nAccessions; i++) {
		cout << "\t" << ExtraInGroupAccession[i];
	}
	cout << endl << "SNPs shared between the source and the accession";
	for (int i = 0; i < nAccessions; i++) {
		cout << "\t" << SharedInGroupAccession[i];
	}
	cout << endl << "SNPs present in the source but absent in the accession";
	for (int i = 0; i < nAccessions; i++) {
		cout << "\t" << AbsentInGroupAccession[i];
	}
	cout << endl << "SNP position shared but different bases in the source and the accession";
	for (int i = 0; i < nAccessions; i++) {
		cout << "\t" << SharedButDifferentInGroupAccession[i];
	}
	
	cout << endl;
*/
}


int main (int argc, char** argv) {
	
	if (argc < 6) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	string Header_File (argv[1]);
	string Source_SNP_File (argv[2]);
	string Group_SNP_List_File (argv[3]);
	string GFF_File (argv[4]);
	string OUT_File (argv[5]);
	
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

	// read the file containing the filenames for SNP lists for the group
	map <string, string> Group_SNP_List;
	ReadFile2(Group_SNP_List_File, Group_SNP_List);
	
	map < string, char** > SNPList;
	map < string, int > RefLength;
	// initialise the SNP list 
	InitialiseSNPList(Header_File, SNPList, RefLength, Group_SNP_List.size());
	// read SNP list of the source
	ReadSNPList(Source_SNP_File, SNPList, 0);
	// read SNP lists for each member of the group
	int idx = 0;
	map <string, string>::iterator itGroup_SNP_List = Group_SNP_List.begin();
	for (; itGroup_SNP_List != Group_SNP_List.end(); itGroup_SNP_List++) {
		string Accession = itGroup_SNP_List->first;
		string Filename =  itGroup_SNP_List->second;
		
		// read SNP list of the current member of the group
		ReadSNPList(Filename, SNPList, ++idx);
		
	}
	
	map < string, map < pair<int, int>, string > > GeneIDs;
	map < string, map < pair<int, int>, int* > > GeneSNPCount;
	ReadGFFFile(GFF_File, GeneIDs, GeneSNPCount, Group_SNP_List.size());


	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		

	// check for common and unique SNPs
	CheckSNPList (SNPList, RefLength, Group_SNP_List, GeneIDs, GeneSNPCount, ofs);
	
//	cout << "Accession\tSource SNPs\tAccession SNPs\tCommon-Shared SNPs\tCommon-Different SNPs" << endl;
	
	
	// close the output file
	ofs.close(); 

	return 0;
}

