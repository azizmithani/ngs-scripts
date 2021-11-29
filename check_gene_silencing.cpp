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

const int MIN_READ_COUNT = 2;
const double MIN_READ_PERCENTAGE = 5.0;
const double GENOME_MEMBERSHIP_THRESHOLD = 0.75;

const int MIN_BASE_QUALITY = 20;
const int OPTIONAL_FIELD_START = 11;

typedef pair<long, char> SNP;

map <char, string> NucleotideMap;

class SNPGroup {
public:
	set< SNP > SNPList; // SNPs in this Group
	map < SNP, int > SNPReadCount; // No. of reads supporting a SNP
//	int ReadCount;
	bool GroupActive;
	int OverlapSize;
	
	bool AGenome;
	bool DGenome;
	double ASNPPercentage;
	double DSNPPercentage;
};

void InitialiseSNPList (string SAM_File, map < string, char* >& SNPList, map < string, long >& RefLength) {
	// read the header from SAM File and initialise the variable to hold SNP List
	string str="";
	// open the sam file stream
	ifstream ifs ( SAM_File.c_str() );
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
			for(long i = 0; i < length; i++) {
				*(RefSNPList+i) = ' ';		
			}
			
			SNPList[reference] = RefSNPList;			
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

int CountSNPsInTheGeneRegion (char* RefSNPList, long GeneStart, long GeneEnd) {
	
	int count = 0;
	for (long l = GeneStart - 1; l < GeneEnd; l++) { // we start from ReadStart - 1 because SNP list is 0-based
		if ( *(RefSNPList+l) != ' ' ) {
			count++;
		}
	}
	
	return count;
}

void CheckGeneSilencing(string GFF_File, map < string, long >& RefLengths, map < string, char* >& CSGenomeSNPList, map < string, char* >& DiploidSNPList, string OUT_File) {

	char* RefCSGenomeSNPList;
	char* RefDiploidSNPList;
	long RefLength;

	// gff file
	ifstream ifs_gff (GFF_File.c_str());
	string strGFF = "";
	string reference = "";
	string previous_reference = "";
	vector< pair< long, long > > RefGFFData;
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		long start = atol(entries[2].c_str());
		long end = atol(entries[3].c_str());
		
		if (previous_reference.compare(reference) != 0) {
			RefCSGenomeSNPList = CSGenomeSNPList[reference];
			RefDiploidSNPList = DiploidSNPList[reference];
			RefLength = RefLengths[reference];
		}
		
		long GeneStart = (start - 100 < 1 ? 1 : start - 100);
		long GeneEnd = (end + 100 > RefLength ? RefLength : end + 100);
		int CSGenomeSNPs = CountSNPsInTheGeneRegion(RefCSGenomeSNPList, GeneStart, GeneEnd);
		int DiploidSNPs = CountSNPsInTheGeneRegion(RefDiploidSNPList, GeneStart, GeneEnd);
		
		if (DiploidSNPs > 0 && CSGenomeSNPs == 0) {
			cout << reference << "\t" << start << "\t" << end << "\t" << DiploidSNPs << "\t" << CSGenomeSNPs << "\tSilenced" << endl;
		//} else if (DiploidSNPs == 0 && CSGenomeSNPs > 0) {
		//	cout << reference << "\t" << start << "\t" << end << "\t" << DiploidSNPs << "\t" << CSGenomeSNPs << "\tUnsilenced" << endl;
		}
	} // end while ifs_gff is good
	
}

int main (int argc, char** argv) {
	
	if (argc < 5) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Header_File (argv[1]);
	string GFF_File (argv[2]);
	string CS_Genome_SNP_File (argv[3]);
	string Diploid_SNP_File (argv[4]);
	string OUT_File (argv[5]);
	
	
	// variable to hold SNP and other related lists 
	map < string, vector< pair<long, long> > > GFFData;

	map < string, char* > CSGenomeSNPList;
	map < string, char* > DiploidSNPList;
	map < string, long > RefLengths;
	
	// initialise SNP lists
	InitialiseSNPList(Header_File, CSGenomeSNPList, RefLengths);	
	InitialiseSNPList(RefLengths, DiploidSNPList);
	
	// Read SNPs
	ReadSNPList(CS_Genome_SNP_File, CSGenomeSNPList);
	ReadSNPList(Diploid_SNP_File, DiploidSNPList);
	
	// process the sam file
	CheckGeneSilencing(GFF_File, RefLengths, CSGenomeSNPList, DiploidSNPList, OUT_File);
	
	return 0;
}
