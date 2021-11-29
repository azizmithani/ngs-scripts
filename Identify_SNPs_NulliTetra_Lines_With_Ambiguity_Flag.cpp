/*
 *
 *  Created by Aziz Mithani on 07/05/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <limits>
#include "Utilities.h"

using namespace std;

const int MIN_COVERAGE = 5;
const int MIN_BASE_COVERAGE = 3;
const double MIN_BASE_PROPORTION = 0.05;
const int OFFSET = 3;

typedef map < pair<string, int>, vector <string> > SNP_TABLE; //[<reference, position>] = list of snps
typedef map < pair<string, int>, int > COVERAGE_LIST; //[<reference, position>] = coverage
typedef map < pair<string, int>, int* > BASE_DIST_LIST; //[<reference, position>] = <,,,,> 

void TabulateSNPPositions (string Position_File, map <string, string>& NT_SNP_Files, SNP_TABLE& NT_SNPs) {
	
	cout << "Tabulating SNP Positions ..." << endl;
	
	ifstream ifs( Position_File.c_str() );
	if (!ifs) { // The file does not exist, move to the next one
		return;
	}
	
	string str = "";
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		if (entries[2].compare("N") == 0) { // ignore positions with "N" in the reference
			continue;
		}
		
		// create the key (chromosome, position)
		pair<string, int> key (entries[0], atoi(entries[1].c_str()));
		
		// find the key in the table
		SNP_TABLE::iterator itNT_SNPs = NT_SNPs.find(key);
		
		vector <string> SNPRow(NT_SNP_Files.size() + OFFSET);
		if (itNT_SNPs != NT_SNPs.end()) {
			SNPRow = itNT_SNPs->second;
		} else {
			SNPRow[0] = entries[0];
			SNPRow[1] = entries[1];
			SNPRow[2] = entries[2];
		}
		
		if (itNT_SNPs != NT_SNPs.end()) {
			itNT_SNPs->second = SNPRow;
		} else {
			NT_SNPs[key] = SNPRow;
		}
		
	} // end while ifs.good()
	
	// close the file
	ifs.close();
	
	
}

void TabulateSNPs (map <string, string>& NT_SNP_Files, SNP_TABLE& NT_SNPs) {
	
	cout << "Tabulating SNPs ..." << endl;
	
	map <string, string>::iterator itNT_SNP_Files = NT_SNP_Files.begin();
	int idxNT = -1;
	for (; itNT_SNP_Files != NT_SNP_Files.end(); itNT_SNP_Files++) {
		string NT_Line = itNT_SNP_Files->first;
		string Filename =  itNT_SNP_Files->second;
		idxNT++;
		
		cout << "\t" << NT_Line << endl;
		
		ifstream ifs( Filename.c_str() );
		if (!ifs) { // The file does not exist, move to the next one
			continue;
		}
		
		string str = "";
		while (ifs.good()) {
			getline(ifs, str);
			
			if (str.empty()) { // ignore emtry lines
				continue;
			}
			
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			if (entries[2].compare("N") == 0) { // ignore positions with "N" in the reference
				continue;
			}
			
			// create the key (chromosome, position)
			pair<string, int> key (entries[0], atoi(entries[1].c_str()));
			
			// find the key in the table
			SNP_TABLE::iterator itNT_SNPs = NT_SNPs.find(key);
			
			vector <string> SNPRow(NT_SNP_Files.size() + OFFSET);
			if (itNT_SNPs != NT_SNPs.end()) {
				SNPRow = itNT_SNPs->second;
			} else {
				SNPRow[0] = entries[0];
				SNPRow[1] = entries[1];
				SNPRow[2] = entries[2];
			}
			
			SNPRow[idxNT + OFFSET] = entries[3];
			
			if (itNT_SNPs != NT_SNPs.end()) {
				itNT_SNPs->second = SNPRow;
			} else {
				NT_SNPs[key] = SNPRow;
			}
			
		} // end while ifs.good()
		
		// close the file
		ifs.close();
		
	} // end for each NT Line
	
	
}

void ReadCoverageFile (string Coverage_File, COVERAGE_LIST& CoverageList) {
	// open the Coverage file stream
	ifstream ifs_Coverage ( Coverage_File.c_str() );
	string str="";
	while (ifs_Coverage.good()) {
		getline(ifs_Coverage, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		// create the key (chromosome, position)
		pair<string, int> key (entries[0], atoi(entries[1].c_str()));
		CoverageList[ key ] = atoi(entries[2].c_str());
		
	}
	ifs_Coverage.close();
}

void CheckCoverages(map <string, string>& CoverageFiles, map <string, string>& NT_SNP_Files, SNP_TABLE& NT_SNPs){
	
	cout << "Checking coverages ..." << endl;
	
	map <string, string>::iterator itNT_SNP_Files = NT_SNP_Files.begin();
	int idxNT = -1;
	for (; itNT_SNP_Files != NT_SNP_Files.end(); itNT_SNP_Files++) {
		set< pair <string, int> > SNPsToBeDeleted;
		
		idxNT++;
		string NTLine = itNT_SNP_Files->first;
		cout << "\t" << NTLine << endl;
		
		// Get the Coverage Filename for this NT Line
		map <string, string>::iterator itCoverageFiles = CoverageFiles.find(NTLine);		
		string Filename =  itCoverageFiles->second;
		// Read the Coverage file
		COVERAGE_LIST CoverageList;
		ReadCoverageFile (Filename, CoverageList);
		
		SNP_TABLE::iterator itNT_SNPs = NT_SNPs.begin();
		for (; itNT_SNPs != NT_SNPs.end(); itNT_SNPs++) {
			vector<string> SNPRow = itNT_SNPs->second;
			
			// create the key (chromosome, position)
			pair<string, int> key (SNPRow[0], atoi(SNPRow[1].c_str()));
			
			// check the coverage
			if (CoverageList[key] < MIN_COVERAGE) {
				SNPsToBeDeleted.insert(itNT_SNPs->first);
				continue;
			}
		} // for each SNP	
		
		set< pair <string, int> >::iterator itSNPsToBeDeleted = SNPsToBeDeleted.begin();
		for (; itSNPsToBeDeleted != SNPsToBeDeleted.end(); itSNPsToBeDeleted++) {
			NT_SNPs.erase(NT_SNPs.find(*itSNPsToBeDeleted));		
		}
		
	} // for each line
}

int CalculateCoverage(int* BaseDist) {
	
	int Coverage = 0;
	for (int b = 0; b < TOTAL_BASES; b++) {
		Coverage += BaseDist[b];
	}
	
	return Coverage;
}

void ReadBaseDistributionFile(string Base_Dist_File, BASE_DIST_LIST& BaseDistList) {
	// open the file stream
	ifstream ifs_base_dist ( Base_Dist_File.c_str() );
	string str="";
	while (ifs_base_dist.good()) {
		getline(ifs_base_dist, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		// create the key (chromosome, position)
		pair<string, int> key (entries[0], atoi(entries[1].c_str()));
		
		int* BaseFreq = new int[TOTAL_BASES];
		for ( int i = 3; i < entries.size(); i++ ) {
			BaseFreq[i-3] = atoi(entries[i].c_str());
		}
		
		BaseDistList[ key ] = BaseFreq;
		
	}
	ifs_base_dist.close();
}

void CheckBaseDistribution(map <string, string>& BaseDistributionsFiles, map <string, string>& NT_SNP_Files, SNP_TABLE& NT_SNPs){
	
	cout << "Checking base distributions ..." << endl;
	
	// Initialise the nucleotide map
	map <string, set<string> > NucleotideMap;
	InitialiseNucleotideMap(NucleotideMap);
	map <string, string > NucleotideMapReverse;
	InitialiseNucleotideMapReverse(NucleotideMapReverse);
	
	map <string, string>::iterator itNT_SNP_Files = NT_SNP_Files.begin();
	int idxNT = -1;
	for (; itNT_SNP_Files != NT_SNP_Files.end(); itNT_SNP_Files++) {
		set< pair <string, int> > SNPsToBeDeleted;
		
		idxNT++;
		string NT_Line = itNT_SNP_Files->first;
		cout << "\t" << NT_Line << endl;
		
		map <string, string>::iterator itBaseDistributionsFiles = BaseDistributionsFiles.find(NT_Line);		
		string Filename =  itBaseDistributionsFiles->second;
		// Read the base distribution file
		BASE_DIST_LIST BaseDistList;
		ReadBaseDistributionFile (Filename, BaseDistList);
		
		SNP_TABLE::iterator itNT_SNPs = NT_SNPs.begin();
		for (; itNT_SNPs != NT_SNPs.end(); itNT_SNPs++) {
			vector<string> SNPRow = itNT_SNPs->second;
			
			
			string Base = (SNPRow[idxNT + OFFSET].length() == 0 ? SNPRow[2] : SNPRow[idxNT + OFFSET]);
			set<string> AllBases = NucleotideMap[Base];
			set<string> NewBases;
			
			// check if a base is present at this position but not called as a SNP (ambiguous or certain) by samtools because of coverage
			// create the key (chromosome, position)
			pair<string, int> key (SNPRow[0], atoi(SNPRow[1].c_str()));
			// get base distributions for this position
			int* BaseDist = BaseDistList[key];
			// calculate coverage for this position
			double Coverage = (double)CalculateCoverage(BaseDist);
			
			// check the coverage
			if (Coverage < MIN_COVERAGE) {
				//SNPsToBeDeleted.insert(itNT_SNPs->first);
				itNT_SNPs->second[idxNT + OFFSET] = "<";
				continue;
			}
			
			NewBases.clear();
			// check for all bases (a,c,g,t,n)
			for (int b = 0; b < TOTAL_BASES; b++) {
				if (BaseDist[b] >= MIN_BASE_COVERAGE && (double)BaseDist[b]/Coverage >= MIN_BASE_PROPORTION) {
					// this base is present in sufficient quantity at this position, add it to the base if not present
					string theBase (1, GetIndexBaseBaseDistFile(b));
					NewBases.insert(theBase);
					//if (SNPRow[1].compare("41") == 0) {
					//	cout << SNPRow[1] << "\t" << Base << "\t" << theBase << endl;
					//	cout << theBase << " added at position " << SNPRow[1] << endl;
					//}
				}
			}
			// Get the new consensus base
			if (NewBases.size() > 0) {
				Base = NucleotideMapReverse[SetToString(NewBases)];
			}
			
			itNT_SNPs->second[idxNT + OFFSET] = Base;
			
		} // for each SNP
		
		set< pair <string, int> >::iterator itSNPsToBeDeleted = SNPsToBeDeleted.begin();
		for (; itSNPsToBeDeleted != SNPsToBeDeleted.end(); itSNPsToBeDeleted++) {
			NT_SNPs.erase(NT_SNPs.find(*itSNPsToBeDeleted));		
		}
		
	} // for each line
	
}

string FindGenomeBase(set<string>& AllBases, set<string>& Nulli, set<string>& NonNulli1, set<string>& NonNulli2 ) {
	
	set<string>::iterator itBase = AllBases.begin();
	for (; itBase != AllBases.end(); itBase++) {
		if (Nulli.find(*itBase) == Nulli.end()) { // base not found in nulli line
			if (NonNulli1.find(*itBase) != NonNulli1.end() && NonNulli2.find(*itBase) != NonNulli2.end()) { //but found in both other non-nulli lines
				return *itBase;
			}
		}
	}
	return "";
	
}

void AssignSNPs(map <string, string>& NT_SNP_Files, SNP_TABLE& NT_SNPs, string OUT_File){
	
	cout << "Assigning SNPs ..." << endl;
	
	// Initialise the nucleotide map
	map <string, set<string> > NucleotideMap;
	InitialiseNucleotideMap(NucleotideMap);
	
	ofstream ofs ( OUT_File.c_str() );
	ofs << "Chromosome\tPosition\tRef Base";
	map <string, string>::iterator itNT_SNP_Files = NT_SNP_Files.begin();
	for (; itNT_SNP_Files != NT_SNP_Files.end(); itNT_SNP_Files++) {
		ofs << "\t" << itNT_SNP_Files->first;
	} // for each line
	ofs << endl;
	
	SNP_TABLE::iterator itNT_SNPs = NT_SNPs.begin();
	for (; itNT_SNPs != NT_SNPs.end(); itNT_SNPs++) {
		vector<string> SNPRow = itNT_SNPs->second;
		
		// get the consensus base for each of the 6 nulli-tetra line
		string baseNXATXB = (SNPRow[OFFSET + 0].length() == 0 ? SNPRow[2] : SNPRow[OFFSET + 0]); 
		string baseNXATXD = (SNPRow[OFFSET + 1].length() == 0 ? SNPRow[2] : SNPRow[OFFSET + 1]); 
		string baseNXBTXA = (SNPRow[OFFSET + 2].length() == 0 ? SNPRow[2] : SNPRow[OFFSET + 2]); 
		string baseNXBTXD = (SNPRow[OFFSET + 3].length() == 0 ? SNPRow[2] : SNPRow[OFFSET + 3]); 
		string baseNXDTXA = (SNPRow[OFFSET + 4].length() == 0 ? SNPRow[2] : SNPRow[OFFSET + 4]); 
		string baseNXDTXB = (SNPRow[OFFSET + 5].length() == 0 ? SNPRow[2] : SNPRow[OFFSET + 5]); 
		
		// if two nulli lines of the same genome don't thav the same base, ignore the position
		// it will also ignore positions where one out of two nulli lines of the same genome has low coverage
		if (baseNXATXB.compare(baseNXATXD) != 0) {
			ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
			continue;
		} else if (baseNXBTXA.compare(baseNXBTXD) != 0) {
			ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
			continue;
		} else if (baseNXDTXA.compare(baseNXDTXB) != 0) {
			ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
			continue;
		}
		
		// check if any genome has low coverage for both nulli lines
		if (baseNXATXB[0] == '<' || baseNXBTXA[0] == '<' || baseNXDTXA[0] == '<') {
			ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
			continue;
		}
		
		// get the corresponding nucleotides (a,c,g,t,n) for each nulli line
		set<string> NXA = NucleotideMap[baseNXATXB];
		set<string> NXB = NucleotideMap[baseNXBTXA];
		set<string> NXD = NucleotideMap[baseNXDTXA];
		
		// get the union of all nucleotides
		set<string> AllBases;
		AllBases.insert(NXA.begin(), NXA.end());
		AllBases.insert(NXB.begin(), NXB.end());
		AllBases.insert(NXD.begin(), NXD.end());
		
		if (AllBases.size() > 3) { // more than 3 bases, must be an error
			ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
			continue;
		}
		
		//Find the base out of all bases found at the position which is missing (if any) from each nulli line 
		string BaseA = FindGenomeBase(AllBases, NXA, NXB, NXD);
		string BaseB = FindGenomeBase(AllBases, NXB, NXA, NXD);
		string BaseD = FindGenomeBase(AllBases, NXD, NXA, NXB);
		
		
		int TotalLength = BaseA.length() + BaseB.length() + BaseD.length();
		if ( TotalLength == 0 ) {
			if (NXA.size() == 1 && NXB.size() == 1 && NXD.size() == 1 && AllBases.size() == 1) { // common snp
				BaseA = *AllBases.begin();
				BaseB = BaseA;
				BaseD = BaseA;
			} else {
				ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
				continue;
			}
			
		} else if ( TotalLength == 1 ) {
			if (AllBases.size() != 2) { // extra bases present, cannot assign properly
				if (AllBases.size() == 1) { 
					BaseA = *AllBases.begin();
					BaseB = BaseA;
					BaseD = BaseA;
				} else {
					ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t*\t*\t*" << endl; 
					continue;
				}
			} else {
				if (BaseA.length() > 0) {
					AllBases.erase(AllBases.find(BaseA));
					BaseB = *AllBases.begin();
					BaseD = BaseB;
				} else if (BaseB.length() > 0) {
					AllBases.erase(AllBases.find(BaseB));
					BaseA = *AllBases.begin();
					BaseD = BaseA;
				} else if (BaseD.length() > 0) {
					AllBases.erase(AllBases.find(BaseD));
					BaseA = *AllBases.begin();
					BaseB = BaseA;
				}
			}
			
		} else if (TotalLength == 2) {
			// One genome is being silenced
			//cout << "Two genomes assigned at " << SNPRow[1] << endl;
			BaseA = (BaseA.length() == 0 ? "." : BaseA);
			BaseB = (BaseB.length() == 0 ? "." : BaseB);
			BaseD = (BaseD.length() == 0 ? "." : BaseD);
		} else {
			// bases were assigned for all three genomes
			// do nothing
		}
		
		ofs << SNPRow[0] << "\t" << SNPRow[1] << "\t" << SNPRow[2] << "\t" << BaseA << "\t" << BaseB << "\t" << BaseD << endl; 
		
	} // for each SNP
	
	ofs.close();
}

void WriteTable (string OUT_File, map <string, string>& NT_SNP_Files, SNP_TABLE& NT_SNPs) {
	
	cout << "Writing output ..." << endl;
	
	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		
	
	// header
	ofs << "Chromosome\tPosition\tRef Base";
	map <string, string>::iterator itNT_SNP_Files = NT_SNP_Files.begin();
	for (; itNT_SNP_Files != NT_SNP_Files.end(); itNT_SNP_Files++) {
		ofs << "\t" << itNT_SNP_Files->first;
	}
	ofs << endl;
	
	SNP_TABLE::iterator itNT_SNPs = NT_SNPs.begin();
	for (; itNT_SNPs != NT_SNPs.end(); itNT_SNPs++) {
		vector<string> SNPRow = itNT_SNPs->second;
		
		vector<string>::iterator itSNPRow = SNPRow.begin();
		ofs << *itSNPRow;
		itSNPRow++;
		for (; itSNPRow != SNPRow.end(); itSNPRow++) {
			ofs << "\t" << ((*itSNPRow).length() == 0 ? SNPRow[2] : *itSNPRow);
		}
		ofs << endl;
	}
	ofs.close();
}

void RemoveCommonSNPs (SNP_TABLE& NT_SNPs) {
	
	cout << "Removing common SNPs ..." << endl;
	
	vector< pair <string, int> > SNPsToBeDeleted;
	SNP_TABLE::iterator itNT_SNPs = NT_SNPs.begin();
	for (; itNT_SNPs != NT_SNPs.end(); itNT_SNPs++) {
		vector<string> SNPRow = itNT_SNPs->second;
		
		int idxCol = -1;
		string FirstBase = "";
		bool isCommonSNP = true;
		vector<string>::iterator itSNPRow = SNPRow.begin();
		for (; itSNPRow != SNPRow.end(); itSNPRow++) {
			idxCol++;
			
			if (idxCol == OFFSET) {
				FirstBase = ((*itSNPRow).length() == 0 ? SNPRow[2] : *itSNPRow);
			} else if (idxCol > OFFSET) {
				string Base = ((*itSNPRow).length() == 0 ? SNPRow[2] : *itSNPRow);
				if (FirstBase[0] != Base[0]) {
					isCommonSNP = false;
					break;
				}
			}
		}	// end for each entry of the variation row
		
		if (isCommonSNP) { // common SNP -> delete it
			SNPsToBeDeleted.push_back(itNT_SNPs->first);
		}
		
	} // for each row of the variation table
	
	vector< pair <string, int> >::iterator itSNPsToBeDeleted = SNPsToBeDeleted.begin();
	for (; itSNPsToBeDeleted != SNPsToBeDeleted.end(); itSNPsToBeDeleted++) {
		NT_SNPs.erase(NT_SNPs.find(*itSNPsToBeDeleted));		
	}
}

int main (int argc, char** argv) {
	
	string File_NT_SNP_Files (argv[1]);
	string File_Base_Dist_Files =  (argv[2]);
	string OUT_File (argv[3]);
	string Position_File (argv[4]);
	
	map <string, string> NT_SNP_Files;
	ReadFile2(File_NT_SNP_Files, NT_SNP_Files);
	
	SNP_TABLE NT_SNPs;
	if (!Position_File.empty()) {	
		TabulateSNPPositions (Position_File, NT_SNP_Files, NT_SNPs);
	}
	TabulateSNPs (NT_SNP_Files, NT_SNPs);
	
	map <string, string> Base_Dist_Files;
	ReadFile2(File_Base_Dist_Files, Base_Dist_Files);
	CheckBaseDistribution(Base_Dist_Files, NT_SNP_Files, NT_SNPs);			
	
	//CheckCoverages(Coverage_Files, NT_SNP_Files, NT_SNPs);			
	
	// remove SNPs that common across all lines (possibly sequencing error/varietal difference in the reference)
	//RemoveCommonSNPs(NT_SNPs);
	
	AssignSNPs(NT_SNP_Files, NT_SNPs, OUT_File);
	
	WriteTable(OUT_File + ".tab", NT_SNP_Files, NT_SNPs);
	
	return 0;
}
