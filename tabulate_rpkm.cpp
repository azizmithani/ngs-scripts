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
#include "Utilities.h"

using namespace std;

typedef map < pair<int, int>, double* > REF_GFF_DATA;
typedef map < string, REF_GFF_DATA > GFF_DATA;
typedef map<int, string > REF_GFF_GENE_NAME;
typedef map<string, REF_GFF_GENE_NAME > GFF_GENE_NAME;

int GetReadEndCoordinates (int ReadStart, string cigar) {
	
	long sum = 0;
	string length = "";
	for (int i = 0; i < cigar.length(); i++) {
		if (isdigit(cigar[i])) {
			length += cigar[i];
		} else if (cigar[i]=='S' || cigar[i]=='H') {
			// ignore
			length="";
		} else if (cigar[i] =='D') {
			sum+= atol(length.c_str());
			length="";
		} else if (cigar[i] =='I') {
			length = "";
		} else if (cigar[i]=='M' || cigar[i]=='N') {
			sum += atol(length.c_str());
			length="";
		}
	} // end for
	
	return ReadStart + sum - 1;
}

int ProcessSAMFile(string SAM_File, GFF_DATA& GFFData, int idxAcc){
	
	// open the sam file stream
	ifstream ifs_sam ( SAM_File.c_str() );
	
	string strSAM = "";
	string Read = "";
	string reference = "";
	string previous_reference = "";
	int TotalReads = 0;
	
	REF_GFF_DATA RefGFFData;
	REF_GFF_DATA::iterator itRefGFFData;
	
	while (ifs_sam.good()) {
		getline(ifs_sam, strSAM);

		if (strSAM.empty()) { // ignore emtry lines
			continue;
		} else if (strSAM[0] == '@') { // ignore the header
			continue;
		}
		
		bool invalidGene = false;
		// get individual fields
		vector<string> entries;
		strsplit(strSAM, entries, "\t");
		
		// get the reference names
		reference = entries[2];
		
		// read coordinates
		int ReadStart = atoi(entries[3].c_str());
		int ReadEnd = GetReadEndCoordinates(ReadStart, entries[5]);
		
		// get the GFF data for this reference
		if (reference.compare(previous_reference) != 0) {
			
			if (previous_reference.length() > 0) {
				// save the GFF data for the previous reference
				GFFData[previous_reference] = RefGFFData;
			}			
			
			// get the new data
			RefGFFData = GFFData[reference];
			itRefGFFData = RefGFFData.begin();
			
			while (itRefGFFData != RefGFFData.end() && ReadStart > itRefGFFData->first.second + 20) {
				itRefGFFData++;
			}
			
			if (itRefGFFData == RefGFFData.end() || ReadEnd < itRefGFFData->first.first - 20) {
				invalidGene = true;
			}
			
			previous_reference = reference;
		}		
		
		if (ReadStart > itRefGFFData->first.second) {
			// read lies after the current gene			
			
			while (itRefGFFData != RefGFFData.end() && ReadStart > itRefGFFData->first.second + 20) {
				itRefGFFData++;
			}
			
			if (itRefGFFData == RefGFFData.end() || ReadEnd < itRefGFFData->first.first - 20) {
				invalidGene = true;
			}
			
		}
		
		if (!invalidGene) 
			itRefGFFData->second[idxAcc]++;
		
		TotalReads++;
		
	} // end if sam file is good
	
	// save the GFF data for the last reference
	GFFData[previous_reference] = RefGFFData;
	
	ifs_sam.close();
	
	return TotalReads;
}

void TabulateRPKM (map <string, string>& Accessions,  GFF_DATA& RPKMTable) {
	
	cout << "Calculating and tabulating RPKM ..." << endl;
	
	map <string, string>::iterator itAccession = Accessions.begin();
	int idxAcc = -1;
	for (; itAccession != Accessions.end(); itAccession++) {
		string Accession = itAccession->first;
		string Filename =  itAccession->second;
		idxAcc++;
		
		cout << "\t" << Accession << endl;
		
		if (!FileExists(Filename) ) { // The file does not exist, move to the next one
			continue;
		}

		// calculate the number of reads per gene
		int TotalReads = ProcessSAMFile(Filename, RPKMTable, idxAcc);

		// calculate the RPKM
		GFF_DATA::iterator itRPKMTable = RPKMTable.begin();
		for (; itRPKMTable != RPKMTable.end(); itRPKMTable++) {
			REF_GFF_DATA RefRPKMTable = itRPKMTable->second;
			
			REF_GFF_DATA::iterator itRefRPKMTable = RefRPKMTable.begin();
			for (; itRefRPKMTable != RefRPKMTable.end(); itRefRPKMTable++) {
				int GeneLength = itRefRPKMTable->first.second - itRefRPKMTable->first.first + 1;
				
				// RPKM = 10^9 C / NL (C: number of mappable reads in the gene/exons, N: total number of reads, L: gene/exons length) 
				itRefRPKMTable->second[idxAcc] = (1000000000.00 / (double)(TotalReads)) * (itRefRPKMTable->second[idxAcc] / (double)(GeneLength));
				//cout << itRefGFFData->second.first << "\t" << GeneLength << "\t" << itRefGFFData->second.second << "\t" << RPKM << endl;
			}
			itRPKMTable->second = RefRPKMTable;
		}
		
	} // end for each accession

	
}

void WriteTable (string OUT_File, map <string, string>& Accessions, GFF_DATA& RPKMTable, GFF_GENE_NAME& GFFGeneName) {
	
	cout << "Writing output ..." << endl;

	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		

	ofs << "Chromosome\tUnigene\tStart\tEnd";
	map <string, string>::iterator itAccession = Accessions.begin();
	for (; itAccession != Accessions.end(); itAccession++) {
		string Accession = itAccession->first;
		ofs << "\t" << Accession;
	}
	ofs << endl;

	GFF_DATA::iterator itRPKMTable = RPKMTable.begin();
	for (; itRPKMTable != RPKMTable.end(); itRPKMTable++) {
		REF_GFF_DATA RefRPKMTable = itRPKMTable->second;
		REF_GFF_GENE_NAME RefGFFGeneName = GFFGeneName[itRPKMTable->first];
		
		REF_GFF_DATA::iterator itRefRPKMTable = RefRPKMTable.begin();
		for (; itRefRPKMTable != RefRPKMTable.end(); itRefRPKMTable++) {
			ofs << itRPKMTable->first << "\t" << RefGFFGeneName[itRefRPKMTable->first.first] << "\t" << itRefRPKMTable->first.first << "\t" << itRefRPKMTable->first.second;
			
			double* RPKM = itRefRPKMTable->second;
			for (int i = 0; i < Accessions.size(); i++) {
				ofs << "\t" << RPKM[i];
			}
			ofs << endl;
		} // end for each RPKM entry for this reference
	} // end for each reference
	
	ofs.close();
}


void ReadGFFFile (string GFF_File, GFF_DATA& RPKMTable, GFF_GENE_NAME& GFFGeneName, int NoOfAccessions) {
	
	REF_GFF_DATA RefGFFData;
	REF_GFF_GENE_NAME RefGFFGeneName;
	// gff file
	ifstream ifs_gff (GFF_File.c_str());
	string strGFF = "";
	string reference = "";
	string previous_reference = "";
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		} else if (strGFF[0] == '#') {
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		// get start and end coordinates
		long start = atol(entries[3].c_str());
		long end = atol(entries[4].c_str());
		
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				RPKMTable[previous_reference] = RefGFFData;
				GFFGeneName[previous_reference] = RefGFFGeneName;
			}
			
			RefGFFData = RPKMTable[reference];
			RefGFFGeneName = GFFGeneName[reference];
			
			previous_reference = reference;
		}
		
		double* RPKM = new double[NoOfAccessions];
		for(int i = 0; i < NoOfAccessions; i++) {
			RPKM[i]= 0.0;
		}
		
		RefGFFData[pair<long, long>(start,end)] = RPKM;
		int StartPos = entries[8].find("ID=") + 3;
		int EndPos = entries[8].find_first_of(";", StartPos);
		RefGFFGeneName[start] = entries[8].substr(StartPos, EndPos - StartPos);
		
	} // end while ifs_gff is good
	RPKMTable[previous_reference] = RefGFFData;
	GFFGeneName[previous_reference] = RefGFFGeneName;
	
	ifs_gff.close();
	
}


int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Accession_File (argv[1]);
	string GFF_File (argv[2]);	
	string OUT_File (argv[3]);	
	
	map <string, string> Accessions;
	ReadFile2(Accession_File, Accessions);

	// variable to store GFF & RPKM Data
	GFF_DATA RPKMTable;
	GFF_GENE_NAME GFFGeneName;
	ReadGFFFile(GFF_File, RPKMTable, GFFGeneName, Accessions.size());
	
	TabulateRPKM (Accessions, RPKMTable);

	
	WriteTable(OUT_File, Accessions, RPKMTable, GFFGeneName);
	
	return 0;
}
