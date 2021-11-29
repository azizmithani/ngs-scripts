/*
 *  process_sam_file.cpp
 *  
 *
 *  Created by Aziz Mithani on 25/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

class Gene {
public:
	int gIdx;
	string GeneName;
	string RefName;
	int Start;
	int End;
	string Strand;
	int Phase;
	int nReadsWithin;
	int nReadsToTheLeft;
	int nReadsToTheRight;
	
	Gene();
};

Gene::Gene() {
	gIdx = -1;
	GeneName = "";
	RefName = "";
	Start = -1;
	End = -1;
	Strand = "";
	Phase = -1;
	nReadsWithin = 0;
	nReadsToTheLeft = 0;
	nReadsToTheRight = 0;
}

const int MQUAL_THRESHOLD = 20;
const int UpdateInterval = 10000;

void strsplit(const string& str, vector<string>& tokens, const string& delimiters = " ") {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
	
    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

int ReadGenes(string GFF_File, map< string, vector< Gene > > *GenesList) {
	
	// open the gff file stream
	ifstream ifs_gff ( GFF_File.c_str() );
	
	int gCount = 0;
	string strGene;
	string GeneName;
	vector<string> entries;
	vector<string> genes;
	string RefName = "";
	string PrevRefName = ""; 
	vector< Gene > RefGenesList;
	while (ifs_gff.good()) {
		getline(ifs_gff, strGene);
		
		entries.clear();
		strsplit(strGene, entries, "\t");
		
		// read the next line if it's not a gene entry
		if (entries[2].compare("gene") != 0) {
			continue;
		}
		
		// get the gene name
		int pos, posColon;
		// find the occurence of the "ID"
		if ((pos = entries[8].find("ID=")) != string::npos) {
			// if found, find the first ";" after it. Everything in between is the id/gene name
			if ((posColon = entries[8].find(";",pos+1)) != string::npos) {
				GeneName = entries[8].substr(pos+3,posColon-3);
			} else { // no semicolon is found then take the entry upto the end of the string
				GeneName = entries[8].substr(pos+3);
			}
			
		} else { // ID not found. Generate a gene name using reference (chromosome) name, start and end positions 
			GeneName = entries[0] + "_" + entries[3] + "_" + entries[3];
		}
		
		Gene gene;
		gene.gIdx = gCount;
		gene.GeneName = GeneName;
		gene.RefName = entries[0];
		gene.Start = atoi(entries[3].c_str());
		gene.End = atoi(entries[4].c_str());
		gene.Strand = entries[6];
		gene.Phase = atoi(entries[7].c_str());
		
		// Get the reference (chromosome) name
		RefName = entries[0];
		
		// store the gene entry in a list		
		if (RefName.compare(PrevRefName) != 0) { // get the gene list for this reference (chromosome), if we already dont have it
			
			if (PrevRefName.compare("") != 0) {
				// put the list for the previous reference in the list of genes
				GenesList->insert( pair<string, vector< Gene > >(PrevRefName, RefGenesList) );
				RefGenesList.clear();
			}
			
			map< string, vector< Gene > >::iterator it = GenesList->find(RefName);
			// check if it is empty
			if (it !=  GenesList->end()) { // if not then get it
				RefGenesList = it->second;				
			} else { // otherwise initialise it
				RefGenesList.clear();
			}
		} 
		// add the gene entry to the 
		RefGenesList.push_back(gene);
		
		// increment the gene count
		gCount++;
		// save the reference name for the next iteration
		PrevRefName = RefName;
	}
	// put the gene list for the last reference name into the list of genes 
	GenesList->insert( pair< string, vector< Gene > >(RefName, RefGenesList) );
	// close the gff file
	ifs_gff.close();
	
/*
	//========================= DEBUG =======================
	 map< string, vector< Gene > >::iterator it = GenesList->begin();
	 while (it != GenesList->end()) {
		 cout << it->first << ": "  << endl;
		 RefGenesList = it->second;
		 vector< Gene >::iterator itGene = RefGenesList.begin();
		 while (itGene != RefGenesList.end()) {
			 cout << "\t" << itGene->GeneName << "\t" << itGene->Start << endl;
			 itGene++;
		 }
		 it++;
	 }
	 //======================================================
*/
	
	return gCount;
}

int GetReadEndCoordinates(int ReadStart, string cigar) {

	vector<string> SegmentLengths;
	strsplit(cigar, SegmentLengths, "MN");
	
	// sum the segment lengths to get the end coordinate
	int sum = 0;
	vector<string>::iterator it = SegmentLengths.begin();
	while (it != SegmentLengths.end()) {
		sum += atoi(it->c_str());
		it++;
	}
	
	return ReadStart + sum - 1;
}

bool UpdateReadsCount(vector< Gene > *RefGenesList, int ReadStart, int ReadEnd) {
/**
	Returns true if the read was mapped either within, to the left or to the right of any gene
 */
	
	bool ReadAssociatedWithAGene = false;
	// find the genes that are associated with this read (read either lies completely within the gene, or extend past one border)
	vector< Gene >::iterator itGene = RefGenesList->begin();
	while (itGene != RefGenesList->end()) {
				
		if (itGene->Start <= ReadStart && itGene->End >= ReadEnd) { // within
			itGene->nReadsWithin++;
			ReadAssociatedWithAGene = true;
		} else if (itGene->Start > ReadStart && itGene->Start <= ReadEnd) { // left
			itGene->nReadsToTheLeft++;
			ReadAssociatedWithAGene = true;
		} else if (itGene->End >= ReadStart && itGene->End < ReadEnd) { // right
			itGene->nReadsToTheRight++;
			ReadAssociatedWithAGene = true;
		}
		
		// move to the next gene
		itGene++;
	}

	return ReadAssociatedWithAGene;
}

int ProcessReads(string SAM_File, map< string, vector< Gene > > *GenesList, int& nPoorQualityReads, int& nReadsMappedOutsideGenes) {
/**
	Returns the total number of reads
 */
	
	// open the read file stream
	ifstream ifs_sam ( SAM_File.c_str() );
	
	int rCount = 0;
	string strRead;
	string ReadRefName;
	string PrevReadRefName = "";
	int ReadStart, ReadEnd;
	vector<string> entries;
	vector< Gene > *RefGenesList;
	while (ifs_sam.good()) {
		getline(ifs_sam, strRead);

		if (strRead.empty()) { // ignore emtry lines
			continue;
		} else if (strRead[0] == '@') { // ignore the header
			continue;
		}

		// increment the read count
		rCount++;
				
		entries.clear();
		strsplit(strRead, entries, "\t");
		
		//filter out the poor quality read
		if ( atoi(entries[4].c_str()) < MQUAL_THRESHOLD ) {
			nPoorQualityReads++;
		} else {
			ReadRefName = entries[2];
			ReadStart = atoi(entries[3].c_str());
			ReadEnd = GetReadEndCoordinates(ReadStart, entries[5]);
			
			if (ReadRefName.compare(PrevReadRefName) != 0) { // get the gene list for this reference (chromosome), if we already dont have it
				// Get the list of genes corresponding to the read's reference name 
				map< string, vector< Gene > >::iterator it = GenesList->find(ReadRefName);
				
				// check if it is empty
				if (it !=  GenesList->end()) { // if not then get it
					RefGenesList = &it->second;				
				} else { // otherwise go to the next read
					continue;
				}
			}
			
			// update the frequency of read lying within, to the left or to the right of the gene
			bool ReadAssociatedWithAGene = UpdateReadsCount(RefGenesList, ReadStart, ReadEnd);

			// the read was not associated with any gene ... increment the number of reads mapped outside genes 
			if ( !ReadAssociatedWithAGene ) {
				nReadsMappedOutsideGenes++;
			}
			
			// save the reference name for the next iteration
			PrevReadRefName = ReadRefName;
			
		} 

		// show the update
		if (UpdateInterval > 0 && rCount % UpdateInterval == 0) {
			cout << rCount << endl;
		}

	}

	// close the gff file
	ifs_sam.close();
		
	return rCount;
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
		
	string GFF_File (argv[1]);
	string SAM_File (argv[2]);

	// Read the genes
	map< string, vector< Gene > > GenesList;	
	int gCount = ReadGenes(GFF_File, &GenesList);
	
/*	// initialise the array of reads count
	int ReadsCount[gCount][3];
	for (int gIdx = 0; gIdx < gCount; gIdx++) {
		for (int i = 0; i < 3; i++) {
			ReadsCount[gIdx][i] = 0;
		}
	}
*/	
	
	int nPoorQualityReads = 0;
	int nReadsMappedOutsideGenes = 0;
	// process the reads
	int nTotalReads = ProcessReads(SAM_File, &GenesList, nPoorQualityReads, nReadsMappedOutsideGenes);

	// open a file for outputting the results
	ofstream OutFile (SAM_File.append(".rc").c_str());//, ios::out | ios::app | ios::trunc );
	
	// output the results
	OutFile << "Total Reads : " << nTotalReads << endl;
	OutFile << "Poor Quality Reads (Threshold < " << MQUAL_THRESHOLD << "): " << nPoorQualityReads << endl;
	OutFile << "Reads Mapped Outside Genes: " << nReadsMappedOutsideGenes << endl;

	OutFile << "Reference\tGene Name\tReads Within\tReads Left\tReads Right" << endl;
	map< string, vector< Gene > >::iterator it = GenesList.begin();
	while (it != GenesList.end()) {
		vector< Gene > RefGenesList = it->second;
		vector< Gene >::iterator itGene = RefGenesList.begin();
		while (itGene != RefGenesList.end()) {
			OutFile << it->first << "\t" << itGene->GeneName << "\t" << itGene->nReadsWithin << "\t" << itGene->nReadsToTheLeft << "\t" << itGene->nReadsToTheRight << endl;
			itGene++;
		}
		it++;
	}
	OutFile.close();
			
	return 0;
}
