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

void ReadGFFFile(string GFF_File, map<string, vector< vector<string> > >& GFF) {

	ifstream ifs_gff ( GFF_File.c_str() );
	string strGFF = "";

	string reference = "";
	string previous_reference = "";
	vector< vector<string> > RefGFF;
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		}

		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFF[previous_reference] = RefGFF;
			}
			
			RefGFF = GFF[reference];
			previous_reference = reference;
		}

		RefGFF.push_back(entries);
		//cout << ID << endl;
	} // end while
	// update the last reference
	GFF[reference] = RefGFF;

}

void ReferenceToMultiFasta(string FASTA_File, string OUT_File) {

	map<string, vector< vector<string> > > GFF;
	ReadGFFFile(FASTA_File + ".gff", GFF);

	//  open the output file
	ofstream ofs (OUT_File.c_str());
	
	ifstream ifs_fasta ( FASTA_File.c_str() );
	string strFASTA = "";
	string reference = "";
	while (ifs_fasta.good()) {
		// read the next line from the fasta file
		getline(ifs_fasta, strFASTA);
		
		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}

		if (strFASTA[0] == '>') {
			reference = strFASTA.substr(1,strFASTA.length());
		}
		
		map<string, vector< vector<string> > >::iterator itGFF = GFF.find(reference);
		if (itGFF == GFF.end()) {
			cout << "Error. No entry found in the GFF file for the reference: " << reference << endl;
			return;
		}

		// get the current reference entry
		vector < vector<string> > RefGFF = itGFF->second;
		vector < vector<string> >::iterator itRefGFF = RefGFF.begin();
		
		int pos = 0;
		int last_pos = 0;
		string sequence = "";
		for (; itRefGFF != RefGFF.end(); itRefGFF++) {
		
			vector<string> entries = *itRefGFF;
			
			int start = atoi(entries[3].c_str());
			int end = atoi(entries[4].c_str());
			
			string header = entries[8].substr(entries[8].find("fasta_tag=",0) + 10, entries[8].length());
			
			//cout << start << "\t" << end << "\t" << last_pos << "\t" << start - last_pos << "\t" << header << endl;
			while (pos < end && ifs_fasta.good()) {
				
				// read the next line from the fasta file
				getline(ifs_fasta, strFASTA);
				
				if (strFASTA.empty()) { // ignore emtry lines
					continue;
				}
				
				if (strFASTA[0] == '>') {
					cout << "Error. Invalid Coordinates in GFF file." << endl;
					return;
				}
				
				RemoveNewLine(strFASTA);
				sequence += strFASTA;
				pos	+= strFASTA.length();
			} // end while
			
			ofs << ">" << header << endl;
			string dna = sequence.substr(start - last_pos - 1, end - start + 1);
			//cout << header << endl << dna << endl;
			int i = 0;
			int RefLength = dna.size();
			while (i < RefLength) {
				string str = dna.substr(i, BLOCK_SIZE);
				//str.erase(remove(str.begin(), str.end(), '\r'), str.end());
				//str.erase(remove(str.begin(), str.end(), '\n'), str.end());
				//transform(str.begin(), str.end(),str.begin(), ::toupper);
				ofs << str << endl;
				i += BLOCK_SIZE;
			}

			sequence = sequence.substr(end - last_pos - 1);//, sequence.length());
			//cout << endl << sequence << endl;
			last_pos = end - 1;
		}		
		
	}
	// close the stream
	ifs_fasta.close();
	ofs.close();
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "FASTA File not specified" << endl;
		exit(0);
	}
	string FASTA_File (argv[1]);
	string OUT_File (argv[2]);
	
	ReferenceToMultiFasta(FASTA_File, OUT_File);
	return 0;
}
