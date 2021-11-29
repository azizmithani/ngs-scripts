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
#include <algorithm>
#include <map>
#include "Utilities.h"

using namespace std;

const int BLOCK_SIZE = 80;
const string CLUSTALW = "/san/001/src/clustalw-2.0.12/src/clustalw2";


void ReadGFFFile(string GFF_File, map<string, vector< vector<string> > >& GFF) {

	ifstream ifs_gff ( GFF_File.c_str() );
	string strGFF = "";

	string reference = "";
	string previous_reference = "";
	vector< vector<string> > RefGFF;
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF[0] == '#') { // ignore header
			continue;
		}

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

void ReadFASTAFile (string FASTA_File, map<string, string>& Sequences) {
	string strFASTA = "";
	ifstream ifs_fasta ( FASTA_File.c_str() );
	string sequence = "";
	string reference = "";
	while (ifs_fasta.good()) {
		// read the next line from the fasta file
		getline(ifs_fasta, strFASTA);
		
		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}
		
		if (strFASTA[0] == '>') {
			if (reference.length() > 0) {
				Sequences[reference] = sequence;
			}
			reference = strFASTA.substr(1,strFASTA.length());
			sequence = "";
			
			continue;
		}
		
		RemoveNewLine(strFASTA);
		sequence += strFASTA;
	
	}		
	Sequences[reference] = sequence;
	
	ifs_fasta.close();

}

void AlignGenes(string FASTA_File_1, string FASTA_File_2, string OUT_File) {

	map<string, vector< vector<string> > > GFF_1;
	ReadGFFFile(FASTA_File_1 + ".gff", GFF_1);

	map<string, vector< vector<string> > > GFF_2;
	ReadGFFFile(FASTA_File_2 + ".gff", GFF_2);

	map<string, string> Sequences_1;
	ReadFASTAFile(FASTA_File_1, Sequences_1);

	map<string, string> Sequences_2;
	ReadFASTAFile(FASTA_File_2, Sequences_2);

	//  open the output file
	ofstream ofs (OUT_File.c_str());
	
	map<string, vector< vector<string> > >::iterator itGFF_1 = GFF_1.begin();
	map<string, vector< vector<string> > >::iterator itGFF_2;
	string RefSequence1, RefSequence2;
	string reference = "";
	string previous_reference = "";
	for (; itGFF_1 != GFF_1.end(); itGFF_1++) {

		// Get the reference names
		reference = itGFF_1->first;

		// find the entry for this reference in GFF2
		itGFF_2 = GFF_2.find(reference);
		if (itGFF_2 == GFF_2.end()) {
			cout << "There is some problem with the GFF files" << endl;
			exit(0);
		}
		
		// get the reference sequences
		RefSequence1 = Sequences_1[reference];
		RefSequence2 = Sequences_2[reference];

		// get the GFF entries for this reference
		vector< vector<string> > RefGFF1, RefGFF2;
		RefGFF1 = itGFF_1->second;
		RefGFF2 = itGFF_2->second;

		vector< vector<string> >::iterator itRefGFF1 = RefGFF1.begin();
		vector< vector<string> >::iterator itRefGFF2 = RefGFF2.begin();
		for (; itRefGFF1 != RefGFF1.end() && itRefGFF2 != RefGFF2.end(); itRefGFF1++, itRefGFF2++) {
			vector<string> entries1 = *itRefGFF1;
			vector<string> entries2 = *itRefGFF2;
			
			int start1 = atoi(entries1[3].c_str());
			int start2 = atoi(entries2[3].c_str());

			int end1 = atoi(entries1[4].c_str());
			int end2 = atoi(entries2[4].c_str());

			cout << reference << "\t" << start1 << "\t" << end1 << "\t" << start2 << "\t" << end2 << endl;

			string TEMP_File = OUT_File + ".tmp";
			ofstream ofs_temp (TEMP_File.c_str());
			ofs_temp << ">" << FASTA_File_1 << endl;
			ofs_temp << RefSequence1.substr(start1 - 1, end2 - start1 + 1) << endl;

			ofs_temp << ">" << FASTA_File_2 << endl;
			ofs_temp << RefSequence2.substr(start2 - 1, end2 - start2 + 1) << endl;
			
			ofs_temp.close();
			
			// run ClustalW
			string command = CLUSTALW + " " + TEMP_File + " -quiet -stats=" + TEMP_File + ".stats";
			system(command.c_str());
			
			// extract the alignement score
			command = "grep \"aln pw-id avg\" " + TEMP_File + ".stats | cut -d ':' -f2 | sed 's/^[ \t]*//'  > " + TEMP_File + ".out";
			system (command.c_str());
			
			ifstream ifs_temp ((TEMP_File + ".out").c_str());
			string strTEMP;
			getline(ifs_temp, strTEMP);
			ofs << reference << "\t" << start1 << "\t" << end1 << "\t" << start2 << "\t" << end2 << "\t" << atof(strTEMP.c_str()) * 100 << endl;
			
			command = "rm " + TEMP_File + "* " + OUT_File + ".dnd " + OUT_File + ".aln";
			system (command.c_str());
			
		//	break;
		}
	}
	// close the stream
	ofs.close();
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "FASTA Files not specified" << endl;
		exit(0);
	}
	string FASTA_File_1 (argv[1]);
	string FASTA_File_2 (argv[2]);
	string OUT_File (argv[3]);
	
	AlignGenes(FASTA_File_1, FASTA_File_2, OUT_File);
	return 0;
}
