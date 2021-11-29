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


void ReferenceToMultiFasta(string FASTA_File, string OUT_File, int Gap) {
	
	string sep (Gap, 'N');
	
	//  open the output file
	ofstream ofs (OUT_File.c_str());
	
	ifstream ifs_fasta ( FASTA_File.c_str() );
	string strFASTA = "";
	string reference = "";
	int idx = 0;
	while (ifs_fasta.good()) {
		// read the next line from the fasta file
		getline(ifs_fasta, strFASTA);

		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}
		
		// match the reference
		if (strFASTA[0] == '>') { // we are starting a new reference .. process the previous one
			if (reference.length() == 0)
				continue;
			int pos = reference.find(sep);
			while (pos != string::npos) {
				stringstream sIdx;
				sIdx << ++idx;
				ofs << ">" + sIdx.str() << endl;

				//cout << pos;
				string dna = reference.substr(0,pos);
				int i = 0;
				int RefLength = dna.size();
				while (i < RefLength) {
					string str = dna.substr(i, 80);
					//str.erase(remove(str.begin(), str.end(), '\r'), str.end());
					//str.erase(remove(str.begin(), str.end(), '\n'), str.end());
					//transform(str.begin(), str.end(),str.begin(), ::toupper);
					ofs << str << endl;
					i += 80;
				}
				reference = reference.substr(pos + 1,reference.length());
				pos = reference.find_first_not_of('N');
				//cout << "\t" << pos << endl;
				reference = reference.substr(pos,reference.length());
				//cout << reference.substr(0,20) << endl;
				pos = reference.find(sep);
			} // end while
			stringstream sIdx;
			sIdx << ++idx;
			ofs << ">" + sIdx.str() << endl;
			
			//cout << pos;
			string dna = reference.substr(0,reference.length());
			int i = 0;
			int RefLength = dna.size();
			while (i < RefLength) {
				string str = dna.substr(i, 80);
				//str.erase(remove(str.begin(), str.end(), '\r'), str.end());
				//str.erase(remove(str.begin(), str.end(), '\n'), str.end());
				//transform(str.begin(), str.end(),str.begin(), ::toupper);
				ofs << str << endl;
				i += 80;
			}
			
			reference = "";
			
		} else {
			RemoveNewLine(strFASTA);
			reference += strFASTA;
		}		
	} // end while ifs_fasta.good()
	int pos = reference.find(sep);
	while (pos != string::npos) {
		stringstream sIdx;
		sIdx << ++idx;
		ofs << ">" + sIdx.str() << endl;
		
		//cout << idx << "\t" << pos;
		string dna = reference.substr(0,pos);
		int i = 0;
		int RefLength = dna.size();
		while (i < RefLength) {
			string str = dna.substr(i, 80);
			//str.erase(remove(str.begin(), str.end(), '\r'), str.end());
			//str.erase(remove(str.begin(), str.end(), '\n'), str.end());
			//transform(str.begin(), str.end(),str.begin(), ::toupper);
			ofs << str << endl;
			i += 80;
		}
		reference = reference.substr(pos + 1,reference.length());
		pos = reference.find_first_not_of('N');
		//cout << "\t" << pos << endl;
		reference = reference.substr(pos,reference.length());
		//cout << reference.substr(0,40) << endl;
		pos = reference.find(sep);
	} // end while
	stringstream sIdx;
	sIdx << ++idx;
	ofs << ">" + sIdx.str() << endl;
	
	//cout << pos;
	string dna = reference.substr(0,reference.length());
	int i = 0;
	int RefLength = dna.size();
	while (i < RefLength) {
		string str = dna.substr(i, 80);
		//str.erase(remove(str.begin(), str.end(), '\r'), str.end());
		//str.erase(remove(str.begin(), str.end(), '\n'), str.end());
		//transform(str.begin(), str.end(),str.begin(), ::toupper);
		ofs << str << endl;
		i += 80;
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
	int Gap(atoi(argv[3]));
	
	ReferenceToMultiFasta(FASTA_File, OUT_File, Gap);
	return 0;
}
