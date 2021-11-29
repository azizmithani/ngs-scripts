/*
 *  dna_translator.cpp
 *  
 *
 *  Created by Aziz Mithani on 21/05/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */




#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "Utilities.h"

using namespace std;


int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "File not specified" << endl;
		exit(0);
	}
	
	string SequenceFile (argv[1]);
	string OuputMultiFastaFormat = "";
	if (argc >= 3) {
		OuputMultiFastaFormat = argv[2];
	}
	bool isOuputMultiFastaFormat = false;
	if (OuputMultiFastaFormat.length() > 0 && (OuputMultiFastaFormat[0] == 'T' || OuputMultiFastaFormat[0] == 't') ) {
		isOuputMultiFastaFormat = true;
	}
	
//	string OutputFile (argv[2]);
	
	int skip = 0; // no lines to skip
	vector< vector< string > > Sequences = ReadFile(SequenceFile, skip);
	
	string dna = "";
	string reference = "";
	string full_reference = "";
	vector< vector< string > >::iterator it = Sequences.begin();	
	if (it != Sequences.end()) {
		vector< string > entry = *it;
		string seq = entry[0];
		if (seq[0] != '>') {
			cout << "Invalid file format." << endl;
			return -1;
		} else {
			vector<string> entries;
			strsplit(seq, entries, " ");
			
			reference = entries[0].substr(1,entries[0].length());
			TrimLeadingSpaces(reference);
			full_reference = entry[0];
		}

		it++;
	}
	
	// process the remaining file
	while (it != Sequences.end()) {
		vector< string > entry = *it;
		string seq = entry[0];
		if (seq[0] == '>') {
			
			// new reference is starting ... translate the previous one 
			string SelectedProtein;
			int frame = TranslateDNA(dna, SelectedProtein);
			if (isOuputMultiFastaFormat) {
				//cout << ">" << reference << endl;
				cout << full_reference << endl;
				int pos = 0;
				int Length = SelectedProtein.size();
				while (pos < Length) {
					string str = SelectedProtein.substr(pos, 80);
					str.erase(remove(str.begin(), str.end(), '\r'), str.end());
					str.erase(remove(str.begin(), str.end(), '\n'), str.end());
					transform(str.begin(), str.end(),str.begin(), ::toupper);
					cout << str << endl;
					pos += 80;
				}
			} else {
				cout << reference << "\t" << dna.length() << "\t" << SelectedProtein.length() << "\t" << SelectedProtein << "\t" << frame << endl;
			}
			dna = "";
			
			// get the reference name
			vector<string> entries;
			strsplit(seq, entries, " ");
			
			reference = entries[0].substr(1,entries[0].length());
			TrimLeadingSpaces(reference);
			full_reference = entry[0];		
		} else{
			RemoveNewLine(seq);
			dna += seq;
		}
		it++;
	}
	// translate the final sequence 
	string SelectedProtein;
	int frame = TranslateDNA(dna, SelectedProtein);
	if (isOuputMultiFastaFormat) {
		//cout << ">" << reference << endl;
		cout << full_reference << endl;
		int pos = 0;
		int Length = SelectedProtein.size();
		while (pos < Length) {
			string str = SelectedProtein.substr(pos, 80);
			str.erase(remove(str.begin(), str.end(), '\r'), str.end());
			str.erase(remove(str.begin(), str.end(), '\n'), str.end());
			transform(str.begin(), str.end(),str.begin(), ::toupper);
			cout << str << endl;
			pos += 80;
		}
	} else {
		cout << reference << "\t" << dna.length() << "\t" << SelectedProtein.length() << "\t" << SelectedProtein << "\t" << frame << endl;
	}
	
	
/*	
//	cout << reference << endl;
	// open the output stream
	ofstream ofs (OutputFile.c_str());
	ofs << ">" + Header << endl;
	int pos = 0;
	int RefLength = reference.size();
	while (pos < RefLength) {
		string str = reference.substr(pos, 80);
		str.erase(remove(str.begin(), str.end(), '\r'), str.end());
		str.erase(remove(str.begin(), str.end(), '\n'), str.end());
		transform(str.begin(), str.end(),str.begin(), ::toupper);
		ofs << str << endl;
		pos += 80;
	}
	
//	ofs << endl;
	
	ofs.close();
*/		
	return 0;
}

