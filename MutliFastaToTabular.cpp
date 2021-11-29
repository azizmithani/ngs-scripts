/*
 *  MultiFastaToReference.cpp
 *  
 *	Converts a multifasta file to tabular format.
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
	string OutputFile (argv[2]);
	string Tag;
	bool TagPresent;
	if (argc >= 4) {
		Tag = argv[3];
		TagPresent = true;
	} else {
		Tag = "";
		TagPresent = false;
	}
		
	int skip = 0; // no lines to skip
	vector< vector< string > > Sequences = ReadFile(SequenceFile, skip);
	
	vector< vector< string > >::iterator it = Sequences.begin();	
	if (it != Sequences.end()) {
		vector< string > entry = *it;
		string str = entry[0];
		if (str[0] != '>') {
			cout << "Invalid file format." << endl;
			return -1;
		} else {
			if (!TagPresent)
				Tag = str.substr(1,str.length()-1); // Tag for the next entry
		}

		it++;
	}
	// process the remaining file
	string seq = "";
	ofstream ofs (OutputFile.c_str());
	while (it != Sequences.end()) {
		vector< string > entry = *it;
		string str = entry[0];
		if (str[0] == '>') {
			ofs << Tag << "\t" << seq << endl;
			seq = ""; // re-initialise the variable
			if (!TagPresent)
				Tag = str.substr(1,str.length()-1); // Tag for the next entry
		} else
			seq += str;
		it++;
	}
	// write the last entry
	ofs << Tag << "\t" << seq << endl;
	
	ofs.close();
		
	return 0;
}

