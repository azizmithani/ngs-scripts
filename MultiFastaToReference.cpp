/*
 *  MultiFastaToReference.cpp
 *  
 *	Creates a reference genome by using all the sequences from the given file.
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
	int Gap (atoi(argv[3]));
	string Header (argv[4]);
	
	string seperator( Gap, 'N' );
	
	int skip = 0; // no lines to skip
	vector< vector< string > > Sequences = ReadFile(SequenceFile, skip);
	
	vector< vector< string > >::iterator it = Sequences.begin();	
	if (it != Sequences.end()) {
		vector< string > entry = *it;
		string seq = entry[0];
		if (seq[0] != '>') {
			cout << "Invalid file format." << endl;
			return -1;
		}
		it++;
	}
	// process the remaining file
	string reference = "";
	while (it != Sequences.end()) {
		vector< string > entry = *it;
		string seq = entry[0];
		if (seq[0] == '>')
			reference += seperator;
		else{
			RemoveNewLine(seq);
			reference += seq;
		}
		it++;
	}
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
		
	return 0;
}

