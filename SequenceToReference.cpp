/*
 *  SequenceToReference.cpp
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
	
	string seperator( Gap, 'N' );
	
	int skip = 0; // no lines to skip
	vector< vector< string > > Sequences = ReadFile(SequenceFile, skip);
	//cout << SequenceData.size() << endl;
	
	// group the sequences according to chromosomes (or group ids)
	map< string, vector<string> > GroupedSequences;
	//vector< vector< string > >::iterator it;
	for (vector< vector< string > >::iterator it = Sequences.begin(); it != Sequences.end(); it++) {
		vector< string > entry = *it;
		GroupedSequences[entry[0]].push_back(entry[1]);
	}
	
	// open the output stream
	ofstream ofs (OutputFile.c_str());
	//map<string, vector< string > >::iterator it;
	// write the first sequence to the output
	for (map<string, vector< string > >::iterator it = GroupedSequences.begin(); it != GroupedSequences.end(); it++ ) { 

		string reference = "";
		vector< string > ChrSequences =  it->second;
		vector< string >::iterator itSeq = ChrSequences.begin();
		if (itSeq != ChrSequences.end()) {
			reference += *itSeq;
			itSeq++;
		}
		// append all the remaining sequence. also add extra 'N'
		while (itSeq != ChrSequences.end()) {
			reference += seperator + *itSeq;
			itSeq++;
		}

		ofs << ">" + it->first << endl;
		int pos = 0;
		int RefLength = reference.size();
		while (pos < RefLength) {
			ofs << reference.substr(pos, 80) << endl;
			pos += 80;
		}
	}
	
	ofs << endl;
	
	ofs.close();
		
	return 0;
}

