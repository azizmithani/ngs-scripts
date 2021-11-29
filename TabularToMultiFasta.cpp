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

	int skip = 0; // no lines to skip
	vector< vector< string > > Sequences = ReadFile(SequenceFile, skip);
	//cout << SequenceData.size() << endl;
	
	// open the output stream
	ofstream ofs (OutputFile.c_str());
	// write the sequence to the output
	//vector< vector< string > >::iterator it;
	for (vector< vector< string > >::iterator it = Sequences.begin(); it != Sequences.end(); it++) {
		vector< string > entry = *it;

		if (entry[0].length() == 0) {
			continue;
		}
		ofs << ">" + entry[0] << endl;
		string sequence = entry[1];
		int pos = 0;
		int SeqLength = sequence.size();
		while (pos < SeqLength) {
			ofs << sequence.substr(pos, 80) << endl;
			pos += 80;
		}
	}
		
	ofs.close();
		
	return 0;
}

