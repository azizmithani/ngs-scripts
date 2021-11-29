/*
 *
 *  Created by Aziz Mithani on 21/05/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */




#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

using namespace std;

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

void ReadIdsToExtract (string Ids_File, set<string>& IdsToExtract) {

	// read the  file
	ifstream ifs ( Ids_File.c_str() );
	string str = "";
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		}
		
		IdsToExtract.insert(str);
		
	}
	
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "File not specified" << endl;
		exit(0);
	}
	
	string Sequence_File (argv[1]);
	string Ids_File (argv[2]);

	set<string> IdsToExtract;
	ReadIdsToExtract(Ids_File, IdsToExtract);
	
	// read the sequence file
	ifstream ifs_seq ( Sequence_File.c_str() );
	string strSeq = "";
	string Id = "";
	bool OutputSequence = false;
	while (ifs_seq.good()) {
		getline(ifs_seq, strSeq);
		
		if (strSeq.empty()) { // ignore emtry lines
			continue;
		}
		
		if (strSeq[0] == '>') {			
			// Extract the Id for the new sequence
			Id = ExtractReferenceName(strSeq);
			
			if (IdsToExtract.find(Id) != IdsToExtract.end()) {
				OutputSequence = true;
				cout << strSeq << endl;
			} else {
				OutputSequence = false;
			}
		} else{
			if (OutputSequence) {
				cout << strSeq << endl;
			}
		}
	}
	
	return 0;
}

