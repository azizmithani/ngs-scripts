/*
 *  Blast_Sequences.cpp
 *  
 *
 *  Created by Aziz Mithani on 20/05/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <string>
#include <vector>
#include <map>
#include "Utilities.h"

using namespace std;
/*
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
*/



int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "File not specified" << endl;
		exit(0);
	}
	
	string MappedLociFile (argv[1]);

	int skip = 1; // skip the first line (headings)
	vector< vector< string > > MappedLoci = ReadFile(MappedLociFile, skip);
	//cout << SequenceData.size() << endl;
	
	
	// process each mapped locus
	vector< vector< string > >::iterator it;
	for (it = MappedLoci.begin(); it != MappedLoci.end(); it++) {
		vector< string > entries = *it;  
		// open the streams for output
		ofstream ofsQuery ("in.txt");

		ofsQuery << entries[7] << endl;

		cout << "Blasting the est: " + entries[0] << endl;
		
		ofsQuery.close();
		string command = "blastn -query in.txt -db Brachypodium_distachyon -out " + entries[0] + ".txt -evalue 1000 -word_size 7 -penalty -3 -reward 1 -gapopen 4 -gapextend 2";
		system(command.c_str());
	}
	
	system("rm in.txt");
	
	return 0;
}
