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

void split(const string& str, vector<string>& tokens, const string& delimiters = " ") {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
	
    if (string::npos != pos && string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        lastPos = str.find_first_not_of(delimiters, pos);
		tokens.push_back(str.substr(lastPos, str.length()));
        // Skip delimiters.  Note the "not_of"
        // Find next "non-delimiter"
        //pos = str.find_first_of(delimiters, lastPos);
    } else {
		tokens.push_back(str);
	}


	//cout << tokens[0] << " : " << tokens[1] << endl;
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "File not specified" << endl;
		exit(0);
	}
	
	string DataFile (argv[1]);
	
	string ID = "-";
	string Gene = "-";
	string Description = "-";
	
	cout << "ID\tGene\tDescription" << endl;
	string str="";
	// open the sam file stream
	ifstream ifs ( DataFile.c_str() );
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else { 
			vector<string> entries;
			split(str, entries, " ");
			
			if (entries[0].compare("ID") == 0) {
				ID = entries[1];
			} else if (entries[0].compare("GENE") == 0) {
				Gene = entries[1];
			} else if (entries[0].compare("TITLE") == 0) {
				Description = entries[1];
			} else if (entries[0].compare("//") == 0) {
				cout << ID << "\t" << Gene << "\t" << Description << endl;
				ID = "-";
				Gene = "-";
				Description = "-";
			}
			
		} // end if str.empty()
	} // while ifs.close()
	
	ifs.close();
	
	return 0;
}

