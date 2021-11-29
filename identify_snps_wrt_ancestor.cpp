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

using namespace std;

const int BLOCK_SIZE = 80;

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "FASTA File not specified" << endl;
		exit(0);
	}
	string AncestorFile (argv[1]);
	string DiploidFile (argv[2]);
	string GenomeFile (argv[3]);
	string OutFile (argv[4]);

	ofstream ofsOut (OutFile.c_str());

	ifstream ifsAncestor ( AncestorFile.c_str() );
	ifstream ifsGenome ( GenomeFile.c_str() );
	ifstream ifsDiploid ( DiploidFile.c_str() );
	string strAncestor = "";
	string strGenome = "";
	string strDiploid = "";
	string Reference = "";
	int fastaTotalReadPos = 0; // position until which the file has been read
	while (ifsAncestor.good()) {
		// read the next line from the fasta files
		getline(ifsAncestor, strAncestor);
		getline(ifsGenome, strGenome);
		getline(ifsDiploid, strDiploid);
		
		if (strAncestor.empty()) { // ignore emtry lines
			continue;
		}
		
		// match the reference
		if (strAncestor[0] == '>') { // we are reading a new reference (chromosome)
			Reference = ExtractReferenceName(strAncestor);
			fastaTotalReadPos = 0;
			continue;
	    }
		
		for (int i = 0; i < strAncestor.length(); i++) {
			if (strAncestor[i] == 'N') {
				continue;
			}
			if (strDiploid[i] == '*' || strDiploid[i] == '0' || strDiploid[i] == ' ') {
				continue;
			}
			if (strGenome[i] == ' ' || strGenome[i] == '0') {
				continue;
			}
			if (strAncestor[i] != strGenome[i] || strAncestor[i] != strDiploid[i]) {
				ofsOut << Reference << "\t" << fastaTotalReadPos + i + 1 << "\t" << strAncestor[i] << "\t" <<  strDiploid[i] << "\t" <<  strGenome[i] << endl;
			}
		}
		fastaTotalReadPos += strAncestor.length();
	} // end while ifsAncestor is good
	
	ofsOut.close();
	return 0;
}
