/*
 *  change_sequence_orientation.cpp
 *  reads a fasta or multi-fasta file and changes the sequence orientation, if required, to the one that gives the largest protein coding sequence  
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
			// new reference is starting ... process the previous one 
		
			// convert to upper case
			dna.erase(remove(dna.begin(), dna.end(), '\r'), dna.end());
			dna.erase(remove(dna.begin(), dna.end(), '\n'), dna.end());
			transform(dna.begin(), dna.end(),dna.begin(), ::toupper);
			string polyA (10, 'A');
			string polyT (10, 'T');
			
			string SelectedSequence = "";
			// see if poly A tail is present and starts within last 25 bases
			int polyA_pos = dna.find(polyA, dna.length() - 25);
			if ( polyA_pos != string::npos) { 
				//cout << "poly A sequence found" << endl << full_reference << endl << dna << endl;
				SelectedSequence = dna;
		    } else { // otherwise see if poly T is present and starts within first 15 bases
				int polyT_pos = dna.find(polyT);				
				if ( polyT_pos != string::npos && polyT_pos < 15 ) {
					//cout << "poly T sequence found" << endl << full_reference << endl << dna << endl;
					SelectedSequence = ReverseComplementDNA(dna);
				}
			}

			if (SelectedSequence.length() == 0) {
				// check the proteins
				string rev_dna = ReverseComplementDNA(dna);
				
				string protein = TranslateDNA(dna);
				string rev_protein = TranslateDNA(rev_dna);
								
				if (protein.size() < rev_protein.size()) {
					//cout << full_reference << endl << dna << endl << full_reference << endl << rev_dna << endl;
					SelectedSequence = rev_dna;
				} else {
					SelectedSequence = dna;
				}
			}
			
			cout << full_reference << endl;
			int pos = 0;
			int Length = SelectedSequence.size();
			while (pos < Length) {
				string str = SelectedSequence.substr(pos, 80);
				cout << str << endl;
				pos += 80;
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
	
	// process the final sequence 
	// convert to upper case
	dna.erase(remove(dna.begin(), dna.end(), '\r'), dna.end());
	dna.erase(remove(dna.begin(), dna.end(), '\n'), dna.end());
	transform(dna.begin(), dna.end(),dna.begin(), ::toupper);
	string polyA (10, 'A');
	string polyT (10, 'T');
	
	string SelectedSequence = "";
	// see if poly A tail is present and starts within last 25 bases
	int polyA_pos = dna.find(polyA, dna.length() - 25);
	if ( polyA_pos != string::npos ) { 
		SelectedSequence = dna;
	} else { // otherwise see if poly T is present and starts within first 15 bases
		int polyT_pos = dna.find(polyT);
		if ( polyT_pos != string::npos && polyT_pos < 15 ) {
			SelectedSequence = ReverseComplementDNA(dna);
		}
	}
	
	if (SelectedSequence.length() == 0) {
		// check the proteins
		string rev_dna = ReverseComplementDNA(dna);
		
		string protein = TranslateDNA(dna);
		string rev_protein = TranslateDNA(rev_dna);
		
		if (protein.size() < rev_protein.size()) {
			//cout << full_reference << endl << dna << endl << full_reference << endl << rev_dna << endl;
			SelectedSequence = rev_dna;
		} else {
			SelectedSequence = dna;
		}
	}
	
	cout << full_reference << endl;
	int pos = 0;
	int Length = SelectedSequence.size();
	while (pos < Length) {
		string str = SelectedSequence.substr(pos, 80);
		cout << str << endl;
		pos += 80;
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

