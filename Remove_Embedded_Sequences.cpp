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

vector< string > ReadSequenceFile(string FileName){
	
	// open the file stream
	ifstream ifs ( FileName.c_str() );
	
	string str="";
	vector< string > v;
	while (ifs.good()) {
		getline(ifs, str);
				
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		v.push_back (str);
		
	}
	
	ifs.close();
	
	return v;
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string SequenceFile (argv[1]);
	string OutputFile (argv[2]);
		
	vector< string > Sequences = ReadSequenceFile(SequenceFile);
	//cout << SequenceData.size() << endl;
	
	// open the output stream
	ofstream ofs (OutputFile.c_str());
	for (vector< string >::iterator it = Sequences.begin(); it != Sequences.end(); it++ ) { 
		string str = *it;
		
		// split the entries
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		string Query = entries[1];
		
		bool isEmbedded = false;
		for (vector< string >::iterator it1 = Sequences.begin(); it1 != Sequences.end(); it1++ ) { 
			
			string Text = *it1;
			
			// cout << "\t" << Text << endl;
			if (it == it1) {// don't match the entry against itself
				continue;
			}
/*			if (str.compare(Text) == 0) {// don't match the entry against itself
				cout << &it << "\t" << &it1 << endl;
				continue;
*/			else if (Query.size() > Text.size())
				continue; // query sequence cannot be embedded in the text if the text is smaller than query
			else if ( Text.find(Query) == string::npos) // query sequence not found in the text
				continue;
			// we come here if query sequence is present in the text. Mark the query sequence to be removed
			isEmbedded = true;
			break;
		}
		// if the query sequence is not embedded in any other sequence, then write it to the output file 
		if (!isEmbedded)
			ofs << str << endl;
	}	
		
	ofs.close();

	
	return 0;
}

/*
int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string SequenceFile (argv[1]);
	string OutputFile (argv[2]);
	
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
	map< string, vector<string> > SequencesToBeRemoved;
	for (map<string, vector< string > >::iterator it = GroupedSequences.begin(); it != GroupedSequences.end(); it++ ) { 
		vector< string > ChrSequences =  it->second;
		
		for (vector< string >::iterator itSeq = ChrSequences.begin(); itSeq != ChrSequences.end(); itSeq++) {
			string Query = *itSeq;
			bool isEmbedded = false;
			//			cout << Query << endl;
			for (map<string, vector< string > >::iterator it1 = GroupedSequences.begin(); it1 != GroupedSequences.end(); it1++ ) { 
				vector< string > ChrSequences1 =  it1->second;
				
				for (vector< string >::iterator itSeq1 = ChrSequences1.begin(); itSeq1 != ChrSequences1.end(); itSeq1++) {
					string Text = *itSeq1;
					
					//					cout << "\t" << Text << endl;
					if (Query.compare(Text) == 0)
						continue;
					else if (Query.size() > Text.size())
						continue; // query sequence cannot be embedded in the text if the text is smaller than query
					else if ( Text.find(Query) == string::npos) // query sequence not found in the text
						continue;
					// we come here if query sequence is present in the text. Mark the query sequence to be removed
					isEmbedded = true;
					break;
				}
				if (isEmbedded) // don't look in other groups if one hit is found
					break;
			}
			// if the query sequence is not embedded in any other sequence, then write it to the output file 
			if (!isEmbedded)
				ofs << it->first << "\t" << Query << endl;
		}
	}	
	
	ofs.close();
	
	
	return 0;
}

*/
