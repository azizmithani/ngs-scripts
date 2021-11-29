/*
 *  find_differences.cpp
 *  
 *
 *  Created by Aziz Mithani on 07/05/2010.
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
#include <set>
#include <stdio.h>

using namespace std;

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

set< string > ReadFile(string FileName){
	
	// open the read file stream
	ifstream ifs ( FileName.c_str() );

	string str="";
	set< string > v;
	while (ifs.good()) {
		getline(ifs, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		v.insert (entries[0] + ":" + entries[1]);
		//v.push_back (pair<string,int> (entries[0],atoi(entries[1].c_str())));
		
	}
	
	ifs.close();
	
	return v;
}


map< string, vector<string> > ReadFileCompletely(string FileName) {
	
	// open the read file stream
	ifstream ifs ( FileName.c_str() );
	
	string str="";
	map< string, vector<string> > m;
	while (ifs.good()) {
		getline(ifs, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		m.insert( pair<string, vector<string> >(entries[0] + ":" + entries[1], entries));
		
	}
	
	ifs.close();
	
	return m;
}

int FindDifferences(string File_1, string File_2, string OUT_File_Diff) {
	/**
	 Reads the two tab-seperated file and performs set minus File_1 - File_2 and writes the output 
	 */
	
	
	set< string > v_1 = ReadFile(File_1);
	set< string > v_2 = ReadFile(File_2);
	// allocate a vector for the differences
	//vector< pair<string,int> > v_diff (v_1.size());

	map< string, vector<string> > m_1 = ReadFileCompletely(File_1);
	
	// find the differences
	//set_difference(v_1.begin(), v_1.end(), v_2.begin(), v_2.end(), v_diff.begin());
	
	
	// open the streams for output
	ofstream ofsDiff (OUT_File_Diff.c_str());
	set<string>::iterator it;
	for ( it=v_1.begin() ; it != v_1.end(); it++ ) {
		string entry = *it;
	 	if (entry.empty()) {
			continue;
		} else if (v_2.find(entry) != v_2.end()) { // if this item is present in the second list, ignore and move to the next one
			continue;
		}
		
		// output the complete entry from the map
		string key = entry;
		vector<string>::iterator itRow;
		for (itRow = m_1[key].begin(); itRow != m_1[key].end(); itRow++) {
			ofsDiff << *itRow << "\t";
		}
		ofsDiff << endl;
	}
	
	// close the file
	ofsDiff.close();
	
	return 0;
}


int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	// read the file names
	int nFiles = argc;
	vector< string > Filenames;
	for (int i = 1; i < argc; i++) {
		string Filename (argv[i]);
		Filenames.push_back(Filename);		
	}
	
	string outFilename = "temp.txt";
	// in a loop, read all files but one, write it to a new file and call the 'find_differences' script, delete the file
	for (vector< string >::iterator it1 = Filenames.begin(); it1 != Filenames.end(); it1++) {
		set< string > UnionSet;
		for (vector< string >::iterator it2 = Filenames.begin(); it2 != Filenames.end(); it2++) {
			if (it1 == it2) {
				continue;
			}
			set< string > Mutations = ReadFile(*it2);
			UnionSet.insert(Mutations.begin(), Mutations.end());
		}
		// write the union set to the file
		ofstream ofs (outFilename.c_str());
		for (set< string >::iterator it = UnionSet.begin(); it != UnionSet.end(); it++) {
			vector<string> entries;
			strsplit(*it, entries, ":");
			
			ofs << entries[0] + "\t" << entries[1] << endl;
		}
		ofs.close();
		FindDifferences(*it1, outFilename, *it1 + ".uniq");
	}
	remove(outFilename.c_str());
	
	return 0;
}


/*
int FindDifferences(string File_1, string File_2, string OUT_File_Diff) {
	/**
	 Reads the two tab-seperated file and performs set minus File_1 - File_2 and writes the output 
	 */
/*	
	
	vector< pair<string,int> > v_1 = ReadFile(File_1);
	vector< pair<string,int> > v_2 = ReadFile(File_2);
	// allocate a vector for the differences
	vector< pair<string,int> > v_diff (v_1.size());
	
	map< string, vector<string> > m_1 = ReadFileCompletely(File_1);
	
	// find the differences
	set_difference(v_1.begin(), v_1.end(), v_2.begin(), v_2.end(), v_diff.begin());
	
	
	// open the streams for output
	ofstream ofsDiff (OUT_File_Diff.c_str());
	vector< pair<string,int> >::iterator it;
	for ( it=v_diff.begin() ; it != v_diff.end(); it++ ) {
		pair<string,int> entry = *it;
		if (entry.first.empty())
			continue;
		
		// output the complete entry from the map
		stringstream ss;
		ss << entry.second;
		string key = entry.first + ":" + ss.str(); 
		vector<string>::iterator itRow;
		for (itRow = m_1[key].begin(); itRow != m_1[key].end(); itRow++) {
			ofsDiff << *itRow << "\t";
		}
		ofsDiff << endl;
	}
	
	// close the file
	ofsDiff.close();
	
	return 0;
}
*/


/*

vector< pair<string,int> > ReadFile(string FileName){
	
	// open the read file stream
	ifstream ifs ( FileName.c_str() );
	
	string str="";
	vector< pair<string,int> > v;
	while (ifs.good()) {
		getline(ifs, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		v.push_back (pair<string,int> (entries[0],atoi(entries[1].c_str())));
		
	}
	
	ifs.close();
	
	return v;
}

int FindDifferences(string File_1, string File_2, string OUT_File_Diff) {
	
	
	vector< pair<string,int> > v_1 = ReadFile(File_1);
	vector< pair<string,int> > v_2 = ReadFile(File_2);
	// allocate a vector for the differences
	vector< pair<string,int> > v_diff (v_1.size());
	
	// find the differences
	set_difference(v_1.begin(), v_1.end(), v_2.begin(), v_2.end(), v_diff.begin());
	
	
	// open the streams for output
	ofstream ofsDiff (OUT_File_Diff.c_str());
	vector< pair<string,int> >::iterator it;
	for ( it=v_diff.begin() ; it != v_diff.end(); it++ ) {
		pair<string,int> entry = *it;
		if (entry.first.empty())
			continue;
		ofsDiff << entry.first << "\t" << entry.second << endl;
	}
	
	// close the file
	ofsDiff.close();
	
	return 0;
}


int main (int argc, char** argv) {
	
	if (argc < 4) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string File_1 (argv[1]);
	string File_2 (argv[2]);
	string File_out (argv[3]);
	
	FindDifferences(File_1, File_2, File_out);
	
	return 0;
}
*/
