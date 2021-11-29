/*
 *  chop_reads.cpp
 *  
 *  Chops the desired number of bases from the .fastq / .fa file  
 * 
 *  Created by Aziz Mithani on 20/04/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

int main (int argc, char** argv) {

	string DirName="/san/001/BGI/Data from BGI/Data from 21st Apr 2010/";
	string MutantName="E72";
	string LeftReadFilename = DirName + MutantName + "_1.fq";
	string RightReadFilename = DirName + MutantName + "_2.fq";
	string OutLeftReadFilename = DirName + MutantName + "_1_f.fq";
	string OutRightReadFilename = DirName + MutantName + "_2_f.fq";
	int LeftReadsLength = 73;
	int RightReadsLength = 75;
	
	// open the read file stream
	ifstream ifsLeft ( LeftReadFilename.c_str() );
	ifstream ifsRight ( RightReadFilename.c_str() );
	
	vector<string> ReadNames;
	map<string, vector<string> > LeftReads;
	string str = "";
	string str1 = "";
	string str2 = "";
	string str3 = "";
	while (ifsLeft.good()) {
		// get the lines containing read
		getline(ifsLeft, str);
		getline(ifsLeft, str1);
		getline(ifsLeft, str2);
		getline(ifsLeft, str3);
		
		if (str1.length() != LeftReadsLength)
			continue;
		vector<string> read;
		read.push_back(str);
		read.push_back(str1);
		read.push_back(str2);
		read.push_back(str3);
		
		string key = str.substr(0, str.length()-2);
		// store the read
		LeftReads.insert(pair<string, vector<string> > (key, read));
		
		// store the read name
//		ReadNames.push_back(key);
	}	
	cout << "Left reads: " << LeftReads.size() << endl;
	
	map<string, vector<string> > RightReads;
	str = "";
	str1 = "";
	str2 = "";
	str3 = "";
	while (ifsRight.good()) {
		// get the lines containing read
		getline(ifsRight, str);
		getline(ifsRight, str1);
		getline(ifsRight, str2);
		getline(ifsRight, str3);

		if (str1.length() != RightReadsLength)
			continue;

		string key = str.substr(0, str.length()-2);
		// move to the next read, if there is no matching left read 
		if (LeftReads.find(key) == LeftReads.end()) {
			//cout << str << endl;
			continue;
		}		
		
		vector <string> read;
		read.push_back(str);
		read.push_back(str1);
		read.push_back(str2);
		read.push_back(str3);
		
		// store the read
		RightReads.insert(pair<string, vector<string> > (key, read));
				
	}	
	cout << "Right reads: " << RightReads.size() << endl;

	
	// go through the left reads and remove the unpaired ones
	map<string, vector<string> > LeftReadsNew;
	LeftReadsNew.clear();
	map<string, vector<string> >::iterator it;
	for (it=LeftReads.begin(); it != LeftReads.end(); it++ ) {
		if ( RightReads.find(it->first) != RightReads.end() ) {		
			LeftReadsNew.insert(pair<string, vector<string> > (it->first, it->second));			
		} else {
			//cout << it->first << endl;
		}
	}
	cout << "Left reads: " << LeftReadsNew.size() << endl;

	// open a file for the output
	ofstream ofsLeft (OutLeftReadFilename.c_str());
	ofstream ofsRight (OutRightReadFilename.c_str());	

	// write the output
	for (it=LeftReadsNew.begin(); it != LeftReadsNew.end(); it++) {
		if (it->first.empty())
			continue;
		vector<string> read;
		vector<string>::iterator itRead;
		read = it->second;
		for ( itRead=read.begin() ; itRead != read.end(); itRead++ ) {
			ofsLeft << *itRead << endl;
		}
		
	}
	for (it=RightReads.begin(); it != RightReads.end(); it++) {
		if (it->first.empty())
			continue;
		vector<string> read;
		vector<string>::iterator itRead;
		read = it->second;
		for ( itRead=read.begin() ; itRead != read.end(); itRead++ ) {
			ofsRight << *itRead << endl;
		}
		
	}
	
	
/*	// write the output
	vector<string>::iterator it1;
	for (it1 = ReadNames.begin(); it1 != ReadNames.end(); it1++ ) {
		
		map<string, vector<string> >::iterator itReadEntry = LeftReadsNew.find(*it1);
		if (itReadEntry == LeftReadsNew.end()) { // this read name doesnot correspond to paired read, move to the next one
			continue;
		}
		
		vector<string> read;
		vector<string>::iterator itRead;

		read = itReadEntry->second;
		for ( itRead=read.begin() ; itRead != read.end(); itRead++ ) {
			ofsLeft << *itRead << endl;
		}

		read = RightReads.find(*it1)->second;
		for ( itRead=read.begin() ; itRead != read.end(); itRead++ ) {
			ofsRight << *itRead << endl;
		}
		
	}
*/	
	ofsLeft.close();
	ofsRight.close();
	
	return 0;
}
