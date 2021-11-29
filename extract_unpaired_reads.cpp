/*
 *  filter_sam_file.cpp
 *  
 *	Removes low quality, unpaired and non-unique reads from a sam file
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
#include <set>
#include <map>
#include "Utilities.h"

using namespace std;


const int FLAG_PAIRED_READ = 2;

set <string> LeftReads, RightReads;
void ExtractReads(string SAM_File, string Left_Read_File, string Right_Read_File, string OUT_Left_Read_File, string OUT_Right_Read_File) {
	// open the sam file stream
	ifstream ifs_sam ( SAM_File.c_str() );

	string strRead="";
	while (ifs_sam.good()) {
		getline(ifs_sam, strRead);
		
		if (strRead.empty()) { // ignore emtry lines
			continue;
		} else if (strRead[0] == '@') { // ignore the header
			continue;
		}

		vector<string> entries;
		strsplit(strRead, entries, "\t");

		if ( !( atoi(entries[1].c_str()) & FLAG_PAIRED_READ ) ) {
			string ReadName = entries[0];
			if (ReadName[ReadName.length() - 1] == '1' ) {
				LeftReads.insert("@" + ReadName);
			} else if (ReadName[ReadName.length() - 1] == '2' ) {
				RightReads.insert("@" + ReadName);
			} else {
				LeftReads.insert("@" + ReadName + "/1");
				RightReads.insert("@" + ReadName + "/2");
			}
		}
				
	} // end while
	
	// close the file
	ifs_sam.close();
	
	// open the left read file stream
	ifstream ifs_left ( Left_Read_File.c_str() );
	ofstream ofs_left ( OUT_Left_Read_File.c_str() );
	string strLeft="";
	while (ifs_left.good()) {
		getline(ifs_left, strLeft);
		
		if (strLeft.empty()) { // ignore emtry lines
			continue;
		}
		
		//cout << strLeft << "\t" << strLeft.substr(0,strLeft.length()-2) << endl;
		if (LeftReads.find(strLeft) != LeftReads.end() ){// || LeftReads.find(strLeft.substr(0,strLeft.length()-2)) != LeftReads.end()) {
			ofs_left << strLeft << endl;
			for (int i = 0; i < 3; i++) {
				getline(ifs_left, strLeft);
				ofs_left << strLeft  << endl;
			}
		} else {
			for (int i = 0; i < 3; i++) {
				getline(ifs_left, strLeft);
			}
		}

	} // end while
	ifs_left.close();
	ofs_left.close();
	
	// open the right read file stream
	ifstream ifs_right ( Right_Read_File.c_str() );
	ofstream ofs_right ( OUT_Right_Read_File.c_str() );
	string strRight="";
	while (ifs_right.good()) {
		getline(ifs_right, strRight);
		
		if (strRight.empty()) { // ignore emtry lines
			continue;
		}
		
		if (RightReads.find(strRight) != RightReads.end() ) { //|| RightReads.find(strRight.substr(0,strRight.length()-2)) != RightReads.end()) {
			ofs_right << strRight << endl;
			for (int i = 0; i < 3; i++) {
				getline(ifs_right, strRight);
				ofs_right << strRight  << endl;
			}
		} else {
			for (int i = 0; i < 3; i++) {
				getline(ifs_right, strRight);
			}
		}
	} // end while
	ifs_right.close();
	ofs_right.close();
	
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string SAM_File (argv[1]);
	string Left_Read_File (argv[2]);
	string Right_Read_File (argv[3]);
	string OUT_Left_Read_File (argv[4]);
	string OUT_Right_Read_File (argv[5]);
	
	ExtractReads(SAM_File, Left_Read_File, Right_Read_File, OUT_Left_Read_File, OUT_Right_Read_File);
			
	return 0;
}
