/*
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

void ReadGeneInfo(string GeneInfoFile, map<string, vector<string> >& GeneInfo){

	string str="";
	// open the sam file stream
	ifstream ifs ( GeneInfoFile.c_str() );
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		GeneInfo[entries[0]] = entries;
	}
}


int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "File not specified" << endl;
		exit(0);
	}
	
	string GeneInfoFile (argv[1]);
	string GFFFile (argv[2]);
	
	// Read the gene info
	map<string, vector<string> > GeneInfo;
	ReadGeneInfo(GeneInfoFile, GeneInfo);
	
	// rename the gff file to keep a backup copy
	string BackupGFFFile = (GFFFile + ".old");
	rename (GFFFile.c_str(), BackupGFFFile.c_str()); 
	
	string str="";
	// open the file stream
	ifstream ifs ( BackupGFFFile.c_str() );
	// open the out file stream
	ofstream ofs ( GFFFile.c_str() );
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str[0] == '#') { // ignore emtry lines
			ofs << str	<< endl;
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		string attributes = entries[8];

		string UG = "";
		int posUG = attributes.find("/ug");
		int endPosUG;
		if (posUG != string::npos) {
			endPosUG = attributes.find(" ", posUG);
			UG = attributes.substr(posUG+4, endPosUG - posUG - 4);
		}
		
		vector<string> info = GeneInfo[UG];
		
		int posID = attributes.find("ID=");
		int endPosID;
		if (posID != string::npos) {
			endPosID = attributes.find(";", posID);
		}
		ofs << entries[0] << "\t" << entries[1] << "\t" << entries[2] << "\t" << entries[3] << "\t" << entries[4] << "\t" << entries[5] << "\t" << entries[6] << "\t" << entries[7] << "\t";
		ofs << attributes.substr(0,endPosID) << ";Name=" << info[1] << ";description=" << info[2] << ";" << entries[8].substr(endPosID+1) << endl;
	} // while ifs.close()
	
	ifs.close();
	ofs.close();
	
	//rename (BackupGFFFile.c_str(), GFFFile.c_str());

	return 0;
}

