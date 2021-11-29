/*
 *  find_indeds_from_sam_file.cpp
 *  
 *	NOTE: Run awk -F"\t" {'print $1"\t"$2"\t"$8'} <in-consensus>.pileup > | uniq > -u <out>.coverage before running the program
 *
 *  Created by Aziz Mithani on 25/03/2010.
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

using namespace std;


const int MQUAL_THRESHOLD = 20;
const int COUNT_THRESHOLD = 1;

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

void FindIndels(string SAM_File, string Coverage_File, string OUT_File_Ins, string OUT_File_Del) {
/**
	Returns the total number of reads
 */
	
	// open the sam file stream
	ifstream ifs_sam ( SAM_File.c_str() );

	string strRead="";
	//set< pair<string,int> > Insertions;
	set< pair<string, string> > Insertions;
	set< pair<string,int> > Deletions;
	map< string, int > InsCount;
	map< string, int > DelCount;
	while (ifs_sam.good()) {
		getline(ifs_sam, strRead);
		
		if (strRead.empty()) { // ignore emtry lines
			continue;
		} else if (strRead[0] == '@') { // ignore the header
			continue;
		}

		vector<string> entries;
		strsplit(strRead, entries, "\t");

		// find the occurence of "D" in 6th column (5th with 0-based indexing)
		int pos, posColon;
		if ((pos = entries[5].find("D")) == string::npos && ((pos = entries[5].find("I")) == string::npos)) {
			continue;
		} else {
			// process this read
			int ReadStart = atoi(entries[3].c_str());

			string cigar = entries[5];
			string length = "";
			int LastNonDigitPosition = -1;
			int sum = 0;
			int i = 0;
			int start = 0;
			int BasesBeforeIns = 0;
			while (i < cigar.length()) { 
				for (i=start; i<cigar.length(); i++) {
					if (isdigit(cigar[i])) {
						length += cigar[i];
					} else if (cigar[i]=='S') {
						BasesBeforeIns += atoi(length.c_str());
						// ignore
						length="";
						LastNonDigitPosition = i;
					} else if (cigar[i] =='D') {
						int DelLength = atoi(cigar.substr(LastNonDigitPosition+1, i-LastNonDigitPosition).c_str());
						LastNonDigitPosition = i;
						int DelStart = ReadStart + sum;
						start = i+1;
						sum+= DelLength;
						length="";
						
						// cout << "Deletion:" << entries[0] << "\t" << DelStart << "\t" << DelLength << endl;
						stringstream sDelStart;
						sDelStart << DelStart;
						string key = entries[2] + ":" + sDelStart.str(); 
						Deletions.insert(pair <string, int> (key, DelLength));
						// increment the count of this deletion
						stringstream sDelLength;
						sDelLength << DelLength;
						DelCount[key + ":" + sDelLength.str()] = DelCount[key + ":" + sDelLength.str()]++;
						break;
					} else if (cigar[i] =='I') {
						int InsLength = atoi(cigar.substr(LastNonDigitPosition+1, i-LastNonDigitPosition).c_str());
						LastNonDigitPosition = i;
						int InsStart = ReadStart + sum;
						start = i+1;
						//sum+= InsLength;

						// cout << "Insertion:" << entries[0] << "\t" << InsStart << "\t" << InsLength << endl;
						stringstream sInsStart;
						sInsStart << InsStart;
						string key = entries[2] + ":" + sInsStart.str();
						Insertions.insert(pair <string, string> (key, entries[9].substr(BasesBeforeIns,InsLength)));
						//cout << "Insertion:" << key << "\t" << entries[9].substr(BasesBeforeIns,InsLength) << "\t" << entries[0] << "\t" << BasesBeforeIns << "\t" << InsLength << endl;
						// increment the count of this insertion
						stringstream sInsLength;
						sInsLength << InsLength;
						InsCount[key + ":" + sInsLength.str() + ":" + entries[9].substr(BasesBeforeIns,InsLength)] = InsCount[key + ":" + sInsLength.str() + ":" + entries[9].substr(BasesBeforeIns,InsLength)]++;
						
						BasesBeforeIns += atoi(length.c_str());
						length="";
		
						break;
					} else if (cigar[i]=='M' || cigar[i]=='N') {
						BasesBeforeIns += atoi(length.c_str());
						sum += atoi(length.c_str());
						length="";
						LastNonDigitPosition = i;
					}
				}
			} // end while
			
		}
	}
	
	
	map< string,int > Coverage;
	if (!Coverage_File.empty()) {
		// open the Coverage file stream
		ifstream ifs_Coverage ( Coverage_File.c_str() );
		string str="";
		while (ifs_Coverage.good()) {
			getline(ifs_Coverage, str);
			
			
			if (str.empty()) { // ignore emtry lines
				continue;
			} 
			
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			stringstream sPosition;
			sPosition << entries[1];
			Coverage.insert( pair<string,int>(entries[0] + ":" + sPosition.str(), atoi(entries[2].c_str())) );
			
		}
		ifs_Coverage.close();
	}	
	
	// open the streams for output
	ofstream ofsIns (OUT_File_Ins.c_str());
	ofstream ofsDel (OUT_File_Del.c_str());

	ofsIns << "chromosome\tstart\tposition\tlength\tnumber of reads\tcoverage pecentage\tPoor Global coverage\tInserted Base(s)" << endl;	
	set< pair<string,string> >::iterator it;
	//ofsIns << "Insertions" << endl;
	for ( it=Insertions.begin() ; it != Insertions.end(); it++ ) {
		pair<string,string> entry = *it;
		
		stringstream sLength;
		sLength << entry.second.length();
		int count = InsCount[entry.first + ":" + sLength.str() + ":" + entry.second];
		//cout << entry.first << "\t" << entry.second << count << endl;
		if (count <= COUNT_THRESHOLD )
			continue;
		
		if (Coverage.empty()) {
			ofsIns << entry.first << "\t" << entry.second << "\t" << count  << endl;
		} else {
			double coverage; 
			coverage = (Coverage[entry.first] == 0 ? 0 : (double)count/(double)Coverage[entry.first] * 100.00);
			stringstream sCoverage;
			sCoverage << coverage;

			vector<string> entries;
            strsplit(entry.first, entries, ":");
			ofsIns << entries[0] << "\t" << entries[1] << "\t" << entry.first << "\t" << entry.second.length() << "\t" << count << "\t" << Coverage[entry.first] << "\t" << sCoverage.str() << "\t" << entry.second << endl;
		}
	}

	ofsDel << "chromosome\tstart\tposition\tlength\tnumber of reads\tcoverage pecentage\tPoor Global coverage" << endl;     
	//ofsDel << "Deletions" << endl;
	set< pair<string,int> >::iterator it1;
	for ( it1=Deletions.begin() ; it1 != Deletions.end(); it1++ ) {
		pair<string,int> entry = *it1;

		stringstream sLength;
		sLength << entry.second;
		int count = DelCount[entry.first + ":" + sLength.str()];
		if (count <= COUNT_THRESHOLD )
			continue;
		
		if (Coverage.empty()) {
			ofsDel << entry.first << "\t" << entry.second << "\t" << count  << endl;
		} else {
			double coverage; 
			coverage = (Coverage[entry.first] == 0 ? 0 : (double)count/(double)Coverage[entry.first] * 100.00);
			stringstream sCoverage;
			sCoverage << coverage;

			vector<string> entries;
            strsplit(entry.first, entries, ":");
            ofsDel << entries[0] << "\t" << entries[1] << "\t" << entry.first << "\t" << entry.second << "\t" << count << "\t" << Coverage[entry.first] << "\t" << sCoverage.str() << endl;
		}
	}
	// close the files
	ifs_sam.close();
	ofsIns.close();
	ofsDel.close();

}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string Coverage_File;
	if (argc < 3) {
		string Coverage_File  = "";
	} else {
		Coverage_File  = argv[2];
	}
	string SAM_File (argv[1]);
	
	string OUT_File_Ins = SAM_File.substr(0,SAM_File.find_last_of('.')) + ".ins";
	string OUT_File_Del = SAM_File.substr(0,SAM_File.find_last_of('.')) + ".del";
	FindIndels(SAM_File, Coverage_File, OUT_File_Ins, OUT_File_Del);
			
	return 0;
}
