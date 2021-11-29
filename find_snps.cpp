/*
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
#include "Utilities.h"

using namespace std;

const int MINIMUM_NO_OF_READS = 5;
const int nCols = 5;

void FindSNPs(string Pileup_File, string Base_Dist_File, string OUT_File) {
	
	map < string, int** > Frequency;
	map < string, int > RefLength;
	
	ifstream ifs_dist ( Base_Dist_File.c_str() );
	string str = "";
	string reference;
	string PreviousReference = "";
	int** RefFrequency;
	while (ifs_dist.good()) {
		getline(ifs_dist, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			int length = atoi(strLength.c_str());
			
			// memory initialisation
			//int **RefFrequency;
			RefFrequency = new int* [length];
			for(int i = 0; i < length; i++) {
				*(RefFrequency+i) = new int[nCols];		
				for (int j = 0; j < nCols; j++) {
					(*(RefFrequency+i))[j] = 0;
				}
			}
			
			RefLength[reference] = length;
			Frequency[reference] = RefFrequency;
		} else {
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[0];
			int pos = atoi(entries[1].c_str());
	
			
			if (PreviousReference.compare(reference) != 0) {
				RefFrequency = Frequency[reference];
				PreviousReference = reference;
			}
			int* PosFreqency = *(RefFrequency + pos - 1);
			
			for ( int i = 3; i < entries.size(); i++ ) {
				PosFreqency[i-3] = atoi(entries[i].c_str());
			}
			
		}
	} // end while
	
	ifs_dist.close();

	// open the snp file stream
	ifstream ifs_pileup ( Pileup_File.c_str() );
	// and the stream for output
	ofstream ofs (OUT_File.c_str());
	
	string strPileup = "";
	while (ifs_pileup.good()) {
		getline(ifs_pileup, strPileup);
		
		if (strPileup.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strPileup, entries, "\t");
		
		reference = entries[0];
		int pos = atoi(entries[1].c_str());

		if (PreviousReference.compare(reference) != 0) {
			RefFrequency = Frequency[reference];
			PreviousReference = reference;
		}
		int* PosFreqency = *(RefFrequency + pos - 1);
	
		string RefBase = entries[2];
		string ConsensusBase = entries[3];
		
		if (RefBase.compare(ConsensusBase) == 0) { // same 
			continue;
		}
		
		if (RefBase[0] == '*') { // not a snp
			continue;
		}

		int ConsensusQuality = atoi(entries[4].c_str());
		if (ConsensusQuality == 0) { // no consensus
			continue;
		}
		
		int Coverage = 0;
		for (int i = 0; i < nCols; i++) {
			Coverage += PosFreqency[i];
		}
		if (Coverage < MINIMUM_NO_OF_READS) { // low coverage 
			continue;
		}
		
		if ((double)Coverage/atof(entries[7].c_str()) < 0.8) { // to get rid of SNPs that lie in a deletion
			continue;
		}
		
		ofs << entries[0] << "\t" << entries[1] << "\t" << entries[2] << "\t" << entries[3] << "\t" << entries[4] << "\t" << entries[5] << "\t" << entries[6] << "\t" << entries[7] << endl;
		
	} // end while

	// close the file
	ofs.close();
	ifs_pileup.close();
	
}


int main (int argc, char** argv) {
	
	if (argc < 4) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Pileup_File (argv[1]);
	string Base_Dist_File (argv[2]);
	string OUT_File (argv[3]);
			
	FindSNPs(Pileup_File, Base_Dist_File, OUT_File);

	return 0;
}



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
