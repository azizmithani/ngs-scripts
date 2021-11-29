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
#include "Utilities.h"

using namespace std;

//const int MINIMUM_NO_OF_READS = 7;
//const int MAXIMUM_NO_OF_READS = 75;

const int nCols = 5;

map <char, string> NucleotideMap;
/*
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
		
		//v.insert (entries[0] + ":" + entries[1] + ":" + entries[2] + ":" + entries[3]);
		v.insert (entries[0] + ":" + entries[1]);
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
*/
void FindMinimumAllele(string SNP_File, string Base_Dist_File, string OUT_File) {
	
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
	ifstream ifs_snp ( SNP_File.c_str() );
	// and the stream for output
	ofstream ofs (OUT_File.c_str());
	
	int offset = -1;
	string strSNP = "";
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		reference = entries[0];
		int pos = atoi(entries[1].c_str());
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		if (PreviousReference.compare(reference) != 0) {
			RefFrequency = Frequency[reference];
			PreviousReference = reference;
		}
		int* PosFreqency = *(RefFrequency + pos - 1);
	
		int Coverage = 0;
		for (int i = 0; i < nCols; i++) {
			Coverage += PosFreqency[i];
		}
		
		string SNPBase = entries[3+offset];
		string strBase = NucleotideMap[SNPBase[0]];
		char min_freq_base;
		int min_freq = INT_MAX;
		for (int i = 0; i < strBase.length(); i++) {			
			int base;
			switch (strBase[i]) {
				case 'A':
					base = 0;
					break;
				case 'C':
					base = 1;
					break;
				case 'G':
					base = 2;
					break;
				case 'T':
					base = 3;
					break;
				default:
					base = 4;
					break;
			}
			
			if (PosFreqency[base] < min_freq) {
				min_freq = PosFreqency[base];
				min_freq_base = strBase[i];
			}
		} // end for
		
		ofs << strSNP << "\t" << min_freq_base << "\t" << min_freq << "\t" << Coverage << endl;
		
	} // end while
	

	// close the file
	ofs.close();
	ifs_snp.close();
	
}


int main (int argc, char** argv) {
	
	if (argc < 4) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string SNP_File (argv[1]);
	string Base_Dist_File (argv[2]);
	string OUT_File (argv[3]);
	
	NucleotideMap['A'] = "A";
	NucleotideMap['C'] = "C";
	NucleotideMap['G'] = "G";
	NucleotideMap['T'] = "T";
	NucleotideMap['R'] = "AG";
	NucleotideMap['Y'] = "CT";
	NucleotideMap['M'] = "AC";
	NucleotideMap['K'] = "GT";
	NucleotideMap['S'] = "CG";
	NucleotideMap['W'] = "AT";
	NucleotideMap['H'] = "ACT";
	NucleotideMap['B'] = "CGT";
	NucleotideMap['D'] = "AGT";
	NucleotideMap['V'] = "ACG";
	NucleotideMap['N'] = "N";
		
	FindMinimumAllele(SNP_File, Base_Dist_File, OUT_File);

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
