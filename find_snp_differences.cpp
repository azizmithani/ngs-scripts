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

using namespace std;

const int MINIMUM_NO_OF_READS = 0;
const int MAXIMUM_NO_OF_READS = 100000;

const int nCols = 5;

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

int FindDifferences(string File_1, string File_2, string OUT_File_Diff) {
	/**
	 Reads the two tab-seperated file and performs set minus File_1 - File_2 and writes the output 
	 */
	
	set<string> Nucleotides;
	Nucleotides.insert("A");
	Nucleotides.insert("C");
	Nucleotides.insert("G");
	Nucleotides.insert("T");
/*	Nucleotides.insert("R");
	Nucleotides.insert("Y");
	Nucleotides.insert("M");
	Nucleotides.insert("K");
	Nucleotides.insert("S");
	Nucleotides.insert("W");
	Nucleotides.insert("T");
*/	
	set<string> v_1 = ReadFile(File_1);
	set<string> v_2 = ReadFile(File_2);

	map< string, vector<string> > m_1 = ReadFileCompletely(File_1);
	
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
/*		vector<string> tokens;
		strsplit(entry,tokens,":");
		string key = tokens[0] + ":" + tokens[1]; 
*/
		string key = entry;
		vector<string> entries = m_1[key];
		
		int number_of_reads = atoi(entries[7].c_str());
		if (atoi(entries[5].c_str()) == 0) { // ignore snp calls with 0 snp-quality
			continue;
		} else if (Nucleotides.find(entries[2]) == Nucleotides.end()) { // ignore if reference is ambiguous
			continue;
		} else if (Nucleotides.find(entries[3]) == Nucleotides.end()) { // ignore if consensus is ambiguous
			continue;
		}
		if (number_of_reads >= MINIMUM_NO_OF_READS && number_of_reads <= MAXIMUM_NO_OF_READS){ // ignore the ones with single read
			// output: Chr:position, chr, position, ref base, consensus base, consensus quality, SNP quality, no. of reads
			// NOTES:	Consensus quality is the Phred-scaled probability that the consensus is wrong. 
			//			SNP quality is the Phred-scaled probability that the consensus is identical to the reference. 
			//			For SNP calling, SNP quality is of more importance.
			ofsDiff << entries[0] << "\t" << entries[1] << "\t" << key << "\t" << entries[2] << "\t" << entries[3]  << "\t" << entries[4]  << "\t" << entries[5]  << "\t" << entries[7] << endl;
		}
	}
	
	// close the file
	ofsDiff.close();
	
	return 0;
}

int FindDifferences(string File_1, string File_2, string OUT_File_Diff, string BASE_DIST_File) {
	/**
	 Reads the two tab-seperated file and performs set minus File_1 - File_2 and writes the output 
	 */
	
	set<string> Nucleotides;
	Nucleotides.insert("A");
	Nucleotides.insert("C");
	Nucleotides.insert("G");
	Nucleotides.insert("T");
/*	Nucleotides.insert("R");
	Nucleotides.insert("Y");
	Nucleotides.insert("M");
	Nucleotides.insert("K");
	Nucleotides.insert("S");
	Nucleotides.insert("W");
	Nucleotides.insert("T");
*/	
	
	map < string, int** > Frequency;
	map < string, int > RefLength;
	
	ifstream ifs_dist ( BASE_DIST_File.c_str() );
	string str = "";
	while (ifs_dist.good()) {
		getline(ifs_dist, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			int length = atoi(strLength.c_str());
			
			// memory initialisation
			int **RefFrequency;
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
			
			int** RefFrequency = Frequency[reference];
			int* PosFreqency = *(RefFrequency + pos - 1);
			
			for ( int i = 3; i < entries.size(); i++ ) {
				PosFreqency[i-3] = atoi(entries[i].c_str());
			}
			
		}
	} // end while
	
	ifs_dist.close();
	
	set<string> v_1 = ReadFile(File_1);
	set<string> v_2 = ReadFile(File_2);
	
	map< string, vector<string> > m_1 = ReadFileCompletely(File_1);
	
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
		/*		vector<string> tokens;
		 strsplit(entry,tokens,":");
		 string key = tokens[0] + ":" + tokens[1]; 
		 */
		string key = entry;
		vector<string> entries = m_1[key];
		
		int number_of_reads = atoi(entries[7].c_str());
		if (atoi(entries[5].c_str()) == 0) { // ignore snp calls with 0 snp-quality
			continue;
		} else if (Nucleotides.find(entries[2]) == Nucleotides.end()) { // ignore if reference is ambiguous
			continue;
		} else if (Nucleotides.find(entries[3]) == Nucleotides.end()) { // ignore if consensus is ambiguous
			continue;
		}
		if (number_of_reads >= MINIMUM_NO_OF_READS && number_of_reads <= MAXIMUM_NO_OF_READS){ // ignore the ones which are less than the cut-off
			// output: Chr:position, chr, position, ref base, consensus base, consensus quality, SNP quality, no. of reads
			// NOTES:	Consensus quality is the Phred-scaled probability that the consensus is wrong. 
			//			SNP quality is the Phred-scaled probability that the consensus is identical to the reference. 
			//			For SNP calling, SNP quality is of more importance.

			if (atoi(entries[1].c_str()) > RefLength[entries[0]] || number_of_reads == 0) { // ideally number of reads should never be zero
			/*	for (int i = 0; i < nCols; i++) {
					ofsDiff << "\t" << 0;
				}
				for (int i = 0; i < nCols; i++) {
					ofsDiff << "\t" << 0;
				}
			 */
				continue;
			} else {
				int** RefFrequency = Frequency[entries[0]];			
				int* PosFreqency = *(RefFrequency + atoi(entries[1].c_str()) - 1);
				int coverage = 0;
				for (int i = 0; i < nCols; i++) {
					coverage += PosFreqency[i];
				}
				if (coverage == 0) {
					continue;
				}
				ofsDiff << entries[0] << "\t" << entries[1] << "\t" << key << "\t" << entries[2] << "\t" << entries[3]  << "\t" << entries[4]  << "\t" << entries[5]  << "\t" << entries[7];
				for (int i = 0; i < nCols; i++) {
					ofsDiff << "\t" << (coverage == 0 ? 0 : (double) PosFreqency[i]/ (double)coverage * 100);
				}
				for (int i = 0; i < nCols; i++) {
					ofsDiff << "\t" << PosFreqency[i];
				}
				ofsDiff << endl;
			}
		}
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
	string File_Base_Dist = "";
	if (argc >= 5) {
		File_Base_Dist = argv[4];
	}
	
	if (File_Base_Dist.length() == 0)
		FindDifferences(File_1, File_2, File_out);
	else
		FindDifferences(File_1, File_2, File_out, File_Base_Dist);

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
