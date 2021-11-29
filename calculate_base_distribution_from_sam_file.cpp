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


const int MQUAL_THRESHOLD = 20;
const bool IGNORE_UNPAIRED_READS = true;
const bool IGNORE_NON_UNIQUE_READS = true;

const int FLAG_PAIRED_READ = 2;
const int nCols = 5;


void ProcessSAMFile(string SAM_File, string OUT_File) {

	// open the sam file stream
	ifstream ifs ( SAM_File.c_str() );
	
	// read the header and initialise the variable to hold frequencies
	string str="";
	map < string, int** > Frequency;
	map < string, int > RefLength;
	while (ifs.good()) {
		getline(ifs, str);
		
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
			break;
		}

	}
	
	// open the sam file stream again
	ifstream ifs_sam ( SAM_File.c_str() );

	// also open the output file
	ofstream ofs (OUT_File.c_str());
	
	string strRead="";
	while (ifs_sam.good()) {
		getline(ifs_sam, strRead);
		
		if (strRead.empty()) { // ignore emtry lines
			continue;
		} else if (strRead[0] == '@') { //  header
			if (strRead.find("@SQ") != string::npos) { // get sequence length from the header
				ofs << strRead << endl;
			}
			//cout << strRead << endl;
			continue;
		}
		
		//cout << strRead << endl;

		vector<string> entries;
		strsplit(strRead, entries, "\t");

		string reference = entries[2];
		int start = atoi(entries[3].c_str());
		string read = entries[9];

		//cout << entries[0] << endl;

		// find the occurence of "D","I" or "S" in 6th column (5th with 0-based indexing)
		int pos, posColon;
		if ((pos = entries[5].find("D")) != string::npos || (pos = entries[5].find("I")) != string::npos || (pos = entries[5].find("S")) != string::npos) {
			// process this read
			
			string NewRead = "";
			string cigar = entries[5];
			string length = "";
			int LastNonDigitPosition = -1;
			int sum = 0;
			int i = 0;
			int start = 0;
			int TotalDeleted = 0;
			while (i < cigar.length()) { 
				for (i=start; i<cigar.length(); i++) {
					if (isdigit(cigar[i])) {
						length += cigar[i];
					} else if (cigar[i]=='S') {
						if (i == cigar.length() - 1) {
							for (int j = read.length() - 1; j >= read.length() - atoi(length.c_str()); j--) {
								read[j] = '-';
							}
						} else
							read = read.substr(atoi(length.c_str()));//,read.length());
						length="";
						LastNonDigitPosition = i;
					} else if (cigar[i] =='D') {
						int DelLength = atoi(cigar.substr(LastNonDigitPosition+1, i-LastNonDigitPosition).c_str());
						LastNonDigitPosition = i;
						int DelStart = sum;
						start = i+1;
						length="";
						//sum+= DelLength;
						
						string gap(DelLength,'.');
						NewRead += read.substr(0,DelStart) + gap;
						//cout  << "\t" << NewRead << endl;
						//cout << "\t" << DelStart << ", " << read.length() << ", " << TotalDeleted << "\t";
						read = read.substr(DelStart,read.length());
						TotalDeleted += DelStart;
						sum  = 0;
						// cout << "Deletion:" << entries[0] << "\t" << DelStart << "\t" << DelLength << endl;
						break;
					} else if (cigar[i] =='I') {
						int InsLength = atoi(cigar.substr(LastNonDigitPosition+1, i-LastNonDigitPosition).c_str());
						LastNonDigitPosition = i;
						int InsStart = sum;
						start = i+1;
						sum = 0;
						length="";

						NewRead += read.substr(0,InsStart);
						read = read.substr(InsStart + InsLength,read.length());
						TotalDeleted += InsStart + InsLength;
						
						break;
					} else if (cigar[i]=='M' || cigar[i]=='N') {
						sum += atoi(length.c_str());
						length="";
						LastNonDigitPosition = i;
					}
				}
			} // end while
			NewRead += read.substr(0,read.length());
		
			read = NewRead;
		}
		
		//cout << read << "\t" << start << endl;
		int** RefFrequency = Frequency[reference];
		for (int i = 0; i < read.length(); i++) {
			if (read[i] == '-') { // soft-clip .. ignore this position in this read
				continue;
			} else if (read[i] == '.') { // gap due to deletion .. ignore
				continue;
			}
			
			if (start + i > RefLength[reference])
				break;
			int* PosFreqency = *(RefFrequency + start + i - 1);
			int base;
			switch (read[i]) {
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
			} // end switch
			//cout << reference << "\t" << start+i << "\t" << read[i] << "\t" << PosFreqency[base] << endl;
			PosFreqency[base]++;
		} // end for (each charachter in the read) 
		
	} // end while
	

	map < string, int** >::iterator itFreq;
	for (itFreq = Frequency.begin(); itFreq != Frequency.end(); itFreq++) {
		int** RefFrequency = itFreq->second;
		int length = RefLength[itFreq->first];
		for (int i = 0; i < length; i++) {
			ofs << itFreq->first << "\t" << i+1 << "\t" << itFreq->first << ":" << i+1;
			int* PosFreqency = RefFrequency[i];
			int p = 0;
			for (int p = 0; p < nCols; p++) {
				ofs << "\t" << PosFreqency[p];
			}
			ofs << endl;
		}
	}
	
	// close the file
	ifs_sam.close();
	ofs.close();
 
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string SAM_File (argv[1]);
	string OUT_File (argv[2]);
	
	ProcessSAMFile(SAM_File, OUT_File);
			
	return 0;
}


/*
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include "Utilities.h"

using namespace std;


const int MQUAL_THRESHOLD = 20;
const bool IGNORE_UNPAIRED_READS = true;
const bool IGNORE_NON_UNIQUE_READS = true;

const int FLAG_PAIRED_READ = 2;

void ProcessSAMFile(string SAM_File, string OUT_File) {
	
	// open the sam file stream
	ifstream ifs ( SAM_File.c_str() );
				
	// open the sam file stream again
	ifstream ifs_sam ( SAM_File.c_str() );
	
	string strRead="";
	map < string, map<int, vector<int> > > Frequency;
	while (ifs_sam.good()) {
		getline(ifs_sam, strRead);
		
		if (strRead.empty()) { // ignore emtry lines
			continue;
		} else if (strRead[0] == '@') { // ignore the header
			//cout << strRead << endl;
			continue;
		}
		
		vector<string> entries;
		strsplit(strRead, entries, "\t");
		
		cout << entries[0] << endl;
		
		string reference = entries[2];
		int start = atoi(entries[3].c_str());
		string read = entries[9];
		map < string, map<int, vector<int> > >::iterator itFreq = Frequency.find(reference);
		if ( itFreq == Frequency.end() ) {
			map<int, vector<int> > RefFrequency;
			Frequency[reference] = RefFrequency;
			itFreq = Frequency.find(reference);
		}		
		map<int, vector<int> > RefFrequency = itFreq->second;
		for (int i = 0; i < read.length(); i++) {
			map<int, vector<int> >::iterator itRefFreq = RefFrequency.find(start+i);
			if ( itRefFreq == RefFrequency.end() ) {
				vector<int> PosFreqency(5,0);
				RefFrequency[start+i] = PosFreqency;
				itRefFreq = RefFrequency.find(start+i);
			}
			vector<int> PosFreqency = itRefFreq->second;
			int base;
			switch (read[i]) {
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
			} // end switch
			//cout << reference << "\t" << start+i << "\t" << read[i] << "\t" << PosFreqency[base] << endl;
			PosFreqency[base]++;
			RefFrequency[start+i] = PosFreqency;
			Frequency[reference] = RefFrequency;
		} // end for (each charachter in the read) 
		
	} // end while
	
	ofstream ofs (OUT_File.c_str());
	
	map < string, map<int, vector<int> > >::iterator itFreq;
	for (itFreq = Frequency.begin(); itFreq != Frequency.end(); itFreq++) {
		map<int, vector<int> > RefFrequency = itFreq->second;
		map<int, vector<int> >::iterator itRefFreq;
		for (itRefFreq = RefFrequency.begin(); itRefFreq != RefFrequency.end(); itRefFreq++) {
			ofs << itFreq->first << "\t" << itRefFreq->first << "\t" << itFreq->first << ":" << itRefFreq->first;
			vector<int> PosFreqency = itRefFreq->second;
			vector<int>::iterator itPosFreq;
			for (itPosFreq = PosFreqency.begin(); itPosFreq != PosFreqency.end(); itPosFreq++) {
				ofs << "\t" << *itPosFreq;
			}
			ofs << endl;
		}
	}
	
	// close the file
	ifs_sam.close();
	ofs.close();
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string SAM_File (argv[1]);
	string OUT_File (argv[2]);
	
	ProcessSAMFile(SAM_File, OUT_File);
	
	return 0;
}
*/
