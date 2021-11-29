/*
 *  fix_mutations.cpp
 *  
 *	fixes the supplied mutation in the given FASTA file
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
#include <algorithm>
#include <map>
#include "Utilities.h"

using namespace std;

const int BLOCK_SIZE = 80;

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

void FixSNPs(string FASTA_File, string SNP_File, string OUT_File) {
	
	cout << "Fixing SNPs ... " << endl;
	map < string, map <int, string> > SNPList;
	
	// open the file streams
	ifstream ifs_snp ( SNP_File.c_str() );
	
	// go through the list of snp file to save snps each referecne (chromosome)
	int offset = 0;
	string strSNP = "";
	string reference = "";
	string currentReference = ""; // reference (chromosome) being currently read 
	map<int, string> RefSNPList;
	while (ifs_snp.good()) {
		
		// read the snp 
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
			offset = 1;
		} else {
			offset = 0;
		}
		
		reference = entries[0];
		
		if (reference.compare(currentReference) != 0) { // new reference
			if (currentReference.length() > 0) {
				SNPList[currentReference] = RefSNPList;
			}
			RefSNPList = SNPList[reference];
			currentReference = reference;
		}
		
		int SNPPos = atoi(entries[1].c_str());
		RefSNPList[SNPPos] = entries[3+offset];
		
		
	} // end while
	SNPList[currentReference] = RefSNPList;
	
	ifs_snp.close();
	
	// open the output file
	ofstream ofs (OUT_File.c_str());
	
	ifstream ifs_fasta ( FASTA_File.c_str() );
	string strFASTA = "";
	int fastaPos = 0; // current position in the fasta file
	int fastaTotalReadPos = 0; // position until which the file has been read
	bool refNotFound = false;
	reference = "";
	currentReference = "";
	map< int, string >::iterator itSNP;
	while (ifs_fasta.good()) {
		// read the next line from the fasta file
		getline(ifs_fasta, strFASTA);
		
		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}
		
		// match the reference
		if (strFASTA[0] == '>') { // we are reading a new reference (chromosome)
			reference = ExtractReferenceName(strFASTA);
			
			if (reference.compare(currentReference) != 0) { // new reference
				if ( SNPList.find(reference) == SNPList.end() )
					refNotFound = true;
				else {		
					RefSNPList = SNPList[reference];
					itSNP = RefSNPList.begin();
					currentReference = reference;
					fastaPos = 0;
					fastaTotalReadPos = 0; 
					refNotFound = false;
				}
			}
			//cout << reference << "\t" << RefSNPList.size() << endl;
			ofs << strFASTA << endl;
			continue;
		}
		
		// sequence data... check if there is any SNP present for this reference. If not, output it as it is
		if (refNotFound) {
			ofs << strFASTA << endl;
			continue;
		}
		
		// sequence data... some SNPs may be present ... process them
		// the position up to which the fasta file has been read  
		fastaTotalReadPos += strFASTA.length();
		
		if (itSNP == RefSNPList.end()) { // no more snp for this reference
			ofs << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}		
		
		if (itSNP->first > fastaTotalReadPos) { // no snp in this line .. output it as it is
			ofs << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}
		
		// one or more SNPs are present in this line
		for (int i = 0; i < strFASTA.length(); i++) {
			fastaPos++;
			if (itSNP->first == fastaPos) {
				strFASTA[i] = itSNP->second[0];				
				itSNP++;
			}
		}
		ofs << strFASTA << endl;
		
	} // end while ifs_fasta.good()
	
	// close the stream
	ifs_fasta.close();
	ofs.close();
	
}

void FixDeletions(string FASTA_File, string Deletion_File, string OUT_File, map < string, map <int, int> >& DelList) {
	//map < string, int >& RefLength, map < string, int* >& PositionMap) {
	
	/* **************************** *
	 *		D E L E T I O N S		*
	 * **************************** */
	cout << "Fixing deletions ... " << endl;
	
	// open the file stream
	ifstream ifs_del ( Deletion_File.c_str() );
	
	// go through the list of deletion file to save deletions for each referecne (chromosome)
	int offset = 0;
	string strDel = "";
	string reference = "";
	string currentReference = ""; // reference (chromosome) being currently read 
	map<int, int> RefDelList;
	while (ifs_del.good()) {
		
		// read the deletion 
		getline(ifs_del, strDel);
		
		if (strDel.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strDel, entries, "\t");
		
		if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
			offset = 1;
		} else {
			offset = 0;
		}
		
		reference = entries[0];
		
		if (reference.compare(currentReference) != 0) { // new reference
			if (currentReference.length() > 0) {
				DelList[currentReference] = RefDelList;
				//ofs << currentReference << "\t" << RefDelList.size() << endl;
			}
			RefDelList = DelList[reference];
			
			currentReference = reference;
		}
		
		//RefDelList[atoi(entries[1].c_str())] = atoi(entries[2+offset].c_str());
		int position = atoi(entries[1].c_str()) + 1;
		int length = entries[2+offset].length();
		RefDelList[position] = length;
		
	} // end while
	DelList[currentReference] = RefDelList;
	
	ifs_del.close();
	
	// open the temporary output file
	ofstream ofs_del (OUT_File.c_str());
	
	ifstream ifs_fasta_del ( FASTA_File.c_str() );
	string strFASTA = "";
	int fastaPos = 0; // current position in the fasta file
	int fastaTotalReadPos = 0; // position until which the file has been read
	bool refNotFound = false;
	reference = "";
	currentReference = "";
	map< int, int >::iterator itDel;
	int BasesToBeDeleted = 0;
	while (ifs_fasta_del.good()) {
		// read the next line from the fasta file
		getline(ifs_fasta_del, strFASTA);
		
		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}
		
		// match the reference
		if (strFASTA[0] == '>') { // we are reading a new reference (chromosome)
			reference = ExtractReferenceName(strFASTA);
			
			if (reference.compare(currentReference) != 0) { // new reference
				if ( DelList.find(reference) == DelList.end() )
					refNotFound = true;
				else {		
					RefDelList = DelList[reference];
					itDel = RefDelList.begin();
					currentReference = reference;
					fastaPos = 0;
					fastaTotalReadPos = 0; 
					refNotFound = false;
				}
			}
			//cout << reference << "\t" << RefDelList.size() << endl;
			ofs_del << strFASTA << endl;
			continue;
		}
		
		// sequence data... check if there is any deletion present for this reference. If not, output it as it is
		if (refNotFound) {
			//			cout << "Reference '" << reference << "' not found" << endl;
			ofs_del << strFASTA << endl;
			continue;
		}
		
		// sequence data... some deletions may be present ... process them
		// the position up to which the fasta file has been read  
		fastaTotalReadPos += strFASTA.length();
		
		if (BasesToBeDeleted > 0) {
			for (int i = 0; i < strFASTA.length() && BasesToBeDeleted > 0; i++) {
				strFASTA[i] = '-';								
				BasesToBeDeleted--;
			}
		}
		if (itDel == RefDelList.end()) { // no more snp for this reference
			//			cout << "No deletion found for the reference '" << reference << "'" << endl;
			ofs_del << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}		
		
		
		if (itDel->first > fastaTotalReadPos) { // no snp in this line .. output it as it is
			//			cout << "Del at " << itDel->first << "\t" << "Total bases Read: " << fastaTotalReadPos << "\t";
			//			cout << "No deletion found in the current line for the reference '" << reference << "'" << endl;
			ofs_del << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}
		
		//cout << reference << "\t" << strFASTA << endl;
		
		// one or more deletions are present in this line
		//		string strFASTANew = "";
		for (int i = 0; i < strFASTA.length(); i++) {
			fastaPos++;
			if (BasesToBeDeleted > 0) {
				strFASTA[i] = '-';								
				BasesToBeDeleted--;
			}
			if (itDel->first == fastaPos) {
				//cout << "Fixing " << itDel->second << "bp deletion at " << reference << ":" << fastaPos << endl;
				strFASTA[i] = '-';								
				BasesToBeDeleted = max(BasesToBeDeleted - 1,itDel->second - 1);
				//if (BasesToBeDeleted == 0) { // all the bases in this deletion have been deleted
				itDel++;
				//	cout << itDel->first << "\t" << fastaPos << endl;
				//}
				//continue;
			}
		}
		ofs_del << strFASTA << endl;
		
	} // end while ifs_fasta.good()
	
	// close the streams
	ifs_fasta_del.close();
	ofs_del.close();
	
}

int FixInsertions(string FASTA_File, string Insertion_File, string OUT_File, map < string, map <int, string> >& InsList) {
	//map < string, int >& RefLength, map < string, int* >& PositionMap) {
	
	/* **************************** *
	 *		I N S E R T I O N S		*
	 * **************************** */
	cout << "Fixing insertions ... " << endl;
	
	// open the file streams
	ifstream ifs_ins ( Insertion_File.c_str() );
	
	// go through the list of insertion file to save insertions for each referecne (chromosome)
	int offset = 0;
	string strIns = "";
	string reference = "";
	string currentReference = ""; // reference (chromosome) being currently read 
	map<int, string> RefInsList;
	while (ifs_ins.good()) {
		
		// read the insertion 
		getline(ifs_ins, strIns);
		
		if (strIns.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strIns, entries, "\t");
		
		if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
			offset = 1;
		} else {
			offset = 0;
		}
		
		reference = entries[0];
		
		if (reference.compare(currentReference) != 0) { // new reference
			if (currentReference.length() > 0) {
				InsList[currentReference] = RefInsList;
			}
			
			RefInsList = InsList[reference];
			
			currentReference = reference;
		}
		
		int position = atoi(entries[1].c_str());
		//RefInsList[atoi(entries[1].c_str())] = entries[6+offset];
		RefInsList[position + 1] = entries[2+offset];
		
		
	} // end while
	InsList[currentReference] = RefInsList;
	//	ofs << currentReference << "\t" << RefDelList.size() << endl;
	
	ifs_ins.close();
	
	// also open the temporary output file
	ofstream ofs_ins (OUT_File.c_str());
	
	ifstream ifs_fasta_ins ( FASTA_File.c_str() );
	string strFASTA = "";
	int fastaPos = 0; // current position in the fasta file
	int fastaTotalReadPos = 0; // position until which the file has been read
	bool refNotFound = false;
	reference = "";
	currentReference = "";
	map< int, string >::iterator itIns;
	int prevInsPos = 0;
	string prevInsBases = "";
	while (ifs_fasta_ins.good()) {
		// read the next line from the fasta file
		getline(ifs_fasta_ins, strFASTA);
		
		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}
		
		// match the reference
		if (strFASTA[0] == '>') { // we are reading a new reference (chromosome)
			reference = ExtractReferenceName(strFASTA);
			
			if (reference.compare(currentReference) != 0) { // new reference
				if ( InsList.find(reference) == InsList.end() )
					refNotFound = true;
				else {		
					RefInsList = InsList[reference];
					itIns = RefInsList.begin();
					currentReference = reference;
					fastaPos = 0;
					fastaTotalReadPos = 0; 
					refNotFound = false;
				}
			}
			//cout << reference << "\t" << RefDelList.size() << endl;
			ofs_ins << strFASTA << endl;
			continue;
		}
		
		// sequence data... check if there is any insertion present for this reference. If not, output it as it is
		if (refNotFound) {
			//			cout << "Reference '" << reference << "' not found" << endl;
			ofs_ins << strFASTA << endl;
			continue;
		}
		
		// sequence data... some insertion may be present ... process them
		// the position up to which the fasta file has been read  
		fastaTotalReadPos += strFASTA.length();
		
		if (itIns == RefInsList.end()) { // no more insertion for this reference
			ofs_ins << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}		
		
		if (itIns->first > fastaTotalReadPos) { // no insetion in this line .. output it as it is
			ofs_ins << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}
		
		// one or more insertions are present in this line
		string strFASTANew = "";
		for (int i = 0; i < strFASTA.length(); i++) {
			fastaPos++;
			if (itIns->first == fastaPos) {
				//cout << "Fixing " << itIns->second.length() << "bp insertion at " << reference << ":" << fastaPos << endl;
				if (itIns->first == prevInsPos) {
					cout << "There is some problemetic data. Aborting " << endl;
					ifs_fasta_ins.close();
					ofs_ins.close();
					return 1;
					
					if (itIns->second.length() <= prevInsBases.length()) {
						// this insertion is embeded in the previous do nothing
					} else {
						int pos = itIns->second.find_first_of(prevInsBases[prevInsBases.length()-1]);
						strFASTANew += itIns->second.substr(pos+1, itIns->second.length());
					}
					
				} else {
					//strFASTANew += itIns->second;
					for (int j = 0; j < itIns->second.length(); j++) {
						strFASTANew += itIns->second.substr(j,1) + "+";
					}
					//string tag(itIns->second.length(), '+');
					//strFASTANew += itIns->second + tag;
				}
				strFASTANew += strFASTA.substr(i,1);							
				prevInsPos = fastaPos;
				prevInsBases = itIns->second;
				//if (BasesToBeDeleted == 0) { // all the bases in this deletion have been deleted
				itIns++;
				//	cout << itDel->first << "\t" << fastaPos << endl;
				//}
				continue;
			}
			strFASTANew += strFASTA.substr(i,1);							
		}
		ofs_ins << strFASTANew << endl;
		
	} // end while ifs_fasta.good()
	
	// close the streams
	ifs_fasta_ins.close();
	ofs_ins.close();
	
}

bool FinaliseOutput(string FASTA_File, string OUT_File, map<string, map<pair<int, int>, pair<int, int> > >& GFF_New_Coord) {
	
	cout << "Finalising and writing output ..." << endl;
	//  delete the '-'s, fix the insertions ('+'s) and write to the output
	ifstream ifs_temp (FASTA_File.c_str());
	ofstream ofs_fasta (OUT_File.c_str());
	
	// update the gff file
	map<pair<int, int>, pair<int, int> > RefGFF_New_Coord;
	map<pair<int, int>, pair<int, int> >::iterator itGFF;
	int seq_start, seq_end;
	
	string reference = "";
	string strFASTA = "";
	string sequence = "";
	int indel_pos = 0;
	int TotalDeleted = 0;
	int TotalInserted = 0;
	while (ifs_temp.good()) {
		// read the next line from the fasta file
		getline(ifs_temp, strFASTA);
		
		if (strFASTA.empty()) { // ignore emtry lines
			continue;
		}
		
		if (strFASTA[0] == '>') {
			// output everything from the previous refrence
			int pos = 0;
			int SeqLength = sequence.size();
			while (pos < SeqLength) {
				ofs_fasta << sequence.substr(pos, BLOCK_SIZE) << endl;
				pos += BLOCK_SIZE;
			}
			sequence = "";
			
			if (reference.length() != 0) {
				GFF_New_Coord[reference] = RefGFF_New_Coord;
			}			
			
			vector<string> entries;
			strsplit(strFASTA, entries, "\t");
			
			reference = ExtractReferenceName(entries[0]);
			
			if (GFF_New_Coord.find(reference) == GFF_New_Coord.end()) {
				cout << "Error: Reference '" << reference << "' not found in the GFF file" << endl;
				ifs_temp.close();
				ofs_fasta.close();
				remove(OUT_File.c_str());
				return false;
			}
			
			RefGFF_New_Coord = GFF_New_Coord[reference];
			itGFF = RefGFF_New_Coord.begin();
			seq_start = itGFF->first.first;
			seq_end = itGFF->first.second;
			indel_pos = 0;
			TotalDeleted = 0;
			TotalInserted = 0;
			// output this reference name
			ofs_fasta << strFASTA << endl;
			continue;
		}
		
		if (strFASTA.find("-") == string::npos && strFASTA.find("+") == string::npos) { //no indels in this line
			sequence += strFASTA;
			indel_pos += strFASTA.length();
		} else {
			char previous_base = ' ';
			for (int i = 0; i < strFASTA.length(); i++) {
				if (strFASTA[i] != '-' && strFASTA[i] != '+') {
					sequence += strFASTA.substr(i,1);
					//cout << sequence << endl;
					indel_pos++;
				} else {
					if (strFASTA[i] == '-') { //previous_base != '-' && 
						indel_pos++;
					}
					
					//cout << previous_base << "\t" << indel_pos << "\t" << strFASTA[i] << endl;
					
					// update the gff entries
					itGFF = RefGFF_New_Coord.begin();
					for (; itGFF != RefGFF_New_Coord.end(); itGFF++) {
						if (indel_pos < itGFF->first.first) { // before this gff entry
							if (strFASTA[i] == '+') { 
								itGFF->second.first++;
								itGFF->second.second++;
							} else {
								itGFF->second.first--;
								itGFF->second.second--;
							}							
						} else if (indel_pos > itGFF->first.second) { //after this gff entry
							// do nothing
						} else  {
							if (strFASTA[i] == '+') { 
								itGFF->second.second++;
							} else {
								itGFF->second.second--;
							}							
						}
					}
					
					if (strFASTA[i] == '+') {
						indel_pos--;
					}
					
					if (previous_base == '-' && strFASTA[i] == '-') {
						//			indel_pos++;
					}
					if (previous_base == '+' && strFASTA[i] == '+') {
						//			indel_pos--;
					}
					
				} // end if no indel at this position
				
				previous_base = strFASTA[i];
			}
		}
		
		if (sequence.length() == BLOCK_SIZE) {
			ofs_fasta << sequence << endl;
			sequence = "";
		} else if (sequence.length() > BLOCK_SIZE) {
			ofs_fasta << sequence.substr(0, BLOCK_SIZE) << endl;
			sequence = sequence.substr(BLOCK_SIZE,sequence.length());
		}
	}
	
	GFF_New_Coord[reference] = RefGFF_New_Coord;
	
	int pos = 0;
	int SeqLength = sequence.size();
	while (pos < SeqLength) {
		ofs_fasta << sequence.substr(pos, BLOCK_SIZE) << endl;
		pos += BLOCK_SIZE;
	}
	sequence = "";
	
	
	// close the streams
	ifs_temp.close();
	ofs_fasta.close();
	
	return true;
	
}

void WriteGFFFile(string OUT_GFF_File, vector< vector<string> > GFF_Data, map<string, map<pair<int, int>, pair<int, int> > > GFF_New_Coord) {
	
	cout << "Writing new GFF file ..." << endl;
	
	ofstream ofs_gff (OUT_GFF_File.c_str());
	ofs_gff << "##gff-version 3" << endl;
	
	string reference = "";
	string previous_reference = "";
	vector< vector<string> >::iterator itGFF = GFF_Data.begin();
	map< pair<int,int>, pair<int, int> > RefGFF_New_Coord;
	for (; itGFF != GFF_Data.end(); itGFF++) {
		vector<string> entries = *itGFF;
		string reference = entries[0];
		int start = atoi(entries[3].c_str());
		int end = atoi(entries[4].c_str());
		
		if (reference.compare(previous_reference) != 0) {
			RefGFF_New_Coord = GFF_New_Coord[reference];
			
			previous_reference = reference;
		}
		
		pair<int, int> new_coord = RefGFF_New_Coord[pair<int, int> (start,end)];
		
		string attributes = entries[8];
		int posOrigCoord = attributes.find("original_coordinates");
		if (posOrigCoord == string::npos && (start != new_coord.first || end != new_coord.second) ) {//only add the original coordinates if not present previously
			if (attributes[attributes.length() - 1] != ';') {
				attributes += ";";
			}
			attributes += "original_coordinates=(" + entries[3] + "," + entries[4] + ")";
		}
		ofs_gff << reference << "\t" << entries[1] << "\t" << entries[2] << "\t" << new_coord.first << "\t" << new_coord.second << "\t" << entries[5] << "\t" << entries[6] << "\t" << entries[7] << "\t" << attributes	<< endl;
	}
	
	ofs_gff.close();
}

void CalculateRefLengths(string FASTA_File, map < string, int >& RefLength) {
	
	ifstream ifs( FASTA_File.c_str() );
	string str = "";
	string reference = "";
	int length = 0;
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str[0] == '>') { // we are reading a new reference (chromosome)
			if (reference.length() > 0) {
				RefLength[reference] = length;
			}
			
			reference = ExtractReferenceName(str);
			length = 0;
		} else {
			length += str.length();			
		}
	} // end while
	
	// save the last reference
	RefLength[reference] = length;
	
	ifs.close();
	
}

void InitialisePositionMap(map < string, int >& RefLength, map< string, int* >& PositionMap) {
	
	map<string, int>::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		string reference = itRefLength->first;
		int length = itRefLength->second;
		
		int *RefPositionMap = new int[length + 1];
		for(int i = 0; i <= length; i++) {
			*(RefPositionMap+i) = i;
		}
		PositionMap[reference] = RefPositionMap;
	}
	
}


void ReadPositionMapFile(string Map_File, map< string, int* >& PositionMap) {
	
	cout << "Reading position map file ..." << endl;

	ifstream ifs( Map_File.c_str() );
	
	string reference = "";
	string previous_reference = "";
	int *RefPositionMap;
	string str = "";
	int Position = 0;
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else {
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[0];
			
			if (reference.compare(previous_reference) != 0) {
				if (previous_reference.length() > 0) {
					PositionMap[previous_reference] = RefPositionMap;
				}
				
				RefPositionMap = PositionMap[reference];
				Position  = 0;
				
				previous_reference = reference;
			}
			
			//int CurrentPosition = atoi(entries[1].c_str());
			int InitialPosition = atoi(entries[2].c_str());
			Position++;
			*(RefPositionMap + Position) = InitialPosition;
				
		}
	} // end while
	PositionMap[reference] = RefPositionMap;
	
	ifs.close();
	
}

//void UpdatePositionMap(map < string, int >& RefLength, map < string, int >& NewRefLength, map < string, int* >& PositionMap,  map < string, map <int, int> >& DelList,  map < string, map <int, string> >& InsList) {
void UpdatePositionMap(string OUT_MAP_File, map < string, int >& RefLength, map < string, int* >& PositionMap,  map < string, map <int, int> >& DelList,  map < string, map <int, string> >& InsList) {
	
	cout << "Writing new position map file ..." << endl;
	
	ofstream ofs (OUT_MAP_File.c_str());
	
	map<string, int>::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		string reference = itRefLength->first;
		int ReferenceLength = itRefLength->second;
		//int NewReferenceLength = NewRefLength[reference];
		//int ReferenceLength = min(OldReferenceLength, NewReferenceLength);

		int* RefPositionMap = PositionMap[reference];
		
		// deletion
		int previous_position = 0;
		//int nDeletedBases = 0;		
		map<int, int> RefDelList = DelList[reference];
		map<int, int>::iterator itRefDelList = RefDelList.begin();
		for (; itRefDelList != RefDelList.end(); itRefDelList++) {
			int position = itRefDelList->first;
			int length = itRefDelList->second;
			
			// mark the deleted bases
			for (int i = 0; i < length; i++) {
				*(RefPositionMap + position + i) = -1;
			}
			
			previous_position = position + length - 1;
			
		}
		
		
		
		// insertion
		previous_position = 0;
		//int nInsertedBases = 0;
		map<int, string> RefInsList = InsList[reference];
		map<int, string>::iterator itRefInsList = RefInsList.begin();
		int NewPosition = 0;
		for (; itRefInsList != RefInsList.end(); itRefInsList++) {
			int position = itRefInsList->first;
			int length = itRefInsList->second.length();
			
			for (int i = previous_position + 1; i < position; i++) {
				if ( *(RefPositionMap+i) != -1 ) {
					//*(RefPositionMap+i) += nInsertedBases;
					ofs << reference << "\t" << ++NewPosition << "\t" << *(RefPositionMap+i) << endl;
				}
			}
			for (int i = 0; i < length; i++) {
				//*(RefPositionMap+position+i) = 0;
				ofs << reference << "\t" << ++NewPosition << "\t" << 0 << endl;
			}
			
			//nInsertedBases += length;
			previous_position = position - 1; // subtract 1 because we added 1 when building the insertion list
		}			
		// process the data from last insertion for the current reference
		int i = previous_position + 1;
		while (NewPosition < ReferenceLength) {
			if ( *(RefPositionMap+i) != -1 ) {
				//*(RefPositionMap+i) -= nInsertedBases;
				ofs << reference << "\t" << ++NewPosition << "\t" << *(RefPositionMap+i) << endl;
			}
			i++;
		}
//		for (int i = previous_position + 1; i <= ReferenceLength; i++) {
//		}
		
/*		if (NewReferenceLength > OldReferenceLength) {
			int Position = *(RefPositionMap + OldReferenceLength);
			for (int i = OldReferenceLength + 1; i <= NewReferenceLength; i++) {
				ofs << reference << "\t" << ++NewPosition << "\t" << ++Position << endl;
			}
			
		}
*/		
//		PositionMap[reference] = RefPositionMap;
	}
	ofs.close();
	
}

int main (int argc, char** argv) {
	
	if (argc < 6) {
		cout <<  "Files not specified" << endl;
		exit(0);
	}
	string FASTA_File (argv[1]);
	string OUT_File (argv[2]);
	string SNP_File (argv[3]);
	string Deletion_File (argv[4]);
	string Insertion_File (argv[5]);
	
	map < string, int > RefLength;
	map < string, int* > PositionMap;
	CalculateRefLengths(FASTA_File, RefLength);
	InitialisePositionMap(RefLength, PositionMap);
	
	string Map_File = FASTA_File + ".map";
	if (FileExists(Map_File)) {
		ReadPositionMapFile(Map_File, PositionMap);
	}
	
	// gff file
	ifstream ifs_gff ((FASTA_File + ".gff").c_str());
	string strGFF = "";
	vector< vector<string> > GFF_Data;
	map<string, map<pair<int, int>, pair<int, int> > > GFF_New_Coord;
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		}
		if (strGFF[0] == '#') { // ignore header
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		GFF_Data.push_back(entries);
		int start = atoi(entries[3].c_str());
		int end = atoi(entries[4].c_str());
		map< string, map< pair<int,int>, pair<int, int> > >::iterator itRefGFF_New_Coord = GFF_New_Coord.find(entries[0]);
		if (itRefGFF_New_Coord == GFF_New_Coord.end()) {
			GFF_New_Coord[entries[0]];
			itRefGFF_New_Coord = GFF_New_Coord.find(entries[0]);
		}
		itRefGFF_New_Coord->second.insert(pair< pair<int, int>, pair<int, int> > ( pair<int, int>(start,end), pair<int, int>(start,end) ));
		//GFF_New_Coord[entries[0]] = RefGFF_New_Coord;
	} // end while ifs_gff is good
	
	string Temp_File;
	map < string, map <int, int> > DelList;
	map < string, map <int, string> > InsList;
	
	if (SNP_File.length() > 0) {
		Temp_File = OUT_File + ".stmp";
		FixSNPs(FASTA_File, SNP_File, Temp_File);
		FASTA_File = Temp_File;
	}
	if (Deletion_File.length() > 0) {
		Temp_File = OUT_File + ".dtmp";
		FixDeletions(FASTA_File, Deletion_File, Temp_File, DelList);
		FASTA_File = Temp_File;
	}
	if (Insertion_File.length() > 0) {
		Temp_File = OUT_File + ".itmp";
		FixInsertions(FASTA_File, Insertion_File, Temp_File, InsList);
		FASTA_File = Temp_File;
	}
	if ( FinaliseOutput(FASTA_File, OUT_File, GFF_New_Coord) ) {
		
		WriteGFFFile((OUT_File + ".gff"), GFF_Data, GFF_New_Coord);
				
		map < string, int > NewRefLength;
		CalculateRefLengths(OUT_File, NewRefLength);

		UpdatePositionMap((OUT_File + ".map"), NewRefLength, PositionMap, DelList, InsList);
		
	}
	
	remove((OUT_File + ".stmp").c_str());
	remove((OUT_File + ".dtmp").c_str());
	remove((OUT_File + ".itmp").c_str());
	
	return 0;
}
