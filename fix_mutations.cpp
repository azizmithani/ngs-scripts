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
				//ofs << currentReference << "\t" << RefSNPList.size() << endl;
			}
			RefSNPList = SNPList[reference];
			currentReference = reference;
		}
		
		RefSNPList[atoi(entries[1].c_str())] = entries[3+offset];
	} // end while
	SNPList[currentReference] = RefSNPList;
//	ofs << currentReference << "\t" << RefSNPList.size() << endl;
	
	ifs_snp.close();
	
	// also open the output file
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
			reference = strFASTA.substr(1,strFASTA.length());
			
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
//			cout << "Reference '" << reference << "' not found" << endl;
			ofs << strFASTA << endl;
			continue;
		}
		
		// sequence data... some SNPs may be present ... process them
		// the position up to which the fasta file has been read  
		fastaTotalReadPos += strFASTA.length();

		if (itSNP == RefSNPList.end()) { // no more snp for this reference
//			cout << "No SNP found for the reference '" << reference << "'" << endl;
			ofs << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}		

		if (itSNP->first > fastaTotalReadPos) { // no snp in this line .. output it as it is
//			cout << "SNP at " << itSNP->first << "\t" << "Total bases Read: " << fastaTotalReadPos << "\t";
//			cout << "No SNP found in the current line for the reference '" << reference << "'" << endl;
			ofs << strFASTA << endl;
			fastaPos += strFASTA.length();
			continue;
		}

		//cout << reference << "\t" << strFASTA << endl;
		
		// one or more SNPs are present in this line
		for (int i = 0; i < strFASTA.length(); i++) {
			fastaPos++;
			if (itSNP->first == fastaPos) {
				//cout << "Fixing SNP at " << reference << ":" << fastaPos << "\t" << strFASTA[i] << " --> " << itSNP->second[0] << endl;
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

void FixDeletions(string FASTA_File, string Deletion_File, string OUT_File) {
	
//	string Temp_File_Del = OUT_File + ".dtmp";
//	string Temp_File_Ins = OUT_File + ".itmp";
	
	
	/* **************************** *
	 *		D E L E T I O N S		*
	 * **************************** */
	cout << "Fixing deletions ... " << endl;
	map < string, map <int, int> > DelList;
	
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
		
		RefDelList[atoi(entries[1].c_str())] = atoi(entries[2+offset].c_str());
	} // end while
	DelList[currentReference] = RefDelList;
	//	ofs << currentReference << "\t" << RefDelList.size() << endl;
	
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
			reference = strFASTA.substr(1,strFASTA.length());
			
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
				//if (BasesToBeDeleted == 0) { // all the bases in this deletion have been deleted
				//	itDel++;
				//	cout << itDel->first << "\t" << fastaPos << endl;
				//}
			//	continue;
			}
			if (itDel->first == fastaPos) {
				//cout << "Fixing " << itDel->second << "bp deletion at " << reference << ":" << fastaPos << endl;
				strFASTA[i] = '-';								
				BasesToBeDeleted = max(BasesToBeDeleted - 1,itDel->second - 1);
				//if (BasesToBeDeleted == 0) { // all the bases in this deletion have been deleted
					itDel++;
				//	cout << itDel->first << "\t" << fastaPos << endl;
				//}
				continue;
			}
		}
		ofs_del << strFASTA << endl;
				
	} // end while ifs_fasta.good()

	// close the streams
	ifs_fasta_del.close();
	ofs_del.close();

}

int FixInsertions(string FASTA_File, string Insertion_File, string OUT_File) {

	/* **************************** *
	 *		I N S E R T I O N S		*
	 * **************************** */
	cout << "Fixing insertions ... " << endl;
	map < string, map <int, string> > InsList;
	
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
		
		RefInsList[atoi(entries[1].c_str())] = entries[6+offset];
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
			reference = strFASTA.substr(1,strFASTA.length());
			
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
		
		// sequence data... some deletions may be present ... process them
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
					
				} else
					strFASTANew += itIns->second;
				
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

void FinaliseOutput(string FASTA_File, string OUT_File) {

	cout << "Writing output ..." << endl;
	//  delete the '-'s, fix the insertions and write to the output 
	ifstream ifs_temp (FASTA_File.c_str());
	ofstream ofs_fasta (OUT_File.c_str());

	string strFASTA = "";
	string sequence = "";
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
			
			// output this reference name
			ofs_fasta << strFASTA << endl;
			continue;
		}

		if (strFASTA.find("-") == string::npos) { //no deletions in this line
			sequence += strFASTA;
		} else {
			for (int i = 0; i < strFASTA.length(); i++) {
				if (strFASTA[i] != '-')
					sequence += strFASTA.substr(i,1);
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
	
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "FASTA File not specified" << endl;
		exit(0);
	}
//	string Mutation_Type (argv[1]);
	string FASTA_File (argv[1]);
	string OUT_File (argv[2]);
	string SNP_File (argv[3]);
	string Deletion_File (argv[4]);
	string Insertion_File (argv[5]);

	string Temp_File;
	if (SNP_File.length() > 0) {
		Temp_File = OUT_File + ".stmp";
		FixSNPs(FASTA_File, SNP_File, Temp_File);
		FASTA_File = Temp_File;
	}
	if (Deletion_File.length() > 0) {
		Temp_File = OUT_File + ".dtmp";
		FixDeletions(FASTA_File, Deletion_File, Temp_File);
		FASTA_File = Temp_File;
	}
	if (Insertion_File.length() > 0) {
		Temp_File = OUT_File + ".itmp";
		FixInsertions(FASTA_File, Insertion_File, Temp_File);
		FASTA_File = Temp_File;
	}
	FinaliseOutput(FASTA_File, OUT_File);
	remove((OUT_File + ".stmp").c_str());
	remove((OUT_File + ".dtmp").c_str());
	remove((OUT_File + ".itmp").c_str());
			
	return 0;
}
