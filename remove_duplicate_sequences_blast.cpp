/*
s *
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

int main (int argc, char** argv) {
	
	if (argc < 7) {
		cout << "Parameters: BLAST_DB(string) Sequence_File(sring) Out_File(string) Self_Blast(true/false) Percent_Identity(double) Align_Coverage_Percentage(double)" << endl;
		exit(0);
	}
	
	string BLAST_DB (argv[1]);
	string Sequence_File (argv[2]);
	string OUT_File (argv[3]);
	string Self_Blast(argv[4]);
	
	double Percent_Identity = atof(argv[5]);//95.00;
	double Align_Coverage = atof(argv[6]);//95.00;

	bool MaskLowComplexity = false;
	if (argc >= 7) {
		string ComplexityFiltering(argv[7]);
		if (ComplexityFiltering[0] == 'T' || ComplexityFiltering[0] == 't' ) {
			MaskLowComplexity = true;
		}
	}
	
	
	string TEMP_File = OUT_File + ".tmp";
	string command = "blastn -query " + Sequence_File + " -db " + BLAST_DB + " -out " + TEMP_File + " -outfmt '6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore'";
	if (MaskLowComplexity) {
		stringstream s;
		s << Percent_Identity;
		command += " -dust no -perc_identity " + s.str() + " -num_threads 4";
	}
	system(command.c_str());
	
	bool isSelfBlast = false;
	if (Self_Blast[0] == 'T' || Self_Blast[0] == 't' ) {
		isSelfBlast = true;
	}
	
	// open the sequence file stream
	ifstream ifs ( Sequence_File.c_str() );
	string strSeq = "";
	string sequence = "";
	string reference = "";
	map<string, int> SeqLengths;
	while (ifs.good()) {
		getline(ifs, strSeq);
		
		if (strSeq.empty()) { // ignore emtry lines
			continue;
		}
		
		if (strSeq[0] == '>') { // sequence header
			if (reference.length() != 0) {
				// save the value for previous sequence
				SeqLengths[reference] = sequence.length();
			}
			
			vector<string> entries;
			strsplit(strSeq, entries, " ");
			
			reference = entries[0].substr(1,entries[0].length());
			TrimLeadingSpaces(reference);
			
			sequence = "";
		} else {
			RemoveNewLine(strSeq);
			sequence += strSeq;
		}
	}
	// save the value for the last sequence
	SeqLengths[reference] = sequence.length();
	ifs.close();
	
	// open the file stream
	ifstream ifs_blast ( TEMP_File.c_str() );
	// read the blast output in a map
	map< string, vector< vector<string> > > SequenceList;
	string strBLAST = "";
	while (ifs_blast.good()) {
		getline(ifs_blast, strBLAST);
		
		if (strBLAST.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strBLAST, entries, "\t");
		
		SequenceList[entries[0]].push_back(entries);
	}
	// close the stream
	ifs_blast.close();
		
	set< string > AcceptedSequences;
	map< string, vector< vector<string> > >::iterator it;
	if (isSelfBlast) {
		// accept the sequence if only one hit is found
		for (it = SequenceList.begin(); it != SequenceList.end(); it++) {
			if (it->second.size() == 1) {
				AcceptedSequences.insert(it->first);
			} 
		}
	}
	//cout << AcceptedSequences.size() << endl;
	for (it = SequenceList.begin(); it != SequenceList.end(); it++) {
		bool isDuplicate = false;
		string SelectedSeq = "";
		if (isSelfBlast && it->second.size() == 1)
			continue;
		vector< vector<string> >::iterator itSeq;
		for (itSeq = it->second.begin(); itSeq != it->second.end(); itSeq++) {
			vector<string> entries = *itSeq;
			int QueryAlignLength = atoi(entries[8].c_str()) - atoi(entries[7].c_str()) + 1; // query end - query start + 1
			if (isSelfBlast && entries[0].compare(entries[1]) == 0 && atof(entries[2].c_str()) == 100.00 && QueryAlignLength == SeqLengths[entries[0]]) { // self blast, this is the same sequence
				//SelectedSeq = entries[0];
				continue;
			} else if (atof(entries[2].c_str()) < Percent_Identity || (double)QueryAlignLength/(double)SeqLengths[entries[0]]*100 < Align_Coverage ) { // ignore entries where there's no good match
				continue;
			//} else if (isSelfBlast && AcceptedSequences.find(entries[1]) != AcceptedSequences.end()) {
			} else if (AcceptedSequences.find(entries[1]) != AcceptedSequences.end()) {
				isDuplicate = true;
				break;
			} else {
				if (SelectedSeq.length() == 0)
					SelectedSeq = ( SeqLengths[entries[1]] > SeqLengths[entries[0]] ? entries[1] : entries[0] );
				else
					SelectedSeq = ( SeqLengths[entries[1]] > SeqLengths[SelectedSeq] ? entries[1] : SelectedSeq );
				
				if (!isSelfBlast) {
					isDuplicate = true;
					break;
				}
			}
		} // end for all sequences
		if (!isDuplicate) {
			if (SelectedSeq.length() > 0 && AcceptedSequences.find(SelectedSeq) == AcceptedSequences.end()) {
				AcceptedSequences.insert(SelectedSeq);
			} else {
				AcceptedSequences.insert(it->first);
			}
		}		
	} // end for all queries

	// open the fasta file stream
	ifstream ifs_fasta ( Sequence_File.c_str() );
	// open the stream for output
	ofstream ofs (OUT_File.c_str());
	strSeq = "";
	bool accepted = false;
	while (ifs_fasta.good()) {
		getline(ifs_fasta, strSeq);
		
		if (strSeq.empty()) { // ignore emtry lines
			continue;
		}
		
		if (strSeq[0] == '>') { // sequence header
			vector<string> entries;
			strsplit(strSeq, entries, " ");
			
			string reference = entries[0].substr(1,entries[0].length());
			TrimLeadingSpaces(reference);
			
			if (AcceptedSequences.find(reference) != AcceptedSequences.end()) {
				accepted = true;
				ofs << strSeq << endl;
			} else if (SequenceList.find(reference) == SequenceList.end()) { // no hit was found for this reference. accept it
				accepted = true;
				//cout << reference << endl;
				ofs << strSeq << endl;
			} else
				accepted = false;
		} else {
			if (accepted) {
				ofs << strSeq << endl;
			}
		}
	}
	
	command = "rm " + TEMP_File;
	system(command.c_str());
	
	return 0;
}
