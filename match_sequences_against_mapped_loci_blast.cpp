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
	
//	if (argc < 7) {
//		cout << "Parameters: BLAST_DB(string) Sequence_File(sring) Out_File(string) Self_Blast(true/false) Percent_Identity(double) Align_Coverage_Percentage(double)" << endl;
//		exit(0);
//	}
	
	string BLAST_DB (argv[1]);
	//string Sequence_File (argv[2]);
	string Mapped_Loci_File (argv[2]);
	string OUT_File (argv[3]);

	double Percent_Identity = atof(argv[4]);//80.00;
	double Percent_Coverage = atof(argv[5]);//80.00;
		
	map< string, int > ChromosomeMap; // map < seq_name, chr_number (1..7) >
	string TEMP_File = OUT_File + ".tmp";
	string TEMP_OUT_File = OUT_File + ".out.tmp";
	
	// read the mapped loci file
	int skip = 1; // skip the header
	vector< vector< string > > MappedLoci = ReadFile(Mapped_Loci_File, skip);
	// go through each mapped locus
	vector< vector< string > >::iterator itMappedLoci = MappedLoci.begin();
	for (; itMappedLoci != MappedLoci.end(); itMappedLoci++) {
		// get the current entry
		vector< string > MappedLocus = *itMappedLoci;
		
		//cout << MappedLocus[0] << "\t" << MappedLocus[1] << "\t" << MappedLocus[2] << "\t" << MappedLocus[3] << endl;
		
		// write the sequence to a temporary file
		ofstream ofs_temp ( TEMP_File.c_str() ); // open the output file stream		
		ofs_temp << ">" << MappedLocus[0] << endl; // write the sequence header
		ofs_temp << MappedLocus[4] << endl; // the sequence
		ofs_temp.close(); // close the output stream
		
		// get the sequence length
		int SeqLength = MappedLocus[4].length();
		
		// blast the sequence against the given database
		string command = "blastn -query " + TEMP_File + " -db " + BLAST_DB + " -out " + TEMP_OUT_File + " -outfmt '6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore'";
		system(command.c_str());
		
		// open the file stream
		ifstream ifs_blast ( TEMP_OUT_File.c_str() );
		string strBLAST = "";
		while (ifs_blast.good()) {
			getline(ifs_blast, strBLAST);
			
			if (strBLAST.empty()) { // ignore emtry lines
				continue;
			}

			//cout << strBLAST << endl;
		
			vector<string> entries;
			strsplit(strBLAST, entries, "\t");
			
			int QueryAlignLength = atoi(entries[8].c_str()) - atoi(entries[7].c_str()) + 1; // query end - query start + 1
			if (atof(entries[2].c_str()) < Percent_Identity || (double)QueryAlignLength/(double)SeqLength*100 < Percent_Coverage ) { // ignore entries where there's no good match
				// do nothing
			} else {
				
				string SeqName = entries[1];
				int Chr = atoi(entries[0].c_str());
				
				// see if this sequence has been assigned a chromosome yet
				map< string, int >::iterator itChromosomeMap = ChromosomeMap.find(SeqName);
				if (itChromosomeMap == ChromosomeMap.end()) { //this sequence hasn't been assigned a chromosome yet
					ChromosomeMap[SeqName] =  Chr;
					
					cout << SeqName << "\t" << Chr << "\t" << Percent_Identity << "\t" << Percent_Coverage << "\t" <<  MappedLocus[3] << endl;
				} else if (itChromosomeMap->second != Chr) { // this sequence was not assigned the same choromosome
					// remove the mapping (we don't want any sequences which have duplicate mapping)
					//cout << "Removing " << SeqName << endl;
					ChromosomeMap.erase(itChromosomeMap);
				}

			}
		} // end while blast file is good
		// close the stream
		ifs_blast.close();
		
	} // end for itMappedLoci

	// remove the temporary files
	string command = "rm -f " + TEMP_File + " " + TEMP_OUT_File;
	system(command.c_str());

	ofstream ofs_out ( OUT_File.c_str() ); // open the output file stream
	map< string, int >::iterator itChromosomeMap = ChromosomeMap.begin();
	for (; itChromosomeMap != ChromosomeMap.end(); itChromosomeMap++) {
		ofs_out << itChromosomeMap->first << "\t" << itChromosomeMap->second << endl;
	}
	ofs_out.close(); // close the output stream
	
	
/*	
		
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

 */
	
	return 0;
}
