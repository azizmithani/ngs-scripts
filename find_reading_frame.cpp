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

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

void ReadBlastOutput (string Blast_Output_File, map< string, vector< string > >& BlastOutput) {
	// open the file stream
	ifstream ifs_blast ( Blast_Output_File.c_str() );

	string strBLAST = "";
	while (ifs_blast.good()) {
		getline(ifs_blast, strBLAST);
		
		if (strBLAST.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strBLAST, entries, "\t");
		
		BlastOutput[entries[0]] = entries;
	}
	// close the stream
	ifs_blast.close();
}	

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Parameters: Sequence_File(string) Blast_Output_File(string) Out_File(string)" << endl;
		exit(0);
	}
	
	string Sequence_File (argv[1]);
	string Blast_Output_File (argv[2]);
	string OUT_File (argv[3]);
	
	// read the blast output in the map
	map< string, vector< string > > BlastOutput;
	ReadBlastOutput(Blast_Output_File, BlastOutput);
	
	ifstream ifs_seq ( Sequence_File.c_str() );
	string strSeq = "";
	string dna = "";
	string Id = "";
	while (ifs_seq.good()) {
		getline(ifs_seq, strSeq);
		
		if (strSeq.empty()) { // ignore emtry lines
			continue;
		}

		if (strSeq[0] == '>') {
			// new reference is starting ... process the previous one 

			if (dna.length() != 0) {			

				map< string, vector< string > >::iterator itBlastOutput = BlastOutput.find(Id);
				
				if (itBlastOutput == BlastOutput.end()) {
					// no frame information was found for this sequence by blast
					cout << Id << endl;
					
				} else {
					// get the blast output for this id
					vector<string> entries = itBlastOutput->second;
					// get the frame
					int frame = atoi(entries[13].c_str());
					// get coding sequence coordinates 
					int qStart = atoi(entries[7].c_str());
					int qEnd = atoi(entries[8].c_str());
										
					// backup of the original values
					int frameOriginal = frame;
					int	qStartOriginal = qStart;
					int qEndOriginal = qEnd;
					
					if (frame < 0) {
						dna = ReverseComplementDNA(dna);
						qStart = dna.length() - qStart + 1;
						qEnd = dna.length() - qEnd + 1;
						frame *= -1;
					}
					int qStartNew = qStart;
					int qEndNew = qEnd;
					
					// 0 Base Indexing
					qStart--;
					qEnd--;
					frame--;
					
					// translate the dna sequence
					string protein = TranslateDNASequence(dna, frame);
					
					// Get the first amino acid position in the translated sequence
					int pStart = qStart / 3;
					int pEnd = qEnd / 3;				
					
					string cds = protein.substr(pStart, pEnd - pStart);
					
					// check if this is an "M" (start codon)
					if (protein[pStart] != 'M') {
						// if not, check to the left of it for the first "M" such that there is no "*" (stop codon) in between.
						string cds_prefix = protein.substr(0, pStart + 1);
						int firstM = cds_prefix.find_last_of('M');
						int firstStop = cds_prefix.find_last_of('*');
						
						//cout << "\t" << firstM << "\t" << firstStop;
						if (firstM > firstStop) { // an "M" is present to the left with no stop codon in the middle
							// change the query start to this position
							qStartNew = firstM * 3 + frame;
						} else if (firstM < 0 && firstStop < 0) {// if neither start nor stop codon is present, then take the sequence from the begining 
							qStartNew = 1;
						}
						
					}
					
					if (frameOriginal < 0) {
						qStartNew = dna.length() - qStartNew + 1;
						qEndNew = dna.length() - qEndNew + 1;
					}
					cout << Id << "\t" << qStartOriginal << "\t" << qEndOriginal << "\t" << frameOriginal << "\t" << protein[pStart] << "\t" <<  qStartNew << "\t" << qEndNew << endl;
					
					
					// Get the amino acid after the CDS
					
					// check if this is  a stop codon "*"
					
					// if not, check to the right of it for the first "*".
					// if no stop codon is found then take the sequence up to the end
					
					
					/*
					 
					 // Get the position of first M in the translated sequence			
					 int firstM = protein.find_first_of('M');
					 
					 if (firstM ) {
					 }
					 
					 // extract the coding sequence
					 string cds = dna.substr(qStart, qEnd - qStart);
					 // translate the dna sequence
					 string protein = TranslateDNASequence(cds);
					 
					 if (protein[0] == 'M') { // protein starts with an 'M'
					 // do nothing
					 } else {
					 }
					 
					 dna = "";
					 
					 // get the reference name
					 vector<string> entries;
					 strsplit(seq, entries, " ");
					 
					 reference = entries[0].substr(1,entries[0].length());
					 TrimLeadingSpaces(reference);
					 full_reference = entry[0];		
					 
					 */
					
				}
			}
			
			Id = ExtractReferenceName(strSeq);
			dna = "";
				

		} else{
			RemoveNewLine(strSeq);
			dna += strSeq;
		}
	}
	// open the stream for output
	ofstream ofs (OUT_File.c_str());

	// close the stream
	ofs.close();
	
	return 0;
}
