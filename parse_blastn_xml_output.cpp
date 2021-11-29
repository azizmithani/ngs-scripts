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
#include <limits>
#include "Utilities.h"


using namespace std;

struct HSP {
	int HspFrom;
	int HspTo;
	int HspFrame;
	double HspEValue;
	int HspIdentity;
	int HspGaps;
	int HspAlignLength;
	int HspFromOld;
	int HspToOld;
};
struct HIT {
	string HitId;
	string HitDef;
	int HitLength;
	vector<HSP> Hsps;
};
struct QUERY_OUTPUT {
	string QueryDef;
	int QueryLength;
	vector<HIT> Hits;
};
typedef map< string, QUERY_OUTPUT > BLAST_OUTPUT;

int nNoHit = 0;
int nSingleHit = 0;
int nMultipleHit = 0;

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

void ReadBlastOutput (string Blast_Output_File, BLAST_OUTPUT& BlastOutput) {
	// open the file stream
	ifstream ifs_blast ( Blast_Output_File.c_str() );

	QUERY_OUTPUT *QueryOutput = NULL;
	HIT *Hit = NULL;
	HSP *Hsp = NULL;
	int pos = 0;
	//QueryOutput.QueryDef = "";
	string strBLAST = "";
	while (ifs_blast.good()) {
		getline(ifs_blast, strBLAST);
		
		if (strBLAST.empty()) { // ignore emtry lines
			continue;
		}
		
		//cout << strBLAST << endl;
		if ((pos = strBLAST.find("<Iteration_query-def>")) != string::npos) {
			// save previous query output
			if (QueryOutput != NULL) {
 				// save the last hit for previous query
				if (Hit != NULL) {
					if (Hsp != NULL) {
						Hit->Hsps.push_back(*Hsp);
						Hsp = NULL;
					}
					QueryOutput->Hits.push_back(*Hit);
					Hit = NULL;
				}
				
				string QueryId = QueryOutput->QueryDef.substr(0, QueryOutput->QueryDef.find_first_of(' '));
				BlastOutput[QueryId] = *QueryOutput;
				
			}
			// allocate a new query output block
			QueryOutput = new QUERY_OUTPUT;
			// save query description
			QueryOutput->QueryDef = strBLAST.substr(pos + 21, strBLAST.find_first_of('<', pos + 1) - pos - 21);
		} else if ((pos = strBLAST.find("<Iteration_query-len>")) != string::npos) {
			// save query length
			QueryOutput->QueryLength = atoi(strBLAST.substr(pos + 21, strBLAST.find_first_of('<', pos + 1) - pos - 21).c_str());
		} else if ((pos = strBLAST.find("<Hit_id>")) != string::npos) {
			// save previous hit for this query
			if (Hit != NULL) {
				// save the last Hsp for previous Hit
				if (Hsp != NULL) {
					Hit->Hsps.push_back(*Hsp);
					Hsp = NULL;
				}
				QueryOutput->Hits.push_back(*Hit);
			}
			// allocate a new hit block
			Hit = new HIT;
			// save Hit Id
			Hit->HitId = strBLAST.substr(pos + 8, strBLAST.find_first_of('<', pos + 1) - pos - 8);
		} else if ((pos = strBLAST.find("<Hit_def>")) != string::npos) {
			// save Hit description
			Hit->HitDef = strBLAST.substr(pos + 9, strBLAST.find_first_of('<', pos + 1) - pos - 9);
		} else if ((pos = strBLAST.find("<Hit_len>")) != string::npos) {
			// save Hit length
			Hit->HitLength = atoi(strBLAST.substr(pos + 9, strBLAST.find_first_of('<', pos + 1) - pos - 9).c_str());
		} else if ((pos = strBLAST.find("<Hsp_num>")) != string::npos) {
			// save previous hsp for this hit
			if (Hsp != NULL) {
				Hit->Hsps.push_back(*Hsp);
			}
			// allocate a new HSP block
			Hsp = new HSP;
		} else if ((pos = strBLAST.find("<Hsp_query-from>")) != string::npos) {
			// save Hsp query from
			Hsp->HspFrom = atoi(strBLAST.substr(pos + 16, strBLAST.find_first_of('<', pos + 1) - pos - 16).c_str());
		} else if ((pos = strBLAST.find("<Hsp_query-to>")) != string::npos) {
			// save Hsp query to
			Hsp->HspTo = atoi(strBLAST.substr(pos + 14, strBLAST.find_first_of('<', pos + 1) - pos - 14).c_str());
		} else if ((pos = strBLAST.find("<Hsp_query-frame>")) != string::npos) {
			// save Hsp query to
			Hsp->HspFrame = atoi(strBLAST.substr(pos + 17, strBLAST.find_first_of('<', pos + 1) - pos - 17).c_str());
		} else if ((pos = strBLAST.find("<Hsp_evalue>")) != string::npos) {
			// save Hsp evalue
			Hsp->HspEValue = atof(strBLAST.substr(pos + 12, strBLAST.find_first_of('<', pos + 1) - pos - 12).c_str());
		} else if ((pos = strBLAST.find("<Hsp_identity>")) != string::npos) {
			// save Hsp identity
			string strIdentity ();
			Hsp->HspIdentity = atoi(strBLAST.substr(pos + 14, strBLAST.find_first_of('<', pos + 1) - pos - 14).c_str());
		} else if ((pos = strBLAST.find("<Hsp_gaps>")) != string::npos) {
			// save Hsp gaps
			string strGaps (); 
			Hsp->HspGaps = atoi(strBLAST.substr(pos + 10, strBLAST.find_first_of('<', pos + 1) - pos - 10).c_str());
		} else if ((pos = strBLAST.find("<Hsp_align-len>")) != string::npos) {
			// save Hsp gaps
			string strGaps (); 
			Hsp->HspAlignLength = atoi(strBLAST.substr(pos + 15, strBLAST.find_first_of('<', pos + 1) - pos - 15).c_str());
		}
		
	} // end while ifs_blast is good
	// save last query output
	if (QueryOutput != NULL) {
		// save the last hit for previous query
		if (Hit != NULL) {
			// save the last Hsp for previous Hit
			if (Hsp != NULL) {
				Hit->Hsps.push_back(*Hsp);
				Hsp = NULL;
			}
			QueryOutput->Hits.push_back(*Hit);
			Hit = NULL;
		}
		
		string QueryId = QueryOutput->QueryDef.substr(0, QueryOutput->QueryDef.find_first_of(' '));
		BlastOutput[QueryId] = *QueryOutput;
	}
	
	// close the stream
	ifs_blast.close();
		
}	

void CheckHSP (string& dna, HSP& Hsp) {

	// get the frame
	int frame = Hsp.HspFrame;
	// get coding sequence coordinates 
	int qStart = Hsp.HspFrom;
	int qEnd = Hsp.HspTo;
	
	//cout << frame << "\t" << qStart << "\t" << qEnd << "\t" << dna.length() << endl;
	
	// backup of the original values
	int frameOriginal = frame;
	
	if (frame < 0) {
		dna = ReverseComplementDNA(dna);
		int	qStartOriginal = qStart;
		qStart = dna.length() - qEnd + 1;
		qEnd = dna.length() - qStartOriginal + 1;
		frame *= -1;
	}
	int qStartNew = qStart;
	int qEndNew = qEnd;
	
	// 0 Base Indexing
	qStart--;
	qEnd--;
	frame--;

	//cout << frame << "\t" << qStart << "\t" << qEnd << endl;
	

	// translate the dna sequence
	string protein = TranslateDNASequence(dna, frame);
	
	// Get the first & last amino acid positions in the translated sequence
	int pStart = qStart / 3;
	int pEnd = qEnd / 3;				
	
	string cds = protein.substr(pStart, pEnd - pStart);
	
	//cout << protein << "\t" << cds << endl;
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
	
	// check if the codon after the last a.a. is a stop codon "*"
	if (protein.length() >= pEnd && protein[pEnd+1] != '*') {
		// if not, check to the right of it for the first "*"
		int firstStop = protein.find_first_of('*', pEnd);
		if (firstStop == string::npos) { // no stop found, take the sequence up to the end
			qEndNew = dna.length();
		} else {
			// otherwise change the query end to this position
			qEndNew = firstStop * 3 + frame;
		}

	}

	if (frameOriginal < 0) {
		int qStartNewOriginal = qStartNew;
		qStartNew = dna.length() - qEndNew + 1;
		qEndNew = dna.length() - qStartNewOriginal + 1;
	}
	
	// save the old coordinates
	Hsp.HspFromOld = Hsp.HspFrom;
	Hsp.HspToOld = Hsp.HspTo;
	// update the coordinates
	Hsp.HspFrom = qStartNew;
	Hsp.HspTo = qEndNew;

}

bool CheckHSPFrames (vector<HSP>& Hsps) {
	vector<HSP>::iterator itHsp = Hsps.begin();
	// get the first frame
	int frame = itHsp->HspFrame;
	itHsp++;
	// check if the remaining frames are the same
	for (; itHsp != Hsps.end(); itHsp++) {
		if (itHsp->HspFrame != frame) {
			return false;
		}
	}
	// no discrepancy found, return true
	return true;
}

bool MergeHsps(vector<HSP>& Hsps, HSP& newHsp, int frame) {

	bool HspFound = false;
	newHsp.HspFrom = numeric_limits<int>::max();
	newHsp.HspTo = numeric_limits<int>::min();
	newHsp.HspFrame = frame;
	vector<HSP>::iterator itHsp = Hsps.begin();
	vector< vector<HSP>::iterator > HspsToBeDeleted;
	// merge the Hsps
	for (; itHsp != Hsps.end(); itHsp++) {
		if (itHsp->HspFrame == frame) { // correct frame, merge it
			newHsp.HspFrom = min(itHsp->HspFrom, newHsp.HspFrom);
			newHsp.HspTo = max(itHsp->HspTo, newHsp.HspTo);

			// mark for deletion
			HspsToBeDeleted.push_back(itHsp); 
			
			// at least one Hsp found, mark the flag
			HspFound = true;
		}
	}
	// remove the Hsps which were marked for deletion
	vector< vector<HSP>::iterator >::reverse_iterator ritHspsToBeDeleted = HspsToBeDeleted.rbegin();
	for (; ritHspsToBeDeleted != HspsToBeDeleted.rend(); ++ritHspsToBeDeleted) {
		Hsps.erase(*ritHspsToBeDeleted);
	}
	
	return HspFound;
}

void WriteHSP(ofstream& ofs, string Id, QUERY_OUTPUT& QueryOutput, HIT& Hit, HSP& Hsp) {
	ofs << Id << "\t" << Hit.HitId << "\t" << Hit.HitDef << "\t" << (double)(Hsp.HspTo - Hsp.HspFrom)/(double)QueryOutput.QueryLength * 100.00 << "\t" << Hsp.HspIdentity << "\t" << Hsp.HspGaps << "\t" << Hsp.HspAlignLength << "\t" << Hsp.HspEValue << endl;
}

void ProcessBlastOutput(BLAST_OUTPUT& BlastOutput, string& dna, string Id, ofstream& ofs) {
	BLAST_OUTPUT::iterator itBlastOutput = BlastOutput.find(Id);
	
	//cout << Id << "\t" << itBlastOutput->second.Hits.size() << endl;
	
	
	//cout << Id << endl;
	if (itBlastOutput != BlastOutput.end()) {
		// get the blast output for this id
		QUERY_OUTPUT QueryOutput = itBlastOutput->second;
		
		if (QueryOutput.Hits.size() == 0) { // no hit found for this query
			string protein = "";
			int qStart, qEnd;
			
//			int frame = FindReadingFrame(dna, protein, qStart, qEnd, true);
//			ofs << Id << "\t-\t-\t" << qStart << "\t" << qEnd << "\t" << frame << endl;
			
			// increment the count
			nNoHit++;
		} else if (QueryOutput.Hits.size() == 1) { // only one hit found for this query
			HIT *Hit  = &(*QueryOutput.Hits.begin()); 
			string HitId = Hit->HitId;
			string HitDef = Hit->HitDef;
			int HitLength = Hit->HitLength;
			
			//cout << QueryOutput.QueryDef << "\t" << QueryOutput.Hits.size() << "\t" << Hit->Hsps.size() << endl;
			if (Hit->Hsps.size() == 1) { // only one hsp found for this hit (majority scenario)
				//CheckHSP(dna, *Hit->Hsps.begin());
				WriteHSP(ofs, Id, QueryOutput, *Hit, *Hit->Hsps.begin());
				
				// increment the count
				nSingleHit++;
			} else {
				vector<HSP>::iterator the_hsp = Hit->Hsps.begin();
				for (; the_hsp != Hit->Hsps.end(); the_hsp++) {
					WriteHSP(ofs, Id, QueryOutput, *Hit, *the_hsp);					
				}
				
				nMultipleHit++;

			}
		}
	}
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Parameters: Sequence_File(string) Blast_Output_File(string) Out_File(string)" << endl;
		exit(0);
	}
	
	//blastx -query /san/001/Aziz/Wheat_UniGenes_Build_58/Ta.seq.uniq -out /san/001/Aziz/Wheat_UniGenes_Build_58/Ta.seq.uniq.blastx.xml -outfmt 5 -db nr -num_threads 8 -max_target_seqs 1
	
	string Sequence_File (argv[1]);
	string Blast_Output_File (argv[2]);
	string OUT_File (argv[3]);
	
	// read the blast output in the map
	BLAST_OUTPUT BlastOutput;
	ReadBlastOutput(Blast_Output_File, BlastOutput);
	
	// open the stream for output
	ofstream ofs (OUT_File.c_str());
	
	// read the sequence file
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
				//cout << Id << endl;
				ProcessBlastOutput(BlastOutput, dna, Id, ofs);
			}
			
			// Extract the Id for the new sequence
			Id = ExtractReferenceName(strSeq);
			// reset dna
			dna = "";
				

		} else{
			RemoveNewLine(strSeq);
			dna += strSeq;
		}
	}
	// process the last sequence
	ProcessBlastOutput(BlastOutput, dna, Id, ofs);

	cout << "No Hits: " << nNoHit << endl;
	cout << "Single Hits: " << nSingleHit << endl;
	cout << "Multiple Hits: " << nMultipleHit << endl;
	// close the stream
	ofs.close();
	
	return 0;
}
