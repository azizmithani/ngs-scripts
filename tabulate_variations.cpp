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
#include <map>
#include <set>
#include <limits>
#include "Utilities.h"

using namespace std;


const int TABULATE_VARIATIONS = 0;
const int INDIVIDUAL_VARIATIONS = 1;
const int GEOGRAPHIC_VARIATIONS = 2;
const int MAJORITY_VARIATIONS = 3;
const int VARIATIONS_PER_GENE = 4;
const int UNIQUE_VARIATIONS_PER_GENE = 5;
const int SYNONYMOUS_VS_NON_SYNONYMOUS_SNP = 6;
const int SYNONYMOUS_VS_NON_SYNONYMOUS_UNIQUE_SNP = 7;

const int MIN_COVERAGE = 3;
const double MIN_CONSENSUS_SNP = 80.0;
const double MIN_CONSENSUS_INDEL = 30.0;
const int FLANKING_BASES = 5; // number of bases on either sides of a variation to be checked 

struct HSP {
	int HspFrom;
	int HspTo;
	int HspFrame;
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


//Global variables
bool isSNPData;
bool isInsData;
int offset;

string ExtractReferenceName(string str) {
	if (str.find_first_of(' ') != string::npos) {
		return str.substr(1,str.find_first_of(' ')-1);
	} else {
		return str.substr(1, str.length());
	}
	
}

void ReadBlastOutput (string Blast_Output_File, BLAST_OUTPUT& BlastOutput) {
	
	cout << "Reading reading frame and coding sequence coordinates ... " << endl;
	
	// open the file stream
	ifstream ifs_blast ( Blast_Output_File.c_str() );
	
	QUERY_OUTPUT *QueryOutput = NULL;
	HIT *Hit = NULL;
	HSP *Hsp = NULL;
	int pos = 0;
	//QueryOutput.QueryDef = "";
	string strBLAST = "";
	string QueryId = "";
	string PreviousQueryId = "";
	while (ifs_blast.good()) {
		getline(ifs_blast, strBLAST);
		
		if (strBLAST.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strBLAST, entries, "\t");
		

		// query description
		QueryId = entries[0];

		BLAST_OUTPUT::iterator itBlastOutput = BlastOutput.find(QueryId);
		if (itBlastOutput == BlastOutput.end()) {
			// allocate a new query output block
			QueryOutput = new QUERY_OUTPUT;
		} else {
			*QueryOutput = itBlastOutput->second;
		}


		// Hit Id
		string HitId = entries[1];
		
		// check if this hit has been seen before (multiple HSPs)
		Hit = NULL;
		vector<HIT>::iterator itHit = QueryOutput->Hits.begin();
		for (; itHit != QueryOutput->Hits.end(); itHit++) {
			if (HitId.compare(itHit->HitId) == 0) { // This hit has been seen before
				Hit = &(*itHit);
				break;
			}
		}
		if (Hit == NULL) { // New Hit
			// allocate a new Hit block
			Hit = new HIT;
			// Hit Id
			Hit->HitId = HitId;
			// Hit description
			Hit->HitDef = entries[2];
		}
		
		
		HSP Hsp;
		Hsp.HspFrom = atoi(entries[3].c_str());
		Hsp.HspTo = atoi(entries[4].c_str());
		Hsp.HspFrame = atoi(entries[5].c_str());
		
		Hit->Hsps.push_back(Hsp);
		
		if (itHit == QueryOutput->Hits.end()) { // this hit was not seen before, add it to the list of Hits
			QueryOutput->Hits.push_back(*Hit);
		} else { // this hit was seen before, replace it.
			*itHit = *Hit;
		}

		BlastOutput[QueryId] = *QueryOutput;
	}
	
	// close the stream
	ifs_blast.close();
	
}	



void TabulateVariation (map <string, string>& Accessions, map < string, vector <string> >& VariationTable) {
	
	cout << "Tabulating varations ..." << endl;
	
	map <string, string>::iterator itAccession = Accessions.begin();
	int idxAcc = -1;
	for (; itAccession != Accessions.end(); itAccession++) {
		string Accession = itAccession->first;
		string Filename =  itAccession->second;
		idxAcc++;
		
		cout << "\t" << Accession << endl;
		
		ifstream ifs( Filename.c_str() );
		if (!ifs) { // The file does not exist, move to the next one
			cout << "\tFile '" << Filename << "' not found" << endl;
			continue;
		}
		
		string str = "";
		while (ifs.good()) {
			getline(ifs, str);
			
			if (str.empty()) { // ignore emtry lines
				continue;
			}
			
			vector<string> entries;
			strsplit(str, entries, "\t");

			// create the key (chromosome:position)
			string key = entries[0] + ":" + entries[1];
			
			// find the key in the table
			map < string, vector <string> >::iterator itVariation = VariationTable.find(key);

			vector <string> VariationRow(Accessions.size() + offset);
			if (itVariation != VariationTable.end()) {
				VariationRow = itVariation->second;
			} else {
				//VariationRow(Accessions.size() + 3);
				VariationRow[0] = entries[0];
				VariationRow[1] = entries[1];
				if (isSNPData) {
					VariationRow[2] = entries[2];
				}
			}
			
			if (isSNPData) { // save the consensus base for SNP
				VariationRow[idxAcc + offset] = entries[3];
			} else if (isInsData) { // save the base(s) inserted if this is the ins data
				VariationRow[idxAcc + offset] = entries[2];
			} else { // length (del)
				string length;
				stringstream ss;
				ss << entries[2].length();
				ss >> length;
				VariationRow[idxAcc + offset] = length;
			}
			
			if (itVariation != VariationTable.end()) {
				itVariation->second = VariationRow;
			} else {
				VariationTable[key] = VariationRow;
			}
			
		} // end while ifs.good()

		// close the file
		ifs.close();

	} // end for each accession

	
}

void ReadBaseDistributionFile(string Base_Dist_File, map < string, int** >& Frequency, map < string, int >& RefLength) {
	
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
				*(RefFrequency+i) = new int[TOTAL_BASES];		
				for (int j = 0; j < TOTAL_BASES; j++) {
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
}

void ReadCoverageFile (string Coverage_File, map< string, int >& Coverage) {
	// open the Coverage file stream
	ifstream ifs_Coverage ( Coverage_File.c_str() );
	string str="";
	while (ifs_Coverage.good()) {
		getline(ifs_Coverage, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		Coverage[ entries[0] + ":" + entries[1] ] = atoi(entries[2].c_str());
		
	}
	ifs_Coverage.close();
}

string VerifyData(string data, bool isSNPData, vector<string>& VariationRow) {
	if (data.length() == 0) {
		if (isSNPData) {
			return VariationRow[2];
		} else {
			return "-";
		}
	}
	return data;
}

int CalculateCoverage(int* BaseDist) {
	int coverage = 0;
	for (int i = 0; i < TOTAL_BASES; i++) {
		coverage += BaseDist[i];
		//cout << "\t" << BaseDist[i];
	}
	
	//cout << "\t" << coverage << endl;

	return coverage;
}

void CheckBaseDistribution(map <string, string>& BaseDistributionsFiles, map <string, string>& Accessions, map < string, vector <string> >& VariationTable){

//	map< string, map < string, int** > > BaseDistributions;
	map < string, int > RefLength;

	cout << "Checking base distributions ..." << endl;
	
	map <string, string>::iterator itAccession = Accessions.begin();
	int idxAcc = -1;
	for (; itAccession != Accessions.end(); itAccession++) {
		set<string> VariationsToBeDeleted;

		idxAcc++;
		string Accession = itAccession->first;
		cout << "\t" << Accession << endl;
				
		map <string, string>::iterator itBaseDistributionsFiles = BaseDistributionsFiles.find(Accession);		
		string Filename =  itBaseDistributionsFiles->second;
		
		map < string, int** > Frequency;
		ReadBaseDistributionFile (Filename, Frequency, RefLength);
		
		string Reference = "";
		string PreviousReference = "";
		int** RefFrequency;
		int ReferenceLength;
		map <string, vector<string> >::iterator itVariation = VariationTable.begin();
		for (; itVariation != VariationTable.end(); itVariation++) {
			vector<string> VariationRow = itVariation->second;
			
			Reference = VariationRow[0];
			long Position = atol(VariationRow[1].c_str());
			
			if (Reference.compare(PreviousReference) != 0) {
				RefFrequency = Frequency[Reference];
				ReferenceLength = RefLength[Reference];
				
				PreviousReference = Reference;
			}

			string data =  VerifyData(VariationRow[offset + idxAcc], isSNPData, VariationRow);
			string RefBase = VariationRow[2];
			
			if (data[0] != RefBase[0]) {
				// if there is a SNP w.r.t the reference, dont check this position
				continue;
			}

			// check the base coverage and the coverages on either sides of the variation, it must be >= the threshold
			bool LowCoverage = false;
			int BaseCoverage = 0;
			for (int p = Position - FLANKING_BASES; p <= Position + FLANKING_BASES; p++) {
				if (p <= 0 || p > ReferenceLength) {
					continue; // safety check
				}
				int Coverage = CalculateCoverage(RefFrequency[p - 1]);
				if (p == Position) {
					BaseCoverage = Coverage;
				}
				if (Coverage < MIN_COVERAGE) { // delete this variation if the coverage is below the threshold for this accession
					LowCoverage = true;
					break;
				}
			}
			if (LowCoverage) {
				VariationsToBeDeleted.insert(itVariation->first);
				continue;
			}
			
			// This might not have been called a SNP because consensus percentage for SNP didn't meet the threshold. 
			// Check if the reference base meets the consensus. if not, delete this variation
			int* BaseDist = RefFrequency[Position - 1];
			int dist = BaseDist[GetBaseIndexBaseDistFile(data[0])];
			if ((double)dist / (double)BaseCoverage * 100.00 < MIN_CONSENSUS_SNP) {
				//cout << itVariation->first << " Low Consensus" << endl;
				VariationsToBeDeleted.insert(itVariation->first);
				continue;
			}
			
			
		}
		
		
		set<string>::iterator itVariationsToBeDeleted = VariationsToBeDeleted.begin();
		for (; itVariationsToBeDeleted != VariationsToBeDeleted.end(); itVariationsToBeDeleted++) {
			VariationTable.erase(VariationTable.find(*itVariationsToBeDeleted));		
		}
		
	} // for each accession
}

void CheckCoverage(map <string, string>& BaseDistributionsFiles, map <string, string>& CoverageFiles, map <string, string>& Accessions, map < string, vector <string> >& VariationTable){
	
//	map< string, map < string, int > > Coverages;
	map < string, int > RefLength;
	
	cout << "Checking coverages ..." << endl;
	
	map <string, string>::iterator itAccession = Accessions.begin();
	int idxAcc = -1;
	for (; itAccession != Accessions.end(); itAccession++) {
		set<string> VariationsToBeDeleted;
		
		idxAcc++;
		string Accession = itAccession->first;
		cout << "\t" << Accession << endl;
		
		// Get the Coverage Filename for this accesssion
		map <string, string>::iterator itCoverageFiles = CoverageFiles.find(Accession);		
		string Filename =  itCoverageFiles->second;
		// Read the Coverage file
		map < string, int > Coverage;
		ReadCoverageFile (Filename, Coverage);
		
		// Get the Base Distribution Filename for this accesssion
		map <string, string>::iterator itBaseDistributionsFiles = BaseDistributionsFiles.find(Accession);		
		Filename =  itBaseDistributionsFiles->second;
		// Read the Base Distribution File
		map < string, int** > Frequency;
		ReadBaseDistributionFile (Filename, Frequency, RefLength);
		
		string Reference = "";
		string PreviousReference = "";
		int** RefFrequency;
		int ReferenceLength;
		map <string, vector<string> >::iterator itVariation = VariationTable.begin();
		for (; itVariation != VariationTable.end(); itVariation++) {
			vector<string> VariationRow = itVariation->second;
			
			Reference = VariationRow[0];
			int Position = atoi(VariationRow[1].c_str());
			
			if (Reference.compare(PreviousReference) != 0) {
				RefFrequency = Frequency[Reference];
				ReferenceLength = RefLength[Reference];
				
				PreviousReference = Reference;
			}
			
			string data =  VerifyData(VariationRow[offset + idxAcc], isSNPData, VariationRow);
			
			// check the base coverage and the coverages on either sides of the variation, it must be >= the threshold
			bool LowCoverage = false;
			int PositionCoverage = 0;
			for (int p = Position - FLANKING_BASES; p <= Position + FLANKING_BASES; p++) {
				if (p < 0 || p >= ReferenceLength) {
				//if (p < 0) {
					continue; // safety check
				}

				string pos;
				stringstream ss;
				ss << p;
				ss >> pos;
				
				string key = Reference + ":" + pos;
				map<string, int>::iterator itCoverage = Coverage.find(key);
				if (itCoverage == Coverage.end()) { // zero coverage
					LowCoverage = true; 
					break;
				}
				int Coverage = itCoverage->second;
				if (p == Position) {
					PositionCoverage = Coverage;
				}
				if (Coverage < MIN_COVERAGE) { // delete this variation if the  coverage is below the threshold for this accession
					LowCoverage = true; 
					break;
				}
			}
			if (LowCoverage) {
				VariationsToBeDeleted.insert(itVariation->first);
				continue;
			}

			if (data.compare("-") == 0) { // check the consensus if this was not called a deletion
				// This might not have been called a INS/DEL because consensus percentage for the INS/DEL didn't meet the threshold. 
				int BaseCoverage = CalculateCoverage(RefFrequency[Position - 1]); // 0-based indexing
				if (100.00 - ((double)BaseCoverage / (double)PositionCoverage * 100.00) > MIN_CONSENSUS_INDEL) { // check if at least X% of reads have an INS/DEL, if so this is an ambiguous position, delete it
					//cout << itVariation->first << " Low Consensus" << endl;
					VariationsToBeDeleted.insert(itVariation->first);
					continue;
				}
			}				
		}		
		
		set<string>::iterator itVariationsToBeDeleted = VariationsToBeDeleted.begin();
		for (; itVariationsToBeDeleted != VariationsToBeDeleted.end(); itVariationsToBeDeleted++) {
			VariationTable.erase(VariationTable.find(*itVariationsToBeDeleted));		
		}
		
	} // for each accession
}

void WriteTable (string OUT_File, map <string, string>& Accessions, map < string, vector <string> >& VariationTable) {
	
	cout << "Writing output ..." << endl;

	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		

	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	//cout << itVariation->second.size() << "\t" << offset << "\t" << Accessions.size() << endl;
	if (itVariation->second.size() - offset >= Accessions.size()) { // print the header only if full table is to be written
		ofs << "Chromosome\tPosition";
		if (isSNPData) {
			ofs << "\tRef Base";
		}
		map <string, string>::iterator itAccession = Accessions.begin();
		for (; itAccession != Accessions.end(); itAccession++) {
			string Accession = itAccession->first;
			ofs << "\t" << Accession;
		}
		ofs << endl;
	}
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		//ofs << itVariation->first;
		vector<string>::iterator itVariationRow = VariationRow.begin();
		ofs << VerifyData(*itVariationRow, isSNPData, VariationRow);
		itVariationRow++;
		for (; itVariationRow != VariationRow.end(); itVariationRow++) {
			//cout << "\t" << *itVariationRow;
			ofs << "\t" << VerifyData(*itVariationRow, isSNPData, VariationRow);
		}	// end for each entry of the variation row
		//cout << endl;
		ofs << endl;
	}
	
	ofs.close();
}

void RemoveCommonVariation (map < string, vector <string> >& VariationTable, map<string, string>& Accessions) {
	
	cout << "Removing common variations ..." << endl;
		
	vector<string> VariationsToBeDeleted;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		int idxCol = -1;
		int BaseCounts[TOTAL_BASES] = {0,0,0,0,0};
		map<string, int> IndelCounts;
		vector<string>::iterator itVariationRow = VariationRow.begin();
		vector<string>::iterator itFirstAccession, itLastAccession;
		for (; itVariationRow != VariationRow.end(); itVariationRow++) {
			idxCol++;
			
			//cout << PreviousValue << "\t" << *itVariationRow << endl;
			if (idxCol >= offset) {
				if (idxCol == offset) {
					itFirstAccession = itVariationRow;
				}
				itLastAccession = itVariationRow;
				
				string data = VerifyData(*itVariationRow, isSNPData, VariationRow);
				if (isSNPData) {
					BaseCounts[GetBaseIndex(data[0])]++;
				} else {
					IndelCounts[data]++;
				}
			}
		}	// end for each entry of the variation row
		
		//cout << VariationRow[0] << "\t" << VariationRow[1] << "\t" <<VariationRow[2] << "\t" << BaseCounts[0] << "\t" << BaseCounts[1] << "\t" << BaseCounts[2] << "\t" << BaseCounts[3] << endl;
		bool isCommonVariation = false;
		string Data = "";
		if (isSNPData) {
			for (int idx = 0; idx < TOTAL_BASES; idx++) {
				if (BaseCounts[idx] == Accessions.size()) {
					isCommonVariation = true;
					break;
				}
			}
		} else {
			map<string, int>::iterator itIndelCounts = IndelCounts.begin();
			for (; itIndelCounts != IndelCounts.end(); itIndelCounts++) {
				if (itIndelCounts->second == Accessions.size()) {
					isCommonVariation = true;
					break;
				}
			}
		}
		
		if (isCommonVariation) { // common variation -> delete it
			VariationsToBeDeleted.push_back(itVariation->first);
		}
		
	} // for each row of the variation table
	
	vector<string>::iterator itVariationsToBeDeleted = VariationsToBeDeleted.begin();
	for (; itVariationsToBeDeleted != VariationsToBeDeleted.end(); itVariationsToBeDeleted++) {
		VariationTable.erase(VariationTable.find(*itVariationsToBeDeleted));		
	}
}

void FilterForIndividualVariation (map < string, vector <string> >& VariationTable, map<string, string>& Accessions) {

	cout << "Filtering for individual variations ..." << endl;

	vector<string> VariationsToBeDeleted;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		if (isSNPData) {
			string data = VariationRow[2];
			if (data[0] == 'N') {
				VariationsToBeDeleted.push_back(itVariation->first);
				continue;
			}
		}
		
		int idxCol = -1;
		int BaseCounts[TOTAL_BASES] = {0,0,0,0,0};
		map<string, int> IndelCounts;
		vector<string>::iterator itVariationRow = VariationRow.begin();
		for (; itVariationRow != VariationRow.end(); itVariationRow++) {
			idxCol++;
			
			//cout << PreviousValue << "\t" << *itVariationRow << endl;
			if (idxCol >= offset) {
				string data = VerifyData(*itVariationRow, isSNPData, VariationRow);
				if (isSNPData) {
					BaseCounts[GetBaseIndex(data[0])]++;
				} else {
					IndelCounts[data]++;
				}
			}
		}	// end for each entry of the variation row
		
		//cout << VariationRow[0] << "\t" << VariationRow[1] << "\t" <<VariationRow[2] << "\t" << BaseCounts[0] << "\t" << BaseCounts[1] << "\t" << BaseCounts[2] << "\t" << BaseCounts[3] << endl;
		bool isIndvidualVariation = false;
		set<string> SelectedAccessions;
		if (isSNPData) {
			for (int idx = 0; idx < TOTAL_BASES; idx++) {
				if (BaseCounts[idx] == 1) {
					isIndvidualVariation = true;
					
					string SelectedData;
					stringstream ss;
					ss << GetIndexBase(idx);
					ss >> SelectedData;
					idxCol = -1;
					itVariationRow = VariationRow.begin();
					map <string, string>::iterator itAccession = Accessions.begin();
					for (; itVariationRow != VariationRow.end(); itVariationRow++) {
						idxCol++;
						
						//cout << PreviousValue << "\t" << *itVariationRow << endl;
						if (idxCol >= offset) {
							string data = VerifyData(*itVariationRow, isSNPData, VariationRow);
							if (data.compare(SelectedData) == 0) {
								SelectedAccessions.insert(itAccession->first);
							}
							itAccession++;
						}
					}	// end for each entry of the variation row
					
					//break;
				}
			}
		} else {
			map<string, int>::iterator itIndelCounts = IndelCounts.begin();
			for (; itIndelCounts != IndelCounts.end(); itIndelCounts++) {
				if (itIndelCounts->second == 1) {
					isIndvidualVariation = true;

					string SelectedData = itIndelCounts->first;
					idxCol = -1;
					itVariationRow = VariationRow.begin();
					map <string, string>::iterator itAccession = Accessions.begin();
					for (; itVariationRow != VariationRow.end(); itVariationRow++) {
						idxCol++;
						
						//cout << PreviousValue << "\t" << *itVariationRow << endl;
						if (idxCol >= offset) {
							string data = VerifyData(*itVariationRow, isSNPData, VariationRow);
							if (data.compare(SelectedData) == 0) {
								SelectedAccessions.insert(itAccession->first);
							}
							itAccession++;
						}
					}	// end for each entry of the variation row

					//break;
				}
			}
		}
		
		//cout << isIndvidualVariation << endl;
		if (!isIndvidualVariation) {
			VariationsToBeDeleted.push_back(itVariation->first);
		} else {
			set<string>::iterator it = SelectedAccessions.begin();
			string strSelectedAccessions = *it;
			it++;
			for (; it != SelectedAccessions.end(); it++) {
				strSelectedAccessions += ", " + *it;
			}
			itVariation->second.push_back(strSelectedAccessions);
		}


			
	} // for each row of the variation table
	
	vector<string>::iterator itVariationsToBeDeleted = VariationsToBeDeleted.begin();
	for (; itVariationsToBeDeleted != VariationsToBeDeleted.end(); itVariationsToBeDeleted++) {
		VariationTable.erase(VariationTable.find(*itVariationsToBeDeleted));		
	}
}

void FilterForGeographicVariation (map < string, vector <string> >& VariationTable, map<string, string>& Accessions, map<string, string>& GeographicLocations) {
		
	cout << "Filtering for geographic variations ..." << endl;

	vector<string> VariationsToBeDeleted;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		map<string, set<string> > LocalisedVariations;
		map<string, int> LocationAccessionCount;

		vector<string> VariationRow = itVariation->second;
		
		if (isSNPData) {
			string data = VariationRow[2]; // ambiguous reference
			if (data[0] == 'N') {
				VariationsToBeDeleted.push_back(itVariation->first);
				continue;
			}
		}
		
//		cout << itVariation->first;
		int idxCol = -1;
		vector<string>::iterator itVariationRow = VariationRow.begin();
		map< string, string >::iterator itAccessions = Accessions.begin();
		set<string> LocationsPresent; // locations present in the data -- used later to identify which location a variation corresponds to
		for (; itVariationRow != VariationRow.end(); itVariationRow++) {
			idxCol++;
			
			//cout << PreviousValue << "\t" << *itVariationRow << endl;
			if (idxCol >= offset) {
				string data = VerifyData(*itVariationRow, isSNPData, VariationRow);
				
				string Accession = itAccessions->first;
				string Location = GeographicLocations[Accession];
				
				// add this location to the list of present locations
				LocationsPresent.insert(Location);
				
//				cout << "\t" << Accession << ":" << Location << ":" << data;
				// increment the number of accession for this location
				map<string, int>::iterator itLocationAccessionCount = LocationAccessionCount.find(Location);
				if (itLocationAccessionCount == LocationAccessionCount.end()) {
					LocationAccessionCount[Location] = 1;
				} else {
					LocationAccessionCount[Location]++;
				}

				LocalisedVariations[Location].insert(data);

				itAccessions++; // move to the next accessions
			}
		} // end for each entry of the variation row
				
//		cout << endl;
		
		//cout << VariationRow[0] << "\t" << VariationRow[1] << "\t" <<VariationRow[2] << "\t" << BaseCounts[0] << "\t" << BaseCounts[1] << "\t" << BaseCounts[2] << "\t" << BaseCounts[3] << endl;
		bool isGeographicVariation = false;
		map<string, set<string> >::iterator itLocalisedVariations = LocalisedVariations.begin();
		set<string>::iterator itLocationsPresent = LocationsPresent.begin();
		set<string> SelectedLocations;
		for (; itLocalisedVariations != LocalisedVariations.end(); itLocalisedVariations++) {
			string Location = *itLocationsPresent;
			itLocationsPresent++;
//			cout << itLocalisedVariations->first << "\t" << itLocalisedVariations->second.size() << endl;
			// check if all accession of this location have the same variation (set size should be one)
			if (itLocalisedVariations->second.size() == 1) {
				// get the variation for all accessions of this location 
				string CurrentLocationData = *itLocalisedVariations->second.begin();
				// check if no other location has this variation
				bool OtherLocationHasThisVariation = false;
				map<string, set<string> >::iterator itLocalisedVariations1 = LocalisedVariations.begin();
				for (; itLocalisedVariations1 != LocalisedVariations.end(); itLocalisedVariations1++) {
					if (itLocalisedVariations == itLocalisedVariations1) { // ignore the same location
						continue;
					}
					
					set<string>::iterator it = itLocalisedVariations1->second.begin();
					for (; it != itLocalisedVariations1->second.end(); it++) {
						if (CurrentLocationData.compare(*it) == 0) {
							OtherLocationHasThisVariation = true;
							break;
						}
					}
					if (OtherLocationHasThisVariation) {
						break;
					}
				} // end for itLocalisedVariations1
				if (!OtherLocationHasThisVariation) {
					isGeographicVariation = true;
					
					SelectedLocations.insert(Location);

					//break;
				}
			} // end if all accession of this location have the same variation
		} // end for each itlocalisedVariation 
		
//		cout << isGeographicVariation << endl;
		if (!isGeographicVariation) {
			VariationsToBeDeleted.push_back(itVariation->first);
		} else {
			set<string>::iterator it = SelectedLocations.begin();
			string strSelectedLocations = *it;
			it++;
			for (; it != SelectedLocations.end(); it++) {
				strSelectedLocations += ", " + *it;
			}
			itVariation->second.push_back(strSelectedLocations);
		}
		
		
	} // for each row of the variation table
	
	vector<string>::iterator itVariationsToBeDeleted = VariationsToBeDeleted.begin();
	for (; itVariationsToBeDeleted != VariationsToBeDeleted.end(); itVariationsToBeDeleted++) {
		VariationTable.erase(VariationTable.find(*itVariationsToBeDeleted));		
	}
}

void GetMajorityVariation (map < string, vector <string> >& VariationTable, map<string, string>& Accessions, int MajorityCutoff) {
	
	cout << "Checking for majority variations ..." << endl;
	
	if (MajorityCutoff < Accessions.size() / 2) {
		cout << "Invalid majority cutoff value" << endl;
		exit(0);
	}
	
	vector<string> VariationsToBeDeleted;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		int idxCol = -1;
		int BaseCounts[TOTAL_BASES] = {0,0,0,0,0};
		map<string, int> IndelCounts;
		vector<string>::iterator itVariationRow = VariationRow.begin();
		vector<string>::iterator itFirstAccession, itLastAccession;
		for (; itVariationRow != VariationRow.end(); itVariationRow++) {
			idxCol++;
			
			//cout << PreviousValue << "\t" << *itVariationRow << endl;
			if (idxCol >= offset) {
				if (idxCol == offset) {
					itFirstAccession = itVariationRow;
				}
				itLastAccession = itVariationRow;

				string data = VerifyData(*itVariationRow, isSNPData, VariationRow);
				if (isSNPData) {
					BaseCounts[GetBaseIndex(data[0])]++;
				} else {
					IndelCounts[data]++;
				}
			}
		}	// end for each entry of the variation row
		
		//cout << VariationRow[0] << "\t" << VariationRow[1] << "\t" <<VariationRow[2] << "\t" << BaseCounts[0] << "\t" << BaseCounts[1] << "\t" << BaseCounts[2] << "\t" << BaseCounts[3] << endl;
		bool isMajorityVariation = false;
		string Data = "";
		if (isSNPData) {
			for (int idx = 0; idx < TOTAL_BASES; idx++) {
				if (BaseCounts[idx] >= MajorityCutoff) {
					isMajorityVariation = true;
					stringstream ss;
					ss << GetIndexBase(idx);
					ss >> Data;
					break;
				}
			}
		} else {
			map<string, int>::iterator itIndelCounts = IndelCounts.begin();
			for (; itIndelCounts != IndelCounts.end(); itIndelCounts++) {
				if (itIndelCounts->second >= MajorityCutoff) {
					isMajorityVariation = true;
					Data = itIndelCounts->first;
					break;
				}
			}
		}
		
		if (!isMajorityVariation) { // not a majority variation -> delete it
			VariationsToBeDeleted.push_back(itVariation->first);
		} else {
			// check if the base is same as the reference (for SNPs) or "-" (i.e. no indel)
			if ((isSNPData && Data.compare(VariationRow[2]) == 0) || (!isSNPData && Data.compare("-") == 0)) { // if so, delete the row
				VariationsToBeDeleted.push_back(itVariation->first);
			} else { // keep it, but remove all the columns with the majority consensus
				VariationRow.erase(itFirstAccession, itLastAccession + 1);
				VariationRow.push_back(Data);
				itVariation->second = VariationRow;
			}
		}
		
	} // for each row of the variation table
	
	vector<string>::iterator itVariationsToBeDeleted = VariationsToBeDeleted.begin();
	for (; itVariationsToBeDeleted != VariationsToBeDeleted.end(); itVariationsToBeDeleted++) {
		VariationTable.erase(VariationTable.find(*itVariationsToBeDeleted));		
	}
}

void ReadGFFFile (string GFF_File, map<string, map<pair<long, long>, int* > >& GFFData, map<string, map<int, string > >& GFFGeneName, int NoOfAccessions) {
	
	map<pair<long, long>, int* > RefGFFData;
	map<int, string > RefGFFGeneName;
	// gff file
	ifstream ifs_gff (GFF_File.c_str());
	string strGFF = "";
	string reference = "";
	string previous_reference = "";
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		} else if (strGFF[0] == '#') {
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		// get start and end coordinates
		long start = atol(entries[3].c_str());
		long end = atol(entries[4].c_str());
		
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFData;
				GFFGeneName[previous_reference] = RefGFFGeneName;
			}
			
			RefGFFData = GFFData[reference];
			RefGFFGeneName = GFFGeneName[reference];
			
			previous_reference = reference;
		}
		
		int* VarFrequency = new int[NoOfAccessions];
		for(int i = 0; i < NoOfAccessions; i++) {
			VarFrequency[i]= 0;
		}
		
		RefGFFData[pair<long, long>(start,end)] = VarFrequency;
		int StartPos = entries[8].find("ID=") + 3;
		int EndPos = entries[8].find_first_of(";", StartPos);
		RefGFFGeneName[start] = entries[8].substr(StartPos, EndPos - StartPos);
		
	} // end while ifs_gff is good
	GFFData[previous_reference] = RefGFFData;
	GFFGeneName[previous_reference] = RefGFFGeneName;
	
	ifs_gff.close();
	
}

void GetVariationsPerGene(string GFF_File, map < string, vector <string> >& VariationTable, string OUT_File,  map<string, string>& Accessions) {

	cout << "Reading the GFF File and Initialising the Variables ... " << endl;
	int NoOfAccessions = Accessions.size();
	//Read the GFF File
	// variable to store GFF Data
	map<string, map<pair<long, long>, int* > > GFFData;
	map<string, map<int, string > > GFFGeneName;
	ReadGFFFile(GFF_File, GFFData, GFFGeneName, NoOfAccessions);
	
	cout << "Counting the variations per gene ... " << endl;
	string reference = "";
	string previous_reference = "";
	map<pair<long, long>, int* > RefGFFData;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		reference = VariationRow[0];
		long position = atol(VariationRow[1].c_str());
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFData;
			}
			
			RefGFFData = GFFData[reference];
			
			previous_reference = reference;
		}
				
		// find the gene corresponding to this position
		map<pair<long, long>, int* >::iterator itRefGFFData = RefGFFData.begin();
		for (; itRefGFFData != RefGFFData.end(); itRefGFFData++) {
			//cout << itRefGFFData->first.first << "\t" << itRefGFFData->first.second << endl;		
			if (itRefGFFData->first.first - 20 <= position && itRefGFFData->first.second + 20 >= position) {
				break; 
			}
		}
		//cout << position << "\t" << itRefGFFData->first.first << "\t" << itRefGFFData->first.second << endl;		

		int idxCol = -1;
		vector<string>::iterator itVariationRow = VariationRow.begin();
		for (; itVariationRow != VariationRow.end(); itVariationRow++) {
			idxCol++;
			
			if (idxCol >= offset) {
				if ((*itVariationRow).length() != 0) { // the current variation is present in this accession  
					// increment the variation count for this accession
					itRefGFFData->second[idxCol - offset]++;
				}				
			}
		}	// end for each entry of the variation row
		
	} // for each row of the variation table
		
	cout << "Writing variation count to the output ..." << endl;
	
	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		
	
	ofs << "Reference\tUnigene\tStart\tEnd";
	map <string, string>::iterator itAccession = Accessions.begin();
	for (; itAccession != Accessions.end(); itAccession++) {
		string Accession = itAccession->first;
		ofs << "\t" << Accession;
	}
	ofs << endl;
	
	map<string, map<pair<long, long>, int* > >::iterator itGFFData = GFFData.begin();
	for (; itGFFData != GFFData.end(); itGFFData++) {
		map<pair<long, long>, int* >::iterator itRefGFFData = itGFFData->second.begin();
		map<int, string > RefGFFGeneName = GFFGeneName[itGFFData->first];
		for (; itRefGFFData != itGFFData->second.end(); itRefGFFData++) {
			ofs << itGFFData->first << "\t" << RefGFFGeneName[itRefGFFData->first.first] << "\t" << itRefGFFData->first.first << "\t" << itRefGFFData->first.second;
			for (int i = 0; i < NoOfAccessions; i++) {
				ofs << "\t" << itRefGFFData->second[i];
			}
			ofs	<< endl;
		}
	}

	ofs.close();
	
}

void GetUniqueVariationsPerGene(string GFF_File, map < string, vector <string> >& VariationTable, string OUT_File, string Reading_Frame_Coordinates_File) {
	
	cout << "Reading the GFF File and Initialising the Variables ... " << endl;
	//Read the GFF File
	// variable to store GFF Data
	map<string, map<pair<long, long>, int* > > GFFData;
	map<string, map<int, string > > GFFGeneName;
	ReadGFFFile(GFF_File, GFFData, GFFGeneName, 6); //CDS-Length, CDS-SNPs, 5'UTR-Length, 5'UTR-SNPs, 3'UTR-Length, 3'UTR-SNPs

	// read the blast output in the map
	BLAST_OUTPUT BlastOutput;
	ReadBlastOutput(Reading_Frame_Coordinates_File, BlastOutput);
	
	cout << "Counting unique variations per gene ... " << endl;
	string reference = "";
	string previous_reference = "";
	map<pair<long, long>, int* > RefGFFData;
	map<int, string > RefGFFGeneName;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		reference = VariationRow[0];
		long position = atol(VariationRow[1].c_str());
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFData;
			}
			
			RefGFFData = GFFData[reference];
			RefGFFGeneName = GFFGeneName[reference];
			
			previous_reference = reference;
		}
		
		// find the gene corresponding to this position
		map<pair<long, long>, int* >::iterator itRefGFFData = RefGFFData.begin();
		for (; itRefGFFData != RefGFFData.end(); itRefGFFData++) {
			//cout << itRefGFFData->first.first << "\t" << itRefGFFData->first.second << endl;		
			if (itRefGFFData->first.first - 20 <= position && itRefGFFData->first.second + 20 >= position) {
				break; 
			}
		}
		
		if (itRefGFFData == RefGFFData.end()) {
			continue;
		}
		
		string GeneName = RefGFFGeneName[itRefGFFData->first.first];
		BLAST_OUTPUT::iterator itBlastOutput = BlastOutput.find(GeneName);
		if (itBlastOutput == BlastOutput.end()) {
			// no reading frames/cds coordinates found
			// increment the variation count for this gene for CDS and mark the CDS length as -1
			itRefGFFData->second[0] = -1;
			itRefGFFData->second[1]++;
		} else {
			int GeneLength = itRefGFFData->first.second - itRefGFFData->first.first + 1;
			QUERY_OUTPUT QueryOutput = itBlastOutput->second;

			// get the hits
			vector<HIT>::iterator itHit = QueryOutput.Hits.begin();
			for (; itHit != QueryOutput.Hits.end(); itHit++) {
				HIT Hit = *itHit;
				HSP Hsp;
				if (Hit.Hsps.size() == 1) { // only one hit
					Hsp = *Hit.Hsps.begin();
					if (Hsp.HspFrame < 0) {
						int	HspFromOriginal = Hsp.HspFrom;
						Hsp.HspFrom = GeneLength - Hsp.HspTo + 1;
						Hsp.HspTo = GeneLength - HspFromOriginal + 1;
					}
					
				} else { // multiple HSPs
					// create a dummny Hsp
					Hsp.HspFrom =  numeric_limits<int>::max();
					Hsp.HspTo =  numeric_limits<int>::min();
					vector<HSP>::iterator itHsp = Hit.Hsps.begin();
					// set the HSP coordinates as follows. Everything to the left of 1st HSP is 5'UTR and everything to the right of last HSP is 3'UTR. 
					for (; itHsp != Hit.Hsps.end(); itHsp++) {
						int qStart = itHsp->HspFrom;
						int qEnd = itHsp->HspTo;
						if (itHsp->HspFrame < 0) {
							int	qStartOriginal = qStart;
							qStart = GeneLength - qEnd + 1;
							qEnd = GeneLength - qStartOriginal + 1;
						}
						
						Hsp.HspFrame = itHsp->HspFrame;
						Hsp.HspFrom = min(Hsp.HspFrom, qStart);
						Hsp.HspTo = max(Hsp.HspTo, qEnd);
					}
				}
				// CDS length
				itRefGFFData->second[0] = Hsp.HspTo - Hsp.HspFrom + 1;
				// 5'UTR length
				itRefGFFData->second[2] = Hsp.HspFrom - 1;
				// 3'UTR length
				itRefGFFData->second[4] = GeneLength - Hsp.HspTo + 1;
				
		/*		if (GeneName.find("17894659") != string::npos) {
					cout << position << "\t" << itRefGFFData->first.first << "\t" << itRefGFFData->first.second << "\t" << itRefGFFData->first.first + Hsp.HspFrom << "\t" << itRefGFFData->first.first + Hsp.HspTo << endl;
				}
				if (position >= 228000 && position <= 228344) {
					cout << position << endl;
				}
		*/		
				
				if (Hsp.HspFrame > 0) { // this gene is in the forward direction 
					if (position < itRefGFFData->first.first + Hsp.HspFrom) { // SNP is before the start of CDS (5'UTR)
						itRefGFFData->second[3]++;
					} else if (position < itRefGFFData->first.first + Hsp.HspTo) { // SNP is in the CDS
						itRefGFFData->second[1]++;
					} else { // SNP is after the CDS (3'UTR)
						itRefGFFData->second[5]++;
					}
				} else { // this gene is in the reverse direction (|-3'UTR-|--CDS--|-5'UTR-|)
					if (position >= itRefGFFData->first.second - itRefGFFData->second[2]) { // SNP is before the start of CDS (5'UTR)
						itRefGFFData->second[3]++;
					} else if (position > itRefGFFData->first.first + itRefGFFData->second[4]) { // SNP is in the CDS
						itRefGFFData->second[1]++;
					} else { // SNP is after the CDS (3'UTR)
						itRefGFFData->second[5]++;
					}
					
				}
				
			}
			
		}
		
	} // for each row of the variation table
	
	cout << "Writing variation count to the output ..." << endl;
	
	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		
	
	ofs << "Reference\tUnigene\tStart\tEnd\tCDS-Length\tCDS-Variations\t5'UTR-Length\t5'UTR-Variations\t3'UTR-Length\t3'UTR-Variations" << endl;
	
	map<string, map<pair<long, long>, int* > >::iterator itGFFData = GFFData.begin();
	for (; itGFFData != GFFData.end(); itGFFData++) {
		map<pair<long, long>, int* >::iterator itRefGFFData = itGFFData->second.begin();
		map<int, string > RefGFFGeneName = GFFGeneName[itGFFData->first];
		for (; itRefGFFData != itGFFData->second.end(); itRefGFFData++) {
			string GeneName = RefGFFGeneName[itRefGFFData->first.first];
			ofs << itGFFData->first << "\t" << GeneName << "\t" << itRefGFFData->first.first << "\t" << itRefGFFData->first.second;
			
			if (itRefGFFData->second[0] == 0) { // No variation was found for this gene, so save the CDS, 5'UTR and 3'UTR lengths
				BLAST_OUTPUT::iterator itBlastOutput = BlastOutput.find(GeneName);
				if (itBlastOutput == BlastOutput.end()) {
					// no reading frames/cds coordinates found
					// increment the variation count for this gene for CDS and mark the CDS length as -1
					itRefGFFData->second[0] = -1;
				} else {
					int GeneLength = itRefGFFData->first.second - itRefGFFData->first.first + 1;
					QUERY_OUTPUT QueryOutput = itBlastOutput->second;
					
					// get the hits
					vector<HIT>::iterator itHit = QueryOutput.Hits.begin();
					for (; itHit != QueryOutput.Hits.end(); itHit++) {
						HIT Hit = *itHit;
						HSP Hsp;
						if (Hit.Hsps.size() == 1) { // only one hit
							Hsp = *Hit.Hsps.begin();
							if (Hsp.HspFrame < 0) {
								int	HspFromOriginal = Hsp.HspFrom;
								Hsp.HspFrom = GeneLength - Hsp.HspTo + 1;
								Hsp.HspTo = GeneLength - HspFromOriginal + 1;
							}
						} else { // multiple HSPs
							// create a dummny Hsp
							Hsp.HspFrom =  numeric_limits<int>::max();
							Hsp.HspTo =  numeric_limits<int>::min();
							vector<HSP>::iterator itHsp = Hit.Hsps.begin();
							// set the HSP coordinates as follows. Everything to the left of 1st HSP is 5'UTR and everything to the right of last HSP is 3'UTR. 
							for (; itHsp != Hit.Hsps.end(); itHsp++) {
								int qStart = itHsp->HspFrom;
								int qEnd = itHsp->HspTo;
								if (itHsp->HspFrame < 0) {
									int	qStartOriginal = qStart;
									qStart = GeneLength - qEnd + 1;
									qEnd = GeneLength - qStartOriginal + 1;
								}
								
								Hsp.HspFrom = min(Hsp.HspFrom, qStart);
								Hsp.HspTo = max(Hsp.HspTo, qEnd);
							}
						}
						// CDS length
						itRefGFFData->second[0] = Hsp.HspTo - Hsp.HspFrom + 1;
						// 5'UTR length
						itRefGFFData->second[2] = Hsp.HspFrom - 1;
						// 3'UTR length
						itRefGFFData->second[4] = GeneLength - Hsp.HspTo + 1;						
					}
					
				}			
			}
			
			for (int i = 0; i < 6; i++) {
				ofs << "\t" << itRefGFFData->second[i];
			}
			ofs << endl;
		}
	}
	
	ofs.close();
	
}

void ReadSequenceFile(string Sequence_File, map<string, string>& Sequences){

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
				Sequences[Id] = dna;
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
	// add the last sequence
	Sequences[Id] = dna;
}

bool CheckSNP(HSP& Hsp, int GeneLength, int& RelativePosition, string& dna, string SNPBase, bool& isSynonymous) {
	//returns true if its a snp in the coding region otherwise returns false
	
	// get the relative position of this snp
	if (Hsp.HspFrame < 0) {
		int	HspFromOriginal = Hsp.HspFrom;
		Hsp.HspFrom = GeneLength - Hsp.HspTo + 1;
		Hsp.HspTo = GeneLength - HspFromOriginal + 1;
		RelativePosition = GeneLength - RelativePosition + 1;
		dna = ReverseComplementDNA(dna);
		Hsp.HspFrame *= -1;
		//cout << "\t" << Hsp.HspFrom << "\t" << Hsp.HspTo << "\t" << RelativePosition << endl;
	}
	if (RelativePosition >= Hsp.HspFrom && RelativePosition <= Hsp.HspTo) { // this snp is in the coding region

		// get the position of the base in the codon
		int PositionInCodon = RelativePosition % 3;
		if (PositionInCodon == 0) {
			PositionInCodon = 3;
		}
		// extract the codon from the dna
		string Codon = dna.substr(RelativePosition - PositionInCodon, 3);
		// the new codon (with snp base)
		string NewCodon = Codon;
		NewCodon[PositionInCodon - 1] = SNPBase[0];
		
		//cout << "\t" << PositionInCodon << "\t" << dna.substr(RelativePosition - PositionInCodon, 3) << "\t" << Codon  << "\t" << TranslateDNASequence(Codon, 0) << "\t" << SNPBase << "\t" << TranslateDNASequence(NewCodon, 0) << "\t" << NewCodon << endl;
		// Check if it is synonymous or non synonymous snp
		if ( TranslateDNASequence(Codon, 0).compare(TranslateDNASequence(NewCodon, 0)) == 0) { // check if both codons give the same aminoacid
			isSynonymous = true;
		} else {
			isSynonymous = false;
		}
		
		return true;
	} else {
		return false;
	}

	
}

// need to update from SynonymousVsNonSynonymousUniqueVariationsPerGene()
void SynonymousVsNonSynonymousVariationsPerGene(string GFF_File, map < string, vector <string> >& VariationTable, string OUT_File, string Reading_Frame_Coordinates_File, map<string, string>& Accessions, string Sequence_File) {
	
	cout << "Reading the GFF File and Initialising the Variables ... " << endl;
	int NoOfAccessions = Accessions.size();
	//Read the GFF File
	// variable to store GFF Data
	map<string, map<pair<long, long>, int* > > GFFData;
	map<string, map<int, string > > GFFGeneName;
	ReadGFFFile(GFF_File, GFFData, GFFGeneName, NoOfAccessions*3); //(Synonymous, Non-Synonymous, Unassigned)

	cout << "Reading the Sequence File ... " << endl;
	map<string, string> Sequences;
	ReadSequenceFile(Sequence_File, Sequences);
	
	// read the blast output in the map
	BLAST_OUTPUT BlastOutput;
	ReadBlastOutput(Reading_Frame_Coordinates_File, BlastOutput);
	
	cout << "Counting synonymous and non synonymous SNPs per gene ... " << endl;
	string reference = "";
	string previous_reference = "";
	map<pair<long, long>, int* > RefGFFData;
	map<int, string > RefGFFGeneName;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		reference = VariationRow[0];
		long position = atol(VariationRow[1].c_str());
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFData;
			}
			
			RefGFFData = GFFData[reference];
			RefGFFGeneName = GFFGeneName[reference];
			
			previous_reference = reference;
		}
		
		// find the gene corresponding to this position
		map<pair<long, long>, int* >::iterator itRefGFFData = RefGFFData.begin();
		for (; itRefGFFData != RefGFFData.end(); itRefGFFData++) {
			//cout << itRefGFFData->first.first << "\t" << itRefGFFData->first.second << endl;		
			if (itRefGFFData->first.first - 20 <= position && itRefGFFData->first.second + 20 >= position) {
				break; 
			}
		}
		
		if (itRefGFFData == RefGFFData.end()) {
			continue;
		}
		
		string GeneName = RefGFFGeneName[itRefGFFData->first.first];
		BLAST_OUTPUT::iterator itBlastOutput = BlastOutput.find(GeneName);

		if (itBlastOutput == BlastOutput.end()) {
			// no reading frames/cds coordinates found
			// go through each accession
		/*	int idxCol = -1;
			vector<string>::iterator itVariationRow = VariationRow.begin();
			for (; itVariationRow != VariationRow.end(); itVariationRow++) {
				idxCol++;
				
				if (idxCol >= offset) {
					if ((*itVariationRow).length() != 0) { // the current variation is present in this accession  
						// increment the unassinged snp count for this accession
						itRefGFFData->second[(idxCol - offset)*3 + 2]++;
						
					} // end if the variation is present in the current accession
					
				} 
			}	// end for each entry of the variation row
		*/	
			//cout << position << ": No Blast Hit found for " << GeneName << endl  ;
		} else {
			// get the query output for this unigene
			QUERY_OUTPUT QueryOutput = itBlastOutput->second;
			// get the hit
			HIT Hit = *QueryOutput.Hits.begin();
			
			// gene length
			int GeneLength = itRefGFFData->first.second - itRefGFFData->first.first + 1;
			// get the dna sequence for this unigene
			string dna = Sequences[itBlastOutput->first];
			// position of the SNP relative to start of the gene
			int RelativePosition = position - itRefGFFData->first.first + 1;

			if (RelativePosition <= 0 || RelativePosition > GeneLength) { // outside the gene coordinates 
				continue;
			}
			
			//cout << position << "\t" << RelativePosition << "\t" << Hit.Hsps.size() << "\t" << GeneName << "\t" << Hit.Hsps.begin()->HspFrame << endl;

			// go through each accession
			int idxCol = -1;
			vector<string>::iterator itVariationRow = VariationRow.begin();
			for (; itVariationRow != VariationRow.end(); itVariationRow++) {
				idxCol++;
				
				if (idxCol >= offset) {
					if ((*itVariationRow).length() != 0) { // the current variation is present in this accession  
						bool isSynonymous = false;
						int nHsps = 0;
						vector<HSP>::iterator itHsp = Hit.Hsps.begin();
						for (; itHsp != Hit.Hsps.end(); itHsp++) {							
							// check if the snp is synonymous or non-synonymous
							nHsps += CheckSNP(*itHsp, GeneLength, RelativePosition, dna, *itVariationRow, isSynonymous);
						} // for each HSP
						if (nHsps > 1) { // this SNP is found in multiple Hsps 
							// increment the unassinged snp count for this accession
							itRefGFFData->second[(idxCol - offset)*3 + 2]++;
							//cout << "\tUnassigned" << endl;
						} else {
							if (isSynonymous) {
								// increment the synonymous snp count for this accession
								itRefGFFData->second[(idxCol - offset)*3]++;
								//cout << "\tSynonymous" << endl;
							} else {
								// increment the non-synonymous snp count for this accession
								itRefGFFData->second[(idxCol - offset)*3 + 1]++;
								//cout << "\tNon-Synonymous" << endl;
							}
						}

					} // end if the variation is present in the current accession
					
				} 
			}	// end for each entry of the variation row
			
			
		} // if there was a blast output for this unigene
		
	} // for each row of the variation table
	
	cout << "Writing variation count to the output ..." << endl;
	
	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		
	
	ofs << "Reference\tUnigene\tStart\tEnd";
	map <string, string>::iterator itAccession = Accessions.begin();
	for (; itAccession != Accessions.end(); itAccession++) {
		string Accession = itAccession->first;
		ofs << "\t" << Accession << "-Synonymous" << "\t" << Accession << "-Non-Synonymous" << "\t" << Accession << "-Unassinged";
	}
	ofs << endl;
	
	map<string, map<pair<long, long>, int* > >::iterator itGFFData = GFFData.begin();
	for (; itGFFData != GFFData.end(); itGFFData++) {
		map<pair<long, long>, int* >::iterator itRefGFFData = itGFFData->second.begin();
		map<int, string > RefGFFGeneName = GFFGeneName[itGFFData->first];
		for (; itRefGFFData != itGFFData->second.end(); itRefGFFData++) {
			string GeneName = RefGFFGeneName[itRefGFFData->first.first];
			ofs << itGFFData->first << "\t" << GeneName << "\t" << itRefGFFData->first.first << "\t" << itRefGFFData->first.second;
			
			for (int i = 0; i < NoOfAccessions * 3; i++) {
				ofs << "\t" << itRefGFFData->second[i];
			}
			ofs << endl;
		}
	}
	
	ofs.close();
	
}

void SynonymousVsNonSynonymousUniqueVariationsPerGene(string GFF_File, map < string, vector <string> >& VariationTable, string OUT_File, string Reading_Frame_Coordinates_File, map<string, string>& Accessions, string Sequence_File) {
	
	cout << "Reading the GFF File and Initialising the Variables ... " << endl;
	int NoOfAccessions = Accessions.size();
	//Read the GFF File
	// variable to store GFF Data
	map<string, map<pair<long, long>, int* > > GFFData;
	map<string, map<int, string > > GFFGeneName;
	ReadGFFFile(GFF_File, GFFData, GFFGeneName, 4); //(Synonymous, Non-Synonymous, Unassigned-Multiple Frame, Unassigned-Inconsistent)
	
	cout << "Reading the Sequence File ... " << endl;
	map<string, string> Sequences;
	ReadSequenceFile(Sequence_File, Sequences);
	
	// read the blast output in the map
	BLAST_OUTPUT BlastOutput;
	ReadBlastOutput(Reading_Frame_Coordinates_File, BlastOutput);
	
	cout << "Counting synonymous and non synonymous unique SNPs per gene ... " << endl;
	string reference = "";
	string previous_reference = "";
	map<pair<long, long>, int* > RefGFFData;
	map<int, string > RefGFFGeneName;
	map <string, vector<string> >::iterator itVariation = VariationTable.begin();
	for (; itVariation != VariationTable.end(); itVariation++) {
		vector<string> VariationRow = itVariation->second;
		
		reference = VariationRow[0];
		long position = atol(VariationRow[1].c_str());
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFData;
			}
			
			RefGFFData = GFFData[reference];
			RefGFFGeneName = GFFGeneName[reference];
			
			previous_reference = reference;
		}
		
		// find the gene corresponding to this position
		map<pair<long, long>, int* >::iterator itRefGFFData = RefGFFData.begin();
		for (; itRefGFFData != RefGFFData.end(); itRefGFFData++) {
			//cout << itRefGFFData->first.first << "\t" << itRefGFFData->first.second << endl;		
			if (itRefGFFData->first.first - 20 <= position && itRefGFFData->first.second + 20 >= position) {
				break; 
			}
		}
		
		if (itRefGFFData == RefGFFData.end()) {
			continue;
		}
		
		string GeneName = RefGFFGeneName[itRefGFFData->first.first];
		BLAST_OUTPUT::iterator itBlastOutput = BlastOutput.find(GeneName);
		
		if (itBlastOutput == BlastOutput.end()) {
			// no reading frames/cds coordinates found
			// increment the unassinged snp count for this accession
			//itRefGFFData->second[2]++;
			
			//cout << position << ": No Blast Hit found for " << GeneName << endl  ;
		} else {
			// get the query output for this unigene
			QUERY_OUTPUT QueryOutput = itBlastOutput->second;
			// get the hit
			HIT Hit = *QueryOutput.Hits.begin();
			
			// gene length
			int GeneLength = itRefGFFData->first.second - itRefGFFData->first.first + 1;
			// get the dna sequence for this unigene
			string dna = Sequences[itBlastOutput->first];
			// position of the SNP relative to start of the gene
			int RelativePosition = position - itRefGFFData->first.first + 1;
			
			if (RelativePosition <= 0 || RelativePosition > GeneLength) { // outside the gene coordinates 
				continue;
			}
			
			//cout << position << "\t" << RelativePosition << "\t" << Hit.Hsps.size() << "\t" << GeneName << "\t" << Hit.Hsps.begin()->HspFrame << endl;
			
			// go through each accession
			int idxCol = -1;
			bool isSynonymous = false;
			int nSynonymous = 0;
			int nNonSynonymous = 0;
			bool isInMultipleFrame = false;
			//bool isInconsisent = false;
			bool isInCDS = false;
			vector<string>::iterator itVariationRow = VariationRow.begin();
			for (; itVariationRow != VariationRow.end(); itVariationRow++) {
				idxCol++;
				if (idxCol >= offset) {
					if ((*itVariationRow).length() != 0) { // the current variation is present in this accession  
						int nHsps = 0;
						vector<HSP>::iterator itHsp = Hit.Hsps.begin();
						for (; itHsp != Hit.Hsps.end(); itHsp++) {							
							// check if the snp is synonymous or non-synonymous
							nHsps += CheckSNP(*itHsp, GeneLength, RelativePosition, dna, *itVariationRow, isSynonymous);
						} // for each HSP
						if (nHsps > 0) { // the snp was found in CDS
							isInCDS = true;
							if (nHsps > 1) { // this SNP is found in multiple Hsps 
								isInMultipleFrame = true;
								break;
								//cout << "\tUnassigned" << endl;
							} else {
								if (isSynonymous) {
									nSynonymous++;
								} else {
									nNonSynonymous++;
								}
							}
						}						
					} // end if the variation is present in the current accession
				} 
			}	// end for each entry of the variation row
			if (isInCDS) {
				if (isInMultipleFrame) { // multiple frame
					// increment the in-multiple-frame snp count for this gene
					itRefGFFData->second[2]++;
				} else if (nSynonymous > 0 && nNonSynonymous > 0) { // this snp is synonymous is one or more lines and also non-synouous in one or more other lines
					// increment the inconsistent snp count for this gene
					itRefGFFData->second[3]++;
				} else if (nSynonymous > 0) {
					//cout << "Synonymous" << endl;
					// increment the synonymous snp count for this gene
					itRefGFFData->second[0]++;
				} else if (nNonSynonymous > 0) {
					//cout << "Non-Synonymous" << endl;
					// increment the non-synonymous snp count for this gene
					itRefGFFData->second[1]++;
				}	
			}
			
		} // if there was a blast output for this unigene
		
	} // for each row of the variation table
	
	cout << "Writing variation count to the output ..." << endl;
	
	ofstream ofs ( OUT_File.c_str() ); // open the output file stream		
	
	ofs << "Reference\tUnigene\tStart\tEnd\tSynonymous SNPs\tNon-synonymous SNPs\tUnassigned-Multiple-Frame SNPs\tUnassigned-Inconsistent SNPs";
	ofs << endl;
	
	map<string, map<pair<long, long>, int* > >::iterator itGFFData = GFFData.begin();
	for (; itGFFData != GFFData.end(); itGFFData++) {
		map<pair<long, long>, int* >::iterator itRefGFFData = itGFFData->second.begin();
		map<int, string > RefGFFGeneName = GFFGeneName[itGFFData->first];
		for (; itRefGFFData != itGFFData->second.end(); itRefGFFData++) {
			string GeneName = RefGFFGeneName[itRefGFFData->first.first];
			ofs << itGFFData->first << "\t" << GeneName << "\t" << itRefGFFData->first.first << "\t" << itRefGFFData->first.second;
			
			for (int i = 0; i < 4; i++) {
				ofs << "\t" << itRefGFFData->second[i];
			}
			ofs << endl;
		}
	}
	
	ofs.close();
	
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Variation (argv[1]);
	// initialisation - Global variables
	offset = 2;
	isSNPData = false;
	isInsData = false;
	int ParameterOffset = 0;
	if (Variation[0] == 'S' || Variation[0] == 's') {
		offset = 3;
		isSNPData = true;
	} else if (Variation[0] == 'I' || Variation[0] == 'i') {
		isInsData = true;
	}
	
	string Accession_File (argv[2]);
	string Base_Dist_File (argv[3]);
	string Coverage_File;
	if (!isSNPData) {
		Coverage_File =  (argv[4]);
		ParameterOffset = 1;
	}
	string OUT_File (argv[4 + ParameterOffset]);
	int AnalysisType (atoi(argv[5 + ParameterOffset]));
	
	
	map <string, string> Accessions;
	ReadFile2(Accession_File, Accessions);

	map < string, vector <string> > VariationTable;
	TabulateVariation (Accessions, VariationTable);
	cout << "Number of Entries: " << VariationTable.size() << endl;

//	if (AnalysisType != MAJORITY_VARIATIONS && AnalysisType != UNIQUE_VARIATIONS_PER_GENE && Base_Dist_File.length() != 0) {
	if (AnalysisType != MAJORITY_VARIATIONS && Base_Dist_File.length() != 0) {
		map <string, string> BaseDistributionsFiles;
		ReadFile2(Base_Dist_File, BaseDistributionsFiles);
		if (isSNPData) {
			CheckBaseDistribution(BaseDistributionsFiles, Accessions, VariationTable);
		} else {
			map <string, string> CoverageFiles;
			ReadFile2(Coverage_File, CoverageFiles);
			CheckCoverage(BaseDistributionsFiles, CoverageFiles, Accessions, VariationTable);			
		}
	}
		
	string Geographic_Location_File;
	string GFF_File;
	map <string, string> GeographicLocations;
	int MajorityCutoff = 0;
	string strRemoveCommonVariations = "";
	string Reading_Frame_Coordinates_File = "";
	string Sequence_File = "";
	switch (AnalysisType) {
		case TABULATE_VARIATIONS:
			strRemoveCommonVariations = argv[6 + ParameterOffset];
			if (strRemoveCommonVariations[0] == 'T' || strRemoveCommonVariations[0] == 't') {
				RemoveCommonVariation(VariationTable, Accessions);
			}
			WriteTable(OUT_File, Accessions, VariationTable);
			break;

		case INDIVIDUAL_VARIATIONS:
			FilterForIndividualVariation(VariationTable, Accessions);
			WriteTable(OUT_File, Accessions, VariationTable);
			break;
			
		case GEOGRAPHIC_VARIATIONS:
			Geographic_Location_File = argv[6 + ParameterOffset];
			ReadFile2(Geographic_Location_File, GeographicLocations);

			FilterForGeographicVariation(VariationTable, Accessions, GeographicLocations);
			WriteTable(OUT_File, Accessions, VariationTable);
			break;

		case MAJORITY_VARIATIONS:
			MajorityCutoff = atoi(argv[6 + ParameterOffset]);
			GetMajorityVariation(VariationTable, Accessions, MajorityCutoff);
			WriteTable(OUT_File, Accessions, VariationTable);
			break;
			
		case VARIATIONS_PER_GENE:
			GFF_File = argv[6 + ParameterOffset];
			GetVariationsPerGene(GFF_File, VariationTable, OUT_File, Accessions);
			break;

		case UNIQUE_VARIATIONS_PER_GENE:
			GFF_File = argv[6 + ParameterOffset];
			strRemoveCommonVariations = argv[7 + ParameterOffset];
			if (strRemoveCommonVariations[0] == 'T' || strRemoveCommonVariations[0] == 't') {
				RemoveCommonVariation(VariationTable, Accessions);
			}
			Reading_Frame_Coordinates_File = argv[8 + ParameterOffset];
			GetUniqueVariationsPerGene(GFF_File, VariationTable, OUT_File, Reading_Frame_Coordinates_File);
			break;
			
		case SYNONYMOUS_VS_NON_SYNONYMOUS_SNP:
			GFF_File = argv[6 + ParameterOffset];
			strRemoveCommonVariations = argv[7 + ParameterOffset];
			if (strRemoveCommonVariations[0] == 'T' || strRemoveCommonVariations[0] == 't') {
				RemoveCommonVariation(VariationTable, Accessions);
			}
			Reading_Frame_Coordinates_File = argv[8 + ParameterOffset];
			Sequence_File = argv[9 + ParameterOffset];
			SynonymousVsNonSynonymousVariationsPerGene(GFF_File, VariationTable, OUT_File, Reading_Frame_Coordinates_File, Accessions, Sequence_File);
			break;

		case SYNONYMOUS_VS_NON_SYNONYMOUS_UNIQUE_SNP:
			GFF_File = argv[6 + ParameterOffset];
			strRemoveCommonVariations = argv[7 + ParameterOffset];
			if (strRemoveCommonVariations[0] == 'T' || strRemoveCommonVariations[0] == 't') {
				RemoveCommonVariation(VariationTable, Accessions);
			}
			Reading_Frame_Coordinates_File = argv[8 + ParameterOffset];
			Sequence_File = argv[9 + ParameterOffset];
			SynonymousVsNonSynonymousUniqueVariationsPerGene(GFF_File, VariationTable, OUT_File, Reading_Frame_Coordinates_File, Accessions, Sequence_File);
			break;
			
		default:
			break;
	}
	
	return 0;
}
