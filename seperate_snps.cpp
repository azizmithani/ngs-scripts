/*
 *
 *  Created by Aziz Mithani on 29/07/2010.
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

const int MIN_READ_COUNT = 2;
const double MIN_READ_PERCENTAGE = 5.0;
const double GENOME_MEMBERSHIP_THRESHOLD = 0.75;
const int MIN_DIPLOID_COVERAGE = 5;

const int MIN_BASE_QUALITY = 20;
const int OPTIONAL_FIELD_START = 11;
const int FLAG_MATE_UNMAPPED = 8;

typedef pair<long, char> SNP;
typedef map < string, char* > SNP_LIST;
typedef map < string, int* > COVERAGE_LIST;

map <char, string> NucleotideMap;

class SNPGroup {
public:
	set< SNP > SNPList; // SNPs in this Group
	map < SNP, int > SNPReadCount; // No. of reads supporting a SNP
//	int ReadCount;
	bool GroupActive;
	int OverlapSize;
	
	bool AGenome;
	bool BGenome;
	bool DGenome;
	double ASNPPercentage;
	double BSNPPercentage;
	double DSNPPercentage;
};

void InitialiseSNPList (string SAM_File, SNP_LIST& SNPList, map < string, char* >& RefBaseList, map < string, long >& RefLength) {
	// read the header from SAM File and initialise the variable to hold SNP List
	string str="";
	// open the sam file stream
	ifstream ifs ( SAM_File.c_str() );
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			long length = atol(strLength.c_str());
			
			// memory initialisation
			char *RefSNPList = new char[length];
			char *RefRefBaseList = new char[length];
			for(long i = 0; i < length; i++) {
				*(RefSNPList+i) = ' ';		
				*(RefRefBaseList+i) = ' ';		
			}
			SNPList[reference] = RefSNPList;
			RefBaseList[reference] = RefRefBaseList;
			
			RefLength[reference] = length;
		} else {
			break;
		}
		
	}
	ifs.close();
	
}

void InitialiseSNPList (map < string, long >& RefLength, SNP_LIST& SNPList) {
	
	map < string, long >::iterator itRefLength = RefLength.begin();
	for (; itRefLength != RefLength.end(); itRefLength++) {
		
		string reference = itRefLength->first;
		long length = itRefLength->second;
		
		// memory initialisation
		char *RefSNPList = new char[length];
		for(long i = 0; i < length; i++) {
			*(RefSNPList+i) = ' ';		
		}
		SNPList[reference] = RefSNPList;
		
	}
	
}

void ReadSNPList (string SNP_File, SNP_LIST& SNPList) {
	
	// open the snp file stream
	ifstream ifs_snp ( SNP_File.c_str() );
	
	int offset = -1;
	string strSNP = "";
	string reference;
	string PreviousReference = "";
	char* RefSNPList;
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		reference = entries[0];
		long pos = atol(entries[1].c_str());
		string ConsensusBase = entries[3+offset];
		//string ReferenceBase = entries[2+offset];
		
		if (PreviousReference.compare(reference) != 0) {
			map< string, char* >::iterator itSNPList = SNPList.find(reference);
			// check if it is empty
			if (itSNPList !=  SNPList.end()) { // if not then get it
				RefSNPList = itSNPList->second;
			}		
			
			PreviousReference = reference;
		}
		
		*(RefSNPList + pos - 1) = ConsensusBase[0];
		
	}
	// close the SNP file
	ifs_snp.close();
	
}

void ReadSNPList (string SNP_File, SNP_LIST& SNPList, map < string, char* >& RefBaseList) {
	
	// open the snp file stream
	ifstream ifs_snp ( SNP_File.c_str() );
	
	int offset = -1;
	string strSNP = "";
	string reference;
	string PreviousReference = "";
	char* RefSNPList;
	char* RefRefBaseList;
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		reference = entries[0];
		long pos = atol(entries[1].c_str());
		string ConsensusBase = entries[3+offset];
		string ReferenceBase = entries[2+offset];
		
		if (PreviousReference.compare(reference) != 0) {
			map< string, char* >::iterator itSNPList = SNPList.find(reference);
			map< string, char* >::iterator itRefBaseList = RefBaseList.find(reference);
			// check if it is empty
			if (itSNPList !=  SNPList.end()) { // if not then get it
				RefSNPList = itSNPList->second;
			}					
			// check if it is empty
			if (itRefBaseList !=  RefBaseList.end()) { // if not then get it
				RefRefBaseList = itRefBaseList->second;
			}		
			PreviousReference = reference;
		}
		
		*(RefSNPList + pos - 1) = ConsensusBase[0];
		*(RefRefBaseList + pos - 1) = ReferenceBase[0];
		
	}
	// close the SNP file
	ifs_snp.close();
	
}

void ReadCoverageFile (string Coverage_File, map < string, long >& RefLength, COVERAGE_LIST& CoverageList) {
	
	// open the coverage file stream
	ifstream ifs_Coverage ( Coverage_File.c_str() );
	
	string reference;
	string PreviousReference = "";
	string str="";
	int* RefCoverageList;
	while (ifs_Coverage.good()) {
		getline(ifs_Coverage, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		reference = entries[0];
		long pos = atol(entries[1].c_str());
		int coverage = atoi(entries[2].c_str());
		
		if (PreviousReference.compare(reference) != 0) {
			if (PreviousReference.length() > 0) {
				CoverageList[PreviousReference] = RefCoverageList;
			}
			
			long length = RefLength[reference];
			// memory initialisation
			RefCoverageList = new int[length];
			for(long i = 0; i < length; i++) {
				*(RefCoverageList+i) = 0;		
			}
			CoverageList[reference] = RefCoverageList;
			
			PreviousReference = reference;
		}

		*(RefCoverageList + pos - 1) = coverage;
		
	}
	CoverageList[PreviousReference] = RefCoverageList;
	// close the SNP file
	ifs_Coverage.close();
		
}
		
int GetReadEndCoordinates (long ReadStart, string cigar) {
	
	long sum = 0;
	string length = "";
	for (int i = 0; i < cigar.length(); i++) {
		if (isdigit(cigar[i])) {
			length += cigar[i];
		} else if (cigar[i]=='S') {
			// ignore
			length="";
		} else if (cigar[i] =='D') {
			sum+= atol(length.c_str());
			length="";
		} else if (cigar[i] =='I') {
			length = "";
		} else if (cigar[i]=='M' || cigar[i]=='N') {
			sum += atol(length.c_str());
			length="";
		}
	} // end for
	
	return ReadStart + sum - 1;
}

void GetOptionalFields (vector<string>& entries, map<string, string>& Fields) {
	for (int i = OPTIONAL_FIELD_START; i < entries.size(); i++) {
		vector<string> f;
		strsplit(entries[i], f, ":");
		Fields[f[0]] = f[2];
	}
	
	//return Fields;
}

void GetSNPsInTheReadRegion (char* RefSNPList, long ReadStart, long ReadEnd, map<long, char>& ReadRegionSNPList) {//, vector<long>& ReadSNPPositions) {
	
	for (long l = ReadStart - 1; l < ReadEnd; l++) { // we start from ReadStart - 1 because SNP list is 0-based
		if ( *(RefSNPList+l) != ' ' ) {
			long pos = l+1; // get the actual position (1-based)
			ReadRegionSNPList[pos] = *(RefSNPList+l);
			//ReadSNPPositions.push_back(pos);
		}
	}
	
}

int GetRelativePosition (long position, long ReadStart, string MD) {
	
	bool inDeletion = false;
	string length = "";
	long sum = ReadStart;
	int RelativePosition = (int)(position - ReadStart);
	int Gap = 0;
	for (int i = 0; i < MD.length(); i++) {
		if (isdigit(MD[i])) {
			if (inDeletion)
				inDeletion = false;
			length += MD[i];
		} else if (MD[i]=='A' || MD[i]=='C' || MD[i]=='G' || MD[i]=='T') {
			if (!inDeletion) {
				sum += atol(length.c_str());
				//RelativePosition += atoi(length.c_str());
				//RelativePosition++;
				if (sum >= position)
					break;
			} else {
				if (sum == position) {
					return -1;
				}
				Gap++;
			}
			
			sum++; // increment the sum by one to include the current position 
			//RelativePosition++; // increment the RelativePosition by one to include the current position 
			length = "";
			
			
		} else if (MD[i] =='^') {
			sum += atol(length.c_str());
			if (sum > position)
				break;
			else if (sum == position) {
				return -1;
			}
			//RelativePosition += atoi(length.c_str());
			length="";
			inDeletion = true;
		}
		
	} // end for
	return RelativePosition - Gap;
}

bool GetReadSNPs (map < long, char >& ReadRegionSNPList, long& ReadStart, string& MD, string& ReadSeq, string& BaseQualities, set< SNP >& ReadSNPList) {
	// returns true if all snp bases are good quality and false if bad quality snps (which are ignored) are present	
	bool AllGoodQualitySNPs = true;
	map < long, char >::iterator itReadRegionSNPList = ReadRegionSNPList.begin();
	for (; itReadRegionSNPList != ReadRegionSNPList.end(); itReadRegionSNPList++) {
		// Get SNP details
		long SNPPosition = itReadRegionSNPList->first;
		char ConsensusBase = itReadRegionSNPList->second;
		
		// Get the relative position of this SNP on the read 
		int SNPRelativePosition = GetRelativePosition(SNPPosition, ReadStart, MD);
		
		// ignore this base if there's is a deletion in this read at SNP position (SNPRelativePosition == -1)
		if (SNPRelativePosition == -1) {
			continue;
		}
		
		// Igore this base if the sequencing quality is low
		int BaseQuality = int(BaseQualities[SNPRelativePosition]) - 33; // base qualities are ascii - 33
		if (BaseQuality < MIN_BASE_QUALITY) {
			AllGoodQualitySNPs = false;
			continue;
		}
		
		// Get the base on the read
		char ReadBase = ReadSeq[SNPRelativePosition];
		//		if (SNPPosition == 2281) {
		//			cout << "\t" << ReadStart << "\t" << ReadSeq << "\t" << MD << "\t" << SNPPosition << "\t" << SNPRelativePosition << "\t" << ConsensusBase << "\t" << ReadBase << endl;
		//		}
		
		// ignore this base if it's not part of the consensus 
		string strBase = NucleotideMap[ConsensusBase];
		if (strBase.find_first_of(ReadBase) == string::npos) { 
			continue;
		}
		
		// Add the SNP to this read's SNP list
		ReadSNPList.insert (pair <long, char> (SNPPosition, ReadBase)); 
		//cout << ReadStart << "\t" << SNPPosition << "\t" << SNPRelativePosition << "\t" << ConsensusBase << "\t" << ReadBase << endl;
		
	}
	return AllGoodQualitySNPs;
	
}

string CreateSNPGroupId (set <SNP>& SNPGroup) {
	
	stringstream id;
	set< SNP >::iterator itSNPGroup = SNPGroup.begin();
	for (; itSNPGroup != SNPGroup.end(); itSNPGroup++) {
		id << itSNPGroup->first << itSNPGroup->second << "-";
	}
	return id.str();
}

void ResetGroupFlags (vector<SNPGroup>& RefSNPGroups) {
	vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		//		itRefSNPGroups->GroupModified = false;
		itRefSNPGroups->GroupActive = true;
		itRefSNPGroups->OverlapSize = 0;
	}
}

void CreateRefSNPGroups (set< set<SNP> >& ReadSNPGroups, map<string, int> ReadSNPReadCount, vector<SNPGroup>& RefSNPGroups, map<SNP, int>& RefSNPReadCount) { 
/* Add the SNP groups such that largest group is at the begining */
	
	int MaxSize = 0;
	map<int, vector<SNPGroup> > HashRefSNPGroups;
	set< set<SNP> >::iterator itReadSNPGroups = ReadSNPGroups.begin();
	for (; itReadSNPGroups != ReadSNPGroups.end(); itReadSNPGroups++) {
		set<SNP> ReadSNPGroup = *itReadSNPGroups;
		string GroupId = CreateSNPGroupId(ReadSNPGroup);
		int ReadCount = ReadSNPReadCount[GroupId];
		
		// create a new group
		SNPGroup NewSNPGroup;
		set< SNP >::iterator itReadSNPGroup = ReadSNPGroup.begin();
		while (itReadSNPGroup != ReadSNPGroup.end()) {
			NewSNPGroup.SNPList.insert(*itReadSNPGroup);
			NewSNPGroup.SNPReadCount[*itReadSNPGroup] = ReadCount;
			RefSNPReadCount[*itReadSNPGroup] += ReadCount;
			
			itReadSNPGroup++;
		}
		
		// Add to the list of groups based on group size
		HashRefSNPGroups[ReadSNPGroup.size()].push_back(NewSNPGroup);		
	}
	
	// create the final list sorted according to the size (descending order)
	map<int, vector<SNPGroup> >::reverse_iterator rit = HashRefSNPGroups.rbegin();
	for (; rit != HashRefSNPGroups.rend(); ++rit) {
		vector<SNPGroup>::iterator itSNPGroup = rit->second.begin();
		for (; itSNPGroup != rit->second.end(); itSNPGroup++) {
			RefSNPGroups.push_back(*itSNPGroup);
		}
	}
	
}

void RemoveEmbeddedSNPGroups (vector<SNPGroup>& RefSNPGroups) {
	
	ResetGroupFlags(RefSNPGroups);
	
	set<vector< SNPGroup >::iterator> SNPGroupsToBeDeleted;
	vector<SNPGroup>::iterator itRefSNPGroups1 = RefSNPGroups.begin();
	for (; itRefSNPGroups1 != RefSNPGroups.end(); itRefSNPGroups1++) {		
		if (!itRefSNPGroups1->GroupActive) {
			continue;
		}
		
		// get the list SNP position in this group
		long LastSNPPositionGroup1 = itRefSNPGroups1->SNPList.rbegin()->first;

		vector<SNPGroup>::iterator itRefSNPGroups2 = itRefSNPGroups1;//RefSNPGroups.begin();
		itRefSNPGroups2++;
		for (; itRefSNPGroups2 != RefSNPGroups.end(); itRefSNPGroups2++) {

			set< SNP >::iterator itSNPList2 = itRefSNPGroups2->SNPList.begin();	
			if (LastSNPPositionGroup1 < itSNPList2->first) { // if the group 1 ends before group 2 starts, no overlap possible
				continue;
			}
			set< SNP >::iterator itSNPList1 = itRefSNPGroups1->SNPList.begin();
			if (itRefSNPGroups2->SNPList.rbegin()->first < itSNPList1->first) { // if the group 2 ends before group 1 starts, no overlap possible
				continue;
			}

			int RefSNPGroup2Size = itRefSNPGroups2->SNPList.size();

			int OverlappingSNPs = 0;
			while (itSNPList1 != itRefSNPGroups1->SNPList.end() && itSNPList2 != itRefSNPGroups2->SNPList.end()) { 
				if (itSNPList1->first == itSNPList2->first) { // same position
					if (itSNPList1->second == itSNPList2->second) {
						OverlappingSNPs++;
						// move to the next SNP
						itSNPList1++;
						itSNPList2++;
					} else {
						//cout << itSNPList1->first << "\t" << itSNPList1->second << "\t" << itSNPList2->first << "\t" << itSNPList2->second << endl;
						break;
					}
				} else if (itSNPList1->first < itSNPList2->first) {
					itSNPList1++;
				} else  {
					itSNPList2++;
				}
				
			}
			
			if (OverlappingSNPs == RefSNPGroup2Size) { // all SNPs are embedded in group 1
/*				cout << "Embedded groups:" << endl;
				set<SNP>::iterator itSNPGroup1 = itRefSNPGroups1->SNPList.begin();
				cout << "\t-> ";
				for (; itSNPGroup1 != itRefSNPGroups1->SNPList.end(); itSNPGroup1++) {
					cout << itSNPGroup1->first << "," << itSNPGroup1->second << "(" << itRefSNPGroups1->SNPReadCount[*itSNPGroup1] << ")" << " ";
				}
				cout << endl;
				set<SNP>::iterator itSNPGroup2 = itRefSNPGroups2->SNPList.begin();
				cout << "\t-> ";
				for (; itSNPGroup2 != itRefSNPGroups2->SNPList.end(); itSNPGroup2++) {
					cout << itSNPGroup2->first << "," << itSNPGroup2->second << "(" << itRefSNPGroups2->SNPReadCount[*itSNPGroup2] << ")" << " ";
				}
				cout << endl;
*/												
				itRefSNPGroups1->SNPList.insert(itRefSNPGroups2->SNPList.begin(), itRefSNPGroups2->SNPList.end());
				map< SNP, int >::iterator itSNPReadCount = itRefSNPGroups2->SNPReadCount.begin();
				for ( ; itSNPReadCount != itRefSNPGroups2->SNPReadCount.end(); itSNPReadCount++) {
					itRefSNPGroups1->SNPReadCount[itSNPReadCount->first] += itSNPReadCount->second;
				}
//				set< SNP >::iterator itSNPList = itRefSNPGroups2->SNPList.begin();
//				for ( ; itSNPList != itRefSNPGroups2->SNPList.end(); itSNPList++) {
//					itRefSNPGroups1->SNPReadCount[*itSNPList] += itRefSNPGroups2->SNPReadCount[*itSNPList];
//				}
				itRefSNPGroups2->GroupActive = false;
				SNPGroupsToBeDeleted.insert(itRefSNPGroups2);
				
/*				cout << "New Groups:" << endl;
				vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
				for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
					if (!itRefSNPGroups->GroupActive) {
						continue;
					}
					set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
					cout << "\t-> ";
					for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
						cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << " ";
					}
					cout << endl;
				}
				cout << endl;
*/				
			}
		} // end for itRefSNPGroups2 
	} // end while itRefSNPGroups1
	
	for (set<vector< SNPGroup >::iterator>::reverse_iterator it = SNPGroupsToBeDeleted.rbegin(); it != SNPGroupsToBeDeleted.rend(); ++it) {
		RefSNPGroups.erase(*it);
	}
}

void RemoveDuplicateSNPGroups (vector<SNPGroup>& RefSNPGroups) {
	
	set<vector< SNPGroup >::iterator> SNPGroupsToBeDeleted;
	vector<SNPGroup>::iterator itRefSNPGroups1 = RefSNPGroups.begin();
	for (; itRefSNPGroups1 != RefSNPGroups.end(); itRefSNPGroups1++) {		
		// get the list SNP position in this group
		long LastSNPPositionGroup1 = itRefSNPGroups1->SNPList.rbegin()->first;
		int RefSNPGroup1Size = itRefSNPGroups1->SNPList.size();
		
		vector<SNPGroup>::iterator itRefSNPGroups2 = itRefSNPGroups1;//RefSNPGroups.begin();
		itRefSNPGroups2++;
		for (; itRefSNPGroups2 != RefSNPGroups.end(); itRefSNPGroups2++) {
			
			int RefSNPGroup2Size = itRefSNPGroups2->SNPList.size();
			
			if (RefSNPGroup1Size != RefSNPGroup2Size) { // duplicate groups have the same size!
				continue;
			}
			
			set< SNP >::iterator itSNPList2 = itRefSNPGroups2->SNPList.begin();	
			if (LastSNPPositionGroup1 < itSNPList2->first) { // if the group 1 ends before group 2 starts, no overlap possible
				continue;
			}
			set< SNP >::iterator itSNPList1 = itRefSNPGroups1->SNPList.begin();
			if (itRefSNPGroups2->SNPList.rbegin()->first < itSNPList1->first) { // if the group 2 ends before group 1 starts, no overlap possible
				continue;
			}
			
			int OverlappingSNPs = 0;
			while (itSNPList1 != itRefSNPGroups1->SNPList.end() && itSNPList2 != itRefSNPGroups2->SNPList.end()) { 
				if (itSNPList1->first == itSNPList2->first) { // same position
					if (itSNPList1->second == itSNPList2->second) {
						OverlappingSNPs++;
						// move to the next SNP
						itSNPList1++;
						itSNPList2++;
					} else {
						//cout << itSNPList1->first << "\t" << itSNPList1->second << "\t" << itSNPList2->first << "\t" << itSNPList2->second << endl;
						break;
					}
				} else if (itSNPList1->first < itSNPList2->first) {
					itSNPList1++;
				} else  {
					itSNPList2++;
				}
				
			}
			
			if (OverlappingSNPs == RefSNPGroup1Size) { // duplicate group
				map< SNP, int >::iterator itSNPReadCount = itRefSNPGroups2->SNPReadCount.begin();
				for ( ; itSNPReadCount != itRefSNPGroups2->SNPReadCount.end(); itSNPReadCount++) {
					itRefSNPGroups1->SNPReadCount[itSNPReadCount->first] += itSNPReadCount->second;
				}
				SNPGroupsToBeDeleted.insert(itRefSNPGroups2);
				
			}
		} // end for itRefSNPGroups2 
	} // end while itRefSNPGroups1
	
	for (set<vector< SNPGroup >::iterator>::reverse_iterator it = SNPGroupsToBeDeleted.rbegin(); it != SNPGroupsToBeDeleted.rend(); ++it) {
		RefSNPGroups.erase(*it);
	}
}

void RemoveSNPGroupsWithLowReadCount (vector<SNPGroup>& RefSNPGroups) {
	
	vector<vector< SNPGroup >::iterator> SNPGroupsToBeDeleted;
	vector<SNPGroup>::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {		
		map<SNP, int >::iterator itSNPReadCount = itRefSNPGroups->SNPReadCount.begin();
		for (; itSNPReadCount != itRefSNPGroups->SNPReadCount.end(); itSNPReadCount++) {
			if (itSNPReadCount->second < MIN_READ_COUNT) {
				SNPGroupsToBeDeleted.push_back(itRefSNPGroups);
				break;
			}
		}
			
	} // end for itRefSNPGroups
	
	for (vector<vector< SNPGroup >::iterator>::reverse_iterator it = SNPGroupsToBeDeleted.rbegin(); it != SNPGroupsToBeDeleted.rend(); ++it) {
		RefSNPGroups.erase(*it);
	}
}

void RemoveSNPGroupsWithLowReadPercentage (vector<SNPGroup>& RefSNPGroups, map<SNP, int>& RefSNPReadCount) {
	
	vector<vector< SNPGroup >::iterator> SNPGroupsToBeDeleted;
	vector<SNPGroup>::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		
		set< SNP >::iterator itSNPList = itRefSNPGroups->SNPList.begin();
		for (; itSNPList != itRefSNPGroups->SNPList.end(); itSNPList++) {
			int ReadCount = itRefSNPGroups->SNPReadCount[*itSNPList];
			if ((double)ReadCount / (double)RefSNPReadCount[*itSNPList] * 100 < MIN_READ_PERCENTAGE) {
				SNPGroupsToBeDeleted.push_back(itRefSNPGroups);
				break;
			}
		} // end for itSNPList
				
	} // end for itRefSNPGroups
	
	for (vector<vector< SNPGroup >::iterator>::reverse_iterator it = SNPGroupsToBeDeleted.rbegin(); it != SNPGroupsToBeDeleted.rend(); ++it) {
		RefSNPGroups.erase(*it);
	}
}

void RemoveSNPGroupsWithLowReadSupport (vector<SNPGroup>& RefSNPGroups, map<SNP, int>& RefSNPReadCount) {
	
	vector<vector< SNPGroup >::iterator> SNPGroupsToBeDeleted;
	vector<SNPGroup>::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		
		set< SNP >::iterator itSNPList = itRefSNPGroups->SNPList.begin();
		for (; itSNPList != itRefSNPGroups->SNPList.end(); itSNPList++) {
			int ReadCount = itRefSNPGroups->SNPReadCount[*itSNPList];
			if (ReadCount < MIN_READ_COUNT) {
				SNPGroupsToBeDeleted.push_back(itRefSNPGroups);
				break;				
			} else if ((double)ReadCount / (double)RefSNPReadCount[*itSNPList] * 100 < MIN_READ_PERCENTAGE) {
				SNPGroupsToBeDeleted.push_back(itRefSNPGroups);
				break;
			}
		} // end for itSNPList
		
	} // end for itRefSNPGroups
	
	for (vector<vector< SNPGroup >::iterator>::reverse_iterator it = SNPGroupsToBeDeleted.rbegin(); it != SNPGroupsToBeDeleted.rend(); ++it) {
		RefSNPGroups.erase(*it);
	}
}

void SortRefSNPGroupsBySize (vector<SNPGroup>& RefSNPGroups) { 
	/* Sort the SNP groups such that largest group is at the begining */
	
	map<int, vector<SNPGroup> > HashRefSNPGroups;
	vector<SNPGroup>::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		// Add to the list of groups based on group size
		HashRefSNPGroups[itRefSNPGroups->SNPList.size()].push_back(*itRefSNPGroups);		
	}
	
	RefSNPGroups.clear();
	// create the final list sorted according to the size (descending order)
	map<int, vector<SNPGroup> >::reverse_iterator rit = HashRefSNPGroups.rbegin();
	for (; rit != HashRefSNPGroups.rend(); ++rit) {
		vector<SNPGroup>::iterator itSNPGroup = rit->second.begin();
		for (; itSNPGroup != rit->second.end(); itSNPGroup++) {
			RefSNPGroups.push_back(*itSNPGroup);
		}
	}
	
}

void MergeOverlappingSNPGroups (vector<SNPGroup>& RefSNPGroups) {
// iteratively merge the groups with the best overlap.
	
	
	ResetGroupFlags(RefSNPGroups);
	while (true) {
		int MaxOverlapSize = 0;
		vector<SNPGroup>::iterator itSelectedGroup1;
		vector< vector<SNPGroup>::iterator > SelectedGroups2; // multiple overlap is possible, so need to save all overlapping groups
		vector<SNPGroup>::iterator itRefSNPGroups1 = RefSNPGroups.begin();
		for (; itRefSNPGroups1 != RefSNPGroups.end(); itRefSNPGroups1++) {
			
			if (!itRefSNPGroups1->GroupActive) { // this group has been merged in a bigger group once or more before, no need to check this
				continue;
			}
			
			// get the last SNP position in this group
			long LastSNPPositionGroup1 = itRefSNPGroups1->SNPList.rbegin()->first;
			
			vector<SNPGroup>::iterator itRefSNPGroups2 = itRefSNPGroups1;
			itRefSNPGroups2++;
			for (; itRefSNPGroups2 != RefSNPGroups.end(); itRefSNPGroups2++) {
				set< SNP >::iterator itSNPList2 = itRefSNPGroups2->SNPList.begin();	
				if (LastSNPPositionGroup1 < itSNPList2->first) { // if the group 1 ends before group 2 starts, no overlap possible
					continue;
				}
				set< SNP >::iterator itSNPList1 = itRefSNPGroups1->SNPList.begin();
				if (itRefSNPGroups2->SNPList.rbegin()->first < itSNPList1->first) { // if the group 2 ends before group 1 starts, no overlap possible
					continue;
				}
				
				if (itSNPList1->first == itSNPList2->first) { // same position, overlap not possible
					continue;
				}
				
				int OverlappingSNPs = 0;
				bool OverlapFound = false;
				while (itSNPList1 != itRefSNPGroups1->SNPList.end() && itSNPList2 != itRefSNPGroups2->SNPList.end()) { 
					if (itSNPList1->first == itSNPList2->first) { // same position
						if (itSNPList1->second == itSNPList2->second) {
							OverlappingSNPs++;
							OverlapFound = true;
							// move to the next SNP
							itSNPList1++;
							itSNPList2++;
						} else {
							OverlapFound = false;
							break;
						}
					} else if (itSNPList1->first < itSNPList2->first) {
						itSNPList1++;
					} else  {
						itSNPList2++;
					}
					
				}
				if (OverlapFound) { 
					if (OverlappingSNPs == itRefSNPGroups2->SNPList.size()) { // group 2 is embedded in group 1 => they have been merged before 
						continue;
					} else if (!itRefSNPGroups2->GroupActive && OverlappingSNPs < itRefSNPGroups2->OverlapSize) { // this group already had a better overlappe before
						continue;
					}
					
					if (OverlappingSNPs > MaxOverlapSize) { 
						itSelectedGroup1 = itRefSNPGroups1;
						SelectedGroups2.clear();
						SelectedGroups2.push_back(itRefSNPGroups2);
						MaxOverlapSize = OverlappingSNPs;
						
/*						cout << endl << "Overlap of size " << OverlappingSNPs << " found:" << endl;
						set<SNP>::iterator itSNPGroup1 = itRefSNPGroups1->SNPList.begin();
						cout << "1. ";
						for (; itSNPGroup1 != itRefSNPGroups1->SNPList.end(); itSNPGroup1++) {
							cout << itSNPGroup1->first << "," << itSNPGroup1->second << "(" << itRefSNPGroups1->SNPReadCount[*itSNPGroup1] << ")" << " ";
						}
						cout << endl;
						set<SNP>::iterator itSNPGroup2 = itRefSNPGroups2->SNPList.begin();
						cout << "2. ";
						for (; itSNPGroup2 != itRefSNPGroups2->SNPList.end(); itSNPGroup2++) {
							cout << itSNPGroup2->first << "," << itSNPGroup2->second << "(" << itRefSNPGroups2->SNPReadCount[*itSNPGroup2] << ")" << " ";
						}
						cout << endl;
*/						
					} else if (OverlappingSNPs == MaxOverlapSize && itSelectedGroup1 == itRefSNPGroups1) { // multiple overlap for the group 1
						//cout << "Multiple overlap of size " << OverlappingSNPs << " found" << endl;
						SelectedGroups2.push_back(itRefSNPGroups2);

/*						set<SNP>::iterator itSNPGroup2 = itRefSNPGroups2->SNPList.begin();
						cout << "2. ";
						for (; itSNPGroup2 != itRefSNPGroups2->SNPList.end(); itSNPGroup2++) {
							cout << itSNPGroup2->first << "," << itSNPGroup2->second << "(" << itRefSNPGroups2->SNPReadCount[*itSNPGroup2] << ")" << " ";
						}
						cout << endl;
*/						
					}
					
				} 
			} // end for itRefSNPGroups2 
			
		} // end while itRefSNPGroups1
				
		if (MaxOverlapSize > 0) {
/*			cout << "Merging groups with overlap size " << MaxOverlapSize << endl;
//			cout << "Merging the following groups:" << endl;
			set<SNP>::iterator itSNPGroup1 = itSelectedGroup1->SNPList.begin();
			cout << "\t" << itSelectedGroup1->GroupActive << " + ";
			for (; itSNPGroup1 != itSelectedGroup1->SNPList.end(); itSNPGroup1++) {
				//cout << itSNPGroup1->first << "," << itSNPGroup1->second << "(" << itSelectedGroup1->SNPReadCount[*itSNPGroup1] << ")" << " ";
				cout << itSNPGroup1->first << "," << itSNPGroup1->second << " ";
			}
			cout << endl;
*/
			// merge all group 2's with group 1
			SNPGroup SNPGroup1 = *itSelectedGroup1;
			vector< SNPGroup > NewRefSNPGroups; // New SNP groups after duplication 
			vector< vector<SNPGroup>::iterator >::iterator itSelectedGroups2 = SelectedGroups2.begin();
			for (; itSelectedGroups2 != SelectedGroups2.end(); itSelectedGroups2++) {
				SNPGroup SNPGroup2 = **itSelectedGroups2;
				
/*				set<SNP>::iterator itSNPGroup2 = SNPGroup2.SNPList.begin();
				cout << "\t" << SNPGroup2.GroupActive << " + ";
				for (; itSNPGroup2 != SNPGroup2.SNPList.end(); itSNPGroup2++) {
					//cout << itSNPGroup2->first << "," << itSNPGroup2->second << "(" << SNPGroup2.SNPReadCount[*itSNPGroup2] << ")" << " ";
					cout << itSNPGroup2->first << "," << itSNPGroup2->second << " ";
				}
				cout << endl;
*/				
				// mark the group as inactive
				(*itSelectedGroups2)->GroupActive = false;
				// set the overlap size that was found for this group
				(*itSelectedGroups2)->OverlapSize = MaxOverlapSize;

				// duplicate group 1
				SNPGroup NewSNPGroup1 (SNPGroup1);
				// merge this copy of group 1 with current group 2
				NewSNPGroup1.SNPList.insert(SNPGroup2.SNPList.begin(), SNPGroup2.SNPList.end());
				map< SNP, int >::iterator itSNPReadCount = SNPGroup2.SNPReadCount.begin();
				for ( ; itSNPReadCount != SNPGroup2.SNPReadCount.end(); itSNPReadCount++) {
					NewSNPGroup1.SNPReadCount[itSNPReadCount->first] += itSNPReadCount->second;
				}
//				set< SNP >::iterator itSNPList = SNPGroup2.SNPList.begin();
//				for ( ; itSNPList != SNPGroup2.SNPList.end(); itSNPList++) {
//					NewSNPGroup1.SNPReadCount[*itSNPList] += SNPGroup2.SNPReadCount[*itSNPList];
//				}
				// add to the list
				NewRefSNPGroups.push_back(NewSNPGroup1);
			} // end for each selected group 2

			// erase all group 2
//AMT			vector< vector<SNPGroup>::iterator >::reverse_iterator ritSelectedGroups2 = SelectedGroups2.rbegin();
//AMT			for (; ritSelectedGroups2 != SelectedGroups2.rend(); ++ritSelectedGroups2) {
//AMT				RefSNPGroups.erase(*ritSelectedGroups2);
//AMT				(*ritSelectedGroups2)->GroupActive = false;
//AMT			}

			// erase group 1
			RefSNPGroups.erase(itSelectedGroup1);
			// add newly created groups
			vector<SNPGroup>::iterator itNewRefSNPGroups = NewRefSNPGroups.begin();
			for (; itNewRefSNPGroups != NewRefSNPGroups.end(); itNewRefSNPGroups++) {
				RefSNPGroups.push_back(*itNewRefSNPGroups);
			}
			
//			RemoveEmbeddedSNPGroups(RefSNPGroups);
			RemoveDuplicateSNPGroups(RefSNPGroups);

			SortRefSNPGroupsBySize(RefSNPGroups);
			
/*			cout << "New Groups:" << endl;
			vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
			for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
				set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
				cout << "\t" << itRefSNPGroups->GroupActive << " + ";
				for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
					//cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << " ";
					cout << itSNPGroup->first << "," << itSNPGroup->second << " ";
				}
				cout << endl;
			}
			cout << endl;
*/			
		} else {
			break;
		}
	} // end while true

	// remove the group which have been marked as inactive
	vector<vector< SNPGroup >::iterator> SNPGroupsToBeDeleted;
	vector<SNPGroup>::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		if (!itRefSNPGroups->GroupActive) { // this group has been merged in a bigger group once or more before, delete this 
			SNPGroupsToBeDeleted.push_back(itRefSNPGroups);
		}
	}
	for (vector<vector< SNPGroup >::iterator>::reverse_iterator it = SNPGroupsToBeDeleted.rbegin(); it != SNPGroupsToBeDeleted.rend(); ++it) {
		RefSNPGroups.erase(*it);
	}
	
}

bool CheckDiscrepancyInReadPair (set< SNP >& Read1SNPList, set< SNP >& Read2SNPList) {
	
	set< SNP >::iterator itRead1SNPList = Read1SNPList.begin();
	set< SNP >::iterator itRead2SNPList = Read2SNPList.begin();
	while (itRead1SNPList != Read1SNPList.end() && itRead2SNPList != Read2SNPList.end()) {
		
		if (itRead1SNPList->first < itRead2SNPList->first) {
			itRead1SNPList++;
		} else if (itRead1SNPList->first == itRead2SNPList->first) { // same position
			if (itRead1SNPList->second == itRead2SNPList->second) { //same SNP
				itRead1SNPList++;
				itRead2SNPList++;
			} else {
				return true;
			}
		} else if (itRead1SNPList->first > itRead2SNPList->first) {
			itRead2SNPList++;
		} 
	} // end while	
	
	return false;
}

void FinaliseSNPGroups(vector<SNPGroup>& RefSNPGroups, map < SNP, int >& RefSNPReadCount, set< set <SNP> >& ReadSNPGroups, map<string, int>& ReadSNPReadCount) {
	
	if (ReadSNPGroups.size() == 0) {
		return;
	}
	
/*	cout << ReadSNPGroups.size() << endl;
	set< set <SNP> >::iterator itReadSNPGroups = ReadSNPGroups.begin();
	for (; itReadSNPGroups != ReadSNPGroups.end(); itReadSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itReadSNPGroups->begin();
		set<SNP> ReadSNPGroup = *itReadSNPGroups;
		string GroupId = CreateSNPGroupId(ReadSNPGroup);
		cout << ReadSNPReadCount[GroupId] << ": ";
		for (; itSNPGroup != itReadSNPGroups->end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "\t";
		}
		cout << endl;
	}
*/	
	//cout << "\tCreating SNP groups" << endl;
	// create SNP groups from these read SNPs
	CreateRefSNPGroups(ReadSNPGroups, ReadSNPReadCount, RefSNPGroups, RefSNPReadCount);
	
/*	vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
		for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << "\t";
		}
		cout << endl;
	}
*/	
	// Remove the embedded Groups
	//cout << "\tRemoving embedded SNP groups" << endl;
	RemoveEmbeddedSNPGroups(RefSNPGroups);
	
/*	itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
		for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << "\t";
		}
		cout << endl;
	}
*/	
	// merge overlapping groups
	//cout << "\tMerging SNP groups" << endl;
	MergeOverlappingSNPGroups(RefSNPGroups);
	
/*	itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
		for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << "\t";
		}
		cout << endl;
	}
*/
	// Re-Remove the embedded Groups (merging might have resulted in embedded groups)
	//cout << "\tRe-Removing embedded SNP groups" << endl;
	RemoveEmbeddedSNPGroups(RefSNPGroups);
	
/*	vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
		for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << "\t";
		}
		cout << endl;
	}
*/	
//	// Remove the groups for which read count falls below the threshold
//	cout << "Removing SNP groups with low read count" << endl;
//	RemoveSNPGroupsWithLowReadCount(RefSNPGroups);	

	// Remove the groups for which read count falls below the threshold or where the percentage of reads supporting a SNP is less than the threshold
//	cout << "\tRemoving SNP groups with low read support" << endl;
	RemoveSNPGroupsWithLowReadSupport(RefSNPGroups, RefSNPReadCount);	
	
/*	itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
		for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << "\t";
		}
		cout << endl;
	}
/*		
	// Remove unlikely groups (groups where the percentage of reads supporting a SNP is less than the threshold) 
	cout << "Removing SNP groups with low read percentage" << endl;
	RemoveSNPGroupsWithLowReadPercentage(RefSNPGroups, RefSNPReadCount);

	itRefSNPGroups = RefSNPGroups.begin();
	for (; itRefSNPGroups != RefSNPGroups.end(); itRefSNPGroups++) {
		set<SNP>::iterator itSNPGroup = itRefSNPGroups->SNPList.begin();
		for (; itSNPGroup != itRefSNPGroups->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "(" << itRefSNPGroups->SNPReadCount[*itSNPGroup] << ")" << "\t";
		}
		cout << endl;
	}
*/	

/*	// Print the groups
	//vector< SNPGroup >::iterator 
	vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
	while (itRefSNPGroups != RefSNPGroups.end()) {
		set< SNP >::iterator itSNPs = itRefSNPGroups->SNPList.begin();
		cout << "Group Size: " << itRefSNPGroups->SNPList.size() << endl;
		while (itSNPs != itRefSNPGroups->SNPList.end()) {
			cout << "\t" << (*itSNPs).first << "\t" << (*itSNPs).second << "\t" << itRefSNPGroups->SNPReadCount[*itSNPs] << endl;
			//cout << "\t" << (*itSNPs).first << "\t" << (*itSNPs).second << "\t" << RefSNPReadCount[*itSNPs] << endl;
			
			itSNPs++;
		}
		
		itRefSNPGroups++;
	}
*/	
}

bool CheckPositionOverlap(vector<SNPGroup>::iterator &itRefSNPGroup, vector<SNPGroup>::iterator &itSNPGroup) {
	
	set< SNP >::iterator itSNPList2 = itSNPGroup->SNPList.begin();	
	if (itRefSNPGroup->SNPList.rbegin()->first < itSNPList2->first) { // if the group 1 ends before group 2 starts, no overlap possible
		return false;
	}
	set< SNP >::iterator itSNPList1 = itRefSNPGroup->SNPList.begin();
	if (itSNPGroup->SNPList.rbegin()->first < itSNPList1->first) { // if the group 2 ends before group 1 starts, no overlap possible
		return false;
	}
	
	while (itSNPList1 != itRefSNPGroup->SNPList.end() && itSNPList2 != itSNPGroup->SNPList.end()) { 
		if (itSNPList1->first == itSNPList2->first) { // same position
			return true;
		} else if (itSNPList1->first < itSNPList2->first) {
			itSNPList1++;
		} else  {
			itSNPList2++;
		}
	} // end while

	return false;
}

void CheckGenomeMembershipA (vector<SNPGroup>::iterator &itRefSNPGroup, vector <vector< SNPGroup >::iterator> &ASNPGroups, int ASNPCount, int ANoCoverage) {
	
	double GroupSize = itRefSNPGroup->SNPList.size();
	double ASNPPercentage = (double)ASNPCount/(double)(GroupSize - ANoCoverage);
	if (ASNPPercentage >= GENOME_MEMBERSHIP_THRESHOLD) {

		vector < vector <vector< SNPGroup >::iterator>::iterator > ASNPGroupsToBeDeleted;

		bool GroupToBeAdded = false;
		bool OverlappingGroupFound = false;
		vector< vector< SNPGroup >::iterator >::iterator itASNPGroup = ASNPGroups.begin();
		for (; itASNPGroup != ASNPGroups.end(); itASNPGroup++) {
			if (CheckPositionOverlap(itRefSNPGroup, *itASNPGroup)) {
				OverlappingGroupFound = true;
				if (GroupSize >= (*itASNPGroup)->SNPList.size() ) { // ideally, GroupSize will never be greater than SNPlist size since groups are sorted in descending order by size
					double CurrentASNPPercentage = (*itASNPGroup)->ASNPPercentage;
					if (ASNPPercentage > CurrentASNPPercentage) { // this is a better group than previously found
						// mark the previous group as false
						(*itASNPGroup)->AGenome = false;
						//ASNPGroups.erase(itASNPGroup);
						
						ASNPGroupsToBeDeleted.push_back(itASNPGroup);
						
						// add this group to the SNPGroups list
						itRefSNPGroup->AGenome = true;
						itRefSNPGroup->ASNPPercentage = ASNPPercentage;
						
						GroupToBeAdded = true;
					} else if (ASNPPercentage == CurrentASNPPercentage) {
						// add this group to the SNPGroups list
						itRefSNPGroup->AGenome = true;
						itRefSNPGroup->ASNPPercentage = ASNPPercentage;

						GroupToBeAdded = true;
					}
				} else {
					double CurrentASNPPercentage = (*itASNPGroup)->ASNPPercentage;
					if (ASNPPercentage >= CurrentASNPPercentage) { 
						// add this group to the SNPGroups list
						itRefSNPGroup->AGenome = true;
						itRefSNPGroup->ASNPPercentage = ASNPPercentage;
						
						GroupToBeAdded = true;
					}
				} // if group size is equal to current ASNPGroup.size
				//break;
			} // if position overlap found
		} // for each ASNPGroup
	
		if (!OverlappingGroupFound) {
			//cout << "\tA-Genome" << endl;
			itRefSNPGroup->AGenome = true;
			itRefSNPGroup->ASNPPercentage = ASNPPercentage;
			ASNPGroups.push_back(itRefSNPGroup);
		} else {
			vector < vector <vector< SNPGroup >::iterator>::iterator >::reverse_iterator ritASNPGroupsToBeDeleted = ASNPGroupsToBeDeleted.rbegin();
			for (; ritASNPGroupsToBeDeleted != ASNPGroupsToBeDeleted.rend(); ++ritASNPGroupsToBeDeleted) {
				ASNPGroups.erase(*ritASNPGroupsToBeDeleted);
			}
			if (GroupToBeAdded) {
				ASNPGroups.push_back(itRefSNPGroup);
			}
			
		}

		
	} // if ASNPPercentage >= threshold
}

void CheckGenomeMembershipB (vector<SNPGroup>::iterator &itRefSNPGroup, vector <vector< SNPGroup >::iterator> &BSNPGroups, int BSNPCount, int BNoCoverage) {
	
	double GroupSize = itRefSNPGroup->SNPList.size();
	double BSNPPercentage = (double)BSNPCount/(double)(GroupSize - BNoCoverage);
	if (BSNPPercentage >= GENOME_MEMBERSHIP_THRESHOLD) {
		
		vector < vector <vector< SNPGroup >::iterator>::iterator > BSNPGroupsToBeDeleted;
		
		bool GroupToBeAdded = false;
		bool OverlappingGroupFound = false;
		vector< vector< SNPGroup >::iterator >::iterator itBSNPGroup = BSNPGroups.begin();
		for (; itBSNPGroup != BSNPGroups.end(); itBSNPGroup++) {
			if (CheckPositionOverlap(itRefSNPGroup, *itBSNPGroup)) {
				OverlappingGroupFound = true;
				if (GroupSize >= (*itBSNPGroup)->SNPList.size() ) { // ideally, GroupSize will never be greater than SNPlist size since groups are sorted in descending order by size
					double CurrentBSNPPercentage = (*itBSNPGroup)->BSNPPercentage;
					if (BSNPPercentage > CurrentBSNPPercentage) { // this is a better group than previously found
						// mark the previous group as false
						(*itBSNPGroup)->BGenome = false;
						//ASNPGroups.erase(itBSNPGroup);
						
						BSNPGroupsToBeDeleted.push_back(itBSNPGroup);
						
						// add this group to the SNPGroups list
						itRefSNPGroup->BGenome = true;
						itRefSNPGroup->BSNPPercentage = BSNPPercentage;
						
						GroupToBeAdded = true;
					} else if (BSNPPercentage == CurrentBSNPPercentage) {
						// add this group to the SNPGroups list
						itRefSNPGroup->BGenome = true;
						itRefSNPGroup->BSNPPercentage = BSNPPercentage;
						
						GroupToBeAdded = true;
					}
				} else {
					double CurrentBSNPPercentage = (*itBSNPGroup)->BSNPPercentage;
					if (BSNPPercentage >= CurrentBSNPPercentage) { 
						// add this group to the SNPGroups list
						itRefSNPGroup->BGenome = true;
						itRefSNPGroup->BSNPPercentage = BSNPPercentage;
						
						GroupToBeAdded = true;
					}
				} // if group size is equal to current BSNPGroup.size
				//break;
			} // if position overlap found
		} // for each BSNPGroup
		
		if (!OverlappingGroupFound) {
			//cout << "\tB-Genome" << endl;
			itRefSNPGroup->BGenome = true;
			itRefSNPGroup->BSNPPercentage = BSNPPercentage;
			BSNPGroups.push_back(itRefSNPGroup);
		} else {
			vector < vector <vector< SNPGroup >::iterator>::iterator >::reverse_iterator ritBSNPGroupsToBeDeleted = BSNPGroupsToBeDeleted.rbegin();
			for (; ritBSNPGroupsToBeDeleted != BSNPGroupsToBeDeleted.rend(); ++ritBSNPGroupsToBeDeleted) {
				BSNPGroups.erase(*ritBSNPGroupsToBeDeleted);
			}
			if (GroupToBeAdded) {
				BSNPGroups.push_back(itRefSNPGroup);
			}
			
		}
		
		
	} // if BSNPPercentage >= threshold
}

void CheckGenomeMembershipD (vector<SNPGroup>::iterator &itRefSNPGroup, vector <vector< SNPGroup >::iterator> &DSNPGroups, int DSNPCount, int DNoCoverage) {

	double GroupSize = itRefSNPGroup->SNPList.size();
	double DSNPPercentage = (double)DSNPCount/(double)(GroupSize - DNoCoverage);
	if (DSNPPercentage >= GENOME_MEMBERSHIP_THRESHOLD) {
		vector < vector <vector< SNPGroup >::iterator>::iterator > DSNPGroupsToBeDeleted;
		
		bool GroupToBeAdded = false;
		bool OverlappingGroupFound = false;
		vector< vector< SNPGroup >::iterator >::iterator itDSNPGroup = DSNPGroups.begin();
		for (; itDSNPGroup != DSNPGroups.end(); itDSNPGroup++) {
			
			if (CheckPositionOverlap(itRefSNPGroup, *itDSNPGroup)) {
				OverlappingGroupFound = true;
				if (GroupSize >= (*itDSNPGroup)->SNPList.size() ) { // ideally, GroupSize will never be greater than SNPlist size since groups are sorted in descending order by size
					double CurrentDSNPPercentage = (*itDSNPGroup)->DSNPPercentage;
					if (DSNPPercentage > CurrentDSNPPercentage) { // this is a better group than previously found
						// remove the previous group
						(*itDSNPGroup)->DGenome = false;
						//DSNPGroups.erase(itDSNPGroup);
						
						DSNPGroupsToBeDeleted.push_back(itDSNPGroup);
						
						// add this group to the SNPGroups list
						itRefSNPGroup->DGenome = true;
						itRefSNPGroup->DSNPPercentage = DSNPPercentage;
						
						GroupToBeAdded = true;

					} else if (DSNPPercentage == CurrentDSNPPercentage) {
						// add this group to the SNPGroups list
						itRefSNPGroup->DGenome = true;
						itRefSNPGroup->DSNPPercentage = DSNPPercentage;
						
						GroupToBeAdded = true;
					}
				} else {
					double CurrentDSNPPercentage = (*itDSNPGroup)->DSNPPercentage;
					if (DSNPPercentage >= CurrentDSNPPercentage) { 
						// add this group to the SNPGroups list
						itRefSNPGroup->DGenome = true;
						itRefSNPGroup->DSNPPercentage = DSNPPercentage;
						
						GroupToBeAdded = true;
					}
				} // if group size is equal to current DSNPGroup.size
				//break;
			} // if position overlap found
		} // for each DSNPGroup

		if (!OverlappingGroupFound) {
			//cout << "\tD-Genome" << endl;
			itRefSNPGroup->DGenome = true;
			itRefSNPGroup->DSNPPercentage = DSNPPercentage;
			DSNPGroups.push_back(itRefSNPGroup);
		} else {
			vector < vector <vector< SNPGroup >::iterator>::iterator >::reverse_iterator ritDSNPGroupsToBeDeleted = DSNPGroupsToBeDeleted.rbegin();
			for (; ritDSNPGroupsToBeDeleted != DSNPGroupsToBeDeleted.rend(); ++ritDSNPGroupsToBeDeleted) {
				DSNPGroups.erase(*ritDSNPGroupsToBeDeleted);
			}
			if (GroupToBeAdded) {
				DSNPGroups.push_back(itRefSNPGroup);
			}
			
		}
	} // if DSNPPercentage >= threshold
	
}

int WriteSNPsToFile(string reference, map <long, char>& SNPList, char* RefBaseList, char* DiploidSNPList, ofstream &ofs) {
	int SNPCount = 0;
	map<long, char>::iterator itSNP = SNPList.begin();
	for (; itSNP != SNPList.end(); itSNP++) {
		// get the SNP
		long SNPPosition = itSNP->first;
		char SNPBase = itSNP->second;
		
		// get the reference base
		char RefBase = *(RefBaseList + SNPPosition - 1);
		
		// get the diploid base
		char DiploidBase = *(DiploidSNPList + SNPPosition - 1);
		if (DiploidBase == ' ') {
			DiploidBase = RefBase;
		}
		
		if (SNPBase != RefBase || DiploidBase != RefBase) {
			SNPCount++;
		}
		ofs << reference << "\t" << SNPPosition << "\t" << RefBase << "\t" << DiploidBase << "\t" << SNPBase << endl;
	} // end for each SNP
	
	return SNPCount;
}

/*
void WriteSNPsToFile(string reference, set <SNP> SNPList, char* RefBaseList, char* DiploidSNPList, ofstream &ofs) {
	set<SNP>::iterator itSNP = SNPList.begin();
	for (; itSNP != SNPList.end(); itSNP++) {
		// get the SNP
		long SNPPosition = itSNP->first;
		char SNPBase = itSNP->second;
		// get the reference base
		char RefBase = *(RefBaseList + SNPPosition - 1);
		// get the diploid base
		char DiploidBase = *(DiploidSNPList + SNPPosition - 1);
		if (DiploidBase == ' ') {
			DiploidBase = RefBase;
		}
		
		ofs << reference << "\t" << SNPPosition << "\t" << RefBase << "\t" << DiploidBase << "\t" << SNPBase << endl;
	} // end for each SNP
}
*/

void AssignSNPsToGenomes(string reference, vector<SNPGroup>& RefSNPGroups, char* RefRefBaseList, char* RefDiploidASNPList, char* RefDiploidBSNPList, char* RefDiploidDSNPList, int* RefDiploidACoverageList, int* RefDiploidBCoverageList, int* RefDiploidDCoverageList, ofstream &ofs_A, ofstream &ofs_B, ofstream &ofs_D) {
	
	//if (RefSNPGroups.size() == 0) {
	//	return;
	//}
	//cout << "\tAssigning SNP to genomes" << endl;
	
	set< long > DistinctPositions; // distinct SNP positions present in all groups (does not include positions ignored because of low base quality)
	set< long > DistinctPositionsConsidered; // distinct SNP positions considered after removing low coverage positions
	vector <vector< SNPGroup >::iterator> ASNPGroups, BSNPGroups, DSNPGroups;
	vector< SNPGroup >::iterator itRefSNPGroup = RefSNPGroups.begin();
	for (; itRefSNPGroup != RefSNPGroups.end(); itRefSNPGroup++) {
		
		// Initialise necessary flags/variables
		itRefSNPGroup->AGenome = false;
		itRefSNPGroup->BGenome = false;
		itRefSNPGroup->DGenome = false;
		itRefSNPGroup->ASNPPercentage = 0.0;
		itRefSNPGroup->BSNPPercentage = 0.0;
		itRefSNPGroup->DSNPPercentage = 0.0;

		int ASNPCount = 0;
		int BSNPCount = 0;
		int DSNPCount = 0;
		int LowCoverage = 0;
		int ALowCoverage = 0;
		int BLowCoverage = 0;
		int DLowCoverage = 0;
		set< SNP >::iterator itSNP = itRefSNPGroup->SNPList.begin();
		for (; itSNP != itRefSNPGroup->SNPList.end(); itSNP++) {
			long SNPPosition = itSNP->first;
			char SNPBase = itSNP->second;
			
			char RefBase = *(RefRefBaseList + SNPPosition - 1);
			char ABase = *(RefDiploidASNPList + SNPPosition - 1);
			char BBase = *(RefDiploidBSNPList + SNPPosition - 1);
			char DBase;
			if (RefDiploidDSNPList != NULL) {
				DBase = *(RefDiploidDSNPList + SNPPosition - 1);
			}
			int ACoverage = *(RefDiploidACoverageList + SNPPosition - 1);
			int BCoverage = *(RefDiploidBCoverageList + SNPPosition - 1);
			int DCoverage = (RefDiploidDCoverageList == NULL ? numeric_limits<int>::max() : *(RefDiploidDCoverageList + SNPPosition - 1));

			//cout << "\t" << SNPPosition << "\t" << SNPBase << "\t" << ABase << "/" << ACoverage << "\t" << BBase << "/" << BCoverage << "\t" << DBase << "/" << DCoverage;
			DistinctPositions.insert(SNPPosition);

			// Added on 22-03-2011. Ignore positions where one or more diploids have low coverage
			if (ACoverage < MIN_DIPLOID_COVERAGE || BCoverage < MIN_DIPLOID_COVERAGE || DCoverage < MIN_DIPLOID_COVERAGE) { // if there is low coverage in one or more diploids, ignore this position
				LowCoverage++;
				
				//cout << "\tLow Coverage";
				
				continue;
			}
			
			//cout << endl;
			
			DistinctPositionsConsidered.insert(SNPPosition);
			
			//cout << SNPPosition << "\t" << RefBase << "\t" << SNPBase << "\t" << ABase << "\t" << DBase << endl;
			if (SNPBase == ABase) {
				ASNPCount++;
			} else if (ABase == ' ') {
				if (ACoverage < MIN_DIPLOID_COVERAGE) { // if there is low coverage in A diploid, ignore this position
					ALowCoverage++;
				} else if (SNPBase == RefBase) { // no snp in A genome => ABase is same as the reference base.
					ASNPCount++;
				}
			}
			
			if (SNPBase == BBase) {
				BSNPCount++;
			} else if (BBase == ' ') {
				if (BCoverage < MIN_DIPLOID_COVERAGE) { // if there is low coverage in B diploid, ignore this position
					BLowCoverage++;
				} else if (SNPBase == RefBase) { // no snp in B genome => BBase is same as the reference base.
					BSNPCount++;
				}
			}
			
			if (RefDiploidDSNPList != NULL) {
				if (SNPBase == DBase) {
					DSNPCount++;
				} else if (DBase == ' ') {
					if (DCoverage < MIN_DIPLOID_COVERAGE) { // if there is low coverage in D diploid, ignore this position
						DLowCoverage++;
					} else if (SNPBase == RefBase) { // no snp in D genome => DBase is same as the reference base.
						DSNPCount++;
					}
				}
			}				
			
		} // end for each SNP in the group
/*
		double GroupSize = itRefSNPGroup->SNPList.size();
		
		vector< SNPGroup >::iterator itRefSNPGroups = RefSNPGroups.begin();
		set<SNP>::iterator itSNPGroup = itRefSNPGroup->SNPList.begin();
		cout << "(" << GroupSize << "): ";
		for (; itSNPGroup != itRefSNPGroup->SNPList.end(); itSNPGroup++) {
			cout << itSNPGroup->first << "," << itSNPGroup->second << "\t";
		}
		cout << endl;

		cout << (double)ASNPCount/GroupSize << " / " << (double)DSNPCount/GroupSize << endl;
*/
		// Check for A Genome Membership and add the group to the list of A Genome SNP Groups if this is a valid group
		CheckGenomeMembershipA(itRefSNPGroup, ASNPGroups, ASNPCount, LowCoverage);
		// Check for B Genome Membership and add the group to the list of B Genome SNP Groups if this is a valid group
		CheckGenomeMembershipB(itRefSNPGroup, BSNPGroups, BSNPCount, LowCoverage);
		if (RefDiploidDSNPList != NULL) {
			// Check for D Genome Membership and add the group to the list of D Genome SNP Groups if this is a valid group
			CheckGenomeMembershipD(itRefSNPGroup, DSNPGroups, DSNPCount, LowCoverage);
		}
		
		//cout << endl;
	
	}// end for each Ref SNP group

	// Now create a new list of SNPs by removing those groups which are invalid (A/B/D Genome flag is marked false in A/B/D SNPGroup lists)
	// Use set to sort the according to SNP positions 
	//set< SNP > ASNPList, DSNPList;
	map< long, char > ASNPList, BSNPList, DSNPList;
	
	// A-Genome
	int ALowCoverage = 0;
	vector< vector< SNPGroup >::iterator >::iterator itASNPGroup = ASNPGroups.begin();
	for (; itASNPGroup != ASNPGroups.end(); itASNPGroup++) {
		if ((*itASNPGroup)->AGenome) {
			//ASNPList.insert((*itASNPGroup)->SNPList.begin(), (*itASNPGroup)->SNPList.end());
			set<SNP>::iterator itSNP = (*itASNPGroup)->SNPList.begin();
			for (; itSNP != (*itASNPGroup)->SNPList.end(); itSNP++) {
				// get the SNP
				long SNPPosition = itSNP->first;
				char SNPBase = itSNP->second;

				int ACoverage = *(RefDiploidACoverageList + SNPPosition - 1);
				int BCoverage = *(RefDiploidBCoverageList + SNPPosition - 1);
				int DCoverage = (RefDiploidDCoverageList == NULL ? numeric_limits<int>::max() : *(RefDiploidDCoverageList + SNPPosition - 1));
				// Added on 22-03-2011. Ignore positions where one or more diploids have low coverage
				if (ACoverage < MIN_DIPLOID_COVERAGE || BCoverage < MIN_DIPLOID_COVERAGE || DCoverage < MIN_DIPLOID_COVERAGE) {
					if (ACoverage < MIN_DIPLOID_COVERAGE) { // Don't keep the SNP if there is no diploid coverage
						ALowCoverage++;
					}
					continue;
				}
				
				map <long, char>::iterator itBSNP = ASNPList.find(SNPPosition);
				if (itBSNP == ASNPList.end()) { // SNP position not seen before
					ASNPList[SNPPosition] = SNPBase;
				} else if (itBSNP->second != SNPBase) { // ambiguous base
					ASNPList.erase(itBSNP);
				}
				
			} // end for each SNP
 
		}
	} // end for each ASNPGroup

	// B-Genome
	int BLowCoverage = 0;
	vector< vector< SNPGroup >::iterator >::iterator itBSNPGroup = BSNPGroups.begin();
	for (; itBSNPGroup != BSNPGroups.end(); itBSNPGroup++) {
		if ((*itBSNPGroup)->BGenome) {
			//BSNPList.insert((*itBSNPGroup)->SNPList.begin(), (*itBSNPGroup)->SNPList.end());
			set<SNP>::iterator itSNP = (*itBSNPGroup)->SNPList.begin();
			for (; itSNP != (*itBSNPGroup)->SNPList.end(); itSNP++) {
				// get the SNP
				long SNPPosition = itSNP->first;
				char SNPBase = itSNP->second;
				
				int ACoverage = *(RefDiploidACoverageList + SNPPosition - 1);
				int BCoverage = *(RefDiploidBCoverageList + SNPPosition - 1);
				int DCoverage = (RefDiploidDCoverageList == NULL ? numeric_limits<int>::max() : *(RefDiploidDCoverageList + SNPPosition - 1));
				// Added on 22-03-2011. Ignore positions where one or more diploids have low coverage
				if (ACoverage < MIN_DIPLOID_COVERAGE || BCoverage < MIN_DIPLOID_COVERAGE || DCoverage < MIN_DIPLOID_COVERAGE) {
					if (BCoverage < MIN_DIPLOID_COVERAGE) { // Don't keep the SNP if there is no diploid coverage
						BLowCoverage++;
					}
					continue;
				}
				
				map <long, char>::iterator itBSNP = BSNPList.find(SNPPosition);
				if (itBSNP == BSNPList.end()) { // SNP position not seen before
					BSNPList[SNPPosition] = SNPBase;
				} else if (itBSNP->second != SNPBase) { // ambiguous base
					BSNPList.erase(itBSNP);
				}
				
			} // end for each SNP
			
		}
	} // end for each BSNPGroup
	
	// D-Genome
	int DLowCoverage = 0;
	if (RefDiploidDSNPList != NULL) {
		vector< vector< SNPGroup >::iterator >::iterator itDSNPGroup = DSNPGroups.begin();
		for (; itDSNPGroup != DSNPGroups.end(); itDSNPGroup++) {
			if ((*itDSNPGroup)->DGenome) {
				//DSNPList.insert((*itDSNPGroup)->SNPList.begin(), (*itDSNPGroup)->SNPList.end());
				set<SNP>::iterator itSNP = (*itDSNPGroup)->SNPList.begin();
				for (; itSNP != (*itDSNPGroup)->SNPList.end(); itSNP++) {
					// get the SNP
					long SNPPosition = itSNP->first;
					char SNPBase = itSNP->second;
					
					int ACoverage = *(RefDiploidACoverageList + SNPPosition - 1);
					int BCoverage = *(RefDiploidBCoverageList + SNPPosition - 1);
					int DCoverage = *(RefDiploidDCoverageList + SNPPosition - 1);
					// Added on 22-03-2011. Ignore positions where one or more diploids have low coverage
					if (ACoverage < MIN_DIPLOID_COVERAGE || BCoverage < MIN_DIPLOID_COVERAGE || DCoverage < MIN_DIPLOID_COVERAGE) {
						if (DCoverage < MIN_DIPLOID_COVERAGE) { // Don't keep the SNP if there is no diploid coverage
							DLowCoverage++;
						}
						continue;
					}				
					
					map <long, char>::iterator itDSNP = DSNPList.find(SNPPosition);
					if (itDSNP == DSNPList.end()) { // SNP position not seen before
						DSNPList[SNPPosition] = SNPBase;
					} else if (itDSNP->second != SNPBase) { // ambiguous base
						DSNPList.erase(itDSNP);
					}
				} // end for each SNP
				
			}
			
		} // end for each DSNPGroup
		
	}
	
	// Write the SNPs (returns the number of SNPs in each genome (either in the Diploid or CS-Genome))
	int ASNPCount = WriteSNPsToFile(reference, ASNPList, RefRefBaseList, RefDiploidASNPList, ofs_A);
	int BSNPCount = WriteSNPsToFile(reference, BSNPList, RefRefBaseList, RefDiploidBSNPList, ofs_B);
	int DSNPCount;
	if (RefDiploidDSNPList != NULL) {
		DSNPCount = WriteSNPsToFile(reference, DSNPList, RefRefBaseList, RefDiploidDSNPList, ofs_D);
	}
	// write the details to the console
	cout << "\t"  << DistinctPositions.size() << "\t"  << DistinctPositionsConsidered.size() << "\t" << ASNPList.size() << "\t" << BSNPList.size();
	if (RefDiploidDSNPList != NULL) {
		cout << "\t" << DSNPList.size();
	}
	cout << "\t" << ASNPCount << "\t" << BSNPCount;
	if (RefDiploidDSNPList != NULL) {
		cout << "\t" << DSNPCount;
	}
	cout << "\t" << ALowCoverage << "\t" << BLowCoverage;
	if (RefDiploidDSNPList != NULL) {
		cout << "\t" << DLowCoverage;
	}
	cout << endl;
	
}

void ReadGFFFile(string GFF_File, map < string, vector< pair<long, long> > >& GFFData) {
	// gff file
	ifstream ifs_gff (GFF_File.c_str());
	string strGFF = "";
	string reference = "";
	string previous_reference = "";
	vector< pair< long, long > > RefGFFData;
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
		
		reference = entries[0];
		long start = atol(entries[3].c_str());
		long end = atol(entries[4].c_str());
		
		if (previous_reference.compare(reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFData;
			}
			
			RefGFFData = GFFData[reference];
			previous_reference = reference;
		}
		
		RefGFFData.push_back(pair<long, long>(start, end));
	} // end while ifs_gff is good
	GFFData[reference] = RefGFFData;
	
}

void ProcessSAMFile(string SAM_File, map < string, vector< pair<long, long> > >& GFFData, SNP_LIST& SNPList, map < string, char* >& RefBaseList, SNP_LIST& DiploidASNPList, SNP_LIST& DiploidBSNPList, SNP_LIST& DiploidDSNPList, COVERAGE_LIST& DiploidACoverageList, COVERAGE_LIST& DiploidBCoverageList, COVERAGE_LIST& DiploidDCoverageList, string OUT_File_A, string OUT_File_B, string OUT_File_D) {
	
	// open the sam file stream
	ifstream ifs_sam ( SAM_File.c_str() );
	
	// open the streams for output
	ofstream ofs_A (OUT_File_A.c_str());
	ofstream ofs_B (OUT_File_B.c_str());
	ofstream ofs_D (OUT_File_D.c_str());
	
	ofs_A << "Reference\tPosition\tReference Base\tA-Diploid Base\tA-CS-Genome Base" << endl;
	ofs_B << "Reference\tPosition\tReference Base\tB-Diploid Base\tB-CS-Genome Base" << endl;
	if (OUT_File_D.length() > 0) {
		ofs_D << "Reference\tPosition\tReference Base\tD-Diploid Base\tD-CS-Genome Base" << endl;
	}
	
	cout << "Start\tEnd\tTotal SNP Positions\tTotal SNP Positions Grouped\tTotal SNP Positions Assigned";
	cout << "\tSNP Positions assigned to A-Genome\tSNP Positions assigned to B-Genome";
	if (OUT_File_D.length() > 0) {
		cout << "\tSNP Positions assigned to D-Genome";
	}
	cout << "\tSNPs in A-Diploid/CS-Genome\tSNPs in B-Diploid/CS-Genome";
	if (OUT_File_D.length() > 0) {
		cout << "\tSNPs in D-Diploid/CS-Genome";
	}
	cout << "\tPositions ignored for A-Genome (Low Coverage)\tPositions ignored for B-Genome (Low Coverage)";
	if (OUT_File_D.length() > 0) {
		cout << "\tPositions ignored for D-Genome (Low Coverage)";
	}	
	cout << endl;

	string strSAM = "";
	string Read1 = "";
	string Read2 = "";
	string reference1 = "";
	string reference2 = "";
	string previous_reference = "";
	map< string, string > PendingReads;
	
	set< set <SNP> > ReadSNPGroups;
	map<string, int> ReadSNPReadCount;
	
	vector< SNPGroup > RefSNPGroups;
	map < SNP, int > RefSNPReadCount; // No of reads supporting a SNP

	char* RefSNPList;
	char* RefASNPList;
	char* RefBSNPList;
	char* RefDSNPList;
	char* RefDiploidASNPList = NULL;
	char* RefDiploidBSNPList = NULL;
	char* RefDiploidDSNPList = NULL;	
	char* RefRefBaseList;
	int* RefDiploidACoverageList = NULL;
	int* RefDiploidBCoverageList = NULL;
	int* RefDiploidDCoverageList = NULL;
	vector< pair<long, long> > RefGFFData;
	vector< pair<long, long> >::iterator itRefGFFData;
	bool GetNextRead = true;
	bool InvalidGene = false;
	bool ReadPairMissing = false;
	// Initialise
	// clear SNP Groups
	RefSNPGroups.clear();
	// clear No. of Reads supporting a SNP
	RefSNPReadCount.clear();
	// clear pending list of reads
	PendingReads.clear();				
	// clear Read SNP Groups
	ReadSNPGroups.clear();
	// clear No. of reads having a SNP group
	ReadSNPReadCount.clear();
	
	while (ifs_sam.good()) {
		if (GetNextRead) {
			getline(ifs_sam, strSAM);
		}
		
		GetNextRead = true;
		if (strSAM.empty()) { // ignore emtry lines
			continue;
		} else if (strSAM[0] == '@') { // ignore the header
			continue;
		}
		
		// get individual fields
		vector<string> entries;
		strsplit(strSAM, entries, "\t");

		// extract the read name
		string ReadName = entries[0];
		
		
		if (atoi(entries[1].c_str()) & FLAG_MATE_UNMAPPED) {
			Read1 = strSAM;
			Read2 = Read1;
			ReadPairMissing = true;
		} else {
			// see if the read's pair has been seen before
			map< string, string >::iterator itPendingReads = PendingReads.find(ReadName);
			if (itPendingReads == PendingReads.end()) { //no, add it to the list of pending reads
				PendingReads[ReadName] = strSAM;
				continue;
			} else { // get the read pair
				Read1 = itPendingReads->second;
				Read2 = strSAM;
				
				// remove the read from the list of pending reads
				PendingReads.erase(itPendingReads);
			}
			ReadPairMissing = false;
		}
			
		// get individual fields
		vector<string> entries1;
		strsplit(Read1, entries1, "\t");
		
		vector<string> entries2;
		strsplit(Read2, entries2, "\t");
		
		string Read1Name = entries1[0];
		string Read2Name = entries2[0];
		
		// get the reference names
		reference1 = entries1[2];
		reference2 = entries2[2];
		
		if (reference1.compare(reference2) != 0) { // the reads in the pair are mapped on difference chromosomes. ignore the pair
			continue;
		}
		
		// read 1 coordinates
		long Read1Start = atol(entries1[3].c_str());
		long Read1End = GetReadEndCoordinates(Read1Start, entries1[5]);
		
		// read 2 coordinates
		long Read2Start = atol(entries2[3].c_str());
		long Read2End = GetReadEndCoordinates(Read2Start, entries2[5]);
		
		// get the SNP list and GFF data for this reference
		if (reference1.compare(previous_reference) != 0) {
			
			if (previous_reference.length() > 0) {
				// process the previous reference groups
				FinaliseSNPGroups(RefSNPGroups, RefSNPReadCount, ReadSNPGroups, ReadSNPReadCount);
				
				// Assign the snps to A, B & D genome
				AssignSNPsToGenomes(previous_reference, RefSNPGroups, RefRefBaseList, RefDiploidASNPList, RefDiploidBSNPList, RefDiploidDSNPList, RefDiploidACoverageList, RefDiploidBCoverageList, RefDiploidDCoverageList, ofs_A, ofs_B, ofs_D);
								
				// clear SNP Groups
				RefSNPGroups.clear();
				// clear No. of Reads supporting a SNP
				RefSNPReadCount.clear();
				// clear pending list of reads
				PendingReads.clear();				
				// clear Read SNP Groups
				ReadSNPGroups.clear();
				// clear No. of reads having a SNP group
				ReadSNPReadCount.clear();
				
			}			
			
			// get the new data
			RefSNPList = SNPList[reference1];
			RefGFFData = GFFData[reference1];
			RefRefBaseList = RefBaseList[reference1];
			itRefGFFData = RefGFFData.begin();
			
			while (itRefGFFData != RefGFFData.end() && (Read1Start > itRefGFFData->second || Read2Start > itRefGFFData->second)) {
				itRefGFFData++;
			}
			
			if (itRefGFFData == RefGFFData.end() || (Read1End < itRefGFFData->first && Read2End < itRefGFFData->first)) {
				InvalidGene = true;
				break;
			}
			
			map < long, char > UnigeneSNPList;
			GetSNPsInTheReadRegion(RefSNPList, itRefGFFData->first - 20, itRefGFFData->second + 20, UnigeneSNPList);
			cout << itRefGFFData->first << "\t" << itRefGFFData->second << "\t" << UnigeneSNPList.size();
			//cout << "Checking reads in the region: " << itRefGFFData->first << " - " << itRefGFFData->second << endl;
			//cout << CurrentDateTime();

			RefDiploidASNPList = DiploidASNPList[reference1];
			RefDiploidBSNPList = DiploidBSNPList[reference1];
			if (DiploidDSNPList.size() > 0) {
				RefDiploidDSNPList = DiploidDSNPList[reference1];
			}

			RefDiploidACoverageList = DiploidACoverageList[reference1];
			RefDiploidBCoverageList = DiploidBCoverageList[reference1];
			if (DiploidDCoverageList.size() > 0) {
				RefDiploidDCoverageList = DiploidDCoverageList[reference1];
			}
			previous_reference = reference1;
		}		
		
		if (Read1Start > itRefGFFData->second || Read2Start > itRefGFFData->second) {
			// read lies after the current gene
			
			// process the previous reference groups
			FinaliseSNPGroups(RefSNPGroups, RefSNPReadCount, ReadSNPGroups, ReadSNPReadCount);
			
			// Assign the snps to A, B & D genome
			AssignSNPsToGenomes(reference1, RefSNPGroups, RefRefBaseList, RefDiploidASNPList, RefDiploidBSNPList, RefDiploidDSNPList, RefDiploidACoverageList, RefDiploidBCoverageList, RefDiploidDCoverageList, ofs_A, ofs_B, ofs_D);
			
			// clear SNP Groups
			RefSNPGroups.clear();
			// clear No. of Reads supporting a SNP
			RefSNPReadCount.clear();
			// clear pending list of reads
			PendingReads.clear();
			// clear Read SNP Groups
			ReadSNPGroups.clear();
			// clear No. of reads having a SNP group
			ReadSNPReadCount.clear();
			
			GetNextRead = false;
			while (itRefGFFData != RefGFFData.end() && (Read1Start > itRefGFFData->second || Read2Start > itRefGFFData->second)) {
				itRefGFFData++;
			}
			
			if (itRefGFFData == RefGFFData.end() || (Read1End < itRefGFFData->first && Read2End < itRefGFFData->first)) {
//				cout << Read1Name << "\t" << Read1Start << "\t" << Read2Start << "\t" << Read1End << "\t" << Read2End << "\t" << (itRefGFFData == RefGFFData.end() ? -1 : itRefGFFData->first) << endl;
				InvalidGene = true;
				break;
			}
			
			map < long, char > UnigeneSNPList;
			GetSNPsInTheReadRegion(RefSNPList, itRefGFFData->first - 20, itRefGFFData->second + 20, UnigeneSNPList);
			cout << itRefGFFData->first << "\t" << itRefGFFData->second << "\t" << UnigeneSNPList.size();
			//cout << "Checking reads in the region: " << itRefGFFData->first << " - " << itRefGFFData->second << endl;
			//cout << CurrentDateTime();

			continue;
		}
		
		if (Read2Start < Read1Start ) { // swap the reads if read1 starts after read2
			// SAM file entries
			swap(Read1, Read2);
			
			// read name
			swap(Read1Name, Read2Name);
			
			// entries
			swap(entries1, entries2);
			
			// read start
			swap(Read1Start, Read2Start);
			
			// read end
			swap(Read1End, Read2End);			
		}
		
		// get the SNPs that lie in the region of these reads
		map < long, char > Read1RegionSNPList, Read2RegionSNPList;
		GetSNPsInTheReadRegion(RefSNPList, Read1Start, Read1End, Read1RegionSNPList);//, Read1RegionSNPPositions);
		if (!ReadPairMissing) {
			GetSNPsInTheReadRegion(RefSNPList, Read2Start, Read2End, Read2RegionSNPList);//, Read2RegionSNPPositions);		
		}
		
		if ( Read1RegionSNPList.size() == 0 &&  Read2RegionSNPList.size() == 0) { // ignore the read pair if both reads don't have any SNP within their boundaries
			continue;
		}
		
		// Get the optional fields
		map<string, string> Fields1;
		GetOptionalFields(entries1, Fields1);
		string MD1 = Fields1["MD"]; // mismatiching positions/bases in Read 1
		string Read1Seq = entries1[9];
		string BaseQualities1 = entries1[10];
		set< SNP > Read1SNPList; 
				
		// get read 1 SNPs
		bool AllGoodQualitySNPs = GetReadSNPs(Read1RegionSNPList, Read1Start, MD1, Read1Seq, BaseQualities1, Read1SNPList);
		if (!AllGoodQualitySNPs) { // poor data present at SNP positions .. ignore
			continue;
		} 
		
		if (!ReadPairMissing) {
			// Get the optional fields
			map<string, string> Fields2;
			GetOptionalFields(entries2, Fields2);
			string Read2Seq = entries2[9];
			string MD2 = Fields2["MD"]; // mismatiching positions/bases in Read 2
			string BaseQualities2 = entries2[10];
			set< SNP > Read2SNPList;
			
			// get read 2 SNPs
			AllGoodQualitySNPs = GetReadSNPs(Read2RegionSNPList, Read2Start, MD2, Read2Seq, BaseQualities2, Read2SNPList);
			if (!AllGoodQualitySNPs) { // poor data present at SNP positions .. ignore
				continue;
			} 
			
			if (Read1End >= Read2Start && CheckDiscrepancyInReadPair(Read1SNPList, Read2SNPList)) {
				//cout << "Discrepancy in Read pair: " << Read1Name << endl;
				continue;
			}

			// add Read 2 SNPs to the Read 1 SNP list to get a combined SNP list
			Read1SNPList.insert(Read2SNPList.begin(), Read2SNPList.end());
		}
		
		
		
		if (Read1SNPList.size() > 0) {
			// we have the list of valid SNPs present in this Read pair. 
			
/*			cout << Read1Name << " " << Read1Start << "/" << Read2Start << ": ";
			set<SNP>::iterator itReadSNPList = Read1SNPList.begin();
			for (; itReadSNPList != Read1SNPList.end(); itReadSNPList++) {
				cout << itReadSNPList->first << "," << itReadSNPList->second << " ";
			}
			cout << endl;
*/			
			// Put them in their correct group
			ReadSNPGroups.insert(Read1SNPList);
			// Generate a Group Id
			string GroupId = CreateSNPGroupId(Read1SNPList);
			// Increment the number of reads containing this SNP Group
			ReadSNPReadCount[GroupId]++;
			
			
		} // if Read SNP list size > 0
		
	} // end if sam file is good
	
	if (!InvalidGene) {
		// process the final reference groups
		FinaliseSNPGroups(RefSNPGroups, RefSNPReadCount, ReadSNPGroups, ReadSNPReadCount);
		
		// Assign the snps to A, B & D genome
		AssignSNPsToGenomes(reference1, RefSNPGroups, RefRefBaseList, RefDiploidASNPList, RefDiploidBSNPList, RefDiploidDSNPList, RefDiploidACoverageList, RefDiploidBCoverageList, RefDiploidDCoverageList, ofs_A, ofs_B, ofs_D);
	}
	
	ifs_sam.close();
	ofs_A.close();
	ofs_B.close();
	if (OUT_File_D.length() > 0) {
		ofs_D.close();
	}
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string SAM_File (argv[1]);
	string GFF_File (argv[2]);
	string Hexaploid_SNP_File (argv[3]);
	string Diploid_SNP_File_A (argv[4]);
	string Diploid_SNP_File_B (argv[5]);
	string Diploid_SNP_File_D (argv[6]);
	string Diploid_Coverage_File_A (argv[7]); // this coverage file needs to be generated from pileup file without any consensus option
	string Diploid_Coverage_File_B (argv[8]);
	string Diploid_Coverage_File_D (argv[9]);
	string OUT_File_A (argv[10]);
	string OUT_File_B (argv[11]);
	string OUT_File_D (argv[12]);
	
	
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
	
	// variable to hold SNP and other related lists 
	SNP_LIST HexaploidSNPList;
	SNP_LIST DiploidSNPList_A;
	SNP_LIST DiploidSNPList_B;
	SNP_LIST DiploidSNPList_D;
	map < string, char* > RefBaseList;
	map < string, long > RefLength;
	map < string, vector< pair<long, long> > > GFFData;
	COVERAGE_LIST DiploidCoverageList_A;
	COVERAGE_LIST DiploidCoverageList_B;
	COVERAGE_LIST DiploidCoverageList_D;

	// Read the genes' coordinates from the GFF file
	ReadGFFFile(GFF_File, GFFData);
	
	cout << 1 << endl;

	// initialise SNP lists
	InitialiseSNPList(SAM_File, HexaploidSNPList, RefBaseList, RefLength);	
	InitialiseSNPList(RefLength, DiploidSNPList_A);
	InitialiseSNPList(RefLength, DiploidSNPList_B);
	if (Diploid_SNP_File_D.length() > 0) {
		InitialiseSNPList(RefLength, DiploidSNPList_D);
	}

	cout << 2 << endl;
	
	// Read SNPs
	ReadSNPList(Hexaploid_SNP_File, HexaploidSNPList, RefBaseList);
	ReadSNPList(Diploid_SNP_File_A, DiploidSNPList_A);
	ReadSNPList(Diploid_SNP_File_B, DiploidSNPList_B);
	if (Diploid_SNP_File_D.length() > 0) {
		ReadSNPList(Diploid_SNP_File_D, DiploidSNPList_D);
	}

	cout << 3 << endl;
	
	// Read Coverage
	ReadCoverageFile(Diploid_Coverage_File_A, RefLength, DiploidCoverageList_A);
	ReadCoverageFile(Diploid_Coverage_File_B, RefLength, DiploidCoverageList_B);
	if (Diploid_Coverage_File_D.length() > 0) {
		ReadCoverageFile(Diploid_Coverage_File_D, RefLength, DiploidCoverageList_D);
	}

	cout << 4 << endl;
	
	// process the sam file
	ProcessSAMFile(SAM_File, GFFData, HexaploidSNPList, RefBaseList, DiploidSNPList_A, DiploidSNPList_B, DiploidSNPList_D, DiploidCoverageList_A, DiploidCoverageList_B, DiploidCoverageList_D, OUT_File_A, OUT_File_B, OUT_File_D);
	
	return 0;
}

