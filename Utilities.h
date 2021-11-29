/*
 *  Utilities.h
 *  
 *
 *  Created by Aziz Mithani on 20/05/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>

using namespace std;

//void strsplit(const string& str, vector<string>& tokens, const string& delimiters);

//vector< vector<string> > ReadFile(string FileName);

const int TOTAL_BASES = 5;

void InitialiseTranslationTable(map<string, string>& Translation) {
	
	Translation["AAA"] = "K";
	Translation["AAC"] = "N";
	Translation["AAG"] = "K";
	Translation["AAT"] = "N";
	Translation["ACA"] = "T";
	Translation["ACC"] = "T";
	Translation["ACG"] = "T";
	Translation["ACT"] = "T";
	Translation["AGA"] = "R";
	Translation["AGC"] = "S";
	Translation["AGG"] = "R";
	Translation["AGT"] = "S";
	Translation["ATA"] = "I";
	Translation["ATC"] = "I";
	Translation["ATG"] = "M";
	Translation["ATT"] = "I";
	Translation["CAA"] = "Q";
	Translation["CAC"] = "H";
	Translation["CAG"] = "Q";
	Translation["CAT"] = "H";
	Translation["CCA"] = "P";
	Translation["CCC"] = "P";
	Translation["CCG"] = "P";
	Translation["CCT"] = "P";
	Translation["CGA"] = "R";
	Translation["CGC"] = "R";
	Translation["CGG"] = "R";
	Translation["CGT"] = "R";
	Translation["CTA"] = "L";
	Translation["CTC"] = "L";
	Translation["CTG"] = "L";
	Translation["CTT"] = "L";
	Translation["GAA"] = "E";
	Translation["GAC"] = "D";
	Translation["GAG"] = "E";
	Translation["GAT"] = "D";
	Translation["GCA"] = "A";
	Translation["GCC"] = "A";
	Translation["GCG"] = "A";
	Translation["GCT"] = "A";
	Translation["GGA"] = "G";
	Translation["GGC"] = "G";
	Translation["GGG"] = "G";
	Translation["GGT"] = "G";
	Translation["GTA"] = "V";
	Translation["GTC"] = "V";
	Translation["GTG"] = "V";
	Translation["GTT"] = "V";
	Translation["TAA"] = "*";
	Translation["TAC"] = "Y";
	Translation["TAG"] = "*";
	Translation["TAT"] = "Y";
	Translation["TCA"] = "S";
	Translation["TCC"] = "S";
	Translation["TCG"] = "S";
	Translation["TCT"] = "S";
	Translation["TGA"] = "*";
	Translation["TGC"] = "C";
	Translation["TGG"] = "W";
	Translation["TGT"] = "C";
	Translation["TTA"] = "L";
	Translation["TTC"] = "F";
	Translation["TTG"] = "L";
	Translation["TTT"] = "F";
	Translation["aaa"] = "K";
	Translation["aac"] = "N";
	Translation["aag"] = "K";
	Translation["aat"] = "N";
	Translation["aca"] = "T";
	Translation["acc"] = "T";
	Translation["acg"] = "T";
	Translation["act"] = "T";
	Translation["aga"] = "R";
	Translation["agc"] = "S";
	Translation["agg"] = "R";
	Translation["agt"] = "S";
	Translation["ata"] = "I";
	Translation["atc"] = "I";
	Translation["atg"] = "M";
	Translation["att"] = "I";
	Translation["caa"] = "Q";
	Translation["cac"] = "H";
	Translation["cag"] = "Q";
	Translation["cat"] = "H";
	Translation["cca"] = "P";
	Translation["ccc"] = "P";
	Translation["ccg"] = "P";
	Translation["cct"] = "P";
	Translation["cga"] = "R";
	Translation["cgc"] = "R";
	Translation["cgg"] = "R";
	Translation["cgt"] = "R";
	Translation["cta"] = "L";
	Translation["ctc"] = "L";
	Translation["ctg"] = "L";
	Translation["ctt"] = "L";
	Translation["gaa"] = "E";
	Translation["gac"] = "D";
	Translation["gag"] = "E";
	Translation["gat"] = "D";
	Translation["gca"] = "A";
	Translation["gcc"] = "A";
	Translation["gcg"] = "A";
	Translation["gct"] = "A";
	Translation["gga"] = "G";
	Translation["ggc"] = "G";
	Translation["ggg"] = "G";
	Translation["ggt"] = "G";
	Translation["gta"] = "V";
	Translation["gtc"] = "V";
	Translation["gtg"] = "V";
	Translation["gtt"] = "V";
	Translation["taa"] = "*";
	Translation["tac"] = "Y";
	Translation["tag"] = "*";
	Translation["tat"] = "Y";
	Translation["tca"] = "S";
	Translation["tcc"] = "S";
	Translation["tcg"] = "S";
	Translation["tct"] = "S";
	Translation["tga"] = "*";
	Translation["tgc"] = "C";
	Translation["tgg"] = "W";
	Translation["tgt"] = "C";
	Translation["tta"] = "L";
	Translation["ttc"] = "F";
	Translation["ttg"] = "L";
	Translation["ttt"] = "F";
	
}

void InitialiseNucleotideMap(map <string, set<string> >& NucleotideMap){
	
	set<string> A; A.insert("A");
	set<string> C; C.insert("C");
	set<string> G; G.insert("G");
	set<string> T; T.insert("T");
	set<string> R; R.insert("A"); R.insert("G");
	set<string> Y; Y.insert("C"); Y.insert("T");
	set<string> M; M.insert("A"); M.insert("C");
	set<string> K; K.insert("G"); K.insert("T");
	set<string> S; S.insert("C"); S.insert("G");
	set<string> W; W.insert("A"); W.insert("T");
	set<string> H; H.insert("A"); H.insert("C"); H.insert("T");
	set<string> B; B.insert("C"); B.insert("G"); B.insert("T");
	set<string> D; D.insert("A"); D.insert("G"); D.insert("T");
	set<string> V; V.insert("A"); V.insert("C"); V.insert("G");
	set<string> N; N.insert("N");
	NucleotideMap["A"] = A;
	NucleotideMap["C"] = C;
	NucleotideMap["G"] = G;
	NucleotideMap["T"] = T;
	NucleotideMap["R"] = R;
	NucleotideMap["Y"] = Y;
	NucleotideMap["M"] = M;
	NucleotideMap["K"] = K;
	NucleotideMap["S"] = S;
	NucleotideMap["W"] = W;
	NucleotideMap["H"] = H;
	NucleotideMap["B"] = B;
	NucleotideMap["D"] = D;
	NucleotideMap["V"] = V;
	NucleotideMap["N"] = N;
	
}

void InitialiseNucleotideMapReverse(map <string, string >& NucleotideMap){
	
	NucleotideMap["A"] = "A";
	NucleotideMap["C"] = "C";
	NucleotideMap["G"] = "G";
	NucleotideMap["T"] = "T";
	NucleotideMap["AG"] = "R";
	NucleotideMap["CT"] = "Y";
	NucleotideMap["AC"] = "M";
	NucleotideMap["GT"] = "K";
	NucleotideMap["CG"] = "S";
	NucleotideMap["AT"] = "W";
	NucleotideMap["ACT"] = "H";
	NucleotideMap["CGT"] = "B";
	NucleotideMap["AGT"] = "D";
	NucleotideMap["ACG"] = "V";
	NucleotideMap["ACGT"] = "N";
	NucleotideMap["N"] = "N";
	
}

string SetToString(set<string> TheSet, const string delimiter = "") {
	
	if (TheSet.size() == 0) {
		return "";
	}
	
	set<string>::iterator itTheSet = TheSet.begin();
	string str = *itTheSet;
	itTheSet++;
	for (; itTheSet != TheSet.end(); itTheSet++) {
		str += delimiter + *itTheSet;
	}
	
	return str;
	
}
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

vector< vector<string> > ReadFile(string FileName, int skip = 0){
	// skip: no. of lines to skip
	
	
	// open the file stream
	ifstream ifs ( FileName.c_str() );
	
	string str="";
	vector< vector<string> > v;
	int line = 0;
	while (line++ < skip && ifs.good()) {
		getline(ifs, str); // throw away the first line
	}
	while (ifs.good()) {
		getline(ifs, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		v.push_back (entries);
		
	}
	
	ifs.close();
	
	return v;
}

void ReadFile2(string FileName, map <string, string>& Data, int skip = 0){
	// skip: no. of lines to skip
	
	
	// open the file stream
	ifstream ifs ( FileName.c_str() );
	
	string str="";
	int line = 0;
	while (line++ < skip && ifs.good()) {
		getline(ifs, str); // throw away the lines to be skipped
	}
	while (ifs.good()) {
		getline(ifs, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
	
		// split the data
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		// add to the map
		Data[entries[0]] = entries[1];

	}
	
	ifs.close();
	
}

void TrimSpaces( string& str) {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af
	
    // if all spaces or empty return an empty string
    if(( string::npos == startpos ) || ( string::npos == endpos)) {
        str = "";
    }
    else
        str = str.substr( startpos, endpos-startpos+1 );
}

void TrimLeadingSpaces( string& str) {
 	 size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
	 if( string::npos != startpos )
	 str = str.substr( startpos );
}

void TrimTrailingSpaces( string& str) {
	 size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse after excluding trailing blank spaces
	 if( string::npos != endpos )
	 str = str.substr( 0, endpos+1 );
}

void RemoveNewLine( string& str) {
	str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
	str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
}

string ReverseComplementDNA(string dna) {
	map<char, char> Complementation;
	Complementation['A'] = 'T';
	Complementation['T'] = 'A';
	Complementation['C'] = 'G';
	Complementation['G'] = 'C';
	Complementation['a'] = 't';
	Complementation['t'] = 'a';
	Complementation['c'] = 'g';
	Complementation['g'] = 'c';
	
	string result (dna.rbegin(), dna.rend());
	for (int i = 0; i < result.length(); i++) {
		map<char, char>::iterator it;
		if((it = Complementation.find(result[i])) != Complementation.end())
			result[i] = it->second;
	}
	
	return result;
}

string TranslateDNASequence(string& dna) {
	
	map<string, string> TranslationTable;
	InitialiseTranslationTable(TranslationTable);
	
	if (dna.length() < 3) {
		return "";
	}

	string protein = "";
	for (int i = 0; i < dna.length() - 3; i += 3) {
		map<string, string>::iterator it = TranslationTable.find(dna.substr(i,3));
		protein += (it == TranslationTable.end() ? "X" : it->second);
	}
		
	return protein;
}

string TranslateDNASequence(string& dna, int frame) {

	if (frame > 2 || frame < -2) {
		return ""; 
	}
	
	map<string, string> TranslationTable;
	InitialiseTranslationTable(TranslationTable);
	
	if (dna.length() < 3) {
		return "";
	}
	
	string the_dna;
	if (frame < 0) {
		the_dna = ReverseComplementDNA(dna);
		frame *= -1;
	} else {
		the_dna = dna;
	}

	
	string protein = "";
	for (int i = frame; i <= the_dna.length() - 3; i += 3) {
		map<string, string>::iterator it = TranslationTable.find(the_dna.substr(i,3));
		protein += (it == TranslationTable.end() ? "X" : it->second);
	}
	
	return protein;
}

string TranslateDNA(string dna) {
	string StartCodon = "M";
	string StopCodon = "*";
	
	map<string,string> TranslationTable;
	InitialiseTranslationTable(TranslationTable);
	
	if (dna.length() < 3) {
		return "";
	}
	
	string SelectedProtein = "";
	for (int f = 0; f < 3; f++) {
		string protein = "";
		bool TranslationStarted = false;
		for (int i = f; i < dna.length() - 3; i += 3) {
			map<string, string>::iterator it = TranslationTable.find(dna.substr(i,3));
			string codon = (it == TranslationTable.end() ? "X" : it->second);

			//if (dna.substr(i,3).find("N") != string::npos) {
			//	cout << "Ambigous\t" << dna.substr(i,3) << "\t" << codon << endl;
			//}
			if (!TranslationStarted && codon.compare(StartCodon) == 0) {
				TranslationStarted = true;
			}
			if (TranslationStarted && codon.compare(StopCodon) == 0) {
				if (protein.length() > SelectedProtein.length())
					SelectedProtein = protein;
				protein = "";
				TranslationStarted = false;
			}
			if (TranslationStarted) {
				protein += codon;//Translation[dna.substr(i,3)];
			}
		}
		//cout << protein.length() << endl;
		if (protein.length() > SelectedProtein.length())
			SelectedProtein = protein;
		
	} // for all reading frames	
	
	return SelectedProtein;
}

int TranslateDNA(string& dna, string& protein) {
	/* returns the frame: 0, 1, or 2*/
	
	string StartCodon = "M";
	string StopCodon = "*";
	
	map<string,string> TranslationTable;
	InitialiseTranslationTable(TranslationTable);
	
	if (dna.length() < 3) {
		protein = "";
		return -1;
	}
	
	string SelectedProtein = "";
	int SelectedFrame = -1;
	for (int f = 0; f < 3; f++) {
		string current_protein = "";
		bool TranslationStarted = false;
		for (int i = f; i < dna.length() - 3; i += 3) {
			map<string, string>::iterator it = TranslationTable.find(dna.substr(i,3));
			string codon = (it == TranslationTable.end() ? "X" : it->second);
			
			//if (dna.substr(i,3).find("N") != string::npos) {
			//	cout << "Ambigous\t" << dna.substr(i,3) << "\t" << codon << endl;
			//}
			if (!TranslationStarted && codon.compare(StartCodon) == 0) {
				TranslationStarted = true;
			}
			if (TranslationStarted && codon.compare(StopCodon) == 0) {
				if (current_protein.length() > SelectedProtein.length()) {
					SelectedProtein = current_protein;
					SelectedFrame = f;
				}
				current_protein = "";
				TranslationStarted = false;
			}
			if (TranslationStarted) {
				current_protein += codon;//Translation[dna.substr(i,3)];
			}
		}
		//cout << protein.length() << endl;
		if (current_protein.length() > SelectedProtein.length()) {
			SelectedProtein = current_protein;
			SelectedFrame = f;
		}
		
	} // for all reading frames	
	
	protein = SelectedProtein;
	return SelectedFrame;
	
}

/*
int FindReadingFrame(string& dna, string& protein, bool UseStartCodon) {
	// returns the frame: 1, 2, 3, -1, -2 or -3 
	
	string StartCodon = "M";
	string StopCodon = "*";
	
	map<string,string> TranslationTable;
	InitialiseTranslationTable(TranslationTable);

	
	if (dna.length() < 3) {
		protein = "";
		return 0;
	}
	
	string SelectedProtein = "";
	int SelectedFrame = 0;
	// forward frames
	for (int f = 0; f < 3; f++) {
		string current_protein = "";
		bool TranslationStarted = (UseStartCodon ? false : true);
		for (int i = f; i < dna.length() - 3; i += 3) {
			map<string, string>::iterator it = TranslationTable.find(dna.substr(i,3));
			string codon = (it == TranslationTable.end() ? "X" : it->second);

			if (!TranslationStarted && codon.compare(StartCodon) == 0) {
				TranslationStarted = true;
			} else if (TranslationStarted) {
				if (codon.compare(StopCodon) == 0) {
					if (current_protein.length() > SelectedProtein.length()) {
						SelectedProtein = current_protein;
						SelectedFrame = f + 1;
					}
					current_protein = "";
					TranslationStarted = (UseStartCodon ? false : true);
				}
				current_protein += codon;
			}

		}
		//cout << protein.length() << endl;
		if (current_protein.length() > SelectedProtein.length()) {
			SelectedProtein = current_protein;
			SelectedFrame = f + 1;
		}
	} // for all reading frames	

	string rdna = ReverseComplementDNA(dna);
	// reverse frames
	for (int f = 0; f < 3; f++) {
		string current_protein = "";
		bool TranslationStarted = (UseStartCodon ? false : true);
		for (int i = f; i < rdna.length() - 3; i += 3) {
			map<string, string>::iterator it = TranslationTable.find(rdna.substr(i,3));
			string codon = (it == TranslationTable.end() ? "X" : it->second);
			
			if (!TranslationStarted && codon.compare(StartCodon) == 0) {
				TranslationStarted = true;
			} else if (TranslationStarted) {
				if (codon.compare(StopCodon) == 0) {
					if (current_protein.length() > SelectedProtein.length()) {
						SelectedProtein = current_protein;
						SelectedFrame = (f + 1) * -1;
					}
					current_protein = "";
					TranslationStarted = (UseStartCodon ? false : true);
				}
				current_protein += codon;
			}
			
		}
		//cout << protein.length() << endl;
		if (current_protein.length() > SelectedProtein.length()) {
			SelectedProtein = current_protein;
			SelectedFrame = (f + 1) * -1;
		}
	} // for all reading frames	
	
	protein = SelectedProtein;
	return SelectedFrame;
	
}
*/

int FindReadingFrame(string& dna, string& protein, int& qStart, int& qEnd, bool UseStartCodon) {
	/* returns the frame: 1, 2, 3, -1, -2 or -3 */
	
	string StartCodon = "M";
	string StopCodon = "*";
	
	map<string,string> TranslationTable;
	InitialiseTranslationTable(TranslationTable);
	
	
	if (dna.length() < 3) {
		protein = "";
		return 0;
	}
	
	string SelectedProtein = "";
	int SelectedStart = 0;
	int SelectedEnd = 0;
	int SelectedFrame = 0;
	// forward frames
	for (int f = 0; f < 3; f++) {
		string current_protein = "";
		bool TranslationStarted = (UseStartCodon ? false : true);
		for (int i = f; i < dna.length() - 3; i += 3) {
			map<string, string>::iterator it = TranslationTable.find(dna.substr(i,3));
			string codon = (it == TranslationTable.end() ? "X" : it->second);
			if (!TranslationStarted && codon.compare(StartCodon) == 0) {
				TranslationStarted = true;
				current_protein = codon;
				qStart = i;
			} else if (TranslationStarted) {
				qEnd = i+3;
				if (codon.compare(StopCodon) == 0) {
					if (current_protein.length() > SelectedProtein.length()) {
						SelectedProtein = current_protein;
						SelectedFrame = f + 1;
						SelectedStart = qStart;
						SelectedEnd = qEnd;
					}
					current_protein = "";
					TranslationStarted = (UseStartCodon ? false : true);
				} else {
					if (current_protein.length() == 0) {
						qStart = i;
					}
					current_protein += codon;
				}
			}
			
		}
		if (current_protein.length() != 0) {
			qEnd = dna.length();
		}
		
		//cout << protein.length() << endl;
		if (current_protein.length() > SelectedProtein.length()) {
			SelectedProtein = current_protein;
			SelectedFrame = f + 1;
			SelectedStart = qStart;
			SelectedEnd = qEnd;
		}
	} // for all reading frames	
	
	string rdna = ReverseComplementDNA(dna);
	// reverse frames
	for (int f = 0; f < 3; f++) {
		string current_protein = "";
		bool TranslationStarted = (UseStartCodon ? false : true);
		for (int i = f; i < rdna.length() - 3; i += 3) {
			map<string, string>::iterator it = TranslationTable.find(rdna.substr(i,3));
			string codon = (it == TranslationTable.end() ? "X" : it->second);
			
			if (!TranslationStarted && codon.compare(StartCodon) == 0) {
				TranslationStarted = true;
				current_protein = codon;
				qStart = i;
			} else if (TranslationStarted) {
				qEnd = i+3;
				if (codon.compare(StopCodon) == 0) {
					if (current_protein.length() > SelectedProtein.length()) {
						SelectedProtein = current_protein;
						SelectedFrame = (f + 1) * -1;
						SelectedStart = qStart;
						SelectedEnd = qEnd;
					}
					current_protein = "";
					TranslationStarted = (UseStartCodon ? false : true);
				} else {
					if (current_protein.length() == 0) {
						qStart = i;
					}

					current_protein += codon;
				}
			}
			
		}
		if (current_protein.length() != 0) {
			qEnd = dna.length();
		}
		//cout << protein.length() << endl;
		if (current_protein.length() > SelectedProtein.length()) {
			SelectedProtein = current_protein;
			SelectedFrame = (f + 1) * -1;
			SelectedStart = qStart;
			SelectedEnd = qEnd;
		}
	} // for all reading frames	
	
	protein = SelectedProtein;
	
	if (SelectedFrame < 0) {
		qStart = dna.length() - SelectedEnd + 1;
		qEnd = dna.length() - SelectedStart + 1;
		//cout << dna.length() << "\t" << SelectedStart << "\t" << SelectedEnd << "\t" << qStart << "\t" << qEnd << endl;
	} else {
		qStart = SelectedStart;
		qEnd = SelectedEnd;
	}

	
	return SelectedFrame;

}

string CurrentDateTime() {
	time_t rawtime;
	struct tm * timeinfo;
	
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	
	return asctime (timeinfo);
}

int GetBaseIndex(char base) {

	switch (base) {
		case 'A':
			return 0;
		case 'T':
			return 1;
		case 'G':
			return 2;
		case 'C':
			return 3;
		case 'N':
			return 4;
		default:
			return -1;
	}
}

char GetIndexBase(int index) {
	
	switch (index) {
		case 0:
			return 'A';
		case 1:
			return 'T';
		case 2:
			return 'G';
		case 3:
			return 'C';
		case 4:
			return 'N';
		default:
			return ' ';
	}
}

int GetBaseIndexBaseDistFile(char base) {
	
	switch (base) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		case 'N':
			return 4;
		default:
			return -1;
	}
}

char GetIndexBaseBaseDistFile(int index) {
	
	switch (index) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		case 4:
			return 'N';
		default:
			return ' ';
	}
}

bool FileExists(string filename) {
	ifstream ifile(filename.c_str());
	return ifile.good();
}

#endif
