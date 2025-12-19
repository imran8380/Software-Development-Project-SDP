#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

//ORF STRUCT
struct ORF
{
    int start;
    int end;
    string sequence;
};

//SEQUENCE CLASS
class Sequence
{
private:
    string seq;
    bool isRNA;

public:
    Sequence(const string &s)
    {
        seq = s;
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        isRNA = detectRNA();
        if(isRNA) convertToDNA();
    }

    bool detectRNA() const
    {
        for(char c: seq)
            if(c=='U') return true;
        return false;
    }

    bool getIsRNA() const
    {
        return isRNA;
    }
    string getSeq() const
    {
        return seq;
    }

    void convertToDNA()
    {
        for(char &c: seq)
            if(c=='U') c='T';
    }

    string convertDNAtoRNA(const string &dna) const
    {
        string rna = dna;
        for(char &c: rna)
            if(c=='T') c='U';
        return rna;
    }

    //GC CONTENT
    string calculateGC() const
    {
        int a=0, t=0, g=0, c=0;
        for(char ch: seq)
        {
            if(ch=='A') a++;
            else if(ch=='T') t++;
            else if(ch=='G') g++;
            else if(ch=='C') c++;
        }

        int total = a+t+g+c;
        if(total==0) return "Error: Invalid sequence!\n";

        double gcPercent = (double)(g+c) * 100.0 / total;

        stringstream out;
        out << "\n===== GC CONTENT RESULT =====\n";
        out << "A: " << a << ", T: " << t
            << ", G: " << g << ", C: " << c << "\n";
        out << fixed << setprecision(2)
            << "GC% = " << gcPercent << "%\n";

        return out.str();
    }

    //REVERSE COMPLEMENT
    string reverseComplement() const
    {
        string r = seq;
        reverse(r.begin(), r.end());
        for(char &b: r)
        {
            switch(b)
            {
            case 'A':
                b='T';
                break;
            case 'T':
                b='A';
                break;
            case 'G':
                b='C';
                break;
            case 'C':
                b='G';
                break;
            default:
                b='N';
            }
        }
        return r;
    }

    bool isStopCodon(const string &codon) const
    {
        return codon=="TAA" || codon=="TAG" || codon=="TGA";
    }

    //ORF FINDER
    vector<ORF> findORFs(bool forward=true) const
    {
        vector<ORF> orfs;
        string sequence = forward ? seq : reverseComplement();
        string startCodon = "ATG";
        int n = sequence.length();

        for(int frame=0; frame<3; frame++)
        {
            for(int i=frame; i+2<n; i+=3)
            {
                if(sequence.substr(i,3)==startCodon)
                {
                    for(int j=i+3; j+2<n; j+=3)
                    {
                        if(isStopCodon(sequence.substr(j,3)))
                        {
                            ORF o;
                            o.start = i+1;
                            o.end = j+3;
                            o.sequence = sequence.substr(i, j-i+3);
                            orfs.push_back(o);
                            break;
                        }
                    }
                }
            }
        }
        return orfs;
    }
};

//FASTA READER
string readFASTA(const string &filename)
{
    ifstream file(filename);
    if(!file.is_open())
    {
        cout << "Error: Cannot open file!\n";
        return "";
    }
    string line, seq="";
    while(getline(file,line))
    {
        if(!line.empty() && line[0]=='>') continue;
        seq += line;
    }
    return seq;
}

//SAVE TO FILE
void saveToFile(const string &text, const string &filename)
{
    ofstream file(filename);
    file << text;
    file.close();
    cout << "\nOutput saved to file: " << filename << endl;
}

//MAIN FUNCTION
int main()
{
    string input;
    int choice;

    cout << "Input Options:\n1. Direct input\n2. FASTA file input\nChoose (1/2): ";
    cin >> choice;

    if(choice==1)
    {
        cout << "Enter DNA/RNA sequence: ";
        cin >> input;
    }
    else if(choice==2)
    {
        string filename;
        cout << "Enter FASTA filename: ";
        cin >> filename;
        input = readFASTA(filename);
    }
    else
    {
        cout << "Invalid option!\n";
        return 0;
    }

    Sequence seqObj(input);

    //REPORT
    string report = "\n== PROCESSED SEQUENCE ==\n";
    report += input + "\n";

    //Show detected sequence type
    report += "\nSequence Type Detected: ";
    report += seqObj.getIsRNA() ? "RNA\n" : "DNA\n";

    report += seqObj.calculateGC();

    //FORWARD ORFs
    report += "\nFORWARD STRAND ORFs\n";
    vector<ORF> fOrfs = seqObj.findORFs(true);

    for(int frame=0; frame<3; frame++)
    {
        report += "Frame +" + to_string(frame+1) + ":\n";
        bool found=false;
        for(auto &o: fOrfs)
        {
            if((o.start-1)%3==frame)
            {
                report += "  Start: " + to_string(o.start) +
                          ", End: " + to_string(o.end) +
                          ", Length: " + to_string(o.sequence.size()) + " nt\n";
                string outSeq = seqObj.getIsRNA()
                                ? seqObj.convertDNAtoRNA(o.sequence)
                                : o.sequence;
                report += "  Sequence: " + outSeq + "\n";
                found=true;
            }
        }
        if(!found) report += "  No ORF found\n";
    }

    //REVERSE ORFs
    report += "\n=REVERSE STRAND ORFs=\n";
    vector<ORF> rOrfs = seqObj.findORFs(false);
    int seqLen = seqObj.getSeq().size();

    for(int frame=0; frame<3; frame++)
    {
        report += "Frame -" + to_string(frame+1) + ":\n";
        bool found=false;
        for(auto &o: rOrfs)
        {
            int origEnd = seqLen - (o.start - 1);
            int origStart = origEnd - o.sequence.size() + 1;
            if((origStart-1)%3==frame)
            {
                report += "  Start: " + to_string(origStart) +
                          ", End: " + to_string(origEnd) +
                          ", Length: " + to_string(o.sequence.size()) + " nt\n";
                string outSeq = seqObj.getIsRNA()
                                ? seqObj.convertDNAtoRNA(o.sequence)
                                : o.sequence;
                report += "  Sequence: " + outSeq + "\n";
                found=true;
            }
        }
        if(!found) report += "  No ORF found\n";
    }

    //OUTPUT
    cout << "\nOutput Options:\n1. Show in console\n2. Save to .txt file\nChoose (1/2): ";
    cin >> choice;

    if(choice==1) cout << report;
    else if(choice==2)
    {
        string filename;
        cout << "Enter output filename: ";
        cin >> filename;
        saveToFile(report, filename);
    }

    return 0;
}
