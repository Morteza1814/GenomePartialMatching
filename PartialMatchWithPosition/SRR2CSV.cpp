#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;
using namespace std;

int main()
{
  ifstream file("/bigtemp/rgq5aw/SRR12464727.fasta");
  ofstream MyExcelFile;
  MyExcelFile.open("/bigtemp/rgq5aw/SRR.csv");
  string str;
  int count = 0;

  string refParts;
  //std::cout << std::getline(file, str) << "\n";
  while (std::getline(file, str) && count < 1000000)
  {
    //cout << str <<"\n";
    if (str.find("SRR") == std::string::npos)
    {
      //		cout << "str = " << str << "\n";
      refParts += str;
    }
    else
    {
      //		cout << "ref = " << refParts << "\n";
      count++;
      MyExcelFile /*<< "ref_" << count << ", " */ << refParts << endl;
      refParts = "";
    }
  }

  MyExcelFile.close();

  file.close();
}
