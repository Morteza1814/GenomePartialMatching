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
	ifstream file("/bigtemp/rgq5aw/RMYN02.fasta");
	ofstream MyExcelFile;
	MyExcelFile.open("/bigtemp/rgq5aw/RMY.csv");
	string str;
	int count = 0, count150 = 0;
	char ch;
	string refParts;
	cout << getline(file, str) << "\n";
	while (file >> ch)
	{
		//cout << str <<"\n";
		if (ch == '\n')
		{
			continue;
		}
		else if (count150 < 150)
		{
			//		cout << "str = " << refParts << "\n";
			refParts.push_back(ch);
			count150++;
		}
		else
		{
			count150 = 0;
			cout << refParts << "\n";
			count++;
			MyExcelFile /*<< "ref_" << count << ", " */ << refParts << endl;
			refParts = "";
			refParts.push_back(ch);
			count150++;
		}
	}

	MyExcelFile.close();

	file.close();
}
