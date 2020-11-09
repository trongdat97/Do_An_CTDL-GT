#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include<iomanip>
using namespace std;

const int N = 10000;

int main(){
	srand(time(NULL));
	ofstream outfile;
	outfile.open("input.txt");
    for (int i = 0; i<N; i++ ) outfile<< setw(7) << left << rand() << " ";
    outfile.close();
    cout<<"Da tao mang gom "<<N<<" phan tu va luu vao file input.txt!";
    return 0;
}
