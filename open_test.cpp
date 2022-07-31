#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(){
  char name[30];
  int i = 0;
  sprintf(name, "pos_vel/pos_vel_%d.dat", i);
  ifstream IN(name);
  if(! IN){
    IN.close();
    cout << "error\n";
  }else
    cout << "success!\n";
  IN.close();
  return 0;
}
