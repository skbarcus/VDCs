#include <string>
#include <vector>
using namespace std;

std::vector<Double_t> v1 = 0.;

void vector_test()
{
  v1.push_back(1.0);
  v1.push_back(2.0);

  for(Int_t i=0;i<v1.size();i++)
    {
      cout<<v1.at(i)<<endl;
    }
  v1.clear();
  for(Int_t i=0;i<v1.size();i++)
    {
      cout<<v1.at(i)<<endl;
    }
}
