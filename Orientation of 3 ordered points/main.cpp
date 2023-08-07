#include <iostream>

using namespace std;

struct vertex
{
    double x;
    double y;
};

int orientation (vertex p1, vertex p2, vertex p3 )
{

    int val = ((p2.y-p1.y)/(p2.x-p1.x)-(p3.y-p2.y)/(p3.x-p2.x));

      if (val == 0)
        return 0;
      else
       return (val > 0)? 1: 2;
        //0 for Linear
        //1 for clockwise
        //2 for counter clockwise

}

int main()
{
    vertex p1 = {0, 0}, p2 = {4, 4}, p3 = {5,3};
    int o = orientation(p1, p2, p3);
    if (o==0)         cout << "Linear";
    else if (o == 1)  cout << "Clockwise";
    else              cout << "Counter-Clockwise";
    return 0;
}
