#include <iostream>

using namespace std;

struct vertex
{

    double x;
    double y;
};

// Given three points p, q, r, the function checks if point q lies on line segment 'pr'
bool onSegment(vertex p, vertex q, vertex r)
{
    if ((q.x<=max(p.x,r.x)&& q.x>=min(p.x,r.x))
       && (q.y<=max(p.y,r.y)&& q.y>=min(p.y,r.y)))

    return true;

    else
        return false;
}

// To find orientation of ordered triplet (p, q, r).
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(vertex p, vertex q, vertex r)

{
   double val = (q.y-p.y)/(q.x-p.x) - (r.y-q.y)/(r.x-q.x);

   if (val==0);
   return 0;

  return  (val>0 ? 1 : 2);

}

bool doIntersect(vertex p1, vertex q1, vertex p2, vertex q2)

{ //General Case

    int o1 = orientation(p1,q1,p2);
    int o2 = orientation(p1,q1,q2);
    int o3 = orientation(p2,q2,p1);
    int o4 = orientation(p2,q2,q1);

    if (o1!=o2 && o3!=o4)
        return  true;

  //Special Cases
    if (o1==0 && onSegment(p1,p2,q1))
        return true;
    if (o2==0 && onSegment(p1,q2,q1))
        return true;
    if (o3==0 && onSegment(p2,p1,q2))
        return true;
    if (o4==0 && onSegment(p2,q1,q2))
        return true;

    else
        return false;
}


int main()
{   struct vertex p1 = {1, 1}, q1 = {10, 1};
    struct vertex p2 = {1, 2}, q2 = {10, 2};

    doIntersect(p1, q1, p2, q2)? cout << "Yes\n": cout << "No\n";

    p1 = {10, 0}, q1 = {0, 10};
    p2 = {0, 0}, q2 = {10, 10};
    doIntersect(p1, q1, p2, q2)? cout << "Yes\n": cout << "No\n";

    p1 = {-5, -5}, q1 = {0, 0};
    p2 = {1, 1}, q2 = {10, 10};
    doIntersect(p1, q1, p2, q2)? cout << "Yes\n": cout << "No\n";

    return 0;
}
