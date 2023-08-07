#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include<iomanip>

using namespace std;

#define INF 100000

ofstream myfile ("Result.csv");

struct vertex
{
    double x;
    double y;
};

double distance1(vertex p1, vertex p2)
{
    return (sqrt(pow((p2.x-p1.x),2)+pow((p2.y-p1.y),2)));
}

double dotProduct(vertex p1, vertex p2)
{

    return p2.x * p1.x + p2.y * p1.y;
}

vertex getVector(vertex p1, vertex p2)
{
    vertex P;
    P.x = p1.x-p2.x;
    P.y = p1.y-p2.y;
    return P;
}

double findAngle(vertex p1, vertex p2)
{
    return atan2(p2.y-p1.y, p2.x-p1.x);
}


// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(double angle[], int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temp arrays
    double L[n1], R[n2];

    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        L[i] = angle[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = angle[m + 1 + j];

    // Merge the temp arrays back into arr[l..r]

    // Initial index of first subarray
    int i = 0;

    // Initial index of second subarray
    int j = 0;

    // Initial index of merged subarray
    int k = l;

    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            angle[k] = L[i];
            i++;
        }
        else
        {
            angle[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of
    // L[], if there are any
    while (i < n1)
    {
        angle[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of
    // R[], if there are any
    while (j < n2)
    {
        angle[k] = R[j];
        j++;
        k++;
    }
}

// l is for left index and r is
// right index of the sub-array
// of arr to be sorted */
void mergeSort(double angle[],int l,int r)
{
    if(l>=r)
    {
        return;//returns recursively
    }
    int m =l+ (r-l)/2;
    mergeSort(angle,l,m);
    mergeSort(angle,m+1,r);
    merge(angle,l,m,r);
}

// UTILITY FUNCTIONS
// Function to print an array
void printArray(double angle[], int size)
{
    for (int i = 0; i < size; i++)
        cout << i+1<<") "<<angle[i] <<endl;
}

double* anglePolygon (vector<vertex>vertices, vertex p, int n)  //Function to calculate all angles of polygon
{
    static double angleArray[] {0.0,0.0,0.0,0.0,0.0,0.0}; //Array to store angle between vertex and horizontal line passing through G

    for (int i=0; i<n; i++)
    {
        angleArray[i] = findAngle (p,vertices[i])* 180 / 3.141592;
        if (angleArray[i]<0)
            angleArray[i] = 360+angleArray[i];

    }

// printArray(angleArray,n);
    return angleArray;
}

double* sortAnglePolygon (double* angleArray, int n)  //Function to sort all angles of polygon
{
    static double angle_Array[10];
    for (int i=0; i<n; i++)
    {
        angle_Array[i] = angleArray[i];
    }

    mergeSort(angle_Array, 0, n-1);

    //  cout << "\n\nSorted array is \n";
    //  printArray(angle_Array, n);

    return angle_Array;
}

int* indexPolygon (double* angleArray, double* sortedAngleArray, int n)  //Function to return index of polygon
{

    static int polygonIndex[] {0,0,0,0,0,0};

    for (int i =0; i<n; i++)
    {
        int j=0;
        if (sortedAngleArray[i] == angleArray[i])
            polygonIndex[i]=i;
        else
        {
            while (sortedAngleArray[i] != angleArray[j])
                j++;
            polygonIndex[i]=j;
        }

    }

//cout << "\n\nPolygon Index is \n";

// for (int i =0; i<n; i++)

//    cout<<i+1<<") "<<polygonIndex[i]+1<<endl;

    return polygonIndex;

}

double minDistance (vertex A, vertex B, vertex P)
{


    vertex AB;
    vertex BP;
    vertex AP;

    AB = getVector(A,B);
    BP = getVector(B,P);
    AP = getVector(A,P);

    // Variables to store dot product
    double AB_BP, AB_AP;

    AB_BP = dotProduct(AB,BP);
    AB_AP = dotProduct(AB,AP);

    double minDistance;

    // Case 1
    if (AB_BP > 0)

        minDistance = distance1 (B,P);


    // Case 2
    else if (AB_AP < 0)

        minDistance = distance1 (A,P);

    // Case 3
    else
    {
        double mod = distance1 (A,B);
        minDistance = abs(AB.x * AP.y - AB.y * AP.x) / mod;
    }

    return minDistance;
}


double findAreaTriangle(vertex v1, vertex v2, vertex v3)

{
    double area;
    area = 0.5 * (v1.x *(v2.y-v3.y)+ v2.x * (v3.y - v1.y) + v3.x * (v1.y-v2.y) );
    return area;
}


///Generalized Mean Value Coordinates
double* meanValue(vector <vertex>Vertices, vertex p)
{
//static double MV[6] {0,0,0,0,0,0};
    vector <vertex>vertices;
    vertices = Vertices;

    vertex P;
    P=p;


    int vSize = vertices.size();

    double w[6] {0,0,0,0,0,0};
    double wt[6] {0,0,0,0,0,0};
    double W=0;

    for (int i=0; i<vSize; i++)
    {
        double prod =1; // Second Big term in MV Coordinates

        int prev = (i + vSize - 1) % vSize;
        int next = (i + 1) % vSize;

        vertex vPrev = vertices[prev];
        vertex vNext = vertices[next];
        vertex vi = vertices[i];

        double rPrev = distance1(vPrev, p);
        double rNext = distance1(vNext, p);

        vertex dPrev = getVector(vPrev,p);
        vertex dNext = getVector(vNext,p);

        wt[i] = sqrt(rPrev * rNext - dotProduct (dPrev,dNext));
        //  cout<<wt[i]<<endl;

        for(int j=0; j<vSize; j++)
        {
            double prod1 =1.0;
            if (j!=prev && j!=i)
            {

                int nextJ = (j + 1) % vSize;

                vertex vNextJ = vertices[nextJ];
                vertex viJ = vertices[j];

                double r = distance1(viJ,p);
                double r1= distance1(vNextJ,p);
                vertex d= getVector(viJ,p);
                vertex dNext= getVector(vNextJ,p);

                double dotJ= dotProduct(d,dNext);

                prod1 = sqrt(r*r1 + dotJ);

                // cout<<r<<"  "<<r1<<"  "<<dotJ<<"   "<<prod1<<endl;
            }
            prod = prod*prod1;
        }

        w[i] = wt[i]* prod;
//   cout<<w[i]<<endl;
        W = W + w[i];
    }

// ofstream myfile ("MV.csv");

    if (myfile.is_open())

    {
        myfile <<"Distance from G:"<<endl;
        double dis;
        for (int i=0; i< vSize; i++)
        {
            dis = distance1(P, vertices[i]);
            myfile<<dis<<endl;
        }

        myfile <<"\n\nMV Coordinates are:"<<endl;
// cout<<"MV Coordinates are:"<<endl;
        static double MV[6] {0,0,0,0,0,0};

        double sum =0;

        for (int i=0; i< vSize; i++)
        {
            MV[i] = w[i]/W;
            //    cout<<"MV "<<i+1<<":  "<<MV[i]*100<<endl;

            //  myfile<<"MV "<<i+1<<":  "<<MV[i]*100<<endl;
            myfile<<MV[i]*100<<endl;

            sum = sum + MV[i]*100;

        }
        return MV;
//cout<<sum;
        myfile<<sum<<endl<<endl;

    }

    //  myfile.close();
}


/// BaryCentric Coordinates
void BarycentricNew(vector <vertex>Vertices, vertex p )
{

    vector <vertex>vertices;
    vertices = Vertices;

    vertex P;
    P=p;

    int vSize = vertices.size();


    double bc[6]= {0,0,0,0,0,0}; //Barycentric Coordinates
    double w[6] = {0,0,0,0,0,0};
    double wTotal =0;

    for (int i=0; i<vSize ; i++)
    {
        int prev = (i + vSize - 1) % vSize;
        int next = (i + 1) % vSize;

        vertex vPrev = vertices[prev];
        vertex vNext = vertices[next];
        vertex vi = vertices[i];

        double A1 = findAreaTriangle(vPrev, vi, vNext);
        double A2 = findAreaTriangle(p, vPrev, vi);
        double A3 = findAreaTriangle(p, vi, vNext);

        w[i] = A1/(A2*A3);
        wTotal = wTotal+w[i];

    }


//  ofstream myfile ("Bary.csv");

    if (myfile.is_open())

    {

        myfile <<"\nBaryCentric Coordinates are:"<<endl;
        cout<<"\n\nBaryCentric Coordinates are:"<<endl;

        double sum =0;
        for (int i=0; i<vSize ; i++)
        {

            bc[i]= w[i]/wTotal;
            cout<<"BC "<<i+1<<":  "<<bc[i]*100<<endl;

            //  myfile<<"BaryCentric Coordinate "<<i+1<<":  "<<bc[i]*100<<endl;
            myfile<<bc[i]*100<<endl;

            sum = sum + bc[i]*100;

        }

        cout<<sum;
        myfile<<sum;

    }
    // myfile.close();
}


///Generalized BaryCentric Coordinates

void genBarycentric(vector <vertex>Vertices, vertex p)
{

    vector <vertex>vertices;
    vertices = Vertices;

    vertex P;
    P=p;


    int vSize = vertices.size();

    double w[6] {0,0,0,0,0,0};
    double W=0;


    for (int i=0; i<vSize; i++)
    {
        double prod =1;
        int prev = (i + vSize - 1) % vSize;
        int next = (i + 1) % vSize;

        vertex vPrev = vertices[prev];
        vertex vNext = vertices[next];
        vertex vi = vertices[i];
        double areaPart1 = findAreaTriangle(vPrev, vi,vNext);

        for(int j=0; j<vSize; j++)
        {
            double prod1 =1.0;
            if (j!=prev && j!=i)
            {

                int nextJ = (j + 1) % vSize;

                vertex vNextJ = vertices[nextJ];
                vertex viJ = vertices[j];
                prod1 = findAreaTriangle (p,viJ,vNextJ);
            }
            prod = prod*prod1;
        }

        w[i] = areaPart1* prod;
//   cout<<w[i]<<endl;
        W = W + w[i];

    }

// ofstream myfile ("GBary.csv");

    if (myfile.is_open())

    {

        myfile <<"\n\nGeneralized BaryCentric Coordinates are:"<<endl;
        cout<<"\n\nGeneralized BaryCentric Coordinates are:"<<endl;
        double B[6] {0,0,0,0,0,0};

        double sum =0;

        for (int i=0; i< vSize; i++)
        {
            B[i] = w[i]/W;
            cout<<"Gen BC "<<i+1<<":  "<<B[i]*100<<endl;

            //myfile<<"BaryCentric Coordinate "<<i+1<<":  "<<B[i]*100<<endl;
            myfile<<B[i]*100<<endl;

            sum = sum + B[i]*100;

        }

        cout<<sum;
        myfile<<sum<<endl<<endl;

    }
    myfile.close();

}

///Checking whether point lies inside or outside polygon

// Given three colinear points p, q, r, the function checks if point q lies on line segment 'pr'

bool onSegment(vertex p, vertex q, vertex r)

{
    if ((q.x<=max(p.x,r.x) && (q.x>=min(p.x,r.x)))
            &&
            (q.y<=max(p.y,r.y) && (q.y>=min(p.y,r.y)))
       )
        return true;

    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(vertex p, vertex q, vertex r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0; // colinear
    return (val > 0)? 1: 2; // clock or counterclock wise
}

// The function that returns true if line segment 'p1q1' and 'p2q2' intersect.
bool doIntersect(vertex p1, vertex q1, vertex p2, vertex q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

// Returns true if the point p lies inside the polygon[] with n vertices

double isInside(vector <vertex>polygon, int n, vertex p)

{

    double distance [6] {0,0,0,0,0,0};

    // There must be at least 3 vertices in polygon[]
    if (n < 3) return false;

    // Create a point for line segment from p to infinite

    vertex extreme;   ///line from p extends to infinite in x direction
    extreme.x = INF;
    extreme.y = p.y;

    //Calculating shortest distance of point P from each line segment of polygon

    int j = 0;
    do
    {
        int next = (j+1)%n;
        distance [j] = minDistance (polygon[j], polygon[next], p);
        j = next;

    }
    while (j!=0);

    cout<<"Minimum Distance of point P from all edges are\n";
    for (int i=0; i<n ; i++)
    {
        cout<<distance[i]<<endl;
    }

    double minDistance = distance[0];

    for (int i=1; i<n; i++)
    {
        if (distance[i]<minDistance)
            minDistance=distance[i];
    }

    cout<<"Minimum among distance is: "<<minDistance<<endl;


    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;

    do
    {
        int next = (i+1)%n;

        // Check if the line segment from 'p' to 'extreme' intersects with the line segment from 'polygon[i]' to 'polygon[next]'

        if (doIntersect(polygon[i], polygon[next], p, extreme))
        {
            // If the point 'p' is colinear with line segment i & 'i-next',
            // then check if it lies on segment. If it lies, return true, otherwise false

            if ( orientation(polygon [i],p, polygon[next] ) == 0)
                return onSegment(polygon[i], p, polygon[next]);

            count++;

        }

        i = next;

    }
    while (i!=0);

    // Return true if count is odd, false otherwise


    if (count%2 ==1)
        return true;
    else
        return false;



}




int main()
{
    vector <vertex>vertices;

    vertex v1;
    v1.x = 3.0;
    v1.y = 8.0;
    vertices.push_back(v1);


    vertex v2;
    v2.x = 1.0;
    v2.y = 6.0;
    vertices.push_back(v2);

    vertex v3;
    v3.x = 2.0;
    v3.y = 3.0;
    vertices.push_back(v3);

    vertex v4;
    v4.x = 5.0;
    v4.y = 3.0;
    vertices.push_back(v4);


    vertex v5;
    v5.x = 6.0;
    v5.y = 6.0;
    vertices.push_back(v5);

    vertex v6;
    v6.x = 5.0;
    v6.y = 8.0;
    vertices.push_back(v6);



    vertex p;
    p.x = 3.0;//3.0;//3.0;   //2.0
    p.y = 1.0;//7.0;//1.0;   //7.0

    int n= vertices.size();

    cout<<fixed << setprecision(1);

    double* angleArray = anglePolygon (vertices, p, n);

//cout << "\n\nAngles are \n";
//printArray(angleArray,n);



    double* sortedAngleArray = sortAnglePolygon (angleArray,n);

//cout << "\n\nSorted Angles are \n";
//printArray(sortedAngleArray,n);


    int* polygonIndex = indexPolygon (angleArray, sortedAngleArray, n);



//cout << "\n\nPolygon Index is \n";
//for (int i =0; i<n; i++)
//   cout<<i+1<<") "<<polygonIndex[i]+1<<endl;

//cout << "\n\nPolygon Index\tAngle\tPolygon Index Post Sorting\tSorted Angle\n";




    vector <vertex>verticesSorted;

    for (int i=0; i<n; i++)
    {
        verticesSorted.push_back (vertices[polygonIndex[i]]);
    }

    for (int i=0; i<n;i++)
        cout<<verticesSorted[i].x<<"   "<<verticesSorted[i].y<<endl;

//if (isInside(vertices, n, p))




        cout << "\n\nPolygon \tAngle\t    Polygon Index \tSorted Angle\n";
        cout << "Index\t     \t            Post Sorting\t\n";

        for (int i=0; i<n; i++)
          {
            cout<<i+1<<")\t\t"<<angleArray[i]<<"\t\t"<<polygonIndex[i]+1<<") \t\t   "<<sortedAngleArray[i]<<endl;
          }

//cout<<endl<<endl;

        cout<<endl<<endl;

        //cout<<"Index\tMV Coordinates\t\tMV Coordinates Post Sorting\tSorted Index:"<<endl;
     //   cout<<"Index\t       MV Coordinates\t\tMV Coordinates       Sorted Index:"<<endl;
     //   cout<<"                                        Post Sorting"<<endl;
        cout<<fixed << setprecision(2);
        double* mvOriginal =  meanValue(vertices,p);
        cout<<"Index and Original MV Coordinates\n";
        for (int i=0; i<n; i++)
            cout<<i+1<<"   "<<mvOriginal[i]*100<<endl;
        cout<<endl<<endl;
        cout<<"Sorted Index & MV Coordinates after Sorting\n";
        double* mvPostSorting = meanValue(verticesSorted,p);
        for (int i=0; i<n; i++)
            cout<<polygonIndex[i]+1<<"    "<<mvPostSorting[i]*100<<endl;



      /*  for (int i=0; i<n; i++)
        {
           cout<<i+1<<")\t\t  "<<mvOriginal[i]*100<<"\t\t\t  "<<mvPostSorting[i]*100<<"\t\t\t"<<polygonIndex[i]+1<<endl;

        }

       */

    // BarycentricNew(vertices,p);
    // genBarycentric(vertices,p);

//  else
//cout<<"Point outside";

    return 0;

}
