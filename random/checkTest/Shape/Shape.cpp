// Shape.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>    
#include <algorithm> 

using namespace std;

struct ElPoint
{
  ElPoint(double _x, double _y) : x(_x), y(_y) { }

  double x;
  double y;
};

class ElShape
{
public:
  ElShape(string s) : m_sName(s) { }
  virtual double dist_origin() const = 0;
  virtual double area() const = 0; 
  virtual ~ElShape() {}
  const string& GetName() const { return m_sName; }

private:
  string m_sName;
};

class ElCircle : public ElShape
{
public:
  ElCircle(string s, ElPoint zCenetr, double r) : ElShape(s), m_zCenter(zCenetr), m_dRadius(r) { }
  virtual ~ElCircle() {}

  double dist_origin() const { 
    return sqrt(pow(m_zCenter.x,2) + pow(m_zCenter.y,2));
  }

  double area() const {
    const double pi = 3.1415926535897;
    return pow(m_dRadius,2)*pi;
  }

private:
    
  ElPoint m_zCenter;
  double  m_dRadius;
};

// Square's width and height are both m_dLength
class ElSquare : public ElShape
{
public:
  ElSquare(string s, ElPoint zCenetr, double dLength) : ElShape(s), m_zCenter(zCenetr),  m_dLength(dLength) { }
  virtual ~ElSquare() {}

  double dist_origin() const {
    return sqrt(pow(m_zCenter.x,2) + pow(m_zCenter.y,2) );
  }

  double area() const  {
    return pow(m_dLength,2);
  }

private:
  ElPoint m_zCenter;
  double  m_dLength;
};

class ElTriangle : public ElShape
{
public:
  ElTriangle(string s, ElPoint zP1, ElPoint zP2, ElPoint zP3) : ElShape(s), m_zP1(zP1), m_zP2(zP2), m_zP3(zP3) { }
  virtual ~ElTriangle() {}

  double dist_origin() const {

    ElPoint m_zCenter(double(m_zP1.x + m_zP2.x + m_zP3.x)/double(3),double(m_zP1.y + m_zP2.y + m_zP3.y)/double(3));

    return sqrt(pow(m_zCenter.x,2) + pow(m_zCenter.y,2) );
  }

  double area() const {
    return (( m_zP1.x*(m_zP2.y-m_zP3.y) +  m_zP2.x*(m_zP3.y-m_zP1.y)+  m_zP3.x*(m_zP1.y-m_zP2.y))/2.0); 
  }

private:
  ElPoint m_zP1;
  ElPoint m_zP2;
  ElPoint m_zP3;
};

bool comp_area(ElShape* lhs, ElShape* rhs)
{
  return lhs->area() < rhs->area();
}

bool comp_dist(ElShape* lhs, ElShape* rhs)
{
  return lhs->dist_origin() < rhs->dist_origin();
}


int main(int argc, char* argv[])
{
  //Samples of shape objects
  /*
    ElCircle A("A", ElPoint(0,0), 10);
    ElCircle B("B", ElPoint(10,0), 12);
    ElSquare C("C", ElPoint(5,5),  8);
    ElSquare D("D", ElPoint(-15,2),  4);
    ElTriangle E("E", ElPoint(0,0), ElPoint(5,0), ElPoint(0,5));
    ElTriangle F("F", ElPoint(-1, -1), ElPoint(4,0), ElPoint(0,3));
  */
  ElShape* objects[6]; 
  objects[0] = new ElCircle("A", ElPoint(0,0), 10);
  objects[1] = new ElCircle("B", ElPoint(10,0), 12);
  objects[2] = new ElSquare("C", ElPoint(5,5),  8);
  objects[3] = new ElSquare("D", ElPoint(-15,2),  4);
  objects[4] = new ElTriangle("E", ElPoint(0,0), ElPoint(5,0), ElPoint(0,5));
  objects[5] = new ElTriangle("F", ElPoint(-1, -1), ElPoint(4,0), ElPoint(0,3));

  cout << "========================================================= " << endl; 
  cout << "Before Soring with respect to area " << endl; 
  cout << endl; 
  for(int i =0; i< 6 ; i++){
    cout << objects[i]->area() << endl; 
  }
  
  sort(objects, objects+6, comp_area);
 
  cout << endl; 
  cout << "========================================================= " << endl; 
  cout << "After Soring with respect to area " << endl; 
  cout << endl; 
  for(int i =0; i< 6 ; i++){
    cout << objects[i]->area() << endl; 
  }

  cout << endl; 
  cout << "========================================================= " << endl; 
  cout << "Before Soring with respect to distance " << endl; 
  cout << endl; 
  for(int i =0; i< 6 ; i++){
    cout << objects[i]->dist_origin() << endl; 
  }
  
  
  sort(objects, objects+6, comp_dist);

  cout << endl; 

  cout << "========================================================= " << endl; 
  cout << "After Soring with respect to distance " << endl; 
  cout << endl; 

  
  for(int i =0; i< 6 ; i++){
    cout << objects[i]->dist_origin() << endl; 
  }
  cout << endl;


  return 0;
}

