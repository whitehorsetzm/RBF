#include <iostream>
#include <math.h>
#include "dataclass.h"
#include <map>
#include <vector>
#include "geomobjects.hpp"
using namespace std;

class Vector
{
public:
    Vector() { x = y = z = 0.; }
    Vector(double v1, double v2, double v3) {
        x = v1; y = v2; z = v3;
    }
    Vector(const Vector& v) {
        x = v.x;  y = v.y; z = v.z;
    }

    /* 矢量求模 calc. magnitude */
    double magnitude() const;

    /* 归一 normalization */
    void normalize();

    //calculate the distance
    double getDistance(const Vector &a);

    /* **********************************
     * 重载操作符 overloaded operator
     * **********************************/

    /* 赋值 assignment */
    const Vector& operator = (const Vector& v);
    const Vector& operator = (double v);

    /* 四则运算 simple mathamatic operator */
    const Vector operator + (const Vector& v) const;
    const Vector operator - (const Vector& v) const;
    const Vector operator * (double s) const;
    const Vector operator / (double s) const;

    const Vector& operator += (const Vector& v);
    const Vector& operator -= (const Vector& v);
    const Vector& operator += (double delta);
    const Vector& operator -= (double delta);

    bool operator==(const Vector& v) const;
    bool operator!=(const Vector& v) const;

    /* 点积 & 叉积 dot & cross */
    double operator * (const Vector& v) const;
    const Vector operator ^ (const Vector& V) const;
    const Vector operator - ();

//	friend const Vector operator - (const Vector& p);
    friend const Vector operator * (double scale, const Vector& v);

public:
    double x, y, z;
};


class Box
{
  public:
  Vector pmin, pmax;
  vector<Vector> point;
  set<int> node;
  Box () { ; }

  Box ( const Vector & p1, const Vector & p2)
  {

          pmin.x = min(p1.x, p2.x);
          pmax.x = max(p1.x, p2.x);

          pmin.y = min(p1.y, p2.y);
          pmax.y = max(p1.y, p2.y);

          pmin.z = min(p1.z, p2.z);
          pmax.z = max(p1.z, p2.z);

          point.push_back(pmax);

          Vector temp = pmax;
          temp.x-(pmax.x-pmin.x);
          point.push_back(temp);

          temp = pmax;
          temp.y-(pmax.y-pmin.y);
          point.push_back(temp);

          temp = pmax;
          temp.z-(pmax.z-pmin.z);
          point.push_back(temp);

          temp = pmax;
          temp.x-(pmax.x-pmin.x);
          temp.y-(pmax.y-pmin.y);
          point.push_back(temp);

          temp = pmax;
          temp.x-(pmax.x-pmin.x);
          temp.z-(pmax.z-pmin.z);
          point.push_back(temp);


          temp = pmax;
          temp.z-(pmax.z-pmin.z);
          temp.y-(pmax.y-pmin.y);
          point.push_back(temp);


          point.push_back(pmin);

  }

  const Vector & PMin () const
  {
      return pmin;
  }
  const Vector & PMax () const
  {
      return pmax;
  }

  void Set (const Vector & p)
  {
      pmin = pmax = p;
  }

  void Add (const Vector & p)
  {
      for (int i = 0; i < D; i++)
      {
          if (p(i) < pmin(i)) pmin(i) = p(i);
          else if (p(i) > pmax(i)) pmax(i) = p(i);
      }
  }

  double & operator[] (int i)
  {
      if(i<=D)
          return pmin[i];
      else
          return pmax[i-D];
  }
  const	double & operator[] (int i) const
  {
      if(i<=D)
          return pmin[i];
      else
          return pmax[i-D];
  }
/*
  template <typename T1, typename T2>
  void Set (const IndirectArray<T1, T2> & points)
  {
    Set (points[points.Begin()]);
    for (int i = points.Begin()+1; i < points.End(); i++)
      Add (points[i]);
  }

  template <typename T1, typename T2>
  void Add (const IndirectArray<T1, T2> & points)
  {
    for (int i = points.Begin(); i < points.End(); i++)
      Add (points[i]);
  }
*/
  //盒子的中心
  Vector Center () const
  {
      Vector c;
      for (int i = 0; i < D; i++)
          c(i) = 0.5 * (pmin(i)+pmax(i));
      return c;
  }
  //对角线的长度
  double Diam () const
  {
      return Abs (pmax-pmin);
  }
  //获取第nr个点
  Vector GetPointNr (int nr) const
  {
      Vector p;
      for (int i = 0; i < D; i++)
      {
          p(i) = (nr & 1) ? pmax(i) : pmin(i);
          nr >>= 1;
      }
      return p;
  }

  bool Intersect (const Box<D> & box2) const
  {
      for (int i = 0; i < D; i++)
          if (pmin(i) > box2.pmax(i) ||
              pmax(i) < box2.pmin(i)) return 0;
      return 1;
  }
  //判断点是否在盒子内
  bool IsIn (const Vector & p) const
  {
      for (int i = 0; i < D; i++)
          if (p(i) < pmin(i) || p(i) > pmax(i)) return 0;
      return 1;
  }

};



int init_box(int &x,int &y,int &z,mesh a){
    x[0]=max;
    x[1]=min;
    y[0]=max;
    y[1]=min;
    z[0]=max;
    z[1]=min;
    for(int i=0;i<a.pointernum;++i){
        if(a.point[i].x<x[0])
            x[0]=a.point[i].x;
        else if(a.point[i].x>x[1])
            x[1]=a.point[i].x;
        if(a.point[i].y<y[0])
            y[0]=a.point[i].y;
        else if(a.point[i].y>y[1])
            y[1]=a.point[i].y;
        if(a.point[i].z<z[0])
            z[0]=a.point[i].z;
        else if(a.point[i].z>z[1])
            z[1]=a.point[i].z;
    }
}

double get_minor(double a,double b){
    if(a>b)
        return b;
    else
        return a;
}

double get_major(double a,double b){
    if(a<b)
        return b;
    else
        return a;
}


//已知3点坐标，求平面ax+by+cz+d=0;

void get_panel(Vector p1,Vector p2,Vector p3,double &a,double &b,double &c,double &d)

{

    a = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y);

    b = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z);

    c = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);

    d = 0 - (a * p1.x + b*p1.y + c*p1.z);

}


// 已知三点坐标，求法向量

Vector get_Normal(Vector p1,Vector p2,Vector p3)

{

    a = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y);

    b = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z);

    c = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);

    return Vector(a, b, c);

}


//点到平面距离

double dis_pt2panel(Vector pt,double a,double b,double c,double d)
{

    return f_abs(a * pt.x + b*pt.y + c*pt.z + d) / sqrt(a * a + b * b + c * c);

}




//// <summary>
//// 求一条直线与平面的交点
//// </summary>
//// <param name="planeVector">平面的法线向量，长度为3</param>
//// <param name="planePoint">平面经过的一点坐标，长度为3</param>
//// <param name="lineVector">直线的方向向量，长度为3</param>
//// <param name="linePoint">直线经过的一点坐标，长度为3</param>
//// <returns>返回交点坐标，长度为3</returns>
//Vector returnResult;
// Vector CalPlaneLineIntersectPoint(Vector planeVector, Vector planePoint, Vector lineVector, Vector linePoint)
//{

//float vp1, vp2, vp3, n1, n2, n3, v1, v2, v3, m1, m2, m3, t,vpt;
//vp1 = planeVector.x;
//vp2 = planeVector.y;
//vp3 = planeVector.z;
//n1 = planePoint.x;
//n2 = planePoint.y;
//n3 = planePoint.z;
//v1 = lineVector.x;
//v2 = lineVector.y;
//v3 = lineVector.z;
//m1 = linePoint.x;
//m2 = linePoint.y;
//m3 = linePoint.z;
//vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;


//t = ((n1 - m1) * vp1 + (n2 - m2) * vp2 + (n3 - m3) * vp3) / vpt;
//returnResult[0] = m1 + v1 * t;
//returnResult[1] = m2 + v2 * t;
//returnResult[2] = m3 + v3 * t;


//}



void get_box(HYBRID_MESH mesh,double &x_min,double &y_min,double &z_min,double &x_max,double &y_max,double &z_max){
    double x[2],y[2],z[2];
    x[0]=999999;
    x[1]=-999999;
    y[0]=999999;
    y[1]=-999999;
    z[0]=999999;
    z[1]=-999999;

    for(int i=0;i<mesh.NumNodes;++i){
        if(mesh.nodes[i].coord.x<x[0])
            x[0]=mesh.nodes[i].coord.x;
        else if(mesh.nodes[i].coord.x>x[1])
            x[1]=mesh.nodes[i].coord.x;

        if(mesh.nodes[i].coord.y<y[0])
            y[0]=mesh.nodes[i].coord.y;
        else if(mesh.nodes[i].coord.y>y[1])
            y[1]=mesh.nodes[i].coord.y;

        if(mesh.nodes[i].coord.z<z[0])
            z[0]=mesh.nodes[i].coord.z;
        else if(mesh.nodes[i].coord.z>z[1])
            z[1]=mesh.nodes[i].coord.z;

    }
    x_min=x[0];
    y_min=y[0];
    z_min=z[0];
    x_max=x[1];
    y_max=y[1];
    z_max=z[1];
}



int main()
{
 Box a;



    HYBRID_MESH mesh;
    for(int j=0;j<mesh.NumNodes;++j){
    double d_inf,max_x_s=0,x_i_s;


    set<int> suf;
    set<int> vol;
    for(int i=0;i<mesh.NumTris;++i){         //表面点集合
       suf.insert(mesh.pTris[i].vertices[0]);
       suf.insert(mesh.pTris[i].vertices[1]);
       suf.insert(mesh.pTris[i].vertices[2]);
    }
    for(int i=0;i<mesh.NumQuads;++i){
        suf.insert(mesh.pQuads[i].vertices[0]);
        suf.insert(mesh.pQuads[i].vertices[1]);
        suf.insert(mesh.pQuads[i].vertices[2]);
        suf.insert(mesh.pQuads[i].vertices[3]);
    }

    set<int>::iterator suf_iter=suf.begin();    //内部点集合
    for(int i=0;i<mesh.NumNodes;++i){
        if(i!=*suf_iter)
            vol.insert(i);
        else
            suf_iter++;
    }


    Vector min,max;
    get_box(mesh,min.x,min.y,min.z,max.x,max.y,max.z);

    vector<*Box> box;
    Box *temp = new Box(min,max);
    box.push_back(temp);


    for(int i=0;i<suf.size();++i){    //i为表面点
        x_i_s=sqrt(pow(mesh.nodes[suf[i]].coord_m.x,2)+pow(mesh.nodes[suf[i]].coord_m.y,2)+pow(mesh.nodes[suf[i]].coord_m.z,2));
        if(max_x_s<x_i_s){
            max_x_s=x_i_s;
        }
    }


    d_inf=3*max_x_s;      //公式（9）


    double d_ij,t_ij,d_j_inf;

    double d_min=99999;//d_j到所有表面点中最近的那个；
    double d;
    for(int i=0;i<suf.size();++i){    //i为表面点
        d=sqrt(pow(mesh.nodes[suf[i]].coord.x-mesh.nodes[j].coord.x,2)+pow(mesh.nodes[suf[i]].coord.y-mesh.nodes[j].coord.y,2)+pow(mesh.nodes[suf[i]].coord.z-mesh.nodes[j].coord.z,2));
        if(d<d_min){
            d_min=d;
        }

    }
    d_j_inf=get_minor(d_inf,2*d_min);//公式 (5)

    double p = 1-(10/9)*max(0,min(d_min,d_inf)/d_inf-0.1); //公式（7）
    double n = p*p*(3-2*p); //公式（6）

    double x_v,w,x_s;

    double t_ij;
    double w_ij;
    double w=0,w_x=0;
    double x_move;
    for(int i=0;i<suf.size();++i){    //i为表面点
        d_ij=sqrt(pow(mesh.nodes[suf[i]].coord.x-mesh.nodes[j].coord.x,2)+pow(mesh.nodes[suf[i]].coord.y-mesh.nodes[j].coord.y,2)+pow(mesh.nodes[suf[i]].coord.z-mesh.nodes[j].coord.z,2));
        t_ij=min(d_j_inf,d_ij)/d_j_inf;//公式（4）
        w_ij=pow((1-t_ij),4)*(4*t_ij+1);//公式（3）
        for(int i=0;i<suf.size();++i){
            w+=w_ij;
            w_x+=w_ij*x_move;
        }

    }
    double x_j_v;
    x_j_v=n*w_x/w;     //公式(8)
}
    return 0;
}
