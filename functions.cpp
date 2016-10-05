#include "functions.h"
double sort_increase(double a, double b)
{
    return a<b;
}

vector<int> absSortIndex(vector<double> const& values, bool isIncrease)
{
    vector<int> indices;
    for(size_t i=0;i<values.size();i++)
        indices.push_back(i);
    if(isIncrease)
        sort ( begin(indices), end(indices), [&](int a, int b) { return fabs(values[a])<fabs(values[b]); });
    else
        sort ( begin(indices), end(indices), [&](int a, int b) { return fabs(values[a])>fabs(values[b]); });

    return indices;
}

vector<vector<RealPoint> > getContours(const vector<vector<RealPoint> > &vecPoint)
{
    vector<vector<RealPoint> > vecContour;
    vector<RealPoint> vecC;
    for(auto it:vecPoint)
    {
        for(auto it_bis:it)
        {
            RealPoint p(it_bis[0],it_bis[1]);
            if(find(vecC.begin(),vecC.end(),p)==vecC.end())
                vecC.push_back(p);
        }
        vecContour.push_back(vecC);
        vecC.clear();
    }
    return vecContour;
}

vector<double> readMeanindfulThicknessFile(const char* filename)
{
    vector<double> P;
    ifstream myfile (filename);
    string lineTmp;
    if (myfile.is_open())
    {
        int count=0;
        //idx noiselvl code x y
        double m=0,x=0,y=0;
        std::getline(myfile, lineTmp);//ignore the first two lines
        std::getline(myfile, lineTmp);
        while(myfile >> x >> y >> m)// X Y noiseLevel
        {
            P.push_back(2.0*m);
            count++;
        }
        myfile.close();
        //cout<<count<<" points are read and globalNoise="<<globalNoise<<endl;
        cout<<count<<" points are read "<<endl;
    }
    else cout << "Unable to open file " <<filename;

    return P;
}

void writeFile(const vector<RealPoint>& v, const char* filename)
{
    ofstream myfile (filename);
    if (myfile.is_open())
    {
        for(vector<RealPoint>::const_iterator it=v.begin(); it != v.end(); it++)
            myfile <<(*it)[0]<<" "<<(*it)[1]<<endl;
        myfile.close();
    }
    else cout << "Unable to open file " <<filename;
}

Point getStartPoint(const AlphaThickSegmentComputer2D s)
{
    if(s.size()==0)
        return *(s.begin()); //add point from iterator
    return *(s.containerBegin());
}

Point getEndPoint(const AlphaThickSegmentComputer2D s)
{
    Point p;
    if(s.size()==0)
    {
        for(vector<Point>::const_iterator it=s.begin();it != s.end();it++)
            p=*(it);
    }
    else
    {
        for(vector<Point>::const_iterator it=s.containerBegin();it != s.containerEnd();it++)
            p=*(it);
    }
    return p;
}

RealPoint getStartPoint(const AlphaThickSegmentComputer2DD s)
{
    return *(s.begin());
}

RealPoint getEndPoint(const AlphaThickSegmentComputer2DD s)
{
    RealPoint p;
    for(vector<RealPoint>::const_iterator it=s.begin();it != s.end();it++)
    {
        p=*(it);
    }
    return p;
}

//(Nx,Ny) : normal vector of the line
double getSlope(double Nx,double Ny)
{
    if(fabs(Ny)<1e-6) //vertical line
    {
        if(Nx>0)
            return -1e6; // -INFINITY ?
        else
            return 1e6; //INFINITY ?
    }
    return (-Nx/Ny);
}

double acuteAngle(RealPoint p1, RealPoint p2, RealPoint p3)
{
    RealPoint v1=RealPoint(p1[0]-p2[0],p1[1]-p2[1]);
    RealPoint v2=RealPoint(p3[0]-p2[0],p3[1]-p2[1]);
    double angle=atan2(v2[0],-v2[1]) - atan2(v1[0],-v1[1]);
    if(angle<0)
        angle += 2*M_PI; //report angle to (0,2pi)
    if(angle>M_PI)
        angle=2*M_PI - angle;
    return angle;
}

double signedAngle(RealPoint v1, RealPoint v2)
{
    double Na=sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
    double Nb=sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
    double C=(v1[0]*v2[0] + v1[1]*v2[1])/(Na*Nb);
    double S=(v1[0]*v2[1] - v1[1]*v2[0]);
    if(S<0)
        return(-acos(C));
    else
        return (acos(C));
}

//cos(THETA)=a*b / (||a|| ||b||)
//This gives us just the relative angle between a and b, in [0,PI].
double relativeAngle(RealPoint v1, RealPoint v2)
{
    //normalize the vector
    PointVector<2,double> n_v1, n_v2;
    n_v1[0]=v1[0]/(sqrt(v1[0]*v1[0] + v1[1]*v1[1]));
    n_v1[1]=v1[1]/(sqrt(v1[0]*v1[0] + v1[1]*v1[1]));

    n_v2[0]=v2[0]/(sqrt(v2[0]*v2[0] + v2[1]*v2[1]));
    n_v2[1]=v2[1]/(sqrt(v2[0]*v2[0] + v2[1]*v2[1]));

    //float dotProduct=ax*bx + ay*by;
    //float theta=acos( dotProduct / (vec1Len*vec2Len) );
    return acos(n_v1[0]*n_v2[0]+n_v1[1]*n_v2[1]);
}

// angle between three points (p1,p2,p3)=(b,c,a)
//c*c=a*a + b*b - 2abcos(alpha)
double relativeAngle(RealPoint p1, RealPoint p2, RealPoint p3)
{
    if(p1 == p2 || p1 == p3 || p2==p3)
        return 0;
    double a, b, c, ac;
    a=distancePoints(p1,p2);
    b=distancePoints(p2,p3);
    c=distancePoints(p1,p3);
    ac=(c*c-a*a-b*b);
    if(fabs(ac)<1e-6)
        return M_PI/2;
    if(fabs(c-a-b)<1e-6)
        return M_PI;
    //ac=(c*c-a*a-b*b)/(-2*a*b);
    return acos(ac/(-2*a*b));
}

double distancePointSegment(RealPoint p, double a, double b, double c) //Segment : a x + by + c=0
{
    if(a==0 && b==0)
        return 0;
    return fabs(a*p[0] + b*p[1] + c)/sqrt(a*a+b*b);
}

double distancePointSegment(RealPoint p, RealPoint s1, RealPoint s2)
{
    double a, b, c;
    //line : ax + by + c=0
    //(y1 – y2)x + (x2 – x1)y + (x1y2 – x2y1)=0
    a=s1[1] - s2[1];
    b=s2[0] - s1[0];
    c=s1[0]*s2[1] - s2[0]*s1[1];
    return distancePointSegment(p, a, b, c);
}

double distancePoints(RealPoint p1, RealPoint p2)
{
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]));
}

double distancePointCircle(RealPoint p, RealPoint center, double radius)
{
    return fabs(distancePoints(p,center) - radius);
}

double signedDistancePointCircle(RealPoint p, RealPoint center, double radius)
{
    return distancePoints(p,center) - radius;
}

int findElement(const vector<RealPoint>& vec, RealPoint p, int start)
{
    int it;
    for(it=start; it<vec.size(); it++)
        if(vec.at(it)==p)
            return it;
    //if(it==vec.size())
    for(it=0; it<start; it++)
        if(vec.at(it)==p)
            //return vec.size()+it;//FIXME
            return (int)vec.size()+it;
    return -1;
}

int findElement(const vector<RealPoint>& vec, RealPoint p)
{
    vector<RealPoint>::const_iterator it=find(vec.begin(),vec.end(),p);
    if(it != vec.end())
        return (it-vec.begin());
    else
        return -1;
}

int findStartElement(const vector<RealPoint>& vec, const AlphaThickSegmentComputer2DD s)
{
    RealPoint pStart=getStartPoint(s);
    vector<RealPoint>::const_iterator it=find(vec.begin(),vec.end(),pStart);
    if(it != vec.end())
        return (it-vec.begin());
    return -1;
}

int findEndElement(const vector<RealPoint>& vec, const AlphaThickSegmentComputer2DD s)
{
    RealPoint pStart=getStartPoint(s);
    RealPoint pEnd=getEndPoint(s);
    int indexStart=-1;
    vector<RealPoint>::const_iterator itStart=find(vec.begin(),vec.end(),pStart);
    if(itStart != vec.end())
        indexStart=(itStart-vec.begin());
    int indexEnd=-1;
    vector<RealPoint>::const_iterator itEnd=find(vec.begin(),vec.end(),pEnd);
    if(itEnd != vec.end())
        indexEnd=(itEnd-vec.begin());
    return indexEnd;
}

/******* Evaluation criteria **********/
double error_CR(const vector<RealPoint>& contour, const vector<RealPoint>& DP)
{
    return (double)contour.size()/DP.size();
}

double error_ISE(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed)
{
    double ise=0.0, d=0.0;
    for(size_t it=0; it+1<DP.size(); it++)
    {
        int indexStart=indexDP.at(it);
        int indexEnd=indexDP.at(it+1);
        if(indexStart<0 || indexEnd<0)
            cout<<"Pb in error_ISE index of DP"<<endl;
        for(int index=indexStart+1; index<indexEnd; index++)
        {
            d=distancePointSegment(contour.at(index),DP.at(it), DP.at(it+1));
            ise += d*d;
        }
        if(isClosed && indexStart>indexEnd)
        {
            //cout<<"it="<<it<<" : indexStart="<<indexStart<<", indexEnd="<<indexEnd<<endl;
            for(int index=indexStart+1; index<indexEnd+contour.size(); index++)
            {
                d=distancePointSegment(contour.at(index%contour.size()),DP.at(it), DP.at(it+1));
                ise += d*d;
            }
        }
    }
    //cout<<"ise="<<ise<<endl;
    if(isClosed)//consider the last seg DP.front() and DP.back()
    {
        int indexStart=indexDP.back();
        int indexEnd=indexDP.front();
        if(indexStart>indexEnd)
            indexEnd=indexEnd+contour.size();
        //cout<<"isClosed"<<" : indexStart="<<indexStart<<", indexEnd="<<indexEnd<<endl;
        if(indexStart<0 || indexEnd<0)
            cout<<"Pb in error_ISE index of DP"<<endl;
        for(int index=indexStart+1; index<indexEnd; index++)
        {
            d=distancePointSegment(contour.at(index%contour.size()),DP.back(), DP.front());
            ise += d*d;
        }
    }
    return ise;
}

double error_L_infini(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed)
{
    double d_max=0.0, d;
    for(size_t it=0; it+1<DP.size(); it++)
    {
        int indexStart=indexDP.at(it);
        int indexEnd=indexDP.at(it+1);
        if(indexStart<0 || indexEnd<0)
            cout<<"Pb in error_ISE index of DP"<<endl;
        for(int index=indexStart+1; index<indexEnd; index++)
        {
            d=distancePointSegment(contour.at(index),DP.at(it), DP.at(it+1));
            if(d>d_max)
                d_max=d;
        }
    }
    if(isClosed)//consider the last seg DP.front() and DP.back()
    {
        int indexStart=indexDP.back();
        int indexEnd=indexDP.front()+contour.size();
        //cout<<"it="<<it<<" : indexStart="<<indexStart<<", indexEnd="<<indexEnd<<endl;
        if(indexStart<0 || indexEnd<0)
            cout<<"Pb in error_ISE index of DP"<<endl;
        for(int index=indexStart+1; index<indexEnd; index++)
        {
            d=distancePointSegment(contour.at(index%contour.size()),DP.back(), DP.front());
            if(d>d_max)
                d_max=d;
        }
    }
    return d_max;
}

void error_All(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed)
{
    double cr=error_CR(contour,DP);
    double ise=error_ISE(contour,DP,indexDP,isClosed);
    double dmax=error_L_infini(contour,DP,indexDP,isClosed);
    double fom=cr/ise;
    double fom_M=(cr*cr)/ise;
    double fom_ND=(cr*cr*cr)/(ise*dmax);

    cout<<"cr="<<cr<<", ise="<<ise<<", dmax="<<dmax<<", FOM="<<fom<<", FOM_M="<<fom_M<<", FOM_ND="<<fom_ND<<endl;
}

/******* Evaluation criteria **********/

/******* Arc/segment decomposition ****/
//verify whether p3 is left of p1p2
int isLeft(RealPoint p3, RealPoint p1, RealPoint p2)
{
    double v1, v2, v3, v4;
    v1=p2[0] - p1[0];
    v2=p3[1] - p1[1];
    v3=p3[0] - p1[0];
    v4=p2[1] - p1[1];

    double t=v1*v2 - v3*v4;
    if(t>0) return 1;//left
    else if(t<0) return -1;//right
    else return 0;//on the line
}

RealPoint determineCenter(RealPoint p1, RealPoint p2, RealPoint p3)
{
    double a1, b1, a2, b2, c1, c2;
    double xA, yA, xB, yB, xC, yC;
    xA=p1[0];
    yA=p1[1];
    xB=p2[0];
    yB=p2[1];
    xC=p3[0];
    yC=p3[1];

    a1=xA-xB;
    b1=yA-yB;
    a2=xA-xC;
    b2=yA-yC;
    c1=(xA*xA-xB*xB+yA*yA-yB*yB)/2;
    c2=(xA*xA-xC*xC+yA*yA-yC*yC)/2;
    double x,y,dentaY;
    dentaY=b1*a2-a1*b2;
    if(dentaY!=0){
        y= (double)(c1*a2-a1*c2)/(double)dentaY;
        if (a1!=0) x=(double)(c1-b1*y)/(double)a1;
        else if(a2!=0) x=(double)(c2-b2*y)/(double)a2;
        else {
            cout<<"Error: 3 points of the arc are colinear."<<endl;
            x=-1;
            y=-1;
        }
    }
    else
    {
        x=-1;
        y=-1;
    }
    return RealPoint(x,y);
}

double determineRadius(RealPoint p1, RealPoint p2, RealPoint p3)
{
    RealPoint center=determineCenter(p1,p2,p3);
    return (distancePoints(p1,center) + distancePoints(p2,center) + distancePoints(p3,center))/3.0;
}

double determineRadius(RealPoint center, RealPoint p)
{
    return distancePoints(center,p);
}

RealPoint determineCenter(Point p1, Point p2, Point p3)
{
    double a1, b1, a2, b2, c1, c2;
    double xA, yA, xB, yB, xC, yC;
    xA=p1[0];
    yA=p1[1];
    xB=p2[0];
    yB=p2[1];
    xC=p3[0];
    yC=p3[1];

    a1=xA-xB;
    b1=yA-yB;
    a2=xA-xC;
    b2=yA-yC;
    c1=(xA*xA-xB*xB+yA*yA-yB*yB)/2;
    c2=(xA*xA-xC*xC+yA*yA-yC*yC)/2;
    double x,y,dentaY;
    dentaY=(b1*a2-a1*b2);
    if(dentaY!=0) {
        y= (c1*a2-a1*c2)/dentaY;
        if (a1!=0) x=(c1-b1*y)/a1;
        else if(a2!=0) x=(c2-b2*y)/a2;
        else {
            cout<<"Error: 3 points of the arc are colinear."<<endl;
            x=-1;
            y=-1;
        }
    }
    else {
        x=-c1;
        y=-1;
    }
    return RealPoint(x,y);
}

RealPoint determineCenter(Point p1, Point p2, double radius, bool negatif)
{
    double a, b, c, q;
    double xA, yA, xB, yB, xC, yC;
    xA=p1[0];
    yA=p1[1];
    xB=p2[0];
    yB=p2[1];
    xC=(p1[0]+p2[0])/2.0;
    yC=(p1[1]+p2[1])/2.0;

    a=xB-xA;
    b=yA-yB;
    q=sqrt(a*a+b*b);
    c=sqrt(radius*radius - (q*q)/4.0);

    //find direction from p1 to p2
    double x,y;
    if(negatif){ //p1 is left and p2 is right
        x=xC - c*b/q;
        y=yC - c*a/q;
    }
    else{ //p1 is left and p2 is right
        x=xC + c*b/q;
        y=yC + c*a/q;
    }
    return RealPoint(x,y);
}

double determineRadius(Point p1, Point p2, Point p3)
{
    RealPoint center=determineCenter(p1,p2,p3);
    return (distancePoints(p1,center) + distancePoints(p2,center) + distancePoints(p3,center))/3.0;
}

double determineRadius(RealPoint center, Point p)
{
    return distancePoints(p,center);
}

double arcLength(Point p1, Point p2, Point p3)
{
    double R=determineRadius(p1,p2,p3);
    double angle=2*(M_PI - acuteAngle(p1,p2,p3));//2*M_PI - 2*angle;
    return R*angle;
}
/******* Arc/segment decomposition ****/

/******* Error Arc/Segment appox ******/
double iseContourCircle(const vector<RealPoint>& contour, Point p1, Point p2, RealPoint center, double radius)
{
    double ise=0, d=0;
    int indexP1=findElement(contour,p1);
    int indexP2=findElement(contour,p2);
    if(indexP1 != -1 && indexP2 != -1)
        for(int it=indexP1+1; it<indexP2; it++)
        {
            d=distancePointCircle(contour.at(it),center,radius);
            ise += d*d;
        }
    return ise;
}

double iseContourCircle(const vector<RealPoint>& contour, int indexP1, int indexP2, RealPoint center, double radius)
{
    double ise=0, d=0;
    for(int it=indexP1+1; it<indexP2; it++)
    {
        d=distancePointCircle(contour.at(it),center,radius);
        ise += d*d;
    }
    return ise;
}

double lmaxContourCircle(const vector<RealPoint>& contour, Point p1, Point p2, RealPoint center, double radius)
{
    double lmax=0, l=0;
    int indexP1=findElement(contour,p1);
    int indexP2=findElement(contour,p2);
    if(indexP1 != -1 && indexP2 != -1)
        for(int it=indexP1+1; it<indexP2; it++)
        {
            l=distancePointCircle(contour.at(it),center,radius);
            if(lmax<l)
                lmax=l;
        }
    return lmax;
}

double lmaxContourCircle(const vector<RealPoint>& contour, int indexP1, int indexP2, RealPoint center, double radius)
{
    double lmax=0, l=0;
    for(int it=indexP1+1; it<indexP2; it++)
    {
        l=distancePointCircle(contour.at(it),center,radius);
        if(lmax<l)
            lmax=l;
    }
    return lmax;
}

double iseContourSegment(const vector<RealPoint>& contour, Point p1, Point p2)
{
    double ise=0, d=0;
    int indexP1=findElement(contour,p1);
    int indexP2=findElement(contour,p2);
    if(indexP1 != -1 && indexP2 != -1)
        for(int it=indexP1+1; it<indexP2; it++)
        {
            d=distancePointSegment(contour.at(it),contour.at(indexP1),contour.at(indexP2));
            ise += d*d;
        }
    return ise;
}

double iseContourSegment(const vector<RealPoint>& contour, int indexP1, int indexP2)
{
    double ise=0, d=0;
    for(int it=indexP1+1; it<indexP2; it++)
    {
        d=distancePointSegment(contour.at(it),contour.at(indexP1),contour.at(indexP2));
        ise += d*d;
    }
    return ise;
}

double lmaxContourSegment(const vector<RealPoint>& contour, Point p1, Point p2)
{
    double lmax=0, l=0;
    int indexP1=findElement(contour,p1);
    int indexP2=findElement(contour,p2);
    if(indexP1 != -1 && indexP2 != -1)
        for(int it=indexP1+1; it<indexP2; it++)
        {
            l=distancePointSegment(contour.at(it),contour.at(indexP1),contour.at(indexP2));
            if(lmax<l)
                lmax=l;
        }
    return lmax;
}

double lmaxContourSegment(const vector<RealPoint>& contour, int indexP1, int indexP2)
{
    double lmax=0, l=0;
    for(int it=indexP1+1; it<indexP2; it++)
    {
        l=distancePointSegment(contour.at(it),contour.at(indexP1),contour.at(indexP2));
        if(lmax<l)
            lmax=l;
    }
    return lmax;
}
/******* Error Arc/Segment appox ******/
