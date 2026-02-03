#ifndef GEOMETRY_h
#define GEOMETRY_h
#include <string>

class Point{
    public:
        double x;
        double y;
        //defualt constructor and constructor
        Point(): x(0), y(0) {}
        Point(double x, double y): x(x), y(y) {}
        std::string display() const;

};

class Line{
    //will define the line with a point and a direction vector L(t) = O + dir(t)
public:
    Line(Point origin, double dx, double dy);
    Point project_point_on_line(const Point& point) const;
    double projection_error(const Point& p1, const Point& p2) const;
    std::string display();
    Point O; // the reference point on the line
    Point dir; //normalised direction vector
};

double euclidean_distance(const Point&, const Point&);
#endif
