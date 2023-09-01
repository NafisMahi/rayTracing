#include <windows.h> // for MS Windows
#include <GL/glut.h> // GLUT, include glu.h and gl.h
#include <cmath>
#include <iostream>
#include <vector>

void drawAxes()
{
    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0); // Red
    // X axis
    glVertex3f(0, 0, 0);
    glVertex3f(1, 0, 0);

    glColor3f(0, 1, 0); // Green
    // Y axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 1, 0);

    glColor3f(0, 0, 1); // Blue
    // Z axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 1);
    glEnd();
}

class Vector3D
{
public:
    double x;
    double y;
    double z;

    Vector3D()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector3D operator+(const Vector3D &v)
    {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    Vector3D operator-(const Vector3D &v)
    {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    Vector3D operator*(const double d)
    {
        return Vector3D(x * d, y * d, z * d);
    }

    Vector3D operator/(const double d)
    {
        return Vector3D(x / d, y / d, z / d);
    }
};

Vector3D normalize(Vector3D v)
{
    double magnitude = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x /= magnitude;
    v.y /= magnitude;
    v.z /= magnitude;
    return v;
}

Vector3D cross(Vector3D a, Vector3D b)
{
    Vector3D c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

double dot(Vector3D a, Vector3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3D reflectedRay(Vector3D incident, Vector3D normal)
{
    // Assume incident and normal are normalized vectors.
    double dotProduct = dot(incident, normal);

    // Calculate the reflected ray using the formula.
    Vector3D reflected = incident - normal * (2.0 * dotProduct);

    return normalize(reflected); // Return the normalized reflected ray
}

bool isPointInTriangle(Vector3D p, Vector3D a, Vector3D b, Vector3D c)
{
    // Vectors
    Vector3D ab = b - a;
    Vector3D ac = c - a;
    Vector3D ap = p - a;

    // Compute barycentric coordinates
    // double d = ab.dot(ab) * ac.dot(ac) - ab.dot(ac) * ab.dot(ac);
    double d = dot(ab, ab) * dot(ac, ac) - dot(ab, ac) * dot(ab, ac);
    double invD = 1.0 / d;

    // double u = (ac.dot(ac) * ap.dot(ab) - ab.dot(ac) * ap.dot(ac)) * invD;
    double u = (dot(ac, ac) * dot(ap, ab) - dot(ab, ac) * dot(ap, ac)) * invD;
    // double v = (ab.dot(ab) * ap.dot(ac) - ab.dot(ac) * ap.dot(ab)) * invD;
    double v = (dot(ab, ab) * dot(ap, ac) - dot(ab, ac) * dot(ap, ab)) * invD;

    // Check if the point is inside the triangle
    if (u >= 0 && v >= 0 && (u + v) <= 1)
    {
        return true;
    }

    return false;
}

class Ray
{
public:
    Vector3D start, dir;
    Ray(Vector3D start, Vector3D dir)
    {
        this->start = start;
        this->dir = normalize(dir);
    }
};

class Color
{
public:
    double r, g, b;
    Color() {}
    Color(double r, double g, double b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    Color operator*(const double d)
    {
        return Color(r * d, g * d, b * d);
    }
};

double raytriangleIntersection(Ray ray, Vector3D p0, Vector3D p1, Vector3D p2)
{
    Vector3D e1 = p1 - p0;
    Vector3D e2 = p2 - p0;
    Vector3D h = cross(ray.dir, e2);
    double a = dot(e1, h);
    if (a > -0.00001 && a < 0.00001)
    {
        return -1;
    }
    double f = 1 / a;
    Vector3D s = ray.start - p0;
    double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0)
    {
        return -1;
    }
    Vector3D q = cross(s, e1);
    double v = f * dot(ray.dir, q);
    if (v < 0.0 || u + v > 1.0)
    {
        return -1;
    }
    double t = f * dot(e2, q);
    if (t > 0.00001)
    {
        return t;
    }
    else
    {
        return -1;
    }
}

double rayQuadIntersection(Ray ray, Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3)
{
    double t1 = raytriangleIntersection(ray, p0, p1, p2);
    double t2 = raytriangleIntersection(ray, p0, p2, p3);
    if (t1 == -1 && t2 == -1)
    {
        return -1;
    }
    else if (t1 == -1)
    {
        return t2;
    }
    else if (t2 == -1)
    {
        return t1;
    }
    else
    {
        return t1 < t2 ? t1 : t2;
    }
}

double rayFloorIntersection(Ray ray)
{
    // Function to find the intersection point of a ray and the checkerboard
    // Define the plane equation for the infinite checkerboard
    // In this case, it is z = 0, which can be represented as (0, 0, 1) . (x, y, z) = 0

    Vector3D color1 = Vector3D(0, 0, 0); // Black
    Vector3D color2 = Vector3D(1, 1, 1); // White

    Vector3D planeNormal = Vector3D(0, 0, 1);
    double planeD = 0;

    // Calculate the intersection using dot products
    double t = -(dot(planeNormal, ray.start) + planeD) / dot(planeNormal, ray.dir);

    // If t is negative, the ray does not intersect the plane
    if (t < 0)
    {
        return -1; // Or any other way to indicate no intersection
    }

    // Calculate the intersection point using the ray equation
    Vector3D intersectionPoint = ray.start + ray.dir * t;

    return t;
}

double rayfloorintersection(Ray ray, Vector3D reference_point)
{
    Vector3D normal = Vector3D(0, 0, 1);
        double dotP = dot(normal, ray.dir);
        
        if (round(dotP * 100) == 0)
			return -1;

        double t = -(dot(normal,ray.start)) / dotP;

        Vector3D p = ray.start + ray.dir * t;

        if(p.x <= reference_point.x || p.x >= abs(reference_point.x) && p.y <= reference_point.y && p.y >= abs(reference_point.y)){
            return -1;
        }
        
        return t;
    }


void drawInfiniteCheckerBoard()
{
      // draw in infinite checker board with given size and coeff
    double width,height;
    width = height = 50;
    for (int i = -100; i < 100; i++)
    {
        for (int j = -100; j < 100; j++)
        {
            if ((i + j) % 2 == 0)
                glColor3f(1, 1, 1);
            else
                glColor3f(0, 0, 0);

            glBegin(GL_QUADS);
            glVertex3f(i * width, j * height, 0);
            glVertex3f((i + 1) * width, j * height, 0);
            glVertex3f((i + 1) * width, (j + 1) * height, 0);
            glVertex3f(i * width, (j + 1) * height, 0);
            glEnd();
        }
    }
}

// Create a base class named Object in the header/src file mentioned in Step 2. You should define separate classes for each object (shape), and all of them should inherit the Object class

// start code here

class Object
{
public:
    Vector3D reference_point;
    double height, width, length;
    Color color;
    double coEfficients[4];
    int shine;

    Object()
    {
        // Initialize members here if needed
    }

    virtual void draw() {}
    virtual double getIntersectingT(Ray ray) { return 0; }
    virtual Vector3D getNormal(Vector3D intersectionPoint) { return Vector3D(); }
    void setColor() {}
    void setShine() {}
    void setCoEfficients() {}
};

class Sphere : public Object
{
public:
    Sphere(Vector3D center, double radius, Color color, double ambient, double diffuse, double specular, double reflection, double shine)
    {
        this->reference_point = center;
        this->length = radius;
        this->width = radius;
        this->height = radius;
        this->color = color;
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = reflection;
        this->shine = shine;
    }

    void draw()
    {
        int slices = 50;
        int stacks = 50;
        glPushMatrix();
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        glColor3f(color.r, color.g, color.b);
        glutSolidSphere(length, slices, stacks);
        glPopMatrix();
    }

    double getIntersectingT(Ray ray)
    {
        ray.start = ray.start - reference_point; // adjust ray start

        double a = 1;
        double b = 2 * dot(ray.dir, ray.start);
        double c = dot(ray.start, ray.start) - (length * length);

        double discriminant = pow(b, 2) - 4 * a * c;
        double t = -1;
        if (discriminant < 0)
        {
            t = -1;
        }
        else
        {

            if (fabs(a) < 1e-5)
            {
                t = -c / b;
                return t;
            }

            double t1 = (-b - sqrt(discriminant)) / (2 * a);
            double t2 = (-b + sqrt(discriminant)) / (2 * a);

            if (t2 < t1)
            {
                double temp = t1;
                t1 = t2;
                t2 = temp;
            }

            if (t1 > 0)
            {
                t = t1;
            }
            else if (t2 > 0)
            {
                t = t2;
            }
            else
            {
                t = -1;
            }
        }

        return t;
    }

    Vector3D getNormal(Vector3D intersectionPoint)
    {
        Vector3D normal = intersectionPoint - reference_point;
        normal = normalize(normal);
        return normal;
    }

    void setColor(double r, double g, double b)
    {
        this->color = Color(r, g, b);
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double reflection)
    {
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = reflection;
    }
};

class Pyramid : public Object
{
public:
    Pyramid(Vector3D lowest_coordinate, double length, double width, double height, Color color, double ambient, double diffuse, double specular, double reflection, double shine)
    {
        this->reference_point = lowest_coordinate;
        this->length = length;
        this->width = width;
        this->height = height;
        this->color = color;
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = reflection;
        this->shine = shine;
    }

    void draw()
    {
        double x = reference_point.x;
        double y = reference_point.y;
        double z = reference_point.z;

        glColor3f(color.r, color.g, color.b);

        // Base (now aligned with XY-plane)
        glBegin(GL_QUADS);
        glVertex3f(x - width / 2, y - width / 2, z);
        glVertex3f(x + width / 2, y - width / 2, z);
        glVertex3f(x + width / 2, y + width / 2, z);
        glVertex3f(x - width / 2, y + width / 2, z);
        glEnd();

        // Sides (triangles)
        glBegin(GL_TRIANGLES);
        glVertex3f(x, y, z + height);
        glVertex3f(x - width / 2, y - width / 2, z);
        glVertex3f(x + width / 2, y - width / 2, z);

        glVertex3f(x, y, z + height);
        glVertex3f(x + width / 2, y - width / 2, z);
        glVertex3f(x + width / 2, y + width / 2, z);

        glVertex3f(x, y, z + height);
        glVertex3f(x + width / 2, y + width / 2, z);
        glVertex3f(x - width / 2, y + width / 2, z);

        glVertex3f(x, y, z + height);
        glVertex3f(x - width / 2, y + width / 2, z);
        glVertex3f(x - width / 2, y - width / 2, z);
        glEnd();
    }

    Vector3D getNormal(Vector3D intersectionPoint)
    {
        double epsilon = 1e-5; // Small value to handle floating-point errors

        Vector3D top = Vector3D(reference_point.x, reference_point.y, reference_point.z + height);
        Vector3D bottomLeft = Vector3D(reference_point.x - width / 2, reference_point.y - width / 2, reference_point.z);
        Vector3D bottomRight = Vector3D(reference_point.x + width / 2, reference_point.y - width / 2, reference_point.z);
        Vector3D topLeft = Vector3D(reference_point.x - width / 2, reference_point.y + width / 2, reference_point.z);
        Vector3D topRight = Vector3D(reference_point.x + width / 2, reference_point.y + width / 2, reference_point.z);

        // Check if intersection is on base
        if (abs(intersectionPoint.z - reference_point.z) < epsilon)
        {
            return Vector3D(0, 0, -1);
        }

        // Else, check which side face is intersected
        Vector3D vertices[4][3] = {
            {top, bottomLeft, bottomRight},
            {top, bottomRight, topRight},
            {top, topRight, topLeft},
            {top, topLeft, bottomLeft},
        };

        for (int i = 0; i < 4; ++i)
        {
            Vector3D A = vertices[i][0];
            Vector3D B = vertices[i][1];
            Vector3D C = vertices[i][2];
            Vector3D normal = normalize(cross(B - A, C - A));

            // Assume you have a function that checks whether intersectionPoint is on triangle ABC
            if (isPointInTriangle(intersectionPoint, A, B, C))
            {
                return normal;
            }
        }

        // Handle error cases or edge cases here
        return Vector3D(0, 0, 0);
    }

    void setColor(double r, double g, double b)
    {
        this->color = Color(r, g, b);
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double reflection)
    {
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = reflection;
    }

    double getIntersectingT(Ray ray)
    {
        // Initialize to maximum double value, as we are looking for the closest intersection
        double closestT = std::numeric_limits<double>::max(); // changed from INT_MAX for consistency

        // Calculate vertices of the base quad for XY plane
        Vector3D bottom_left = Vector3D(reference_point.x - length / 2, reference_point.y - width / 2, reference_point.z);
        Vector3D bottom_right = Vector3D(reference_point.x + length / 2, reference_point.y - width / 2, reference_point.z);
        Vector3D top_left = Vector3D(reference_point.x - length / 2, reference_point.y + width / 2, reference_point.z);
        Vector3D top_right = Vector3D(reference_point.x + length / 2, reference_point.y + width / 2, reference_point.z);

        // Calculate apex of the pyramid
        Vector3D apex = Vector3D(reference_point.x, reference_point.y, reference_point.z + height);

        // Check intersection with the base
        double t = rayQuadIntersection(ray, bottom_left, bottom_right, top_right, top_left);
        if (t > 0)
        {
            closestT = std::min(closestT, t);
        }

        // Check intersections with the 4 triangular faces
        t = raytriangleIntersection(ray, apex, bottom_left, bottom_right);
        if (t > 0)
        {
            closestT = std::min(closestT, t);
        }
        t = raytriangleIntersection(ray, apex, bottom_right, top_right);
        if (t > 0)
        {
            closestT = std::min(closestT, t);
        }
        t = raytriangleIntersection(ray, apex, top_right, top_left);
        if (t > 0)
        {
            closestT = std::min(closestT, t);
        }
        t = raytriangleIntersection(ray, apex, top_left, bottom_left);
        if (t > 0)
        {
            closestT = std::min(closestT, t);
        }

        // If no intersection has been found, return a negative number
        if (closestT == std::numeric_limits<double>::max())
        {
            return -1;
        }
        return closestT;
    }
};

// now create a class named Cube which inherits object class,you are given the bottom left corner of the cube and the length of the side of the cube
// start code here
class Cube : public Object
{
public:
    Cube(Vector3D bottom_left, double length, Color color, double ambient, double diffuse, double specular, double reflection, double shine)
    {
        this->reference_point = bottom_left;
        this->length = length;
        this->width = length;
        this->height = length;
        this->color = color;
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = reflection;
        this->shine = shine;
    }

    // Function to draw a cube
    void drawCube(float x, float y, float z, float sideLength)
    {

        glColor3f(color.r, color.g, color.b);
        glBegin(GL_QUADS);

        // Front face
        glVertex3f(x, y, z);
        glVertex3f(x + sideLength, y, z);
        glVertex3f(x + sideLength, y + sideLength, z);
        glVertex3f(x, y + sideLength, z);

        // Back face
        glVertex3f(x, y, z + sideLength);
        glVertex3f(x + sideLength, y, z + sideLength);
        glVertex3f(x + sideLength, y + sideLength, z + sideLength);
        glVertex3f(x, y + sideLength, z + sideLength);

        // Left face
        glVertex3f(x, y, z);
        glVertex3f(x, y, z + sideLength);
        glVertex3f(x, y + sideLength, z + sideLength);
        glVertex3f(x, y + sideLength, z);

        // Right face
        glVertex3f(x + sideLength, y, z);
        glVertex3f(x + sideLength, y, z + sideLength);
        glVertex3f(x + sideLength, y + sideLength, z + sideLength);
        glVertex3f(x + sideLength, y + sideLength, z);

        // Top face
        glVertex3f(x, y + sideLength, z);
        glVertex3f(x + sideLength, y + sideLength, z);
        glVertex3f(x + sideLength, y + sideLength, z + sideLength);
        glVertex3f(x, y + sideLength, z + sideLength);

        // Bottom face
        glVertex3f(x, y, z);
        glVertex3f(x + sideLength, y, z);
        glVertex3f(x + sideLength, y, z + sideLength);
        glVertex3f(x, y, z + sideLength);

        glEnd();
    }

    void draw()
    {
        drawCube(reference_point.x, reference_point.y, reference_point.z, length);
    }

    Vector3D getNormal(Vector3D intersectionPoint)
    {
        double epsilon = 1e-5; // Small value to handle floating-point errors

        if (abs(intersectionPoint.x - reference_point.x) < epsilon)
            return Vector3D(-1, 0, 0);
        if (abs(intersectionPoint.x - (reference_point.x + length)) < epsilon)
            return Vector3D(1, 0, 0);

        if (abs(intersectionPoint.y - reference_point.y) < epsilon)
            return Vector3D(0, -1, 0);
        if (abs(intersectionPoint.y - (reference_point.y + height)) < epsilon)
            return Vector3D(0, 1, 0);

        if (abs(intersectionPoint.z - reference_point.z) < epsilon)
            return Vector3D(0, 0, -1);
        if (abs(intersectionPoint.z - (reference_point.z + width)) < epsilon)
            return Vector3D(0, 0, 1);
        return Vector3D(0, 0, 0);
    }

    void setColor(double r, double g, double b)
    {
        this->color = Color(r, g, b);
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double reflection)
    {
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = reflection;
    }

    double getIntersectingT(Ray ray)
    {
        double closestT = std::numeric_limits<double>::max();

        Vector3D A = reference_point;
        Vector3D B = reference_point + Vector3D(length, 0, 0);
        Vector3D C = reference_point + Vector3D(length, 0, length);
        Vector3D D = reference_point + Vector3D(0, 0, length);
        Vector3D E = reference_point + Vector3D(0, length, 0);
        Vector3D F = reference_point + Vector3D(length, length, 0);
        Vector3D G = reference_point + Vector3D(length, length, length);
        Vector3D H = reference_point + Vector3D(0, length, length);

        // Check each face of the cube for intersection
        double t;
        t = rayQuadIntersection(ray, A, B, C, D); // Bottom face
        if (t > 0)
            closestT = std::min(closestT, t);

        t = rayQuadIntersection(ray, E, F, G, H); // Top face
        if (t > 0)
            closestT = std::min(closestT, t);

        t = rayQuadIntersection(ray, A, E, H, D); // Left face
        if (t > 0)
            closestT = std::min(closestT, t);

        t = rayQuadIntersection(ray, B, F, G, C); // Right face
        if (t > 0)
            closestT = std::min(closestT, t);

        t = rayQuadIntersection(ray, A, B, F, E); // Front face
        if (t > 0)
            closestT = std::min(closestT, t);

        t = rayQuadIntersection(ray, D, C, G, H); // Back face
        if (t > 0)
            closestT = std::min(closestT, t);

        if (closestT == std::numeric_limits<double>::max())
        {
            return -1;
        }
        return closestT;
    }
};

struct PointLight
{
    Vector3D pos;
    Color color;
    double falloff;

public:
    PointLight(Vector3D pos, Color color, double falloff)
    {
        this->pos = pos;
        this->color = color;
        this->falloff = falloff;
    }

    void draw()
    {
        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(color.r, color.g, color.b);
        glVertex3f(pos.x, pos.y, pos.z);
        glEnd();
    }
};

struct SpotLight
{
    Vector3D pos;
    Vector3D direction;
    Color color;
    double falloff; // this is different from the spotlight
    double cutoff;

public:
    SpotLight(Vector3D pos, Color color, double falloff, Vector3D direction, double cutoff)
    {
        this->pos = pos;
        this->color = color;
        this->falloff = falloff;
        this->direction = direction;
        this->cutoff = cutoff;
    }
    void draw()
    {
        glPointSize(15);
        glBegin(GL_POINTS);
        glColor3f(color.r, color.g, color.b);
        glVertex3f(pos.x, pos.y, pos.z);
        glEnd();
    }
};