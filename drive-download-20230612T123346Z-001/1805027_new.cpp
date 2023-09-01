#include <GL/glut.h> // GLUT, include glu.h and gl.h
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "bitmap_image.hpp"
#include "1805027.h"
using namespace std;

struct point
{
    GLfloat x, y, z;
    point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    point operator+(const point &p)
    {
        return point(x + p.x, y + p.y, z + p.z);
    }

    point operator-(const point &p)
    {
        return point(x - p.x, y - p.y, z - p.z);
    }

    point operator*(const double d)
    {
        return point(x * d, y * d, z * d);
    }
};

Vector3D pos; // position of the eye
Vector3D l;   // look/forward direction
Vector3D r;   // right direction
Vector3D u;   // up direction
Vector3D c;
double near_plane, far_plane, fov_y, fov_x, aspect_ratio;
int recursion_level, pixels;
int screenWidth, screenHeight;
double checkerBoardWidth = 50;
double checkerBoard_width;
Color checkerBoard_color;
vector<Object *> objects;
vector<PointLight *> pointLights;
vector<SpotLight *> spotLights;

float cameraX = 5.0f;
float cameraY = 5.0f;
float cameraZ = 5.0f;
float centerX, centerY, centerZ;
GLfloat reference = 0.0, triangle = 0.0, stepsize = 16.0, stepno = 0, cylinderRadius = 0.0, cylinderHeight = sqrt(2.0), scaleFactor = 16.0, cylinderScalefactor = 0.0, sphereScaleFactor = 0.0, angle = 0.0;

/* Initialize OpenGL Graphics */
void initGL()
{
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}

// Global variables
GLfloat eyex = 4, eyey = 4, eyez = 4;
GLfloat centerx = 0, centery = 0, centerz = 0;
GLfloat upx = 0, upy = 1, upz = 0;
bool isAxes = true, isCube = false, isPyramid = false;

vector<vector<Vector3D>> createPointBuffer()
{
    // Vector3D pointBuffer[screenWidth][screenHeight];
    vector<vector<Vector3D>> pointBuffer(screenWidth, vector<Vector3D>(screenHeight));

    // Calculate the midpoint
    Vector3D midpoint = pos + (l * near_plane);

    // Calculate dimensions of the screen in world space
    double height = 2 * near_plane * tan(fov_y / 2.0 * M_PI / 180.0);
    // double width = height * aspect_ratio;
    double width = 2 * near_plane * tan(((fov_y * aspect_ratio) / 2.0) * M_PI / 180.0);

    double dx = height / screenWidth;
    double dy = width / screenHeight;

    int midIndex = screenWidth / 2;

    Vector3D topLeft = midpoint + u * midIndex * dy - r * midIndex * dx;

    // Calculate the lower-left corner of the screen in world coordinates
    // Vector3D lowerLeft = midpoint - (u * (height / 2.0)) - (r * (width / 2.0));

    // Generate the pointBuffer array
    for (int i = 0; i < screenWidth; ++i)
    {
        for (int j = 0; j < screenHeight; ++j)
        {

            // Interpolation factors
            // double alpha = (double)i / (screenWidth - 1);
            // double beta = (double)j / (screenHeight - 1);

            // Calculate the 3D point that corresponds to this pixel
            Vector3D Point = topLeft + (r * (j * dx)) - (u * (i * dy));

            pointBuffer[i][j] = Point;
        }
    }

    for (int i = 0; i < screenWidth; i++)
    {
        for (int j = 0; j < screenHeight; j++)
        {
            // cout<<pointBuffer[i][j].x<<" "<<pointBuffer[i][j].y<<" "<<pointBuffer[i][j].z<<endl;
        }
    }
    return pointBuffer;
}

bitmap_image createBitmapObject()
{
    bitmap_image image(screenWidth, screenHeight);
    image.set_all_channels(0, 0, 0);
    image.save_image("output.bmp");
    return image;
}

bool isObscuredPointLight(PointLight *pointLight, Vector3D intersectionPoint)
{
    Vector3D direction = normalize(pointLight->pos - intersectionPoint);
    Vector3D source = intersectionPoint + (direction) * 0.0001;
    Ray ray = Ray(source, direction);
    double distance = sqrt(dot((pointLight->pos - intersectionPoint), ((pointLight->pos - intersectionPoint))));
    for(int i =0; i<objects.size(); i++){
        double t = objects[i]->getIntersectingT(ray);
        if(t > 0 && t < distance){
           // cout << " in condition " << endl;
            return true;
        }
    }
    double t = rayFloorIntersection(ray);
    if(t > 0 && t < distance){
        // cout << " in condition of floor" << endl;
        return true;
    }
   // cout << "falseeee" << endl;
    return false;
}

bool isObscuredSpotLight(SpotLight* spotLight, Vector3D intersectionPoint)
{
    Vector3D lightStart;
    Vector3D lightDirection;
    lightDirection = intersectionPoint - lightStart;
    Vector3D direction = normalize(spotLight->pos - intersectionPoint);
    Vector3D source = intersectionPoint + (direction) * 0.0001;

    double angle = acos(dot(lightDirection, normalize(spotLight->direction))) * 180.0 / M_PI;
    if (angle > spotLight->cutoff * 180.0 / M_PI)
    {
        return true;
    }

    Ray ray = Ray(source, direction);
    double distance = sqrt(dot((spotLight->pos - intersectionPoint), ((spotLight->pos - intersectionPoint))));
    for(int i =0; i<objects.size(); i++){
        double t = objects[i]->getIntersectingT(ray);
        if(t > 0 && t < distance){
            return true;
        }
    }
    double t = rayFloorIntersection(ray);
    if(t > 0 && t < distance){
        return true;
    }
    return false;
    //t > 0 and t <distance between intersection point and light source
}

// boolean isObscuredPointLight(PointLight *pointLight, Vector3D intersectionPoint)
// {
//     Vector3D lightStart;
//     Vector3D lightDirection;
//     double t;
//     lightStart = pointLight->pos;
//     lightDirection = intersectionPoint - lightStart;
//     normalize(lightDirection);
//     double t_object = (intersectionPoint - lightStart).x / lightDirection.x;
//     Ray ray = Ray(lightStart, lightDirection);
//     double min_t = std::numeric_limits<double>::max();

//     for (int i = 0; i < objects.size(); ++i)
//     {
//         t = objects[i]->getIntersectingT(ray);
//         if (t > -0.00001 && t < min_t)
//         {
//             min_t = t;
//             if (min_t < t_object)
//             {
//                 return true;
//             }
//         }
//     }
//     return false;
// }

// boolean isObscuredSpotLight(SpotLight *spotLight, Vector3D intersectionPoint)
// {
//     Vector3D lightStart;
//     Vector3D lightDirection;
//     double t;
//     lightStart = spotLight->pos;
//     lightDirection = intersectionPoint - lightStart;
//     normalize(lightDirection);

//     // check 1
//     double angle = acos(dot(lightDirection, normalize(spotLight->direction))) * 180.0 / M_PI;
//     if (angle > spotLight->cutoff * 180.0 / M_PI)
//     {
//         return true;
//     }

//     // check 2
//     double t_object = (intersectionPoint - lightStart).x / lightDirection.x;
//     Ray ray = Ray(lightStart, lightDirection);
//     double min_t = std::numeric_limits<double>::max();

//     for (int i = 0; i < objects.size(); ++i)
//     {
//         t = objects[i]->getIntersectingT(ray);
//         if (t > -0.00001 && t < min_t)
//         {
//             min_t = t;
//             if (min_t < t_object)
//             {
//                 return true;
//             }
//         }
//     }
//     return false;
// }

Color calculateColorCoefficients(Vector3D intersectionPoint, Color intersectionPointColor, Vector3D normal, Object *object,int level)
{
    Color color;
    normal = normalize(normal);
    // Ambient
    if (object == NULL)
    {
        color = intersectionPointColor * checkerBoard_color.r; //ambient of floor * intersection point color
    }

    else
    {
        color = intersectionPointColor * object->coEfficients[0];
    }

    double lambert = 0.0;
    double phong = 0.0;
    for (int i = 0; i < pointLights.size(); i++)
    {
        if (!isObscuredPointLight(pointLights[i], intersectionPoint))
        {
            // cout << " in point light" << endl;
            // cout << "trueeeeee " <<endl;
            Vector3D toSource = pointLights[i]->pos - intersectionPoint;
            toSource = normalize(toSource);
            double distance = sqrt(dot((pointLights[i]->pos - intersectionPoint), ((pointLights[i]->pos - intersectionPoint))));
            double scaling_factor = std::exp(-distance * distance * pointLights[i]->falloff);
            if (dot(normal, toSource) > 0)
            {
                lambert += dot(normal, toSource) * scaling_factor;
            }

            // cout << "dot  " << dot(normal, toSource) << endl;
            Vector3D reflected = reflectedRay(normalize(intersectionPoint - pos), normal);
            if (object != NULL)
            {
                if (dot(reflected, toSource) > 0)
                {
                    phong += pow(dot(reflected, toSource), object->shine ) * scaling_factor;
                }
                // phong += std::pow(dot(reflected, toSource), object->shine) * scaling_factor;
            }
        }

        else{
            // cout << "falseeeeeeeeeeee" <<endl;
        }
    }

    for (int i = 0; i < spotLights.size(); i++)
    {
        if (!isObscuredSpotLight(spotLights[i], intersectionPoint ))
        {
            // cout << " FALSEE " << endl;
            Vector3D toSource = spotLights[i]->pos - intersectionPoint;
            toSource = normalize(toSource);
            double distance = sqrt(dot((spotLights[i]->pos - intersectionPoint), ((spotLights[i]->pos - intersectionPoint))));
            double scaling_factor = std::exp(-distance * distance * spotLights[i]->falloff);
            if (dot(normal, toSource) > 0)
            {
                lambert += dot(normal, toSource) * scaling_factor;
            }
            // cout << "lambert " << lambert << endl;
            Vector3D reflected = reflectedRay(normalize(intersectionPoint - pos), normal);
            if (object != NULL)
            {
                if (dot(reflected, toSource) > 0)
                {
                    phong += pow(dot(reflected, toSource), object->shine ) * scaling_factor;
                }
            }
        }
        else{
            //  cout << "TRUEEE" <<endl;
        }
    }
    
    
    
    if (object == NULL)
    {
        color.r = color.r + (lambert * checkerBoard_color.g) * color.r;
        color.g = color.g + (lambert * checkerBoard_color.g) * color.g;
        color.b = color.b + (lambert * checkerBoard_color.g) * color.b;
        // cout << "color " << color.r << " " << color.g << " " << color.b << endl;
    }
    else
    {
        color.r = color.r + (lambert * object->coEfficients[1] * intersectionPointColor.r) + (phong * object->coEfficients[2] * intersectionPointColor.r);
        color.g = color.g + (lambert * object->coEfficients[1] * intersectionPointColor.g) + (phong * object->coEfficients[2] * intersectionPointColor.g);
        color.b = color.b + (lambert * object->coEfficients[1] * intersectionPointColor.b) + (phong * object->coEfficients[2] * intersectionPointColor.b);
    }

    if (level < recursion_level)
    {
        Vector3D reflected = reflectedRay(normalize(intersectionPoint - pos), normal);
        reflected = normalize(reflected);
        Ray ray = Ray(intersectionPoint + reflected * 0.0001, reflected);
        double t;
        double index = -1;
        double min_t = std::numeric_limits<double>::max();
        Color color_;
        for (int i = 0; i < objects.size(); i++)
        {
            t = objects[i]->getIntersectingT(ray);
            if (t > 0 && t < min_t)
            {
                min_t = t;
                object = objects[i];
                index = i;
                color_ = object->color;
            }
        }
        
        t = rayFloorIntersection(ray);
        if(t < min_t && t > 0)
        {
            int x = static_cast<int>(floor(intersectionPoint.x) / checkerBoardWidth);
            int y = static_cast<int>(floor(intersectionPoint.y) / checkerBoardWidth);
             if ((x + y) % 2 == 0)
                {
                    color_ = Color(1, 1, 1);
                }
                else
                {
                    color_ = Color(0, 0, 0);
                }
            index = 100000;
        }
        if (index != -1)
        {
            Color colortemp(0, 0, 0);
            colortemp = calculateColorCoefficients(intersectionPoint + reflected * 0.0001 + reflected * min_t, color_, reflected, object, level + 1);
            color.r += colortemp.r * object->coEfficients[3];
            color.g += colortemp.g * object->coEfficients[3];
            color.b += colortemp.b * object->coEfficients[3];
        }
    }

    return color;
}

void capture()
{
    vector<vector<Vector3D>> pointBuffer = createPointBuffer();
    bitmap_image image = createBitmapObject();

    // Calculate the color of each pixel
    for (int i = 0; i < screenWidth; ++i)
    {
        for (int j = 0; j < screenHeight; ++j)
        {
            // Calculate the 3D point that corresponds to this pixel
            Vector3D Point = pointBuffer[i][j];
            bool isfloor = false;

            // Calculate the ray direction
            Vector3D rayDir = Point - pos;
            normalize(rayDir);

            // Create the ray
            Ray ray(Point, rayDir);

            // for each object,find the intersecting point and save the least positive value of t
            // start code here
            double t;
            double min_t = 99999;
            Object *nearestObject = NULL;
            for (int k = 0; k < objects.size(); ++k)
            {
                t = objects[k]->getIntersectingT(ray);
                // cout << "t" << t << endl;
                if (t > 0 && t < min_t)
                {
                    min_t = t;
                    nearestObject = objects[k];
                }
            }

            double floor_t = rayFloorIntersection(ray);
            // double floor_t = rayfloorintersection(ray,Vector3D(pos.x,pos.y,pos.z));
            if (floor_t > 0 && floor_t < min_t)
            {
                min_t = floor_t;
                nearestObject = NULL;
                isfloor = true;
            }

            Color color;

            if (nearestObject != NULL)
            {
                color = nearestObject->color;
                Vector3D normal = nearestObject->getNormal(ray.start + ray.dir * min_t);
                if (dot(rayDir, normal) > -0.00001)
                {
                    normal = normal * -1;
                }
                color = calculateColorCoefficients(ray.start + ray.dir * min_t, color, normal, nearestObject,0);
                image.set_pixel(j, i, 255 * color.r, 255 * color.g, 255 * color.b);
            }

            else if (isfloor)
            {
                Vector3D intersectionPoint = ray.start + ray.dir * min_t;
                // Find which square of the checkerboard the intersection point is in
                int x = static_cast<int>(floor(intersectionPoint.x) / checkerBoardWidth);
                int y = static_cast<int>(floor(intersectionPoint.y) / checkerBoardWidth);

                Vector3D normal = Vector3D(0, 0, 1);
                // if (dot(rayDir, normal) > -0.00001)
                // {
                //     normal = normal * -1;
                // }
                // Determine the color based on the square's coordinates

                if ((x + y) % 2 == 0)
                {
                    color = Color(1, 1, 1);
                }
                else
                {
                    color = Color(0, 0, 0);
                }
                color = calculateColorCoefficients(ray.start + ray.dir * min_t, color, normal, NULL,0);
                image.set_pixel(j, i, 255 * color.r, 255 * color.g, 255 * color.b);
            }
            else
            {
                // Set the pixel color to black
                image.set_pixel(j, i, 0, 0, 0);
            }
        }
    }
    image.save_image("output.bmp");
}

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */

void display()
{
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); // To operate on Model-View matrix
    glLoadIdentity();           // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);

    // control viewing (or camera)
    // gluLookAt(eyex,eyey,eyez,
    //           centerx,centery,centerz,
    //           upx,upy,upz);

    gluLookAt(pos.x, pos.y, pos.z,
              pos.x + l.x, pos.y + l.y, pos.z + l.z,
              u.x, u.y, u.z);

    //  drawAxes();
    drawInfiniteCheckerBoard();

    for (int i = 0; i < objects.size(); i++)
    {
        // std ::cout << objects.size() << endl;
        objects[i]->draw();
    }

    for (int i = 0; i < pointLights.size(); i++)
    {
        pointLights[i]->draw();
    }

    for (int i = 0; i < spotLights.size(); i++)
    {
        spotLights[i]->draw();
    }

    // glRotatef(angle,0,1,0);

    glutSwapBuffers(); // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(80.0f, 1, 1, 1000);
}

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y)
{

    float v = 0.25;
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;
    double rate = 0.1;
    c.x = pos.x + l.x;
    c.y = pos.y + l.y;
    c.z = pos.z + l.z;

    switch (key)
    {
    case '1':
        r.x = r.x * cos(rate) + l.x * sin(rate);
        r.y = r.y * cos(rate) + l.y * sin(rate);
        r.z = r.z * cos(rate) + l.z * sin(rate);

        l.x = l.x * cos(rate) - r.x * sin(rate);
        l.y = l.y * cos(rate) - r.y * sin(rate);
        l.z = l.z * cos(rate) - r.z * sin(rate);
        break;

    case '2':
        r.x = r.x * cos(-rate) + l.x * sin(-rate);
        r.y = r.y * cos(-rate) + l.y * sin(-rate);
        r.z = r.z * cos(-rate) + l.z * sin(-rate);

        l.x = l.x * cos(-rate) - r.x * sin(-rate);
        l.y = l.y * cos(-rate) - r.y * sin(-rate);
        l.z = l.z * cos(-rate) - r.z * sin(-rate);
        break;

    case '3':
        l.x = l.x * cos(rate) + u.x * sin(rate);
        l.y = l.y * cos(rate) + u.y * sin(rate);
        l.z = l.z * cos(rate) + u.z * sin(rate);

        u.x = u.x * cos(rate) - l.x * sin(rate);
        u.y = u.y * cos(rate) - l.y * sin(rate);
        u.z = u.z * cos(rate) - l.z * sin(rate);
        break;

    case '4':
        l.x = l.x * cos(-rate) + u.x * sin(-rate);
        l.y = l.y * cos(-rate) + u.y * sin(-rate);
        l.z = l.z * cos(-rate) + u.z * sin(-rate);

        u.x = u.x * cos(-rate) - l.x * sin(-rate);
        u.y = u.y * cos(-rate) - l.y * sin(-rate);
        u.z = u.z * cos(-rate) - l.z * sin(-rate);
        break;

    case '5':
        u.x = u.x * cos(rate) + r.x * sin(rate);
        u.y = u.y * cos(rate) + r.y * sin(rate);
        u.z = u.z * cos(rate) + r.z * sin(rate);

        r.x = r.x * cos(rate) - u.x * sin(rate);
        r.y = r.y * cos(rate) - u.y * sin(rate);
        r.z = r.z * cos(rate) - u.z * sin(rate);
        break;

    case '6':
        u.x = u.x * cos(-rate) + r.x * sin(-rate);
        u.y = u.y * cos(-rate) + r.y * sin(-rate);
        u.z = u.z * cos(-rate) + r.z * sin(-rate);

        r.x = r.x * cos(-rate) - u.x * sin(-rate);
        r.y = r.y * cos(-rate) - u.y * sin(-rate);
        r.z = r.z * cos(-rate) - u.z * sin(-rate);
        break;

    case ',':
        if (stepno < 16)
        {
            scaleFactor -= 1.0;
            sphereScaleFactor += 1.0;
            cylinderScalefactor += 1.0;
            cylinderRadius += 1 / (sqrt(3) * 16.0);
            cylinderHeight -= sqrt(2) / 16.0;
            stepno += 1.0;
        }
        break;
    case '.':
        if (stepno > 0)
        {
            scaleFactor += 1.0;
            sphereScaleFactor -= 1.0;
            cylinderScalefactor -= 1.0;
            cylinderRadius -= 1 / (sqrt(3) * 16.0);
            cylinderHeight += sqrt(2) / 16.0;
            stepno -= 1.0;
        }
        break;

    case 'o':
        capture();
    case 'a':
        // std::cout << "hello" << std ::endl;
        angle -= 5.0;
        // glRotatef(angle,0,1,0);
        break;
    case 'd':
        angle += 5.0;
        // glRotatef(angle,0,1,0);
        break;
    case 'w':
        pos.x += v * u.x;
        pos.y += v * u.y;
        pos.z += v * u.z;

        l.x = c.x - pos.x;
        l.y = c.y - pos.y;
        l.z = c.z - pos.z;
        break;
    case 's':
        pos.x -= v * u.x;
        pos.y -= v * u.y;
        pos.z -= v * u.z;

        l.x = c.x - pos.x;
        l.y = c.y - pos.y;
        l.z = c.z - pos.z;
        break;
    // Control exit
    case 27:     // ESC key
        exit(0); // Exit window
        break;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP: // down arrow key
        pos.x += l.x;
        pos.y += l.y;
        pos.z += l.z;
        break;
    case GLUT_KEY_DOWN: // up arrow key
        pos.x -= l.x;
        pos.y -= l.y;
        pos.z -= l.z;
        break;

    case GLUT_KEY_RIGHT:
        pos.x += r.x;
        pos.y += r.y;
        pos.z += r.z;
        break;
    case GLUT_KEY_LEFT:
        pos.x -= r.x;
        pos.y -= r.y;
        pos.z -= r.z;
        break;

    case GLUT_KEY_PAGE_UP:
        pos.x += u.x;
        pos.y += u.y;
        pos.z += u.z;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos.x -= u.x;
        pos.y -= u.y;
        pos.z -= u.z;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:
        break;
    case GLUT_KEY_END:
        break;

    default:
        break;
    }
    glutPostRedisplay();
}

void loaddata()
{

    ifstream file("description.txt");

    file >> near_plane >> far_plane;
    file >> fov_y;
    file >> aspect_ratio;
    file >> recursion_level;
    file >> pixels;
    fov_x = fov_y * aspect_ratio;
    screenHeight = screenWidth = pixels;
    double ambient, diffuse, reflection, specular, shininess;
    // std::cout << near_distance << " " << far_distance << " " << fovY << " " << aspectRatio << " " << recursion_level << " " << screen_height << std::endl;

    double cell_width, a_coeff, d_coeff, r_coeff;
    int num_objects;
    file >> checkerBoardWidth;
    file >> ambient >> diffuse >> reflection;
    checkerBoard_width = checkerBoardWidth;
    checkerBoard_color.r = ambient;
    checkerBoard_color.g = diffuse;
    checkerBoard_color.b = reflection;

    file >> num_objects;

    double px_1, py_1, pz_1;
    double side, radius, width, height;
    double color_r, color_g, color_b;

    for (int i = 0; i < num_objects; i++)
    {
        std::cout << "object " << i << std::endl;
        std ::string s;
        file >> s;
        if (s == "cube")
        {

            file >> px_1 >> py_1 >> pz_1;
            file >> side;
            file >> color_r >> color_g >> color_b;

            file >> ambient >> diffuse >> specular >> reflection;
            file >> shininess;
            Object *obj = new Cube(Vector3D(px_1, py_1, pz_1), side, Color(color_r, color_g, color_b), ambient, diffuse, specular, reflection, shininess);
            objects.push_back(obj);
        }
        else if (s == "sphere")
        {
            file >> px_1 >> py_1 >> pz_1;
            file >> radius;
            file >> color_r >> color_g >> color_b;
            file >> ambient >> diffuse >> specular >> reflection;
            file >> shininess;
            Object *obj = new Sphere(Vector3D(px_1, py_1, pz_1), radius, Color(color_r, color_g, color_b), ambient, diffuse, specular, reflection, shininess);
            objects.push_back(obj);
        }
        else
        {
            file >> px_1 >> py_1 >> pz_1;
            file >> width >> height;
            file >> color_r >> color_g >> color_b;
            file >> ambient >> diffuse >> specular >> reflection;
            file >> shininess;
            Object *obj = new Pyramid(Vector3D(px_1, py_1, pz_1), width, width, height, Color(color_r, color_g, color_b), ambient, diffuse, specular, reflection, shininess);
            objects.push_back(obj);
        }
    }

    int number_of_normal_lights;
    double x, y, z, fall_off_parameter;
    file >> number_of_normal_lights;
    for (int i = 0; i < number_of_normal_lights; i++)
    {

        file >> x >> y >> z;
        file >> fall_off_parameter;

        PointLight *light = new PointLight(Vector3D(x, y, z), Color(255, 255, 255), fall_off_parameter);
        pointLights.push_back(light);
    }

    int number_of_spot_lights;
    double cut_off_angle;
    file >> number_of_spot_lights;
    for (int i = 0; i < number_of_spot_lights; i++)
    {

        file >> x >> y >> z;
        Vector3D p(x, y, z);
        file >> fall_off_parameter;
        file >> x >> y >> z;
        Vector3D direction(x, y, z);
        file >> cut_off_angle;

        SpotLight *light = new SpotLight(p, Color(255, 255, 255), fall_off_parameter, direction, cut_off_angle);
        spotLights.push_back(light);
    }

    file.close();
}

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char **argv)
{

    pos.x = 0;
    pos.y = -160;
    pos.z = 60;

    l.x = 0, l.y = 1, l.z = 0;
    u.x = 0, u.y = 0, u.z = 1;

    r.x = 1;
    r.y = 0;
    r.z = 0;
    loaddata();
    createPointBuffer();
    capture();
    // Sphere *sphere = new Sphere(Vector3D(3, 3, 3), 3, Color(1, 0, 0), 0.4, 0.2, 0.2, 0.2, 10);
    //  Ray ray = Ray(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
    //  cout << sphere->getIntersectingT(ray) << endl;
    //  createBitmapObject();
    //  bitmap_image image = createBitmapObject();
    //  for(int i =0 ;i<screenHeight;i++){
    //      for(int j =0 ;j<screenWidth;j++){
    //          image.set_pixel(j,i,1,1,1);
    //      }
    //  }
    //  image.save_image("output.bmp");

    for (int i = 0; i < objects.size(); i++)
    {
        cout << "reference point " << objects[i]->reference_point.x << " " << objects[i]->reference_point.y << " " << objects[i]->reference_point.z << endl;
        cout << "color " << objects[i]->color.r << " " << objects[i]->color.g << " " << objects[i]->color.b << endl;
        cout << "coEfficients " << objects[i]->coEfficients[0] << " " << objects[i]->coEfficients[1] << " " << objects[i]->coEfficients[2] << " " << objects[i]->coEfficients[3] << " " << objects[i]->shine << endl;
        cout << "length " << objects[i]->length << endl;
        cout << "width " << objects[i]->width << endl;
        cout << "height " << objects[i]->height << endl;
        cout << "shine" << objects[i]->shine << endl;
    }

    glutInit(&argc, argv);                                    // Initialize GLUT
    glutInitWindowSize(768, 768);                             // Set the window's initial width & height
    glutInitWindowPosition(50, 50);                           // Position the window's initial top-left corner
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color
    glutCreateWindow("Ray tracing");                          // Create a window with the given title
    glutDisplayFunc(display);                                 // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);                         // Register callback handler for window re-shape
    // glutKeyboardFunc(keyboard);
    glutKeyboardFunc(keyboardListener);  // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener); // Register callback handler for special-key event
    initGL();                            // Our own OpenGL initialization
    glutMainLoop();                      // Enter the event-processing loop
    return 0;
}
