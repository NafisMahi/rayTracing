#include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

GLfloat triangle=0.0,stepsize=16.0,stepno=0,cylinderRadius=0.0,cylinderHeight=sqrt(2.0),scaleFactor=16.0,cylinderScalefactor=0.0,sphereScaleFactor=0.0,angle=0.0;
/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
}

// Global variables
GLfloat eyex = 4, eyey = 4, eyez = 4;
GLfloat centerx = 0, centery = 0, centerz = 0;
GLfloat upx = 0, upy = 1, upz = 0;
bool isAxes = true, isCube = false, isPyramid = false;

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
        glColor3f(1,0,0);   // Red
        // X axis
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);
  
        glColor3f(0,1,0);   // Green
        // Y axis
        glVertex3f(0,0,0);
        glVertex3f(0,1,0);

        glColor3f(0,0,1);   // Blue
        // Z axis
        glVertex3f(0,0,0);
        glVertex3f(0,0,1);
    glEnd();
}

void glDraw3Dtriangle()
{
    glBegin(GL_TRIANGLES);

       glVertex3f(1,0,0);
       glVertex3f(0,1,0);
       glVertex3f(0,0,1);
       
    glEnd();

}

void drawScaledTriangle()
{
    glPushMatrix();
        glTranslated(stepno/48.0,stepno/48.0,stepno/48.0);
        glScaled(scaleFactor/16.0,scaleFactor/16.0,scaleFactor/16.0);
        glDraw3Dtriangle();
    glPopMatrix();        
}

void glDrawTetrahedron(float r,float g,float b,float x)
{
    //glDraw3Dtriangle();
    glColor3f(r,g,b);
    //glDraw3Dtriangle();
    drawScaledTriangle();

    glPushMatrix();
        glRotatef(90,0,1,0);
        glColor3f(1.0-r,1.0-g,0.0f);
       // glDraw3Dtriangle();
        drawScaledTriangle();
        glRotatef(90,0,1,0);
        glColor3f(r,g,b);
        //glDraw3Dtriangle();
        drawScaledTriangle();
        glRotatef(90,0,1,0);
        glColor3f(1.0-r,1.0-g,0.0f);
        //glDraw3Dtriangle();
        drawScaledTriangle();
    
    glPopMatrix();
}

void drawOctahedron(float x)
{
    glDrawTetrahedron(1.0f,0.0f,0.0f,x);
    
    glPushMatrix();
       glRotatef(180,1,0,0);
       glDrawTetrahedron(1.0f,0.0f,0.0f,x);
    glPopMatrix();
}

double diffAmount()
{
   return (1.0-1.0/3.0)/16.0;
}

std::vector<float> verticesCubeSphere(int subdivision)
{
  const float DEG2RAD = acos(-1) / 180.0f;

    std::vector<float> vertices;
    float n1[3];        // normal of longitudinal plane rotating along Y-axis
    float n2[3];        // normal of latitudinal plane rotating along Z-axis
    float v[3];         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = 1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            // add a vertex into array
            vertices.push_back(v[0]);
            vertices.push_back(v[1]);
            vertices.push_back(v[2]);
        }
    }

    return vertices;
}

void drawOneSixthSphere()
{
    glBegin(GL_POLYGON);
       std::vector<float> vertices = verticesCubeSphere(8);
       for(int i = 0;i < vertices.size(); i+=3)
       {
         glVertex3f(vertices[i],vertices[i+1],vertices[i+2]);
       }
    glEnd();
}

void drawScaledFaceOfSphere()
{
    glPushMatrix();
        glTranslated((16.0 - stepno)/16.0,0.0,0.0);
        glScaled(sphereScaleFactor/16.0,sphereScaleFactor/16.0,sphereScaleFactor/16.0);
        glScaled(1/sqrt(3.0),1/sqrt(3.0),1/sqrt(3.0));
        drawOneSixthSphere();
    glPopMatrix();   
}

void glDrawSphere()
{
    glColor3f(0.0f,1.0f,0.0f);
    //drawOneSixthSphere();
    // glPushMatrix();
    //     glScaled(1/3.0,1/3.0,1/3.0);
    //     drawScaledFaceOfSphere();
    // glPopMatrix();
    drawScaledFaceOfSphere();

    glPushMatrix();
        //glScaled(1/3.0,1/3.0,1/3.0);
        glRotatef(90,0,1,0);
        glColor3f(0.0f,0.0f,1.0f);
        // drawOneSixthSphere();
        drawScaledFaceOfSphere();
        
        glRotatef(90,0,1,0);
        glColor3f(0.0f,1.0f,0.0f);
        // drawOneSixthSphere();
        drawScaledFaceOfSphere();

        glRotatef(90,0,1,0);
        glColor3f(0.0f,0.0f,1.0f);
        // drawOneSixthSphere();
        drawScaledFaceOfSphere();
    glPopMatrix();

    glPushMatrix();
        glRotatef(90,0,0,1);
        glColor3f(1.0f,0.0f,0.0f);
        // drawOneSixthSphere();
        drawScaledFaceOfSphere();
        
        glRotatef(180,0,0,1);
        glColor3f(1.0f,0.0f,0.0f);
        // drawOneSixthSphere();
        drawScaledFaceOfSphere();
    glPopMatrix();
}

void draw_cylinder(GLfloat radius,
                   GLfloat height,
                   GLubyte R,
                   GLubyte G,
                   GLubyte B)
{
    GLfloat x              = 0.0;
    GLfloat y              = 0.0;
    GLfloat angle          = 0.0;
    GLfloat angle_stepsize = 0.1;

    /** Draw the tube */
    glColor3ub(R-40,G-40,B-40);
    glBegin(GL_QUAD_STRIP);
    angle = 0;
        while( angle <  71 ) {
            x = radius * cos((angle*M_PI)/180.0);
            y = radius * sin((angle*M_PI)/180.0);
            glVertex3f(x, height/2.0 , y);
            glVertex3f(x, -height/2.0 , y);
            angle = angle + angle_stepsize;
        }
        glVertex3f(radius, height/2.0, 0.0);
        glVertex3f(radius, -height/2.0, 0.0);
    glEnd();

}

void drawScaledCylinder()
{
    glPushMatrix();
        //glTranslated((16 - stepno)/32.0,(16 - stepno)/32.0,0.0);
        glTranslated((16 - stepno)/(2*16),(16 - stepno)/(2*16),0.0);
        //glScaled(cylinderScalefactor/16.0,cylinderScalefactor/16.0,cylinderScalefactor/16.0);
        glRotatef(45.0,0,0,1);
        glRotatef(35.3,0,1,0);
        draw_cylinder(cylinderRadius,cylinderHeight,255,140,200);
    glPopMatrix();
}

void drawCylinders()
{
   drawScaledCylinder();

    // glPushMatrix();
    //     glRotatef(90,0,1,0);
    //     drawScaledCylinder();
        
    //     glRotatef(90,0,1,0);
    //     drawScaledCylinder();

    //     glRotatef(90,0,1,0);
    //     drawScaledCylinder();

    //     glRotatef(90,1,0,0);
    //     drawScaledCylinder();

    //     drawScaledCylinder();
        
    //     glRotatef(90,0,0,1);
    //     drawScaledCylinder();

    //     glRotatef(90,0,0,1);
    //     drawScaledCylinder();

    //     glRotatef(90,0,0,1);
    //     drawScaledCylinder();

    //     glRotatef(90,1,0,0);
    //     drawScaledCylinder();

    //      glRotatef(90,0,1,0);
    //      drawScaledCylinder();

    //      glRotatef(90,0,1,0);
    //      drawScaledCylinder();

    //      glRotatef(90,0,1,0);
    //      drawScaledCylinder();
    
    glPopMatrix();
   //drawScaledCylinder();

    //  glPushMatrix();
    //    glRotatef(90,1,0,0);
    //    drawScaledCylinder();
       
    //    //glRotatef(-90,1,0,0);
    //    glRotatef(90,0,1,0);
    // //    drawScaledCylinder();

    // //    glRotatef(9g0,0,1,0);
    // //    drawScaledCylinder();
    //  glPopMatrix();

    //    glPushMatrix();
    //        glRotatef(90,0,0,1);
    //        drawScaledCylinder();
    // //   glPopMatrix();

}

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display() {
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);

    // control viewing (or camera)
    gluLookAt(eyex,eyey,eyez,
              centerx,centery,centerz,
              upx,upy,upz);
    drawAxes();
    //glDraw3Dtriangle();
    //drawOctahedron();
   // drawOctahedron();
    //drawOneSixthSphere();
    //   glDrawSphere();
    //   drawOctahedron(triangle);
    // glPushMatrix();
    //      glTranslated(1/2.0,1/2.0,0);
    //      glRotatef(45.0,0,0,1);
    //     //  glRotatef(35.3,0,1,0);
    //      //glRotatef(45,0,0,1);   
         
    //      draw_cylinder(1/sqrt(3.0), sqrt(2.0), 255, 160, 100);
    // glPopMatrix();
    //drawScaledCylinder();
    glPushMatrix();
       drawCylinders();
    glPopMatrix();
    //glDrawTetrahedron(1.0f,0.0f,0.0f,triangle);
    glutSwapBuffers();  // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

void keyboard(unsigned char key, int x, int y) {
    // key is the char pressed, e.g., 'a' or 27 for ESC
    // (x, y) is the mouse location in Windows' coordinates
    std :: cout << "helloooo" << std :: endl;
    switch (key) {
    case 'u':
        triangle -= diffAmount();
        std::cout << triangle << std :: endl;
        std::cout << "heloooo" << std :: endl;
        break;
    case 'i':
        triangle += diffAmount();
    
    default:
        break;
    }
    glutPostRedisplay();
}

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y) {
    std :: cout << "helloooo" << std :: endl;
    float v = 0.1;
    switch (key) {
    // Control eye (location of the eye)
    // control eyex
    case '1':
        eyex += v;
        break;
    case '2':
        eyex -= v;
        break;
    // control eyey
    case '3':
        eyey += v;
        break;
    case '4':
        eyey -= v;
        break;
    // control eyez
    case '5':
        eyez += v;
        break;
    case '6':
        eyez -= v;
        break;

    // Control center (location where the eye is looking at)
    // control centerx
    case 'q':
        centerx += v;
        break;
    case 'w':
        centerx -= v;
        break;
    // control centery
    case 'e':
        centery += v;
        break;
    case 'r':
        centery -= v;
        break;
    // control centerz
    case 't':
        centerz += v;
        break;
    case 'y':
        centerz -= v;
        break;

    // Control what is shown
    case 'a':
        isAxes = !isAxes;   // show/hide Axes if 'a' is pressed
        break;
    case 'c':
        isCube = !isCube;   // show/hide Cube if 'c' is pressed
        break;
    case 'p':
        isPyramid = !isPyramid; // show/hide Pyramid if 'p' is pressed
        break;
    case 'u':
        if(stepno < 16)
        {
          scaleFactor -= 1.0;
          sphereScaleFactor += 1.0;
          cylinderScalefactor +=1.0;
          cylinderRadius += 1/(sqrt(3)*16.0);
          cylinderHeight -= sqrt(2)/16.0;
          std :: cout << cylinderHeight << std :: endl;
          stepno += 1.0;
        }
        break;
    case 'i':
        if(stepno > 0)
        {
          scaleFactor += 1.0;
          sphereScaleFactor -= 1.0;
          cylinderScalefactor -= 1.0;
          cylinderRadius -= 1/(sqrt(3)*16.0);
          cylinderHeight += sqrt(2)/16.0;
          stepno -= 1.0;
        }
         break;

    // Control exit
    case 27:    // ESC key
        exit(0);    // Exit window
        break;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}

/* Callback handler for special-key event */
void specialKeyListener(int key, int x,int y) {
    double v = 0.25;
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;
    switch (key) {
    case GLUT_KEY_LEFT:
        eyex += v * (upy*lz);
        eyez += v * (-lx*upy);
        s = sqrt(eyex*eyex + eyez*eyez) / (4 * sqrt(2));
        eyex /= s;
        eyez /= s;
        break;
    case GLUT_KEY_RIGHT:
        eyex += v * (-upy*lz);
        eyez += v * (lx*upy);
        s = sqrt(eyex*eyex + eyez*eyez) / (4 * sqrt(2));
        eyex /= s;
        eyez /= s;
        break;
    case GLUT_KEY_UP:
        eyey += v;
        break;
    case GLUT_KEY_DOWN:
        eyey -= v;
        break;
    
    default:
        return;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    glutInit(&argc, argv);                      // Initialize GLUT
    glutInitWindowSize(640, 640);               // Set the window's initial width & height
    glutInitWindowPosition(50, 50);             // Position the window's initial top-left corner
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("OpenGL 3D Drawing");      // Create a window with the given title
    glutDisplayFunc(display);                   // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape
    //glutKeyboardFunc(keyboard);
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);        // Register callback handler for special-key event
    initGL();                                   // Our own OpenGL initialization
    glutMainLoop();                             // Enter the event-processing loop
    return 0;
}
