#include<GL/glut.h>

void start_animation(int argc, char **argv)
{
     glutInit(&argc,argv);
     glutInitWindowSize(512,512);
     glutCreateWindow("SpatCorr");
}

void anim_loop()
{
     glutMainLoop();
}

void display()
{
     glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
     glFlush();
}
