#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

#include <GL/freeglut.h>
#include <Windows.h>

#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <iostream>
#include <ctime>
#include <bits/stdc++.h>

#include "TriMesh.h"
#include "util.h"

#include "Camera.h"

enum MouseMode {
    MM_CAMERA
};

struct RayHit {
    double t;
    double alpha;
    double beta;
    int mesh_idx;
};

struct RectangleObj {
    Eigen::Vector3d centerPos;
    Eigen::Vector3d widthVec;
    Eigen::Vector3d heightVec;
    Eigen::Vector3d color;
    Eigen::Vector3d n;
    bool is_light;
    double kd;
};

const double __FAR__ = 1.0e33;

const int g_FilmWidth = 640;
const int g_FilmHeight = 480;
float *g_FilmBuffer = nullptr;
GLuint g_FilmTexture = 0;

bool g_DrawFilm = true;

int width = 640;
int height = 480;

MouseMode g_MouseMode = MM_CAMERA;
int mx, my;

double g_FrameSize_WindowSize_Scale_x = 1.0;
double g_FrameSize_WindowSize_Scale_y = 1.0;

Camera g_Camera;

RectangleObj rects[2];

TriMesh *triMeshes;

const int sample_num = 100;

double intensity = 10;

void initFilm() {
    g_FilmBuffer = (float *) malloc(
            sizeof(float) * g_FilmWidth * g_FilmHeight * 3);
    memset(g_FilmBuffer, 0, sizeof(float) * g_FilmWidth * g_FilmHeight * 3);

    glGenTextures(1, &g_FilmTexture);
    glBindTexture(GL_TEXTURE_2D, g_FilmTexture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, g_FilmWidth, g_FilmHeight, 0, GL_RGB,
                 GL_FLOAT, g_FilmBuffer);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

void updateFilm() {
    glBindTexture(GL_TEXTURE_2D, g_FilmTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, g_FilmWidth, g_FilmHeight, GL_RGB,
                    GL_FLOAT, g_FilmBuffer);
}

void drawFilm() {
    Eigen::Vector3d screen_center = g_Camera.getEyePoint() -
                                    g_Camera.getZVector() *
                                    g_Camera.getFocalLength();
    Eigen::Vector3d p1 = screen_center -
                         g_Camera.getXVector() * g_Camera.getScreenWidth() *
                         0.5 -
                         g_Camera.getYVector() * g_Camera.getScreenHeight() *
                         0.5;
    Eigen::Vector3d p2 = screen_center +
                         g_Camera.getXVector() * g_Camera.getScreenWidth() *
                         0.5 -
                         g_Camera.getYVector() * g_Camera.getScreenHeight() *
                         0.5;
    Eigen::Vector3d p3 = screen_center +
                         g_Camera.getXVector() * g_Camera.getScreenWidth() *
                         0.5 +
                         g_Camera.getYVector() * g_Camera.getScreenHeight() *
                         0.5;
    Eigen::Vector3d p4 = screen_center -
                         g_Camera.getXVector() * g_Camera.getScreenWidth() *
                         0.5 +
                         g_Camera.getYVector() * g_Camera.getScreenHeight() *
                         0.5;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, g_FilmTexture);

    glBegin(GL_TRIANGLES);
    glColor3f(1.0, 1.0, 1.0);

    glTexCoord2f(0.0, 1.0);
    glVertex3f(p1.x(), p1.y(), p1.z());
    glTexCoord2f(1.0, 1.0);
    glVertex3f(p2.x(), p2.y(), p2.z());
    glTexCoord2f(1.0, 0.0);
    glVertex3f(p3.x(), p3.y(), p3.z());

    glTexCoord2f(0.0, 1.0);
    glVertex3f(p1.x(), p1.y(), p1.z());
    glTexCoord2f(1.0, 0.0);
    glVertex3f(p3.x(), p3.y(), p3.z());
    glTexCoord2f(0.0, 0.0);
    glVertex3f(p4.x(), p4.y(), p4.z());

    glEnd();

    glDisable(GL_TEXTURE_2D);
}

void clearRayTracedResult() {
    memset(g_FilmBuffer, 0, sizeof(float) * g_FilmWidth * g_FilmHeight * 3);
}

void
makeRectangle(const Eigen::Vector3d &centerPos, const Eigen::Vector3d &widthVec,
              const Eigen::Vector3d &heightVec, const Eigen::Vector3d &color,
              RectangleObj &out_Rect) {
    out_Rect.centerPos = centerPos;
    out_Rect.widthVec = widthVec;
    out_Rect.heightVec = heightVec;
    out_Rect.color = rgbNormalize(color);
    const Eigen::Vector3d v1 = centerPos - widthVec - heightVec;
    const Eigen::Vector3d v2 = centerPos + widthVec + heightVec;
    const Eigen::Vector3d v3 = centerPos - widthVec + heightVec;
    out_Rect.n = ((v1 - v3).cross(v2 - v3)).normalized();
}

void raySquareIntersection(const RectangleObj &rect, const int in_Triangle_idx,
                           const Ray &in_Ray, RayHit &out_Result) {
    out_Result.t = __FAR__;

    const Eigen::Vector3d v1 = rect.centerPos - rect.widthVec - rect.heightVec;
    const Eigen::Vector3d v2 = rect.centerPos + rect.widthVec + rect.heightVec;
    const Eigen::Vector3d v3 = rect.centerPos - rect.widthVec + rect.heightVec;

    Eigen::Vector3d triangle_normal = (v1 - v3).cross(v2 - v3);
    triangle_normal.normalize();

    const double denominator = triangle_normal.dot(in_Ray.d);
    if (denominator >= 0.0)
        return;

    const double t = triangle_normal.dot(v3 - in_Ray.o) / denominator;
    if (t <= 0.0)
        return;

    const Eigen::Vector3d x = in_Ray.o + t * in_Ray.d;

    Eigen::Matrix<double, 3, 2> A;
    A.col(0) = v1 - v3;
    A.col(1) = v2 - v3;

    Eigen::Matrix2d ATA = A.transpose() * A;
    const Eigen::Vector2d b = A.transpose() * (x - v3);

    const Eigen::Vector2d alpha_beta = ATA.inverse() * b;

    if (alpha_beta.x() < 0.0 || 1.0 < alpha_beta.x() || alpha_beta.y() < 0.0 ||
        1.0 < alpha_beta.y())
        return;

    out_Result.t = t;
    out_Result.alpha = alpha_beta.x();
    out_Result.beta = alpha_beta.y();
}

void rayTriangleIntersection(const TriMesh &triMesh, const Ray &in_Ray,
                             RayHit &out_Result) {
    out_Result.t = __FAR__;

    const double denominator = triMesh.n.dot(in_Ray.d);
    if (denominator >= 0.0)
        return;

    const double t = triMesh.n.dot(triMesh.v3 - in_Ray.o) / denominator;
    if (t <= 0.0)
        return;

    const Eigen::Vector3d x = in_Ray.o + t * in_Ray.d;

    Eigen::Matrix<double, 3, 2> A;
    A.col(0) = triMesh.v1 - triMesh.v3;
    A.col(1) = triMesh.v2 - triMesh.v3;

    Eigen::Matrix2d ATA = A.transpose() * A;
    const Eigen::Vector2d b = A.transpose() * (x - triMesh.v3);

    const Eigen::Vector2d alpha_beta = ATA.inverse() * b;

    /**
     * 0 < alpha < 1 && 0 < 0 < beta < 1
     */
    if (alpha_beta.x() < 0.0 || 1.0 < alpha_beta.x() || alpha_beta.y() < 0.0
        || 1.0 < alpha_beta.y())
        return;

    /**
     * 0 < 1 - alpha - beta < 1
     */
    if (1 - alpha_beta.x() - alpha_beta.y() < 0.0
        || 1.0 < 1 - alpha_beta.x() - alpha_beta.y())
        return;

    out_Result.t = t;
    out_Result.alpha = alpha_beta.x();
    out_Result.beta = alpha_beta.y();
}


void rayTracing(const Ray &in_Ray, RayHit &io_Hit) {
    double t_min = __FAR__;
    double alpha_I = 0.0, beta_I = 0.0;
    int mesh_idx = -1;

    for (int i = 0; i < (int) (sizeof(triMeshes) / sizeof(TriMesh)); i++) {
        RayHit temp_hit{};
        rayTriangleIntersection(triMeshes[i], in_Ray, temp_hit);
        if (temp_hit.t < t_min) {
            t_min = temp_hit.t;
            alpha_I = temp_hit.alpha;
            beta_I = temp_hit.beta;
            mesh_idx = i;
        }
    }

    io_Hit.t = t_min;
    io_Hit.alpha = alpha_I;
    io_Hit.beta = beta_I;
    io_Hit.mesh_idx = mesh_idx;
}

void Learn() {
    for (int y = 0; y < g_FilmHeight; y++) {
        for (int x = 0; x < g_FilmWidth; x++) {
            const int pixel_idx = y * g_FilmWidth + x;

            const double p_x = (x + 0.5) / g_FilmWidth;
            const double p_y = (y + 0.5) / g_FilmHeight;

            Ray ray;
            g_Camera.screenView(p_x, p_y, ray);

            RayHit ray_hit;
            rayTracing(ray, ray_hit);
            if (ray_hit.mesh_idx >= 0) {
                if (triMeshes[ray_hit.mesh_idx].is_light == true) {
                    // 当たった四角形が光源ならば白を返す
                    g_FilmBuffer[pixel_idx * 3]
                            = triMeshes[ray_hit.mesh_idx].color.x();
                    g_FilmBuffer[pixel_idx * 3 + 1]
                            = triMeshes[ray_hit.mesh_idx].color.y();
                    g_FilmBuffer[pixel_idx * 3 + 2]
                            = triMeshes[ray_hit.mesh_idx].color.z();
                } else {
                    // 四方八方にレイを飛ばして、色を決める
                    Eigen::Vector3d sum{0.0, 0.0, 0.0};
                    for (int n = 0; n < sample_num; n++) {
                        const Eigen::Vector3d x = ray.o + ray_hit.t * ray.d;
                        const double theta = asin(sqrt(drand48()));
                        const double phi = 2 * M_PI * drand48();

                        const Eigen::Vector3d omega{
                                sin(theta) * cos(phi),
                                cos(theta),
                                sin(theta) * sin(phi)
                        };

                        Ray _ray;
                        _ray.o = x;
                        _ray.d = omega;

                        RayHit _ray_hit;

                        rayTracing(_ray, _ray_hit);

                        if (_ray_hit.mesh_idx == -1 ||
                            !triMeshes[_ray_hit.mesh_idx].is_light) {
                            sum = sum + Eigen::Vector3d{0.0, 0.0, 0.0};
                        } else {
                            const Eigen::Vector3d l = x + _ray_hit.t * omega;
                            sum = sum +
                                  intensity *
                                  triMeshes[_ray_hit.mesh_idx].color;
                        }
                    }
                    const Eigen::Vector3d i_n = sum / sample_num;
                    // kd は実際はr,g,b各成分に対しての反射率を持つ。
                    const Eigen::Vector3d pixel_color =
                            triMeshes[ray_hit.mesh_idx].kd *
                            triMeshes[ray_hit.mesh_idx].color.cwiseProduct(i_n);
                    g_FilmBuffer[pixel_idx * 3] = pixel_color.x();
                    g_FilmBuffer[pixel_idx * 3 + 1] = pixel_color.y();
                    g_FilmBuffer[pixel_idx * 3 + 2] = pixel_color.z();
                }
            } else {
                const Eigen::Vector3d pixel_color = rgbNormalize(
                        Eigen::Vector3d{173, 216, 230});

                g_FilmBuffer[pixel_idx * 3] = pixel_color.x();
                g_FilmBuffer[pixel_idx * 3 + 1] = pixel_color.y();
                g_FilmBuffer[pixel_idx * 3 + 2] = pixel_color.z();
            }
        }
    }
    updateFilm();
    glutPostRedisplay();
}

void drawRectangleObj(const RectangleObj &rect) {
    glBegin(GL_TRIANGLES);
    glColor3f(rect.color(0), rect.color(1), rect.color(2));
    Eigen::Vector3d vertice1 = rect.centerPos + rect.widthVec + rect.heightVec;
    Eigen::Vector3d vertice2 = rect.centerPos - rect.widthVec + rect.heightVec;
    Eigen::Vector3d vertice3 = rect.centerPos - rect.widthVec - rect.heightVec;
    Eigen::Vector3d vertice4 = rect.centerPos + rect.widthVec - rect.heightVec;

    glVertex3f(vertice1(0), vertice1(1), vertice1(2));
    glVertex3f(vertice2(0), vertice2(1), vertice2(2));
    glVertex3f(vertice3(0), vertice3(1), vertice3(2));

    glVertex3f(vertice1(0), vertice1(1), vertice1(2));
    glVertex3f(vertice3(0), vertice3(1), vertice3(2));
    glVertex3f(vertice4(0), vertice4(1), vertice4(2));

    glEnd();
}

void drawTriMesh(const TriMesh &triMesh) {
    glBegin(GL_TRIANGLES);
    glColor3f(triMesh.color(0), triMesh.color(1), triMesh.color(2));

    glVertex3f(triMesh.v1.x(), triMesh.v1.y(), triMesh.v1.z());
    glVertex3f(triMesh.v2.x(), triMesh.v2.y(), triMesh.v2.z());
    glVertex3f(triMesh.v3.x(), triMesh.v3.y(), triMesh.v3.z());

    glEnd();
}

void mouseDrag(int x, int y) {
    int _dx = x - mx, _dy = y - my;
    mx = x;
    my = y;

    double dx = double(_dx) / double(width);
    double dy = -double(_dy) / double(height);

    if (g_MouseMode == MM_CAMERA) {
        double scale = 2.0;

        g_Camera.rotateCameraInLocalFrameFixLookAt(dx * scale);
        glutPostRedisplay();
    }
}

void mouseDown(int x, int y) {
    mx = x;
    my = y;
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
        mouseDown(x, y);
}

void key(unsigned char key, int x, int y) {
    switch (key) {
        case 'C':
        case 'c':
            g_MouseMode = MM_CAMERA;
            break;
        case 'f':
        case 'F':
            g_DrawFilm = !g_DrawFilm;
            glutPostRedisplay();
            break;
    }
}

void projection_and_modelview(const Camera &in_Camera) {
    const double fovy_deg = (2.0 * 180.0 / M_PI) *
                            atan(0.024 * 0.5 / in_Camera.getFocalLength());

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovy_deg, double(width) / double(height),
                   0.01 * in_Camera.getFocalLength(), 1000.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    const Eigen::Vector3d lookAtPoint = in_Camera.getLookAtPoint();
    gluLookAt(in_Camera.getEyePoint().x(), in_Camera.getEyePoint().y(),
              in_Camera.getEyePoint().z(), lookAtPoint.x(), lookAtPoint.y(),
              lookAtPoint.z(), in_Camera.getYVector().x(),
              in_Camera.getYVector().y(), in_Camera.getYVector().z());
}

void drawFloor() {
    glBegin(GL_TRIANGLES);
    for (int j = -20; j < 20; j++) {
        for (int i = -20; i < 20; i++) {
            int checker_bw = (i + j) % 2;
            if (checker_bw == 0) {
                glColor3f(0.3, 0.3, 0.3);

                glVertex3f(i * 0.5, 0.0, j * 0.5);
                glVertex3f(i * 0.5, 0.0, (j + 1) * 0.5);
                glVertex3f((i + 1) * 0.5, 0.0, j * 0.5);

                glVertex3f(i * 0.5, 0.0, (j + 1) * 0.5);
                glVertex3f((i + 1) * 0.5, 0.0, (j + 1) * 0.5);
                glVertex3f((i + 1) * 0.5, 0.0, j * 0.5);
            }
        }
    }
    glEnd();
}

void display() {
    glViewport(0, 0, width * g_FrameSize_WindowSize_Scale_x,
               height * g_FrameSize_WindowSize_Scale_y);

    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    projection_and_modelview(g_Camera);

    glEnable(GL_DEPTH_TEST);

    Learn();

    if (g_DrawFilm)
        drawFilm();

    drawRectangleObj(rects[0]);
    drawRectangleObj(rects[1]);

    glDisable(GL_DEPTH_TEST);

    glutSwapBuffers();
}

void resize(int w, int h) {
    width = w;
    height = h;
}

int main(int argc, char *argv[]) {
    g_Camera.setEyePoint(Eigen::Vector3d{0.0, 1.0, 4.0});
    g_Camera.lookAt(Eigen::Vector3d{0.0, 0.5, 0.0},
                    Eigen::Vector3d{0.0, 1.0, 0.0});

    glutInit(&argc, argv);
    glutInitWindowSize(width, height);
    glutInitDisplayMode(
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);

    glutCreateWindow("Hello world!!");

    // With retina display, frame buffer size is twice the window size.
    // Viewport size should be set on the basis of the frame buffer size, rather than the window size.
    // g_FrameSize_WindowSize_Scale_x and g_FrameSize_WindowSize_Scale_y account for this factor.
    GLint dims[4] = {0};
    glGetIntegerv(GL_VIEWPORT, dims);
    g_FrameSize_WindowSize_Scale_x = double(dims[2]) / double(width);
    g_FrameSize_WindowSize_Scale_y = double(dims[3]) / double(height);

    /**
     * メッシュの配置
     */
    const int num_of_elements = 2;
    /**
     * vectorで要素数を確認できるようにすれば解決？？
     * もしくは要素数をグローバルで管理する？？
     */
    triMeshes = new TriMesh[num_of_elements];
    memset(triMeshes, 0, sizeof(TriMesh) * num_of_elements);
    triMeshes[0].setTriMesh(Eigen::Vector3d{0, 1, 0},
                            Eigen::Vector3d{-1, 0, 0},
                            Eigen::Vector3d{1, 0, 0},
                            Eigen::Vector3d{255, 255, 255}, true, 1);

    triMeshes[1].setTriMesh(Eigen::Vector3d{0, 0, 1},
                            Eigen::Vector3d{1, 0, 0},
                            Eigen::Vector3d{-1, 0, 0},
                            Eigen::Vector3d{100, 100, 255}, false, 1);

    std::cout << "各TriMeshの情報を表示します" << std::endl;
    std::cout << sizeof(triMeshes) << std::endl;
    std::cout << sizeof(TriMesh) << std::endl;


    for (int i = 0; i < (int) (sizeof(triMeshes) / sizeof(TriMesh)); i++) {
        triMeshes[i].printInfo();
    }
    /**
     * 乱数シード値の設定
     */
    randInit();

    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMouseFunc(mouse);
    glutMotionFunc(mouseDrag);
    glutKeyboardFunc(key);

    initFilm();
    clearRayTracedResult();


    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    glutMainLoop();
    return 0;
}
