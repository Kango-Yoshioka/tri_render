#include <GL/freeglut.h>

#include <math.h>
#include <vector>
#include <iostream>
#include <bits/stdc++.h>

#include "TriMesh.h"
#include "util.h"

#include "Camera.h"
#include "TriMeshObj.h"

enum MouseMode {
    MM_CAMERA
};

struct RayHit {
    double t;
    double alpha;
    double beta;
    int mesh_idx;
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

std::vector<TriMesh> triMeshes;
std::vector<int> lightTriMeshIdxes;

std::vector<double> sample_lights_cdf;

const int sample_num = 20;

double intensity = 100;

float *g_AccumulationBuffer = nullptr;
int *g_CountBuffer = nullptr;

bool save_flag = false;

constexpr int save_border = 1000;

const Eigen::Vector3d background_color{215, 230, 250};

constexpr int method_idx = 2;

double S_A = 0; // 光源面の面積の合計

void initFilm() {
    g_FilmBuffer = (float *) malloc(
            sizeof(float) * g_FilmWidth * g_FilmHeight * 3);
    memset(g_FilmBuffer, 0, sizeof(float) * g_FilmWidth * g_FilmHeight * 3);

    g_AccumulationBuffer = (float *) malloc(sizeof(float) * g_FilmWidth *
                                            g_FilmHeight * 3);
    g_CountBuffer = (int *) malloc(sizeof(int) * g_FilmWidth * g_FilmHeight);
    glGenTextures(1, &g_FilmTexture);
    glBindTexture(GL_TEXTURE_2D, g_FilmTexture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, g_FilmWidth, g_FilmHeight, 0,
                 GL_RGB, GL_FLOAT, g_FilmBuffer);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

void saveFilm() {
    std::cout << "Save start" << std::endl;
    std::ofstream writing_file;
    std::string filename = "method_" + std::to_string(method_idx) + ".ppm";

    const int scale = 10000;

    writing_file.open(filename, std::ios::out);
    const std::string header = "P3\n" + std::to_string(g_FilmWidth) + " " +
                               std::to_string(g_FilmHeight) + "\n" +
                               std::to_string(scale);
    writing_file << header << std::endl;

    for (int i = 0; i < g_FilmHeight * g_FilmWidth; i++) {
        std::string ppm_text = "";

        ppm_text += std::to_string(int(g_FilmBuffer[i * 3] * scale)) + " "
                    + std::to_string(int(g_FilmBuffer[i * 3 + 1] * scale)) +
                    " " + std::to_string(int(g_FilmBuffer[i * 3 + 2] * scale));

        if (i % g_FilmWidth == (g_FilmWidth - 1)) {
            ppm_text += "\n";
        } else {
            ppm_text += " ";
        }

        writing_file << ppm_text;
    }
    writing_file.close();

    std::cout << "Save complete" << std::endl;
}

void updateFilm() {
    for (int i = 0; i < g_FilmWidth * g_FilmHeight; i++) {
        if (g_CountBuffer[i] > 0) {
            g_FilmBuffer[i * 3] = dClamp(g_AccumulationBuffer[i * 3]
                                         / g_CountBuffer[i], 0, 1);
            g_FilmBuffer[i * 3 + 1] = dClamp(g_AccumulationBuffer[i * 3 + 1]
                                             / g_CountBuffer[i], 0, 1);
            g_FilmBuffer[i * 3 + 2] = dClamp(g_AccumulationBuffer[i * 3 + 2]
                                             / g_CountBuffer[i], 0, 1);

            if (g_CountBuffer[i] == save_border && !save_flag) {
                // 保存処理
                saveFilm();
                save_flag = true;
            }
        } else {
            g_FilmBuffer[i * 3] = 0.0;
            g_FilmBuffer[i * 3 + 1] = 0.0;
            g_FilmBuffer[i * 3 + 2] = 0.0;
        }
    }
    glBindTexture(GL_TEXTURE_2D, g_FilmTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, g_FilmWidth, g_FilmHeight, GL_RGB,
                    GL_FLOAT, g_FilmBuffer);

}

void resetFilm() {
    memset(g_AccumulationBuffer, 0, sizeof(float) * g_FilmWidth *
                                    g_FilmHeight * 3);
    memset(g_CountBuffer, 0, sizeof(int) * g_FilmWidth * g_FilmHeight);

    save_flag = false;
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

    for (int i = 0; i < triMeshes.size(); i++) {
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

void directionalSampling(const int &pixel_idx, Eigen::Vector3d &pixel_color,
                         const Ray &ray, const RayHit &ray_hit) {
    // 四方八方にレイを飛ばして、色を決める
    Eigen::Vector3d sum{0.0, 0.0, 0.0};
    const Eigen::Vector3d x = ray.o + ray_hit.t * ray.d;
    const TriMesh diffuseTriMesh = triMeshes[ray_hit.mesh_idx];
    for (int n = 0; n < sample_num; n++) {
        const double theta = asin(sqrt(drand48()));
        const double phi = 2 * M_PI * drand48();

        /**
         * 拡散面上の半球でサンプリングをする
         * 拡散面で新たに基底ベクトルを設定しランダムな方向omegaを決定する
         * p = (v2 - v1).normalized(), q = n, r = (q × p).normalized()
         * というp, q, rの三本のベクトルを基底ベクトルとしてランダムな方向omega
         * 決めたとすると、デカルト座標系においては,
         * x = ox * px + oy * qx + oz * rx
         * y = ox * py + oy * qy + oz * ry
         * x = ox * pz + oy * qz + oz * rz と表せる
         */

        const Eigen::Vector3d p = (diffuseTriMesh.v2 -
                                   diffuseTriMesh.v1).normalized();
        const Eigen::Vector3d q = diffuseTriMesh.n;
        const Eigen::Vector3d r = (q.cross(p)).normalized();

        const Eigen::Vector3d omega{
                sin(theta) * cos(phi),
                cos(theta),
                sin(theta) * sin(phi)
        };

        const Eigen::Vector3d sampledDirection = {
                omega.x() * p.x() + omega.y() * q.x() + omega.z() * r.x(),
                omega.x() * p.y() + omega.y() * q.y() + omega.z() * r.y(),
                omega.x() * p.z() + omega.y() * q.z() + omega.z() * r.z()
        };
        Ray _ray;
        _ray.o = x;
        _ray.d = sampledDirection;

        RayHit _ray_hit;
        rayTracing(_ray, _ray_hit);

        // 飛ばしたレイが面に当たり、かつ光源である時、
        if (_ray_hit.mesh_idx != -1 &&
            triMeshes[_ray_hit.mesh_idx].is_light) {
            const Eigen::Vector3d L_in = intensity *
                                         triMeshes[_ray_hit.mesh_idx].color;
            sum = sum + L_in * diffuseTriMesh.kd;
        }
    }
    const Eigen::Vector3d I_n = sum;
    pixel_color = diffuseTriMesh.color.cwiseProduct(I_n);
    g_CountBuffer[pixel_idx] += sample_num;
}

void sampleLightTriMesh(int &out_sampled_light_triMesh_idx) {
    const double xi = drand48();
    for (int i = 0; i < sample_lights_cdf.size(); i++) {
        if (xi <= sample_lights_cdf[i]) {
            out_sampled_light_triMesh_idx = lightTriMeshIdxes[i];
            return;
        }
    }
    std::cerr << "Not selected Light" << std::endl;
    exit(EXIT_FAILURE);
}

void isVisible(const Eigen::Vector3d &x, const Eigen::Vector3d &y,
               const int &x_triMesh_idx, const int &y_triMesh_idx,
               int &out_is_visible) {
    int is_visible_x_y;
    int is_visible_y_x;
    /**
     * x->yまで光が届くか？
     */
    Ray ray_x_y;
    ray_x_y.o = x;
    ray_x_y.d = (y - x).normalized();
    RayHit ray_hit_x_y;
    rayTracing(ray_x_y, ray_hit_x_y);
    /**
     * レイが面に当たり、その面がyが存在するメッシュであるか？
     */
    if (ray_hit_x_y.mesh_idx != -1 && ray_hit_x_y.mesh_idx == y_triMesh_idx) {
        is_visible_x_y = 1;
    } else {
        is_visible_x_y = 0;
    }

    /**
     * y->まで光が届くか？
     */
    Ray _ray_y_x;
    _ray_y_x.o = y;
    _ray_y_x.d = (x - y).normalized();
    RayHit _ray_hit_y_x;
    rayTracing(_ray_y_x, _ray_hit_y_x);
    /**
     * レイが面に当たり、その面がXが存在するメッシュであるか？
     */
    if (_ray_hit_y_x.mesh_idx != -1 && _ray_hit_y_x.mesh_idx == x_triMesh_idx) {
        is_visible_y_x = 1;
    } else {
        is_visible_y_x = 0;
    }
    out_is_visible = is_visible_y_x * is_visible_x_y;
}

void areaSampling(const int &pixel_idx, Eigen::Vector3d &pixel_color,
                  const Ray &ray, const RayHit &ray_hit) {
    Eigen::Vector3d sum{0.0, 0.0, 0.0};
    const TriMesh diffuseTriMesh = triMeshes[ray_hit.mesh_idx];
    for (int i = 0; i < sample_num; i++) {
        /**
         * 複数光源の中から1つをランダムに選ぶ。
         */
        int sampled_light_triMesh_idx;
        sampleLightTriMesh(sampled_light_triMesh_idx);
        TriMesh sampledLightTriMesh = triMeshes[sampled_light_triMesh_idx];
        /**
         * 選ばれた光源面上のランダムな一点をサンプルする。
         * y = v1 + beta * (v2 - v1) + gamma * (v3 - v1)
         */
        const double gamma = 1 - sqrt(drand48());
        const double beta = (1 - gamma) * drand48();
        const Eigen::Vector3d y =
                sampledLightTriMesh.v1 +
                beta * (sampledLightTriMesh.v2 - sampledLightTriMesh.v1) +
                gamma * (sampledLightTriMesh.v3 - sampledLightTriMesh.v1);
        const Eigen::Vector3d x = ray.o + ray_hit.t * ray.d;
        const Eigen::Vector3d y_x = x - y;
        const Eigen::Vector3d w = y_x.normalized();

        /**
         * 光源から拡散面へと光が届くかどうかを判定する。
         */
        int is_visible;
        isVisible(x, y, ray_hit.mesh_idx, sampled_light_triMesh_idx,
                  is_visible);
        const Eigen::Vector3d ny = sampledLightTriMesh.n;
        const Eigen::Vector3d nx = diffuseTriMesh.n;

        const double cosx = (-w).dot(nx);
        const double cosy = w.dot(ny);

        const Eigen::Vector3d L_in = intensity * sampledLightTriMesh.color;
        const double fr = diffuseTriMesh.kd / M_PI;
        const double geometry = cosx * cosy / y_x.squaredNorm() * is_visible;

        sum = sum + S_A * L_in * fr * geometry;
    }
    const Eigen::Vector3d I_n = sum;
    pixel_color = diffuseTriMesh.color.cwiseProduct(I_n);
    g_CountBuffer[pixel_idx] += sample_num;
}

void methods(const int &method_idx, const int &pixel_idx,
             Eigen::Vector3d &pixel_color,
             const Ray &ray, const RayHit &ray_hit) {
    switch (method_idx) {
        case 1:
            directionalSampling(pixel_idx, pixel_color, ray, ray_hit);
            break;
        case 2:
            areaSampling(pixel_idx, pixel_color, ray, ray_hit);
            break;
    }
}

void render() {
    for (int y = 0; y < g_FilmHeight; y++) {
        for (int x = 0; x < g_FilmWidth; x++) {
            const int pixel_idx = y * g_FilmWidth + x;

            const double p_x = (x + 0.5) / g_FilmWidth;
            const double p_y = (y + 0.5) / g_FilmHeight;

            Ray ray;
            g_Camera.screenView(p_x, p_y, ray);

            RayHit ray_hit;
            rayTracing(ray, ray_hit);
            Eigen::Vector3d pixel_color;
            if (ray_hit.mesh_idx >= 0) {
                if (triMeshes[ray_hit.mesh_idx].is_light == true) {
                    // 当たった四角形が光源ならば光源の色を返す
                    pixel_color = triMeshes[ray_hit.mesh_idx].color;
                    g_CountBuffer[pixel_idx] += 1;
                } else {
                    methods(method_idx, pixel_idx, pixel_color, ray, ray_hit);
                }
            } else {
                pixel_color = rgbNormalize(background_color);
                g_CountBuffer[pixel_idx] += 1;
            }

            g_AccumulationBuffer[pixel_idx * 3] += pixel_color.x();
            g_AccumulationBuffer[pixel_idx * 3 + 1] += pixel_color.y();
            g_AccumulationBuffer[pixel_idx * 3 + 2] += pixel_color.z();
        }
    }
    updateFilm();
    glutPostRedisplay();
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
        resetFilm();
        updateFilm();
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

    render();

    if (g_DrawFilm)
        drawFilm();

    for (int i = 0; i < triMeshes.size(); i++) {
        drawTriMesh(triMeshes[i]);
    }

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

    /**
    triMeshes.push_back(
            TriMesh(Eigen::Vector3d{-2.0, 2.25, -2.0},
                    Eigen::Vector3d{-3.0, 0.0, -2.0},
                    Eigen::Vector3d{0.0, 0.0, -2.0},
                    Eigen::Vector3d{215.0, 14.0, 74.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(Eigen::Vector3d{0.25, 0.0, -2.0},
                    Eigen::Vector3d{0.5, 1.5, -2.0},
                    Eigen::Vector3d{-1.5, 2.25, -2.0},
                    Eigen::Vector3d{1.0, 195.0, 215.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(Eigen::Vector3d{1.5, 2.5, -2.0},
                    Eigen::Vector3d{0.75, 1.0, -2.0},
                    Eigen::Vector3d{1.75, 0.0, -2.0},
                    Eigen::Vector3d{50.0, 205.0, 50.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(Eigen::Vector3d{0.75, 0.75, -2.0},
                    Eigen::Vector3d{0.5, 0.0, -2.0},
                    Eigen::Vector3d{1.5, 0.0, -2.0},
                    Eigen::Vector3d{255.0, 244.0, 1.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{2.0, 1.5, -2.0},
                    Eigen::Vector3d{2.0, 0.0, -2.0},
                    Eigen::Vector3d{3.0, 0.0, -2.0},
                    Eigen::Vector3d{147.0, 112.0, 219.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{0.5, 0.25, 0.5},
                    Eigen::Vector3d{0.25, 0.0, 0.5},
                    Eigen::Vector3d{0.75, 0.0, 0.25},
                    Eigen::Vector3d{255.0, 255.0, 255.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{-1.0, 1.5, 1.25},
                    Eigen::Vector3d{-1.0, 0.0, 1.25},
                    Eigen::Vector3d{-1.75, 0.0, 0.5},
                    Eigen::Vector3d{210.0, 105.0, 30.0}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{0.0, 0.0, 1.5},
                    Eigen::Vector3d{3.0, 0.0, -1.5},
                    Eigen::Vector3d{-3.0, 0.0, -1.5},
                    Eigen::Vector3d{255.0, 255.0, 255.0}, false, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{-1.0, 0.0, -1.0},
                    Eigen::Vector3d{-0.5, 0.5, -1.0},
                    Eigen::Vector3d{-1.5, 0.5, -1.0},
                    Eigen::Vector3d{255.0, 255.0, 255.0}, false, 1.0
            ));
    **/
    TriMeshObj triMeshObjs[4];
    triMeshObjs[0] = TriMeshObj(Eigen::Vector3d{255, 255, 255},
                                false, 1.0);
    triMeshObjs[0].setRectangleMesh(Eigen::Vector3d{0.0, -1.5, 0.0},
                                    Eigen::Vector3d{2.5, 0.0, 0.0},
                                    Eigen::Vector3d{0.0, 0.0, -3.0});

    triMeshObjs[1] = TriMeshObj(Eigen::Vector3d{255, 255, 255},
                                false, 1.0);
    triMeshObjs[1].setRectangleMesh(Eigen::Vector3d{0.0, 0.75, -3.0},
                                    Eigen::Vector3d{2.5, 0.0, 0.0},
                                    Eigen::Vector3d{0.0, 2.25, 0.0});

    triMeshObjs[2] = TriMeshObj(Eigen::Vector3d{255, 255, 255},
                                false, 1.0);
    triMeshObjs[2].setRectangleMesh(Eigen::Vector3d{-2.5, 0.75, -0.75},
                                    Eigen::Vector3d{0.0, 0.0, -2.25},
                                    Eigen::Vector3d{0.0, 2.25, 0.0});

    triMeshObjs[3] = TriMeshObj(Eigen::Vector3d{255, 255, 255},
                                false, 1.0);
    triMeshObjs[3].setRectangleMesh(Eigen::Vector3d{2.5, 0.75, -0.75},
                                    Eigen::Vector3d{0.0, 0.0, 2.25},
                                    Eigen::Vector3d{0.0, 2.25, 0.0});

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < triMeshObjs[i].triMeshes.size(); j++) {
            triMeshes.push_back(triMeshObjs[i].triMeshes[j]);
        }
    }

    /**
     * light1
     */
    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{0.0, 0.3, -1.45},
                    Eigen::Vector3d{0.15, -0.3, -1.6},
                    Eigen::Vector3d{-0.15, -0.3, -1.6},
                    Eigen::Vector3d{255, 105, 180}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{0.0, 0.3, -1.45},
                    Eigen::Vector3d{0.0, -0.3, -1.3},
                    Eigen::Vector3d{0.15, -0.3, -1.6},
                    Eigen::Vector3d{255, 105, 180}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{0.0, 0.3, -1.45},
                    Eigen::Vector3d{-0.15, -0.3, -1.6},
                    Eigen::Vector3d{0.0, -0.3, -1.3},
                    Eigen::Vector3d{255, 105, 180}, true, 1.0
            ));

    /**
     * light2
     */
    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{1.4, 0.3, 0.0},
                    Eigen::Vector3d{1.55, -0.3, -0.15},
                    Eigen::Vector3d{1.25, -0.3, -0.15},
                    Eigen::Vector3d{105, 255, 255}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{1.4, 0.3, 0.0},
                    Eigen::Vector3d{1.4, -0.3, 0.15},
                    Eigen::Vector3d{1.55, -0.3, -0.15},
                    Eigen::Vector3d{105, 255, 255}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{1.4, 0.3, 0.0},
                    Eigen::Vector3d{1.25, -0.3, -0.15},
                    Eigen::Vector3d{1.4, -0.3, 0.15},
                    Eigen::Vector3d{105, 255, 255}, true, 1.0
            ));

    /**
     * light3 105 G:255 B:105
     */
    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{-1.4, 0.3, 0.0},
                    Eigen::Vector3d{-1.25, -0.3, -0.15},
                    Eigen::Vector3d{-1.55, -0.3, -0.15},
                    Eigen::Vector3d{105, 255, 105}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{-1.4, 0.3, 0.0},
                    Eigen::Vector3d{-1.4, -0.3, 0.15},
                    Eigen::Vector3d{-1.25, -0.3, -0.15},
                    Eigen::Vector3d{105, 255, 105}, true, 1.0
            ));

    triMeshes.push_back(
            TriMesh(
                    Eigen::Vector3d{-1.4, 0.3, 0.0},
                    Eigen::Vector3d{-1.55, -0.3, -0.15},
                    Eigen::Vector3d{-1.4, -0.3, 0.15},
                    Eigen::Vector3d{105, 255, 105}, true, 1.0
            ));

    /**
     * 配置メッシュの情報の表示
     */
    std::cout << "Print triMeshes information" << std::endl;
    std::cout << triMeshes.size() << " elements in triMeshes" << std::endl;

    for (int i = 0; i < triMeshes.size(); i++) {
        triMeshes[i].printInfo();
    }

    /**
     * 配置メッシュの中で光源のものをリストに格納する
     */
    for (int i = 0; i < triMeshes.size(); i++) {
        if (triMeshes[i].is_light) {
            lightTriMeshIdxes.push_back(i);
        }
    }

    /**
     * 光源を選択する累積分布関数を作る。
     */
    for (int i = 0; i < lightTriMeshIdxes.size(); i++) {
        if (i == 0) {
            sample_lights_cdf.push_back(triMeshes[lightTriMeshIdxes[i]].A);
        } else {
            sample_lights_cdf.push_back(sample_lights_cdf[i - 1] +
                                        triMeshes[lightTriMeshIdxes[i]].A);
        }
        S_A += triMeshes[lightTriMeshIdxes[i]].A;
    }

    std::cout << "sample_lights_cdf" << std::endl;
    for (int i = 0; i < lightTriMeshIdxes.size(); i++) {
        sample_lights_cdf[i] = sample_lights_cdf[i] / S_A;
        std::cout << i << ":\t" << sample_lights_cdf[i] << std::endl;
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
