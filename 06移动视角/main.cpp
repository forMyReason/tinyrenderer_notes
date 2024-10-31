#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <algorithm>
const int width  = 800;
const int height = 800;
const int depth  = 255;

Model *model = NULL;
int *zbuffer = NULL;
Vec3f light_dir = Vec3f(0,-1,-1).normalize();
//摄像机位置
Vec3f eye(2,1,3);
//焦点位置
Vec3f center(0,0,1);

//视角矩阵，用于将(-1,1),(-1,1),(-1,1)映射到(1/8w,7/8w),(1/8h,7/8h),(0,255)
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

//TODO:怪不得之前笔试的时候让推导lookat矩阵
//朝向矩阵，变换矩阵
//更改摄像机视角=更改物体位置和角度，操作为互逆矩阵
//摄像机变换是先旋转再平移，所以物体需要先平移后旋转，且都是逆矩阵
//正好对应之前GAMES中的view矩阵
Matrix lookat(Vec3f eye, Vec3f center, Vec3f up)
{
    Vec3f z = (eye - center).normalize();
    Vec3f x = (up ^ z).normalize();
    Vec3f y = (z ^ x).normalize();
    Matrix rotation = Matrix::identity(4);
    Matrix translation = Matrix::identity(4);
    //矩阵的第四列是用于平移的。因为观察位置从原点变为了center，所以需要将物体平移-center
    for (int i = 0; i < 3; i++) {
        rotation[i][3] = -center[i];
    }
    //正交矩阵的逆 = 正交矩阵的转置
    //矩阵的第一行即是现在的x
    //矩阵的第二行即是现在的y
    //矩阵的第三行即是现在的z
    //***矩阵的三阶子矩阵是当前视线旋转矩阵的逆矩阵***
    for (int i = 0; i < 3; i++) {
        rotation[0][i] = x[i];
        rotation[1][i] = y[i];
        rotation[2][i] = z[i];
    }
    //这样乘法的效果是先平移物体，再旋转
    Matrix res = rotation*translation;
    return res;
}

/**
 * @brief 三角形光栅化 (屏幕空间顶点坐标，灯光强度，uv坐标，顶点到摄像机距离，tga指针，zbuffer)
 * @param t0/t1/t2 三角形各个点屏幕空间坐标
 * @param ity0/ity1/ity2 三角形各个点灯光强度
 * @param uv0/uv1/uv2 三角形各个点uv坐标
 * @param dis0/dis1/dis2 三角形各个点到摄像机距离
 */
void triangle(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, Vec2i uv0, Vec2i uv1, Vec2i uv2,float dis0, float dis1, float dis2, TGAImage &image, int *zbuffer) {
    if (t0.y==t1.y && t0.y==t2.y) return;

    if (t0.y>t1.y)   { std::swap(t0, t1); std::swap(ity0, ity1); std::swap(uv0, uv1); }
    if (t0.y > t2.y) { std::swap(t0, t2); std::swap(ity0, ity2); std::swap(uv0, uv2); }
    if (t1.y > t2.y) { std::swap(t1, t2); std::swap(ity1, ity2); std::swap(uv1, uv2); }

    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++)
    {
        // 是否属于上半部分
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;

        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height;

        // 屏幕空间坐标插值，灯光强度插值，uv坐标插值，顶点到摄像机距离插值
        Vec3i A    =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B    = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0)*beta;
        // 因为加入了对于光照强度的插值，所以最终渲染的结果会更加真实
        float ityA =               ity0 +   (ity2-ity0)*alpha;
        float ityB = second_half ? ity1 +   (ity2-ity1)*beta : ity0 + (ity1-ity0)*beta;
        Vec2i uvA = uv0 + (uv2 - uv0) * alpha;
        Vec2i uvB = second_half ? uv1 + (uv2 - uv1) * beta : uv0 + (uv1 - uv0) * beta;
        float disA = dis0 + (dis2 - dis0) * alpha;
        float disB = second_half ? dis1 + (dis2 - dis1) * beta : dis0 + (dis1 - dis0) * beta;

        // 确保A在B的左边
        if (A.x>B.x) { std::swap(A, B); std::swap(ityA, ityB); }

        for (int j=A.x; j<=B.x; j++) {
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(B.x-A.x);

            // 插值计算当前绘制点的屏幕空间坐标，灯光强度，uv坐标，顶点到摄像机距离
            // 这样每一个点就不是同一个光照强度了，而是根据插值计算出来的，于是消除了光照的跳变
            Vec3i    P = Vec3f(A) +  Vec3f(B-A)*phi;
            float ityP =    ityA  + (ityB-ityA)*phi;
            ityP = std::min(1.f, std::abs(ityP)+0.01f);     // 当前插值点的亮度需要保证在0-1之间
            Vec2i uvP = uvA + (uvB - uvA) * phi;
            float disP = disA + (disB - disA) * phi;        // 计算当前插值点的深度disP

            // 当前插值点在zbuffer中的索引，如果不在，跳过该点
            int idx = P.x+P.y*width;
            if (P.x>=width||P.y>=height||P.x<0||P.y<0) continue;

            // 如果当前插值点的深度比zbuffer中的深度小，更新zbuffer和image
            if (zbuffer[idx]<P.z)
            {
                zbuffer[idx] = P.z;
                TGAColor color = model->diffuse(uvP);
                // 根据亮度和深入调整颜色
                image.set(P.x, P.y, TGAColor(color.bgra[2], color.bgra[1], color.bgra[0])*ityP*(20.f/std::pow(disP,2.f)));
                //image.set(P.x, P.y, TGAColor(255,255,255)* ityP);
            }
        }
    }
}

int main(int argc, char** argv) {
    //读取模型
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    //构造zbuffer并初始化
    zbuffer = new int[width*height];
    for (int i=0; i<width*height; i++) {
        zbuffer[i] = std::numeric_limits<int>::min();
    }
    //绘制模型
    {
        //模型变换矩阵
        Matrix ModelView  = lookat(eye, center, Vec3f(0,1,0));

        //透视矩阵
        Matrix Projection = Matrix::identity(4);
        Projection[3][2] = -1.f / (eye - center).norm();

        //视角矩阵
        Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);

        /*std::cerr << ModelView << std::endl;
        std::cerr << Projection << std::endl;
        std::cerr << ViewPort << std::endl;
        Matrix z = (ViewPort*Projection*ModelView);
        std::cerr << z << std::endl;*/

        TGAImage image(width, height, TGAImage::RGB);
        for (int i=0; i<model->nfaces(); i++) {
            std::vector<int> face = model->face(i);
            Vec3i screen_coords[3];
            float intensity[3];
            float distance[3];
            
            for (int j=0; j<3; j++) {
                Vec3f v = model->vert(face[j]);
                Matrix m_v = ModelView* Matrix(v);
                screen_coords[j] =  Vec3f(ViewPort*Projection* m_v);
                intensity[j] = model->norm(i, j)*light_dir;
                Vec3f new_v = Vec3f(m_v);
                distance[j] = std::pow((std::pow(new_v.x - eye.x,2.0f)+ std::pow(new_v.y - eye.y, 2.0f)+ std::pow(new_v.z - eye.z, 2.0f)),0.5f);
            }
            Vec2i uv[3];
            for (int k = 0; k < 3; k++) {
                uv[k] = model->uv(i, k);
            }
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0], intensity[1], intensity[2], uv[0], uv[1], uv[2], distance[0], distance[1], distance[2], image, zbuffer);
        }
        image.flip_vertically();
        image.write_tga_file("output.tga");
    }
    //输出zbuffer图像
    {
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width]));
            }
        }
        zbimage.flip_vertically();
        zbimage.write_tga_file("zbuffer.tga");
    }
    delete model;
    delete [] zbuffer;
    return 0;
}

