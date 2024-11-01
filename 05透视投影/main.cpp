#include <vector>
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
Vec3f light_dir(0.2,0.15,-1);
Vec3f camera(0,0,3);

// 该函数将一个 4x1 齐次坐标矩阵转换为一个 3x1 向量
// 通常用于将齐次坐标转为三维坐标，通过除以w分量进行归一化
// TODO:当最后一个分量为0，表示向量，不为0，表示坐标
Vec3f m2v(Matrix m) {
    // 列向量罢了
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

//3d-->4d
//添加一个1表示坐标
Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

//将物体x，y坐标(-1,1)转换到屏幕坐标(100,700)    1/8width~7/8width
//zbuffer(-1,1)转换到0~255
// 视口变换矩阵，用于从NDC坐标系到屏幕坐标系的转换
// xy是视口左下角的坐标，也就是最终渲染图像的原点位置；视口左下角通常对应屏幕坐标系原点(0,0)
/**
 * @param x 定义视口左下角在屏幕上的 x 轴向位置
 * @param y 定义视口左下角在屏幕上的 y 轴向位置
 * @attention xy通常是屏幕坐标的起始位置
 * @param w 屏幕宽度
 * @param h 屏幕高度
 * @brief 视口变换矩阵，用于从NDC坐标系到屏幕坐标系的转换
 * @return 视口变换矩阵
 */
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    //第4列表示平移信息
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;        // NDC在Z轴方向上[-1,1]，而深度缓冲区在[0,255]，所以要将Z轴坐标映射到[0,255]范围
    //对角线表示缩放信息
    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

/**
* @brief 分段绘制三角形_uv (屏幕空间顶点坐标，uv坐标，tga指针，亮度，zbuffer)
* @param t0/t1/t2 三角形屏幕空间顶点坐标
*/
void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage &image, float intensity, int *zbuffer) {
    if (t0.y==t1.y && t0.y==t2.y) return;
    //分割成两个三角形，t0<t1<t2
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
    if (t0.y>t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
    if (t1.y>t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }

    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++)
    {
        // 是否属于上半部分
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        //alpha表示当前点在t0与t2之间的比例，beta表示当前点在每一部分中间的比例
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
        //A表示t0与t2之间的点，B表示t0与t1之间的点
        //A是三角形左边的点，B是右边的点
        Vec3i A   =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B   = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
        Vec2i uvA =               uv0 +      (uv2-uv0)*alpha;
        Vec2i uvB = second_half ? uv1 +      (uv2-uv1)*beta : uv0 +      (uv1-uv0)*beta;
        //保证B在A的右边
        if (A.x > B.x) { std::swap(A, B); }// std::swap(uvA, uvB);}
        //用横坐标作为循环控制，对这一行进行着色
        for (int j=A.x; j<=B.x; j++) {
            //在AB之间横向的比例
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(float)(B.x-A.x);
            //插值计算当前点的坐标
            Vec3i   P = Vec3f(A) + Vec3f(B-A)*phi;
            Vec2i uvP =     uvA +   (uvB-uvA)*phi;
            if (P.x < width && P.y < height)
            {
                //计算当前zbuffer下标=P.x+P.y*width
                int idx = P.x+P.y*width;
                //当前点的z大于zbuffer信息，覆盖掉，并更新zbuffer
                if (zbuffer[idx]<P.z) {
                    zbuffer[idx] = P.z;
                    TGAColor color = model->diffuse(uvP);   // 从贴图中获取纹理颜色
                    image.set(P.x, P.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity));
                }
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

    //构造zbuffer
    zbuffer = new int[width*height];
    for (int i=0; i<width*height; i++) {
        //初始化zbuffer
        zbuffer[i] = std::numeric_limits<int>::min();
    }

    // TODO:面的填色有两种方法，这里只实现了三角形横线扫描填充绘制，还可以通过质心坐标插值进行填充
    //绘制
    {
        Matrix Projection = Matrix::identity(4);
        Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
        //投影矩阵[3][2]=-1/c，c为相机z坐标
        Projection[3][2] = -1.f/camera.z;

        TGAImage image(width, height, TGAImage::RGB);

        for (int i=0; i<model->nfaces(); i++) {
            std::vector<int> face = model->face(i);     // 返回第i个面的三个顶点索引
            Vec3i screen_coords[3];
            Vec3f world_coords[3];
            for (int j=0; j<3; j++) {
                Vec3f v = model->vert(face[j]);         // 对于第i个面的第j个顶点，返回顶点坐标，这里是世界坐标
                // 矩阵乘法顺序不可修改
                screen_coords[j] = m2v(ViewPort * Projection * v2m(v));
                world_coords[j]  = v;
            }
            //TODO:计算法向量的时候向量的顺序
            Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
            n.normalize();
            //计算光照
            float intensity = n*light_dir;
            intensity = std::min(std::abs(intensity),1.f);
            if (intensity>0) {
                Vec2i uv[3];
                for (int k=0; k<3; k++) {
                    uv[k] = model->uv(i, k);
                }
                //绘制三角形
                triangle(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
            }
        }
        //tga默认原点在左上，现改为左下
        image.flip_vertically();
        image.write_tga_file("output.tga");
    }
    //输出zbuffer
    {
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i + j * width], 1));
            }
        }
        zbimage.flip_vertically();
        zbimage.write_tga_file("zbuffer.tga");
    }
    delete model;
    delete [] zbuffer;
    return 0;
}

