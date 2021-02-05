#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <thread>
#include <mutex>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <dosl>
#include <armadillo>
#include "0.Hungarian_Algorithm.h"
#define PI 3.14159265

using namespace std;
using namespace cv;

/**************************************************************************************************/

bool SE2 = true;
string image_name = "BC.png";
string text_name = "BC.txt";
string nice_image_name = "BC(original).png";
Mat nice_image = imread(nice_image_name, CV_LOAD_IMAGE_COLOR);
bool nice_image_show = false;
double alpha_start(180), alpha_final(60), step_alpha(5), delta_theta(10), rotation(-90);
int Rs(120), Rf(20), step_R(20), shift(12);
Vec3b color = Vec3b(255, 200, 200);
Vec3b obs_clr = Vec3b(0, 0, 0);
Vec3b vis_clr = Vec3b(100, 100, 100);

/**************************************************************************************************/

Mat ref_image = imread(nice_image_name);
Mat PthHmlgy_image = imread(nice_image_name);
Mat ifl_image = imread(image_name);
Mat falseHolesIMG = imread(nice_image_name);
Mat sc_image(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
Mat rbt_image(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
Mat clusters_img(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
Mat homology_image(nice_image.rows, nice_image.cols, CV_8UC3, Scalar (255, 255, 255));
Mat simplex_image(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));

int rbt(4);
double sensor_radius(20.00), sensor_angle(PI/3);
double max_r(20.00), min_r(10.00), max_d(10.00), min_d(5.00);
int iteration(0), sub_iteration(0);
int expl_count(0), sc_count(0), ISW_counter(0);
int sc_number(1439);
		// M2: 545190 (75093), L457: 130128 (25945/NR:1531), BC: 236254 (39574/NR:1439) 
int factor, num_ISW, num_RW;
int jump_pos(0);
int animation_wait(10);
bool live_simplical_comlex_image = false; //4.Exploration_Checks.h:168
bool Show_image = true; //5.Robot_Path_Planning.h:131-150-255
int homology_step(0);
int C2_size(0);

double first_percent_complete(0.85);
double second_percent_complete(0.85);
double third_percent_complete(0.98);
double rt(1e6);

arma::Mat<double> x00;

vector<bool> HIW_finished(rbt, false);
vector<bool> RWISW_finished(rbt, false);

mutex sc_mtx;
mutex img_mtx;
mutex landmarks_mtx;
mutex observations_mtx;

ofstream myfile;
ofstream file;
ofstream file_ratio;
