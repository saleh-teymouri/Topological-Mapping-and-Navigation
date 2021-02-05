/** Landmark Complex Construction with Search Problem Using Raster Scan **/
/* c++ -std=c++11 -g -I. -Idosl RS.cpp -lopencv_core -lopencv_highgui -lopencv_imgproc */

#include "1.Libraries.h"
#include "2.User_Defined_Classes.h"
#include "3.User_Defined_Functions.h"
#include "4.Exploration_Checks.h"
#include "5.Robot_Path_Planning.h"
#include "6.Landmark_Placement.h"

int main() {

//	vector<thread> threads;

	for (int p = 0; p < rbt; ++p) {
		A:
		int x = rand()%ref_image.cols;
		int y = rand()%ref_image.rows;
		if (ref_image.at<Vec3b> (y, x) == Vec3b(255, 255, 255))
			observations.push_back(dpoint (x, y, (((double)rand())/((double)RAND_MAX))*2*PI-PI,
								sensor_radius, sensor_angle, p));
		else
			goto A;
	}

	read_landmarks (text_name, landmarks);
	for (int a = 0; a < landmarks.size(); ++a) {
		obs_list[0].insert(landmarks[a].id);
	}

	Raster_scan ();

	cout << endl;

	cout << "C2_count: " << sc.C[2].size() << ", " << C2_size << endl;

	cout << endl;

	Mat IMG1(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
	Mat image1 = imread (image_name);
	sc.draw_simplicial_complex (IMG1, landmarks);

	imshow("IMG1", IMG1);
	waitKey (0);

return 0; }
