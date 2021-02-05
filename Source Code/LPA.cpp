/**Landmark Placement Algorithm**/
/* c++ -std=c++11 -g -I. -Idosl LPA.cpp -lopencv_core -lopencv_highgui -lopencv_imgproc */

#include "1.Libraries.h"
#include "2.User_Defined_Classes.h"
#include "3.User_Defined_Functions.h"
#include "4.Exploration_Checks.h"
#include "5.Robot_Path_Planning.h"
#include "6.Landmark_Placement.h"

int main() {

	vector<Landmark> landmarks;

	double delta_R = (double)(Rs - Rf) / (double)step_R;
	double delta_alpha = (double)(alpha_start - alpha_final) / (double)(2 * step_alpha);

	Mat image0 = imread(image_name, CV_LOAD_IMAGE_COLOR);

	int rows, cols, dir(45);

	rows = image0.rows;
	cols = image0.cols;

	Mat image1(rows, cols, CV_8UC3, Scalar (255, 255, 255));
	Mat image2(rows, cols, CV_8UC3, Scalar (255, 255, 255));

//	landmarks.push_back (Landmark(160, 280, 0));

	R2_LPA (landmarks);
	if (SE2) SE2_LPA (landmarks);

	ofstream myfile;
	char filename [50];
	sprintf (filename, "LPA/results/textfiles/BC(shift-%d, radius-%d).txt", shift, Rs);
	myfile.open (filename);

	for (int a = 0; a < landmarks.size(); ++a) {
		disp (landmarks[a].coord, image1, Vec3b (0, 0, 255), 0);
		LOS_circle_fill (image0, landmarks[a].coord, Rf, 180,
			(dir - (alpha_final/2)), (dir + (alpha_final/2)), obs_clr, Vec3b(10, 10, 10));
		myfile << landmarks[a].coord << endl;
	}

	myfile.close();

	cout << landmarks.size() << endl;

	imshow(filename, image1);
	imshow("Image0", image0);
	waitKey(0);

return 0; }
