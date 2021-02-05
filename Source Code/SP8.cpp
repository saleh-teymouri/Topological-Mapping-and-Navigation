/* c++ -std=c++11 -g -I. -Idosl SP8.cpp -lopencv_core -lopencv_highgui -lopencv_imgproc -lpthread -DARMA_DONT_USE_WRAPPER -lopenblas -llapack */

#include "1.Libraries.h"
#include "2.User_Defined_Classes.h"
#include "3.User_Defined_Functions.h"
#include "4.Exploration_Checks.h"
#include "5.Robot_Path_Planning.h"
#include "6.Landmark_Placement.h"

int main(int argc, char* argv[]) {

	factor = 1;//atof(argv[1]);
	num_ISW = 2;//atof(argv[2]);
	num_RW = 1;//atof(argv[3]);

	cout << factor << ", " << num_ISW << ", " << num_RW << endl;

	vector<thread> threads;

	for (int p = 0; p < rbt; ++p) {
		A:
		int x = rand()%ref_image.cols;
		int y = rand()%ref_image.rows;
		if (ifl_image.at<Vec3b> (y, x) == Vec3b(255, 255, 255))
			observations.push_back(dpoint (x, y, (((double)rand())/((double)RAND_MAX))*2*PI-PI,
								sensor_radius, sensor_angle, p));
		else
			goto A;
	}

//	observations.push_back(dpoint (300, 200, (((double)rand())/((double)RAND_MAX))*2*PI-PI,
//		sensor_radius, sensor_angle, 0));

	read_landmarks (text_name, landmarks);
	for (int a = 0; a < landmarks.size(); ++a) {
		obs_list[0].insert(landmarks[a].id);
	}

	for (int r = 0; r < rbt; ++r)
		threads.push_back (thread (move_robot, observations[r].rbt_id));

	sleep(1);
//	file.open("Results/file.txt");
//	file_ratio.open("Results/file_ratio.txt");
	int it = 0;

	Mat IMG1(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
	Mat IMG2(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
//	Mat IMG1 = imread(image_name, CV_LOAD_IMAGE_COLOR);
//	Mat IMG2 = imread(image_name, CV_LOAD_IMAGE_COLOR);

	int previous_expl_count(0);
	int previous_sc_count(0);

	while (!Thread_finished_check(RWISW_finished)) {

		img_mtx.lock();
		sc_image = 0.95*sc_image + 0.05*PthHmlgy_image;
		PthHmlgy_image = imread(nice_image_name);
		img_mtx.unlock();

		cout<<"sc.C[2].size(): "<<sc.C[2].size()<<", "<<C2_size<<", expl_count: "<<expl_count
			<<", percentage: "<<(double(C2_size)/sc_number)*100<<"%"<<"\r"<<flush;

		if (it%1000 == 0) {

//			file<<expl_count<<"\t"<<(double(sc.C[2].size())/sc_number)*100<<endl;

//			sc.draw_simplicial_complex (IMG1, landmarks);

			rt = double(C2_size-previous_sc_count)/double(expl_count-previous_expl_count);
			cout << "\nit: " << it << "\trt: " << rt << endl;
				//"\t" << white_pixel_count(IMG1) <<endl;

			previous_expl_count = expl_count;
			previous_sc_count = C2_size;
		}

		it++;

		imshow ("img0", sc_image);
		imshow ("simplex_image", simplex_image);
		waitKey(animation_wait);

	}

	RWISW_finished.assign(rbt, false);

	cout << endl;
	cout << "RW&ISW sc/obs, " << it << ": " <<
		double(C2_size-previous_sc_count)/double(expl_count-previous_expl_count)
		<< endl;
	sc.draw_simplicial_complex (IMG1, landmarks);
	imshow("IMG1", IMG1);
	waitKey (2e3);

//	IMG2 = IMG1.clone();
//	cout << "white pixel left: " << white_pixel_count(IMG2) << endl;

//	file<<endl;
//	file_ratio<<endl;

	while ((double(C2_size)/sc_number) < third_percent_complete) { //white_pixel_count(IMG2) > 100

		homology_step++;
		sc.Homology ();

//		myfile.open("Results/res.txt", ios_base::app);
//		myfile<<factor<<"\t"<<num_ISW<<"\t"<<num_RW<<"\t"<<expl_count<<endl;

		while (!Thread_finished_check(HIW_finished)) {

			for (int r = 0; r < rbt; ++r) {
				if (!HIW_finished[r] && !observations[r].Assigned) {
					threads.push_back (thread (
						updated_Homology_Informed_Walk, observations[r].rbt_id));
				}
			}

			cout<<"sc.C[2].size(): "<<sc.C[2].size()<<", "<<C2_size<<", expl_count: "
				<<expl_count<<", percentage: "<<(double(C2_size)/sc_number)*100
				<<"%"<<"\r"<<flush;

			if (it%1000 == 0) {

/*				file_ratio<<it<<"\t"<<
					double(sc.C[2].size()-previous_sc_count)/
					double(expl_count-previous_expl_count)*100
					<<endl;
*/
				previous_expl_count = expl_count;
				previous_sc_count = C2_size;
			}

			it++;

			img_mtx.lock();
			sc_image = 0.95*sc_image + 0.05*PthHmlgy_image;
			PthHmlgy_image = imread(nice_image_name);
			img_mtx.unlock();
			imshow ("img0", sc_image);
			imshow ("simplex_image", simplex_image);
			waitKey (animation_wait);

		}

		HIW_finished.assign(rbt, false);

		for (int r = 0; r < rbt; ++r)
			observations[r].Assigned = false;

		cout << endl;
		cout << "HIW sc/obs, " << it << ": " <<
			double(C2_size-previous_sc_count)/double(expl_count-previous_expl_count) << endl; 
				//"\t" << white_pixel_count(IMG2) << endl;

		previous_expl_count = expl_count;
		previous_sc_count = C2_size;

		sc.draw_simplicial_complex (IMG2, landmarks);
		imshow("IMG2", IMG2);
		waitKey (2e3);

	}

//	cout << "white pixel left: " << white_pixel_count(IMG2) << endl;
	waitKey (0);

return 0; }
