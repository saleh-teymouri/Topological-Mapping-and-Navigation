void triangular_landmark_generating (vector<Landmark>& landmarks, Mat image, int d, Vec3b c) {

	int id(0);
	double r = pow((pow(d, 2) + ((pow(d, 2))/4)), 0.5);

	for (int n1 = 0; n1 < (image.cols/d); ++n1) {
		for (int m1 = 0; m1 < (image.rows/r); ++m1) {
			if (image.at<Vec3b> (m1, n1) != c) {
				landmarks.push_back (Landmark(0+(d*n1), 0+(2*r*m1), id));
				id++;
			}
		}
	}

	for (int n2 = 0; n2 < (image.cols/d); ++n2) {
		for (int m2 = 0; m2 < (image.rows/r); ++m2) {
			if (image.at<Vec3b> (m2, n2) != c) {
				landmarks.push_back (Landmark((d/2)+(d*n2), (r)+(2*r*m2), id));
				id++;
			}
		}
	}
}

void read_landmarks (string str, vector<Landmark>& landmarks) {

	ifstream myfile(str);
	int x(0), y(0), id(0);
	char c;

	if (myfile.is_open()) {
		while (myfile >> c >> x >> c >> y >> c) {
			landmarks.push_back (Landmark(x, y, id));
			++id;
		}
		myfile.close();
	}
}

void center_of_blob(Mat & img, Point & pt, vector<Point> & contours) {

	bool do_break = false;
	for (int i = 0; i < img.cols; ++i) {
		for (int j = 0; j < img.rows; ++j) {
			if (img.at<Vec3b> (j, i) == Vec3b (255, 255, 255)) {
				floodFill (img, Point(i, j), Scalar (0, 150, 0), 0,
					Scalar (0, 0, 0), Scalar (255, 255, 255), 4);
				for (int i2 = 0; i2 < img.cols; ++i2) {
					for (int j2 = 0; j2 < img.rows; ++j2) {
						if (img.at<Vec3b> (j2, i2) == Vec3b (0, 150, 0)) {
							contours.push_back (Point(i2, j2));
						}
					}
				}

				int x(0), y(0), xbar(0), ybar(0), shifted_xbar(0), shifted_ybar(0);
				for (int it = 0; it < contours.size(); ++it) {				
					x += contours[it].x;
					y += contours[it].y;
				}

				xbar = (x / contours.size());
				ybar = (y / contours.size());
				pt = Point(xbar, ybar);
				do_break = true;
				break;
			}

			if (do_break) break;
		}

		if (do_break) break;
	}
}

void move_to_white_placement (Point pt, Mat & image, vector<Landmark> & landmarks,
		vector<Point> contours, int R, int theta, double start_angle, double end_angle) {

	Point Possible_temp;
	double d(10e6);
	for (int it = 0; it < contours.size(); ++it) {
		double distance = dist (pt, contours[it]);
		if (distance < d) {
		d = distance;
		Possible_temp = contours[it];
		}
	}

	landmarks.push_back (Landmark(Possible_temp.x, Possible_temp.y, landmarks.size()));
	LOS_circle_fill (image, Point(Possible_temp.x, Possible_temp.y),
		R, theta, start_angle, end_angle, obs_clr, vis_clr);
	if (nice_image_show)
		LOS_circle_fill3 (nice_image, Point(Possible_temp.x, Possible_temp.y),
			R, theta, start_angle, end_angle, obs_clr, color);
}

void make_valid_placement (Point pt, Point temp, Mat & image,
		vector<Landmark> & landmarks, int R, int theta, double start_angle, double end_angle) {

	Point Counter_Landmark(pt.x, pt.y);
	Point Possible_Landmark = Counter_Landmark;
	Point Possible_Landmark2 = Possible_Landmark;
	double big(1);

	if (abs(pt.x - temp.x) > abs(pt.y - temp.y)) 
		{ big = abs(pt.x - temp.x); }
	else { big = abs(pt.y - temp.y); }

	for (int n = 0; n <= abs(pt.x - temp.x) || n <= abs(pt.y - temp.y); ++n) {
		Counter_Landmark = (Point(
		(pt.x + round((double)(((double)n)*((double)(temp.x - pt.x)))/((double)big))),
		(pt.y + round((double)(((double)n)*((double)(temp.y - pt.y)))/((double)big)))));

		if (!boundary_check (Counter_Landmark, image, obs_clr)) {
			landmarks.push_back (Landmark(Possible_Landmark2.x, Possible_Landmark2.y,
				landmarks.size()));
			LOS_circle_fill (image, Point(Possible_Landmark2.x, Possible_Landmark2.y),
				R, theta, start_angle, end_angle, obs_clr, vis_clr);
			if (nice_image_show)
				LOS_circle_fill3 (nice_image, Point(Possible_Landmark2.x, Possible_Landmark2.y),
					R, theta, start_angle, end_angle, obs_clr, color);
			break;
		}

		else {
			Possible_Landmark2 = Possible_Landmark;
			Possible_Landmark = Counter_Landmark;
		}
	}
}

void Place_Landmark (Mat & image, vector<Landmark> & landmarks, double shift, double R,
					double theta, double start_angle, double end_angle) {

	vector<Point> contours;
	do {
		Point pt;
		Mat img = image.clone();
		contours.clear();
		center_of_blob(img, pt, contours);
		if (contours.size() > 0) {
			double rad_theta = theta * PI / 180.00;
			int shifted_xbar = pt.x - (shift * cos(rad_theta));
			int shifted_ybar = pt.y - (shift * sin(rad_theta));
			Point temp = Point(shifted_xbar, shifted_ybar);
			if (image.at<Vec3b> (temp.y, temp.x) != Vec3b (255, 255, 255)
				&& image.at<Vec3b> (pt.y, pt.x) != Vec3b (255, 255, 255))
				move_to_white_placement(pt, image, landmarks, contours,
					R, theta, start_angle, end_angle);

			else {

				if (validation_check (pt, temp, image, obs_clr)) {
					landmarks.push_back (Landmark(temp.x, temp.y,
						landmarks.size()));
					LOS_circle_fill (image, Point(temp.x, temp.y),
						R, theta, start_angle, end_angle, obs_clr, vis_clr);
					if (nice_image_show)
						LOS_circle_fill3 (nice_image, Point(temp.x, temp.y),
							R, theta, start_angle, end_angle, obs_clr, color);
				}

				else
					make_valid_placement (pt, temp, image,
						landmarks, R, theta, start_angle, end_angle);
			}
		}

		for (int it3 = 0; it3 < landmarks.size(); ++it3) {
			disp (img, landmarks[it3].coord, Vec3b (255, 0, 0), 1);
			if (nice_image_show)
				disp (nice_image, landmarks[it3].coord, Vec3b (0, 0, 255), 2);
		}

		if (nice_image_show) {
			imshow("nice", nice_image);
			char filename [50];
			sprintf (filename, "LPA/images/temp/%d_%d.png",iteration, sub_iteration);
			imwrite(filename, nice_image);
		}

//		imshow("Img", img);
//		waitKey(0);
		sub_iteration++;

	} while (contours.size() > 0);
}

void R2_LPA (vector<Landmark> & landmarks) {

	double delta_R = (double)(Rs - Rf) / (double)step_R;
	double shift(0), theta(0), start_angle(0), end_angle(360);

	for (int m = 0; m <= step_R; m++) {	cout << landmarks.size() << endl;
		double R = Rs - (m * delta_R);
		sub_iteration = 0;
		iteration++;
		Mat image = imread(image_name, CV_LOAD_IMAGE_COLOR);
		if (nice_image_show)
			nice_image = imread(nice_image_name, CV_LOAD_IMAGE_COLOR);
		for (int p = 0; p < landmarks.size(); ++p) {
			LOS_circle_fill (image, landmarks[p].coord, R, theta,
				start_angle, end_angle, obs_clr, vis_clr);
			if (nice_image_show)
				LOS_circle_fill3 (nice_image, landmarks[p].coord, R, theta,
					start_angle, end_angle, obs_clr, color);
		}

		if (nice_image_show) {
			for (int it3 = 0; it3 < landmarks.size(); ++it3)
				disp (nice_image, landmarks[it3].coord, Vec3b (0, 0, 255), 2);
		}

		Place_Landmark (image, landmarks, shift, R, theta, start_angle, end_angle);
	}
}

void SE2_LPA (vector<Landmark> & landmarks) {

	double delta_R = (double)(Rs - Rf) / (double)step_R;
	double delta_alpha = (double)(alpha_start - alpha_final) / (double)(2 * step_alpha);

	for (int h = 0; h <= (step_R + step_alpha) ; ++h) {	cout << landmarks.size() << endl;
		for (int n = 0; n <= h && n <= step_alpha; ++n) {
			int m = h - n;
			if (m <= step_R) {
				double R = Rs - (m * delta_R);
				double start_angle = (0 + rotation) + (n * delta_alpha);
				double end_angle = (alpha_start + rotation) - (n * delta_alpha);
				for (int theta = -180; theta <= 180; theta += delta_theta) {
					sub_iteration = 0;
					iteration++;
					double rad_theta = theta * PI / 180.00;
					Mat image = imread(image_name, CV_LOAD_IMAGE_COLOR);
					if (nice_image_show)
						nice_image = imread(nice_image_name, CV_LOAD_IMAGE_COLOR);
					for (int p = 0; p < landmarks.size(); ++p) {
						LOS_circle_fill (image, landmarks[p].coord, R, theta,
							start_angle, end_angle, obs_clr, vis_clr);
						if (nice_image_show)
							LOS_circle_fill3 (nice_image, landmarks[p].coord, R, theta,
								start_angle, end_angle, obs_clr, color);
					}

					if (nice_image_show) {
						for (int it3 = 0; it3 < landmarks.size(); ++it3)
							disp (nice_image, landmarks[it3].coord, Vec3b (0, 0, 255), 2);
					}

					Place_Landmark (image, landmarks, shift, R, theta, start_angle, end_angle);
				}
			}
		}
	}
}
