double dist (Point P1, Point P2) {

	double xd = (P1.x)-(P2.x);
	double yd = (P1.y)-(P2.y);

return sqrt (xd*xd+yd*yd); }

double dist2 (Point P1, Point P2, double d) {

	double xd = (P1.x)-(P2.x);
	double yd = (P1.y)-(P2.y);

	if (xd > d || yd > d)
		return 10e6;

	else
		return sqrt (xd*xd+yd*yd);
}

int fact(int n) {
 
    int res = 1; 
    for (int i = 2; i <= n; i++) 
        res = res * i; 
    return res; 
}

int nCr (int n, int r) { 
return fact(n) / (fact(r) * fact(n - r)); }

void disp (Mat & image, Point n, Vec3b m, int r=1) {

	for (int i = -r; i <= r; ++i) {
		for (int i2 = -r; i2 <= r; ++i2)
			if (n.y+i >= 0 && n.y+i < image.rows && n.x+i2 >= 0 && n.x+i2 < image.cols)
				image.at<Vec3b> ((n.y)+i, (n.x)+i2) = m;
	}
}

void circle_coords (Point center, int r, double start_angle, double end_angle,
					vector<Point> & circle_pts) {

	double step_theta = 1.00 / ((double)r);
	double rad_start_angle = start_angle * PI / 180.00;
	double rad_end_angle = end_angle * PI / 180.00;

	if (start_angle < end_angle) {
		for (double theta = rad_start_angle; theta <= rad_end_angle; theta += step_theta)
			circle_pts.push_back(Point(round(center.x + (r*cos(theta))),
				round(center.y + (r*sin(theta)))));
	}

	if (start_angle > end_angle) {
		for (double theta = rad_start_angle; theta <= rad_end_angle
		|| theta <= (rad_end_angle + 2*PI) || theta <= (rad_end_angle - 2*PI);theta += step_theta)
			circle_pts.push_back(Point(round(center.x + (r*cos(theta))),
				round(center.y + (r*sin(theta)))));
	}
}

bool quadrant_check (double theta_start, double theta_end, double theta_counter) {

	while (theta_counter > theta_end)
		theta_counter -= 2*PI;

	while (theta_counter < theta_start)
		theta_counter += 2*PI;

	if (theta_counter >= theta_start && theta_counter < theta_end)
		return true;

	else
		return false;
}

void midPoint_circle_coords (Point center, int r, double start_angle, double end_angle,
							vector<Point> & circle_pts) {

	Point first_midPoint = Point(center.x+r*cos(PI/4), center.y+r*sin(PI/4));
	Point second_midPoint = Point(center.x+r*cos(3*PI/4), center.y+r*sin(3*PI/4));
	Point third_midPoint = Point(center.x+r*cos(-3*PI/4), center.y+r*sin(-3*PI/4));
	Point fourth_midPoint = Point(center.x+r*cos(-PI/4), center.y+r*sin(-PI/4));

	double start_theta = start_angle * PI / 180.00;
	double end_theta = end_angle * PI / 180.00;

	Point start_Point = Point(center.x+r*cos(start_theta), center.y+r*sin(start_theta));
	Point end_Point = Point(center.x+r*cos(end_theta), center.y+r*sin(end_theta));

	double rad_start_angle = atan2(start_Point.y - center.y, start_Point.x - center.x);
	double rad_end_angle = atan2(end_Point.y - center.y, end_Point.x - center.x);

	Point counter_Point = Point(start_Point.x, start_Point.y);

	double counter_angle = atan2(counter_Point.y - center.y, counter_Point.x - center.x);

	if (rad_start_angle > rad_end_angle)
		rad_end_angle += 2*PI;

	if (rad_start_angle < 0) {
		rad_end_angle += 2*PI;
		rad_start_angle += 2*PI; 
	}

	while (quadrant_check(rad_start_angle, rad_end_angle, counter_angle)) {

		if (quadrant_check(-PI/4, PI/4, counter_angle)) {
			for (double yi = counter_Point.y; yi <= first_midPoint.y+1; ++yi) {
				double x = abs(pow((pow(r, 2) - pow(yi-center.y, 2)), 0.5)) + center.x;
				circle_pts.push_back (Point(x, yi));
				counter_Point = Point (x, yi);
				counter_angle = atan2(counter_Point.y - center.y, counter_Point.x - center.x);
				if (!quadrant_check(rad_start_angle, rad_end_angle, counter_angle))
					break;
			}
		}

		if (quadrant_check(PI/4, 3*PI/4, counter_angle)) {
			for (double xi = counter_Point.x; xi >= third_midPoint.x-1; --xi) {
				double y = abs(pow((pow(r, 2) - pow(xi-center.x, 2)), 0.5)) + center.y;
				circle_pts.push_back (Point(xi, y));
				counter_Point = Point (xi, y);
				counter_angle = atan2(counter_Point.y - center.y, counter_Point.x - center.x);
				if (!quadrant_check(rad_start_angle, rad_end_angle, counter_angle))
					break;
			}
		}

		if (quadrant_check(3*PI/4, 5*PI/4, counter_angle)) {
			for (double yi = counter_Point.y; yi >= third_midPoint.y-1; --yi) {
				double x = -abs(pow((pow(r, 2) - pow(yi-center.y, 2)), 0.5)) + center.x;
				circle_pts.push_back (Point(x, yi));
				counter_Point = Point (x, yi);
				counter_angle = atan2(counter_Point.y - center.y, counter_Point.x - center.x);
				if (!quadrant_check(rad_start_angle, rad_end_angle, counter_angle))
					break;
			}
		}
			
		if (quadrant_check(5*PI/4, 7*PI/4, counter_angle)) {
			for (double xi = counter_Point.x; xi <= first_midPoint.x+1; ++xi) {
				double y = -abs(pow((pow(r, 2) - pow(xi-center.x, 2)), 0.5)) + center.y;
				circle_pts.push_back (Point(xi, y));
				counter_Point = Point (xi, y);
				counter_angle = atan2(counter_Point.y - center.y, counter_Point.x - center.x);
				if (!quadrant_check(rad_start_angle, rad_end_angle, counter_angle))
					break;
			}
		}
	}
}

void line_coords (Point pt1, Point pt2, vector<Point> & line_pts) {

	vector<Point> pts;
	double big(1);

	if (abs(pt2.x - pt1.x) > abs(pt2.y - pt1.y)) { big = abs(pt2.x - pt1.x); }
	else { big = abs(pt2.y - pt1.y); }

	for (double n = 0; n <= abs(pt2.x - pt1.x) || n <= abs(pt2.y - pt1.y); ++n)
		pts.push_back (Point(
		(pt1.x + round((double)(((double)n)*((double)(pt2.x - pt1.x)))/((double)big))),
		(pt1.y + round((double)(((double)n)*((double)(pt2.y - pt1.y)))/((double)big)))));
}

void myellipse (Mat &image, Point pt, Size axes, double angle, double startangle, double endangle, 
		const Scalar & color, int t) {

	if (startangle > endangle) 
		ellipse (image, pt, axes, angle + 180,  -180 - endangle, 180 - startangle, color, t);
	else
		ellipse (image, pt, axes, angle, startangle, endangle, color, t);

}

void change_color (Mat &image, Vec3b c1, Vec3b c2) {

	for (int i = 0; i < image.cols; ++i)
		for (int j = 0; j < image.rows; ++j)
			if (image.at<Vec3b> (j, i) == c1)
				image.at<Vec3b> (j, i) = c2;
}

void draw_robot (Mat &image0, dpoint observations) {

	double d = observations.sensor_radius;
	double angle_range = observations.sensor_angle;

	disp(image0, observations.coord, Vec3b(0, 0, 255));
	Point line_point = Point(observations.coord.x + (d/2)*cos(observations.angle),
		observations.coord.y + (d/2)*sin(observations.angle));
	Point line_point1 = Point(
		observations.coord.x + d*cos(observations.angle-angle_range/2),
		observations.coord.y + d*sin(observations.angle-angle_range/2));
	Point line_point2 = Point(
		observations.coord.x + d*cos(observations.angle+angle_range/2),
		observations.coord.y + d*sin(observations.angle+angle_range/2));
	line (image0, observations.coord, line_point, Scalar(0, 0, 255), 2);
	line (image0, observations.coord, line_point1, Scalar(0, 255, 0), 1);
	line (image0, observations.coord, line_point2, Scalar(0, 255, 0), 1);
	myellipse (image0, observations.coord, Size(d, d), 0,
		(observations.angle-angle_range/2)*180/PI,
		(observations.angle+angle_range/2)*180/PI, Scalar(0, 255, 0), 1);
}

void LOS_circle_fill (Mat &image1, Point center, int r, double theta,
					double start_angle, double end_angle, Vec3b v1, Vec3b v2) {

	start_angle += theta;
	end_angle += theta;
	vector<Point> circle_pts;
	circle_coords (center, r, start_angle, end_angle, circle_pts);

	for (int i = 0; i < circle_pts.size(); ++i) {
		Point line_pts;
		int big(1);

		if (abs(center.x - circle_pts[i].x) > abs(center.y - circle_pts[i].y)) 
			{ big = abs(center.x - circle_pts[i].x); }
		else { big = abs(center.y - circle_pts[i].y); }

		bool do_break = false;
		for (int n = 0; n <= abs(center.x - circle_pts[i].x)
				|| n <= abs(center.y - circle_pts[i].y); ++n) {
			line_pts = (Point(
				(center.x + round((double)(((double)n)*((double)
				(circle_pts[i].x - center.x)))/((double)big))),
				(center.y + round((double)(((double)n)*((double)
				(circle_pts[i].y - center.y)))/((double)big)))));
			for (int r = 0; r <= 1; ++r) {
				for (int r2 = 0; r2 <= 1; ++r2) {
					if (line_pts.y+r < image1.rows && line_pts.y+r >= 0
						&& line_pts.x+r2 < image1.cols && line_pts.x+r2 >= 0
						&& image1.at<Vec3b> (line_pts.y+r, line_pts.x+r2) != v1) {
						image1.at<Vec3b> (line_pts.y+r, line_pts.x+r2) = v2;
					}

					else {
						do_break = true;
						break;
					}
				}
				if (do_break) break;
			}
			if (do_break) break;
		}
	}
}

Vec3b Combined_Pixels (Mat image0, Point pts, Vec3b c1) {

	Vec3b c2 = image0.at<Vec3b> (pts.y, pts.x);

	if (c1 == Vec3b(0, 0, 0) && c2 == Vec3b(0, 0, 0)) {
		return Vec3b(0, 0, 0);
	}

	if (c1 == Vec3b(0, 0, 0) && c2 == Vec3b(128, 128, 128)) {
		return Vec3b(0, 0, 0);
	}

	if (c1 == Vec3b(0, 0, 0) && c2 == Vec3b(255, 255, 255)) {
		return Vec3b(0, 0, 0);
	}

	if (c1 == Vec3b(128, 128, 128) && c2 == Vec3b(0, 0, 0)) {
		return Vec3b(0, 0, 0);
	}

	if (c1 == Vec3b(128, 128, 128) && c2 == Vec3b(128, 128, 128)) {
		return Vec3b(0, 0, 0);
	}

	if (c1 == Vec3b(128, 128, 128) && c2 == Vec3b(255, 255, 255)) {
		return Vec3b(128, 128, 128);
	}

	if (c1 == Vec3b(255, 255, 255) && c2 == Vec3b(0, 0, 0)) {
		return Vec3b(0, 0, 0);
	}

	if (c1 == Vec3b(255, 255, 255) && c2 == Vec3b(128, 128, 128)) {
		return Vec3b(128, 128, 128);
	}

	if (c1 == Vec3b(255, 255, 255) && c2 == Vec3b(255, 255, 255)) {
		return Vec3b(255, 255, 255);
	}
}

void LOS_circle_fill2 (Mat& image1, Point center, int r, double theta,
					double start_angle, double end_angle) {

	Mat image0 = image1.clone();
	start_angle += theta;
	end_angle += theta;
	vector<Point> circle_pts;
	circle_coords (center, r, start_angle, end_angle, circle_pts);

	for (int i = 0; i < circle_pts.size(); ++i) {
		Point line_pts;
		int big(1);

		if (abs(center.x - circle_pts[i].x) > abs(center.y - circle_pts[i].y)) 
			{ big = abs(center.x - circle_pts[i].x); }
		else { big = abs(center.y - circle_pts[i].y); }

		bool do_break = false;
		for (int n = 0; n <= abs(center.x - circle_pts[i].x)
				|| n <= abs(center.y - circle_pts[i].y); ++n) {
			line_pts = (Point(
				(center.x + round((double)(((double)n)*((double)
				(circle_pts[i].x - center.x)))/((double)big))),
				(center.y + round((double)(((double)n)*((double)
				(circle_pts[i].y - center.y)))/((double)big)))));
			for (int r = 0; r <= 1; ++r) {
				for (int r2 = 0; r2 <= 1; ++r2) {
					if (line_pts.y+r < image1.rows && line_pts.y+r >= 0
					&& line_pts.x+r2 < image1.cols && line_pts.x+r2 >= 0
					&& image1.at<Vec3b> (line_pts.y+r, line_pts.x+r2) != obs_clr) {
						image1.at<Vec3b> (line_pts.y+r, line_pts.x+r2) =
						Combined_Pixels (image0, line_pts, Vec3b(128, 128, 128));
					}

					else {
						do_break = true;
						break;
					}			
				}
				if (do_break) break;
			}	
			if (do_break) break;
		}
	}
}

Vec3b Combined_Pixels2 (Mat image0, Point pts, Vec3b c1) {

	Vec3b c2 = image0.at<Vec3b> (pts.y, pts.x);

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 255, 255)) {
		return Vec3b(255, 200, 200);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 200, 200)) {
		return Vec3b(255, 150, 150);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 150, 150)) {
		return Vec3b(255, 100, 100);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 100, 100)) {
		return Vec3b(255, 50, 50);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 50, 50)) {
		return Vec3b(255, 0, 0);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 0, 0)) {
		return Vec3b(255, 0, 0);
	}

/*	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 255, 255)) {
		return Vec3b(255, 200, 200);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 200, 200)) {
		return Vec3b(255, 180, 180);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 180, 180)) {
		return Vec3b(255, 160, 160);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 160, 160)) {
		return Vec3b(255, 140, 140);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 140, 140)) {
		return Vec3b(255, 120, 120);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 120, 120)) {
		return Vec3b(255, 100, 100);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 100, 100)) {
		return Vec3b(255, 80, 80);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 80, 80)) {
		return Vec3b(255, 60, 60);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 60, 60)) {
		return Vec3b(255, 40, 40);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 40, 40)) {
		return Vec3b(255, 20, 20);
	}

	if (c1 == Vec3b(255, 200, 200) && c2 == Vec3b(255, 20, 20)) {
		return Vec3b(255, 0, 0);
	}*/

	if (c2 == Vec3b(10, 10, 10)) {
		return Vec3b(10, 10, 10);
	}
}

void LOS_circle_fill3 (Mat &image1, Point center, int r, double theta,
					double start_angle, double end_angle, Vec3b v1, Vec3b v2) {

	Mat image0 = image1.clone();
	start_angle += theta;
	end_angle += theta;
	vector<Point> circle_pts;
	circle_coords (center, r, start_angle, end_angle, circle_pts);

	Mat contour = imread(nice_image_name, CV_LOAD_IMAGE_COLOR);

	for (int i = 0; i < circle_pts.size(); ++i) {
		Point line_pts;
		int big(1);

		if (abs(center.x - circle_pts[i].x) > abs(center.y - circle_pts[i].y)) 
			{ big = abs(center.x - circle_pts[i].x); }
		else { big = abs(center.y - circle_pts[i].y); }

		bool do_break = false;
		for (int n = 0; n <= abs(center.x - circle_pts[i].x)
				|| n <= abs(center.y - circle_pts[i].y); ++n) {
			line_pts = (Point(
				(center.x + round((double)(((double)n)*((double)
				(circle_pts[i].x - center.x)))/((double)big))),
				(center.y + round((double)(((double)n)*((double)
				(circle_pts[i].y - center.y)))/((double)big)))));
			for (int r = 0; r <= 0; ++r) {
				for (int r2 = 0; r2 <= 0; ++r2) {
					if (line_pts.y+r < image1.rows && line_pts.y+r >= 0
					&& line_pts.x+r2 < image1.cols && line_pts.x+r2 >= 0
					&& image1.at<Vec3b> (line_pts.y+r, line_pts.x+r2) != v1) {
						image1.at<Vec3b> (line_pts.y+r, line_pts.x+r2) =
							Combined_Pixels2 (image0, line_pts, v2);
						contour.at<Vec3b> (line_pts.y+r, line_pts.x+r2) = v2;
					}

					else {
						do_break = true;
						break;
					}		
				}

				if (do_break) break;
			}

			if (do_break) break;
		}
	}

	for (int j0 = 0; j0 < contour.rows; ++j0) {
		for(int i0 = 0; i0 < contour.cols; ++i0) {
			Point pt = Point(i0, j0);
			if (contour.at<Vec3b> (pt.y, pt.x) == Vec3b(255, 255, 255)) {
				bool do_break2 = false;
				int counter(0);
				for (int r = -1; r <= 1; ++r) {
					for (int r2 = -1; r2 <= 1; ++r2) {
						if (pt.y+r < contour.rows && pt.y+r >= 0
						&& pt.x+r2 < contour.cols && pt.x+r2 >= 0
						&& contour.at<Vec3b> (pt.y+r, pt.x+r2) != v1) {
							if (contour.at<Vec3b> (pt.y+r, pt.x+r2) == v2) {
								counter++;
							}
						}
					}
				}

				if (counter >= 5) {
					contour.at<Vec3b> (pt.y, pt.x) = v2;
					image1.at<Vec3b> (pt.y, pt.x) =
						Combined_Pixels2 (image0, pt, v2);
				}
			}
		}
	}

	for (int j2 = 0; j2 < contour.rows; ++j2) {
		for(int i2 = 0; i2 < contour.cols; ++i2) {
			Point pt = Point(i2, j2);
			if (contour.at<Vec3b> (pt.y, pt.x) == v2) {
				for (int r = -1; r <= 1; ++r) {
					for (int r2 = -1; r2 <= 1; ++r2) {
						if (pt.y+r < contour.rows && pt.y+r >= 0
						&& pt.x+r2 < contour.cols && pt.x+r2 >= 0
						&& contour.at<Vec3b> (pt.y+r, pt.x+r2) != v1) {
							if (contour.at<Vec3b> (pt.y+r, pt.x+r2) == Vec3b(255, 255, 255)) {
								nice_image.at<Vec3b> (pt.y, pt.x) = Vec3b(10, 10, 10);
							}
						}

						else
							nice_image.at<Vec3b> (pt.y, pt.x) = Vec3b(10, 10, 10);
					}
				}
			}
		}
	}
}

void left_or_right (dpoint observations, Landmark landmark, int &m) {

	double betha = atan2 ((landmark.coord.y - observations.coord.y),
						(landmark.coord.x - observations.coord.x));

	/** m = 1 : Right **/
	/** m = -1 : Left **/

	while (observations.angle > PI)
		observations.angle -= 2*PI;
	while (observations.angle < -PI)
		observations.angle += 2*PI;

	if (observations.angle * betha >= 0) {
		if (observations.angle >= betha)
			m = 1;
		else
			m = -1;
	}

	else {

		if ((observations.angle >= 0 && observations.angle < (PI/2))
		|| (observations.angle < 0 && observations.angle < (-PI/2)))
			m = 1;

		else
			m = -1;
	}
}

int locate_Jump (vector<double> data) {

	double sum(0), mean(0), standardDeviation(0);

	for(int i = 0; i < data.size(); ++i) 
		sum += data[i];

	mean = sum/data.size();

	for(int i = 0; i < data.size(); ++i)
		standardDeviation += pow(data[i] - mean, 2);

	double std_dvn = sqrt(standardDeviation / data.size());

	vector<double> Diff;
	for (int it = 0; it < data.size(); ++it) {
		if (it != 0) 
			Diff.push_back(pow(data[it], 2) - pow(data[it-1], 2));
		else
			Diff.push_back(0);
	}

	double data_max(-1e6), data_min(1e6);

	for (int it = 0; it < data.size(); ++it) {
		if (data[it] < data_min)
			data_min = data[it];

		if (data[it] > data_max)
			data_max = data[it];
	}

	for (int it = Diff.size()-1; it >= 0; --it) {
		vector<double> avg;
		double sum(0), mean(0);
		for (int it2 = it; it2 > it - (data.size()/1000) && it2 >= 0; --it2) {
			avg.push_back(Diff[it2]);
			sum += Diff[it2];
		}

		mean = sum/avg.size();
		if (mean < (std_dvn/1000)) return it;
	}

return 0; }

bool Thread_finished_check (vector<bool> finished) {

	for (int i = 0; i < finished.size(); ++i)
		if (finished[i] == false) return false;

return true;}

int white_pixel_count (Mat IMG) {

	int white_pixel(0);
	for (int j = 0; j < IMG.rows; ++j)
		for (int i = 0; i < IMG.cols; ++i)
			if (IMG.at<Vec3b> (j, i) == Vec3b(255, 255, 255))
				++white_pixel;

return white_pixel; }

void draw_simplex (Mat &image1, vector<Landmark> landmarks) {

	for (int it = 0; it < landmarks.size(); ++it) {
		for (int it2 = 0; it2 < landmarks.size(); ++it2) {
			if (it2 == it) continue;
			line (image1, landmarks[it].coord, landmarks[it2].coord, Scalar (0, 0, 0), 1, 1);
		}
	}

	for (int it = 0; it < landmarks.size(); ++it) {
		for (int it2 = 0; it2 < landmarks.size(); ++it2) {
			if (it2 == it) continue;
			for (int it3 = 0; it3 < landmarks.size(); ++it3) {
				if (it3 == it2 || it3 == it) continue;
				vector<Point> elm;
				elm.push_back (landmarks[it].coord);
				elm.push_back (landmarks[it2].coord);
				elm.push_back (landmarks[it3].coord);
				const Point* ppt[1] = { &elm[0] };
				int npt[] = { 3 };
				Mat image2 = image1.clone();
				fillPoly (image1, ppt, npt, 1, Scalar (rand()%256, rand()%256, rand()%256), 8);
				image1 = 0.6*image2 + 0.4*image1;
			}
		}
	}
}

void draw_path (Mat &PathIMG, int id, vector<Landmark> PthLndmrks, Scalar S, int thickness = 1) {

	if (PthLndmrks.size() > 1) {
		for (int it = 0; it < PthLndmrks.size()-1; ++it)
			line (PathIMG, PthLndmrks[it].coord, PthLndmrks[it+1].coord, S, thickness, 1);
	}

	else
		return;
}

void draw_falseHoles (Mat &falseHolesIMG, Scalar S, int thickness = 1) {

	for (int it = 0; it < clusters.size(); ++it) {
		for (int it2 = 0; it2 < clusters[it].size(); ++it2) {
			for (int it3 = 0; it3 < clusters[it].size(); ++it3) {
				line (falseHolesIMG, clusters[it][it2].coord,
						clusters[it][it3].coord, S, thickness, 1);
			}
		}
	}
}

void draw_homology_image (vector<int> assignment) {

	if (assignment.size() > 0) {
		for (int id = 0; id < rbt; ++id) {
			int AsgnClstr = assignment[id];
			int Nearest_Landmark_id = observations[id].nearest_landmark_ids[AsgnClstr];
			Landmark Nearest_Landmark = landmarks[Nearest_Landmark_id];

			circle (falseHolesIMG, observations[id].coord, 2, Scalar(100, 50*id, 0), 2);
			draw_robot (falseHolesIMG, observations[id]);
			vector<Landmark> asgn_clster = clusters[assignment[id]];
			int rand_it = ((((double)rand())/((double)RAND_MAX))*(asgn_clster.size()));

			Navigation HIW_exploration(Nearest_Landmark);
			HIW_exploration.start_landmarks = &observations[id].observing_landmarks;
			HIW_exploration.search();

			vector<Landmark> HIW_path_landmarks;
			std::vector<Landmark*> HIW_path =
				HIW_exploration.reconstructPointerPath(Nearest_Landmark);
			for (int p2 = HIW_path.size()-1; p2 >= 0; --p2) {
				HIW_path_landmarks.push_back (Landmark (
					(*HIW_path[p2]).coord.x, (*HIW_path[p2]).coord.y, (*HIW_path[p2]).id));
			}

			draw_path (falseHolesIMG, id, HIW_path_landmarks, Scalar(128, 0, 0), 2);
			draw_falseHoles (falseHolesIMG, Scalar(0, 0, 128), 1);

	//		for (int i = 0; i < asgn_clster.size(); ++i)
	//			disp (homology_image, asgn_clster[i].coord, Vec3b(0, 0, 255));
		}

		imshow("falseHolesIMG", falseHolesIMG);
		waitKey (0);
	}
}

void identify_clusters () {

	unordered_set<int> cc_ids_copy = cc_ids;
	clusters.clear();
	clusters_ids.clear();
	while (cc_ids_copy.size() > 0) {
		connected_components cc;
		cc.start_landmark = &landmarks[*cc_ids_copy.begin()];
		cc.search();

		clusters.push_back(cc.branch_landmarks);
		clusters_ids.push_back(cc.branch_landmarks_ids);

		for (int it = 0; it < cc.branch_landmarks.size(); ++it)
			cc_ids_copy.erase(cc.branch_landmarks[it].id);

	}

	clusters_img = nice_image.clone();

	for (int it = 0; it < clusters.size(); ++it)
		for (int it2 = 0; it2 < clusters[it].size(); ++it2)
			disp (clusters_img, clusters[it][it2].coord, Vec3b(it*50, it*20, it*10), 1);

//	imshow ("clusters", clusters_img);
//	waitKey (2e3);
}

void draw_vector (arma::Mat<double> x, Mat & img) {

	arma::Mat<double> xx;
	for (int it = 0; it < x.n_elem; ++it)
		xx = x / x.max();

	for (int it = 0; it < xx.n_elem; ++it) {
		vector<int> C0_ids;
		for(auto p = sc.id_to_simplex[it].S.begin();
		p != sc.id_to_simplex[it].S.end(); ++p)
			C0_ids.push_back(*p);

		if (xx(it,0) > 0) {
			line (img, landmarks[C0_ids[0]].coord,
				landmarks[C0_ids[1]].coord,
				Scalar (255, 255*(1-xx(it,0)), 255*(1-xx(it,0))), 1, 1);
		}

		if (xx(it,0) <= 0) {
			line (img, landmarks[C0_ids[0]].coord,
				landmarks[C0_ids[1]].coord,
				Scalar (255*(1-abs(xx(it,0))), 255*(1-abs(xx(it,0))), 255), 1, 1);
		}
	}
}

void assign_directions (arma::Mat<double> & x, unordered_set<int> edge_ids) {

	vector<simplex<int>> C1_simplices;
	for (int it = 0; it < x.n_elem; ++it) {
		if (edge_ids.find(it) != edge_ids.end()) {
			C1_simplices.push_back(sc.id_to_simplex[it]);
		}
	}

	vector<int> C0_ids;
	int remain_id;
	for (auto k = C1_simplices[0].S.begin(); k != C1_simplices[0].S.end(); ++k) 
		C0_ids.push_back(*k);

	if (C0_ids[0] > C0_ids[1]) {
		x(C1_simplices[0].id,0) = 1;
		remain_id = C0_ids[0];
	}

	if (C0_ids[1] > C0_ids[0]) {
		x(C1_simplices[0].id,0) = 1;
		remain_id = C0_ids[1];
	}

	C1_simplices.erase(C1_simplices.begin());

	while (C1_simplices.size() > 0) {
		int pos(0);
		for (int it = 0; it < C1_simplices.size(); ++it) {
			if (C1_simplices[it].S.find(remain_id) != C1_simplices[it].S.end()) {
				pos = it;
				C0_ids.clear();
				for (auto p = C1_simplices[it].S.begin(); p != C1_simplices[it].S.end(); ++p)
					if (*p != remain_id)
						C0_ids.push_back(*p);

				if (C0_ids[0] > remain_id) {
					x(C1_simplices[it].id,0) = 1;
					remain_id = C0_ids[0];
				}

				if (C0_ids[0] < remain_id) {
					x(C1_simplices[it].id,0) = -1;
					remain_id = C0_ids[0];
				}
			}
		}

		C1_simplices.erase(C1_simplices.begin()+pos);
	}
}
